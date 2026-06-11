%% =========================================================================
%  STUDY_GAMMA_SENSITIVITY.m
%  Detectability study set for gamma-analysis sensitivity to the CBCT1-vs-CBCT3
%  dose difference, on a SINGLE beam segment.
%
%  Headline number is the gamma pass rate, but every claim is anchored to a
%  discriminability metric (d' / ROC AUC) between a "no real change" (null)
%  distribution and a "real change" (signal) distribution. Two questions are
%  kept strictly separate:
%     (1) DETECTION  - can we tell something changed vs. measurement noise?
%     (2) FIDELITY   - does the detected difference match the true dose change?
%
%  Truth = dose-to-tissue from RayStation (rs_dose), per CT image, as loaded by
%  load_recon_dose_data. ETHOS truth (ethos_truth) is a single global field that
%  common-mode-cancels in the differential, so it is used only as an auxiliary
%  statistic.
%
%  THREE studies (each mirrors the existing pipeline's structure):
%    1. study_noise_null      - detection baseline. Compare a CT_1 dose to
%                               ITSELF under two independent noise seeds
%                               (recon[bundle1,seedA] vs recon[bundle1,seedC]),
%                               NOT CT_1 vs CT_3. Null gamma pass-rate band with
%                               error bars, overlaid with the signal band.
%    2. study_recon_fidelity  - fidelity floor. Per-scenario error R_s - D_s and
%                               the differential (R3-D3)-(R1-D1) = the part that
%                               does NOT common-mode-cancel (reconstruction
%                               artifact). Error maps + summary stats.
%    3. study_detectability_curve - gamma pass rate vs criterion (n%/n mm) with
%                               three series: truth CT1-vs-CT3, recon CT1-vs-CT3,
%                               and the null pass rate.
%
%  Criteria selection (post-hoc, cheap): gamma on already-reconstructed doses is
%  free, so the full pct x dta grid is swept once and the criterion that
%  MAXIMIZES d'/AUC between null and signal (not the raw pass-rate gap) is
%  reported, plus the n%-fixed-mm and n-mm-fixed-% slices.
%
%  The expensive k-Wave forward simulation runs ONCE per scenario (CT_1, CT_3).
%  Only the noise draw + Wiener deconvolution + time-reversal reconstruction are
%  repeated across the ensemble, distributed across the available GPUs.
%
%  Outputs (keyed by the recon CONFIG hash, hash8):
%    AnalysisResults/<pid>/<session>/gamma_sensitivity/<hash8>/
%        detectability_ensemble_<hash8>.mat   (raw ensemble; reused if present)
%        gamma_sensitivity_results_<hash8>.mat (all study outputs)
%        figures/*.png                         (every plot)
%
%  Follows CLAUDE.md and CLAUDE-SIMULATION_CONTEXT.md. Build on
%  study_gamma_pass_rates_chart.m, run_standalone_analysis.m,
%  load_recon_dose_data.m and CalcGamma.m.
%
%  NOTE: HIPAA / remote-execution - this file is WRITTEN here but must be RUN on
%  the remote device. Do not execute locally.
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';
CONFIG.treatment_site = 'Pancreas';
CONFIG.gruneisen_method = 'threshold_2';

% --- Single beam segment under study (the CBCT1-vs-CBCT3 pair) -------------
% Leave beam empty to auto-pick the first beam/segment present on BOTH CTs.
CONFIG.beam    = [];        % scalar beam number, or [] to auto-pick a paired field
CONFIG.segment = [];        % scalar segment number, or [] for no segment filter
CONFIG.ct_pair = [1, 3];    % CT labels compared: ct_pair(1) = reference arm

% Explicit recon config-hash override ('' => auto-discover on disk via loader).
CONFIG.config_hash = '';

% --- Gamma sweep -----------------------------------------------------------
% Diagonal sweep (n%/n mm) used by the null band and the detectability curve.
CONFIG.gamma_sweep_n = (2:1:10)';     % criterion n: (2,2) .. (10,10)
% Full grid (percent x DTA) used by the post-hoc criteria selection. Must
% include every value in gamma_sweep_n so the diagonal is a sub-grid of it.
CONFIG.gamma_grid_pct = (2:1:10)';    % percent axis
CONFIG.gamma_grid_dta = (2:1:10)';    % DTA (mm) axis
CONFIG.low_dose_cutoff_frac = 0.10;   % 10% of reference max -> gamma eval mask
CONFIG.selection_metric = 'auc';      % 'auc' | 'dprime' for the criterion pick

% --- Noise ensemble --------------------------------------------------------
% Per realization: 3 reconstructions (recon[b1,seedA], recon[b1,seedC],
% recon[b2,seedB]). The CT_1 reference arm recon[b1,seedA] is SHARED by both the
% null and signal comparisons so the only difference between the two
% distributions is whether the target is a fresh CT_1 recon (null) or a CT_3
% recon (signal).
CONFIG.ensemble_N         = 8;        % number of realizations
CONFIG.ensemble_base_seed = 42;       % RNG base; per-realization seeds derive from it
CONFIG.recompute_ensemble = true;     % false -> load saved ensemble if present

% --- Tissue model (create_medium) -----------------------------------------
CONFIG.force_uniform_density     = false;
CONFIG.force_uniform_sound_speed = false;
CONFIG.force_uniform_attenuation = false;
CONFIG.force_uniform_gruneisen   = false;
CONFIG.uniform_density      = 1000;
CONFIG.uniform_sound_speed  = 1540;
CONFIG.uniform_alpha_coeff  = 0;
CONFIG.uniform_alpha_power  = 1.1;
CONFIG.uniform_gruneisen    = 1.0;

% --- Sensor geometry (determine_sensor_mask) ------------------------------
CONFIG.sensor_placement_method = 'determine_sensor_mask';
CONFIG.sensor_x_index    = 2;
CONFIG.sensor_y_index    = 4;
CONFIG.elements_per_side = 32;
CONFIG.element_pitch_mm  = 4.35;     % kerf derived inside determine_sensor_mask as (pitch - size)
CONFIG.element_size_mm   = 3.65;
CONFIG.sensor_standoff_mm = 5;
CONFIG.jaw_margin_mm      = 10;
CONFIG.sensor_placement   = 'anterior';
CONFIG.aim_at_iso         = true;
CONFIG.force_turn_angle   = 290;     % forced turn must be 290 deg (rotation-bug workaround)

% --- Physics / forward + reconstruction -----------------------------------
CONFIG.dose_per_pulse_cGy = 0.16;
CONFIG.meterset           = 140;     % overridden per-scenario from rtplan.meterset when available
CONFIG.pml_size           = 10;
CONFIG.cfl_number         = 0.3;
CONFIG.use_gpu            = true;
CONFIG.correction_factor  = 1.9;
CONFIG.use_pressure_scale_correction = false;
CONFIG.mask_recon_to_dose_region     = true;

CONFIG.reconstruction_method  = 'tr';   % 'tr' | 'DAS' | 'hybrid'
CONFIG.num_time_reversal_iter = 10;
CONFIG.convergence_tol        = 1e-3;

% --- Pulse convolution / noise / deconvolution ----------------------------
% Electronic noise enters here; the ensemble is meaningless without it, so
% convolution_kernel MUST be > 0.
CONFIG.convolution_kernel  = 4e-6;    % Gaussian sigma in seconds (4 us)
CONFIG.conv_noise_level    = 0.125;   % noise amplitude as fraction of peak sensor signal
CONFIG.conv_deconv_lambda  = 1e-4;    % Wiener regularization

CONFIG.downscale_factor = 1;
CONFIG.use_grid_padding = true;

% --- Visualization ---------------------------------------------------------
CONFIG.viz_smooth_sigma = 1.0;        % display-only smoothing (voxels); 0 disables

assert(CONFIG.convolution_kernel > 0, ...
    'study_gamma_sensitivity:NoNoise', ...
    ['Detectability study requires CONFIG.convolution_kernel > 0 (electronic ' ...
     'noise is injected in the pulse model).']);

%% ========================= PRINT CONFIGURATION ===========================

fprintf('=========================================================\n');
fprintf('  Gamma Detectability Sensitivity Study\n');
fprintf('=========================================================\n');
fprintf('  Patient:        %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Tissue model:   %s\n', CONFIG.gruneisen_method);
fprintf('  Recon method:   %s (%d iter)\n', CONFIG.reconstruction_method, CONFIG.num_time_reversal_iter);
fprintf('  CT pair:        CT_%d (ref) vs CT_%d (signal)\n', CONFIG.ct_pair(1), CONFIG.ct_pair(2));
fprintf('  Ensemble N:     %d (base seed %d)\n', CONFIG.ensemble_N, CONFIG.ensemble_base_seed);
fprintf('  Gamma diagonal: %g..%g %%/mm\n', CONFIG.gamma_sweep_n(1), CONFIG.gamma_sweep_n(end));
fprintf('  Selection:      %s\n', CONFIG.selection_metric);
fprintf('=========================================================\n\n');

%% ===================== LOAD THE CT_1 / CT_3 PAIR =========================
% Load the single beam segment on BOTH CT images. Mode='set' does not require a
% unique match, so the same beam/segment present on CT_1 and CT_3 loads as two
% fields; each is selected by its own CT label.

ct1 = sprintf('CT_%d', CONFIG.ct_pair(1));
ct2 = sprintf('CT_%d', CONFIG.ct_pair(2));

if isempty(CONFIG.beam)
    [beam, seg] = pick_first_paired_field(CONFIG, ct1, ct2);
    fprintf('[Load] Auto-picked paired field: beam %d, segment %s\n', ...
        beam, mat2str(seg));
else
    beam = CONFIG.beam;
    seg  = CONFIG.segment;
end

out = load_field_set(CONFIG, beam, seg);

iA = find_field_by_ct(out.fields, ct1);
iB = find_field_by_ct(out.fields, ct2);
if isempty(iA) || isempty(iB)
    found = strjoin(arrayfun(@(f) field_ct_label(f), out.fields, ...
        'UniformOutput', false), ', ');
    error('study_gamma_sensitivity:NoPair', ...
        ['Beam %s / segment %s does not have both %s and %s (found CT labels: %s).'], ...
        mat2str(beam), mat2str(seg), ct1, ct2, found);
end
fldA = out.fields(iA);   % CT_1 (reference arm)
fldB = out.fields(iB);   % CT_3 (signal arm)

metadata   = out.metadata;
spacing_mm = metadata.spacing(:)';
hash8      = out.config_hash;

beam_meta = [];
if isfield(metadata, 'beam_metadata') && ~isempty(metadata.beam_metadata)
    beam_meta = metadata.beam_metadata;
end

% Truth doses (dose-to-tissue from RayStation), one per CT image.
D1_truth = double(fldA.rs_dose);
D3_truth = double(fldB.rs_dose);

% On-disk pipeline reconstructions (canonical output) kept as a faint reference
% line on the detectability curve.
recon_ondisk_A = [];
recon_ondisk_B = [];
if isfield(fldA, 'recon_dose') && ~isempty(fldA.recon_dose)
    recon_ondisk_A = double(fldA.recon_dose);
end
if isfield(fldB, 'recon_dose') && ~isempty(fldB.recon_dose)
    recon_ondisk_B = double(fldB.recon_dose);
end

% Optional ETHOS global truth (auxiliary only).
ethos_truth = [];
if isfield(out, 'ethos_truth') && ~isempty(out.ethos_truth)
    ethos_truth = double(out.ethos_truth);
end

gantry_A = getfield_def(fldA.rtplan, 'gantry_angle', 0);
gantry_B = getfield_def(fldB.rtplan, 'gantry_angle', 0);
meterset_A = getfield_def(fldA.rtplan, 'meterset', CONFIG.meterset);
meterset_B = getfield_def(fldB.rtplan, 'meterset', CONFIG.meterset);

fprintf('[Load] Grid [%d x %d x %d], spacing [%.2f %.2f %.2f] mm, hash %s\n', ...
    size(D1_truth,1), size(D1_truth,2), size(D1_truth,3), ...
    spacing_mm(1), spacing_mm(2), spacing_mm(3), hash8);
fprintf('       CT_%d: gantry %.1f deg, meterset %.2f MU\n', CONFIG.ct_pair(1), gantry_A, meterset_A);
fprintf('       CT_%d: gantry %.1f deg, meterset %.2f MU\n', CONFIG.ct_pair(2), gantry_B, meterset_B);

%% ===================== OUTPUT PATHS ======================================

out_dir   = fullfile(CONFIG.working_dir, 'AnalysisResults', CONFIG.patient_id, ...
    CONFIG.session, 'gamma_sensitivity', hash8);
plots_dir = fullfile(out_dir, 'figures');
if ~isfolder(out_dir),   mkdir(out_dir);   end
if ~isfolder(plots_dir), mkdir(plots_dir); end
fprintf('[Out]  %s\n', out_dir);

% Deterministic sensor placement: use the summed plan dose so the array
% geometry (and any grid expansion) is IDENTICAL for CT_1 and CT_3.
processed_dir   = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
    CONFIG.patient_id, CONFIG.session, 'processed');
total_dose_file = fullfile(processed_dir, 'total_rs_dose.mat');
if isfile(total_dose_file)
    CONFIG.total_dose_file = total_dose_file;
    fprintf('[Out]  Sensor placement dose: total_rs_dose.mat (deterministic)\n');
else
    CONFIG.total_dose_file = '';
    fprintf('[Out]  [WARN] total_rs_dose.mat missing; sensor placement uses per-field dose.\n');
end

%% ===================== BUILD / LOAD ENSEMBLE =============================
% The forward simulation is the expensive part: it runs ONCE per scenario when
% the ensemble is (re)computed, then the noise + reconstruction loop runs N
% times across the GPUs. The raw ensemble is cached so the (cheap) studies and
% plots can be re-run without re-simulating.

ens_file = fullfile(out_dir, sprintf('detectability_ensemble_%s.mat', hash8));

if ~CONFIG.recompute_ensemble && isfile(ens_file)
    fprintf('\n[Ensemble] Loading cached ensemble from %s\n', ens_file);
    L = load(ens_file, 'ENS');
    ENS = L.ENS;
else
    fprintf('\n[Ensemble] Building forward bundles (one k-Wave forward per scenario)...\n');

    cfgA = CONFIG; cfgA.meterset = meterset_A;
    cfgB = CONFIG; cfgB.meterset = meterset_B;

    bundle1 = build_forward_bundle(D1_truth, fldA.cbct, gantry_A, beam_meta, cfgA, ...
        sprintf('CT_%d', CONFIG.ct_pair(1)));
    bundle2 = build_forward_bundle(D3_truth, fldB.cbct, gantry_B, beam_meta, cfgB, ...
        sprintf('CT_%d', CONFIG.ct_pair(2)));

    if ~isequal(bundle1.gridSize_orig, bundle2.gridSize_orig)
        error('study_gamma_sensitivity:GridMismatch', ...
            ['The two scenario bundles ended on different grids ([%s] vs [%s]).\n' ...
             'They must share a grid for the gamma comparison and recon averaging.\n' ...
             'Ensure deterministic sensor placement (total_rs_dose.mat present).'], ...
            num2str(bundle1.gridSize_orig), num2str(bundle2.gridSize_orig));
    end

    ENS = run_detectability_ensemble(bundle1, bundle2, CONFIG);

    % Attach context needed by the studies / plots.
    ENS.spacing_mm     = spacing_mm;
    ENS.hash8          = hash8;
    ENS.patient_id     = CONFIG.patient_id;
    ENS.session        = CONFIG.session;
    ENS.ct_pair        = CONFIG.ct_pair;
    ENS.D1             = bundle1.truth_dose;     % truth on the recon grid
    ENS.D3             = bundle2.truth_dose;
    ENS.doseMask       = logical(bundle1.doseMask) | logical(bundle2.doseMask);
    ENS.gamma_sweep_n  = CONFIG.gamma_sweep_n(:);
    ENS.pct_vec        = CONFIG.gamma_grid_pct(:);
    ENS.dta_vec        = CONFIG.gamma_grid_dta(:);

    % On-disk recons / ETHOS truth embedded onto the recon grid (alignment-aware
    % when sensor grid-expansion shifted the content).
    ENS.recon_ondisk_A = embed_to_recon_grid(recon_ondisk_A, bundle1);
    ENS.recon_ondisk_B = embed_to_recon_grid(recon_ondisk_B, bundle1);
    ENS.ethos_truth    = embed_to_recon_grid(ethos_truth,    bundle1);

    save(ens_file, 'ENS', '-v7.3');
    fprintf('[Ensemble] Saved raw ensemble to %s\n', ens_file);
end

%% ===================== RUN THE THREE STUDIES ============================

fprintf('\n[Study 1/3] Noise null (detection baseline)...\n');
S1 = study_noise_null(ENS, CONFIG, plots_dir);

fprintf('\n[Study 2/3] Reconstruction fidelity (differential artifact)...\n');
S2 = study_recon_fidelity(ENS, CONFIG, plots_dir);

fprintf('\n[Study 3/3] Detectability curve...\n');
S3 = study_detectability_curve(ENS, CONFIG, plots_dir);

%% ===================== POST-HOC CRITERIA SELECTION ======================

fprintf('\n[Select] Post-hoc criterion selection (%s over the full grid)...\n', ...
    CONFIG.selection_metric);
SEL = select_criteria(ENS.null_pr, ENS.signal_pr, ENS.pct_vec, ENS.dta_vec, ...
    CONFIG.selection_metric);
plot_criteria_selection(SEL, ENS, CONFIG, plots_dir);

fprintf('  Best criterion (%s): %g%% / %g mm   d''=%.3f  AUC=%.3f\n', ...
    CONFIG.selection_metric, SEL.best_pct, SEL.best_dta, SEL.best_dprime, SEL.best_auc);

%% ===================== SAVE COMBINED RESULTS ============================

results = struct();
results.config        = CONFIG;
results.hash8         = hash8;
results.spacing_mm    = spacing_mm;
results.beam          = beam;
results.segment       = seg;
results.noise_null    = S1;
results.fidelity      = S2;
results.detectability = S3;
results.selection     = SEL;
results.timestamp     = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>

results_file = fullfile(out_dir, sprintf('gamma_sensitivity_results_%s.mat', hash8));
save(results_file, 'results', '-v7.3');
fprintf('\n[Done] Saved combined results to %s\n', results_file);
fprintf('[Done] Figures in %s\n', plots_dir);


% =========================================================================
% =============================  LOCAL FUNCTIONS  =========================
% =========================================================================

% -------------------------------------------------------------------------
%  LOADING HELPERS (mirrors run_standalone_analysis.m)
% -------------------------------------------------------------------------

function out = load_field_set(CONFIG, beam, seg)
%LOAD_FIELD_SET Load matched fields via load_recon_dose_data (Mode='set').
%  IncludeEthos=true so the ETHOS global truth is available as an auxiliary stat.
    args = {'Mode', 'set', 'IncludeEthos', true};
    if ~isempty(beam), args = [args, {'Beam', beam}]; end
    if ~isempty(seg),  args = [args, {'Segment', seg}]; end
    if isfield(CONFIG, 'config_hash') && ~isempty(CONFIG.config_hash)
        args = [args, {'Hash', CONFIG.config_hash}];
    end
    out = load_recon_dose_data(CONFIG.patient_id, CONFIG.session, CONFIG, args{:});
end

function [beam, seg] = pick_first_paired_field(CONFIG, ct1, ct2)
%PICK_FIRST_PAIRED_FIELD First beam/segment present on BOTH ct1 and ct2.
    ldcfg.working_dir = CONFIG.working_dir;
    fidx = list_processed_field_doses(CONFIG.patient_id, CONFIG.session, ldcfg);

    n = numel(fidx);
    keys = strings(1, n);
    cts  = strings(1, n);
    for i = 1:n
        keys(i) = sprintf('%d|%d', fidx(i).beam_index, fidx(i).segment);
        cts(i)  = string(ct_label_from_name(fidx(i).source_mat_filename));
    end

    seen = unique(keys, 'stable');
    for k = 1:numel(seen)
        m = keys == seen(k);
        if any(cts(m) == ct1) && any(cts(m) == ct2)
            idx  = find(m, 1);
            beam = fidx(idx).beam_index;
            seg  = fidx(idx).segment;
            return;
        end
    end

    error('study_gamma_sensitivity:NoPairedField', ...
        ['No field is reconstructed on both %s and %s. Set CONFIG.beam / ' ...
         'CONFIG.segment, or adjust CONFIG.ct_pair.'], ct1, ct2);
end

function idx = find_field_by_ct(fields, want_ct)
%FIND_FIELD_BY_CT Index of the first loaded field whose CT label is want_ct.
    idx = [];
    for i = 1:numel(fields)
        if strcmpi(field_ct_label(fields(i)), want_ct)
            idx = i;
            return;
        end
    end
end

function lbl = field_ct_label(fld)
%FIELD_CT_LABEL CT label of a loaded field: rtplan.ct_label, else filename.
    lbl = '';
    if isfield(fld, 'rtplan') && isfield(fld.rtplan, 'ct_label') ...
            && ~isempty(fld.rtplan.ct_label)
        lbl = strrep(char(fld.rtplan.ct_label), '-', '_');
    end
    if isempty(lbl) && isfield(fld, 'source_mat_filename')
        lbl = ct_label_from_name(fld.source_mat_filename);
    end
    if isempty(lbl), lbl = 'unknown'; end
end

function ct_label = ct_label_from_name(name)
%CT_LABEL_FROM_NAME Lenient CT-label parse: matches CT_1, CT-1, CT1, ...
    ct_label = '';
    tok = regexp(char(name), 'CT[_-]?(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        ct_label = sprintf('CT_%s', tok{1});
    end
end

function v = getfield_def(s, f, default)
%GETFIELD_DEF Struct field value, or default when missing/empty/NaN.
    v = default;
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        candidate = s.(f);
        if isnumeric(candidate) && all(isnan(candidate(:)))
            return;
        end
        v = candidate;
    end
end

% -------------------------------------------------------------------------
%  FORWARD BUNDLE  (faithful port of study_gamma_pass_rates_chart.m forward
%  pipeline, lines ~438-1011, minus plotting). Runs ONE k-Wave forward sim and
%  captures everything needed to re-reconstruct from a fresh noise draw.
% -------------------------------------------------------------------------

function B = build_forward_bundle(truthDose, sct, gantry_angle, beam_meta, CONFIG, label)
%BUILD_FORWARD_BUNDLE  Dose -> medium -> p0 -> sensor -> k-Wave forward -> pulse
%  model, returning the reconstruction bundle B (see reconstruct_recon_dose).
%  B.truth_dose is the (body-masked, possibly grid-expanded) dose on the recon
%  grid, so it aligns voxel-for-voxel with any recon produced from B.

    fprintf('  [Bundle %s] forward simulation...\n', label);

    doseGrid = double(truthDose);
    if ~isfield(sct, 'cubeHU')
        error('build_forward_bundle:NoHU', 'CBCT struct missing cubeHU.');
    end
    spacing_mm = sct.spacing(:)';
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;

    gridSize = size(doseGrid);
    Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);

    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    % --- Optional downscaling ---
    if CONFIG.downscale_factor ~= 1
        df     = CONFIG.downscale_factor;
        new_Nx = max(1, round(Nx / df));
        new_Ny = max(1, round(Ny / df));
        new_Nz = max(1, round(Nz / df));
        doseGrid   = max(imresize3(doseGrid, [new_Nx, new_Ny, new_Nz]), 0);
        sct.cubeHU = imresize3(sct.cubeHU, [new_Nx, new_Ny, new_Nz]);
        if isfield(sct, 'bodyMask')
            sct.bodyMask = imresize3(single(sct.bodyMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        if isfield(sct, 'couchMask') && ~isempty(sct.couchMask)
            sct.couchMask = imresize3(single(sct.couchMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        spacing_mm = spacing_mm .* ([Nx, Ny, Nz] ./ [new_Nx, new_Ny, new_Nz]);
        dx = spacing_mm(1) / 1000; dy = spacing_mm(2) / 1000; dz = spacing_mm(3) / 1000;
        sct.spacing = spacing_mm;
        Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
        gridSize = [Nx, Ny, Nz];
    end

    % --- Acoustic medium ---
    medium = create_medium(sct, CONFIG);

    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    % --- Initial pressure p0 ---
    meterset       = CONFIG.meterset;
    num_pulses     = ceil(meterset / CONFIG.dose_per_pulse_cGy);
    dose_per_pulse = doseGrid / num_pulses;
    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    p0 = smooth(p0);

    doseThreshold = 0.01 * max(doseGrid(:));
    doseMask      = doseGrid > doseThreshold;
    if ~any(doseMask(:)) || max(p0(:)) == 0
        error('build_forward_bundle:NoDose', 'No significant dose / zero p0 for %s.', label);
    end

    % --- FFT-optimal grid padding ---
    Nx_orig = Nx; Ny_orig = Ny; Nz_orig = Nz;
    gridSize_orig = gridSize;
    medium_orig   = medium;

    if CONFIG.use_grid_padding
        Nx_pad = find_optimal_kwave_size(Nx, CONFIG.pml_size);
        Ny_pad = find_optimal_kwave_size(Ny, CONFIG.pml_size);
        Nz_pad = find_optimal_kwave_size(Nz, CONFIG.pml_size);
    else
        Nx_pad = Nx; Ny_pad = Ny; Nz_pad = Nz;
    end

    did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
    if did_pad
        [medium, p0] = pad_medium_p0(medium, p0, Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);
        Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
        gridSize = [Nx, Ny, Nz];
    end

    % Track the pre-expansion recon-grid dims and the offset that sensor
    % grid-expansion introduces, so external-grid arrays (ETHOS truth, on-disk
    % recon) can be embedded into the recon grid in alignment.
    dims_pre_expand = [Nx_orig, Ny_orig, Nz_orig];
    embed_offset    = [0, 0, 0];

    % --- Sensor placement ---
    sensor      = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    switch CONFIG.sensor_placement_method
        case 'full_plane_anterior'
            sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
            sensor_info_orig = struct('num_elements', 0);

        case 'determine_sensor_mask'
            sct_for_sensor = sct;
            if ~isfield(sct_for_sensor, 'couchMask') || isempty(sct_for_sensor.couchMask)
                sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
            end
            if ~isfield(sct_for_sensor, 'origin') || isempty(sct_for_sensor.origin)
                sct_for_sensor.origin = [0, 0, 0];
            end
            sct_for_sensor.spacing = spacing_mm;

            field_dose_for_sensor             = struct();
            field_dose_for_sensor.dose_Gy     = doseGrid;
            field_dose_for_sensor.gantry_angle = gantry_angle;
            field_dose_for_sensor.origin      = sct_for_sensor.origin;
            field_dose_for_sensor.spacing     = spacing_mm;
            field_dose_for_sensor.dimensions  = [Nx_orig, Ny_orig, Nz_orig];

            [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask( ...
                sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);

            % --- Grid expansion handling (water padding to clear the exclusion
            % zone), then re-run FFT-optimal padding. Faithful to the template.
            gp = sensor_info_orig.grid_pad;
            if gp.expanded
                [medium, p0, doseGrid, doseMask, sct, medium_orig, ...
                 Nx_orig, Ny_orig, Nz_orig, Nx, Ny, Nz] = ...
                    apply_grid_expansion(medium, p0, doseGrid, doseMask, sct, ...
                        gp, Nx_orig, Ny_orig, Nz_orig, CONFIG);
                gridSize_orig = [Nx_orig, Ny_orig, Nz_orig];
                gridSize      = [Nx, Ny, Nz];
                sensor.mask   = zeros(Nx, Ny, Nz);
                embed_offset  = [gp.y_pre, gp.x_pre, gp.z_pre];
            end

            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));

        otherwise
            error('build_forward_bundle:Sensor', ...
                'Unsupported sensor_placement_method: %s', CONFIG.sensor_placement_method);
    end

    if sum(sensor.mask(:)) == 0
        error('build_forward_bundle:EmptySensor', 'Sensor mask is empty for %s.', label);
    end

    % --- k-Wave grid & medium ---
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    maxC  = max(medium.sound_speed(:));
    minC  = min(medium.sound_speed(medium.sound_speed > 0));
    dt    = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;
    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);
    kgrid.dt = dt;
    kgrid.Nt = Nt;

    kmedium             = struct();
    kmedium.density     = medium.density;
    kmedium.sound_speed = medium.sound_speed;
    kmedium.alpha_coeff = medium.alpha_coeff;
    kmedium.alpha_power = 1.1;

    if CONFIG.use_gpu
        try
            gpuDevice;
            dataCast = 'gpuArray-single';
        catch
            dataCast = 'single';
        end
    else
        dataCast = 'single';
    end
    inputArgs = {'Smooth', false, 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
                 'DataCast', dataCast, 'PlotSim', false};

    % --- Forward simulation ---
    source_fwd    = struct();
    source_fwd.p0 = p0;
    sensorData    = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
    sensorData    = smooth(sensorData);
    FS            = 1 / kgrid.dt;

    % --- Pulse model: convolve, frequency response, noise, Wiener deconv ---
    conv_kernel_sigma  = CONFIG.convolution_kernel;
    conv_deconv_lambda = CONFIG.conv_deconv_lambda;

    sigma_samples = conv_kernel_sigma / dt;
    kernel_half   = ceil(4 * sigma_samples);
    t_kernel      = (-kernel_half : kernel_half)';
    gauss_kernel  = exp(-t_kernel.^2 / (2 * sigma_samples^2));
    gauss_kernel  = gauss_kernel / sum(gauss_kernel);

    sensorData_cpu = double(gather(sensorData));
    Nt_data        = size(sensorData_cpu, 2);

    H       = fft(gauss_kernel, Nt_data).';
    H_conj  = conj(H);
    H_power = abs(H).^2;

    sensorData_conv = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));
    sensorData_resp = gaussianFilter(sensorData_conv, FS, 0.35e6, 100, true);
    noise_amp       = CONFIG.conv_noise_level * max(abs(sensorData_resp(:)));

    fprintf('  [Bundle %s] forward done. Sensor [%d x %d], noise amp %.3e Pa\n', ...
        label, size(sensorData,1), size(sensorData,2), noise_amp);

    % --- Capture the bundle ---
    bodyMask = [];
    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), gridSize_orig)
        bodyMask = double(sct.bodyMask);
    end

    B = struct();
    B.kgrid                = kgrid;
    B.kmedium              = kmedium;
    B.sensor               = sensor;
    B.inputArgs            = inputArgs;     % cell of k-Wave name/value args
    B.p0                   = p0;
    B.gridSize             = gridSize;
    B.gridSize_orig        = gridSize_orig;
    B.did_pad              = did_pad;
    B.Nx_orig              = Nx_orig;
    B.Ny_orig              = Ny_orig;
    B.Nz_orig              = Nz_orig;
    B.dx = dx; B.dy = dy; B.dz = dz; B.dt = dt;
    B.medium_pad           = medium;
    B.medium_orig          = medium_orig;
    B.doseMask             = doseMask;
    B.num_pulses           = num_pulses;
    B.sensor_info_orig     = sensor_info_orig;
    B.das_opts             = struct();
    B.bodyMask             = bodyMask;
    B.reconstruction_method         = CONFIG.reconstruction_method;
    B.num_time_reversal_iter        = CONFIG.num_time_reversal_iter;
    B.convergence_tol               = CONFIG.convergence_tol;
    B.correction_factor             = CONFIG.correction_factor;
    B.use_pressure_scale_correction = CONFIG.use_pressure_scale_correction;
    B.mask_recon_to_dose_region     = CONFIG.mask_recon_to_dose_region;
    B.sensorData_resp      = double(gather(sensorData_resp));
    B.H_conj               = H_conj;
    B.H_power              = H_power;
    B.conv_deconv_lambda   = conv_deconv_lambda;
    B.noise_amp            = noise_amp;
    B.truth_dose           = doseGrid;     % on gridSize_orig, aligns with recon
    B.spacing_mm           = spacing_mm;
    B.gantry_angle         = gantry_angle;
    B.pre_expand_dims      = dims_pre_expand;
    B.embed_offset         = embed_offset;
    B.label                = label;
end

function out = embed_to_recon_grid(A, B)
%EMBED_TO_RECON_GRID  Place an external-grid array A (on B.pre_expand_dims) into
%  the recon grid B.gridSize_orig, applying the sensor grid-expansion offset so
%  it aligns voxel-for-voxel with reconstructions. [] passes through.
    if isempty(A)
        out = [];
        return;
    end
    gs = B.gridSize_orig;
    if isequal(size(A), gs)
        out = A;       % no expansion: already aligned
        return;
    end
    out = zeros(gs);
    d = B.pre_expand_dims; o = B.embed_offset;
    m1 = min(d(1), size(A, 1));
    m2 = min(d(2), size(A, 2));
    m3 = min(d(3), size(A, 3));
    out(o(1) + (1:m1), o(2) + (1:m2), o(3) + (1:m3)) = A(1:m1, 1:m2, 1:m3);
end

function [medium, p0] = pad_medium_p0(medium, p0, Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad)
%PAD_MEDIUM_P0 Zero/water-pad medium + p0 from [Nx Ny Nz] to the padded size.
    density_pad    = ones(Nx_pad, Ny_pad, Nz_pad) * 1000;
    soundSpeed_pad = ones(Nx_pad, Ny_pad, Nz_pad) * 1540;
    alphaCoeff_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
    gruneisen_pad  = zeros(Nx_pad, Ny_pad, Nz_pad);

    density_pad(1:Nx, 1:Ny, 1:Nz)    = medium.density;
    soundSpeed_pad(1:Nx, 1:Ny, 1:Nz) = medium.sound_speed;
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_pad(1:Nx, 1:Ny, 1:Nz) = medium.alpha_coeff;
    else
        alphaCoeff_pad(:) = medium.alpha_coeff;
    end
    gruneisen_pad(1:Nx, 1:Ny, 1:Nz) = medium.gruneisen;

    medium.density     = density_pad;
    medium.sound_speed = soundSpeed_pad;
    medium.alpha_coeff = alphaCoeff_pad;
    medium.gruneisen   = gruneisen_pad;

    p0_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
    p0_pad(1:Nx, 1:Ny, 1:Nz) = p0;
    p0 = p0_pad;
end

function [medium, p0, doseGrid, doseMask, sct, medium_orig, ...
          Nx_orig, Ny_orig, Nz_orig, Nx, Ny, Nz] = ...
        apply_grid_expansion(medium, p0, doseGrid, doseMask, sct, gp, ...
            Nx_orig, Ny_orig, Nz_orig, CONFIG)
%APPLY_GRID_EXPANSION Water-pad the medium/p0/dose to the sensor-expanded grid,
%  then re-run FFT-optimal padding. Faithful to study_gamma_pass_rates_chart.m.

    density_unp    = medium.density(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    soundSpeed_unp = medium.sound_speed(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_unp = medium.alpha_coeff(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    else
        alphaCoeff_unp = medium.alpha_coeff;
    end
    gruneisen_unp  = medium.gruneisen(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    p0_unp         = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

    % Caller dim order: dim1=Nx, dim2=Ny, dim3=Nz -> gp.y_*->dim1, x_*->dim2, z_*->dim3.
    Nx_exp = Nx_orig + gp.y_pre + gp.y_post;
    Ny_exp = Ny_orig + gp.x_pre + gp.x_post;
    Nz_exp = Nz_orig + gp.z_pre + gp.z_post;

    density_exp    = ones(Nx_exp, Ny_exp, Nz_exp)  * 1000;
    soundSpeed_exp = ones(Nx_exp, Ny_exp, Nz_exp)  * 1540;
    alphaCoeff_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
    gruneisen_exp  = zeros(Nx_exp, Ny_exp, Nz_exp);
    p0_exp         = zeros(Nx_exp, Ny_exp, Nz_exp);

    xr = gp.y_pre + (1:Nx_orig);
    yr = gp.x_pre + (1:Ny_orig);
    zr = gp.z_pre + (1:Nz_orig);

    density_exp(xr, yr, zr)    = density_unp;
    soundSpeed_exp(xr, yr, zr) = soundSpeed_unp;
    if numel(alphaCoeff_unp) > 1
        alphaCoeff_exp(xr, yr, zr) = alphaCoeff_unp;
    else
        alphaCoeff_exp(:) = alphaCoeff_unp;
    end
    gruneisen_exp(xr, yr, zr)  = gruneisen_unp;
    p0_exp(xr, yr, zr)         = p0_unp;

    medium.density     = density_exp;
    medium.sound_speed = soundSpeed_exp;
    medium.alpha_coeff = alphaCoeff_exp;
    medium.gruneisen   = gruneisen_exp;
    p0 = p0_exp;

    medium_orig = struct( ...
        'density',     density_exp, ...
        'sound_speed', soundSpeed_exp, ...
        'alpha_coeff', alphaCoeff_exp, ...
        'gruneisen',   gruneisen_exp);

    doseGrid_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
    doseGrid_exp(xr, yr, zr) = doseGrid;
    doseGrid = doseGrid_exp;

    doseMask_exp = false(Nx_exp, Ny_exp, Nz_exp);
    doseMask_exp(xr, yr, zr) = doseMask;
    doseMask = doseMask_exp;

    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), [Nx_orig, Ny_orig, Nz_orig])
        body_exp = false(Nx_exp, Ny_exp, Nz_exp);
        body_exp(xr, yr, zr) = sct.bodyMask;
        sct.bodyMask = body_exp;
    end
    if isfield(sct, 'couchMask') && ~isempty(sct.couchMask) && ...
            isequal(size(sct.couchMask), [Nx_orig, Ny_orig, Nz_orig])
        couch_exp = false(Nx_exp, Ny_exp, Nz_exp);
        couch_exp(xr, yr, zr) = sct.couchMask;
        sct.couchMask = couch_exp;
    end

    Nx_orig = Nx_exp; Ny_orig = Ny_exp; Nz_orig = Nz_exp;

    if CONFIG.use_grid_padding
        Nx_pad2 = find_optimal_kwave_size(Nx_orig, CONFIG.pml_size);
        Ny_pad2 = find_optimal_kwave_size(Ny_orig, CONFIG.pml_size);
        Nz_pad2 = find_optimal_kwave_size(Nz_orig, CONFIG.pml_size);
    else
        Nx_pad2 = Nx_orig; Ny_pad2 = Ny_orig; Nz_pad2 = Nz_orig;
    end

    if ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], [Nx_orig, Ny_orig, Nz_orig])
        [medium, p0] = pad_medium_p0(medium, p0, Nx_orig, Ny_orig, Nz_orig, ...
            Nx_pad2, Ny_pad2, Nz_pad2);
    end

    Nx = Nx_pad2; Ny = Ny_pad2; Nz = Nz_pad2;
end


% -------------------------------------------------------------------------
%  DETECTABILITY ENSEMBLE  (3 reconstructions per realization, GPU-pinned)
% -------------------------------------------------------------------------

function ENS = run_detectability_ensemble(bundle1, bundle2, CONFIG)
%RUN_DETECTABILITY_ENSEMBLE  For each realization r, reconstruct:
%    recA1 = recon(bundle1, seedA)   reference arm (SHARED by null & signal)
%    recA2 = recon(bundle1, seedC)   null target   (same scenario, fresh noise)
%    recB  = recon(bundle2, seedB)   signal target (CT_3 scenario)
%  and record the full pct x dta gamma pass-rate grid for both
%  null   = gamma(recA1, recA2)   and   signal = gamma(recA1, recB).
%  The ensemble-mean recons (R1_mean, R3_mean) feed the fidelity study.

    N         = CONFIG.ensemble_N;
    base_seed = CONFIG.ensemble_base_seed;
    pct_vec   = CONFIG.gamma_grid_pct(:);
    dta_vec   = CONFIG.gamma_grid_dta(:);
    cutoff_f  = CONFIG.low_dose_cutoff_frac;
    Kp = numel(pct_vec); Kd = numel(dta_vec);
    spacing_mm = bundle1.spacing_mm;
    gridc      = bundle1.gridSize_orig;

    % --- Multi-GPU pool: one worker per GPU, each pinned to a device ---
    ngpu = 0;
    try, ngpu = gpuDeviceCount; catch, ngpu = 0; end
    desired_workers = max(1, ngpu);
    pool = gcp('nocreate');
    if isempty(pool) && exist('parpool', 'file') == 2
        try
            parpool('local', desired_workers);
            pool = gcp('nocreate');
        catch ME
            fprintf('   [WARN] parpool failed (%s); running serially.\n', ME.message);
            pool = [];
        end
    end
    if ngpu >= 1 && ~isempty(pool)
        spmd
            gpuDevice(mod(labindex - 1, ngpu) + 1);
        end
        fprintf('   Pinned %d worker(s) across %d GPU(s).\n', pool.NumWorkers, ngpu);
    elseif ngpu == 0
        fprintf('   No GPU detected; reconstructions run on CPU worker(s).\n');
    end

    null_pr   = nan(N, Kp, Kd);
    signal_pr = nan(N, Kp, Kd);
    R1_sum    = zeros(gridc);
    R3_sum    = zeros(gridc);

    ens_tic = tic;
    parfor r = 1:N
        seedA = base_seed + 3*(r-1) + 1;
        seedB = base_seed + 3*(r-1) + 2;
        seedC = base_seed + 3*(r-1) + 3;

        sdA1 = redraw_noisy_deconv(bundle1, seedA);
        sdA2 = redraw_noisy_deconv(bundle1, seedC);
        sdB  = redraw_noisy_deconv(bundle2, seedB);

        recA1 = reconstruct_recon_dose(bundle1, sdA1);
        recA2 = reconstruct_recon_dose(bundle1, sdA2);
        recB  = reconstruct_recon_dose(bundle2, sdB);

        % Both comparisons use recA1 as the reference (same low-dose mask).
        null_pr(r, :, :)   = gamma_pass_grid(recA1, recA2, pct_vec, dta_vec, spacing_mm, cutoff_f);
        signal_pr(r, :, :) = gamma_pass_grid(recA1, recB,  pct_vec, dta_vec, spacing_mm, cutoff_f);

        R1_sum = R1_sum + recA1;
        R3_sum = R3_sum + recB;

        fprintf('   [Ensemble] realization %d/%d complete.\n', r, N);
    end
    fprintf('   Ensemble compute time: %.1f s\n', toc(ens_tic));

    ENS = struct();
    ENS.null_pr   = null_pr;
    ENS.signal_pr = signal_pr;
    ENS.R1_mean   = R1_sum / N;
    ENS.R3_mean   = R3_sum / N;
    ENS.num_realizations = N;
    ENS.base_seed = base_seed;
end


function sd = redraw_noisy_deconv(B, seed)
%REDRAW_NOISY_DECONV  Fresh electronic-noise realization on the cached pre-noise
%  sensor response, then Wiener-deconvolve the pulse kernel (steps 3-4 of the
%  nominal pulse model with a seeded RNG). Copied from
%  study_gamma_pass_rates_chart.m.
    rng(seed);
    noisy = B.sensorData_resp + B.noise_amp * randn(size(B.sensorData_resp));
    deconv = real(ifft( ...
        fft(noisy, [], 2) .* B.H_conj ./ (B.H_power + B.conv_deconv_lambda), [], 2));
    sd = single(deconv);
end


function recon_dose = reconstruct_recon_dose(B, sensorData_measured)
%RECONSTRUCT_RECON_DOSE  Reconstruct a dose from a (noisy, deconvolved) sensor
%  trace using the cached bundle B. Copied verbatim from
%  study_gamma_pass_rates_chart.m (parfor/GPU-safe; no plotting).

    gridSize = B.gridSize;
    Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);
    inputArgs = B.inputArgs;
    N_iter = B.num_time_reversal_iter;

    switch lower(B.reconstruction_method)
    case 'tr'
        reconPressure      = zeros(gridSize);
        reconPressure_prev = zeros(gridSize);
        sd      = sensorData_measured;
        sd_meas = sensorData_measured;
        for it = 1:N_iter
            source_tr        = struct();
            source_tr.p_mask = B.sensor.mask;
            source_tr.p      = fliplr(sd);
            source_tr.p_mode = 'dirichlet';

            sensor_tr        = struct();
            sensor_tr.mask   = ones(Nx, Ny, Nz);
            sensor_tr.record = {'p_final'};

            pr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_tr, sensor_tr, inputArgs{:});
            if isstruct(pr) && isfield(pr, 'p_final')
                reconPressure = reshape(pr.p_final, [Nx, Ny, Nz]);
            else
                reconPressure = reshape(pr, [Nx, Ny, Nz]);
            end
            reconPressure = max(reconPressure, 0);

            if it > 1
                np = norm(reconPressure_prev(:));
                if np > 0
                    rc = norm(reconPressure(:) - reconPressure_prev(:)) / np;
                else
                    rc = Inf;
                end
                if rc < B.convergence_tol
                    break;
                end
            end
            reconPressure_prev = reconPressure;

            if it < N_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sdr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_resid, B.sensor, inputArgs{:});
                sd  = sd + (sd_meas - sdr);
            end
        end
        reconPressure = gather(reconPressure) * B.correction_factor;

    case 'das'
        reconPressure = das_reconstruct(sensorData_measured, B.sensor, B.sensor_info_orig, ...
            B.medium_pad, Nx, Ny, Nz, B.dx, B.dy, B.dz, B.dt, B.das_opts);
        reconPressure = reconPressure * B.correction_factor;

    case 'hybrid'
        reconPressure = das_reconstruct(sensorData_measured, B.sensor, B.sensor_info_orig, ...
            B.medium_pad, Nx, Ny, Nz, B.dx, B.dy, B.dz, B.dt, B.das_opts);
        reconPressure_prev = reconPressure;
        sd      = sensorData_measured;
        sd_meas = sensorData_measured;
        if N_iter > 1
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sdr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_resid, B.sensor, inputArgs{:});
            sd  = sd + (sd_meas - sdr);
        end
        for it = 2:N_iter
            source_tr        = struct();
            source_tr.p_mask = B.sensor.mask;
            source_tr.p      = fliplr(sd);
            source_tr.p_mode = 'dirichlet';

            sensor_tr        = struct();
            sensor_tr.mask   = ones(Nx, Ny, Nz);
            sensor_tr.record = {'p_final'};

            pr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_tr, sensor_tr, inputArgs{:});
            if isstruct(pr) && isfield(pr, 'p_final')
                reconPressure = reshape(pr.p_final, [Nx, Ny, Nz]);
            else
                reconPressure = reshape(pr, [Nx, Ny, Nz]);
            end
            reconPressure = max(reconPressure, 0);

            np = norm(reconPressure_prev(:));
            if np > 0
                rc = norm(reconPressure(:) - reconPressure_prev(:)) / np;
            else
                rc = Inf;
            end
            if rc < B.convergence_tol
                break;
            end
            reconPressure_prev = reconPressure;

            if it < N_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sdr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_resid, B.sensor, inputArgs{:});
                sd  = sd + (sd_meas - sdr);
            end
        end
        reconPressure = gather(reconPressure) * B.correction_factor;

    otherwise
        error('reconstruct_recon_dose:UnknownMethod', ...
            'Unknown reconstruction_method: "%s"', B.reconstruction_method);
    end

    if B.did_pad
        reconPressure = reconPressure(1:B.Nx_orig, 1:B.Ny_orig, 1:B.Nz_orig);
    end

    if B.use_pressure_scale_correction
        p0_max_orig = max(B.p0(1:B.Nx_orig, 1:B.Ny_orig, 1:B.Nz_orig), [], 'all');
        recon_max   = max(reconPressure(:));
        if recon_max > 0
            reconPressure = reconPressure * (p0_max_orig / recon_max);
        end
    end

    conversionFactor = B.medium_orig.gruneisen .* B.medium_orig.density;
    conversionFactor(conversionFactor == 0) = 1;
    reconDosePerPulse = reconPressure ./ conversionFactor;

    cropSize = B.gridSize_orig;
    if ~isempty(B.bodyMask) && isequal(size(B.bodyMask), cropSize)
        body_mask_plot = B.bodyMask;
    else
        body_mask_plot = ones(cropSize);
    end

    if ~B.mask_recon_to_dose_region
        recon_dose = reconDosePerPulse * B.num_pulses .* body_mask_plot;
    else
        recon_dose = reconDosePerPulse * B.num_pulses .* double(B.doseMask) .* body_mask_plot;
    end
end


% -------------------------------------------------------------------------
%  GAMMA + DISCRIMINABILITY HELPERS
% -------------------------------------------------------------------------

function P = gamma_pass_grid(refDose, tgtDose, pct_vec, dta_vec, spacing_mm, cutoff_frac)
%GAMMA_PASS_GRID  Pass rate (%) for every (pct, dta) pair. Reference = refDose
%  with a cutoff_frac low-dose cutoff. Serial over criteria (the realization
%  loop is the parallel one).
    ref_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(refDose));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(tgtDose));

    cutoff    = cutoff_frac * max(refDose(:));
    eval_mask = refDose >= cutoff;

    Kp = numel(pct_vec); Kd = numel(dta_vec);
    P  = nan(Kp, Kd);
    for ip = 1:Kp
        for id = 1:Kd
            pct = pct_vec(ip); dta = dta_vec(id);
            try
                gmap = CalcGamma(ref_struct, tgt_struct, pct, dta, ...
                    'local', 0, 'limit', dta * 2, 'restrict', 1);
                P(ip, id) = 100 * mean(gmap(eval_mask) <= 1);
            catch
                P(ip, id) = NaN;
            end
        end
    end
end

function pr = gamma_diag_pass(refDose, tgtDose, n_vec, spacing_mm, cutoff_frac)
%GAMMA_DIAG_PASS  Pass rate (%) along the diagonal n%/n mm for n in n_vec.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(refDose));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(tgtDose));
    cutoff    = cutoff_frac * max(refDose(:));
    eval_mask = refDose >= cutoff;

    pr = nan(numel(n_vec), 1);
    for k = 1:numel(n_vec)
        c = n_vec(k);
        try
            gmap  = CalcGamma(ref_struct, tgt_struct, c, c, ...
                'local', 0, 'limit', c * 2, 'restrict', 1);
            pr(k) = 100 * mean(gmap(eval_mask) <= 1);
        catch
            pr(k) = NaN;
        end
    end
end

function d = compute_dprime(null_samps, sig_samps)
%COMPUTE_DPRIME  d' = (mean_null - mean_signal) / pooled_std. Null passes higher
%  than signal, so a larger positive d' means better separability.
    n = null_samps(:); s = sig_samps(:);
    n = n(~isnan(n));  s = s(~isnan(s));
    if isempty(n) || isempty(s)
        d = NaN; return;
    end
    pooled = sqrt((var(n, 1) + var(s, 1)) / 2);   % population variance
    if pooled < eps
        if mean(n) == mean(s)
            d = 0;
        else
            d = sign(mean(n) - mean(s)) * Inf;
        end
    else
        d = (mean(n) - mean(s)) / pooled;
    end
end

function a = compute_auc(null_samps, sig_samps)
%COMPUTE_AUC  ROC AUC = P(null_sample > signal_sample) (Mann-Whitney), ties 0.5.
%  AUC near 1 => the null reliably passes higher than the signal (detectable).
    n = null_samps(:); s = sig_samps(:);
    n = n(~isnan(n));  s = s(~isnan(s));
    if isempty(n) || isempty(s)
        a = NaN; return;
    end
    c = 0;
    for i = 1:numel(n)
        c = c + sum(n(i) > s) + 0.5 * sum(n(i) == s);
    end
    a = c / (numel(n) * numel(s));
end

function sel = select_criteria(null_pr, signal_pr, pct_vec, dta_vec, metric)
%SELECT_CRITERIA  d'/AUC for every (pct, dta), the argmax, and the fixed-mm /
%  fixed-% slices through it. null_pr/signal_pr are N x Kp x Kd.
    Kp = numel(pct_vec); Kd = numel(dta_vec);
    D = nan(Kp, Kd); A = nan(Kp, Kd);
    for ip = 1:Kp
        for id = 1:Kd
            nn = null_pr(:, ip, id);
            ss = signal_pr(:, ip, id);
            D(ip, id) = compute_dprime(nn, ss);
            A(ip, id) = compute_auc(nn, ss);
        end
    end

    switch lower(metric)
        case 'dprime', M = D;
        otherwise,     M = A;   % default AUC
    end

    Mfinite = M; Mfinite(~isfinite(Mfinite)) = -Inf;
    [~, lin] = max(Mfinite(:));
    [bp, bd] = ind2sub(size(Mfinite), lin);

    sel = struct();
    sel.metric        = metric;
    sel.pct_vec       = pct_vec;
    sel.dta_vec       = dta_vec;
    sel.dprime_grid   = D;
    sel.auc_grid      = A;
    sel.metric_grid   = M;
    sel.best_ip       = bp;
    sel.best_id       = bd;
    sel.best_pct      = pct_vec(bp);
    sel.best_dta      = dta_vec(bd);
    sel.best_dprime   = D(bp, bd);
    sel.best_auc      = A(bp, bd);
    % Fixed-mm slice: vary % at the best DTA. Fixed-% slice: vary mm at best %.
    sel.fixed_mm_dta      = dta_vec(bd);
    sel.fixed_mm_metric   = M(:, bd);     % over pct_vec
    sel.fixed_pct_pct     = pct_vec(bp);
    sel.fixed_pct_metric  = M(bp, :).';   % over dta_vec
end


% -------------------------------------------------------------------------
%  STUDY 1: NOISE NULL  (detection baseline)
% -------------------------------------------------------------------------

function S = study_noise_null(ENS, CONFIG, plots_dir)
%STUDY_NOISE_NULL  Null = recon vs itself under independent noise; signal = the
%  CT_1-vs-CT_3 comparison. Reports the pass-rate band of each plus d'/AUC.

    n_vec   = ENS.gamma_sweep_n(:);
    pct_vec = ENS.pct_vec(:);
    dta_vec = ENS.dta_vec(:);
    K = numel(n_vec);

    null_diag   = nan(ENS.num_realizations, K);
    signal_diag = nan(ENS.num_realizations, K);
    for k = 1:K
        ip = find(abs(pct_vec - n_vec(k)) < 1e-9, 1);
        id = find(abs(dta_vec - n_vec(k)) < 1e-9, 1);
        if isempty(ip) || isempty(id)
            error('study_noise_null:NotOnGrid', ...
                'Diagonal value %g is not present in gamma_grid_pct/dta.', n_vec(k));
        end
        null_diag(:, k)   = ENS.null_pr(:, ip, id);
        signal_diag(:, k) = ENS.signal_pr(:, ip, id);
    end

    S = struct();
    S.n_vec       = n_vec;
    S.null_pr     = null_diag;
    S.signal_pr   = signal_diag;
    S.null_mean   = mean(null_diag, 1, 'omitnan');
    S.null_std    = std(null_diag, 0, 1, 'omitnan');
    S.signal_mean = mean(signal_diag, 1, 'omitnan');
    S.signal_std  = std(signal_diag, 0, 1, 'omitnan');
    S.dprime      = arrayfun(@(k) compute_dprime(null_diag(:,k), signal_diag(:,k)), 1:K);
    S.auc         = arrayfun(@(k) compute_auc(null_diag(:,k), signal_diag(:,k)), 1:K);

    fprintf('   crit   null%%(mean+/-std)   signal%%(mean+/-std)    d''     AUC\n');
    for k = 1:K
        fprintf('   %4g   %6.2f +/- %4.2f     %6.2f +/- %4.2f     %5.2f  %5.3f\n', ...
            n_vec(k), S.null_mean(k), S.null_std(k), ...
            S.signal_mean(k), S.signal_std(k), S.dprime(k), S.auc(k));
    end

    % --- Plot: null vs signal bands + d'/AUC ---
    fig = figure('Name', 'Noise Null', 'Color', 'w', 'Position', [80, 80, 1100, 460]);

    subplot(1, 2, 1);
    hold on;
    fill_band(n_vec, S.null_mean, S.null_std, [0.2 0.5 0.9]);
    fill_band(n_vec, S.signal_mean, S.signal_std, [0.85 0.33 0.10]);
    hN = plot(n_vec, S.null_mean, '-o', 'Color', [0.1 0.35 0.8], 'LineWidth', 1.8, ...
        'MarkerFaceColor', [0.1 0.35 0.8], 'MarkerSize', 5);
    hS = plot(n_vec, S.signal_mean, '-s', 'Color', [0.75 0.2 0.05], 'LineWidth', 1.8, ...
        'MarkerFaceColor', [0.75 0.2 0.05], 'MarkerSize', 5);
    yline(90, 'r--', '90%', 'LineWidth', 1.2);
    grid on; box on;
    xlabel('Gamma criterion  n%/n mm'); ylabel('Pass rate (%)');
    ylim([0 101]);
    legend([hN hS], {'Null (recon vs itself)', 'Signal (CT_1 vs CT_3)'}, ...
        'Location', 'southwest');
    title('Pass-rate bands (mean \pm std)', 'FontWeight', 'bold');

    subplot(1, 2, 2);
    yyaxis left;
    plot(n_vec, S.dprime, '-o', 'LineWidth', 1.8, 'MarkerFaceColor', [0 0.447 0.741]);
    ylabel('d'' (null vs signal)');
    yyaxis right;
    plot(n_vec, S.auc, '-s', 'LineWidth', 1.8, 'MarkerFaceColor', [0.851 0.325 0.098]);
    ylabel('ROC AUC'); ylim([0.4 1.02]);
    grid on; box on; xlabel('Gamma criterion  n%/n mm');
    title('Discriminability vs criterion', 'FontWeight', 'bold');

    sgtitle(sprintf('Study 1: Detection baseline  |  %s/%s  hash %s', ...
        ENS.patient_id, ENS.session, ENS.hash8), 'FontWeight', 'bold');
    save_fig(fig, plots_dir, 'study1_noise_null');
end


% -------------------------------------------------------------------------
%  STUDY 2: RECONSTRUCTION FIDELITY  (differential artifact)
% -------------------------------------------------------------------------

function S = study_recon_fidelity(ENS, CONFIG, plots_dir)
%STUDY_RECON_FIDELITY  Per-scenario error R_s - D_s and the differential
%  (R3-D3)-(R1-D1), the part that does NOT common-mode-cancel.

    D1 = ENS.D1; D3 = ENS.D3;
    R1 = ENS.R1_mean; R3 = ENS.R3_mean;
    mask = logical(ENS.doseMask);

    E1 = R1 - D1;            % CT_1 reconstruction error
    E3 = R3 - D3;            % CT_3 reconstruction error
    Ed = E3 - E1;            % = (R3-R1) - (D3-D1): recon change minus truth change
    truth_change = D3 - D1;
    recon_change = R3 - R1;

    S = struct();
    S.stats_E1 = error_stats(E1, mask);
    S.stats_E3 = error_stats(E3, mask);
    S.stats_Ed = error_stats(Ed, mask);
    S.stats_truth_change = error_stats(truth_change, mask);
    S.stats_recon_change = error_stats(recon_change, mask);

    % Common-mode cancellation factor: how much the differential suppresses the
    % per-scenario error RMS. >1 means cancellation helped.
    rms_perscenario = 0.5 * (S.stats_E1.rms + S.stats_E3.rms);
    if S.stats_Ed.rms > 0
        S.common_mode_suppression = rms_perscenario / S.stats_Ed.rms;
    else
        S.common_mode_suppression = Inf;
    end

    % ETHOS auxiliary (global truth common-mode-cancels in the differential).
    if ~isempty(ENS.ethos_truth)
        Ee1 = R1 - ENS.ethos_truth;
        Ee3 = R3 - ENS.ethos_truth;
        S.stats_E1_ethos = error_stats(Ee1, mask);
        S.stats_E3_ethos = error_stats(Ee3, mask);
        S.stats_Ed_ethos = error_stats(Ee3 - Ee1, mask);   % == Ed (truth cancels)
    end

    fprintf('   E1  (R1-D1)  : mean %+.4f  rms %.4f  max|.| %.4f Gy\n', ...
        S.stats_E1.mean, S.stats_E1.rms, S.stats_E1.maxabs);
    fprintf('   E3  (R3-D3)  : mean %+.4f  rms %.4f  max|.| %.4f Gy\n', ...
        S.stats_E3.mean, S.stats_E3.rms, S.stats_E3.maxabs);
    fprintf('   Ediff        : mean %+.4f  rms %.4f  max|.| %.4f Gy\n', ...
        S.stats_Ed.mean, S.stats_Ed.rms, S.stats_Ed.maxabs);
    fprintf('   common-mode suppression (perscenario RMS / diff RMS): %.2fx\n', ...
        S.common_mode_suppression);

    % --- Center on the combined max-dose voxel ---
    Dsum = (D1 + D3) .* double(mask);
    [~, mi] = max(Dsum(:));
    [cy, cx, cz] = ind2sub(size(D1), mi);   % (Y, X, Z)
    S.center = [cy, cx, cz];

    % Symmetric color limits from a robust percentile of the differential.
    edv = abs(Ed(mask));
    cl  = prctile(edv(~isnan(edv)), 99);
    if ~(cl > 0), cl = max(edv(:)) + eps; end

    rows = {E1, 'R_1 - D_1  (CT_1 error)'; ...
            E3, 'R_3 - D_3  (CT_3 error)'; ...
            Ed, '(R_3-D_3)-(R_1-D_1)  differential'};
    plot_fidelity_maps(rows, [cy, cx, cz], cl, ENS, plots_dir);

    S.color_limit = cl;
end


% -------------------------------------------------------------------------
%  STUDY 3: DETECTABILITY CURVE
% -------------------------------------------------------------------------

function S = study_detectability_curve(ENS, CONFIG, plots_dir)
%STUDY_DETECTABILITY_CURVE  Gamma pass rate vs criterion for three series:
%  truth CT_1-vs-CT_3, recon CT_1-vs-CT_3 (ensemble mean), and the null floor.

    n_vec   = ENS.gamma_sweep_n(:);
    pct_vec = ENS.pct_vec(:);
    dta_vec = ENS.dta_vec(:);
    spacing = ENS.spacing_mm;
    cutoff  = CONFIG.low_dose_cutoff_frac;
    K = numel(n_vec);

    fprintf('   truth CT_1 vs CT_3...\n');
    truth_pr = gamma_diag_pass(ENS.D1, ENS.D3, n_vec, spacing, cutoff);
    fprintf('   recon CT_1 vs CT_3 (ensemble mean)...\n');
    recon_pr = gamma_diag_pass(ENS.R1_mean, ENS.R3_mean, n_vec, spacing, cutoff);

    % Null floor: mean +/- std of the null diagonal across realizations.
    null_diag = nan(ENS.num_realizations, K);
    for k = 1:K
        ip = find(abs(pct_vec - n_vec(k)) < 1e-9, 1);
        id = find(abs(dta_vec - n_vec(k)) < 1e-9, 1);
        null_diag(:, k) = ENS.null_pr(:, ip, id);
    end
    null_mean = mean(null_diag, 1, 'omitnan');
    null_std  = std(null_diag, 0, 1, 'omitnan');

    ondisk_pr = [];
    if ~isempty(ENS.recon_ondisk_A) && ~isempty(ENS.recon_ondisk_B)
        fprintf('   recon CT_1 vs CT_3 (on-disk pipeline)...\n');
        ondisk_pr = gamma_diag_pass(ENS.recon_ondisk_A, ENS.recon_ondisk_B, ...
            n_vec, spacing, cutoff);
    end

    S = struct();
    S.n_vec      = n_vec;
    S.truth_pr   = truth_pr(:)';
    S.recon_pr   = recon_pr(:)';
    S.null_mean  = null_mean;
    S.null_std   = null_std;
    S.ondisk_pr  = ondisk_pr;

    fig = figure('Name', 'Detectability Curve', 'Color', 'w', ...
        'Position', [90, 90, 760, 520]);
    hold on;
    fill_band(n_vec, null_mean, null_std, [0.5 0.5 0.5]);
    h1 = plot(n_vec, truth_pr, '-^', 'Color', [0 0.5 0], 'LineWidth', 1.9, ...
        'MarkerFaceColor', [0 0.5 0], 'MarkerSize', 6);
    h2 = plot(n_vec, recon_pr, '-s', 'Color', [0.75 0.2 0.05], 'LineWidth', 1.9, ...
        'MarkerFaceColor', [0.75 0.2 0.05], 'MarkerSize', 6);
    h3 = plot(n_vec, null_mean, '-o', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.6, ...
        'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 5);
    hh = [h1 h2 h3];
    labels = {'Truth: CT_1 vs CT_3 (RayStation)', ...
              'Recon: CT_1 vs CT_3 (ensemble mean)', ...
              'Null floor (recon vs itself)'};
    if ~isempty(ondisk_pr)
        h4 = plot(n_vec, ondisk_pr, ':d', 'Color', [0.85 0.5 0.05], 'LineWidth', 1.4, ...
            'MarkerSize', 5);
        hh(end+1) = h4;
        labels{end+1} = 'Recon: CT_1 vs CT_3 (on-disk)';
    end
    yline(90, 'r--', '90%', 'LineWidth', 1.2);
    grid on; box on;
    xlabel('Gamma criterion  n%/n mm'); ylabel('Pass rate (%)');
    ylim([0 101]);
    legend(hh, labels, 'Location', 'southeast');
    title(sprintf('Study 3: Detectability curve  |  %s/%s  hash %s', ...
        ENS.patient_id, ENS.session, ENS.hash8), 'FontWeight', 'bold');
    save_fig(fig, plots_dir, 'study3_detectability_curve');
end


% -------------------------------------------------------------------------
%  PLOTTING HELPERS
% -------------------------------------------------------------------------

function plot_fidelity_maps(rows, center, cl, ENS, plots_dir)
%PLOT_FIDELITY_MAPS  3 rows (E1, E3, Ediff) x 3 ortho columns at CENTER.
    cy = center(1); cx = center(2); cz = center(3);
    cmap = diverging_colormap();
    colNames = {'Sagittal (X fixed)', 'Coronal (Y fixed)', 'Transverse (Z fixed)'};

    fig = figure('Name', 'Recon Fidelity', 'Color', 'w', ...
        'Position', [60, 60, 1080, 940]);
    for ri = 1:size(rows, 1)
        M   = rows{ri, 1};
        lbl = rows{ri, 2};
        sag = squeeze(M(:, cx, :))';    % (Z, Y)
        cor = squeeze(M(cy, :, :))';    % (Z, X)
        tra = squeeze(M(:, :, cz));     % (Y, X)
        panels = {sag, cor, tra};
        for ci = 1:3
            ax = subplot(3, 3, (ri-1)*3 + ci);
            imagesc(ax, panels{ci});
            axis(ax, 'image'); axis(ax, 'xy');
            colormap(ax, cmap); caxis(ax, [-cl, cl]);
            colorbar(ax);
            if ri == 1, title(ax, colNames{ci}, 'FontWeight', 'bold'); end
            if ci == 1, ylabel(ax, lbl, 'FontWeight', 'bold'); end
        end
    end
    sgtitle(fig, sprintf('Study 2: Reconstruction fidelity (Gy)  |  %s/%s  hash %s', ...
        ENS.patient_id, ENS.session, ENS.hash8), 'FontWeight', 'bold');
    save_fig(fig, plots_dir, 'study2_recon_fidelity');
end

function plot_criteria_selection(SEL, ENS, CONFIG, plots_dir)
%PLOT_CRITERIA_SELECTION  Metric heatmap over the (pct, dta) grid with the
%  argmax marked, plus the fixed-mm and fixed-% slices through it.
    pct_vec = SEL.pct_vec; dta_vec = SEL.dta_vec;

    fig = figure('Name', 'Criteria Selection', 'Color', 'w', ...
        'Position', [70, 70, 1180, 460]);

    subplot(1, 3, 1);
    imagesc(dta_vec, pct_vec, SEL.metric_grid);
    axis xy; colorbar; colormap(gca, parula);
    hold on;
    plot(SEL.best_dta, SEL.best_pct, 'kp', 'MarkerSize', 16, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.5);
    xlabel('DTA (mm)'); ylabel('Dose diff (%)');
    title(sprintf('%s grid (\\star best %g%%/%gmm)', upper(SEL.metric), ...
        SEL.best_pct, SEL.best_dta), 'FontWeight', 'bold');

    subplot(1, 3, 2);
    plot(pct_vec, SEL.fixed_mm_metric, '-o', 'LineWidth', 1.8, ...
        'Color', [0.2 0.5 0.9], 'MarkerFaceColor', [0.2 0.5 0.9]);
    hold on; xline(SEL.best_pct, 'k--');
    grid on; box on;
    xlabel('Dose diff (%)'); ylabel(upper(SEL.metric));
    title(sprintf('Fixed DTA = %g mm', SEL.fixed_mm_dta), 'FontWeight', 'bold');

    subplot(1, 3, 3);
    plot(dta_vec, SEL.fixed_pct_metric, '-s', 'LineWidth', 1.8, ...
        'Color', [0.75 0.2 0.05], 'MarkerFaceColor', [0.75 0.2 0.05]);
    hold on; xline(SEL.best_dta, 'k--');
    grid on; box on;
    xlabel('DTA (mm)'); ylabel(upper(SEL.metric));
    title(sprintf('Fixed dose diff = %g%%', SEL.fixed_pct_pct), 'FontWeight', 'bold');

    sgtitle(fig, sprintf('Criterion selection  |  best d''=%.3f  AUC=%.3f  |  hash %s', ...
        SEL.best_dprime, SEL.best_auc, ENS.hash8), 'FontWeight', 'bold');
    save_fig(fig, plots_dir, 'criteria_selection');
end

function fill_band(x, m, s, rgb)
%FILL_BAND  Shaded mean +/- std band (NaN-safe), drawn behind the mean line.
    x = x(:)'; m = m(:)'; s = s(:)';
    good = ~isnan(m) & ~isnan(s);
    x = x(good); m = m(good); s = s(good);
    if numel(x) < 2, return; end
    xb = [x, fliplr(x)];
    yb = [m + s, fliplr(m - s)];
    fill(xb, yb, rgb, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

function st = error_stats(E, mask)
%ERROR_STATS  Summary stats of E within MASK.
    v = E(mask);
    v = v(~isnan(v));
    if isempty(v)
        st = struct('mean', NaN, 'std', NaN, 'rms', NaN, 'maxabs', NaN, 'n', 0);
        return;
    end
    st = struct();
    st.mean   = mean(v);
    st.std    = std(v);
    st.rms    = sqrt(mean(v.^2));
    st.maxabs = max(abs(v));
    st.n      = numel(v);
end

function save_fig(fig, plots_dir, name)
%SAVE_FIG  Save a figure to plots_dir/<name>.png (and .fig for later edits).
    try
        exportgraphics(fig, fullfile(plots_dir, [name, '.png']), 'Resolution', 150);
    catch
        print(fig, fullfile(plots_dir, [name, '.png']), '-dpng', '-r150');
    end
    try
        savefig(fig, fullfile(plots_dir, [name, '.fig']));
    catch
    end
end

function cmap = diverging_colormap()
%DIVERGING_COLORMAP  Blue-white-red diverging map (negative=blue, positive=red).
    n = 256; half = n / 2;
    up   = linspace(0, 1, half)';
    down = linspace(1, 0, half)';
    r = [up;            ones(half, 1)];
    g = [up;            down];
    b = [ones(half, 1); down];
    cmap = [r, g, b];
end

function cmap = gamma_colormap_gyr() %#ok<DEFNU>
%GAMMA_COLORMAP_GYR Green-yellow-red colormap for gamma index (kept for reuse).
    n    = 256;
    half = round(n / 2);
    rest = n - half;
    r = [linspace(0, 1, half)'; ones(rest, 1)];
    g = [linspace(0.85, 1, half)'; linspace(1, 0, rest)'];
    b = zeros(n, 1);
    cmap = [r, g, b];
end


% -------------------------------------------------------------------------
%  ACOUSTIC MEDIUM  (copied verbatim from study_gamma_pass_rates_chart.m)
% -------------------------------------------------------------------------

function medium = create_medium(sct, config)
%CREATE_MEDIUM Build acoustic medium from SCT data and tissue model config.

    HU       = double(sct.cubeHU);
    gridSize = size(HU);
    tables   = define_tissue_tables();

    switch lower(config.gruneisen_method)
        case 'uniform'
            medium.density     = ones(gridSize) * config.uniform_density;
            medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
            medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
            medium.alpha_power = 1.1;
            medium.gruneisen   = ones(gridSize) * config.uniform_gruneisen;

        case {'threshold_1', 'threshold_2'}
            T          = tables.(config.gruneisen_method);
            nTissues   = length(T.tissue_names);
            boundaries = T.hu_boundaries;

            medium.density     = ones(gridSize) * 1000;
            medium.sound_speed = ones(gridSize) * 1540;
            medium.alpha_coeff = ones(gridSize) * 0.5;
            medium.alpha_power = 1.1;
            medium.gruneisen   = ones(gridSize) * 0.11;

            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                medium.density(mask)     = T.density(t);
                medium.sound_speed(mask) = T.sound_speed(t);
                medium.alpha_coeff(mask) = T.alpha_coeff(t);
                medium.gruneisen(mask)   = T.gruneisen(t);
            end

        otherwise
            error('Unknown gruneisen_method: %s', config.gruneisen_method);
    end

    if config.force_uniform_density
        medium.density = ones(gridSize) * config.uniform_density;
    end
    if config.force_uniform_sound_speed
        medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
    end
    if config.force_uniform_attenuation
        medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
        medium.alpha_power = 1.1;
    end
    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gridSize) * config.uniform_gruneisen;
    end

    if isfield(sct, 'bodyMask')
        outside_body = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask')
            outside_body = outside_body | logical(sct.couchMask);
        end
        medium.density(outside_body)     = config.uniform_density;
        medium.sound_speed(outside_body) = config.uniform_sound_speed;
        medium.alpha_coeff(outside_body) = config.uniform_alpha_coeff;
        medium.gruneisen(outside_body)   = config.uniform_gruneisen;
    end

    medium.grid_size = gridSize;
end

function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding.

    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air','Lung','Fat','Water','Blood','Muscle','SoftTissue','Bone','Metal'};
    tables.threshold_1.density       = [1.2,  400,  920, 1000, 1060, 1050, 1040, 1900, 7800];
    tables.threshold_1.sound_speed   = [343,  600, 1450, 1480, 1575, 1580, 1540, 3200, 5900];
    tables.threshold_1.alpha_coeff   = [0,   0.5, 0.48, 0.0022, 0.2, 0.5, 0.5,  10,   0];
    tables.threshold_1.alpha_power   = [1.0, 1.5, 1.5,  2.0,   1.3, 1.0, 1.1,  1.0,  1.0];
    tables.threshold_1.gruneisen     = [0,   0.5, 0.7,  0.11,  0.15, 0.2, 1.0,  0,    0];

    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water','Fat','SoftTissue','Bone'};
    tables.threshold_2.density       = [1000,   920, 1040,         1900];
    tables.threshold_2.sound_speed   = [1480,  1450, 1540,         3200];
    tables.threshold_2.alpha_coeff   = [0.0022, 0.48, 0.5,         10];
    tables.threshold_2.alpha_power   = [2.0,    1.5,  1.1,         1.0];
    tables.threshold_2.gruneisen     = [0.11,   0.7,  1.0,         1.0];
end
