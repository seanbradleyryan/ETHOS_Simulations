%% =========================================================================
%  STUDY_GAMMA_PASS_RATES_CHART.m
%  Two-dose photoacoustic comparison + gamma pass-rate sweep, driven by the
%  ALREADY-RECONSTRUCTED doses on disk (no per-run k-Wave reconstruction).
%
%  The listed dose file selects a beam/segment on one CT image. This script
%  loads that field's pre-computed reconstruction AND its counterpart on the
%  OTHER CT image via load_recon_dose_data (Mode='set'), then sweeps the gamma
%  criteria n%/n mm (n = 1..5 in half-integer steps) between the two
%  reconstructions, computing every pass rate in PARALLEL across the CPU cores
%  (parfor; CalcGamma is not natively parallel). A red reference line marks the
%  90% pass rate.
%
%  The "original" dose used for the display panels and the gamma low-dose
%  cutoff mask is the RayStation field dose rs_dose (a step15 output), loaded
%  alongside each reconstruction.
%
%  When the noise ensemble is enabled, the expensive part (k-Wave forward
%  simulation) is run ONCE per dose via build_forward_bundle to obtain the
%  sensor response; only the noise draw + Wiener deconvolution + reconstruction
%  are repeated across realizations (GPU-pinned) to produce the error bars.
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

CONFIG.dose_filename = 'dose_1194203_Session_1_reference_CT_3_B15_112.mat';

% The two CT image indices to compare. The listed dose file above carries a
% _CT_k token selecting one of these; this script automatically locates that
% field's counterpart on the OTHER CT image among the loaded set, and
% gamma-compares the two reconstructions.
CONFIG.ct_pair = [1, 3];

% Explicit recon config-hash override ('' => auto-discover on disk via loader).
CONFIG.config_hash = '';

CONFIG.sensor_placement_method = 'determine_sensor_mask';
CONFIG.sensor_x_index = 2;
CONFIG.sensor_y_index = 4;

% Physical 2D ultrasound array geometry (sparse element mask).
% Kerf is derived inside determine_sensor_mask as (pitch - size).
CONFIG.elements_per_side = 32;
CONFIG.element_pitch_mm  = 4.35;
CONFIG.element_size_mm   = 3.64;
CONFIG.sensor_standoff_mm = 5;
CONFIG.jaw_margin_mm      = 10;
CONFIG.sensor_placement   = 'anterior';
CONFIG.aim_at_iso         = true;
CONFIG.force_turn_angle   = 290;     % forced turn must be 290 deg (rotation-bug workaround)

CONFIG.gruneisen_method = 'threshold_2';

CONFIG.force_uniform_density     = false;
CONFIG.force_uniform_sound_speed = false;
CONFIG.force_uniform_attenuation = false;
CONFIG.force_uniform_gruneisen   = false;

CONFIG.uniform_density      = 1000;
CONFIG.uniform_sound_speed  = 1540;
CONFIG.uniform_alpha_coeff  = 0;
CONFIG.uniform_alpha_power  = 1.1;
CONFIG.uniform_gruneisen    = 1.0;

CONFIG.dose_per_pulse_cGy     = 0.16;
CONFIG.meterset               = 140;   % overridden per-dose from rtplan.meterset when available
CONFIG.pml_size               = 10;
CONFIG.cfl_number             = 0.3;
CONFIG.use_gpu                = true;
CONFIG.correction_factor           = 1.9;
CONFIG.use_pressure_scale_correction = false;   % divide max(p0) / max(recon_pressure) before dose conversion
CONFIG.mask_recon_to_dose_region     = true;    % zero recon dose outside the dose mask (>1% of original max).

% --- Reconstruction method (noise ensemble only) ---
%   'tr'     : iterative time-reversal (k-Wave back-propagation)
%   'DAS'    : Delay-And-Sum back-projection (homogeneous c, non-iterative)
%   'hybrid' : DAS for iter 1, k-Wave TR with residual correction for iters 2..N
CONFIG.reconstruction_method = 'tr';

CONFIG.num_time_reversal_iter = 10;
CONFIG.convergence_tol        = 1e-3;

% --- Pulse Convolution / Noise / Deconvolution ---
% Mimics a finite transducer impulse response applied to forward sensor data.
% Set convolution_kernel to 0 to disable the entire block (also disables the
% noise ensemble, which has nothing to vary without it).
CONFIG.convolution_kernel  = 4e-6;   % Gaussian sigma in seconds (4 us)
CONFIG.conv_noise_level    = 0.125;  % Noise amplitude as fraction of peak sensor signal
CONFIG.conv_deconv_lambda  = 1e-4;   % Wiener regularization for deconvolution

CONFIG.downscale_factor = 1;
CONFIG.use_grid_padding = true;

CONFIG.save_results = true;
CONFIG.output_file  = 'standalone_recon_results.mat';
CONFIG.plot_results = true;

% Gaussian sigma (in voxels) applied to dose slices for DISPLAY ONLY in the
% dose comparison panels. Set to 0 to disable display smoothing.
CONFIG.viz_smooth_sigma = 1.0;

% Colour scaling for the final reconstructed-vs-original dose panels:
%   'relative' - each dose (original, reconstruction) is normalized to its OWN
%                maximum, so every panel peaks at 100% and the spatial shapes
%                are compared independent of absolute magnitude.
%   'absolute' - both rows share the absolute Gy range [0, max(original)] so
%                magnitudes are directly comparable.
CONFIG.dose_panel_scale = 'relative';

% Normalize: divide both reconstructed doses by their own max before
% comparison / gamma so each peaks at 1.0.
CONFIG.normalize = false;

% Gamma logging: append CONFIG + gamma pass rates to gamma_log.mat (in
% working_dir) after each run.
CONFIG.log_gamma = false;
CONFIG.gamma_log_file = 'gamma_log.mat';

% --- Noise ensemble (gamma pass-rate error bars) ---------------------------
% Estimates how much measurement noise perturbs the gamma pass-rate curve.
% The expensive k-Wave forward simulation runs ONCE per dose; only the noise
% draw + Wiener deconvolution + reconstruction are repeated across
% num_realizations realizations, each with an independent noise seed. The
% pass-rate spread (std across realizations) becomes the error bar.
%
% The realizations are distributed across all available GPUs (one worker per
% GPU, each pinned to a distinct device). Requires the pulse model to be active
% (CONFIG.convolution_kernel > 0), since that is where electronic noise enters.
CONFIG.noise_ensemble.enable           = true;
CONFIG.noise_ensemble.num_realizations = 8;        % ensemble size N
CONFIG.noise_ensemble.recompute        = true;     % false -> load saved results if present
CONFIG.noise_ensemble.results_file     = 'gamma_noise_ensemble.mat';
CONFIG.noise_ensemble.num_iters        = [];       % [] -> use CONFIG.num_time_reversal_iter
CONFIG.noise_ensemble.base_seed        = 42;       % RNG base; per-realization seeds derive from it
%
% Null-hypothesis baseline: in addition to the signal comparison (reference CT
% recon vs the counterpart CT recon), also reconstruct the REFERENCE dose a
% second time under independent noise and gamma-compare it to itself. That
% "reference vs itself" pass-rate band is the noise floor expected when there is
% NO real dose change, so the signal band dropping below it means the CT
% difference is detectable above noise. Adds one extra reconstruction per
% realization. Meaningful only with the noise ensemble (no noise -> trivially
% 100%).
CONFIG.noise_ensemble.include_null     = true;

%% ===================== LOAD THE RECON DOSE PAIR =========================
%  Parse the beam/segment/plan-type/CT from the listed dose filename, load the
%  matched fields on BOTH CT images via load_recon_dose_data (Mode='set'), and
%  pick the listed CT (A) and its counterpart (B). Each field carries its
%  pre-computed reconstruction (recon_dose), the RayStation truth (rs_dose),
%  the CBCT geometry (cbct) and the RTPLAN stats (rtplan).

sel = parse_dose_selection(CONFIG.dose_filename);
if isempty(sel.beam)
    error('study_gamma_pass_rates:NoBeamToken', ...
        'Dose filename "%s" has no _B<beam> token; cannot select a field.', ...
        CONFIG.dose_filename);
end
if isempty(sel.ct) || strcmpi(sel.ct, 'any')
    error('study_gamma_pass_rates:NoCTtoken', ...
        ['Dose filename "%s" has no _CT_k token; cannot locate its ' ...
         'counterpart on the other CT image.'], CONFIG.dose_filename);
end

ct_a = sscanf(sel.ct, 'CT_%d');
if ~ismember(ct_a, CONFIG.ct_pair)
    error('study_gamma_pass_rates:CTnotInPair', ...
        'Dose CT index %d is not in CONFIG.ct_pair = [%s].', ct_a, num2str(CONFIG.ct_pair));
end
ct_b_all = CONFIG.ct_pair(CONFIG.ct_pair ~= ct_a);
ct_b     = ct_b_all(1);
ct_a_str = sprintf('CT_%d', ct_a);
ct_b_str = sprintf('CT_%d', ct_b);

out = load_field_set(CONFIG, sel);

iA = find_field_by_ct(out.fields, ct_a_str);
iB = find_field_by_ct(out.fields, ct_b_str);
if isempty(iA) || isempty(iB)
    found = strjoin(arrayfun(@(f) field_ct_label(f), out.fields, ...
        'UniformOutput', false), ', ');
    error('study_gamma_pass_rates:NoPair', ...
        'Beam %s / segment %s does not have both %s and %s (found CT labels: %s).', ...
        mat2str(sel.beam), mat2str(sel.segment), ct_a_str, ct_b_str, found);
end
fldA = out.fields(iA);   % listed CT
fldB = out.fields(iB);   % counterpart CT

metadata   = out.metadata;
spacing_mm = metadata.spacing(:)';
hash8      = out.config_hash;

beam_meta = [];
if isfield(metadata, 'beam_metadata') && ~isempty(metadata.beam_metadata)
    beam_meta = metadata.beam_metadata;
end

recon_A = double(fldA.recon_dose);
recon_B = double(fldB.recon_dose);
rs_A    = double(fldA.rs_dose);     % step15 truth (display "original" + cutoff)
rs_B    = double(fldB.rs_dose);

gantry_A   = getfield_def(fldA.rtplan, 'gantry_angle', 0);
gantry_B   = getfield_def(fldB.rtplan, 'gantry_angle', 0);
meterset_A = getfield_def(fldA.rtplan, 'meterset', CONFIG.meterset);
meterset_B = getfield_def(fldB.rtplan, 'meterset', CONFIG.meterset);

label_A = sprintf('%s (listed)', ct_a_str);
label_B = sprintf('%s (counterpart)', ct_b_str);

%% ========================= PRINT CONFIGURATION ===========================

fprintf('=========================================================\n');
fprintf('  Gamma Pass-Rate Chart (loaded reconstructions)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', CONFIG.dose_filename);
fprintf('  Beam/segment:    %s / %s\n', mat2str(sel.beam), mat2str(sel.segment));
fprintf('  Recon A:         %s  (%s)\n', label_A, fldA.source_mat_filename);
fprintf('  Recon B:         %s  (%s)\n', label_B, fldB.source_mat_filename);
fprintf('  Config hash:     %s\n', hash8);
fprintf('  Grid:            [%d x %d x %d]  spacing [%.2f %.2f %.2f] mm\n', ...
    size(recon_A,1), size(recon_A,2), size(recon_A,3), ...
    spacing_mm(1), spacing_mm(2), spacing_mm(3));
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('=========================================================\n\n');

if ~isequal(size(recon_A), size(recon_B))
    error('study_gamma_pass_rates:GridMismatch', ...
        ['The two loaded reconstructions are on different grids ([%s] vs [%s]).\n' ...
         'A common grid is required for the gamma comparison.'], ...
        num2str(size(recon_A)), num2str(size(recon_B)));
end

%% ===================== RESOLVE NOISE-ENSEMBLE PLAN ======================
%  Decide whether the noise ensemble is computed fresh or loaded from a saved
%  results file. The forward bundles are built only when recomputing.

ne          = CONFIG.noise_ensemble;
ne_file     = ne.results_file;
if isempty(fileparts(ne_file))
    ne_file = fullfile(CONFIG.working_dir, ne_file);
end

do_ensemble      = ne.enable;
load_ensemble    = do_ensemble && ~ne.recompute && isfile(ne_file);
compute_ensemble = do_ensemble && ~load_ensemble;

% Noise only enters the chain through the pulse model; without it the ensemble
% has nothing to vary, so disable the (re)computation in that case.
if compute_ensemble && CONFIG.convolution_kernel <= 0
    warning(['Noise ensemble requires CONFIG.convolution_kernel > 0 (electronic ' ...
        'noise is injected in the pulse model). Disabling ensemble computation.']);
    compute_ensemble = false;
    do_ensemble      = load_ensemble;   % only a prior saved file could still be shown
end

if do_ensemble
    if load_ensemble
        ne_mode_str = 'LOAD from file';
    else
        ne_mode_str = 'COMPUTE fresh';
    end
    fprintf('[NoiseEnsemble] %s  (N=%d, file=%s)\n', ...
        ne_mode_str, ne.num_realizations, ne_file);
end

%% ========================= NORMALIZE DOSES ===============================
% When enabled, divide each reconstructed dose by its own max so both peak
% at 1.0 before the comparison / gamma.

if isfield(CONFIG, 'normalize') && CONFIG.normalize
    a_max = max(recon_A(:));
    b_max = max(recon_B(:));
    fprintf('\n[NORM] Normalizing both reconstructed doses by their max:\n');
    fprintf('       %s max: %.4f Gy\n', label_A, a_max);
    fprintf('       %s max: %.4f Gy\n', label_B, b_max);
    if a_max > 0, recon_A = recon_A / a_max; end
    if b_max > 0, recon_B = recon_B / b_max; end
end

%% ========================= RESULTS SUMMARY ===============================

fprintf('\n========= COMPARISON RESULTS =========\n');
fprintf('  %-22s [%.6f, %.4f] Gy\n', [label_A ':'], min(recon_A(:)), max(recon_A(:)));
fprintf('  %-22s [%.6f, %.4f] Gy\n', [label_B ':'], min(recon_B(:)), max(recon_B(:)));

cmp_threshold = 0.01 * max(recon_A(:));
cmp_region    = recon_A > cmp_threshold;
if any(cmp_region(:))
    abs_error = abs(recon_B(cmp_region) - recon_A(cmp_region));
    rel_error = abs_error ./ max(recon_A(cmp_region), 1e-10) * 100;
    fprintf('  Mean abs diff: %.6f Gy\n', mean(abs_error));
    fprintf('  Mean rel diff: %.2f%%\n',  mean(rel_error));
    fprintf('  Max  rel diff: %.2f%%\n',  max(rel_error));
end
fprintf('=====================================\n');

%% ===================== GAMMA PASS-RATE SWEEP (PARALLEL) ==================
% Sweep the gamma criteria n%/n mm for n = 1..5 (half-integer steps) between
% the two reconstructed doses and record the pass rate at each. Each criterion
% is independent, so the CalcGamma evaluations are distributed across the CPU
% cores with parfor.
%   Reference = recon A (listed),  Target = recon B (counterpart).
%   Low-dose cutoff mask comes from the RayStation truth rs_dose (step15).

gamma_results = struct();

gamma_n        = (1:0.5:5)';                  % criterion value n
gamma_criteria = cell(numel(gamma_n), 3);       % {pct, dta, label}
for gc = 1:numel(gamma_n)
    gamma_criteria{gc, 1} = gamma_n(gc);
    gamma_criteria{gc, 2} = gamma_n(gc);
    gamma_criteria{gc, 3} = sprintf('%g%%/%gmm', gamma_n(gc), gamma_n(gc));
end

if exist('CalcGamma', 'file') == 2

    fprintf('\n[Gamma] Sweeping gamma pass rate over %d criteria (1/1 .. 5/5)...\n', numel(gamma_n));
    fprintf('        Reference: %s   |   Target: %s\n', label_A, label_B);

    % Broadcast inputs (shared, read-only across all parfor workers).
    ref_struct.start = [0, 0, 0];
    ref_struct.width = spacing_mm;
    ref_struct.data  = double(recon_A);

    tgt_struct.start = [0, 0, 0];
    tgt_struct.width = spacing_mm;
    tgt_struct.data  = double(recon_B);

    % Low-dose cutoff from the RayStation truth (step15 output). Fall back to
    % the reference recon if the truth is degenerate.
    if max(rs_A(:)) > 0
        low_dose_cutoff = 0.10 * max(rs_A(:));
        gamma_eval_mask = rs_A >= low_dose_cutoff;
    else
        low_dose_cutoff = 0.10 * max(recon_A(:));
        gamma_eval_mask = recon_A >= low_dose_cutoff;
    end

    % Ensure a parallel pool is running so the sweep uses the CPU cores.
    % If the Parallel Computing Toolbox is unavailable, parfor falls back to a
    % serial loop automatically, so failure here is non-fatal.
    if exist('parpool', 'file') == 2
        try
            if isempty(gcp('nocreate'))
                parpool('local');
            end
        catch ME
            fprintf('  [WARN] Could not start parallel pool (%s). Running serially.\n', ME.message);
        end
    else
        fprintf('  [WARN] Parallel Computing Toolbox not found; gamma sweep runs serially.\n');
    end

    pass_rates  = nan(numel(gamma_n), 1);
    crit_vals   = gamma_n;                       % local numeric copy for parfor slicing

    parfor gc = 1:numel(crit_vals)
        crit = crit_vals(gc);
        try
            gmap = CalcGamma(ref_struct, tgt_struct, crit, crit, ...
                'local', 0, 'limit', crit * 2, 'restrict', 1);
            eval_vals      = gmap(gamma_eval_mask);
            pass_rates(gc) = 100 * mean(eval_vals <= 1);
        catch ME
            warning('Gamma [%g%%/%gmm] failed: %s', crit, crit, ME.message);
            pass_rates(gc) = NaN;
        end
    end

    gamma_results.pass_rates    = pass_rates;
    gamma_results.criteria      = gamma_criteria;
    gamma_results.criterion_n   = gamma_n;
    gamma_results.cutoff_Gy     = low_dose_cutoff;
    gamma_results.eval_mask     = gamma_eval_mask;

    fprintf('\n  ------ Gamma Pass Rates: %s vs %s (10%% rs_dose cutoff) ------\n', ...
        label_A, label_B);
    for gc = 1:numel(gamma_n)
        if isnan(pass_rates(gc))
            fprintf('  %-12s  FAILED\n', gamma_criteria{gc, 3});
        else
            fprintf('  %-12s  %.2f%%\n', gamma_criteria{gc, 3}, pass_rates(gc));
        end
    end
else
    warning('CalcGamma not found. Skipping gamma analysis.');
    gamma_results = [];
end

%% ===================== NOISE ENSEMBLE (ERROR BARS) =====================
% Estimate the noise-driven spread of the gamma pass-rate curve. The forward
% simulation runs ONCE per dose (build_forward_bundle); here we re-draw noise,
% deconvolve, reconstruct, and gamma-sweep, repeated over num_realizations
% realizations distributed across all available GPUs. Results are saved so a
% later run can load them instead of recomputing.

ens = struct('available', false);

if load_ensemble
    fprintf('\n[NoiseEnsemble] Loading saved ensemble from %s ...\n', ne_file);
    S = load(ne_file);
    if isfield(S, 'pass_rates_ens') && isfield(S, 'criterion_n')
        if isequal(S.criterion_n(:), gamma_n(:))
            ens.pass_rates  = S.pass_rates_ens;
            ens.criterion_n = S.criterion_n(:);
            ens.available   = true;
            if isfield(S, 'null_pass_rates_ens') && ...
                    ~isempty(S.null_pass_rates_ens) && ...
                    ~all(isnan(S.null_pass_rates_ens(:)))
                ens.null_pass_rates = S.null_pass_rates_ens;
                ens.has_null        = true;
            else
                ens.null_pass_rates = [];
                ens.has_null        = false;
            end
            fprintf('   Loaded %d realizations x %d criteria (null: %s).\n', ...
                size(ens.pass_rates, 1), size(ens.pass_rates, 2), ...
                mat2str(ens.has_null));
        else
            warning('Saved ensemble criteria do not match current sweep; ignoring file.');
        end
    else
        warning('Saved ensemble file missing expected variables; ignoring.');
    end

elseif compute_ensemble

    % Deterministic sensor placement: use the summed plan dose so the array
    % geometry (and any grid expansion) is IDENTICAL for both CT images.
    processed_dir   = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    total_dose_file = fullfile(processed_dir, 'total_rs_dose.mat');
    if isfile(total_dose_file)
        CONFIG.total_dose_file = total_dose_file;
        fprintf('\n[NoiseEnsemble] Sensor placement dose: total_rs_dose.mat (deterministic)\n');
    else
        CONFIG.total_dose_file = '';
        fprintf('\n[NoiseEnsemble] [WARN] total_rs_dose.mat missing; placement uses per-field dose.\n');
    end

    cfgA = CONFIG; cfgA.meterset = meterset_A;
    cfgB = CONFIG; cfgB.meterset = meterset_B;
    if ~isempty(ne.num_iters)
        cfgA.num_time_reversal_iter = ne.num_iters;
        cfgB.num_time_reversal_iter = ne.num_iters;
    end

    fprintf('[NoiseEnsemble] Building forward bundles (one k-Wave forward per dose)...\n');
    bundleA = build_forward_bundle(rs_A, fldA.cbct, gantry_A, beam_meta, cfgA, label_A);
    bundleB = build_forward_bundle(rs_B, fldB.cbct, gantry_B, beam_meta, cfgB, label_B);

    if ~isequal(bundleA.gridSize_orig, bundleB.gridSize_orig)
        error('study_gamma_pass_rates:BundleGridMismatch', ...
            ['The two forward bundles ended on different grids ([%s] vs [%s]).\n' ...
             'They must share a grid for the gamma comparison.\n' ...
             'Ensure deterministic sensor placement (total_rs_dose.mat present).'], ...
            num2str(bundleA.gridSize_orig), num2str(bundleB.gridSize_orig));
    end

    N         = CONFIG.noise_ensemble.num_realizations;
    K         = numel(gamma_n);
    base_seed = CONFIG.noise_ensemble.base_seed;
    fprintf('\n[NoiseEnsemble] Computing %d realizations over %d criteria ...\n', N, K);

    % --- Multi-GPU pool setup: one worker per GPU, each pinned to a device ---
    ngpu = 0;
    try
        ngpu = gpuDeviceCount;
    catch
        ngpu = 0;
    end
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
        % Pin each worker to a distinct GPU. The selection persists on the
        % worker for the subsequent parfor, so each realization reconstructs on
        % its own device.
        spmd
            gpuDevice(mod(labindex - 1, ngpu) + 1);
        end
        fprintf('   Pinned %d worker(s) across %d GPU(s).\n', pool.NumWorkers, ngpu);
    elseif ngpu == 0
        fprintf('   No GPU detected; reconstructions run on CPU worker(s).\n');
    end

    crit_vec     = gamma_n(:);
    spacing_loc  = spacing_mm;
    include_null = CONFIG.noise_ensemble.include_null;

    pass_rates_ens      = nan(N, K);
    null_pass_rates_ens = nan(N, K);
    ens_tic = tic;
    parfor r = 1:N
        seedA = base_seed + 3*(r-1) + 1;   % reference arm (signal + null share this)
        seedB = base_seed + 3*(r-1) + 2;   % counterpart arm (signal)
        seedC = base_seed + 3*(r-1) + 3;   % reference arm, independent draw (null)

        sdA = redraw_noisy_deconv(bundleA, seedA);
        sdB = redraw_noisy_deconv(bundleB, seedB);

        recA = reconstruct_recon_dose(bundleA, sdA);
        recB = reconstruct_recon_dose(bundleB, sdB);

        % Signal: reference CT recon vs counterpart CT recon.
        pass_rates_ens(r, :) = gamma_sweep_pass_rates(recA, recB, crit_vec, spacing_loc);

        % Null: reference CT recon vs a second independent-noise recon of the
        % SAME reference dose. This is the pass-rate band expected when there is
        % no real dose change (noise floor).
        if include_null
            sdC  = redraw_noisy_deconv(bundleA, seedC);
            recC = reconstruct_recon_dose(bundleA, sdC);
            null_pass_rates_ens(r, :) = gamma_sweep_pass_rates(recA, recC, crit_vec, spacing_loc);
        end

        fprintf('   [NoiseEnsemble] realization %d/%d complete.\n', r, N);
    end
    fprintf('   Ensemble compute time: %.1f s\n', toc(ens_tic));

    ens.pass_rates      = pass_rates_ens;
    ens.null_pass_rates = null_pass_rates_ens;
    ens.has_null        = include_null;
    ens.criterion_n     = crit_vec;
    ens.available       = true;

    % --- Save the noise calculation for later reuse ---
    criterion_n   = crit_vec;   %#ok<NASGU>
    ensemble_meta = struct( ...
        'num_realizations',      N, ...
        'base_seed',             base_seed, ...
        'include_null',          include_null, ...
        'num_iters',             bundleA.num_time_reversal_iter, ...
        'reconstruction_method', CONFIG.reconstruction_method, ...
        'conv_noise_level',      CONFIG.conv_noise_level, ...
        'convolution_kernel',    CONFIG.convolution_kernel, ...
        'conv_deconv_lambda',    CONFIG.conv_deconv_lambda, ...
        'label_A',               label_A, ...
        'label_B',               label_B, ...
        'spacing_mm',            spacing_mm, ...
        'timestamp',             datestr(now, 'yyyy-mm-dd HH:MM:SS'));   %#ok<NASGU>
    save(ne_file, 'pass_rates_ens', 'null_pass_rates_ens', 'criterion_n', 'ensemble_meta', '-v7.3');
    fprintf('   Saved ensemble results to %s\n', ne_file);
end

% --- Summary statistics (mean / std across realizations) ---
if ens.available
    ens.mean = mean(ens.pass_rates, 1, 'omitnan');
    ens.std  = std(ens.pass_rates, 0, 1, 'omitnan');
    has_null = isfield(ens, 'has_null') && ens.has_null && ...
        isfield(ens, 'null_pass_rates') && ~isempty(ens.null_pass_rates);
    if has_null
        ens.null_mean = mean(ens.null_pass_rates, 1, 'omitnan');
        ens.null_std  = std(ens.null_pass_rates, 0, 1, 'omitnan');
    else
        ens.null_mean = [];
        ens.null_std  = [];
    end

    if has_null
        fprintf(['\n  ----- Noise-ensemble pass rates (mean +/- std, %d realizations) -----\n' ...
                 '  %-12s  %18s   %18s\n'], ...
            size(ens.pass_rates, 1), 'criterion', 'signal', 'null (ref vs itself)');
        for k = 1:numel(gamma_n)
            fprintf('  %-12s  %7.2f%% +/- %5.2f%%   %7.2f%% +/- %5.2f%%\n', ...
                sprintf('%g%%/%gmm', gamma_n(k), gamma_n(k)), ...
                ens.mean(k), ens.std(k), ens.null_mean(k), ens.null_std(k));
        end
    else
        fprintf('\n  ----- Noise-ensemble pass rates (mean +/- std, %d realizations) -----\n', ...
            size(ens.pass_rates, 1));
        for k = 1:numel(gamma_n)
            fprintf('  %-12s  %.2f%%  +/- %.2f%%\n', ...
                sprintf('%g%%/%gmm', gamma_n(k), gamma_n(k)), ens.mean(k), ens.std(k));
        end
    end
end

%% ===================== GAMMA PASS-RATE CHART ============================
% Pass rate as a function of the gamma criterion (1D), with a red reference
% line at the 90% pass rate. When the noise ensemble is available, error bars
% (mean +/- std over realizations) are overlaid for the signal comparison
% (reference CT vs counterpart CT). If the null hypothesis was computed, a
% second band (reference CT recon vs an independent-noise recon of itself) is
% overlaid as the no-change noise floor: the signal band dropping below the
% null band marks where the CT difference is detectable above noise.

if ~isempty(gamma_results) && isfield(gamma_results, 'pass_rates')
    if ens.available
        ens_mean_in = ens.mean;
        ens_std_in  = ens.std;
    else
        ens_mean_in = [];
        ens_std_in  = [];
    end
    if ens.available && isfield(ens, 'null_mean') && ~isempty(ens.null_mean)
        null_mean_in = ens.null_mean;
        null_std_in  = ens.null_std;
    else
        null_mean_in = [];
        null_std_in  = [];
    end
    plot_gamma_pass_rate_curve(gamma_results.criterion_n, gamma_results.pass_rates, ...
        label_A, label_B, CONFIG.patient_id, ens_mean_in, ens_std_in, ...
        null_mean_in, null_std_in);
end

%% ========================= GAMMA LOG ====================================
% Append CONFIG + gamma pass rates to a running .mat log. Maintains a
% per-criterion best record (highest pass rate and the CONFIG that produced
% it). Skipped when gamma analysis failed or produced no results.

if isfield(CONFIG, 'log_gamma') && CONFIG.log_gamma && ...
        ~isempty(gamma_results) && isfield(gamma_results, 'pass_rates')

    if isfield(CONFIG, 'gamma_log_file') && ~isempty(CONFIG.gamma_log_file)
        log_path = CONFIG.gamma_log_file;
    else
        log_path = 'gamma_log.mat';
    end
    if ~isfolder(fileparts(log_path)) && ~isempty(fileparts(log_path))
        log_path = fullfile(CONFIG.working_dir, log_path);
    elseif isempty(fileparts(log_path))
        log_path = fullfile(CONFIG.working_dir, log_path);
    end

    entry = struct();
    entry.timestamp     = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    entry.config        = CONFIG;
    entry.dose_filename = CONFIG.dose_filename;
    entry.criteria      = gamma_results.criteria(:, 3);
    entry.pass_rates    = gamma_results.pass_rates(:);

    n_crit = numel(entry.pass_rates);

    if isfile(log_path)
        L = load(log_path);
        if isfield(L, 'log_entries')
            log_entries = L.log_entries;
        else
            log_entries = struct([]);
        end
        if isfield(L, 'best')
            best = L.best;
        else
            best = repmat(struct('criterion', '', 'pass_rate', -Inf, ...
                'config', [], 'timestamp', ''), n_crit, 1);
            for gc = 1:n_crit
                best(gc).criterion = entry.criteria{gc};
            end
        end
    else
        log_entries = struct([]);
        best = repmat(struct('criterion', '', 'pass_rate', -Inf, ...
            'config', [], 'timestamp', ''), n_crit, 1);
        for gc = 1:n_crit
            best(gc).criterion = entry.criteria{gc};
        end
    end

    if isempty(log_entries)
        log_entries = entry;
    else
        log_entries(end+1) = entry;
    end

    for gc = 1:n_crit
        pr = entry.pass_rates(gc);
        if ~isnan(pr) && pr > best(gc).pass_rate
            best(gc).criterion = entry.criteria{gc};
            best(gc).pass_rate = pr;
            best(gc).config    = CONFIG;
            best(gc).timestamp = entry.timestamp;
        end
    end

    save(log_path, 'log_entries', 'best', '-v7.3');

    fprintf('\n[GammaLog] Appended run to %s (%d total entries).\n', ...
        log_path, numel(log_entries));
    fprintf('           Best pass rates so far:\n');
    for gc = 1:n_crit
        if isfinite(best(gc).pass_rate)
            fprintf('             %-12s  %.2f%%   (%s)\n', ...
                best(gc).criterion, best(gc).pass_rate, best(gc).timestamp);
        else
            fprintf('             %-12s  (none)\n', best(gc).criterion);
        end
    end
end

%% ========================= SAVE RESULTS =================================

if CONFIG.save_results
    output_fname = sprintf('standalone_comparison_%s_%s_%s_%d.mat', ...
        CONFIG.reconstruction_method, ...
        CONFIG.sensor_placement_method, CONFIG.gruneisen_method, ...
        CONFIG.num_time_reversal_iter);

    results.recon_dose_A    = recon_A;
    results.recon_dose_B    = recon_B;
    results.label_A         = label_A;
    results.label_B         = label_B;
    results.rs_dose_A       = rs_A;
    results.rs_dose_B       = rs_B;
    results.config          = CONFIG;
    results.spacing_mm      = spacing_mm;
    results.grid_size       = size(recon_A);
    results.config_hash     = hash8;
    results.source_mat_A    = fldA.source_mat_filename;
    results.source_mat_B    = fldB.source_mat_filename;
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    if ens.available
        results.noise_ensemble = ens;
    end
    save(output_fname, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', output_fname);
end

%% ========================= POST-SIMULATION VISUALIZATION ================

if CONFIG.plot_results

    % 2x3 dose comparison: RayStation truth (rs_dose, "original") top row,
    % reconstruction bottom row, for the listed CT. Coronal | Sagittal | Axial.
    density_disp = get_display_density(fldA.cbct, CONFIG);
    sensor_disp  = zeros(size(rs_A));
    plot_dose_panels(rs_A, recon_A, sensor_disp, density_disp, spacing_mm, ...
        sprintf('RayStation Truth vs Reconstruction  |  %s', label_A), ...
        {'RayStation truth', 'Reconstruction'}, CONFIG.viz_smooth_sigma, ...
        CONFIG.dose_panel_scale);

    % (gamma pass-rate chart is produced unconditionally above)

end

fprintf('\nGamma pass-rate study complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

% -------------------------------------------------------------------------
%  LOADING HELPERS
% -------------------------------------------------------------------------

function out = load_field_set(CONFIG, sel)
%LOAD_FIELD_SET Load matched fields via load_recon_dose_data (Mode='set').
%  Filters by beam/segment/plan-type parsed from the listed dose filename, but
%  NOT by CT label, so both CT images load and can be picked individually.
    args = {'Mode', 'set', 'Beam', sel.beam, 'IncludeEthos', false};
    if ~isempty(sel.segment)
        args = [args, {'Segment', sel.segment}];
    end
    if ~isempty(sel.plan_type) && ~strcmpi(sel.plan_type, 'any')
        args = [args, {'PlanType', sel.plan_type}];
    end
    if isfield(CONFIG, 'config_hash') && ~isempty(CONFIG.config_hash)
        args = [args, {'Hash', CONFIG.config_hash}];
    end
    out = load_recon_dose_data(CONFIG.patient_id, CONFIG.session, CONFIG, args{:});
end

function sel = parse_dose_selection(name)
%PARSE_DOSE_SELECTION Parse plan_type / CT / beam / segment from a dose filename.
%  Pattern: dose_<pid>_<session>_<plantype>_CT_<k>_B<beam>_<seg>.mat
    name = char(name);
    sel = struct('plan_type', 'any', 'ct', 'any', 'beam', [], 'segment', []);

    t = regexp(name, '_(reference|adapted)_', 'tokens', 'once');
    if ~isempty(t), sel.plan_type = t{1}; end

    t = regexp(name, '_CT[_-]?(\d+)', 'tokens', 'once');
    if ~isempty(t), sel.ct = sprintf('CT_%s', t{1}); end

    t = regexp(name, '_B(\d+)_', 'tokens', 'once');
    if ~isempty(t), sel.beam = str2double(t{1}); end

    t = regexp(name, '_B\d+_(\d+)\.mat$', 'tokens', 'once');
    if ~isempty(t), sel.segment = str2double(t{1}); end
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

function density = get_display_density(sct, config)
%GET_DISPLAY_DENSITY Density map (kg/m^3) from a CBCT struct for display only.
%  Built on the same grid as the loaded recon/rs doses (size(sct.cubeHU)).
    density = [];
    if ~isstruct(sct) || ~isfield(sct, 'cubeHU') || isempty(sct.cubeHU)
        return;
    end
    medium  = create_medium(sct, config);
    density = medium.density;
end

% -------------------------------------------------------------------------
%  FORWARD BUNDLE  (one k-Wave forward sim; captures everything needed to
%  re-reconstruct from a fresh noise draw). Ported from study_gamma_sensitivity.
% -------------------------------------------------------------------------

function B = build_forward_bundle(truthDose, sct, gantry_angle, beam_meta, CONFIG, label)
%BUILD_FORWARD_BUNDLE  Dose -> medium -> p0 -> sensor -> k-Wave forward -> pulse
%  model, returning the reconstruction bundle B (see reconstruct_recon_dose).

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
    % grid-expansion introduces.
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
            % zone), then re-run FFT-optimal padding.
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
%  then re-run FFT-optimal padding.

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
%  NOISE REALIZATION + RECONSTRUCTION (parfor-safe; GPU-pinned by caller)
% -------------------------------------------------------------------------

function sd = redraw_noisy_deconv(B, seed)
%REDRAW_NOISY_DECONV  Re-apply a fresh electronic-noise realization to the
%  cached pre-noise sensor response, then Wiener-deconvolve the pulse kernel.
%    B.sensorData_resp   pre-noise, post-frequency-response signal (CPU double)
%    B.noise_amp         noise amplitude (= conv_noise_level * max|resp|)
%    B.H_conj, B.H_power kernel transfer function terms
%    B.conv_deconv_lambda Wiener regularization
    rng(seed);
    noisy = B.sensorData_resp + B.noise_amp * randn(size(B.sensorData_resp));
    deconv = real(ifft( ...
        fft(noisy, [], 2) .* B.H_conj ./ (B.H_power + B.conv_deconv_lambda), [], 2));
    sd = single(deconv);
end

function recon_dose = reconstruct_recon_dose(B, sensorData_measured)
%RECONSTRUCT_RECON_DOSE  Reconstruct a dose from a (noisy, deconvolved) sensor
%  trace using the cached forward bundle B. Mirrors the nominal pipeline's
%  reconstruction -> crop -> pressure-scale -> pressure->dose stages WITHOUT any
%  plotting or convergence bookkeeping, so it is safe to run inside parfor on a
%  worker pinned to a GPU.

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
        % Iteration 1: DAS seed, then TR residual iterations 2..N.
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

    % --- Crop to original size ---
    if B.did_pad
        reconPressure = reconPressure(1:B.Nx_orig, 1:B.Ny_orig, 1:B.Nz_orig);
    end

    % --- Pressure scale correction (optional) ---
    if B.use_pressure_scale_correction
        p0_max_orig = max(B.p0(1:B.Nx_orig, 1:B.Ny_orig, 1:B.Nz_orig), [], 'all');
        recon_max   = max(reconPressure(:));
        if recon_max > 0
            reconPressure = reconPressure * (p0_max_orig / recon_max);
        end
    end

    % --- Pressure -> dose ---
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

function pr = gamma_sweep_pass_rates(reconA, reconB, crit_vec, spacing_mm)
%GAMMA_SWEEP_PASS_RATES  Pass rate for each n%/n mm criterion between two doses.
%  Serial over criteria (the realization loop is the parallel one). Reference =
%  reconA, target = reconB, with a 10% low-dose cutoff on the reference.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(reconA));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(reconB));

    cutoff    = 0.10 * max(reconA(:));
    eval_mask = reconA >= cutoff;

    pr = nan(numel(crit_vec), 1);
    for k = 1:numel(crit_vec)
        c = crit_vec(k);
        try
            gmap  = CalcGamma(ref_struct, tgt_struct, c, c, ...
                'local', 0, 'limit', c * 2, 'restrict', 1);
            pr(k) = 100 * mean(gmap(eval_mask) <= 1);
        catch
            pr(k) = NaN;
        end
    end
end

% -------------------------------------------------------------------------
%  PLOTTING
% -------------------------------------------------------------------------

function plot_gamma_pass_rate_curve(crit_n, pass_rates, label_ref, label_tgt, patient_id, ens_mean, ens_std, null_mean, null_std)
%PLOT_GAMMA_PASS_RATE_CURVE  Pass rate vs gamma criterion (n%/n mm).
%  crit_n      column vector of criterion values
%  pass_rates  nominal (single-draw) pass rates in percent (NaN entries skipped)
%  ens_mean,ens_std  optional noise-ensemble mean and std per criterion. When
%                    provided (non-empty), error bars are overlaid as the
%                    primary (signal) series and the nominal curve is shown
%                    faintly.
%  null_mean,null_std  optional null-hypothesis band (reference CT recon vs an
%                    independent-noise recon of itself). Plotted in green as the
%                    no-change noise floor: where the signal band falls below
%                    the null band, the CT difference is detectable above noise.
%  A red horizontal reference line is drawn at the 90% pass rate.

    if nargin < 6, ens_mean  = []; end
    if nargin < 7, ens_std   = []; end
    if nargin < 8, null_mean = []; end
    if nargin < 9, null_std  = []; end

    crit_n     = crit_n(:);
    pass_rates = pass_rates(:);
    have_ens   = ~isempty(ens_mean)  && ~isempty(ens_std);
    have_null  = ~isempty(null_mean) && ~isempty(null_std);

    figure('Name', 'Gamma Pass Rate vs Criterion', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [120, 120, 760, 470]);
    hold on;

    legend_handles = gobjects(0);
    legend_labels  = {};

    valid = ~isnan(pass_rates);
    if have_ens
        % Nominal single-draw curve shown faintly for reference.
        h_nom = plot(crit_n(valid), pass_rates(valid), '--o', ...
            'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.0, 'MarkerSize', 4);
        legend_handles(end+1) = h_nom; %#ok<AGROW>
        legend_labels{end+1}  = 'Nominal (single draw)';

        % Ensemble mean +/- std error bars (primary series).
        em = ens_mean(:); es = ens_std(:);
        vv = ~isnan(em);
        h_ens = errorbar(crit_n(vv), em(vv), es(vv), 'b-o', 'LineWidth', 1.8, ...
            'MarkerSize', 6, 'MarkerFaceColor', [0.2, 0.4, 1.0], 'CapSize', 8);
        legend_handles(end+1) = h_ens; %#ok<AGROW>
        if have_null
            legend_labels{end+1} = sprintf('Signal: %s vs %s (mean \\pm std)', ...
                label_ref, label_tgt);
        else
            legend_labels{end+1} = 'Ensemble mean \pm std';
        end

        % Annotate ensemble points with mean +/- std
        for k = 1:numel(crit_n)
            if vv(k)
                text(crit_n(k), em(k) + es(k) + 1.5, ...
                    sprintf('%.1f\\pm%.1f', em(k), es(k)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    else
        h_nom = plot(crit_n(valid), pass_rates(valid), 'b-o', 'LineWidth', 1.8, ...
            'MarkerSize', 6, 'MarkerFaceColor', [0.2, 0.4, 1.0]);
        legend_handles(end+1) = h_nom; %#ok<AGROW>
        legend_labels{end+1}  = 'Pass rate';

        for k = 1:numel(crit_n)
            if valid(k)
                text(crit_n(k), pass_rates(k) + 1.5, sprintf('%.1f%%', pass_rates(k)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    end

    % Null-hypothesis band: reference CT recon vs itself (independent noise).
    % This is the no-change noise floor; the signal series dropping below it
    % marks where the CT difference becomes detectable above noise.
    if have_null
        nm = null_mean(:); ns = null_std(:);
        nv = ~isnan(nm);
        null_color = [0.1, 0.6, 0.2];
        h_null = errorbar(crit_n(nv), nm(nv), ns(nv), '-s', ...
            'Color', null_color, 'LineWidth', 1.6, 'MarkerSize', 6, ...
            'MarkerFaceColor', null_color, 'CapSize', 8);
        legend_handles(end+1) = h_null; %#ok<AGROW>
        legend_labels{end+1}  = sprintf('Null: %s vs itself (mean \\pm std)', label_ref);
    end

    % 90% pass-rate reference line (red)
    yline(90, 'r-', '90% pass rate', 'LineWidth', 1.8, ...
        'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'bottom', 'FontSize', 9);

    hold off;
    grid on;
    xlim([min(crit_n) - 0.5, max(crit_n) + 0.5]);
    ylim([0, 105]);
    xticks(crit_n);
    xlabel('Gamma criterion  (n%/n mm)');
    ylabel('Gamma pass rate (%)');
    title(sprintf('Gamma Pass Rate vs Criterion   |   Patient %s\n%s vs %s', ...
        patient_id, label_ref, label_tgt), 'FontWeight', 'bold', 'FontSize', 11);
    legend(legend_handles, legend_labels, 'Location', 'southeast', 'FontSize', 8);
    drawnow;
end

function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr, rowLabels, smooth_sigma, scale_mode)
%PLOT_DOSE_PANELS 2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = original dose,  Row 2 = reconstructed dose.
%  Dose (jet, semi-transparent) is overlaid on the density map (grayscale).
%  scale_mode controls the dose colour range (and the low-dose transparency
%  threshold, which is always 10% of the relevant reference max):
%    'relative' (default) - each row is normalized to its OWN max, so every
%                panel peaks at 100% and shapes are compared independent of
%                absolute magnitude. Colorbars read in % of that row's max.
%    'absolute' - both rows share [0, max(original)] in Gy so magnitudes are
%                directly comparable. Colorbars read in Gy.
%  Isocenter at max-dose voxel.  Sensor contour in red on every panel.
%  smooth_sigma (voxels) Gaussian-smooths each slice for display only, to
%  fill speckle gaps in a spotty recon. Pass 0 to disable.
%
%  Pass density=[] to show dose only with a black background.

    if nargin < 8 || isempty(smooth_sigma), smooth_sigma = 0; end
    if nargin < 9 || isempty(scale_mode), scale_mode = 'relative'; end
    use_relative = strcmpi(scale_mode, 'relative');

    gridSize = size(original);
    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end

    have_density = ~isempty(density) && isequal(size(density), gridSize);

    [~, max_idx] = max(original(:));
    [cx, cy, cz] = ind2sub(gridSize, max_idx);

    % Absolute-mode colour scale anchored to original (reference) dose. In
    % relative mode each row instead uses its own max (computed in the loop).
    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    x_ax = (1:gridSize(1)) * spacing_mm(1);
    y_ax = (1:gridSize(2)) * spacing_mm(2);
    z_ax = (1:gridSize(3)) * spacing_mm(3);

    % Jet LUT for manual RGB blending (avoids per-axes colormap conflict)
    cmap_jet = jet(256);

    % Soft-tissue window/level for the density background (kg/m^3).
    wl_center = 1050;          % kg/m^3  (soft tissue ~1000-1080)
    wl_width  = 350;           % kg/m^3
    wl_min    = wl_center - wl_width / 2;   % 875  kg/m^3
    wl_max    = wl_center + wl_width / 2;   % 1225 kg/m^3

    if use_relative
        clim_str = 'each row scaled to its own max (relative %)';
    else
        clim_str = sprintf('Dose clim [0, %.4f] Gy', max_dose);
    end

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [50, 50, 1380, 700]);
    sgtitle(sprintf('%s\nIsocenter (max dose): X=%d  Y=%d  Z=%d voxel  |  %s', ...
        titleStr, cx, cy, cz, clim_str), 'FontWeight', 'bold', 'FontSize', 11);

    if nargin < 7 || isempty(rowLabels)
        row_labels = {'Original', 'Reconstructed'};
    else
        row_labels = rowLabels;
    end
    doses      = {original, recon};

    for row = 1:2
        d   = doses{row};
        lbl = row_labels{row};

        % Reference max for this row's colour scale and 10% transparency
        % threshold: own max in relative mode, shared original max otherwise.
        if use_relative
            row_max = max(d(:));
            if row_max == 0, row_max = 1; end
        else
            row_max = max_dose;
        end

        % Slice data for the three views
        dose_slices   = { squeeze(d(:, cy, :))',  squeeze(d(cx, :, :))',  squeeze(d(:, :, cz))' };
        sensor_slices = { squeeze(sensor_mask(:, cy, :))', ...
                          squeeze(sensor_mask(cx, :, :))', ...
                          squeeze(sensor_mask(:, :, cz))' };
        xvecs  = { x_ax, y_ax, x_ax };
        yvecs  = { z_ax, z_ax, y_ax };
        xlbls  = { 'X (mm)', 'Y (mm)', 'X (mm)' };
        ylbls  = { 'Z (mm)', 'Z (mm)', 'Y (mm)' };
        tsuffs = { sprintf('Coronal (Y=%d)', cy), ...
                   sprintf('Sagittal (X=%d)', cx), ...
                   sprintf('Axial (Z=%d)', cz) };
        if have_density
            density_slices = { squeeze(density(:, cy, :))', ...
                               squeeze(density(cx, :, :))', ...
                               squeeze(density(:, :, cz))' };
        end

        for col = 1:3
            ax         = subplot(2, 3, (row-1)*3 + col);
            dose_slice = double(dose_slices{col});
            if smooth_sigma > 0
                dose_slice = imgaussfilt(dose_slice, smooth_sigma);
            end
            xv         = xvecs{col};
            yv         = yvecs{col};

            % --- Background: density as grayscale RGB ---
            if have_density
                dn     = (density_slices{col} - wl_min) / wl_width;
                dn     = max(0, min(1, dn));   % clip to [0,1]
                bg_rgb = repmat(dn, [1, 1, 3]);
            else
                bg_rgb = zeros([size(dose_slice), 3]);   % black
            end
            image(ax, xv, yv, bg_rgb);
            hold(ax, 'on');

            % --- Foreground: dose as jet RGB with alpha mask ---
            %   mask: dose >= 10% of row_max
            norm_d   = max(0, min(1, dose_slice / row_max));
            idx      = max(1, min(256, round(norm_d * 255) + 1));
            sz       = size(dose_slice);
            dose_rgb = reshape(cmap_jet(idx(:), :), [sz, 3]);

            above      = dose_slice >= 0.10 * row_max;
            ramp       = 0.45 + 0.35 * min(1, ...
                (dose_slice - 0.10*row_max) / max(0.40*row_max, 1e-12));
            dose_alpha = above .* ramp;

            h_dose           = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = dose_alpha;

            % --- Sensor contour ---
            s = sensor_slices{col};
            if any(s(:))
                contour(ax, xv, yv, s, [0.5, 0.5], 'r-', 'LineWidth', 2);
            end
            hold(ax, 'off');

            axis(ax, 'xy'); axis(ax, 'image');

            % Colorbar: attach jet LUT. In relative mode the manual RGB blend
            % maps [0, row_max] onto [0, 100]% (norm_d), so the colorbar reads
            % in percent; in absolute mode it reads in Gy over [0, row_max].
            colormap(ax, cmap_jet);
            if use_relative
                caxis(ax, [0, 100]);
                cb = colorbar(ax); cb.Label.String = 'Relative dose (% of max)';
            else
                caxis(ax, [0, row_max]);
                cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
            end

            xlabel(ax, xlbls{col}); ylabel(ax, ylbls{col});
            title(ax, sprintf('%s  %s', lbl, tsuffs{col}));
        end
    end
    drawnow;
end

% -------------------------------------------------------------------------
%  ACOUSTIC MEDIUM
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

            fprintf('       Tissue model: %s (%d tissues)\n', config.gruneisen_method, nTissues);
            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                fprintf('         %-12s: %7d voxels (%.1f%%)\n', ...
                    T.tissue_names{t}, sum(mask(:)), 100 * sum(mask(:)) / numel(HU));
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
