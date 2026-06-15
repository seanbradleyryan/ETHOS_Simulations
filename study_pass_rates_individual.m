%% =========================================================================
%  STUDY_PASS_RATES_INDIVIDUAL.m
%  Batch, per-dose photoacoustic gamma analysis driven by the ALREADY-
%  RECONSTRUCTED doses on disk (no per-run k-Wave reconstruction, except for the
%  optional A2 noise ensemble).
%
%  For an ARBITRARY list of dose files (CONFIG.dose_filenames), each file selects
%  a beam/segment; this script loads that field's pre-computed reconstruction on
%  BOTH CT images (CT_1 and CT_3) via load_recon_dose_data (Mode='set'), together
%  with the RayStation truth (rs_dose), CBCT geometry and RTPLAN stats. It then
%  runs, per dose, two independently toggleable analyses:
%
%    A1  recon vs its own RayStation truth (CT_1)   - per-dose detector accuracy
%    A2  CT_1 recon vs CT_3 recon                    - adapted-vs-reference change
%
%  Each analysis produces (all individually gated by CONFIG.enable.*):
%    - a gamma pass-rate vs n%/n mm chart (global gamma, 10% low-dose cutoff),
%    - a 3-panel axial figure (volA | volB | signed difference) at the axial
%      (constant-Z) slice containing the truth max, and
%    - a 10%-dose-area + ultrasound-sensor visualization.
%
%  A2 optionally overlays a NOISE-ENSEMBLE NULL hypothesis: for the CBCT-1 dose,
%  the k-Wave forward simulation runs ONCE (build_forward_bundle); then noise is
%  redrawn + Wiener-deconvolved + reconstructed num_realizations times, and each
%  noisy recalculation is gamma-compared to the ORIGINAL on-disk CT_1 recon. That
%  "same anatomy, noise only" band is the noise floor: where the A2 signal curve
%  sits above the null, the CT_1-vs-CT_3 change is real rather than noise.
%
%  PARALLELIZATION
%    - Gamma sweeps: ONE flattened parfor over every (dose x analysis x
%      criterion) evaluation, so the CPU pool stays maximally busy.
%    - Noise ensemble: parfor over realizations, one GPU per worker (spmd pin).
%      The forward sim runs once per dose before the loop.
%
%  NOTE: the noise ensemble is inherently per-field (one gantry/sensor geometry),
%  so it applies to single-field dose specs, not summed totals.
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

% Arbitrary list of dose files to analyze. Each carries a _B<beam>_<seg> token
% (and a _CT_k token, which is ignored here since both CTs are loaded). Add as
% many as desired; every enabled analysis is run for each.
CONFIG.dose_filenames = { ...
    'dose_1194203_Session_1_reference_CT_3_B15_112.mat' ...
};

% The two CT image indices. A1 analyzes the LOWER index (CT_1); A2 compares the
% lower (reference) against the higher (CT_3) reconstruction.
CONFIG.ct_pair = [1, 3];

% Explicit recon config-hash override ('' => auto-discover on disk via loader).
CONFIG.config_hash = '505ae853';

% --- Per-analysis enable flags (requirement: each individually toggleable) ---
CONFIG.enable.A1             = true;   % recon vs RayStation truth (CT_1)
CONFIG.enable.A2             = true;   % CT_1 recon vs CT_3 recon
CONFIG.enable.dose_panels    = true;   % volA | volB | difference axial figure
CONFIG.enable.sensor_view    = true;   % 10% dose area + sensor model figure
CONFIG.enable.noise_ensemble = false;  % A2 null hypothesis (k-Wave; expensive)

% Gamma sweep criteria n (each evaluated as n%/n mm).
CONFIG.gamma_n = (1:0.5:5)';

% --- Sensor geometry (placement + display) ---
CONFIG.sensor_placement_method = 'determine_sensor_mask';
CONFIG.sensor_x_index = 2;
CONFIG.sensor_y_index = 4;
CONFIG.elements_per_side  = 32;
CONFIG.element_pitch_mm   = 4.35;
CONFIG.element_size_mm    = 3.64;
CONFIG.sensor_standoff_mm = 5;
CONFIG.jaw_margin_mm      = 10;
CONFIG.sensor_placement   = 'anterior';
CONFIG.aim_at_iso         = true;
CONFIG.force_turn_angle   = 290;     % forced turn must be 290 deg (rotation-bug workaround)

% --- Tissue / acoustic medium model ---
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

% --- Forward / reconstruction parameters (noise ensemble only) ---
CONFIG.dose_per_pulse_cGy     = 0.16;
CONFIG.meterset               = 140;   % overridden per-dose from rtplan.meterset when available
CONFIG.pml_size               = 10;
CONFIG.cfl_number             = 0.3;
CONFIG.use_gpu                = true;
CONFIG.correction_factor             = 1.9;
CONFIG.use_pressure_scale_correction = false;
CONFIG.mask_recon_to_dose_region     = true;

CONFIG.reconstruction_method = 'tr';   % 'tr' | 'DAS' | 'hybrid'
CONFIG.num_time_reversal_iter = 5;
CONFIG.convergence_tol        = 1e-3;

% --- Pulse convolution / noise / deconvolution (drives the noise ensemble) ---
CONFIG.convolution_kernel  = 4e-6;   % Gaussian sigma in seconds (4 us); 0 disables noise
CONFIG.conv_noise_level    = 0.125;  % noise amplitude as fraction of peak sensor signal
CONFIG.conv_deconv_lambda  = 1e-4;   % Wiener regularization for deconvolution

CONFIG.downscale_factor = 1;
CONFIG.use_grid_padding = true;

% --- Noise ensemble (A2 null hypothesis) ---
CONFIG.noise_ensemble.num_realizations = 10;   % ensemble size N (configurable)
CONFIG.noise_ensemble.base_seed        = 42;   % RNG base; per-realization seeds derive from it
CONFIG.noise_ensemble.num_iters        = [];   % [] -> use CONFIG.num_time_reversal_iter

% --- Display options ---
CONFIG.viz_smooth_sigma   = 1.0;        % Gaussian sigma (voxels) for display-only smoothing
CONFIG.dose_panel_scale   = 'relative'; % 'relative' | 'absolute'
CONFIG.dose_panel_clip_pct = 99.5;      % percentile anchoring the colour-scale top

% --- Output ---
CONFIG.save_results = true;
CONFIG.output_file  = 'pass_rates_individual_results.mat';

%% ===================== RESOLVE THE DOSE LIST ============================
%  For each dose filename: parse the beam/segment, load the matched fields on
%  BOTH CT images, and pick CT_1 (A1 + A2 reference) and CT_3 (A2 target). Cache
%  per-dose recon/truth volumes, geometry, the 10% eval mask, the display
%  density, and the (possibly grid-expanded) sensor mask for the visualizations.

ct_lo   = min(CONFIG.ct_pair);
ct_hi   = max(CONFIG.ct_pair);
ct1_str = sprintf('CT_%d', ct_lo);
ct3_str = sprintf('CT_%d', ct_hi);

nDose = numel(CONFIG.dose_filenames);
if nDose == 0
    error('study_pass_rates_individual:NoDoses', 'CONFIG.dose_filenames is empty.');
end

doses = struct([]);

for i = 1:nDose
    fn  = CONFIG.dose_filenames{i};
    sel = parse_dose_selection(fn);
    if isempty(sel.beam)
        error('study_pass_rates_individual:NoBeamToken', ...
            'Dose filename "%s" has no _B<beam> token; cannot select a field.', fn);
    end

    fprintf('\n[Load %d/%d] %s\n', i, nDose, fn);
    out = load_field_set(CONFIG, sel);

    iA = find_field_by_ct(out.fields, ct1_str);
    iB = find_field_by_ct(out.fields, ct3_str);
    if isempty(iA) || isempty(iB)
        found = strjoin(arrayfun(@(f) field_ct_label(f), out.fields, ...
            'UniformOutput', false), ', ');
        error('study_pass_rates_individual:NoPair', ...
            'Beam %s / segment %s lacks both %s and %s (found: %s).', ...
            mat2str(sel.beam), mat2str(sel.segment), ct1_str, ct3_str, found);
    end
    fA = out.fields(iA);   % CT_1 (A1 + A2 reference)
    fB = out.fields(iB);   % CT_3 (A2 target)

    bm = [];
    if isfield(out.metadata, 'beam_metadata') && ~isempty(out.metadata.beam_metadata)
        bm = out.metadata.beam_metadata;
    end

    D = struct();
    D.label       = make_dose_label(sel);
    D.filename    = fn;
    D.spacing_mm  = out.metadata.spacing(:)';
    D.beam_meta   = bm;
    D.recon_CT1   = double(fA.recon_dose);
    D.recon_CT3   = double(fB.recon_dose);
    D.rs_CT1      = double(fA.rs_dose);
    D.rs_CT3      = double(fB.rs_dose);
    D.cbct_CT1    = fA.cbct;
    D.gantry_CT1  = getfield_def(fA.rtplan, 'gantry_angle', 0);
    D.meterset_CT1 = getfield_def(fA.rtplan, 'meterset', CONFIG.meterset);
    D.source_CT1  = fA.source_mat_filename;
    D.source_CT3  = fB.source_mat_filename;

    if ~isequal(size(D.recon_CT1), size(D.recon_CT3))
        error('study_pass_rates_individual:GridMismatch', ...
            'Dose "%s": CT_1 and CT_3 recons are on different grids ([%s] vs [%s]).', ...
            D.label, num2str(size(D.recon_CT1)), num2str(size(D.recon_CT3)));
    end

    % 10% low-dose cutoff mask from the RayStation truth (consistent across A1/A2).
    if max(D.rs_CT1(:)) > 0
        D.cutoff    = 0.10 * max(D.rs_CT1(:));
        D.eval_mask = D.rs_CT1 >= D.cutoff;
    else
        D.cutoff    = 0.10 * max(D.recon_CT1(:));
        D.eval_mask = D.recon_CT1 >= D.cutoff;
    end

    D.density_CT1 = get_display_density(D.cbct_CT1, CONFIG);

    % Sensor placement for display (may grid-expand; capture the embed offset).
    D.disp_dims   = size(D.rs_CT1);
    D.embed_off   = [0, 0, 0];
    D.sensor_disp = zeros(size(D.rs_CT1));
    if CONFIG.enable.sensor_view || CONFIG.enable.dose_panels
        cfg_sensor = CONFIG;
        tdf = fullfile(CONFIG.working_dir, 'RayStationFiles', CONFIG.patient_id, ...
            CONFIG.session, 'processed', 'total_rs_dose.mat');
        if isfile(tdf), cfg_sensor.total_dose_file = tdf; end
        try
            [sm, gp]      = compute_display_sensor(D.cbct_CT1, D.rs_CT1, ...
                D.gantry_CT1, bm, cfg_sensor);
            D.sensor_disp = sm;
            D.disp_dims   = size(sm);
            D.embed_off   = [gp.y_pre, gp.x_pre, gp.z_pre];
            fprintf('  [Sensor] %d voxels on grid [%s] (embed offset [%s]).\n', ...
                nnz(sm), num2str(D.disp_dims), num2str(D.embed_off));
        catch ME
            warning('Sensor visualization unavailable for "%s" (%s); using empty mask.', ...
                D.label, ME.message);
        end
    end

    D.cA1_idx   = [];
    D.cA2_idx   = [];
    D.null_mean = [];
    D.null_std  = [];

    if i == 1
        doses = D;
    else
        doses(i) = D; %#ok<SAGROW>
    end
end

gamma_n   = CONFIG.gamma_n(:);
K         = numel(gamma_n);
has_gamma = (exist('CalcGamma', 'file') == 2);
if ~has_gamma
    warning('CalcGamma not found; gamma sweeps and ensemble will be skipped.');
end

%% ===================== GAMMA SWEEPS (FLATTENED PARFOR) ==================
%  Build the list of comparisons (A1 and/or A2 per dose), then run ONE parfor
%  over every (comparison x criterion) job so the CPU pool stays busy. Each job
%  is an independent CalcGamma evaluation with the dose's 10% truth cutoff mask.

cref = {}; ctgt = {}; cmask = {}; cspacing = {}; ckind = {}; cdose = [];
nc = 0;
for i = 1:nDose
    if CONFIG.enable.A1
        nc = nc + 1;
        cref{nc}     = single(doses(i).recon_CT1);  %#ok<SAGROW>
        ctgt{nc}     = single(doses(i).rs_CT1);     %#ok<SAGROW>
        cmask{nc}    = doses(i).eval_mask;          %#ok<SAGROW>
        cspacing{nc} = doses(i).spacing_mm;         %#ok<SAGROW>
        ckind{nc}    = 'A1';                         %#ok<SAGROW>
        cdose(nc)    = i;                            %#ok<SAGROW>
        doses(i).cA1_idx = nc;
    end
    if CONFIG.enable.A2
        nc = nc + 1;
        cref{nc}     = single(doses(i).recon_CT1);  %#ok<SAGROW>
        ctgt{nc}     = single(doses(i).recon_CT3);  %#ok<SAGROW>
        cmask{nc}    = doses(i).eval_mask;          %#ok<SAGROW>
        cspacing{nc} = doses(i).spacing_mm;         %#ok<SAGROW>
        ckind{nc}    = 'A2';                         %#ok<SAGROW>
        cdose(nc)    = i;                            %#ok<SAGROW>
        doses(i).cA2_idx = nc;
    end
end
nComp     = nc;
comp_pass = nan(nComp, K);

if has_gamma && nComp > 0
    nJobs = nComp * K;
    job_c = zeros(nJobs, 1);
    job_k = zeros(nJobs, 1);
    jj = 0;
    for c = 1:nComp
        for k = 1:K
            jj = jj + 1;
            job_c(jj) = c;
            job_k(jj) = k;
        end
    end

    ensure_pool();

    fprintf('\n[Gamma] %d comparison(s) x %d criteria = %d parallel jobs...\n', ...
        nComp, K, nJobs);

    job_pass = nan(nJobs, 1);
    critv    = gamma_n;
    parfor j = 1:nJobs
        c    = job_c(j);
        k    = job_k(j);
        crit = critv(k);
        sp   = cspacing{c};
        ref_struct = struct('start', [0, 0, 0], 'width', sp, 'data', double(cref{c}));
        tgt_struct = struct('start', [0, 0, 0], 'width', sp, 'data', double(ctgt{c}));
        try
            gmap = CalcGamma(ref_struct, tgt_struct, crit, crit, ...
                'local', 0, 'limit', crit * 2, 'restrict', 1);
            m = cmask{c};
            job_pass(j) = 100 * mean(gmap(m) <= 1);
        catch ME
            warning('Gamma job %d (%s) failed: %s', j, ckind{c}, ME.message);
            job_pass(j) = NaN;
        end
    end

    for j = 1:nJobs
        comp_pass(job_c(j), job_k(j)) = job_pass(j);
    end

    % Console summary per comparison.
    for c = 1:nComp
        fprintf('\n  ----- %s | %s (10%% truth cutoff) -----\n', ...
            ckind{c}, doses(cdose(c)).label);
        for k = 1:K
            if isnan(comp_pass(c, k))
                fprintf('    %g%%/%gmm  FAILED\n', gamma_n(k), gamma_n(k));
            else
                fprintf('    %g%%/%gmm  %.2f%%\n', gamma_n(k), gamma_n(k), comp_pass(c, k));
            end
        end
    end
end

%% ===================== A2 NOISE-ENSEMBLE NULL (OPTIONAL) ================
%  For each dose's CT_1: run the k-Wave forward simulation once, then repeatedly
%  redraw noise + deconvolve + reconstruct and gamma-compare each noisy recalc to
%  the ORIGINAL on-disk CT_1 recon. The mean +/- std band is the no-change noise
%  floor overlaid on the A2 chart. Realizations are distributed across GPUs.

do_ensemble = CONFIG.enable.A2 && CONFIG.enable.noise_ensemble && has_gamma;
if do_ensemble && CONFIG.convolution_kernel <= 0
    warning(['Noise ensemble requires CONFIG.convolution_kernel > 0 (noise enters ' ...
        'in the pulse model). Disabling ensemble.']);
    do_ensemble = false;
end

if do_ensemble
    % Deterministic sensor placement (identical geometry per CT) when available.
    tdf = fullfile(CONFIG.working_dir, 'RayStationFiles', CONFIG.patient_id, ...
        CONFIG.session, 'processed', 'total_rs_dose.mat');

    ngpu = 0;
    try, ngpu = gpuDeviceCount; catch, ngpu = 0; end
    ensure_pool(max(1, ngpu));
    pool = gcp('nocreate');
    if ngpu >= 1 && ~isempty(pool)
        spmd
            gpuDevice(mod(labindex - 1, ngpu) + 1);
        end
        fprintf('\n[NoiseEnsemble] Pinned %d worker(s) across %d GPU(s).\n', ...
            pool.NumWorkers, ngpu);
    elseif ngpu == 0
        fprintf('\n[NoiseEnsemble] No GPU detected; reconstructions run on CPU worker(s).\n');
    end

    N        = CONFIG.noise_ensemble.num_realizations;
    base     = CONFIG.noise_ensemble.base_seed;
    critv    = gamma_n;

    for i = 1:nDose
        D   = doses(i);
        cfg = CONFIG;
        cfg.meterset = D.meterset_CT1;
        if isfile(tdf), cfg.total_dose_file = tdf; end
        if ~isempty(CONFIG.noise_ensemble.num_iters)
            cfg.num_time_reversal_iter = CONFIG.noise_ensemble.num_iters;
        end

        fprintf('\n[NoiseEnsemble] Dose %d/%d (%s): forward sim once, %d realizations...\n', ...
            i, nDose, D.label, N);
        bundle = build_forward_bundle(D.rs_CT1, D.cbct_CT1, D.gantry_CT1, ...
            D.beam_meta, cfg, D.label);

        % Original on-disk recon embedded onto the bundle grid (no-op when the
        % sensor placement did not expand the grid).
        ref_recon = embed_on_grid(double(D.recon_CT1), bundle.gridSize_orig, ...
            bundle.embed_offset, 0);

        spacing_loc = bundle.spacing_mm;
        pr_ens      = nan(N, K);
        parfor r = 1:N
            seed = base + (i - 1) * 1000 + r;
            sd   = redraw_noisy_deconv(bundle, seed);
            rec  = reconstruct_recon_dose(bundle, sd);
            pr_ens(r, :) = gamma_sweep_pass_rates(ref_recon, rec, critv, spacing_loc);
            fprintf('   [NoiseEnsemble] dose %d realization %d/%d complete.\n', i, r, N);
        end

        doses(i).null_pass_rates = pr_ens;
        doses(i).null_mean       = mean(pr_ens, 1, 'omitnan');
        doses(i).null_std        = std(pr_ens, 0, 1, 'omitnan');
    end
end

%% ===================== PER-DOSE OUTPUTS =================================
%  For each dose, emit the enabled A1 / A2 figures. Doses, truths and densities
%  are embedded onto the sensor display grid so the red sensor contour
%  co-registers with the dose (no-op when the grid was not expanded).

RESULTS = struct();
RESULTS.config         = CONFIG;
RESULTS.gamma_n        = gamma_n;
RESULTS.dose_filenames = CONFIG.dose_filenames;

for i = 1:nDose
    D         = doses(i);
    disp_dims = D.disp_dims;
    embed_off = D.embed_off;

    rs1_d  = embed_on_grid(D.rs_CT1,      disp_dims, embed_off, 0);
    rec1_d = embed_on_grid(D.recon_CT1,   disp_dims, embed_off, 0);
    rec3_d = embed_on_grid(D.recon_CT3,   disp_dims, embed_off, 0);
    dens_d = embed_on_grid(D.density_CT1, disp_dims, embed_off, CONFIG.uniform_density);
    sens_d = D.sensor_disp;

    RESULTS.doses(i).label      = D.label;
    RESULTS.doses(i).filename   = D.filename;
    RESULTS.doses(i).spacing_mm = D.spacing_mm;
    RESULTS.doses(i).cutoff_Gy  = D.cutoff;

    % ----- A1: recon CT_1 vs RayStation truth -----
    if CONFIG.enable.A1
        a1 = comp_pass(D.cA1_idx, :);
        RESULTS.doses(i).A1.pass_rates = a1;

        if has_gamma
            plot_gamma_pass_rate_curve(gamma_n, a1(:), 'Recon CT\_1', ...
                'RayStation truth', sprintf('%s | %s [A1]', CONFIG.patient_id, D.label), ...
                [], [], [], []);
        end
        if CONFIG.enable.dose_panels
            plot_truth_recon_diff_axial(rs1_d, rec1_d, rs1_d, dens_d, sens_d, ...
                D.spacing_mm, {'RayStation truth', 'Recon CT\_1'}, ...
                sprintf('A1  Truth vs Recon  |  %s', D.label), ...
                CONFIG.viz_smooth_sigma, CONFIG.dose_panel_scale, CONFIG.dose_panel_clip_pct);
        end
        if CONFIG.enable.sensor_view
            dose_mask = rs1_d >= 0.10 * max(rs1_d(:));
            plot_sensor_dose_planes(dose_mask, sens_d, D.spacing_mm, dens_d, CONFIG);
        end
    end

    % ----- A2: recon CT_1 vs recon CT_3 (+ optional null) -----
    if CONFIG.enable.A2
        a2 = comp_pass(D.cA2_idx, :);
        RESULTS.doses(i).A2.pass_rates = a2;
        RESULTS.doses(i).A2.null_mean  = D.null_mean;
        RESULTS.doses(i).A2.null_std   = D.null_std;

        null_mean = D.null_mean;
        null_std  = D.null_std;

        if has_gamma
            plot_gamma_pass_rate_curve(gamma_n, a2(:), 'Recon CT\_1', ...
                'Recon CT\_3', sprintf('%s | %s [A2]', CONFIG.patient_id, D.label), ...
                [], [], null_mean, null_std);
        end
        if CONFIG.enable.dose_panels
            plot_truth_recon_diff_axial(rec1_d, rec3_d, rs1_d, dens_d, sens_d, ...
                D.spacing_mm, {'Recon CT\_1', 'Recon CT\_3'}, ...
                sprintf('A2  CT_1 vs CT_3 Recon  |  %s', D.label), ...
                CONFIG.viz_smooth_sigma, CONFIG.dose_panel_scale, CONFIG.dose_panel_clip_pct);
        end
        if CONFIG.enable.sensor_view
            dose_mask = rs1_d >= 0.10 * max(rs1_d(:));
            plot_sensor_dose_planes(dose_mask, sens_d, D.spacing_mm, dens_d, CONFIG);
        end
    end
end

%% ========================= SAVE RESULTS =================================

if CONFIG.save_results
    out_path = CONFIG.output_file;
    if isempty(fileparts(out_path))
        out_path = fullfile(CONFIG.working_dir, out_path);
    end
    save(out_path, '-struct', 'RESULTS', '-v7.3');
    fprintf('\nResults saved to: %s\n', out_path);
end

fprintf('\nPer-dose gamma pass-rate study complete (%d dose(s)).\n', nDose);


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

% -------------------------------------------------------------------------
%  POOL / LABEL HELPERS
% -------------------------------------------------------------------------

function ensure_pool(desired_workers)
%ENSURE_POOL Start a parallel pool if none exists (non-fatal if PCT absent).
    if exist('parpool', 'file') ~= 2
        fprintf('  [WARN] Parallel Computing Toolbox not found; loops run serially.\n');
        return;
    end
    try
        if isempty(gcp('nocreate'))
            if nargin >= 1 && ~isempty(desired_workers)
                parpool('local', desired_workers);
            else
                parpool('local');
            end
        end
    catch ME
        fprintf('  [WARN] Could not start parallel pool (%s). Running serially.\n', ME.message);
    end
end

function lbl = make_dose_label(sel)
%MAKE_DOSE_LABEL Compact label from a parsed dose selection (CT-independent).
    parts = {};
    if ~isempty(sel.beam),    parts{end+1} = sprintf('B%d', sel.beam);    end %#ok<AGROW>
    if ~isempty(sel.segment), parts{end+1} = sprintf('seg%d', sel.segment); end %#ok<AGROW>
    if ~isempty(sel.plan_type) && ~strcmpi(sel.plan_type, 'any')
        parts{end+1} = sel.plan_type; %#ok<AGROW>
    end
    if isempty(parts), lbl = 'dose'; else, lbl = strjoin(parts, ' '); end
end

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

function [sensor_mask, gp] = compute_display_sensor(sct, doseGrid, gantry_angle, beam_meta, config)
%COMPUTE_DISPLAY_SENSOR  Ultrasound array mask for display via determine_sensor_mask.
%  Mirrors build_forward_bundle's placement call. Returns the logical mask on
%  the (possibly grid-expanded) display grid and gp = sensor_info.grid_pad, whose
%  .y_pre/.x_pre/.z_pre give the offset to embed an original-grid dose so it
%  co-registers with the mask.
    doseGrid   = double(doseGrid);
    spacing_mm = sct.spacing(:)';
    gridSize   = size(doseGrid);

    if isfield(sct, 'bodyMask') && ~isempty(sct.bodyMask)
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    sct_for_sensor = sct;
    if ~isfield(sct_for_sensor, 'couchMask') || isempty(sct_for_sensor.couchMask)
        sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
    end
    if ~isfield(sct_for_sensor, 'origin') || isempty(sct_for_sensor.origin)
        sct_for_sensor.origin = [0, 0, 0];
    end
    sct_for_sensor.spacing    = spacing_mm;
    sct_for_sensor.dimensions = gridSize;

    field_dose              = struct();
    field_dose.dose_Gy      = doseGrid;
    field_dose.gantry_angle = gantry_angle;
    field_dose.origin       = sct_for_sensor.origin;
    field_dose.spacing      = spacing_mm;
    field_dose.dimensions   = gridSize;

    [sensor_mask, sensor_info] = determine_sensor_mask( ...
        sct_for_sensor, field_dose, beam_meta, config);
    sensor_mask = logical(sensor_mask);
    gp = sensor_info.grid_pad;
end

function v = embed_on_grid(vol, dims, off, fillval)
%EMBED_ON_GRID  Place an original-grid volume into a (larger) display grid.
%  dims = [d1 d2 d3] display grid; off = [o1 o2 o3] pre-pad per dim. Region
%  outside the embedded volume is set to fillval (default 0). Empty/already-
%  matching inputs are passed through unchanged.
    if nargin < 4 || isempty(fillval), fillval = 0; end
    if isempty(vol) || (isequal(size(vol), dims) && all(off == 0))
        v = vol;
        return;
    end
    v = fillval * ones(dims, 'like', double(vol));
    s = size(vol);
    v(off(1) + (1:s(1)), off(2) + (1:s(2)), off(3) + (1:s(3))) = vol;
end

% -------------------------------------------------------------------------
%  FORWARD BUNDLE  (one k-Wave forward sim; captures everything needed to
%  re-reconstruct from a fresh noise draw).
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

            field_dose_for_sensor              = struct();
            field_dose_for_sensor.dose_Gy      = doseGrid;
            field_dose_for_sensor.gantry_angle = gantry_angle;
            field_dose_for_sensor.origin       = sct_for_sensor.origin;
            field_dose_for_sensor.spacing      = spacing_mm;
            field_dose_for_sensor.dimensions   = [Nx_orig, Ny_orig, Nz_orig];

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
    B.inputArgs            = inputArgs;
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
    B.truth_dose           = doseGrid;
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
    rng(seed);
    noisy = B.sensorData_resp + B.noise_amp * randn(size(B.sensorData_resp));
    deconv = real(ifft( ...
        fft(noisy, [], 2) .* B.H_conj ./ (B.H_power + B.conv_deconv_lambda), [], 2));
    sd = single(deconv);
end

function recon_dose = reconstruct_recon_dose(B, sensorData_measured)
%RECONSTRUCT_RECON_DOSE  Reconstruct a dose from a (noisy, deconvolved) sensor
%  trace using the cached forward bundle B. parfor/GPU-safe (no plotting).

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
%  Serial over criteria. Reference = reconA, target = reconB, 10% reconA cutoff.
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
%  ens_mean,ens_std  optional ensemble band for the signal series. When provided
%                    (non-empty), error bars are overlaid and the nominal curve
%                    is shown faintly.
%  null_mean,null_std  optional null-hypothesis band (recon vs an independent-
%                    noise recon of itself). Plotted in green as the no-change
%                    noise floor. A red reference line marks the 90% pass rate.

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
        h_nom = plot(crit_n(valid), pass_rates(valid), '--o', ...
            'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.0, 'MarkerSize', 4);
        legend_handles(end+1) = h_nom; %#ok<AGROW>
        legend_labels{end+1}  = 'Nominal (single draw)';

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
        legend_labels{end+1}  = sprintf('Signal: %s vs %s', label_ref, label_tgt);

        for k = 1:numel(crit_n)
            if valid(k)
                text(crit_n(k), pass_rates(k) + 1.5, sprintf('%.1f%%', pass_rates(k)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    end

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
    title(sprintf('Gamma Pass Rate vs Criterion   |   %s\n%s vs %s', ...
        patient_id, label_ref, label_tgt), 'FontWeight', 'bold', 'FontSize', 11);
    legend(legend_handles, legend_labels, 'Location', 'southeast', 'FontSize', 8);
    drawnow;
end

function plot_truth_recon_diff_axial(volA, volB, slice_ref_vol, density, sensor_mask, ...
        spacing_mm, labels, titleStr, smooth_sigma, scale_mode, clip_pct)
%PLOT_TRUTH_RECON_DIFF_AXIAL  1x3 axial figure: volA | volB | (volA - volB).
%  All inputs are on the SAME (display) grid. The axial (constant-Z, dim 3)
%  slice is taken at the max of slice_ref_vol (the truth). Panels 1 and 2 use the
%  shared jet dose rendering (render_dose_panel); panel 3 shows the signed
%  difference with a symmetric diverging colormap (blue<0<red), alpha-ramped by
%  magnitude. Density is the grayscale background; the sensor contour overlays
%  every panel.
%    labels {1x2} panel titles, e.g. {'RayStation truth','Recon CT\_1'}
%    scale_mode / clip_pct  as in render_dose_panel (relative vs absolute).

    if nargin < 9  || isempty(smooth_sigma), smooth_sigma = 0;          end
    if nargin < 10 || isempty(scale_mode),   scale_mode   = 'relative'; end
    if nargin < 11 || isempty(clip_pct),     clip_pct     = 100;        end
    use_relative = strcmpi(scale_mode, 'relative');

    gridSize = size(volA);
    [~, midx]    = max(slice_ref_vol(:));
    [~, ~, cz]   = ind2sub(gridSize, midx);

    x_ax = (1:gridSize(1)) * spacing_mm(1);
    y_ax = (1:gridSize(2)) * spacing_mm(2);

    cmap_jet  = jet(256);
    wl_center = 1050; wl_width = 350;
    wl_min    = wl_center - wl_width / 2;

    if clip_pct >= 100
        clip_str = 'max';
    else
        clip_str = sprintf('%g%%ile', clip_pct);
    end

    have_density = ~isempty(density) && isequal(size(density), gridSize);
    if have_density
        dens_sl = squeeze(density(:, :, cz))';
    else
        dens_sl = [];
    end
    if ~isempty(sensor_mask) && isequal(size(sensor_mask), gridSize)
        sensor_sl = squeeze(sensor_mask(:, :, cz))';
    else
        sensor_sl = [];
    end

    a_sl = squeeze(volA(:, :, cz))';
    b_sl = squeeze(volB(:, :, cz))';

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [60, 80, 1380, 460]);
    if use_relative
        scale_note = sprintf('panels 1-2 scaled to own %s (relative %%)', clip_str);
    else
        scale_note = sprintf('panels 1-2 shared clim (clip @ %s)', clip_str);
    end
    sgtitle(sprintf('%s\nAxial slice at truth max-dose voxel (Z=%d)  |  %s', ...
        titleStr, cz, scale_note), 'FontWeight', 'bold', 'FontSize', 11);

    % Shared reference for absolute mode (clipped max of both volumes).
    shared_max = max(robust_dose_max(volA, clip_pct), robust_dose_max(volB, clip_pct));

    % --- Panel 1: volA ---
    ax1 = subplot(1, 3, 1);
    if use_relative, rmA = robust_dose_max(volA, clip_pct); else, rmA = shared_max; end
    render_dose_panel(ax1, a_sl, dens_sl, sensor_sl, x_ax, y_ax, rmA, cmap_jet, ...
        wl_min, wl_width, use_relative, clip_str, smooth_sigma, 'X (mm)', 'Y (mm)', ...
        sprintf('%s  Axial (Z=%d)', labels{1}, cz));

    % --- Panel 2: volB ---
    ax2 = subplot(1, 3, 2);
    if use_relative, rmB = robust_dose_max(volB, clip_pct); else, rmB = shared_max; end
    render_dose_panel(ax2, b_sl, dens_sl, sensor_sl, x_ax, y_ax, rmB, cmap_jet, ...
        wl_min, wl_width, use_relative, clip_str, smooth_sigma, 'X (mm)', 'Y (mm)', ...
        sprintf('%s  Axial (Z=%d)', labels{2}, cz));

    % --- Panel 3: signed difference (volA - volB) ---
    ax3 = subplot(1, 3, 3);
    diff_sl = a_sl - b_sl;
    if smooth_sigma > 0, diff_sl = imgaussfilt(diff_sl, smooth_sigma); end
    dmax = robust_dose_max(abs(diff_sl), clip_pct);
    if ~(dmax > 0), dmax = max(abs(diff_sl(:))); end
    if ~(dmax > 0), dmax = 1; end

    if have_density
        dn     = (dens_sl - wl_min) / wl_width;
        dn     = max(0, min(1, dn));
        bg_rgb = repmat(dn, [1, 1, 3]);
    else
        bg_rgb = zeros([size(diff_sl), 3]);
    end
    image(ax3, x_ax, y_ax, bg_rgb);
    hold(ax3, 'on');

    cmap_div = blue_white_red(256);
    norm_d   = max(-1, min(1, diff_sl / dmax));
    idx      = max(1, min(256, round((norm_d + 1) / 2 * 255) + 1));
    sz       = size(diff_sl);
    diff_rgb = reshape(cmap_div(idx(:), :), [sz, 3]);
    alpha    = min(1, abs(diff_sl) / (0.30 * dmax));

    h_diff           = image(ax3, x_ax, y_ax, diff_rgb);
    h_diff.AlphaData = alpha;

    if ~isempty(sensor_sl) && any(sensor_sl(:))
        contour(ax3, x_ax, y_ax, double(sensor_sl), [0.5, 0.5], 'k-', 'LineWidth', 1.5);
    end
    hold(ax3, 'off');
    axis(ax3, 'xy'); axis(ax3, 'image');
    colormap(ax3, cmap_div);
    caxis(ax3, [-dmax, dmax]);
    cb = colorbar(ax3); cb.Label.String = 'Difference (Gy)';
    xlabel(ax3, 'X (mm)'); ylabel(ax3, 'Y (mm)');
    title(ax3, sprintf('%s - %s  Axial (Z=%d)', labels{1}, labels{2}, cz));
    drawnow;
end

function cmap = blue_white_red(n)
%BLUE_WHITE_RED  Diverging colormap: blue (negative) - white (0) - red (positive).
    if nargin < 1 || isempty(n), n = 256; end
    bot = [0.230, 0.299, 0.754];   % blue
    mid = [1.000, 1.000, 1.000];   % white
    top = [0.706, 0.016, 0.150];   % red
    h   = floor(n / 2);
    half1 = [linspace(bot(1), mid(1), h)', linspace(bot(2), mid(2), h)', ...
             linspace(bot(3), mid(3), h)'];
    half2 = [linspace(mid(1), top(1), n - h)', linspace(mid(2), top(2), n - h)', ...
             linspace(mid(3), top(3), n - h)'];
    cmap = [half1; half2];
end

function render_dose_panel(ax, dose_slice, density_slice, sensor_slice, xv, yv, ...
        row_max, cmap_jet, wl_min, wl_width, use_relative, clip_str, ...
        smooth_sigma, xlbl, ylbl, ttl)
%RENDER_DOSE_PANEL  Draw one dose panel: density grayscale background, jet dose
%  overlay (alpha-masked at 10% of row_max with a ramp), optional sensor
%  contour, and a colorbar in % (relative) or Gy (absolute).

    dose_slice = double(dose_slice);
    if smooth_sigma > 0
        dose_slice = imgaussfilt(dose_slice, smooth_sigma);
    end

    have_density = ~isempty(density_slice) && isequal(size(density_slice), size(dose_slice));
    if have_density
        dn     = (density_slice - wl_min) / wl_width;
        dn     = max(0, min(1, dn));
        bg_rgb = repmat(dn, [1, 1, 3]);
    else
        bg_rgb = zeros([size(dose_slice), 3]);
    end
    image(ax, xv, yv, bg_rgb);
    hold(ax, 'on');

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

    if ~isempty(sensor_slice) && any(sensor_slice(:))
        contour(ax, xv, yv, double(sensor_slice), [0.5, 0.5], 'r-', 'LineWidth', 2);
    end
    hold(ax, 'off');

    axis(ax, 'xy'); axis(ax, 'image');

    colormap(ax, cmap_jet);
    if use_relative
        caxis(ax, [0, 100]);
        cb = colorbar(ax);
        cb.Label.String = sprintf('Relative dose (%% of %s)', clip_str);
    else
        caxis(ax, [0, row_max]);
        cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
    end

    xlabel(ax, xlbl); ylabel(ax, ylbl);
    title(ax, ttl);
end

function plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, density, config)
%PLOT_SENSOR_DOSE_PLANES  1x3 anatomical view of sensor geometry vs dose mask.
%  Three orthogonal max-projections (coronal, sagittal, axial). CT density is the
%  grayscale background (mean-projection); the dose mask (dose >= 10% max) is a
%  semi-transparent blue region; the sensor is a solid red region.

    [Nx3, Ny3, Nz3] = size(dose_mask);

    x_ax = (1:Nx3) * spacing_mm(1);
    y_ax = (1:Ny3) * spacing_mm(2);
    z_ax = (1:Nz3) * spacing_mm(3);

    dose_cor  = squeeze(any(dose_mask,   2));
    dose_sag  = squeeze(any(dose_mask,   1));
    dose_axi  = squeeze(any(dose_mask,   3));

    sens_cor  = squeeze(any(sensor_mask, 2));
    sens_sag  = squeeze(any(sensor_mask, 1));
    sens_axi  = squeeze(any(sensor_mask, 3));

    have_density = ~isempty(density) && isequal(size(density), size(dose_mask));
    wl_center = 1050; wl_width = 350;
    wl_min    = wl_center - wl_width / 2;
    if have_density
        ct_projs = {
            squeeze(mean(double(density), 2))', ...
            squeeze(mean(double(density), 1))', ...
            squeeze(mean(double(density), 3))'
        };
    end

    figure('Name', 'Sensor Placement vs Dose Mask', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 1300, 420]);
    sgtitle(sprintf('Sensor Placement vs Dose Mask  (\\geq10%% max)   |   Sensor: %s', ...
        config.sensor_placement_method), 'FontWeight', 'bold', 'FontSize', 11);

    view_data = {
        dose_cor', sens_cor', x_ax, z_ax, 'X (mm)', 'Z (mm)', 'Coronal  (mean-proj along Y)';
        dose_sag', sens_sag', y_ax, z_ax, 'Y (mm)', 'Z (mm)', 'Sagittal  (mean-proj along X)';
        dose_axi', sens_axi', x_ax, y_ax, 'X (mm)', 'Y (mm)', 'Axial  (mean-proj along Z)';
    };

    dose_color   = [0.20, 0.50, 0.90];
    sensor_color = [0.90, 0.10, 0.10];

    for col = 1:3
        ax  = subplot(1, 3, col);
        d2d = double(view_data{col, 1});
        s2d = double(view_data{col, 2});
        xv  = view_data{col, 3};
        yv  = view_data{col, 4};

        if have_density
            dn     = (ct_projs{col} - wl_min) / wl_width;
            dn     = max(0, min(1, dn));
            bg_rgb = repmat(dn, [1, 1, 3]);
        else
            bg_rgb = ones([size(d2d), 3]);
        end
        image(ax, xv, yv, bg_rgb);
        hold(ax, 'on');

        if any(d2d(:))
            dose_rgb = repmat(reshape(dose_color, [1,1,3]), [size(d2d), 1]);
            h_dose   = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = d2d * 0.45;
        end

        if any(s2d(:))
            sens_rgb = repmat(reshape(sensor_color, [1,1,3]), [size(s2d), 1]);
            h_sens   = image(ax, xv, yv, sens_rgb);
            h_sens.AlphaData = s2d * 0.85;
        end

        hold(ax, 'off');
        axis(ax, 'xy'); axis(ax, 'image');
        xlabel(ax, view_data{col, 5}); ylabel(ax, view_data{col, 6});
        title(ax, view_data{col, 7}, 'FontSize', 10);

        hold(ax, 'on');
        p1 = patch(ax, NaN, NaN, dose_color,   'FaceAlpha', 0.45, 'EdgeColor', 'none');
        p2 = patch(ax, NaN, NaN, sensor_color, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
        hold(ax, 'off');
        if col == 3
            legend(ax, [p1, p2], {'Dose mask (>=10%)', 'Sensor'}, ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
        end
    end
    drawnow;
end

function m = robust_dose_max(d, clip_pct)
%ROBUST_DOSE_MAX  Colour-scale reference: the clip_pct percentile of the nonzero
%  dose. clip_pct >= 100 (or too few voxels) falls back to the true max.
    if nargin < 2 || isempty(clip_pct), clip_pct = 100; end
    vals = sort(double(d(d > 0)));
    if isempty(vals)
        m = 1;
        return;
    end
    if clip_pct >= 100
        m = vals(end);
    else
        k = max(1, min(numel(vals), ceil(clip_pct / 100 * numel(vals))));
        m = vals(k);
    end
    if ~(m > 0), m = vals(end); end
    if ~(m > 0), m = 1; end
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
