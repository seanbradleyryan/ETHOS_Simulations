%% =========================================================================
%  STUDY_OPTIMIZATION_SWEEPS.m
%  One-factor-at-a-time (OAT) optimization sweeps for a SINGLE beam/segment
%  pair, built on the run_standalone_comparison.m two-dose k-Wave pipeline.
%
%  WHAT IT DOES
%  ------------
%  1. Runs the full two-dose (reference-CT + counterpart-CT) forward +
%     reconstruction pipeline ONCE with all-default parameters. This is the
%     REFERENCE run used to score every sweep.
%  2. For each optimization variable, sweeps that variable over a list of
%     values while holding EVERY other variable at its default (classic OAT
%     sensitivity sweep). Each (variable, value) combination is one "run".
%  3. Every run is executed inside a try/catch so one failure never aborts the
%     sweep, and all runs are dispatched through a parfor so independent
%     configurations run in parallel (one GPU per worker when available).
%  4. Every run is scored against the reference and ALL gamma analyses from
%     run_standalone_comparison are recorded (even for unsuccessful runs, for
%     insight) into a human-readable .txt report plus a companion .mat.
%
%  SUCCESS CRITERIA (evaluated at the primary 3%/3mm gamma criterion)
%  ------------------------------------------------------------------
%  A run is SUCCESSFUL if EITHER:
%    (1) gamma(default recon vs swept recon) >= success_stability_pct  AND
%        runtime <= default runtime            (reproduces the trusted default
%                                                reconstruction, no slower); OR
%    (2) gamma(swept recon vs its truth) > gamma(default recon vs its truth) AND
%        runtime <= default runtime            (actually MORE accurate, no
%                                                slower).
%  "recon" here is recon1 = the reference-CT, full-access reconstruction.
%
%  SWEPT VARIABLES
%  ---------------
%    downscale_factor        : scales sim dx,dy,dz up by 2^p for
%                              p = 1/4,1/3,1/2,1,2,3. The coarse-grid recon is
%                              resampled BACK to the native dose dimensions
%                              before gamma so pass rates stay comparable.
%    pml_size                : reduced pointwise 10-1..10-5  (9,8,7,6,5).
%    cfl_number              : 0.1,0.2,0.3,0.4,0.5.
%    num_time_reversal_iter  : 3,5,7,10.
%    Nt_scaling  (ADDED)     : 4,6,8,10,12. Air (~343 m/s) sets minC and
%                              inflates the recording length Nt; a larger
%                              divisor trims the recording -> less compute per
%                              k-Wave call with little accuracy loss as long as
%                              all signal reaching the sensor is still captured.
%    convergence_tol (ADDED) : 1e-3,3e-3,1e-2,3e-2,1e-1. The residual-corrected
%                              time reversal often converges well before the
%                              iteration cap; a looser tolerance lets it stop
%                              early -> fewer k-Wave calls, minimal accuracy hit.
%
%  NOTE: to keep the sweep clean, the electronic-noise draw is seeded
%  (CONFIG.rng_seed) so the stability comparison reflects the swept parameter
%  and not run-to-run noise variance.
%
%  Remote/HIPAA execution only. Do not run locally.
%% =========================================================================

clear; clc; close all;

%% ========================= BASE CONFIGURATION ===========================
%  These are the DEFAULT values. Every sweep changes exactly one field and
%  leaves the rest at these values. The default run uses all of them.

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';

% Single beam/segment pair to study.
CONFIG.dose_filename = 'dose_1194203_Session_1_reference_CT_3_B15_112.mat';

CONFIG.ct_pair      = [1, 3];
CONFIG.reference_ct = 1;

CONFIG.dose_file_override = '';
CONFIG.cbct_file_override = '';

CONFIG.sensor_placement_method = 'determine_sensor_mask';
CONFIG.sensor_x_index = 2;
CONFIG.sensor_y_index = 4;

CONFIG.elements_per_side = 32;
CONFIG.element_pitch_mm  = 4.35;
CONFIG.element_size_mm   = 3.65;

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
CONFIG.meterset               = 140;
CONFIG.pml_size               = 10;      % DEFAULT  (swept: 9 8 7 6 5)
CONFIG.cfl_number             = 0.3;     % DEFAULT  (swept: .1 .2 .3 .4 .5)
CONFIG.Nt_scaling             = 6;       % DEFAULT  (swept: 4 6 8 10 12)
CONFIG.use_gpu                = true;
CONFIG.correction_factor      = 20;
CONFIG.use_pressure_scale_correction = false;
CONFIG.mask_recon_to_dose_region     = true;

CONFIG.reconstruction_method = 'tr';     % this study only sweeps the TR pipeline
CONFIG.num_time_reversal_iter = 10;      % DEFAULT  (swept: 3 5 7 10)
CONFIG.convergence_tol        = 1e-3;    % DEFAULT  (swept: 1e-3 3e-3 1e-2 3e-2 1e-1)

CONFIG.convolution_kernel  = 4e-6;
CONFIG.conv_noise_level    = 0.125;
CONFIG.conv_deconv_lambda  = 1e-4;

CONFIG.downscale_factor = 1;             % DEFAULT  (swept: 2^[1/4 1/3 1/2 1 2 3])
CONFIG.use_grid_padding = true;

% Study-specific settings (not part of the standalone pipeline).
CONFIG.rng_seed             = 42;        % deterministic electronic-noise draw
CONFIG.success_stability_pct = 98;       % criterion-1 threshold (%)
CONFIG.runtime_tol          = 0.02;      % allow 2% timing jitter for "equal or lower"

% Least-squares relative normalization (scheme from study_pass_rates_allsegments):
% before gamma, each recon is rescaled by the scalar gain that best matches it to
% its OWN-CT truth over that truth's 10% region (recon1->dose1, recon2->dose2),
% cancelling the stored absolute scale (CONFIG.correction_factor). The truths are
% never rescaled. This makes gamma a consistent, scale-invariant metric across
% every swept configuration.
CONFIG.normalize = true;

% Plotting / saving OFF: this runs headless inside a parfor.
CONFIG.plot_results = false;
CONFIG.save_results = false;
CONFIG.log_gamma    = false;

%% ========================= SWEEP / OUTPUT SETTINGS =======================

% Parallelization. Set use_parallel=false to run the sweep serially (useful
% for a single-GPU box, where full k-Wave sims contend for GPU memory). When
% true, num_workers defaults to the number of visible GPUs.
SWEEP.use_parallel = true;
SWEEP.num_workers  = [];   % [] -> #GPUs (GPU) or feature('numcores') (CPU)

% Primary gamma criterion used for the pass/fail decision.
SWEEP.gamma_criteria = {10, 10, '10%/10mm'; 5, 5, '5%/5mm'; 3, 3, '3%/3mm'};
SWEEP.primary_crit_idx = 3;   % 3%/3mm

% Output files (written to working_dir).
run_stamp = datestr(now, 'yyyymmdd_HHMMSS');
SWEEP.report_txt = sprintf('optimization_sweep_report_%s.txt', run_stamp);
SWEEP.report_mat = sprintf('optimization_sweep_results_%s.mat', run_stamp);

%% ===================== RESOLVE DOSE PAIR & CBCT PATHS ====================
%  (Identical resolution logic to run_standalone_comparison: locate dose A,
%   derive counterpart B on the other CT, pair each with its own CBCT.)

if ~isempty(CONFIG.dose_file_override)
    dose_filepath_A = CONFIG.dose_file_override;
    [processed_dir, baseA, extA] = fileparts(dose_filepath_A);
    dose_basename_A = [baseA, extA];
else
    processed_dir   = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    dose_basename_A = CONFIG.dose_filename;
    dose_filepath_A = fullfile(processed_dir, dose_basename_A);
end

ct_tok = regexp(dose_basename_A, '_CT_(\d+)_', 'tokens', 'once');
if isempty(ct_tok)
    error('study_optimization_sweeps:NoCTtoken', ...
        'Dose filename "%s" has no _CT_k token.', dose_basename_A);
end
ct_a = str2double(ct_tok{1});
if ~ismember(ct_a, CONFIG.ct_pair)
    error('study_optimization_sweeps:CTnotInPair', ...
        'Dose CT index %d is not in CONFIG.ct_pair = [%s].', ct_a, num2str(CONFIG.ct_pair));
end
ct_b_all = CONFIG.ct_pair(CONFIG.ct_pair ~= ct_a);
ct_b     = ct_b_all(1);
ct_list  = [ct_a, ct_b];

if ~ismember(CONFIG.reference_ct, CONFIG.ct_pair)
    error('study_optimization_sweeps:BadReferenceCT', ...
        'CONFIG.reference_ct = %d must be one of CONFIG.ct_pair = [%s].', ...
        CONFIG.reference_ct, num2str(CONFIG.ct_pair));
end
recon_cbct_filepath = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', CONFIG.reference_ct));

dose_basename_B = regexprep(dose_basename_A, '_CT_\d+_', sprintf('_CT_%d_', ct_b));
dose_filepath_B = fullfile(processed_dir, dose_basename_B);

if ~isempty(CONFIG.cbct_file_override)
    cbct_filepath_A = CONFIG.cbct_file_override;
else
    cbct_filepath_A = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_a));
end
cbct_filepath_B = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_b));

PATHS = struct();
PATHS.dose        = {dose_filepath_A, dose_filepath_B};
PATHS.cbct        = {cbct_filepath_A, cbct_filepath_B};
PATHS.ct_list     = ct_list;
PATHS.reference_ct = CONFIG.reference_ct;
PATHS.recon_cbct  = recon_cbct_filepath;
PATHS.criteria    = SWEEP.gamma_criteria;
PATHS.primary     = SWEEP.primary_crit_idx;

%% ========================= LOAD PLAN BEAM METADATA =======================

if ~isfield(CONFIG, 'beam_metadata') || isempty(CONFIG.beam_metadata)
    metadata_filepath = fullfile(processed_dir, 'metadata.mat');
    if isfile(metadata_filepath)
        try
            md = load(metadata_filepath, 'metadata');
            if isfield(md, 'metadata') && isfield(md.metadata, 'beam_metadata') ...
                    && ~isempty(md.metadata.beam_metadata)
                CONFIG.beam_metadata = md.metadata.beam_metadata;
                fprintf('  Loaded beam_metadata for %d beams.\n', numel(CONFIG.beam_metadata));
            else
                fprintf('  [WARN] metadata.mat present but no beam_metadata field.\n');
                CONFIG.beam_metadata = [];
            end
        catch ME
            fprintf('  [WARN] Failed to load metadata.mat: %s\n', ME.message);
            CONFIG.beam_metadata = [];
        end
    else
        fprintf('  [WARN] No metadata.mat in processed_dir; exclusion zone will be empty.\n');
        CONFIG.beam_metadata = [];
    end
end

% Shared tissue tables (same builder the pipeline uses).
CONFIG.tissue_tables = define_tissue_tables();
CONFIG.tissue_tables.uniform = struct( ...
    'density',     CONFIG.uniform_density, ...
    'sound_speed', CONFIG.uniform_sound_speed, ...
    'alpha_coeff', CONFIG.uniform_alpha_coeff, ...
    'alpha_power', 1.1, ...
    'gruneisen',   CONFIG.uniform_gruneisen);

%% ========================= DEFINE THE SWEEPS ============================
%  Each sweep: a name, the CONFIG field it changes, and the list of values.
%  Everything else stays at the default CONFIG above.

SWEEPS = {};
SWEEPS{end+1} = make_sweep('voxel_scale (dx,dy,dz x 2^p)', 'downscale_factor', ...
                           2.^[1/4, 1/3, 1/2, 1, 2, 3]);
SWEEPS{end+1} = make_sweep('pml_size (reduce)',            'pml_size', ...
                           10 + [-1, -2, -3, -4, -5]);
SWEEPS{end+1} = make_sweep('cfl_number',                   'cfl_number', ...
                           [0.1, 0.2, 0.3, 0.4, 0.5]);
SWEEPS{end+1} = make_sweep('num_time_reversal_iter',       'num_time_reversal_iter', ...
                           [3, 5, 7, 10]);
SWEEPS{end+1} = make_sweep('Nt_scaling (recording length)','Nt_scaling', ...
                           [4, 6, 8, 10, 12]);
SWEEPS{end+1} = make_sweep('convergence_tol (TR early stop)','convergence_tol', ...
                           [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]);

% Flatten into a job list (one job per value).
jobs = struct('sweep', {}, 'field', {}, 'value', {}, 'label', {});
for s = 1:numel(SWEEPS)
    sw = SWEEPS{s};
    for v = 1:numel(sw.values)
        jobs(end+1) = struct('sweep', sw.name, 'field', sw.field, ...
            'value', sw.values(v), 'label', sw.labels{v}); %#ok<SAGROW>
    end
end
nJobs = numel(jobs);

fprintf('=========================================================\n');
fprintf('  Optimization sweep harness\n');
fprintf('  Patient %s / %s   dose: %s\n', CONFIG.patient_id, CONFIG.session, CONFIG.dose_filename);
fprintf('  %d sweep variables, %d total sweep runs (+1 default).\n', numel(SWEEPS), nJobs);
fprintf('=========================================================\n\n');

%% ========================= DEFAULT / REFERENCE RUN ======================
%  Run the all-default configuration first (serially). Its runtime and recon
%  become the reference every sweep run is scored against.

fprintf('[REFERENCE] Running all-default configuration...\n');
try
    R_default = run_one_config(CONFIG, PATHS, []);
catch ME
    error('study_optimization_sweeps:DefaultFailed', ...
        'The default reference run failed (%s): %s\nCannot score sweeps without a baseline.', ...
        ME.identifier, ME.message);
end
fprintf('[REFERENCE] runtime = %.1f s, recon1-vs-truth (%s) = %.2f%%\n\n', ...
    R_default.runtime_sec, SWEEP.gamma_criteria{SWEEP.primary_crit_idx, 3}, ...
    R_default.pass.recon1_vs_dose1(SWEEP.primary_crit_idx));

% Reference bundle broadcast into the parfor.
DEFAULT_REF = struct();
DEFAULT_REF.recon1          = single(R_default.recon1);
DEFAULT_REF.spacing         = R_default.spacing;
DEFAULT_REF.runtime         = R_default.runtime_sec;
DEFAULT_REF.recon1_vs_dose1 = R_default.pass.recon1_vs_dose1;
DEFAULT_REF.primary         = SWEEP.primary_crit_idx;
DEFAULT_REF.success_stability_pct = CONFIG.success_stability_pct;
DEFAULT_REF.runtime_tol     = CONFIG.runtime_tol;

%% ========================= PARALLEL POOL / GPUs =========================

num_gpus = 0;
if CONFIG.use_gpu
    try, num_gpus = gpuDeviceCount; catch, num_gpus = 0; end
end

if SWEEP.use_parallel && nJobs > 1
    if isempty(SWEEP.num_workers)
        if num_gpus > 0
            SWEEP.num_workers = num_gpus;
        else
            SWEEP.num_workers = feature('numcores');
        end
    end
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool(SWEEP.num_workers);
    end
    % Pin each worker to a distinct GPU (round-robin) when GPUs are present.
    if num_gpus > 1
        spmd
            gpuDevice(mod(labindex - 1, num_gpus) + 1);
        end
    end
    par_M = pool.NumWorkers;   % cap the parfor at the pool size
    fprintf('[PARALLEL] Pool of %d workers, %d GPU(s).\n\n', pool.NumWorkers, num_gpus);
else
    % par_M = 0 forces the parfor to run as a serial loop on the client with
    % NO pool auto-started (important for single-GPU boxes to avoid contention).
    par_M = 0;
    fprintf('[SERIAL] Running sweeps serially (no pool).\n\n');
end

%% ========================= RUN ALL SWEEP JOBS ==========================
%  Every job is wrapped in try/catch: a failed configuration is recorded as
%  ok=false with its error message, never aborting the sweep.

job_results = cell(nJobs, 1);

parfor (j = 1:nJobs, par_M)
    cfg = CONFIG;
    cfg.(jobs(j).field) = jobs(j).value;

    Rj = blank_result();
    Rj.sweep = jobs(j).sweep;
    Rj.field = jobs(j).field;
    Rj.value = jobs(j).value;
    Rj.label = jobs(j).label;

    t_wall = tic;
    try
        Rrun = run_one_config(cfg, PATHS, DEFAULT_REF); %#ok<PFBNS>
        f = fieldnames(Rrun);
        for k = 1:numel(f)
            Rj.(f{k}) = Rrun.(f{k});
        end
        Rj.ok = true;
    catch ME
        Rj.ok = false;
        Rj.error = sprintf('%s | %s', ME.identifier, ME.message);
        fprintf('   [FAIL] %-28s = %-10s : %s\n', jobs(j).field, jobs(j).label, ME.message);
    end
    Rj.wall_sec = toc(t_wall);
    job_results{j} = Rj;
end

% Convert the cell array back into a struct array.
RESULTS = [job_results{:}];

fprintf('\n[DONE] %d/%d sweep runs completed successfully.\n', ...
    sum([RESULTS.ok]), nJobs);

%% ========================= WRITE REPORT ================================

report_path = SWEEP.report_txt;
if isempty(fileparts(report_path))
    report_path = fullfile(CONFIG.working_dir, report_path);
end
mat_path = SWEEP.report_mat;
if isempty(fileparts(mat_path))
    mat_path = fullfile(CONFIG.working_dir, mat_path);
end

write_report(report_path, CONFIG, SWEEP, PATHS, R_default, RESULTS, SWEEPS);
fprintf('\n[REPORT] Written: %s\n', report_path);

% Companion .mat (drop the heavy reference recon volume from the saved copy).
R_default_light = R_default;
R_default_light.recon1 = [];
save(mat_path, 'CONFIG', 'SWEEP', 'PATHS', 'R_default_light', 'RESULTS', '-v7.3');
fprintf('[REPORT] Saved:   %s\n', mat_path);

fprintf('\nOptimization sweep complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function sw = make_sweep(name, field, values)
%MAKE_SWEEP  Build one sweep spec with auto value labels.
    sw = struct();
    sw.name   = name;
    sw.field  = field;
    sw.values = values;
    sw.labels = cell(1, numel(values));
    for i = 1:numel(values)
        v = values(i);
        if v ~= 0 && (abs(v) < 1e-2 || abs(v) >= 1e4)
            sw.labels{i} = sprintf('%.3g', v);
        elseif v == round(v)
            sw.labels{i} = sprintf('%d', v);
        else
            sw.labels{i} = sprintf('%.4g', v);
        end
    end
end


function R = blank_result()
%BLANK_RESULT  Struct skeleton so the parfor sliced output stays homogeneous.
    R = struct();
    R.ok          = false;
    R.error       = '';
    R.sweep       = '';
    R.field       = '';
    R.value       = NaN;
    R.label       = '';
    R.wall_sec    = NaN;
    R.runtime_sec = NaN;
    R.fwd_time    = [NaN NaN];
    R.tr_time     = [NaN NaN];
    R.iters       = [NaN NaN];
    R.dims_native = [NaN NaN NaN];
    R.pass        = struct('dose1_vs_dose2', [], 'recon1_vs_dose1', [], ...
                           'recon2_vs_dose2', [], 'recon2_vs_dose1', []);
    R.spacing     = [];
    R.stability   = [];
    R.crit1       = false;
    R.crit2       = false;
    R.success     = false;
    R.recon1      = [];
end


function R = run_one_config(cfg, paths, default_ref)
%RUN_ONE_CONFIG  Full two-dose forward+recon pipeline + all gamma analyses.
%  Runs simulate_field_once on both the reference-CT dose and the counterpart
%  dose, then computes the four run_standalone_comparison gamma comparisons on
%  the NATIVE dose grid. When default_ref is supplied (sweep runs) it also
%  computes the default-vs-swept stability gamma and the pass/fail decision.

    criteria = paths.criteria;
    primary  = paths.primary;

    % --- run both doses ---
    o1 = simulate_field_once(paths.dose{1}, paths.cbct{1}, paths.ct_list(1), ...
                             paths.reference_ct, paths.recon_cbct, cfg);
    o2 = simulate_field_once(paths.dose{2}, paths.cbct{2}, paths.ct_list(2), ...
                             paths.reference_ct, paths.recon_cbct, cfg);

    % --- map to reference-CT (idx1) vs counterpart-CT (idx2) ---
    idx1 = find(paths.ct_list == paths.reference_ct, 1);
    idx2 = find(paths.ct_list ~= paths.reference_ct, 1);
    if isempty(idx2)
        tmp = setdiff([1, 2], idx1); idx2 = tmp(1);
    end
    outs = {o1, o2};
    oref = outs{idx1};
    ooth = outs{idx2};

    dose1  = double(oref.dose_native);
    recon1 = double(oref.recon_native);
    dose2  = double(ooth.dose_native);
    recon2 = double(ooth.recon_native);
    spacing = oref.spacing_native;

    if ~isequal(size(dose1), size(dose2))
        error('run_one_config:GridMismatch', ...
            'Native grids differ: [%s] vs [%s].', ...
            num2str(size(dose1)), num2str(size(dose2)));
    end

    % --- Least-squares relative normalization (recon -> own-CT truth) ---
    % Applied BEFORE any gamma so every comparison (and the stored reference
    % recon used for the stability comparison) is on the same scale-invariant
    % footing. Cancels CONFIG.correction_factor; truths (dose1/dose2) untouched.
    if isfield(cfg, 'normalize') && cfg.normalize
        recon1 = recon1 * least_squares_gain(dose1, recon1);
        recon2 = recon2 * least_squares_gain(dose2, recon2);
    end

    runtime = sum([o1.fwd_time, o1.tr_time, o2.fwd_time, o2.tr_time]);

    % --- the four standalone-comparison gamma analyses (native grid) ---
    pass = struct();
    pass.dose1_vs_dose2  = gamma_pass_rates(dose1, dose2,  spacing, criteria);
    pass.recon1_vs_dose1 = gamma_pass_rates(dose1, recon1, spacing, criteria);
    pass.recon2_vs_dose2 = gamma_pass_rates(dose2, recon2, spacing, criteria);
    pass.recon2_vs_dose1 = gamma_pass_rates(dose1, recon2, spacing, criteria);

    R = struct();
    R.ok          = true;
    R.runtime_sec = runtime;
    R.fwd_time    = [o1.fwd_time, o2.fwd_time];
    R.tr_time     = [o1.tr_time,  o2.tr_time];
    R.iters       = [o1.num_iters_done, o2.num_iters_done];
    R.dims_native = size(dose1);
    R.pass        = pass;
    R.spacing     = spacing;
    R.recon1      = single(recon1);   % reference run keeps this; sweep runs drop it below
    R.stability   = [];
    R.crit1       = false;
    R.crit2       = false;
    R.success     = false;

    % --- scoring (sweep runs only) ---
    if ~isempty(default_ref)
        r1 = recon1;
        if ~isequal(size(r1), size(default_ref.recon1))
            r1 = max(imresize3(r1, size(default_ref.recon1)), 0);
        end
        R.stability = gamma_pass_rates(double(default_ref.recon1), r1, ...
                                       default_ref.spacing, criteria);

        rt_ok = runtime <= default_ref.runtime * (1 + default_ref.runtime_tol);
        R.crit1 = (R.stability(primary) >= default_ref.success_stability_pct) && rt_ok;
        R.crit2 = (pass.recon1_vs_dose1(primary) > default_ref.recon1_vs_dose1(primary)) && rt_ok;
        R.success = R.crit1 || R.crit2;

        R.recon1 = [];   % free memory: only the reference recon is retained
    end
end


function pr = gamma_pass_rates(ref_vol, tgt_vol, spacing_mm, criteria)
%GAMMA_PASS_RATES  Lean gamma pass-rate vector (no maps retained).
%  ref_vol sets the 10% low-dose evaluation mask. Mirrors run_gamma_pair's
%  CalcGamma call exactly (global gamma, DTA limit 2*dta, restrict=1).
    nCrit = size(criteria, 1);
    pr = nan(1, nCrit);

    if isempty(ref_vol) || isempty(tgt_vol) || max(ref_vol(:)) <= 0
        return;
    end

    ref_struct.start = [0, 0, 0]; ref_struct.width = spacing_mm; ref_struct.data = double(ref_vol);
    tgt_struct.start = [0, 0, 0]; tgt_struct.width = spacing_mm; tgt_struct.data = double(tgt_vol);

    eval_mask = ref_vol >= 0.10 * max(ref_vol(:));
    if ~any(eval_mask(:)), return; end

    for gc = 1:nCrit
        pct = criteria{gc, 1};
        dta = criteria{gc, 2};
        try
            gmap = CalcGamma(ref_struct, tgt_struct, pct, dta, ...
                'local', 0, 'limit', dta * 2, 'restrict', 1);
            pr(gc) = 100 * mean(gmap(eval_mask) <= 1);
        catch
            pr(gc) = NaN;
        end
    end
end


function g = least_squares_gain(rs_truth, recon)
%LEAST_SQUARES_GAIN  Scalar gain aligning a recon to its RS truth (relative norm).
%  g = sum(rs.*recon)/sum(recon.^2) over the truth's 10% low-dose region; the g
%  minimizing ||rs - g*recon||^2 there. Falls back to 1 for empty/zero inputs.
%  (Scheme from study_pass_rates_allsegments.m.)
    rs_truth = double(rs_truth);
    recon    = double(recon);
    if max(rs_truth(:)) > 0
        mask = rs_truth >= 0.10 * max(rs_truth(:));
    else
        mask = true(size(rs_truth));
    end
    r = recon(mask);
    denom = sum(r .^ 2);
    if denom > 0
        g = sum(rs_truth(mask) .* r) / denom;
    else
        g = 1;
    end
end


function out = simulate_field_once(dose_filepath, cbct_filepath, ct_this, ...
        reference_ct, recon_cbct_filepath, CONFIG)
%SIMULATE_FIELD_ONCE  Single-dose forward + TR reconstruction.
%  Runs the identical forward + iterative time-reversal pipeline as
%  run_standalone_comparison for ONE dose file, headless (no plots), and
%  returns the reconstructed dose RESAMPLED BACK to the native dose
%  dimensions (with native spacing) so gamma analyses are comparable across
%  every sweep configuration, including the dx/dy/dz (downscale) sweep.
%
%  OUTPUT struct fields:
%    .recon_native    reconstructed dose at native dims [Nx0 Ny0 Nz0]
%    .dose_native     body-masked truth dose at native dims (gamma reference)
%    .spacing_native  native voxel spacing [dx dy dz] mm (for gamma DTA)
%    .fwd_time        forward-sim wall time (s)
%    .tr_time         time-reversal wall time (s)
%    .num_iters_done  TR iterations actually run

    out = struct();
    blind_recon = (ct_this ~= reference_ct);

    %% ---- LOAD DOSE ----
    if ~isfile(dose_filepath)
        error('simulate_field_once:NoDose', 'Dose file not found: %s', dose_filepath);
    end
    dose_data   = load(dose_filepath);
    dose_fields = fieldnames(dose_data);

    fd_gantry_angle   = [];
    step15_spacing_mm = [];

    if isfield(dose_data, 'field_dose')
        fd = dose_data.field_dose;
        if ~isfield(fd, 'dose_Gy'), error('field_dose struct missing dose_Gy field.'); end
        if (isfield(fd, 'is_sparse') && fd.is_sparse) || issparse(fd.dose_Gy)
            if ~isfield(fd, 'dose_dims'), error('field_dose.dose_dims missing.'); end
            doseGrid = reshape(full(fd.dose_Gy), fd.dose_dims);
        else
            doseGrid = double(fd.dose_Gy);
        end
        if isfield(fd, 'spacing') && ~isempty(fd.spacing)
            step15_spacing_mm = fd.spacing(:)';
        end
        if isfield(fd, 'meterset') && ~isempty(fd.meterset) && fd.meterset > 0
            CONFIG.meterset = fd.meterset;
        end
        if isfield(fd, 'gantry_angle')
            fd_gantry_angle = fd.gantry_angle;
        end
    elseif isfield(dose_data, 'total_rs_dose_sparse')
        if ~isfield(dose_data, 'total_rs_dose_dims'), error('total_rs_dose_dims missing.'); end
        doseGrid = reshape(full(dose_data.total_rs_dose_sparse), dose_data.total_rs_dose_dims);
    elseif isfield(dose_data, 'total_rs_dose')
        doseGrid = dose_data.total_rs_dose;
    elseif isfield(dose_data, 'dose_Gy')
        doseGrid = dose_data.dose_Gy;
    elseif numel(dose_fields) == 1
        doseGrid = dose_data.(dose_fields{1});
    else
        error('Cannot auto-detect dose variable. Fields: %s', strjoin(dose_fields, ', '));
    end
    doseGrid = double(doseGrid);
    if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
        error('Dose data must be a 3D numeric array.');
    end

    gridSize = size(doseGrid);
    Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);

    %% ---- LOAD CBCT ----
    if ~isfile(cbct_filepath)
        error('simulate_field_once:NoCBCT', 'CBCT file not found: %s', cbct_filepath);
    end
    cbct_data = load(cbct_filepath);
    if isfield(cbct_data, 'CBCT1_resampled')
        sct = cbct_data.CBCT1_resampled;
    elseif isfield(cbct_data, 'CBCT3_resampled')
        sct = cbct_data.CBCT3_resampled;
    else
        error('CBCT1_resampled / CBCT3_resampled not found in %s', cbct_filepath);
    end
    if ~isfield(sct, 'cubeHU'), error('CBCT struct missing cubeHU.'); end

    if isfield(sct, 'spacing') && ~isempty(sct.spacing)
        spacing_mm = sct.spacing(:)';
    elseif ~isempty(step15_spacing_mm)
        spacing_mm = step15_spacing_mm;
    else
        error('No spacing available (CBCT or field_dose).');
    end

    if ~isequal(gridSize, size(sct.cubeHU))
        error('Dose grid [%s] != CBCT grid [%s].', ...
            num2str(gridSize), num2str(size(sct.cubeHU)));
    end

    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    % Native reference (BEFORE any downscale) used for gamma comparisons.
    dose_native  = doseGrid;
    spacing_native = spacing_mm;
    nativeDims   = [Nx, Ny, Nz];

    %% ---- GRID DOWNSCALING (scales dx,dy,dz up by downscale_factor) ----
    if CONFIG.downscale_factor ~= 1
        df     = CONFIG.downscale_factor;
        new_Nx = max(2, round(Nx / df));
        new_Ny = max(2, round(Ny / df));
        new_Nz = max(2, round(Nz / df));

        doseGrid   = max(imresize3(doseGrid, [new_Nx, new_Ny, new_Nz]), 0);
        sct.cubeHU = imresize3(sct.cubeHU, [new_Nx, new_Ny, new_Nz]);
        if isfield(sct, 'bodyMask')
            sct.bodyMask = imresize3(single(sct.bodyMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        if isfield(sct, 'couchMask')
            sct.couchMask = imresize3(single(sct.couchMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        if isfield(sct, 'cubeDensity')
            sct.cubeDensity = imresize3(sct.cubeDensity, [new_Nx, new_Ny, new_Nz]);
        end

        spacing_mm = spacing_mm .* ([Nx, Ny, Nz] ./ [new_Nx, new_Ny, new_Nz]);
        sct.spacing = spacing_mm;
        Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
        gridSize = [Nx, Ny, Nz];
    end
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;

    % Downscaled-native dims: the anatomy block size before FFT padding /
    % sensor grid expansion (equals nativeDims when downscale_factor == 1).
    dsDims = [Nx, Ny, Nz];

    %% ---- ACOUSTIC MEDIUM (forward = TRUE geometry) ----
    medium = create_acoustic_medium(sct, CONFIG);
    medium = apply_standalone_medium_overrides(medium, sct, CONFIG);

    if blind_recon
        if ~isfile(recon_cbct_filepath)
            error('Reference CBCT for recon not found: %s', recon_cbct_filepath);
        end
        rec_cbct = load(recon_cbct_filepath);
        if isfield(rec_cbct, 'CBCT1_resampled')
            sct_recon = rec_cbct.CBCT1_resampled;
        elseif isfield(rec_cbct, 'CBCT3_resampled')
            sct_recon = rec_cbct.CBCT3_resampled;
        else
            error('CBCT1_resampled / CBCT3_resampled not found in %s', recon_cbct_filepath);
        end
        if ~isfield(sct_recon, 'cubeHU'), error('Reference CBCT missing cubeHU.'); end
        if ~isequal(size(sct_recon.cubeHU), nativeDims) && CONFIG.downscale_factor == 1
            error('Reference CBCT grid [%s] != dose grid [%s].', ...
                num2str(size(sct_recon.cubeHU)), num2str(nativeDims));
        end
        if CONFIG.downscale_factor ~= 1
            sct_recon.cubeHU = imresize3(sct_recon.cubeHU, size(sct.cubeHU));
            if isfield(sct_recon, 'bodyMask')
                sct_recon.bodyMask = imresize3(single(sct_recon.bodyMask), size(sct.cubeHU), 'nearest') > 0.5;
            end
            if isfield(sct_recon, 'couchMask')
                sct_recon.couchMask = imresize3(single(sct_recon.couchMask), size(sct.cubeHU), 'nearest') > 0.5;
            end
        end
        sct_recon.spacing = spacing_mm;
        medium_recon = create_acoustic_medium(sct_recon, CONFIG);
        medium_recon = apply_standalone_medium_overrides(medium_recon, sct_recon, CONFIG);
    else
        medium_recon = medium;
    end

    %% ---- INITIAL PRESSURE p0 ----
    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end
    meterset       = CONFIG.meterset;
    num_pulses     = ceil(meterset / CONFIG.dose_per_pulse_cGy);
    dose_per_pulse = doseGrid / num_pulses;

    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    p0 = smooth(p0);

    doseThreshold = 0.01 * max(doseGrid(:));
    doseMask      = doseGrid > doseThreshold;
    if ~any(doseMask(:)) || max(p0(:)) == 0
        error('No significant dose or zero initial pressure.');
    end

    %% ---- OPTIMAL GRID PADDING ----
    Nx_orig = Nx; Ny_orig = Ny; Nz_orig = Nz;
    gridSize_orig = gridSize;
    medium_orig   = medium;

    % Core index ranges that locate the downscaled-native anatomy block within
    % the final cropped recon grid. Overwritten if the sensor expands the grid.
    core_xr = 1:dsDims(1);
    core_yr = 1:dsDims(2);
    core_zr = 1:dsDims(3);

    if CONFIG.use_grid_padding
        Nx_pad = find_optimal_kwave_size(Nx, CONFIG.pml_size);
        Ny_pad = find_optimal_kwave_size(Ny, CONFIG.pml_size);
        Nz_pad = find_optimal_kwave_size(Nz, CONFIG.pml_size);
    else
        Nx_pad = Nx; Ny_pad = Ny; Nz_pad = Nz;
    end

    did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
    if did_pad
        [medium, p0] = pad_forward(medium, p0, Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);
        if blind_recon
            medium_recon = pad_medium_to(medium_recon, [Nx_pad, Ny_pad, Nz_pad]);
        else
            medium_recon = medium;
        end
        Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
        gridSize = [Nx, Ny, Nz];
    end

    %% ---- SENSOR PLACEMENT ----
    sensor      = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    switch CONFIG.sensor_placement_method
        case 'full_plane_anterior'
            sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
        case 'full_plane_lateral'
            sensor.mask(:, CONFIG.sensor_y_index, :) = 1;
        case 'determine_sensor_mask'
            sct_for_sensor = sct;
            if ~isfield(sct_for_sensor, 'couchMask')
                sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
            end
            if ~isfield(sct_for_sensor, 'origin')
                sct_for_sensor.origin = [0, 0, 0];
            end
            sct_for_sensor.spacing = spacing_mm;

            field_dose_for_sensor = struct();
            field_dose_for_sensor.dose_Gy = doseGrid;
            if ~isempty(fd_gantry_angle)
                field_dose_for_sensor.gantry_angle = fd_gantry_angle;
            else
                field_dose_for_sensor.gantry_angle = 0;
            end
            field_dose_for_sensor.origin     = sct_for_sensor.origin;
            field_dose_for_sensor.spacing    = spacing_mm;
            field_dose_for_sensor.dimensions = [Nx_orig, Ny_orig, Nz_orig];

            beam_meta = [];
            if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
                beam_meta = CONFIG.beam_metadata;
            end

            [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask( ...
                sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);

            % --- GRID EXPANSION HANDLING (water-pad + re-FFT-pad) ---
            gp = sensor_info_orig.grid_pad;
            if gp.expanded
                density_unp    = medium.density(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                soundSpeed_unp = medium.sound_speed(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                if numel(medium.alpha_coeff) > 1
                    alphaCoeff_unp = medium.alpha_coeff(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                else
                    alphaCoeff_unp = medium.alpha_coeff;
                end
                gruneisen_unp  = medium.gruneisen(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                p0_unp         = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

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

                medium_orig = struct('density', density_exp, 'sound_speed', soundSpeed_exp, ...
                    'alpha_coeff', alphaCoeff_exp, 'gruneisen', gruneisen_exp);

                if blind_recon
                    medium_recon = crop_medium(medium_recon, [Nx_orig, Ny_orig, Nz_orig]);
                    medium_recon = expand_medium(medium_recon, [Nx_exp, Ny_exp, Nz_exp], xr, yr, zr);
                end

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

                % The anatomy block now lives at (xr,yr,zr) in the expanded grid,
                % which survives the FFT-pad -> crop step unchanged. Record it as
                % the core to extract from the final recon.
                core_xr = xr; core_yr = yr; core_zr = zr;

                Nx_orig = Nx_exp; Ny_orig = Ny_exp; Nz_orig = Nz_exp;
                gridSize_orig = [Nx_orig, Ny_orig, Nz_orig];

                if CONFIG.use_grid_padding
                    Nx_pad2 = find_optimal_kwave_size(Nx_orig, CONFIG.pml_size);
                    Ny_pad2 = find_optimal_kwave_size(Ny_orig, CONFIG.pml_size);
                    Nz_pad2 = find_optimal_kwave_size(Nz_orig, CONFIG.pml_size);
                else
                    Nx_pad2 = Nx_orig; Ny_pad2 = Ny_orig; Nz_pad2 = Nz_orig;
                end

                if ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], [Nx_orig, Ny_orig, Nz_orig])
                    [medium, p0] = pad_forward(medium, p0, Nx_orig, Ny_orig, Nz_orig, ...
                                               Nx_pad2, Ny_pad2, Nz_pad2);
                    if blind_recon
                        medium_recon = pad_medium_to(medium_recon, [Nx_pad2, Ny_pad2, Nz_pad2]);
                    end
                end
                if ~blind_recon
                    medium_recon = medium;
                end

                Nx = Nx_pad2; Ny = Ny_pad2; Nz = Nz_pad2;
                gridSize = [Nx, Ny, Nz];
                sensor.mask = zeros(Nx, Ny, Nz);
            end

            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
        otherwise
            error('simulate_field_once:UnknownSensor', ...
                'This study supports determine_sensor_mask / full_plane_* only (got "%s").', ...
                CONFIG.sensor_placement_method);
    end

    if sum(sensor.mask(:)) == 0
        error('Sensor mask is empty.');
    end

    %% ---- k-WAVE GRID & MEDIUM SETUP ----
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    maxC = max(medium.sound_speed(:));
    minC = min(medium.sound_speed(medium.sound_speed > 0));
    dt   = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;

    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);

    if isfield(CONFIG, 'Nt_scaling') && CONFIG.Nt_scaling ~= 0
        air_c = 343;
        if isfield(CONFIG, 'tissue_tables') && isfield(CONFIG.tissue_tables, 'threshold_2') ...
                && isfield(CONFIG.tissue_tables.threshold_2, 'air')
            air_c = CONFIG.tissue_tables.threshold_2.air.sound_speed;
        end
        if minC <= air_c
            Nt = max(1, ceil(Nt / CONFIG.Nt_scaling));
        end
    end

    kgrid.dt = dt;
    kgrid.Nt = Nt;

    kmedium_fwd             = struct();
    kmedium_fwd.density     = medium.density;
    kmedium_fwd.sound_speed = medium.sound_speed;
    kmedium_fwd.alpha_coeff = medium.alpha_coeff;
    kmedium_fwd.alpha_power = 1.1;

    kmedium_recon             = struct();
    kmedium_recon.density     = medium_recon.density;
    kmedium_recon.sound_speed = medium_recon.sound_speed;
    kmedium_recon.alpha_coeff = medium_recon.alpha_coeff;
    kmedium_recon.alpha_power = 1.1;

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

    %% ---- FORWARD SIMULATION ----
    source_fwd    = struct();
    source_fwd.p0 = p0;

    fwd_tic    = tic;
    sensorData = kspaceFirstOrder3D(kgrid, kmedium_fwd, source_fwd, sensor, inputArgs{:});
    fwd_time   = toc(fwd_tic);

    sensorData          = smooth(sensorData);
    FS                  = 1 / kgrid.dt;
    sensorData_measured = sensorData;

    %% ---- PULSE / FREQUENCY RESPONSE / NOISE / DECONV ----
    if CONFIG.convolution_kernel > 0
        conv_kernel_sigma  = CONFIG.convolution_kernel;
        conv_noise_level   = CONFIG.conv_noise_level;
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

        % Deterministic electronic-noise draw so the stability comparison
        % isolates the swept parameter (not run-to-run RNG variance).
        if isfield(CONFIG, 'rng_seed') && ~isempty(CONFIG.rng_seed)
            rng(CONFIG.rng_seed);
        end
        noise_amp        = conv_noise_level * max(abs(sensorData_resp(:)));
        sensorData_noisy = sensorData_resp + noise_amp * randn(size(sensorData_resp));

        sensorData_deconv = real(ifft( ...
            fft(sensorData_noisy, [], 2) .* H_conj ./ (H_power + conv_deconv_lambda), [], 2));

        sensorData          = single(sensorData_deconv);
        sensorData_measured = single(sensorData_deconv);
    else
        sensorData          = gaussianFilter(sensorData, FS, 0.35e6, 100, true);
        sensorData_measured = sensorData;
    end

    %% ---- ITERATIVE TIME REVERSAL ----
    reconPressure      = zeros(gridSize);
    reconPressure_prev = zeros(gridSize);
    num_iters_done     = 0;

    tr_total_tic = tic;
    for tr_iter = 1:CONFIG.num_time_reversal_iter
        source_tr        = struct();
        source_tr.p_mask = sensor.mask;
        source_tr.p      = fliplr(sensorData);
        source_tr.p_mode = 'dirichlet';

        sensor_tr        = struct();
        sensor_tr.mask   = ones(Nx, Ny, Nz);
        sensor_tr.record = {'p_final'};

        p0_recon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_tr, sensor_tr, inputArgs{:});
        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
        else
            reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
        end
        reconPressure  = max(reconPressure, 0);
        num_iters_done = tr_iter;

        converged = false;
        if tr_iter > 1
            norm_prev = norm(reconPressure_prev(:));
            if norm_prev > 0
                rel_change = norm(reconPressure(:) - reconPressure_prev(:)) / norm_prev;
            else
                rel_change = Inf;
            end
            if rel_change < CONFIG.convergence_tol
                converged = true;
            end
        end
        reconPressure_prev = reconPressure;

        if converged, break; end

        if tr_iter < CONFIG.num_time_reversal_iter
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_resid, sensor, inputArgs{:});
            sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
        end
    end
    reconPressure = gather(reconPressure) * CONFIG.correction_factor;
    tr_time = toc(tr_total_tic);

    %% ---- CROP TO ORIGINAL (padded -> gridSize_orig) ----
    if did_pad || ~isequal(size(reconPressure), gridSize_orig)
        reconPressure = reconPressure(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
        Nx = Nx_orig; Ny = Ny_orig; Nz = Nz_orig;
        gridSize = gridSize_orig;
        medium   = medium_orig;
        sensor.mask = sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    end

    if blind_recon
        medium = crop_medium(medium_recon, gridSize);
    end

    %% ---- PRESSURE SCALE CORRECTION ----
    if CONFIG.use_pressure_scale_correction
        p0_max_orig = max(p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig), [], 'all');
        recon_max   = max(reconPressure(:));
        if recon_max > 0
            reconPressure = reconPressure * (p0_max_orig / recon_max);
        end
    end

    %% ---- PRESSURE -> DOSE ----
    conversionFactor = medium.gruneisen .* medium.density;
    conversionFactor(conversionFactor == 0) = 1;
    reconDosePerPulse = reconPressure ./ conversionFactor;

    body_mask_plot = ones(gridSize);
    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), gridSize)
        body_mask_plot = double(sct.bodyMask);
    end
    if isfield(CONFIG, 'mask_recon_to_dose_region') && ~CONFIG.mask_recon_to_dose_region
        recon_dose = reconDosePerPulse * num_pulses .* body_mask_plot;
    else
        recon_dose = reconDosePerPulse * num_pulses .* double(doseMask) .* body_mask_plot;
    end
    recon_dose(medium.density < 100) = 0;

    %% ---- EXTRACT ANATOMY CORE & RESAMPLE BACK TO NATIVE DIMS ----
    % recon_dose is at gridSize_orig (downscaled-native, possibly expanded).
    % Crop out the anatomy core (dsDims), then resample to native dims so all
    % gamma analyses run on the original resolution regardless of the sweep.
    recon_core = recon_dose(core_xr, core_yr, core_zr);
    if ~isequal(size(recon_core), nativeDims)
        recon_native = max(imresize3(recon_core, nativeDims), 0);
    else
        recon_native = recon_core;
    end

    out.recon_native   = recon_native;
    out.dose_native    = dose_native;
    out.spacing_native = spacing_native;
    out.fwd_time       = fwd_time;
    out.tr_time        = tr_time;
    out.num_iters_done = num_iters_done;
end


function [medium, p0] = pad_forward(medium, p0, Nx, Ny, Nz, Nxp, Nyp, Nzp)
%PAD_FORWARD  Water-pad the forward medium + p0 into an [Nxp Nyp Nzp] grid.
    density_pad    = ones(Nxp, Nyp, Nzp) * 1000;
    soundSpeed_pad = ones(Nxp, Nyp, Nzp) * 1540;
    alphaCoeff_pad = zeros(Nxp, Nyp, Nzp);
    gruneisen_pad  = zeros(Nxp, Nyp, Nzp);
    p0_pad         = zeros(Nxp, Nyp, Nzp);

    density_pad(1:Nx, 1:Ny, 1:Nz)    = medium.density;
    soundSpeed_pad(1:Nx, 1:Ny, 1:Nz) = medium.sound_speed;
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_pad(1:Nx, 1:Ny, 1:Nz) = medium.alpha_coeff;
    else
        alphaCoeff_pad(:) = medium.alpha_coeff;
    end
    gruneisen_pad(1:Nx, 1:Ny, 1:Nz) = medium.gruneisen;
    p0_pad(1:Nx, 1:Ny, 1:Nz)        = p0;

    medium.density     = density_pad;
    medium.sound_speed = soundSpeed_pad;
    medium.alpha_coeff = alphaCoeff_pad;
    medium.gruneisen   = gruneisen_pad;
    p0 = p0_pad;
end


function medium = apply_standalone_medium_overrides(medium, sct, config)
%APPLY_STANDALONE_MEDIUM_OVERRIDES  Force-uniform toggles + coupling bath.
    gs = medium.grid_size;
    if config.force_uniform_density
        medium.density = ones(gs) * config.uniform_density;
    end
    if config.force_uniform_sound_speed
        medium.sound_speed = ones(gs) * config.uniform_sound_speed;
    end
    if config.force_uniform_attenuation
        medium.alpha_coeff = ones(gs) * config.uniform_alpha_coeff;
        medium.alpha_power = 1.1;
    end
    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gs) * config.uniform_gruneisen;
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
end


function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES  Tissue property lookup tables for HU thresholding.
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
    tables.threshold_2.air_hu_threshold = -300;
    tables.threshold_2.air = struct( ...
        'density',     1.2, ...
        'sound_speed', 343, ...
        'alpha_coeff', 0, ...
        'alpha_power', 1.0, ...
        'gruneisen',   0);
end


function m = pad_medium_to(m, sz)
%PAD_MEDIUM_TO  Water-pad the four medium fields into an sz grid.
    d = ones(sz) * 1000;
    d(1:size(m.density,1),     1:size(m.density,2),     1:size(m.density,3))     = m.density;
    c = ones(sz) * 1540;
    c(1:size(m.sound_speed,1), 1:size(m.sound_speed,2), 1:size(m.sound_speed,3)) = m.sound_speed;
    a = zeros(sz);
    if numel(m.alpha_coeff) > 1
        a(1:size(m.alpha_coeff,1), 1:size(m.alpha_coeff,2), 1:size(m.alpha_coeff,3)) = m.alpha_coeff;
    else
        a(:) = m.alpha_coeff;
    end
    g = zeros(sz);
    g(1:size(m.gruneisen,1),   1:size(m.gruneisen,2),   1:size(m.gruneisen,3))   = m.gruneisen;
    m.density = d; m.sound_speed = c; m.alpha_coeff = a; m.gruneisen = g;
end


function m = expand_medium(m_unp, sz, xr, yr, zr)
%EXPAND_MEDIUM  Place cropped medium fields into a water-filled sz grid.
    d = ones(sz) * 1000;  d(xr, yr, zr) = m_unp.density;
    c = ones(sz) * 1540;  c(xr, yr, zr) = m_unp.sound_speed;
    a = zeros(sz);
    if numel(m_unp.alpha_coeff) > 1
        a(xr, yr, zr) = m_unp.alpha_coeff;
    else
        a(:) = m_unp.alpha_coeff;
    end
    g = zeros(sz);        g(xr, yr, zr) = m_unp.gruneisen;
    m = struct('density', d, 'sound_speed', c, 'alpha_coeff', a, 'gruneisen', g);
end


function m = crop_medium(m, sz)
%CROP_MEDIUM  Crop the four medium fields to the leading sz sub-block.
    m.density     = m.density(1:sz(1),     1:sz(2),     1:sz(3));
    m.sound_speed = m.sound_speed(1:sz(1), 1:sz(2),     1:sz(3));
    if numel(m.alpha_coeff) > 1
        m.alpha_coeff = m.alpha_coeff(1:sz(1), 1:sz(2), 1:sz(3));
    end
    m.gruneisen   = m.gruneisen(1:sz(1),   1:sz(2),     1:sz(3));
end


function write_report(report_path, CONFIG, SWEEP, PATHS, R_default, RESULTS, SWEEPS)
%WRITE_REPORT  Human-readable .txt sweep report.
    crit    = SWEEP.gamma_criteria;
    nCrit   = size(crit, 1);
    primary = SWEEP.primary_crit_idx;
    plabel  = crit{primary, 3};

    fid = fopen(report_path, 'w');
    if fid < 0
        warning('Could not open report file %s; printing to console.', report_path);
        fid = 1;
    end
    cl = onCleanup(@() fid > 2 && fclose(fid));

    fprintf(fid, '================================================================\n');
    fprintf(fid, '  ETHOS OPTIMIZATION SWEEP REPORT\n');
    fprintf(fid, '  Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '================================================================\n');
    fprintf(fid, '  Patient / session : %s / %s\n', CONFIG.patient_id, CONFIG.session);
    fprintf(fid, '  Dose file         : %s\n', CONFIG.dose_filename);
    fprintf(fid, '  Reference CT      : CT_%d   (CT pair [%s])\n', CONFIG.reference_ct, num2str(CONFIG.ct_pair));
    fprintf(fid, '  Reconstruction    : %s, sensor=%s, tissue=%s\n', ...
        CONFIG.reconstruction_method, CONFIG.sensor_placement_method, CONFIG.gruneisen_method);
    fprintf(fid, '  Native dims       : [%s]\n', num2str(R_default.dims_native));
    if isfield(CONFIG, 'normalize') && CONFIG.normalize
        fprintf(fid, '  Normalization     : ON  (least-squares gain, recon -> own-CT truth 10%% region)\n');
    else
        fprintf(fid, '  Normalization     : OFF (absolute Gy, correction_factor scale)\n');
    end
    fprintf(fid, '  RNG seed (noise)  : %d  (deterministic noise draw)\n', CONFIG.rng_seed);
    fprintf(fid, '\n');
    fprintf(fid, '  METHOD: one-factor-at-a-time. Each run changes ONE variable and\n');
    fprintf(fid, '  holds all others at their default. The all-default run below is the\n');
    fprintf(fid, '  reference. recon = recon1 (reference-CT, full-access reconstruction).\n');
    fprintf(fid, '  Runtime = sum of forward + time-reversal wall time over BOTH doses.\n');
    fprintf(fid, '\n');
    fprintf(fid, '  SUCCESS if EITHER (at %s):\n', plabel);
    fprintf(fid, '    (1) gamma(default recon vs swept recon) >= %g%%  AND  runtime <= default\n', ...
        CONFIG.success_stability_pct);
    fprintf(fid, '    (2) gamma(swept recon vs truth) > default''s  AND  runtime <= default\n');
    fprintf(fid, '    (runtime tolerance for "equal or lower": %.0f%%)\n', 100 * CONFIG.runtime_tol);
    fprintf(fid, '\n');

    % --- default baseline ---
    fprintf(fid, '----------------------------------------------------------------\n');
    fprintf(fid, '  DEFAULT / REFERENCE RUN\n');
    fprintf(fid, '----------------------------------------------------------------\n');
    fprintf(fid, '  downscale_factor=%g  pml_size=%d  cfl=%.2f  tr_iter=%d  Nt_scaling=%g  conv_tol=%.1e\n', ...
        CONFIG.downscale_factor, CONFIG.pml_size, CONFIG.cfl_number, ...
        CONFIG.num_time_reversal_iter, CONFIG.Nt_scaling, CONFIG.convergence_tol);
    fprintf(fid, '  runtime = %.1f s   (fwd %s | tr %s ; iters %s)\n', ...
        R_default.runtime_sec, mat2str(round(R_default.fwd_time,1)), ...
        mat2str(round(R_default.tr_time,1)), mat2str(R_default.iters));
    print_pass_block(fid, R_default.pass, crit);
    fprintf(fid, '\n');

    % --- each sweep ---
    for s = 1:numel(SWEEPS)
        sw = SWEEPS{s};
        mask = strcmp({RESULTS.sweep}, sw.name);
        rs = RESULTS(mask);

        fprintf(fid, '================================================================\n');
        fprintf(fid, '  SWEEP: %s   (field: %s)\n', sw.name, sw.field);
        fprintf(fid, '================================================================\n');

        % Compact summary table for this sweep.
        fprintf(fid, '  %-10s %8s %8s %9s %10s %10s %6s %6s %8s\n', ...
            'value', 'runtime', 'dRun%', 'iters', ['stab-' plabel], ['R1vT-' plabel], 'c1', 'c2', 'RESULT');
        fprintf(fid, '  %s\n', repmat('-', 1, 86));
        for i = 1:numel(rs)
            r = rs(i);
            if ~r.ok
                fprintf(fid, '  %-10s %8s   %-60s\n', r.label, 'FAIL', ['ERROR: ' r.error]);
                continue;
            end
            dRun = 100 * (r.runtime_sec - R_default.runtime_sec) / max(R_default.runtime_sec, eps);
            stab_p = NaN; if ~isempty(r.stability), stab_p = r.stability(primary); end
            r1vt_p = r.pass.recon1_vs_dose1(primary);
            res_str = 'no';
            if r.success, res_str = 'SUCCESS'; end
            fprintf(fid, '  %-10s %7.1fs %+7.1f %9s %9.2f%% %9.2f%% %6s %6s %8s\n', ...
                r.label, r.runtime_sec, dRun, mat2str(r.iters), ...
                stab_p, r1vt_p, yn(r.crit1), yn(r.crit2), res_str);
        end
        fprintf(fid, '\n');

        % Full per-value gamma detail (all four comparisons, all criteria).
        for i = 1:numel(rs)
            r = rs(i);
            fprintf(fid, '  --- %s = %s ---\n', sw.field, r.label);
            if ~r.ok
                fprintf(fid, '      FAILED: %s\n\n', r.error);
                continue;
            end
            fprintf(fid, '      runtime %.1f s (default %.1f s), iters %s\n', ...
                r.runtime_sec, R_default.runtime_sec, mat2str(r.iters));
            if ~isempty(r.stability)
                fprintf(fid, '      stability (default recon vs this recon):%s\n', ...
                    fmt_crit(r.stability, crit));
            end
            print_pass_block(fid, r.pass, crit);
            fprintf(fid, '\n');
        end
    end

    % --- overall successes ---
    okmask = [RESULTS.ok];
    succ = RESULTS(okmask & [RESULTS.success]);
    fprintf(fid, '================================================================\n');
    fprintf(fid, '  SUCCESSFUL RUNS (%d)\n', numel(succ));
    fprintf(fid, '================================================================\n');
    if isempty(succ)
        fprintf(fid, '  (none met either success criterion)\n');
    else
        for i = 1:numel(succ)
            r = succ(i);
            why = '';
            if r.crit1, why = [why '[crit1: stable+faster] ']; end
            if r.crit2, why = [why '[crit2: more accurate+faster]']; end
            fprintf(fid, '  %-30s = %-10s  runtime %.1fs (default %.1fs)  %s\n', ...
                r.sweep, r.label, r.runtime_sec, R_default.runtime_sec, why);
        end
    end
    fprintf(fid, '\n(End of report)\n');
end


function print_pass_block(fid, pass, crit)
%PRINT_PASS_BLOCK  Print the four standalone-comparison gamma pass-rate rows.
    names = {'dose1_vs_dose2 (truth change)  ', ...
             'recon1_vs_dose1 (ref detector) ', ...
             'recon2_vs_dose2 (blind vs own) ', ...
             'recon2_vs_dose1 (blind vs ref) '};
    fields = {'dose1_vs_dose2', 'recon1_vs_dose1', 'recon2_vs_dose2', 'recon2_vs_dose1'};
    for k = 1:numel(fields)
        if isfield(pass, fields{k}) && ~isempty(pass.(fields{k}))
            fprintf(fid, '      %s:%s\n', names{k}, fmt_crit(pass.(fields{k}), crit));
        else
            fprintf(fid, '      %s: (n/a)\n', names{k});
        end
    end
end


function s = fmt_crit(vals, crit)
%FMT_CRIT  Format a pass-rate vector as "  10%/10mm=xx.x%  5%/5mm=..".
    s = '';
    for gc = 1:size(crit, 1)
        if gc <= numel(vals) && ~isnan(vals(gc))
            s = sprintf('%s   %s=%.2f%%', s, crit{gc, 3}, vals(gc));
        else
            s = sprintf('%s   %s=FAIL', s, crit{gc, 3});
        end
    end
end


function s = yn(tf)
    if tf, s = 'Y'; else, s = 'n'; end
end
