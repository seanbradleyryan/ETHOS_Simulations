%% =========================================================================
%  DIAGNOSTIC_RECON_COMPARISON.m
%  Reconstruction-setting comparison over one or more fields
%
%  Sweeps every combination of:
%    - dose file    CONFIG.dose_filenames (CT_1 field doses)
%    - method       CONFIG.recon_methods
%    - TV weight    CONFIG.tv_lambdas (adjoint_tv* methods only)
%    - iterations   CONFIG.num_iters  (iterative methods; 'das'/'ubp' run once)
%  One combination = one job. Each job goes through run_standalone_field
%  (load dose+CBCT -> build medium -> run_single_field_simulation). Every
%  setting on a dose file sees the SAME measured data (same dose, sensor and
%  noise draw) and runs exactly its iteration count (no early stop). Each
%  recon is normalized and scored against its truth with gamma (every
%  CONFIG.gamma_criteria row) and local SSIM, both over the truth's >=10%
%  region. Tables give each setting's score per dose file and the mean +/- std
%  over dose files; the best setting per method is printed.
%
%  CONFIG.eval_ct3_blind = true also runs each job's CT_3 partner field (same
%  name, _CT_1_ -> _CT_3_) blind: forward on CBCT3, reconstruction on CBCT1,
%  with the CT_1 run's sensor, as in Step 2. It then scores the Step 2.5
%  cross-CT comparisons: truth3 vs recon3, truth1 vs recon3 and truth1 vs
%  truth3. Recon differential = (truth1 vs recon1) - (truth1 vs recon3), as in
%  Step 3. It should stay near the truth differential, 100 - (truth1 vs
%  truth3). A setting that raises scores by smearing out the CT_1/CT_3
%  difference shrinks it.
%
%  CONFIG.use_parallel = true runs the jobs in a parfor, each worker pinned to
%  a GPU round-robin (as in pipeline_simulate). With adjoint_tv* methods the
%  first TV job of each dose file runs serially first, so the adjoint power
%  cache is on disk before the workers start.
%
%  Note: every job re-runs the forward simulation; rng(CONFIG.rng_seed) before
%  each run makes the noise identical. Every recon and truth1-vs-recon1 gamma
%  map is kept for the plots and the saved .mat, so memory grows with the
%  number of jobs.
%  =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'pipeline')));  % step* functions live here

%% ========================= CONFIGURATION ================================

CONFIG = get_default_config();

% --- Sweep settings: every dose file x method x lambda x iteration count ---
CONFIG.dose_filenames = {
    CONFIG.dose_filename
    % 'dose_1194203_Session_1_reference_CT_1_B15_112.mat'
    };   % CT_1 field doses in processed/; the CT_3 partners are derived
CONFIG.recon_methods  = {
    'tr',
    %'adjoint',
    'adjoint_tv',
    'adjoint_tv_tr'};   % any of 'tr','das','ubp','adjoint','adjoint_tv','adjoint_tv_tr'
CONFIG.num_iters      = 10;     % iteration counts to sweep, e.g. [1, 5, 10]
CONFIG.rng_seed       = 42;     % same electronic-noise draw for every method
CONFIG.gamma_criteria = {3, 3, '3%/3mm'};
CONFIG.eval_ct3_blind    = true;                   % also run the CT_3 partner field blind (doubles the k-Wave runs)
CONFIG.cbct_filename_ct3 = 'CBCT3_resampled.mat';  % CT_3 forward geometry; recon uses CONFIG.cbct_filename (CT_1)
CONFIG.Nt_scaling = 5;
CONFIG.conv_noise_level = .125;
CONFIG.conv_use_deconv = false;
CONFIG.gruneisen_method = 'continuous_1'

% --- Parallel: the jobs are independent, so they can run in a parfor ---
CONFIG.use_parallel         = false;   % true: parfor over the jobs
CONFIG.num_parallel_workers = 2;       % pool size; workers share the GPUs round-robin


% CONFIG.sensor_placement_method = 'box';
% CONFIG.gruneisen_method = 'uniform';
% 50, 75

%CONFIG.num_iter = 1;
%CONFIG.Nt_scaling = 3;
%33, 30.8, 30

% CONFIG.adjoint_tv_lambda = 0;

% Pass rate of 60 percent with lambda .1!!!! Note this



% Each job sets num_time_reversal_iter, adjoint_num_iter and
% adjoint_tv_inner_iter to its iteration count; tol 0 disables the early stop.
CONFIG.convergence_tol        = 0;

% --- adjoint_tv knobs (see run_single_field_simulation) ---
CONFIG.tv_lambdas            = [
    %0,
    %0.001,
    0.01,
    .03,
    0.1,
    %0.3,
    %0.5,
    % 1
    ];   % adjoint_tv_lambda sweep (adjoint_tv* methods only)
CONFIG.adjoint_tv_force_power_recalc = false;   % true = ignore processed/simulated_doses/adjoint_power_cache_<hash>.mat, recompute + overwrite
CONFIG.adjoint_tv_power_iter = 5;

% --- Example overrides ---
% CONFIG.sensor_placement_method = 'spherical';

% Gamma runs inside the jobs; fail now rather than after the k-Wave runs.
if exist('CalcGamma', 'file') ~= 2
    error('diagnostic_recon_comparison:NoCalcGamma', 'CalcGamma must be on the MATLAB path.');
end

% CT_3 partner of each CT_1 dose file: same beam/segment, other CT image.
numDoses     = numel(CONFIG.dose_filenames);
doseFilesCT3 = repmat({''}, numDoses, 1);
if CONFIG.eval_ct3_blind
    doseFilesCT3 = strrep(CONFIG.dose_filenames, '_CT_1_', '_CT_3_');
    noToken      = strcmp(doseFilesCT3, CONFIG.dose_filenames);
    if any(noToken)
        error('diagnostic_recon_comparison:NoCT3Name', ...
            'Dose file(s) with no _CT_1_ token to swap for _CT_3_: %s', ...
            strjoin(CONFIG.dose_filenames(noToken), ', '));
    end
end

norm_method = 'none';
if CONFIG.normalize
    norm_method = lower(CONFIG.normalize_method);
end

%% ========================= SETTINGS AND JOBS ============================
% A setting = method x lambda x iteration count. A job = one setting on one
% dose file. Jobs go dose file by dose file: job j = setting s on dose file d,
% with j = (d-1)*numSettings + s.

isTV        = startsWith(CONFIG.recon_methods, 'adjoint_tv');   % only these use adjoint_tv_lambda
isIterative = ~ismember(CONFIG.recon_methods, {'das', 'ubp'});  % only these use num_iters

SETTINGS = struct([]);
for mi = 1:numel(CONFIG.recon_methods)
    lambdas = NaN;   % NaN = unused: run once
    if isTV(mi)
        lambdas = CONFIG.tv_lambdas;
    end
    iters = NaN;     % NaN = unused: run once
    if isIterative(mi)
        iters = CONFIG.num_iters;
    end
    for li = 1:numel(lambdas)
        for ii = 1:numel(iters)
            s = numel(SETTINGS) + 1;
            SETTINGS(s).method    = CONFIG.recon_methods{mi};
            SETTINGS(s).tv_lambda = lambdas(li);
            SETTINGS(s).num_iter  = iters(ii);
            SETTINGS(s).label     = CONFIG.recon_methods{mi};
            if ~isnan(lambdas(li))
                SETTINGS(s).label = sprintf('%s, lambda=%g', SETTINGS(s).label, lambdas(li));
            end
            if ~isnan(iters(ii))
                SETTINGS(s).label = sprintf('%s, iter=%d', SETTINGS(s).label, iters(ii));
            end
        end
    end
end
numSettings = numel(SETTINGS);
numJobs     = numDoses * numSettings;

JOBS = struct([]);
for d = 1:numDoses
    for s = 1:numSettings
        j = (d - 1) * numSettings + s;
        JOBS(j).dose_idx          = d;
        JOBS(j).setting_idx       = s;
        JOBS(j).dose_filename     = CONFIG.dose_filenames{d};
        JOBS(j).dose_filename_ct3 = doseFilesCT3{d};
        JOBS(j).method            = SETTINGS(s).method;
        JOBS(j).tv_lambda         = SETTINGS(s).tv_lambda;
        JOBS(j).num_iter          = SETTINGS(s).num_iter;
        JOBS(j).label             = SETTINGS(s).label;
        JOBS(j).score_truth       = (s == 1);   % one job per dose file also keeps the truths + scores truth1 vs truth3
    end
end

fprintf('=========================================================\n');
fprintf('  Reconstruction setting comparison\n');
fprintf('=========================================================\n');
fprintf('  Patient:      %s / %s\n', CONFIG.patient_id, CONFIG.session);
for d = 1:numDoses
    fprintf('  Dose %d:       %s\n', d, CONFIG.dose_filenames{d});
    if CONFIG.eval_ct3_blind
        fprintf('    CT_3 blind: %s\n', doseFilesCT3{d});
    end
end
if CONFIG.eval_ct3_blind
    fprintf('  CT_3 blind:   forward %s, recon %s\n', CONFIG.cbct_filename_ct3, CONFIG.cbct_filename);
end
fprintf('  Sensor:       %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model: %s\n', CONFIG.gruneisen_method);
fprintf('  Methods:      %s\n', strjoin(CONFIG.recon_methods, ', '));
fprintf('  Iterations:   %s\n', mat2str(CONFIG.num_iters));
fprintf('  TV lambdas:   %s\n', mat2str(CONFIG.tv_lambdas));
fprintf('  Jobs:         %d settings x %d dose file(s) = %d (parallel: %d)\n', ...
    numSettings, numDoses, numJobs, CONFIG.use_parallel);
fprintf('=========================================================\n\n');

%% ========================= RUN THE JOBS =================================
% Serial mode: every job on this client, in order. Parallel mode: the
% power-cache warm-up jobs on this client, then the rest in a parfor.

serialJobs = 1:numJobs;
parJobs    = [];
if CONFIG.use_parallel && numJobs > 1
    % The first adjoint_tv* job of each dose file runs first, so its adjoint
    % power cache (run_single_field_simulation) is on disk before the workers
    % start. Otherwise every worker computes it at once and saves the same
    % file (same warm-up as pipeline_simulate).
    serialJobs = [];
    for d = 1:numDoses
        firstTV    = find([JOBS.dose_idx] == d & startsWith({JOBS.method}, 'adjoint_tv'), 1);
        serialJobs = [serialJobs, firstTV]; %#ok<AGROW>
    end
    parJobs = setdiff(1:numJobs, serialJobs);
end

JOB_OUT = cell(numJobs, 1);
for j = serialJobs
    JOB_OUT{j} = run_job(CONFIG, JOBS(j), norm_method);
end

if ~isempty(parJobs)
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= CONFIG.num_parallel_workers
        if ~isempty(pool), delete(pool); end
        parpool(CONFIG.num_parallel_workers);
    end
    % Pin each worker to a GPU round-robin (worker 1 -> GPU 1, worker 2 ->
    % GPU 2, worker 3 -> GPU 1, ...), as in pipeline_simulate.
    if CONFIG.use_gpu && gpuDeviceCount > 0
        numGpus = gpuDeviceCount;
        spmd
            gpu = gpuDevice(mod(labindex - 1, numGpus) + 1);
            fprintf('Worker %d -> GPU %d (%s)\n', labindex, gpu.Index, gpu.Name);
        end
    end

    parConfig = CONFIG;
    parConfig.adjoint_tv_force_power_recalc = false;   % the warm-up jobs already recomputed it
    parJobList = JOBS(parJobs);
    parOut     = cell(numel(parJobs), 1);
    parfor k = 1:numel(parJobs)
        parOut{k} = run_job(parConfig, parJobList(k), norm_method);
    end
    JOB_OUT(parJobs) = parOut;
end

% Every job returns the same fields, so the outputs form one struct array.
RES        = [JOB_OUT{:}];
truthJobs  = find([RES.score_truth]);   % RES(truthJobs(d)) holds dose file d's truths, sensor, density
spacing_mm = RES(1).spacing_mm;

%% ========================= SCORE TABLES =================================
% Score arrays are (setting, dose file, metric). Metrics: each gamma
% criterion's pass rate (%), then the mean local SSIM (%).

critLabels   = CONFIG.gamma_criteria(:, 3)';
metricLabels = [strcat({'gamma '}, critLabels), {'SSIM'}];
numMetrics   = numel(metricLabels);
scoreShape   = [numSettings, numDoses, numMetrics];

T1R1      = reshape([vertcat(RES.pass_t1r1), [RES.ssim_t1r1]'], scoreShape);
reconTime = mean(reshape([RES.tr_time], numSettings, numDoses), 2);   % mean over dose files
if CONFIG.eval_ct3_blind
    T3R3 = reshape([vertcat(RES.pass_t3r3), [RES.ssim_t3r3]'], scoreShape);
    T1R3 = reshape([vertcat(RES.pass_t1r3), [RES.ssim_t1r3]'], scoreShape);
    T1T3 = [vertcat(RES(truthJobs).pass_t1t3), [RES(truthJobs).ssim_t1t3]'];   % (dose file, metric)
end

fprintf('\n========= TRUTH1 vs RECON1 (normalization: %s) =========\n', norm_method);
for m = 1:numMetrics
    fprintf('\n  %s (%%)\n', metricLabels{m});
    fprintf('  %-36s %10s', 'Setting', 'Recon (s)');
    for d = 1:numDoses
        fprintf(' %10s', sprintf('Dose %d', d));
    end
    if numDoses > 1
        fprintf(' %16s', 'Mean +/- std');
    end
    fprintf('\n');
    for s = 1:numSettings
        fprintf('  %-36s %10.1f', SETTINGS(s).label, reconTime(s));
        fprintf(' %10.2f', T1R1(s, :, m));
        if numDoses > 1
            fprintf(' %16s', mean_std_str(T1R1(s, :, m)));
        end
        fprintf('\n');
    end
end

% CT_3 blind: one table per metric. Recon diff = T1vR1 - T1vR3 (Step 3's
% recon differential, paired per dose file); it should stay near the truth
% differential.
if CONFIG.eval_ct3_blind
    fprintf('\n========= CT_3 BLIND (mean +/- std over %d dose file(s)) =========\n', numDoses);
    for m = 1:numMetrics
        fprintf('\n  %s: truth1 vs truth3 %s, truth differential %s\n', metricLabels{m}, ...
            mean_std_str(T1T3(:, m)), mean_std_str(100 - T1T3(:, m)));
        fprintf('  %-36s %16s %16s %16s %16s\n', 'Setting', 'T1vR1', 'T3vR3', 'T1vR3', 'Recon diff');
        for s = 1:numSettings
            fprintf('  %-36s %16s %16s %16s %16s\n', SETTINGS(s).label, ...
                mean_std_str(T1R1(s, :, m)), mean_std_str(T3R3(s, :, m)), ...
                mean_std_str(T1R3(s, :, m)), mean_std_str(T1R1(s, :, m) - T1R3(s, :, m)));
        end
    end
end

% Best setting per method (over its lambdas and iteration counts), ranked by
% the first metric's mean over dose files.
fprintf('\n');
for mi = 1:numel(CONFIG.recon_methods)
    idx = find(strcmp({SETTINGS.method}, CONFIG.recon_methods{mi}));
    if numel(idx) < 2, continue; end
    [bestScore, b] = max(mean(T1R1(idx, :, 1), 2));
    fprintf('  Best %s setting: %s (%s mean %.2f%%)\n', CONFIG.recon_methods{mi}, ...
        SETTINGS(idx(b)).label, metricLabels{1}, bestScore);
end
fprintf('=========================================================\n');

%% ========================= SAVE RESULTS ================================

if CONFIG.save_results
    output_fname = sprintf('recon_comparison_%s_%s_%s.mat', CONFIG.sensor_placement_method, ...
        CONFIG.gruneisen_method, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    results = struct();
    results.settings      = SETTINGS;               % one per row of the score arrays
    results.dose_files    = CONFIG.dose_filenames;  % one per column
    results.metric_labels = metricLabels;           % one per page
    results.t1r1          = T1R1;
    results.res           = RES;                    % one per job: recons, gamma maps, timings
    results.truth_jobs    = truthJobs;              % res(truth_jobs(d)): dose file d's truths, sensor, truth1 vs truth3
    results.spacing_mm    = spacing_mm;
    results.config        = CONFIG;
    if CONFIG.eval_ct3_blind
        results.t3r3 = T3R3;
        results.t1r3 = T1R3;
        results.t1t3 = T1T3;   % (dose file, metric)
    end
    save(output_fname, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', output_fname);
end

%% ========================= PLOTS =======================================

if CONFIG.plot_results
    % Truth1 vs recon1 per setting (mean over dose files), one bar per metric; red line = 90%.
    settingLabels = {SETTINGS.label};
    figure('Name', 'Truth1 vs Recon1', 'Color', 'w', 'NumberTitle', 'off');
    bar(categorical(settingLabels, settingLabels), reshape(mean(T1R1, 2), numSettings, numMetrics));
    yline(90, 'r-', 'LineWidth', 1.5);
    ylim([0, 100]); grid on;
    ylabel('Gamma pass rate / SSIM (%)');
    legend([metricLabels, {'90%'}], 'Location', 'southoutside', 'Orientation', 'horizontal');
    title(sprintf('Truth1 vs recon1, mean over %d dose file(s)', numDoses));

    if CONFIG.eval_ct3_blind
        % First gamma criterion and SSIM. Dashed line = truth1 vs truth3: a
        % recon that keeps the CT_1/CT_3 difference has its T1vR3 bar near it.
        figure('Name', 'CT_3 Blind', 'Color', 'w', 'NumberTitle', 'off', 'Position', [100, 100, 1300, 450]);
        tiledlayout(1, 2);
        for m = [1, numMetrics]
            nexttile;
            bar(categorical(settingLabels, settingLabels), ...
                [mean(T1R1(:, :, m), 2), mean(T3R3(:, :, m), 2), mean(T1R3(:, :, m), 2)]);
            yline(mean(T1T3(:, m)), 'k--', 'LineWidth', 1.5);
            ylim([0, 100]); grid on;
            ylabel(sprintf('%s (%%), mean over %d dose file(s)', metricLabels{m}, numDoses));
            legend({'truth1 vs recon1', 'truth3 vs recon3 (blind)', 'truth1 vs recon3 (blind)', ...
                'truth1 vs truth3'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
            title(sprintf('%s: T1vR3 should stay near truth1 vs truth3', metricLabels{m}));
        end
    end

    % Dose panels and gamma maps (truth1 vs recon1): one window each, one tab per job.
    figDose  = figure('Name', 'Dose Panels', 'Color', 'w', 'NumberTitle', 'off', 'Position', [50, 50, 1380, 720]);
    tgDose   = uitabgroup(figDose);
    figGamma = figure('Name', 'Gamma Maps', 'Color', 'w', 'NumberTitle', 'off', 'Position', [50, 100, 1450, 400]);
    tgGamma  = uitabgroup(figGamma);
    for j = 1:numJobs
        T = RES(truthJobs(RES(j).dose_idx));   % this job's dose file: truth, sensor, density
        [~, maxIdx]     = max(T.truth1(:));
        [~, ~, czGamma] = ind2sub(size(T.truth1), maxIdx);
        tabTitle = sprintf('D%d | %s', RES(j).dose_idx, RES(j).label);
        plot_dose_panels(T.truth1, RES(j).recon_dose, T.sensor_mask, T.density, spacing_mm, ...
            sprintf('Truth vs %s', tabTitle), CONFIG.viz_smooth_sigma, ...
            {'Truth', RES(j).label}, uitab(tgDose, 'Title', tabTitle));
        plot_gamma_and_error_axial(RES(j).gamma, T.truth1, RES(j).recon_dose, T.sensor_mask, ...
            czGamma, tabTitle, uitab(tgGamma, 'Title', tabTitle));
    end
end

fprintf('\nReconstruction comparison complete.\n');


%% ========================= LOCAL HELPERS ================================

function R = run_job(CONFIG, job, norm_method)
%RUN_JOB One setting on one dose file: simulate the CT_1 field and, when
%  CONFIG.eval_ct3_blind, its CT_3 partner blind with the CT_1 sensor; score
%  each recon against its truth (gamma + SSIM). Self-contained so it can run
%  on a parfor worker. R has the same fields for every job (empty / NaN when
%  unused), so the job outputs concatenate into one struct array.
    cfg = CONFIG;
    cfg.dose_filename         = job.dose_filename;
    cfg.reconstruction_method = job.method;
    if ~isnan(job.tv_lambda)
        cfg.adjoint_tv_lambda = job.tv_lambda;
    end
    if ~isnan(job.num_iter)
        cfg.num_time_reversal_iter = job.num_iter;
        cfg.adjoint_num_iter       = job.num_iter;
        cfg.adjoint_tv_inner_iter  = job.num_iter;
    end

    R = job;
    R.tr_time        = 0;
    R.spacing_mm     = [];
    R.recon_dose     = [];
    R.gamma          = [];    % truth1 vs recon1, with maps (gamma-map plot)
    R.pass_t1r1      = [];
    R.ssim_t1r1      = NaN;
    R.recon_dose_ct3 = [];
    R.pass_t3r3      = [];
    R.ssim_t3r3      = NaN;
    R.pass_t1r3      = [];
    R.ssim_t1r3      = NaN;
    % score_truth job only: the dose file's truths, sensor, density, truth1 vs truth3
    R.truth1         = [];
    R.truth3         = [];
    R.sensor_mask    = [];
    R.density        = [];
    R.pass_t1t3      = [];
    R.ssim_t1t3      = NaN;

    fprintf('\n########## %s | %s ##########\n', job.dose_filename, job.label);
    rng(CONFIG.rng_seed);   % identical noise -> identical measured data
    out = run_standalone_field(cfg);
    R.tr_time        = getf(out.sim_results, 'tr_time_s', 0);
    R.spacing_mm     = out.spacing;

    truth1 = out.doseGrid;
    if strcmp(norm_method, 'max')
        truth1 = truth1 / max(truth1(:));
    end
    R.recon_dose = normalize_recon(out.recon_dose, truth1, norm_method);
    [R.gamma, R.ssim_t1r1] = score_pair(truth1, R.recon_dose, R.spacing_mm, ...
        CONFIG.gamma_criteria, 'truth1 vs recon1');
    R.pass_t1r1 = R.gamma.pass_rates(:)';
    if job.score_truth
        R.truth1      = truth1;
        R.sensor_mask = out.sensor_mask;
        R.density     = out.density;
    end

    if CONFIG.eval_ct3_blind
        fprintf('\n########## %s | %s, CT_3 blind ##########\n', job.dose_filename_ct3, job.label);
        cfg3                     = cfg;
        cfg3.dose_filename       = job.dose_filename_ct3;
        cfg3.cbct_filename       = CONFIG.cbct_filename_ct3;   % forward on CT_3
        cfg3.blind_recon         = true;
        cfg3.recon_cbct_filename = CONFIG.cbct_filename;       % reconstruct on CT_1
        cfg3.precomputed_sensor  = getf(out.sim_results, 'precomputed_sensor', []);   % CT_1 sensor ([] = geometric methods, placed identically)
        rng(CONFIG.rng_seed);
        out3 = run_standalone_field(cfg3);

        truth3 = out3.doseGrid;
        if strcmp(norm_method, 'max')
            truth3 = truth3 / max(truth3(:));
        end
        % Scaled to its own-CT truth (truth3), as in Step 2.5.
        R.recon_dose_ct3 = normalize_recon(out3.recon_dose, truth3, norm_method);
        [g, R.ssim_t3r3] = score_pair(truth3, R.recon_dose_ct3, R.spacing_mm, ...
            CONFIG.gamma_criteria, 'truth3 vs recon3 (blind)');
        R.pass_t3r3 = g.pass_rates(:)';
        [g, R.ssim_t1r3] = score_pair(truth1, R.recon_dose_ct3, R.spacing_mm, ...
            CONFIG.gamma_criteria, 'truth1 vs recon3 (blind)');
        R.pass_t1r3 = g.pass_rates(:)';
        if job.score_truth
            R.truth3 = truth3;
            [g, R.ssim_t1t3] = score_pair(truth1, truth3, R.spacing_mm, ...
                CONFIG.gamma_criteria, 'truth1 vs truth3');
            R.pass_t1t3 = g.pass_rates(:)';
        end
    end
end

function [g, ssimPct] = score_pair(ref, tgt, spacing_mm, criteria, name)
%SCORE_PAIR Gamma pass rates (every criterion) and mean local SSIM (%) of tgt
%  against ref, both over ref's >=10% region (compute_gamma's eval_mask).
    fprintf('\n[Gamma + SSIM] %s\n', name);
    g = compute_gamma(ref, tgt, spacing_mm, 'Criteria', criteria);
    [~, ssimMean] = compute_local_ssim(ref, tgt, g.eval_mask);
    ssimPct = 100 * ssimMean;
    fprintf('  [SSIM] %.2f%%\n', ssimPct);
end

function recon = normalize_recon(recon, truth, norm_method)
%NORMALIZE_RECON 'max' = divide by its own max; 'lsq' = least-squares gain to
%  the (unscaled) truth; 'none' = unchanged.
    switch norm_method
        case 'max'
            if max(recon(:)) > 0, recon = recon / max(recon(:)); end
        case 'lsq'
            recon = recon * least_squares_gain(truth, recon);
        case 'none'
        otherwise
            error('diagnostic_recon_comparison:BadNormalizeMethod', ...
                'Unknown CONFIG.normalize_method "%s" (use ''max'' or ''lsq'').', norm_method);
    end
end

function s = mean_std_str(v)
%MEAN_STD_STR '62.10' for one value, '62.10 +/- 3.40' (mean +/- std) for several.
    if numel(v) == 1
        s = sprintf('%.2f', v);
    else
        s = sprintf('%.2f +/- %.2f', mean(v), std(v));
    end
end

function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
