function results = step3_analysis(patient_id, session, config)
%% STEP3_ANALYSIS - Total-dose + per-segment analysis, figures, statistics
%
%   results = step3_analysis(patient_id, session, config)
%
%   PURPOSE:
%   Final analysis step. Runs after Step 2 (k-Wave recons) and Step 2.5
%   (per-segment gamma/SSIM folded into the CT_1 recon files). Two parts:
%
%   PART A - TOTAL-DOSE COMPARISONS (gamma + local SSIM, orthogonal-view figures)
%     1. ETHOS RTDOSE truth vs RayStation total for ONE CT label
%        (config.analysis_ethos_ct_label). Both are a single-plan dose.
%     2. RayStation total vs reconstructed total (both = sum of every field).
%
%   PART B - PER-SEGMENT STUDY (migrated from study_pass_rates_allsegments.m)
%     Reads the Step 2.5 summary (runs Step 2.5 if it is missing). For BOTH
%     metrics (gamma pass rate, mean local SSIM):
%       - per-beam mean +/- SE table for the four comparisons
%         (truth1_vs_truth3, truth1_vs_recon1, truth1_vs_recon3, truth3_vs_recon3)
%       - pass-rate, differential and own-CT fidelity plots by beam
%       - segment statistics: paired t-test recon1 vs recon3, proportion above
%         the noise floor, correlation of (recon1 - recon3) with the true change
%         and with the CT_1 SNR (+ scatter plot)
%     Plus per-beam electronic-noise SNR and N random segment panels
%     (truth | recon | gamma map | SSIM map).
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct (the pipeline CONFIG works as-is):
%           .working_dir              - Base directory (REQUIRED)
%           .treatment_site           - default 'Pancreas'
%           .gruneisen_method         - default 'threshold_2'
%           .config_hash              - 8-char sim hash; '' = auto (default '')
%           .gamma_dose_pct           - default 3.0
%           .gamma_dist_mm            - default 3.0
%           .gamma_dose_cutoff_pct    - eval mask, % of reference max (default 10)
%           .metrics_plan_type        - Step 2.5 plan type (default 'reference')
%           .metrics_ct_pair          - default [1, 3]
%           .metrics_normalize        - LS gain recon -> truth (default true)
%           .analysis_ethos_ct_label  - CT whose RS total is compared with ETHOS
%                                       (default 'CT_1')
%           .analysis_beams           - beams for Part B; [] = all (default [])
%           .analysis_noise_floor_pct - floor for the "above floor" test;
%                                       [] = Step 2.5 noise-floor mean (default [])
%           .analysis_plot_results    - save figures (default true)
%           .analysis_n_random        - random segment panels (default 5)
%           .analysis_random_seed     - default 42
%
%   OUTPUTS:
%       results - Struct:
%           .ethos_vs_rs / .rs_vs_recon - .gamma (pass_rate, mean_gamma,
%               max_gamma, num_evaluated, num_passed, gamma_map) and
%               .ssim (mean_pct, ssim_map)
%           .recon_gain       - LS gain applied to the total recon
%           .segments.gamma / .segments.ssim - per-beam table + segment stats
%           .snr              - per-beam SNR (mean, se, n)
%           .noise_floor      - Step 2.5 noise-only null (or [])
%           .random_segments  - beam/segment/scores of the random panels
%           .metadata
%
%   FILES CREATED (AnalysisResults/[PatientID]/[Session]/[method]/[hash]/):
%       step3_results.mat         - the results struct
%       beam_summary.csv          - per-beam means/SE, both metrics, + SNR
%       segment_statistics.csv    - t-test / floor / correlation rows
%       step3_console_log.txt     - full console output of this run
%       figures/*.png             - every figure
%
%   ALGORITHM:
%       1. Load total recon, RS total, ETHOS truth (position-resampled) and CBCT
%          via load_recon_dose_data; RS per-CT total from processed/.
%       2. Gamma (CalcGamma, global) + local SSIM over the 10% reference mask.
%       3. Load/run Step 2.5 summary; build per-beam tables and statistics.
%       4. Plot, save .mat/.csv, move the console log into the results folder.
%
%   EXAMPLE:
%       config.working_dir      = get_repo_root();
%       config.gruneisen_method = 'threshold_2';
%       config.config_hash      = 'a9a3e1e6';
%       results = step3_analysis('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       CalcGamma, compute_local_ssim, load_recon_dose_data,
%       step25_segment_metrics, Image Processing Toolbox (ssim).
%       Statistics are base MATLAB (no Statistics Toolbox).
%
%   See also: step25_segment_metrics, load_recon_dose_data, compute_local_ssim

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('step3_analysis:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

if ~ischar(session) && ~isstring(session)
    error('step3_analysis:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

if ~isstruct(config) || ~isfield(config, 'working_dir')
    error('step3_analysis:MissingConfig', 'config must be a struct containing working_dir.');
end

%% ======================== SET DEFAULTS ========================

if ~isfield(config, 'treatment_site'),           config.treatment_site = 'Pancreas'; end
if ~isfield(config, 'gruneisen_method'),         config.gruneisen_method = 'threshold_2'; end
if ~isfield(config, 'config_hash'),              config.config_hash = ''; end
if ~isfield(config, 'gamma_dose_pct'),           config.gamma_dose_pct = 3.0; end
if ~isfield(config, 'gamma_dist_mm'),            config.gamma_dist_mm = 3.0; end
if ~isfield(config, 'gamma_dose_cutoff_pct'),    config.gamma_dose_cutoff_pct = 10.0; end
if ~isfield(config, 'metrics_plan_type'),        config.metrics_plan_type = 'reference'; end
if ~isfield(config, 'metrics_ct_pair'),          config.metrics_ct_pair = [1, 3]; end
if ~isfield(config, 'metrics_normalize'),        config.metrics_normalize = true; end
if ~isfield(config, 'analysis_ethos_ct_label'),  config.analysis_ethos_ct_label = 'CT_1'; end
if ~isfield(config, 'analysis_beams'),           config.analysis_beams = []; end
if ~isfield(config, 'analysis_noise_floor_pct'), config.analysis_noise_floor_pct = []; end
if ~isfield(config, 'analysis_plot_results'),    config.analysis_plot_results = true; end
if ~isfield(config, 'analysis_n_random'),        config.analysis_n_random = 5; end
if ~isfield(config, 'analysis_random_seed'),     config.analysis_random_seed = 42; end

ct1_str = sprintf('CT_%d', min(config.metrics_ct_pair));
ct3_str = sprintf('CT_%d', max(config.metrics_ct_pair));
crit    = config.gamma_dose_pct;

%% ======================== CONSOLE LOG ========================
% The results folder depends on the resolved hash, so the log starts in the
% session folder and is moved into the results folder at the end.

session_dir = fullfile(config.working_dir, 'AnalysisResults', patient_id, session);
if ~isfolder(session_dir), mkdir(session_dir); end
tmp_log = fullfile(session_dir, 'step3_console_log_inprogress.txt');
if isfile(tmp_log), delete(tmp_log); end
diary off;
diary(tmp_log);
log_guard = onCleanup(@() diary('off'));   % stop logging even if an error is thrown

analysis_timer = tic;
fprintf('\n=========================================================\n');
fprintf('  [STEP 3] Total-Dose + Per-Segment Analysis\n');
fprintf('  Patient: %s | Session: %s | %s\n', patient_id, session, datetime('now'));
fprintf('  Gamma: %.1f%% / %.1f mm, eval mask >= %.0f%% of reference max\n', ...
    config.gamma_dose_pct, config.gamma_dist_mm, config.gamma_dose_cutoff_pct);
fprintf('=========================================================\n\n');

%% ======================== PART A: LOAD TOTAL DOSES ========================

fprintf('[STEP 3] Part A: total-dose comparisons\n');
load_args = {'Mode', 'total'};
if ~isempty(config.config_hash)
    load_args = [load_args, {'Hash', config.config_hash}];
end
T = load_recon_dose_data(patient_id, session, config, load_args{:});

hash8              = T.config_hash;
config.config_hash = hash8;   % pin for Step 2.5 / random-panel loads
sim_dir       = fullfile(config.working_dir, 'SimulationResults', patient_id, session, ...
    T.gruneisen_method);
processed_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session, 'processed');
out_dir       = fullfile(session_dir, T.gruneisen_method, hash8);
fig_dir       = fullfile(out_dir, 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end
fprintf('  Results folder: %s\n', out_dir);

spacing     = T.metadata.spacing(:)';
recon_total = double(gather(T.recon_dose));
rs_total    = double(T.rs_dose);
ethos_truth = double(T.ethos_truth);
rs_ct_total = load_ct_total_dose(processed_dir, config.analysis_ethos_ct_label);
dims        = size(rs_total);

if ~isequal(size(recon_total), dims) || ~isequal(size(rs_ct_total), dims) ...
        || ~isequal(size(ethos_truth), dims)
    error('step3_analysis:SizeMismatch', ...
        'Grid mismatch: RS %s, recon %s, RS %s total %s, ETHOS %s.', mat2str(dims), ...
        mat2str(size(recon_total)), config.analysis_ethos_ct_label, ...
        mat2str(size(rs_ct_total)), mat2str(size(ethos_truth)));
end
fprintf('  Grid [%d x %d x %d], spacing [%.2f %.2f %.2f] mm\n', dims, spacing);
fprintf('  Max dose: ETHOS %.4f | RS %s %.4f | RS all %.4f | recon %.4f Gy\n', ...
    max(ethos_truth(:)), config.analysis_ethos_ct_label, max(rs_ct_total(:)), ...
    max(rs_total(:)), max(recon_total(:)));

% Same least-squares scaling Step 2.5 applies to each segment recon.
recon_gain = 1;
if config.metrics_normalize
    recon_gain  = least_squares_gain(rs_total, recon_total, config.gamma_dose_cutoff_pct / 100);
    recon_total = recon_total * recon_gain;
    fprintf('  Recon LS gain (recon -> RS total): %.4g\n', recon_gain);
end

% Overlays: body contour from the CT_1 CBCT, sensor from Step 2's saved mask.
body_mask = [];
if isfield(T, 'cbct') && isfield(T.cbct, ct1_str) && isfield(T.cbct.(ct1_str), 'bodyMask') ...
        && isequal(size(T.cbct.(ct1_str).bodyMask), dims)
    body_mask = logical(T.cbct.(ct1_str).bodyMask);
else
    fprintf('  [NOTE] No %s body mask on the dose grid; contour overlay skipped.\n', ct1_str);
end
sensor_mask = load_sensor_mask(sim_dir, hash8, dims);
clear T;

%% ======================== PART A: GAMMA + SSIM ========================

fprintf('\n  ETHOS truth vs RayStation %s total...\n', config.analysis_ethos_ct_label);
results.ethos_vs_rs = compare_total_doses(ethos_truth, rs_ct_total, spacing, config);

fprintf('\n  RayStation total vs reconstructed total...\n');
results.rs_vs_recon = compare_total_doses(rs_total, recon_total, spacing, config);
results.recon_gain  = recon_gain;

if config.analysis_plot_results
    fprintf('\n  Saving total-dose figures...\n');
    ethos_rs_label = sprintf('RS %s', strrep(config.analysis_ethos_ct_label, '_', '\_'));
    plot_total_comparison(ethos_truth, rs_ct_total, results.ethos_vs_rs, ...
        'ETHOS truth', ethos_rs_label, 'ethos_vs_rs', spacing, body_mask, sensor_mask, fig_dir);
    plot_total_comparison(rs_total, recon_total, results.rs_vs_recon, ...
        'RS total', 'Recon total', 'rs_vs_recon', spacing, body_mask, sensor_mask, fig_dir);
end
clear ethos_truth rs_ct_total rs_total recon_total body_mask sensor_mask;

%% ======================== PART B: STEP 2.5 SUMMARY ========================

fprintf('\n[STEP 3] Part B: per-segment study (Step 2.5 results)\n');
SM = load_step25_summary(patient_id, session, config, sim_dir, hash8);

if ~strcmpi(char(SM.plan_type), char(config.metrics_plan_type))
    warning('step3_analysis:PlanTypeMismatch', ...
        'config.metrics_plan_type=%s but the Step 2.5 summary is %s; using the summary.', ...
        config.metrics_plan_type, char(SM.plan_type));
    config.metrics_plan_type = char(SM.plan_type);
end
if isfield(SM, 'gamma_dose_pct') && ~isequal(SM.gamma_dose_pct, config.gamma_dose_pct)
    warning('step3_analysis:CritMismatch', ...
        'Step 2.5 used %g%%/%g mm; per-segment gamma is labelled with that.', ...
        SM.gamma_dose_pct, SM.gamma_dist_mm);
    crit = SM.gamma_dose_pct;
end

c_tt = find(strcmp(SM.comparisons, 'truth1_vs_truth3'));
c_r1 = find(strcmp(SM.comparisons, 'truth1_vs_recon1'));
c_r3 = find(strcmp(SM.comparisons, 'truth1_vs_recon3'));
c_33 = find(strcmp(SM.comparisons, 'truth3_vs_recon3'));
if isempty(c_tt) || isempty(c_r1) || isempty(c_r3) || isempty(c_33)
    error('step3_analysis:MissingComparison', ...
        'Step 2.5 summary must contain the four study comparisons.');
end

% Beam selection (sorted), restricted to config.analysis_beams when given.
if isempty(config.analysis_beams)
    want = SM.beams(:)';
else
    want = config.analysis_beams(:)';
end
[found, loc] = ismember(want, SM.beams(:)');
if ~all(found)
    error('step3_analysis:BeamsNotInSummary', ...
        'Beam(s) %s are not in the Step 2.5 summary.', mat2str(want(~found)));
end
[beam_list, ord] = sort(SM.beams(loc));
beam_list = beam_list(:)';   % row, so "for b = beam_list" visits each beam
loc       = loc(ord);
n_seg     = SM.n_segments(loc);
nB        = numel(beam_list);
x_labels  = [arrayfun(@(b) sprintf('%d', b), beam_list, 'UniformOutput', false), {'All'}];

noise_floor = [];
if isfield(SM, 'noise_floor') && ~isempty(SM.noise_floor)
    noise_floor = SM.noise_floor;
end
floor_pct = config.analysis_noise_floor_pct;
if isempty(floor_pct)
    if ~isempty(noise_floor)
        floor_pct = noise_floor.mean_pass_rate;
    else
        floor_pct = 17;
        warning('step3_analysis:NoNoiseFloor', ...
            'No Step 2.5 noise floor; using a placeholder floor of %g%%.', floor_pct);
    end
end

fprintf('  Plan %s | hash %s | beams %s (%d segments)\n', config.metrics_plan_type, ...
    hash8, mat2str(beam_list), sum(n_seg));
if ~isempty(noise_floor)
    fprintf('  Noise-only floor: %.2f +/- %.2f %% (n=%d)\n', noise_floor.mean_pass_rate, ...
        noise_floor.std_pass_rate, noise_floor.num_samples);
end

%% ======================== PART B: SNR ========================

snr = gather_snr(sim_dir, hash8, config.metrics_plan_type, ct1_str, beam_list, n_seg);
fprintf('\n==================== BEAM SNR (mean +/- SE over fields) ====================\n');
for n = 1:nB + 1
    fprintf('  beam %-4s  SNR %7.2f +/- %6.2f  (%d field(s))\n', x_labels{n}, ...
        snr.mean(n), snr.se(n), snr.n(n));
end
results.snr = rmfield(snr, 'seg_ct1');
if config.analysis_plot_results
    plot_beam_snr(x_labels, snr, patient_id, session, fullfile(fig_dir, 'beam_snr.png'));
end

%% ======================== PART B: BOTH METRICS ========================

metric_keys = {'gamma', 'ssim'};
for m = 1:2
    key = metric_keys{m};
    if strcmp(key, 'gamma')
        seg_by_beam = SM.gamma.seg_pass(loc);
        mlabel = sprintf('Gamma pass rate (%%) @ %g%%/%g mm', crit, crit);
        mword  = 'Gamma pass rate';
        nf     = noise_floor;      % the noise-only null is a gamma pass rate
        fl     = floor_pct;
    else
        seg_by_beam = SM.ssim.seg_mean(loc);
        mlabel = 'Mean local SSIM (%)';
        mword  = 'Local SSIM';
        nf     = [];
        fl     = NaN;              % no SSIM floor exists
    end

    B = beam_stats(seg_by_beam);   % rows 1..nB = beams, row nB+1 = All
    fprintf('\n==================== BEAM %s (mean +/- SE over segments) ====================\n', ...
        upper(mword));
    for n = 1:nB + 1
        fprintf('\n----- [beam %s]  (%d segments) -----\n', x_labels{n}, B.n(n));
        for c = 1:numel(SM.comparisons)
            fprintf('  %-18s   %6.2f%% +/- %5.2f%%\n', SM.comparisons{c}, B.mean(n, c), B.se(n, c));
        end
    end

    D     = differential_stats(seg_by_beam, c_r1, c_r3, c_tt);
    stats = segment_statistics(seg_by_beam, c_r1, c_r3, c_tt, snr.seg_ct1, fl, x_labels);

    fprintf('\n==================== %s SEGMENT STATISTICS ====================\n', upper(mword));
    fprintf('A = recon_CT1 vs truth_CT1 | B = recon_CT3 vs truth_CT1 | diff = A - B\n');
    for n = 1:nB + 1
        s = stats(n);
        fprintf('\n----- [beam %s] -----\n', x_labels{n});
        fprintf('  (1) Paired t: A > B in %5.1f%% of %d seg | t = %7.3f, one-sided p = %.3g\n', ...
            s.pct_a_gt_b, s.n_pairs, s.t_stat, s.p_one_sided);
        if isfinite(fl)
            fprintf('  (2) Above %.1f%% floor: A %5.1f%% +/- %4.1f%% | B %5.1f%% +/- %4.1f%%\n', ...
                fl, s.above_floor_a_pct, s.above_floor_a_se, s.above_floor_b_pct, s.above_floor_b_se);
        end
        fprintf('  (3) diff vs truth: R^2 = %.3f, rho = %+.3f | diff vs SNR: R^2 = %.3f, rho = %+.3f\n', ...
            s.r2_vs_truth, s.rho_vs_truth, s.r2_vs_snr, s.rho_vs_snr);
    end

    if config.analysis_plot_results
        title_tag = sprintf('%s / %s', strrep(patient_id, '_', '\_'), strrep(session, '_', '\_'));
        plot_pass_rates(x_labels, B, c_r1, c_r3, c_tt, nf, mlabel, mword, title_tag, ...
            fullfile(fig_dir, sprintf('%s_pass_rates_by_beam.png', key)));
        plot_differentials(x_labels, D, B.mean(:, c_r1), nf, mword, title_tag, ...
            fullfile(fig_dir, sprintf('%s_differentials_by_beam.png', key)));
        plot_fidelity(x_labels, B, c_r1, c_33, mlabel, title_tag, ...
            fullfile(fig_dir, sprintf('%s_recon_fidelity_by_beam.png', key)));
        plot_correlations(stats(end).pooled_diff, stats(end).pooled_truth, ...
            stats(end).pooled_snr, mword, title_tag, ...
            fullfile(fig_dir, sprintf('%s_diff_correlations.png', key)));
    end

    results.segments.(key).comparisons  = SM.comparisons;
    results.segments.(key).beam_labels  = x_labels;
    results.segments.(key).beam_stats   = B;
    results.segments.(key).differential = D;
    results.segments.(key).statistics   = rmfield(stats, {'pooled_diff', 'pooled_truth', 'pooled_snr'});
    results.segments.(key).seg_by_beam  = seg_by_beam;
    results.segments.(key).floor_pct    = fl;
end
results.noise_floor = noise_floor;

%% ======================== PART B: RANDOM SEGMENT PANELS ========================

results.random_segments = [];
if config.analysis_plot_results && config.analysis_n_random > 0
    fprintf('\n----------- RANDOM SEGMENT PANELS -----------\n');
    picks = select_random_segments(patient_id, session, config, beam_list, ct1_str, ct3_str);
    for i = 1:numel(picks)
        plot_random_segment(picks(i), crit, fullfile(fig_dir, ...
            sprintf('random_segment_B%d_S%d.png', picks(i).beam, picks(i).seg)));
    end
    results.random_segments = rmfield(picks, 'rows');
end

%% ======================== SAVE RESULTS ========================

fprintf('\n  Saving results...\n');
results.metadata.patient_id   = patient_id;
results.metadata.session      = session;
results.metadata.config_hash  = hash8;
results.metadata.timestamp    = datetime('now');
results.metadata.config       = config;
results.metadata.spacing_mm   = spacing;
results.metadata.grid_size    = dims;
results.metadata.output_dir   = out_dir;

save(fullfile(out_dir, 'step3_results.mat'), 'results', '-v7.3');
fprintf('    Saved: step3_results.mat\n');
write_beam_csv(fullfile(out_dir, 'beam_summary.csv'), results, beam_list);
write_stats_csv(fullfile(out_dir, 'segment_statistics.csv'), results);

%% ======================== SUMMARY ========================

fprintf('\n=========================================================\n');
fprintf('  [STEP 3] Analysis Complete (%.1f sec)\n', toc(analysis_timer));
fprintf('  ETHOS vs RS %s: gamma %.1f%% | local SSIM %.1f%%\n', config.analysis_ethos_ct_label, ...
    results.ethos_vs_rs.gamma.pass_rate, results.ethos_vs_rs.ssim.mean_pct);
fprintf('  RS vs Recon total: gamma %.1f%% | local SSIM %.1f%%\n', ...
    results.rs_vs_recon.gamma.pass_rate, results.rs_vs_recon.ssim.mean_pct);
fprintf('  All-segment gamma: recon1 %.1f%% | recon3 %.1f%% | paired p = %.3g\n', ...
    results.segments.gamma.beam_stats.mean(end, c_r1), ...
    results.segments.gamma.beam_stats.mean(end, c_r3), ...
    results.segments.gamma.statistics(end).p_one_sided);
fprintf('  Output: %s\n', out_dir);
fprintf('=========================================================\n\n');

diary off;
movefile(tmp_log, fullfile(out_dir, 'step3_console_log.txt'), 'f');
clear log_guard;

end


%% =========================================================================
%  LOADING
%  =========================================================================

function d = load_ct_total_dose(processed_dir, ct_label)
% Per-CT RayStation total (sum of that CT's fields) written by Step 1.5.
    p = fullfile(processed_dir, sprintf('total_dose_%s.mat', ct_label));
    if ~isfile(p)
        error('step3_analysis:FileNotFound', ...
            '%s not found. Step 1.5 writes it when field filenames carry a CT label.', p);
    end
    L = load(p);
    if isfield(L, 'ct_total_sparse')
        d = reshape(full(L.ct_total_sparse), L.ct_total_dims);
    else
        d = L.ct_total;
    end
    d = double(d);
end


function sensor_mask = load_sensor_mask(sim_dir, hash8, dims)
% Plan-level sensor mask saved by pipeline_simulate. [] when absent or not on
% the dose grid (e.g. when the sensor method expanded the k-Wave grid).
    sensor_mask = [];
    p = fullfile(sim_dir, sprintf('sensor_mask_%s.mat', hash8));
    if ~isfile(p)
        fprintf('  [NOTE] No sensor_mask_%s.mat; sensor overlay skipped.\n', hash8);
        return;
    end
    L = load(p, 'precomputed_sensor');
    m = L.precomputed_sensor.sensor_mask;
    if ~isequal(size(m), dims)
        fprintf('  [NOTE] Sensor mask grid %s differs from dose grid %s; overlay skipped.\n', ...
            mat2str(size(m)), mat2str(dims));
        return;
    end
    sensor_mask = logical(m);
end


function SM = load_step25_summary(patient_id, session, config, sim_dir, hash8)
% Step 2.5 rollup for this hash; runs Step 2.5 (all beams) when it is missing.
    p = fullfile(sim_dir, sprintf('segment_metrics_summary_%s.mat', hash8));
    if ~isfile(p)
        fprintf('  No Step 2.5 summary for hash %s; running step25_segment_metrics...\n', hash8);
        c = config;
        c.metrics_beams         = [];
        c.metrics_write_summary = true;
        step25_segment_metrics(patient_id, session, c);
        if ~isfile(p)
            error('step3_analysis:NoSummary', ...
                'Step 2.5 ran but %s was not written. Check Step 2 recons exist.', p);
        end
    end
    fprintf('  Loading %s\n', p);
    SM = load(p);
end


%% =========================================================================
%  PART A: TOTAL-DOSE COMPARISON
%  =========================================================================

function cmp = compare_total_doses(reference, evaluated, spacing, config)
% Global gamma (CalcGamma) + local SSIM over the reference's 10% eval mask.
% Same mask, gamma search limit and SSIM as Step 2.5 so totals and segments
% are directly comparable.
    threshold = max(reference(:)) * config.gamma_dose_cutoff_pct / 100;
    mask      = reference >= threshold;
    fprintf('    Eval mask: %d voxels >= %.4f Gy\n', nnz(mask), threshold);

    % CalcGamma widths follow array dims: (row, col, slice) = (Y, X, Z).
    width   = [spacing(2), spacing(1), spacing(3)];
    ref_str = struct('start', [0, 0, 0], 'width', width, 'data', reference);
    tgt_str = struct('start', [0, 0, 0], 'width', width, 'data', evaluated);

    gamma_timer = tic;
    gamma_map = CalcGamma(ref_str, tgt_str, config.gamma_dose_pct, config.gamma_dist_mm, ...
        'local', 0, 'restrict', 1, 'limit', 2 * config.gamma_dist_mm);
    gamma_map = double(gather(gamma_map));
    g = gamma_map(mask);
    gamma_map(~mask) = NaN;

    cmp.gamma.pass_rate     = 100 * mean(g <= 1);
    cmp.gamma.mean_gamma    = mean(g, 'omitnan');
    cmp.gamma.max_gamma     = max(g);
    cmp.gamma.num_evaluated = numel(g);
    cmp.gamma.num_passed    = sum(g <= 1);
    cmp.gamma.gamma_map     = single(gamma_map);
    cmp.gamma.time_sec      = toc(gamma_timer);

    [ssim_map, ssim_mean] = compute_local_ssim(reference, evaluated, mask);
    if ~isempty(ssim_map)
        ssim_map(~mask) = NaN;
    end
    cmp.ssim.mean_pct = 100 * ssim_mean;
    cmp.ssim.ssim_map = single(ssim_map);

    cmp.threshold_Gy = threshold;
    fprintf('    Gamma pass %.1f%% (%d/%d), mean %.3f, max %.3f (%.0f s) | local SSIM %.1f%%\n', ...
        cmp.gamma.pass_rate, cmp.gamma.num_passed, cmp.gamma.num_evaluated, ...
        cmp.gamma.mean_gamma, cmp.gamma.max_gamma, cmp.gamma.time_sec, cmp.ssim.mean_pct);
end


function g = least_squares_gain(rs_truth, recon, cutoff_fr)
% Scalar g minimising ||rs - g*recon||^2 over the truth's cutoff region (1 if undefined).
    mask = rs_truth >= cutoff_fr * max(rs_truth(:));
    r = recon(mask);
    if sum(r .^ 2) > 0
        g = sum(rs_truth(mask) .* r) / sum(r .^ 2);
    else
        g = 1;
    end
end


%% =========================================================================
%  PART A: ORTHOGONAL-VIEW FIGURES
%  =========================================================================

function plot_total_comparison(reference, evaluated, cmp, ref_label, tgt_label, ...
        file_tag, spacing, body_mask, sensor_mask, fig_dir)
% Two figures at the reference max-dose voxel, views = transverse/coronal/sagittal.
%   dose_<tag>.png : rows reference | evaluated | difference
%   maps_<tag>.png : rows gamma map | local SSIM map | histograms + stats
    [~, lin] = max(reference(:));
    [iy, ix, iz] = ind2sub(size(reference), lin);
    ijk   = [iy, ix, iz];
    views = {'Transverse', 'Coronal', 'Sagittal'};

    dose_max = max(max(reference(:)), max(evaluated(:)));
    diff_vol = evaluated - reference;
    diff_lim = max(abs(diff_vol(:)));
    if diff_lim == 0, diff_lim = 1; end

    % ---- Figure 1: doses ----
    vols   = {reference, evaluated, diff_vol};
    labels = {ref_label, tgt_label, 'Difference (eval - ref)'};
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50, 50, 1350, 1050]);
    tl  = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for r = 1:3
        for v = 1:3
            ax = nexttile(tl);
            draw_view(ax, vols{r}, views{v}, ijk, spacing);
            if r < 3
                colormap(ax, jet(256)); clim(ax, [0, dose_max]);
            else
                colormap(ax, rdbu_colormap()); clim(ax, [-diff_lim, diff_lim]);
            end
            colorbar(ax);
            overlay_masks(ax, body_mask, sensor_mask, views{v}, ijk, spacing);
            title(ax, sprintf('%s | %s', labels{r}, views{v}), 'Interpreter', 'tex');
        end
    end
    title(tl, {sprintf('%s vs %s (Gy) at max-dose voxel [%d %d %d]', ref_label, tgt_label, ijk), ...
        'green = body contour | purple = sensor footprint'}, 'FontWeight', 'bold', 'Interpreter', 'tex');
    save_figure(fig, fullfile(fig_dir, sprintf('dose_%s.png', file_tag)));

    % ---- Figure 2: gamma + SSIM maps ----
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50, 50, 1350, 1050]);
    tl  = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for v = 1:3
        ax = nexttile(tl);
        draw_view(ax, double(cmp.gamma.gamma_map), views{v}, ijk, spacing);
        colormap(ax, gamma_colormap()); clim(ax, [0, 2]); colorbar(ax);
        overlay_masks(ax, body_mask, sensor_mask, views{v}, ijk, spacing);
        title(ax, sprintf('\\gamma | %s', views{v}));
    end
    for v = 1:3
        ax = nexttile(tl);
        if ~isempty(cmp.ssim.ssim_map)
            draw_view(ax, double(cmp.ssim.ssim_map), views{v}, ijk, spacing);
            colormap(ax, parula(256)); clim(ax, [0, 1]); colorbar(ax);
            overlay_masks(ax, body_mask, sensor_mask, views{v}, ijk, spacing);
        end
        title(ax, sprintf('Local SSIM | %s', views{v}));
    end

    ax = nexttile(tl);
    g = cmp.gamma.gamma_map(~isnan(cmp.gamma.gamma_map));
    histogram(ax, g, 0:0.05:3, 'FaceColor', [0.3, 0.5, 0.8], 'EdgeColor', 'none');
    xline(ax, 1, 'r--', 'LineWidth', 2);
    xlabel(ax, 'Gamma index'); ylabel(ax, 'Voxels'); xlim(ax, [0, 3]); grid(ax, 'on');
    title(ax, sprintf('Gamma pass %.1f%%', cmp.gamma.pass_rate));

    ax = nexttile(tl);
    if ~isempty(cmp.ssim.ssim_map)
        s = cmp.ssim.ssim_map(~isnan(cmp.ssim.ssim_map));
        histogram(ax, s, 0:0.02:1, 'FaceColor', [0.3, 0.6, 0.4], 'EdgeColor', 'none');
    end
    xlabel(ax, 'Local SSIM'); ylabel(ax, 'Voxels'); xlim(ax, [0, 1]); grid(ax, 'on');
    title(ax, sprintf('Mean local SSIM %.1f%%', cmp.ssim.mean_pct));

    ax = nexttile(tl); axis(ax, 'off');
    text(ax, 0, 0.5, { ...
        sprintf('Eval voxels: %d', cmp.gamma.num_evaluated), ...
        sprintf('Mask threshold: %.4f Gy', cmp.threshold_Gy), ...
        sprintf('Mean \\gamma: %.3f', cmp.gamma.mean_gamma), ...
        sprintf('Max \\gamma: %.3f', cmp.gamma.max_gamma)}, 'FontSize', 11);

    title(tl, sprintf('%s vs %s: \\gamma and local SSIM', ref_label, tgt_label), ...
        'FontWeight', 'bold', 'Interpreter', 'tex');
    save_figure(fig, fullfile(fig_dir, sprintf('maps_%s.png', file_tag)));
end


function [img, xv, yv] = ortho_slice(vol, view_name, ijk, spacing)
% 2D slice through voxel ijk = [iy, ix, iz]. Axes in mm from the grid corner.
% Transverse: rows Y (anterior top). Coronal / sagittal: rows Z (plotted
% with YDir normal so superior = higher iz is at the top).
    [ny, nx, nz] = size(vol);
    x = (0:nx - 1) * spacing(1);
    y = (0:ny - 1) * spacing(2);
    z = (0:nz - 1) * spacing(3);
    switch view_name
        case 'Transverse'
            img = vol(:, :, ijk(3));             xv = x; yv = y;
        case 'Coronal'
            img = squeeze(vol(ijk(1), :, :))';   xv = x; yv = z;
        case 'Sagittal'
            img = squeeze(vol(:, ijk(2), :))';   xv = y; yv = z;
    end
end


function draw_view(ax, vol, view_name, ijk, spacing)
    [img, xv, yv] = ortho_slice(vol, view_name, ijk, spacing);
    h = imagesc(ax, xv, yv, img);
    set(h, 'AlphaData', ~isnan(img));        % NaN (outside eval mask) = blank
    axis(ax, 'image');
    if ~strcmp(view_name, 'Transverse')
        set(ax, 'YDir', 'normal');
    end
    xlabel(ax, 'mm'); ylabel(ax, 'mm');
end


function overlay_masks(ax, body_mask, sensor_mask, view_name, ijk, spacing)
% Body contour on the displayed slice (green); sensor footprint projected
% along the viewing axis (purple) so it shows even when the slice misses it.
    hold(ax, 'on');
    if ~isempty(body_mask)
        [m, xv, yv] = ortho_slice(body_mask, view_name, ijk, spacing);
        if any(m(:))
            contour(ax, xv, yv, double(m), [0.5, 0.5], 'Color', [0, 0.8, 0], 'LineWidth', 1.2);
        end
    end
    if ~isempty(sensor_mask)
        switch view_name
            case 'Transverse', proj = any(sensor_mask, 3);
            case 'Coronal',    proj = squeeze(any(sensor_mask, 1))';
            case 'Sagittal',   proj = squeeze(any(sensor_mask, 2))';
        end
        [~, xv, yv] = ortho_slice(sensor_mask, view_name, ijk, spacing);
        [r, c] = find(proj);
        plot(ax, xv(c), yv(r), '.', 'Color', [0.7, 0, 1], 'MarkerSize', 4);
    end
    hold(ax, 'off');
end


%% =========================================================================
%  PART B: STATISTICS (base MATLAB only)
%  =========================================================================

function B = beam_stats(seg_by_beam)
% Per-beam mean/SE over segments; the last row is all segments pooled.
    nB = numel(seg_by_beam);
    nC = size(seg_by_beam{1}, 2);
    B.mean = nan(nB + 1, nC);
    B.std  = nan(nB + 1, nC);
    B.se   = nan(nB + 1, nC);
    B.n    = zeros(nB + 1, 1);
    for n = 1:nB + 1
        if n <= nB
            M = seg_by_beam{n};
        else
            M = vertcat(seg_by_beam{:});
        end
        B.n(n)       = size(M, 1);
        B.mean(n, :) = mean(M, 1, 'omitnan');
        B.std(n, :)  = std(M, 0, 1, 'omitnan');
        B.se(n, :)   = B.std(n, :) ./ sqrt(max(sum(isfinite(M), 1), 1));
    end
end


function D = differential_stats(seg_by_beam, c_r1, c_r3, c_tt)
% Paired per-segment differentials, mean/SE per beam (+ pooled last row):
%   recon = (recon1 vs truth1) - (recon3 vs truth1);  truth = 100 - (truth1 vs truth3).
    nB = numel(seg_by_beam);
    D.recon_mean = nan(nB + 1, 1); D.recon_se = nan(nB + 1, 1);
    D.truth_mean = nan(nB + 1, 1); D.truth_se = nan(nB + 1, 1);
    for n = 1:nB + 1
        if n <= nB
            M = seg_by_beam{n};
        else
            M = vertcat(seg_by_beam{:});
        end
        rd = M(:, c_r1) - M(:, c_r3);
        td = 100 - M(:, c_tt);
        D.recon_mean(n) = mean(rd, 'omitnan');
        D.recon_se(n)   = std(rd, 'omitnan') / sqrt(max(sum(isfinite(rd)), 1));
        D.truth_mean(n) = mean(td, 'omitnan');
        D.truth_se(n)   = std(td, 'omitnan') / sqrt(max(sum(isfinite(td)), 1));
    end
end


function stats = segment_statistics(seg_by_beam, c_r1, c_r3, c_tt, snr_seg, floor_pct, labels)
% Per beam (+ pooled last element):
%   (1) paired one-sided t-test, H1: mean(A - B) > 0  (A = recon1, B = recon3)
%   (2) % of segments above floor_pct for A and B, with SE of the proportion
%   (3) Pearson R^2 / Spearman rho of (A - B) vs truth1_vs_truth3 and vs CT_1 SNR
    nB = numel(seg_by_beam);
    pooled_A = []; pooled_B = []; pooled_T = []; pooled_S = [];
    for n = 1:nB + 1
        if n <= nB
            A = seg_by_beam{n}(:, c_r1);
            B = seg_by_beam{n}(:, c_r3);
            T = seg_by_beam{n}(:, c_tt);
            S = snr_seg{n}(:);
            pooled_A = [pooled_A; A]; pooled_B = [pooled_B; B]; %#ok<AGROW>
            pooled_T = [pooled_T; T]; pooled_S = [pooled_S; S]; %#ok<AGROW>
        else
            A = pooled_A; B = pooled_B; T = pooled_T; S = pooled_S;
        end
        d = A - B;

        s.beam = labels{n};
        [s.t_stat, s.p_one_sided, s.pct_a_gt_b, s.n_pairs] = paired_t_one_sided(A, B);
        [s.above_floor_a_pct, s.above_floor_a_se] = proportion_above(A, floor_pct);
        [s.above_floor_b_pct, s.above_floor_b_se] = proportion_above(B, floor_pct);
        s.r2_vs_truth  = pearson_r(d, T) ^ 2;
        s.rho_vs_truth = pearson_r(rank_average(d, T), rank_average(T, d));
        s.r2_vs_snr    = pearson_r(d, S) ^ 2;
        s.rho_vs_snr   = pearson_r(rank_average(d, S), rank_average(S, d));
        s.pooled_diff  = d;
        s.pooled_truth = T;
        s.pooled_snr   = S;
        stats(n) = s; %#ok<AGROW>
    end
end


function [t_stat, p, pct_gt, n] = paired_t_one_sided(A, B)
% Paired t-test, one-sided upper tail, via the Student-t CDF written with
% betainc:  P(T > t) = 0.5 * betainc(df / (df + t^2), df/2, 1/2)  for t >= 0.
    ok = isfinite(A) & isfinite(B);
    d  = A(ok) - B(ok);
    n  = numel(d);
    t_stat = NaN; p = NaN; pct_gt = NaN;
    if n >= 1, pct_gt = 100 * mean(d > 0); end
    if n < 2 || std(d) == 0, return; end
    df     = n - 1;
    t_stat = mean(d) / (std(d) / sqrt(n));
    tail   = 0.5 * betainc(df / (df + t_stat ^ 2), df / 2, 0.5);
    if t_stat >= 0
        p = tail;
    else
        p = 1 - tail;
    end
end


function [pct, se] = proportion_above(x, floor_pct)
    x = x(isfinite(x));
    if isempty(x) || ~isfinite(floor_pct)
        pct = NaN; se = NaN;
        return;
    end
    p   = mean(x > floor_pct);
    pct = 100 * p;
    se  = 100 * sqrt(p * (1 - p) / numel(x));
end


function r = pearson_r(x, y)
% Pearson r over pairwise-finite entries; NaN with < 3 pairs or no spread.
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    r = NaN;
    if numel(x) < 3 || std(x) == 0 || std(y) == 0, return; end
    C = corrcoef(x, y);
    r = C(1, 2);
end


function rk = rank_average(x, y)
% Ranks of x (ties averaged) over the entries where BOTH x and y are finite;
% NaN elsewhere. Pearson r of two such rank vectors = Spearman rho.
    rk = nan(size(x));
    ok = find(isfinite(x) & isfinite(y));
    v  = x(ok);
    [sorted, order] = sort(v);
    r  = zeros(size(v));
    r(order) = 1:numel(v);
    for u = unique(sorted)'
        tie = (v == u);
        r(tie) = mean(r(tie));
    end
    rk(ok) = r;
end


function snr = gather_snr(sim_dir, hash8, plan_type, ct1_str, beam_list, n_seg)
% noise_stats.snr from each recon file (header + small variable only).
%   .mean/.se/.n : per beam over all its fields (both CTs), last row pooled
%   .seg_ct1     : {beam} per-segment CT_1 SNR, sorted by segment, aligned
%                  to the Step 2.5 rows (NaN column when counts disagree)
    nB       = numel(beam_list);
    all_vals = cell(1, nB);
    ct1_seg  = cell(1, nB);
    ct1_val  = cell(1, nB);
    ct_pat   = ['_', regexprep(ct1_str, '_', '[_-]?'), '_'];

    listing = dir(fullfile(sim_dir, sprintf('*_recon_%s.mat', hash8)));
    for k = 1:numel(listing)
        name = listing(k).name;
        tok  = regexp(name, '_B(\d+)_(\d+)_recon_', 'tokens', 'once');
        if isempty(tok), continue; end
        slot = find(beam_list == str2double(tok{1}), 1);
        if isempty(slot), continue; end
        if ~strcmpi(plan_type, 'any')
            ptok = regexp(name, '_(adapted|reference)_', 'tokens', 'once');
            if isempty(ptok) || ~strcmpi(ptok{1}, plan_type), continue; end
        end

        fpath = fullfile(sim_dir, name);
        vars  = who('-file', fpath);
        s = NaN;
        if ismember('noise_stats', vars)
            L = load(fpath, 'noise_stats');
            if isfield(L.noise_stats, 'snr') && ~isempty(L.noise_stats.snr)
                s = double(L.noise_stats.snr);
            end
        end
        if isfinite(s)
            all_vals{slot}(end + 1) = s;
        end
        % Step 2.5 rows = CT_1 recon files carrying a segment_metrics fold.
        if ~isempty(regexp(name, ct_pat, 'once')) && ismember('segment_metrics', vars)
            ct1_seg{slot}(end + 1) = str2double(tok{2});
            ct1_val{slot}(end + 1) = s;
        end
    end

    snr.mean = nan(nB + 1, 1); snr.se = nan(nB + 1, 1); snr.n = zeros(nB + 1, 1);
    for n = 1:nB + 1
        if n <= nB
            v = all_vals{n};
        else
            v = [all_vals{:}];
        end
        snr.n(n) = numel(v);
        if ~isempty(v)
            snr.mean(n) = mean(v);
            snr.se(n)   = std(v) / sqrt(numel(v));
        end
    end

    snr.seg_ct1 = cell(1, nB);
    for n = 1:nB
        [~, ord] = sort(ct1_seg{n});
        if numel(ord) == n_seg(n)
            snr.seg_ct1{n} = ct1_val{n}(ord)';
        else
            warning('step3_analysis:SNRAlignSkip', ...
                'Beam #%d: %d CT_1 SNR value(s) vs %d Step 2.5 segment(s); SNR correlation skipped.', ...
                beam_list(n), numel(ord), n_seg(n));
            snr.seg_ct1{n} = nan(n_seg(n), 1);
        end
    end
end


%% =========================================================================
%  PART B: BEAM-AXIS FIGURES
%  =========================================================================

function plot_pass_rates(x_labels, B, c_r1, c_r3, c_tt, noise_floor, mlabel, mword, title_tag, fpath)
% recon1 vs truth1 (green), recon3 vs truth1 (red), truth1 vs truth3 (blue),
% connector between the recon means, noise-only null band (gamma only).
    green = [0.15, 0.60, 0.20]; red = [0.80, 0.15, 0.15];
    blue  = [0.20, 0.40, 0.80]; grey = [0.35, 0.35, 0.35];
    nX = numel(x_labels);
    x  = (1:nX)';   % column, same shape as the B.mean columns

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, max(760, 55 * nX + 240), 500]);
    ax = axes(fig); hold(ax, 'on');
    h = gobjects(0); leg = {};
    if ~isempty(noise_floor)
        m = noise_floor.mean_pass_rate; s = noise_floor.std_pass_rate;
        patch(ax, [0.5, nX + 0.5, nX + 0.5, 0.5], [m - s, m - s, m + s, m + s], grey, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
        h(end + 1) = plot(ax, [0.5, nX + 0.5], [m, m], ':', 'Color', grey, 'LineWidth', 1.5);
        leg{end + 1} = sprintf('Noise-only %.1f \\pm %.1f%%', m, s);
    end
    for n = 1:nX
        a = B.mean(n, c_r1); b = B.mean(n, c_r3);
        if a >= b, lc = green; else, lc = red; end
        plot(ax, [n, n], [a, b], '-', 'Color', lc, 'LineWidth', 1.5);
    end
    h(end + 1) = errorbar(ax, x, B.mean(:, c_r1), B.se(:, c_r1), 'o', 'Color', green, ...
        'MarkerFaceColor', green, 'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7);
    h(end + 1) = errorbar(ax, x, B.mean(:, c_r3), B.se(:, c_r3), 'o', 'Color', red, ...
        'MarkerFaceColor', red, 'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7);
    h(end + 1) = errorbar(ax, x, B.mean(:, c_tt), B.se(:, c_tt), 'o', 'Color', blue, ...
        'MarkerFaceColor', blue, 'MarkerSize', 7, 'LineStyle', 'none', 'CapSize', 7);
    leg = [leg, {'Recon CT\_1 vs Truth CT\_1', 'Recon CT\_3 vs Truth CT\_1', 'Truth CT\_1 vs Truth CT\_3'}];
    yline(ax, 90, 'k--', '90%');
    beam_axis(ax, x_labels);
    ylim(ax, [0, 105]); ylabel(ax, mlabel);
    title(ax, sprintf('%s by beam (mean \\pm SE over segments) | %s', mword, title_tag));
    legend(ax, h, leg, 'Location', 'best', 'FontSize', 9);
    save_figure(fig, fpath);
end


function plot_differentials(x_labels, D, m_r1, noise_floor, mword, title_tag, fpath)
% recon1 - recon3 (green >= 0, red < 0), 100 - truth1_vs_truth3 (blue),
% recon1 - noise-only null (grey, gamma only).
    green = [0.15, 0.60, 0.20]; red = [0.80, 0.15, 0.15];
    blue  = [0.20, 0.40, 0.80]; grey = [0.35, 0.35, 0.35];
    nX = numel(x_labels);
    x  = (1:nX)';

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, max(760, 55 * nX + 240), 500]);
    ax = axes(fig); hold(ax, 'on');
    yline(ax, 0, 'k-');
    for n = 1:nX
        if D.recon_mean(n) >= 0, c = green; else, c = red; end
        errorbar(ax, n, D.recon_mean(n), D.recon_se(n), 'o', 'Color', c, 'MarkerFaceColor', c, ...
            'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7);
    end
    h = [plot(ax, nan, nan, 'o', 'Color', green, 'MarkerFaceColor', green), ...
         plot(ax, nan, nan, 'o', 'Color', red, 'MarkerFaceColor', red)];
    h(end + 1) = errorbar(ax, x, D.truth_mean, D.truth_se, 's', 'Color', blue, ...
        'MarkerFaceColor', blue, 'MarkerSize', 7, 'LineStyle', 'none', 'CapSize', 7);
    leg = {'Recon CT\_1 - CT\_3 \geq 0', 'Recon CT\_1 - CT\_3 < 0', '100 - (Truth CT\_1 vs CT\_3)'};
    if ~isempty(noise_floor)
        h(end + 1) = plot(ax, x, m_r1 - noise_floor.mean_pass_rate, 'd', 'Color', grey, ...
            'MarkerFaceColor', grey, 'LineStyle', 'none');
        leg{end + 1} = 'Recon CT\_1 - noise-only';
    end
    beam_axis(ax, x_labels);
    ylabel(ax, sprintf('%s differential (%%)', mword));
    title(ax, sprintf('Change-detection differentials (mean \\pm SE) | %s', title_tag));
    legend(ax, h, leg, 'Location', 'best', 'FontSize', 9);
    save_figure(fig, fpath);
end


function plot_fidelity(x_labels, B, c_r1, c_33, mlabel, title_tag, fpath)
% Each recon vs its OWN-CT truth: truth1/recon1 (green), truth3/recon3 (red).
    green = [0.15, 0.60, 0.20]; red = [0.80, 0.15, 0.15];
    nX  = numel(x_labels);
    gap = @(v) [v(1:end - 1); NaN; v(end)];    % NaN gap: no line into "All"
    xg  = gap((1:nX)');

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, max(720, 55 * nX + 220), 480]);
    ax = axes(fig); hold(ax, 'on');
    h1 = errorbar(ax, xg, gap(B.mean(:, c_r1)), gap(B.se(:, c_r1)), '-o', 'Color', green, ...
        'MarkerFaceColor', green, 'MarkerSize', 8, 'CapSize', 7);
    h2 = errorbar(ax, xg, gap(B.mean(:, c_33)), gap(B.se(:, c_33)), '-o', 'Color', red, ...
        'MarkerFaceColor', red, 'MarkerSize', 8, 'CapSize', 7);
    yline(ax, 90, 'k--', '90%');
    beam_axis(ax, x_labels);
    ylim(ax, [0, 105]); ylabel(ax, mlabel);
    title(ax, sprintf('Reconstruction fidelity vs own-CT truth (mean \\pm SE) | %s', title_tag));
    legend(ax, [h1, h2], {'Truth CT\_1 vs Recon CT\_1', 'Truth CT\_3 vs Recon CT\_3'}, ...
        'Location', 'best', 'FontSize', 9);
    save_figure(fig, fpath);
end


function plot_beam_snr(x_labels, snr, patient_id, session, fpath)
    purple = [0.45, 0.25, 0.65];
    x = (1:numel(x_labels))';
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, max(760, 55 * numel(x) + 240), 480]);
    ax = axes(fig); hold(ax, 'on');
    if any(isfinite(snr.mean))
        errorbar(ax, x, snr.mean, snr.se, 'o', 'Color', purple, 'MarkerFaceColor', purple, ...
            'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7);
    else
        text(ax, 0.5, 0.5, 'No noise\_stats SNR in the recon files', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
    end
    beam_axis(ax, x_labels);
    ylabel(ax, 'SNR (signal peak / noise amplitude)');
    title(ax, sprintf('Electronic-noise SNR by beam (mean \\pm SE over fields) | %s / %s', ...
        strrep(patient_id, '_', '\_'), strrep(session, '_', '\_')));
    save_figure(fig, fpath);
end


function plot_correlations(d, T, S, mword, title_tag, fpath)
% Pooled scatter of the recon differential vs the true change and vs CT_1 SNR.
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1000, 440]);
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: recon differential correlations | %s', mword, title_tag), ...
        'FontWeight', 'bold');
    xs   = {T, S};
    xlab = {'Truth CT\_1 vs Truth CT\_3 (%)', 'CT\_1 SNR'};
    for k = 1:2
        ax = nexttile(tl); hold(ax, 'on');
        x  = xs{k};
        ok = isfinite(x) & isfinite(d);
        scatter(ax, x(ok), d(ok), 18, [0.20, 0.40, 0.80], 'filled', 'MarkerFaceAlpha', 0.5);
        if nnz(ok) >= 3 && max(x(ok)) > min(x(ok))
            coef = polyfit(x(ok), d(ok), 1);
            xl   = [min(x(ok)), max(x(ok))];
            plot(ax, xl, polyval(coef, xl), 'r-', 'LineWidth', 1.5);
        end
        rho = pearson_r(rank_average(d, x), rank_average(x, d));
        title(ax, sprintf('R^2 = %.3f | Spearman \\rho = %+.3f (n = %d)', ...
            pearson_r(d, x) ^ 2, rho, nnz(ok)), 'FontWeight', 'normal');
        xlabel(ax, xlab{k}); ylabel(ax, 'Recon CT\_1 - CT\_3 (%)');
        grid(ax, 'on'); box(ax, 'on');
    end
    save_figure(fig, fpath);
end


function beam_axis(ax, x_labels)
% Beam tick labels plus a dotted rule before the pooled "All" slot.
    n = numel(x_labels);
    xlim(ax, [0.5, n + 0.5]);
    xticks(ax, 1:n);
    xticklabels(ax, x_labels);
    xline(ax, n - 0.5, ':', 'Color', [0.4, 0.4, 0.4]);
    xlabel(ax, 'Beam number');
    grid(ax, 'on'); box(ax, 'on');
end


%% =========================================================================
%  PART B: RANDOM SEGMENT PANELS
%  =========================================================================

function picks = select_random_segments(patient_id, session, config, beam_list, ct1_str, ct3_str)
% Seeded random segments. For each: axial slice (reference max dose) of
% truth CT_1, recon CT_1 / CT_3, and the Step 2.5 folded gamma + SSIM maps
% for truth1_vs_recon1 and truth1_vs_recon3. Loads one beam at a time.
    picks = struct('beam', {}, 'seg', {}, 'iz', {}, 'scores', {}, 'rows', {});
    rng(config.analysis_random_seed);
    beam_order = beam_list(randperm(numel(beam_list)));
    comps = {'truth1_vs_recon1', 'Recon CT\_1'; 'truth1_vs_recon3', 'Recon CT\_3'};

    for b = beam_order
        if numel(picks) >= config.analysis_n_random, break; end
        args = {'Mode', 'set', 'Beam', b, 'IncludeEthos', false, 'IncludeCBCT', false, ...
                'Hash', config.config_hash};
        if ~strcmpi(config.metrics_plan_type, 'any')
            args = [args, {'PlanType', config.metrics_plan_type}]; %#ok<AGROW>
        end
        try
            S = load_recon_dose_data(patient_id, session, config, args{:});
        catch ME
            warning('step3_analysis:RandomSkipBeam', 'Random panels: beam #%d skipped: %s', b, ME.message);
            continue;
        end

        segs = arrayfun(@(f) double(f.rtplan.seg_num), S.fields);
        cts  = arrayfun(@(f) strrep(char(f.rtplan.ct_label), '-', '_'), S.fields, 'UniformOutput', false);
        useg = unique(segs(~isnan(segs)));
        useg = useg(randperm(numel(useg)));

        for s = useg(:)'
            if numel(picks) >= config.analysis_n_random, break; end
            i1 = find(segs == s & strcmp(cts, ct1_str), 1);
            i3 = find(segs == s & strcmp(cts, ct3_str), 1);
            if isempty(i1) || isempty(i3), continue; end

            rs1    = double(S.fields(i1).rs_dose);
            recons = {double(S.fields(i1).recon_dose), double(S.fields(i3).recon_dose)};
            if config.metrics_normalize
                cut = config.gamma_dose_cutoff_pct / 100;
                recons{1} = recons{1} * least_squares_gain(rs1, recons{1}, cut);
                recons{2} = recons{2} * least_squares_gain(double(S.fields(i3).rs_dose), recons{2}, cut);
            end
            [~, lin] = max(rs1(:));
            [~, ~, iz] = ind2sub(size(rs1), lin);

            pk.beam   = b;
            pk.seg    = s;
            pk.iz     = iz;
            pk.scores = struct('comparison', {}, 'gamma_pass_pct', {}, 'ssim_mean_pct', {});
            pk.rows   = struct('label', {}, 'truth', {}, 'recon', {}, 'gamma', {}, 'ssim', {});
            for c = 1:2
                [gmap, gpass, smap, smean] = folded_maps(S.fields(i1).recon_file, comps{c, 1}, ...
                    config.config_hash, size(rs1));
                pk.rows(c).label = comps{c, 2};
                pk.rows(c).truth = rs1(:, :, iz);
                pk.rows(c).recon = recons{c}(:, :, iz);
                pk.rows(c).gamma = gmap(:, :, iz);
                pk.rows(c).ssim  = smap(:, :, iz);
                pk.scores(c).comparison     = comps{c, 1};
                pk.scores(c).gamma_pass_pct = gpass;
                pk.scores(c).ssim_mean_pct  = smean;
            end
            picks(end + 1) = pk; %#ok<AGROW>
        end
        clear S;
    end
    fprintf('  Selected %d random segment(s) (requested %d).\n', numel(picks), config.analysis_n_random);
end


function [gmap, gpass, smap, smean] = folded_maps(recon_file, comp_name, hash8, vol_size)
% Re-expand the Step 2.5 folded gamma / SSIM values (NaN outside the eval
% mask). All-NaN maps when the fold is missing or from another hash.
    gmap = nan(vol_size); smap = nan(vol_size); gpass = NaN; smean = NaN;
    vars = who('-file', recon_file);
    if ~ismember('segment_metrics', vars)
        warning('step3_analysis:NoFold', 'No segment_metrics in %s.', recon_file);
        return;
    end
    L  = load(recon_file, 'segment_metrics');
    sm = L.segment_metrics;
    if ~strcmpi(char(sm.config_hash), hash8) || ~isequal(sm.vol_size(:)', vol_size)
        return;
    end
    k = find(strcmp({sm.comparison.name}, comp_name), 1);
    if isempty(k), return; end
    c = sm.comparison(k);
    if ~isempty(c.gamma_vals), gmap(c.mask_idx) = double(c.gamma_vals); end
    if ~isempty(c.ssim_vals),  smap(c.mask_idx) = double(c.ssim_vals);  end
    gpass = c.gamma_pass_rate;
    smean = c.ssim_mean;
end


function plot_random_segment(pk, crit, fpath)
% 2 rows (recon CT_1, recon CT_3) x 4 columns (truth | recon | gamma | SSIM).
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1400, 650]);
    tl  = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Beam %d, Segment %d | axial slice %d (truth CT\\_1 max)', ...
        pk.beam, pk.seg, pk.iz), 'FontWeight', 'bold');
    for r = 1:numel(pk.rows)
        row  = pk.rows(r);
        cmax = max(row.truth(:));
        if ~(cmax > 0), cmax = 1; end

        ax = nexttile(tl); imagesc(ax, row.truth); axis(ax, 'image', 'off');
        colormap(ax, hot(256)); clim(ax, [0, cmax]); colorbar(ax);
        title(ax, 'Truth CT\_1');

        ax = nexttile(tl); imagesc(ax, row.recon); axis(ax, 'image', 'off');
        colormap(ax, hot(256)); clim(ax, [0, cmax]); colorbar(ax);
        title(ax, row.label);

        ax = nexttile(tl); h = imagesc(ax, row.gamma); axis(ax, 'image', 'off');
        set(h, 'AlphaData', ~isnan(row.gamma));
        colormap(ax, gamma_colormap()); clim(ax, [0, 2]); colorbar(ax);
        title(ax, sprintf('\\gamma %g%%/%gmm (%.1f%% \\leq 1)', crit, crit, pk.scores(r).gamma_pass_pct));

        ax = nexttile(tl); h = imagesc(ax, row.ssim); axis(ax, 'image', 'off');
        set(h, 'AlphaData', ~isnan(row.ssim));
        colormap(ax, parula(256)); clim(ax, [0, 1]); colorbar(ax);
        title(ax, sprintf('Local SSIM (mean %.1f%%)', pk.scores(r).ssim_mean_pct));
    end
    save_figure(fig, fpath);
end


%% =========================================================================
%  OUTPUT HELPERS
%  =========================================================================

function save_figure(fig, fpath)
    exportgraphics(fig, fpath, 'Resolution', 150);
    close(fig);
    fprintf('    Saved figure: %s\n', fpath);
end


function write_beam_csv(fpath, results, beam_list)
% One row per beam (+ All): segment count, mean/SE of each comparison for
% both metrics, SNR mean/SE.
    labels = results.segments.gamma.beam_labels(:);
    comps  = results.segments.gamma.comparisons;
    T = table(labels, results.segments.gamma.beam_stats.n, 'VariableNames', {'Beam', 'NSegments'});
    keys = {'gamma', 'ssim'};
    for m = 1:2
        Bs = results.segments.(keys{m}).beam_stats;
        for c = 1:numel(comps)
            T.(sprintf('%s_%s_mean', keys{m}, comps{c})) = Bs.mean(:, c);
            T.(sprintf('%s_%s_se',   keys{m}, comps{c})) = Bs.se(:, c);
        end
    end
    T.snr_mean = results.snr.mean;
    T.snr_se   = results.snr.se;
    T.snr_n    = results.snr.n;
    writetable(T, fpath);
    fprintf('    Saved: %s (%d beams + All)\n', fpath, numel(beam_list));
end


function write_stats_csv(fpath, results)
    keys = {'gamma', 'ssim'};
    T = table();
    for m = 1:2
        S  = struct2table(results.segments.(keys{m}).statistics(:));
        S.metric    = repmat(keys(m), height(S), 1);
        S.floor_pct = repmat(results.segments.(keys{m}).floor_pct, height(S), 1);
        T = [T; S]; %#ok<AGROW>
    end
    T = movevars(T, 'metric', 'Before', 1);
    writetable(T, fpath);
    fprintf('    Saved: %s\n', fpath);
end


function cmap = rdbu_colormap()
    % Blue (negative) -> white -> red (positive) for dose differences.
    half = 128;
    r = [linspace(0.1, 1, half)'; linspace(1, 0.8, half)'];
    g = [linspace(0.2, 1, half)'; linspace(1, 0.1, half)'];
    b = [linspace(0.7, 1, half)'; linspace(1, 0.1, half)'];
    cmap = [r, g, b];
end


function cmap = gamma_colormap()
    % Solid green for gamma <= 1 (pass), green -> red over 1 < gamma <= 2.
    x     = linspace(0, 2, 256)';
    green = [0.15, 0.60, 0.20];
    red   = [0.80, 0.15, 0.15];
    f     = min(max(x - 1, 0), 1);
    cmap  = (1 - f) * green + f * red;
end
