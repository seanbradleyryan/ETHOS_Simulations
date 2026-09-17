%% =========================================================================
%  STUDY_PASS_RATES_ALLSEGMENTS.m
%  Per-BEAM photoacoustic gamma / SSIM summary, averaged over ALL segments.
%
%  This script is now a PLOTTER: it DEFERS every gamma-index and SSIM
%  calculation to Step 2.5 (step25_segment_metrics.m), which runs after the
%  k-Wave reconstruction and stores the four per-segment comparisons keyed to
%  the simulation config hash. Here we only LOAD those precomputed results and
%  draw the summary figures + optional random per-segment panels. Nothing is
%  recomputed -- a segment with no Step-2.5 result on disk is skipped / blank.
%
%  Where Step 2.5 stores its results (read here):
%    - segment_metrics_summary_<hash>.mat  (beside the recon doses in
%      SimulationResults/[PatientID]/[Session]/[method]/) -- per-beam & pooled
%      mean/std of the gamma pass rate (%) and mean local SSIM (%) per
%      comparison, the per-segment matrices, and the noise-only null floor.
%      This drives the summary figures and the console table.
%    - each segment's CT_1 recon .mat carries a 'segment_metrics' variable with
%      the masked-region gamma-index / SSIM values (mask_idx + *_vals). The
%      random per-segment panels re-expand these maps on demand.
%
%  EVALUATION METRIC (CONFIG.eval_method):
%    'gamma_index' (default) - per-segment global gamma pass rate (%).
%    'ssim'                  - per-segment mean local SSIM (%) over the 10%
%                              reference eval mask.
%  Both are read from the SAME summary file; the criterion is whatever Step 2.5
%  used (CONFIG.gamma_n is validated against it, never used to recompute).
%
%  Comparisons (per segment, in the summary's stored order):
%    truth1_vs_truth3  RayStation truth CT_1 vs truth CT_3
%    truth1_vs_recon1  RayStation truth CT_1 vs recon CT_1
%    truth1_vs_recon3  RayStation truth CT_1 vs recon CT_3
%    truth3_vs_recon3  RayStation truth CT_3 vs recon CT_3
%
%  FIGURE 1 - change detection, a TABBED window (x = beam number):
%    Tab "Pass rates"    - recon_CT1 (green), recon_CT3 (red) and truth_CT1-vs-
%                          truth_CT3 (blue), std error bars, a per-beam connector
%                          between the two recon means, and (gamma only, when
%                          CONFIG.include_noise_floor) the noise-only null band.
%    Tab "Differentials" - recon_CT1 - recon_CT3 (green >=0 / red <0), the true
%                          change 100 - (truth_CT1 vs truth_CT3) in blue, and
%                          (gamma only) recon_CT1 - noise.
%  FIGURE 2 - reconstruction fidelity (each recon vs its OWN-CT truth).
%  Both carry a trailing pooled "All" entry over every plotted segment.
%
%  OPTIONAL random per-segment panels: N random segments, each a tab of
%  truth | recon | metric-map rows. The metric map is the Step-2.5 folded
%  gamma/SSIM map re-expanded from disk (never recomputed); recon volumes are
%  loaded only for the truth/recon image columns.
%
%  NOTE: HIPAA / remote-execution - this file is WRITTEN here but must be RUN on
%  the remote device. Do not execute locally.
%  =========================================================================

clear; clc; close all;

% Script lives one level below the repo root; add utils/ and pipeline/ from there.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'utils')));
addpath(genpath(fullfile(repoRoot, 'pipeline')));
run_timer = tic;   % program runtime record

%% ========================= CONFIGURATION ================================

% Which precomputed per-segment metric drives the charts / console table:
%   'gamma_index' - global gamma pass rate (%);  'ssim' - mean local SSIM (%).
% Both are read from the Step-2.5 summary; the SSIM run uses *_ssim-suffixed
% log/cache filenames so the two metrics never clobber each other.
CONFIG.eval_method = 'gamma_index';   % 'gamma_index' | 'ssim'

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';
CONFIG.treatment_site = 'Pancreas';

% Beams to summarize (must be a subset of the beams Step 2.5 processed).
CONFIG.beams = 1:17;

% Restrict to a single plan type (validated against the summary).
CONFIG.plan_type = 'reference';   % 'reference' | 'adapted' | 'any'

% The two CT image indices (lower -> *_CT1 volumes, higher -> *_CT3 volumes).
CONFIG.ct_pair = [1, 3];

% Recon config-hash. '' => auto-discover the single summary on disk.
CONFIG.config_hash = 'a9a3e1e6';

% Needed by load_recon_dose_data (random panels) to resolve the method folder.
CONFIG.gruneisen_method = 'threshold_2';

% Comparison display labels {name, ref, tgt, ref_label, tgt_label}. The column
% order here only labels the random panels; the summary's own column order
% drives the summary figures.
CONFIG.comparisons = { ...
    'truth1_vs_truth3', 'rs_CT1', 'rs_CT3',    'Truth CT\_1', 'Truth CT\_3'; ...
    'truth1_vs_recon1', 'rs_CT1', 'recon_CT1', 'Truth CT\_1', 'Recon CT\_1'; ...
    'truth1_vs_recon3', 'rs_CT1', 'recon_CT3', 'Truth CT\_1', 'Recon CT\_3'; ...
    'truth3_vs_recon3', 'rs_CT3', 'recon_CT3', 'Truth CT\_3', 'Recon CT\_3'  ...
};

% Gamma criterion n (n%/n mm) EXPECTED from Step 2.5. Validated against the
% summary; a mismatch warns and the summary's value is used (never recomputed).
CONFIG.gamma_n = 3;

% Least-squares recon->own-truth normalization. Must match what Step 2.5 used
% (validated); also applied to the recon volumes shown in the random panels.
CONFIG.normalize = true;

% Draw the Step-2.5 noise-only null floor (gamma metric only) when present.
CONFIG.include_noise_floor = true;

% --- Output / logging ---
% The console log is written beside the recon doses unless output_dir is set.
CONFIG.output_dir = '';   % '' => the recon-dose directory
CONFIG.log_file   = 'pass_rates_allsegments_log.txt';

% --- Random per-segment panel visualization ---
% N random segments, each a tab of truth | recon | metric-map rows. The metric
% map is the Step-2.5 folded gamma/SSIM map re-expanded from disk (never
% recomputed); recon volumes are loaded only for the truth/recon image columns.
% Selected panels are cached so a matching re-run replots without reloading.
CONFIG.plot_random_segments = true;
CONFIG.n_random             = 5;     % number of random segments to display
CONFIG.random_seed          = 42;    % reproducible selection (in the cache key)
CONFIG.random_comparisons   = {'truth1_vs_recon1', 'truth1_vs_recon3'};
CONFIG.random_cache_file    = 'pass_rates_allsegments_random_viz.mat';
CONFIG.use_cache            = true;  % reuse the random-panel cache when it matches
CONFIG.save_results         = true;  % save the selected random panels to that cache

%% ===================== SETUP ============================================

% Normalize / validate the evaluation metric, and give the SSIM run its own
% cache/log filenames so the two metrics do not overwrite each other on disk.
CONFIG.eval_method = lower(char(CONFIG.eval_method));
if ~ismember(CONFIG.eval_method, {'gamma_index', 'ssim'})
    error('study_pass_rates_allsegments:BadEvalMethod', ...
        'CONFIG.eval_method must be ''gamma_index'' or ''ssim'' (got ''%s'').', ...
        CONFIG.eval_method);
end
if strcmp(CONFIG.eval_method, 'ssim')
    CONFIG.log_file          = add_name_suffix(CONFIG.log_file, '_ssim');
    CONFIG.random_cache_file = add_name_suffix(CONFIG.random_cache_file, '_ssim');
end

% Cache/log directory: the recon-dose folder unless explicitly overridden.
if isempty(CONFIG.output_dir)
    CONFIG.output_dir = fullfile(CONFIG.working_dir, 'SimulationResults', ...
        CONFIG.patient_id, CONFIG.session, CONFIG.gruneisen_method);
end
if ~isfolder(CONFIG.output_dir)
    mkdir(CONFIG.output_dir);
end

% Mirror the whole console session to the log file (overwrite any prior run).
log_path = CONFIG.log_file;
if isempty(fileparts(log_path))
    log_path = fullfile(CONFIG.output_dir, log_path);
end
diary off;
if isfile(log_path), delete(log_path); end
diary(log_path);

ct_lo   = min(CONFIG.ct_pair);
ct_hi   = max(CONFIG.ct_pair);
ct1_str = sprintf('CT_%d', ct_lo);
ct3_str = sprintf('CT_%d', ct_hi);

%% ===================== LOAD PRECOMPUTED STEP 2.5 SUMMARY ================
% Everything the summary figures need is precomputed by step25_segment_metrics
% and stored in segment_metrics_summary_<hash>.mat beside the recon doses.

[summary_path, hash_used] = find_summary_file(CONFIG.output_dir, CONFIG.config_hash);
CONFIG.config_hash = hash_used;   % pin the resolved hash for the random panels
fprintf('[STEP 2.5] Loading precomputed metrics: %s\n', summary_path);
SM = load(summary_path);

req  = {'comparisons', 'beams', 'n_segments', 'gamma', 'ssim'};
miss = req(~isfield(SM, req));
if ~isempty(miss)
    error('study_pass_rates_allsegments:BadSummary', ...
        'Summary %s missing field(s): %s', summary_path, strjoin(miss, ', '));
end

% Validate the settings that must match what Step 2.5 computed. We cannot
% recompute here, so a mismatch is either an error or adopt-the-summary + warn.
if isfield(SM, 'normalize') && ~isequal(logical(SM.normalize), logical(CONFIG.normalize))
    warning('study_pass_rates_allsegments:NormalizeMismatch', ...
        ['CONFIG.normalize=%d but the summary was computed with normalize=%d; ' ...
         'using the summary as-is.'], logical(CONFIG.normalize), logical(SM.normalize));
    CONFIG.normalize = logical(SM.normalize);
end
if isfield(SM, 'plan_type') && ~strcmpi(char(SM.plan_type), char(CONFIG.plan_type))
    warning('study_pass_rates_allsegments:PlanTypeMismatch', ...
        'CONFIG.plan_type=%s but the summary is %s; using the summary.', ...
        char(CONFIG.plan_type), char(SM.plan_type));
    CONFIG.plan_type = char(SM.plan_type);
end

% Criterion: the gamma metric is labelled with whatever Step 2.5 used.
crit = CONFIG.gamma_n(1);
if strcmp(CONFIG.eval_method, 'gamma_index') ...
        && isfield(SM, 'gamma_dose_pct') && ~isempty(SM.gamma_dose_pct)
    if ~isequal(SM.gamma_dose_pct, crit)
        warning('study_pass_rates_allsegments:CritMismatch', ...
            ['CONFIG.gamma_n=%g but Step 2.5 used %g%%/%g mm; ' ...
             'labelling with the summary value.'], ...
            crit, SM.gamma_dose_pct, SM.gamma_dist_mm);
    end
    crit = SM.gamma_dose_pct;
end

% Metric-specific labels (axis text, chart-title word, console header) resolved
% once the criterion is known, then threaded to every plotting / printing site.
metric = build_metric_descriptor(CONFIG.eval_method, crit);

% Column indices INTO THE SUMMARY for the comparisons the summary figures need.
d_tt = comp_col_names(SM.comparisons, 'truth1_vs_truth3');   % fig 1, blue
d_r1 = comp_col_names(SM.comparisons, 'truth1_vs_recon1');   % fig 1 green, fig 2 green
d_r3 = comp_col_names(SM.comparisons, 'truth1_vs_recon3');   % fig 1, red
d_33 = comp_col_names(SM.comparisons, 'truth3_vs_recon3');   % fig 2, red
if any(isnan([d_tt, d_r1, d_r3, d_33]))
    error('study_pass_rates_allsegments:MissingComparison', ...
        ['The Step-2.5 summary must define truth1_vs_truth3, truth1_vs_recon1, ' ...
         'truth1_vs_recon3 and truth3_vs_recon3 for the summary figures.']);
end

% Pick the metric block (per-beam mean/std + per-segment matrices) that matches
% CONFIG.eval_method. Both live in the same summary.
if metric.is_ssim
    sm_mean = SM.ssim.mean;        sm_std = SM.ssim.std;        sm_seg = SM.ssim.seg_mean;
else
    sm_mean = SM.gamma.mean_pass;  sm_std = SM.gamma.std_pass;  sm_seg = SM.gamma.seg_pass;
end

% Restrict to the requested beams (all must be present in the summary).
[tf, loc] = ismember(CONFIG.beams(:)', SM.beams(:)');
if ~all(tf)
    error('study_pass_rates_allsegments:BeamsNotInSummary', ...
        'Requested beam(s) %s are not in the Step-2.5 summary. Re-run Step 2.5 or edit CONFIG.beams.', ...
        mat2str(CONFIG.beams(~tf)));
end
beam_list        = SM.beams(loc);
mean_pass        = sm_mean(loc, :);
std_pass         = sm_std(loc, :);
nseg_used        = SM.n_segments(loc);
seg_pass_by_beam = sm_seg(loc);      % {1 x nBeam}, each [nSeg x nComp]
nComp            = numel(SM.comparisons);
nproc            = numel(beam_list);

% Noise-only null floor (gamma metric only): taken straight from the summary.
noise_floor = [];
if CONFIG.include_noise_floor && ~metric.is_ssim && isfield(SM, 'noise_floor')
    noise_floor = SM.noise_floor;
end

fprintf('============================================================\n');
fprintf(' STUDY_PASS_RATES_ALLSEGMENTS (plots Step-2.5 results)\n');
fprintf(' Patient %s | %s | plan=%s | hash=%s\n', ...
    CONFIG.patient_id, CONFIG.session, CONFIG.plan_type, CONFIG.config_hash);
fprintf(' Beams: %s  (%d)\n', mat2str(beam_list), nproc);
fprintf(' Metric: %s%s | normalize=%d\n', ...
    metric.name, metric.crit_suffix, CONFIG.normalize);
fprintf(' Summary : %s\n', summary_path);
fprintf(' Log file: %s\n', log_path);
fprintf('============================================================\n');

% Sort so the axis is well-defined regardless of how CONFIG.beams was written.
[beam_list, ord]  = sort(beam_list);
mean_pass         = mean_pass(ord, :);
std_pass          = std_pass(ord, :);
nseg_used         = nseg_used(ord);
seg_pass_by_beam  = seg_pass_by_beam(ord);

%% ===================== POOLED "ALL BEAMS" AGGREGATE =====================
%  Mean +/- std over every plotted segment (each segment weighted equally),
%  recomputed from the selected beams so the trailing "All" entry reflects
%  exactly what is drawn (rather than the summary's all-beam pool).

all_seg_pass = vertcat(seg_pass_by_beam{:});   % [sum(nseg_used) x nComp]
all_mean     = mean(all_seg_pass, 1, 'omitnan');
all_std      = std(all_seg_pass, 0, 1, 'omitnan');
all_nseg     = size(all_seg_pass, 1);

%% ===================== CONSOLE SUMMARY (BY BEAM) =======================

comp_names = SM.comparisons(:)';
fprintf('\n==================== BEAM %s (mean over segments) ====================\n', ...
    upper(metric.title));
fprintf('(%s, mean +/- std over segments)\n', metric.header_detail);
for n = 1:nproc
    fprintf('\n----- [beam #%d]  (%d segments) -----\n', beam_list(n), nseg_used(n));
    for d = 1:nComp
        fprintf('  %-18s   %6.2f%% +/- %5.2f%%\n', ...
            comp_names{d}, mean_pass(n, d), std_pass(n, d));
    end
end

fprintf('\n----- [ALL BEAMS]  (%d segments over %d beams) -----\n', all_nseg, nproc);
for d = 1:nComp
    fprintf('  %-18s   %6.2f%% +/- %5.2f%%\n', comp_names{d}, all_mean(d), all_std(d));
end
fprintf('\n=============================================================\n');

%% ===================== SUMMARY PLOTS ====================================
%  x positions are 1..nproc for the beams plus one trailing slot for "All", so
%  the pooled entry sits on the same axis without colliding with a beam number.

x_pos    = 1:(nproc + 1);
x_labels = [arrayfun(@(b) sprintf('%d', b), beam_list, 'UniformOutput', false), {'All'}];

% Figure 1: tabbed change-detection figure.
%   Tab 1 "Pass rates"    - both recons vs the CT_1 truth plus truth-vs-truth,
%                           with the noise floor drawn as a horizontal band.
%   Tab 2 "Differentials" - per-beam recon1-recon3 (green if >0, else red) and
%                           the true change 100-(truth1 vs truth3) (blue).
plot_change_detection_tabs(x_pos, x_labels, mean_pass, std_pass, ...
    all_mean, all_std, seg_pass_by_beam, all_seg_pass, ...
    d_r1, d_r3, d_tt, noise_floor, metric, CONFIG.patient_id, CONFIG.session);

% Figure 2: reconstruction fidelity (each recon against its OWN-CT truth).
plot_recon_fidelity_summary(x_pos, x_labels, ...
    [mean_pass(:, d_r1); all_mean(d_r1)], [std_pass(:, d_r1); all_std(d_r1)], ...
    [mean_pass(:, d_33); all_mean(d_33)], [std_pass(:, d_33); all_std(d_33)], ...
    metric, CONFIG.patient_id, CONFIG.session);

%% ============= RANDOM PER-SEGMENT PANEL VISUALIZATION ==================
% N random segments, each in its own tab: one row of truth | recon | metric per
% requested comparison. The metric map is the Step-2.5 folded map re-expanded
% from disk; only the recon/truth volumes are loaded here (no metric recompute).

if CONFIG.plot_random_segments
    show_random_segment_panels(CONFIG, crit, ct1_str, ct3_str, metric);
end

%% ========================= WRAP UP =====================================

total_runtime = toc(run_timer);
fprintf('\nTotal runtime: %.1f s (%.2f min) | %d beam(s), %d segment(s).\n', ...
    total_runtime, total_runtime/60, nproc, all_nseg);
fprintf('Console log written to: %s\n', log_path);
diary off;


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function name = add_name_suffix(name, suffix)
%ADD_NAME_SUFFIX Insert suffix before the extension of a filename (keeps any dir).
%  e.g. add_name_suffix('results.mat','_ssim') -> 'results_ssim.mat'.
    [d, base, ext] = fileparts(char(name));
    name = fullfile(d, [base, suffix, ext]);
    if isempty(d)
        name = [base, suffix, ext];   % keep it relative when no directory was given
    end
end

function [summary_path, hash] = find_summary_file(output_dir, config_hash)
%FIND_SUMMARY_FILE Locate the Step-2.5 rollup segment_metrics_summary_<hash>.mat.
%  With an explicit config_hash the matching file is required. Otherwise the
%  directory is scanned; a single match is used and its hash returned, while zero
%  or several matches raise (listing what was found) so the caller can disambiguate.
    if ~isempty(config_hash)
        hash = char(config_hash);
        summary_path = fullfile(output_dir, ...
            sprintf('segment_metrics_summary_%s.mat', hash));
        if exist(summary_path, 'file') ~= 2
            error('study_pass_rates_allsegments:NoSummary', ...
                ['Step-2.5 summary not found:\n  %s\n' ...
                 'Run step25_segment_metrics first (or clear CONFIG.config_hash to auto-discover).'], ...
                summary_path);
        end
        return;
    end

    d = dir(fullfile(output_dir, 'segment_metrics_summary_*.mat'));
    if isempty(d)
        error('study_pass_rates_allsegments:NoSummary', ...
            ['No segment_metrics_summary_*.mat in\n  %s\n' ...
             'Run step25_segment_metrics first.'], output_dir);
    end
    if numel(d) > 1
        hashes = regexprep({d.name}, '^segment_metrics_summary_(.+)\.mat$', '$1');
        error('study_pass_rates_allsegments:AmbiguousSummary', ...
            'Multiple Step-2.5 summaries in %s (hashes: %s). Set CONFIG.config_hash.', ...
            output_dir, strjoin(hashes, ', '));
    end
    summary_path = fullfile(output_dir, d(1).name);
    tok  = regexp(d(1).name, '^segment_metrics_summary_(.+)\.mat$', 'tokens', 'once');
    hash = tok{1};
end

function metric = build_metric_descriptor(eval_method, crit)
%BUILD_METRIC_DESCRIPTOR Labels/titles for the active evaluation metric.
%  Fields: .name (legend/console word), .title (chart-title word), .axis_label
%  (y-axis text), .crit_suffix (header criterion text), .header_detail (console
%  sub-header), .is_ssim (logical).
    switch eval_method
        case 'ssim'
            metric.is_ssim     = true;
            metric.name        = 'Mean local SSIM';
            metric.title       = 'SSIM';
            metric.axis_label  = 'Mean local SSIM (%)';
            metric.crit_suffix = '';
            metric.header_detail = 'mean local SSIM over the 10% region';
        otherwise   % 'gamma_index'
            metric.is_ssim     = false;
            metric.name        = 'Gamma pass rate';
            metric.title       = 'Pass-Rate';
            metric.axis_label  = sprintf('Gamma pass rate (%%)  @ %g%%/%g mm', crit, crit);
            metric.crit_suffix = sprintf(' @ %g%%/%g mm', crit, crit);
            metric.header_detail = sprintf('global gamma %g%%/%gmm', crit, crit);
    end
end

function d = comp_col_names(names, want)
%COMP_COL_NAMES Column index of comparison `want` within a cell vector of names.
    d = NaN;
    for r = 1:numel(names)
        if strcmpi(names{r}, want)
            d = r;
            return;
        end
    end
end

function d = comp_col(comparisons, name)
%COMP_COL Column index (in a {name,...} cell matrix) for a comparison name.
    d = NaN;
    for r = 1:size(comparisons, 1)
        if strcmpi(comparisons{r, 1}, name)
            d = r;
            return;
        end
    end
end

%% =========================================================================
%  CHANGE-DETECTION + FIDELITY PLOTS
%% =========================================================================

function plot_change_detection_tabs(x_pos, x_labels, mean_pass, std_pass, ...
        all_mean, all_std, seg_pass_by_beam, all_seg_pass, ...
        d_r1, d_r3, d_tt, noise_floor, metric, patient_id, session)
%PLOT_CHANGE_DETECTION_TABS One figure, two tabs on the shared beam x-axis:
%    Tab 1 "Pass rates"    - the per-beam pass-rate summary (recon1/recon3/truth-
%                            truth) with the noise-only null floor as a band.
%    Tab 2 "Differentials" - the per-beam recon1-recon3 and 100-(truth1_vs_truth3)
%                            differentials.
%  The final x slot on both tabs is the pooled "All" entry.
    fig = figure('Name', 'Change Detection (all segments)', 'Color', 'w', ...
        'NumberTitle', 'off', ...
        'Position', [100, 100, max(760, 55 * numel(x_pos) + 240), 500]);
    tg = uitabgroup(fig);

    tab1 = uitab(tg, 'Title', 'Pass rates');
    ax1  = axes('Parent', tab1); %#ok<LAXES>
    render_beam_pass_rate(ax1, x_pos, x_labels, ...
        [mean_pass(:, d_r1); all_mean(d_r1)], [std_pass(:, d_r1); all_std(d_r1)], ...
        [mean_pass(:, d_r3); all_mean(d_r3)], [std_pass(:, d_r3); all_std(d_r3)], ...
        [mean_pass(:, d_tt); all_mean(d_tt)], [std_pass(:, d_tt); all_std(d_tt)], ...
        noise_floor, metric, patient_id, session);

    tab2 = uitab(tg, 'Title', 'Differentials');
    ax2  = axes('Parent', tab2); %#ok<LAXES>
    [rd_m, rd_s, td_m, td_s] = beam_differentials(seg_pass_by_beam, all_seg_pass, ...
        d_r1, d_r3, d_tt);
    m_r1_all = [mean_pass(:, d_r1); all_mean(d_r1)];   % CT_1 pass rate per x slot
    render_differentials(ax2, x_pos, x_labels, rd_m, rd_s, td_m, td_s, ...
        m_r1_all, noise_floor, metric, patient_id, session);

    drawnow;
end

function render_beam_pass_rate(ax, x_pos, x_labels, m_r1, s_r1, m_r3, s_r3, ...
        m_tt, s_tt, noise_floor, metric, patient_id, session)
%RENDER_BEAM_PASS_RATE Per-beam mean metric (%) with std error bars, into axes ax.
%  Three series (recon1 green, recon3 red, truth-truth blue), a per-beam connector
%  between the two recon means (colour = greater one), and -- when noise_floor is
%  supplied -- the noise-only null as a horizontal mean +/- std band behind them.
    green = [0.15, 0.60, 0.20];
    red   = [0.80, 0.15, 0.15];
    blue  = [0.20, 0.40, 0.80];
    grey  = [0.35, 0.35, 0.35];

    x_pos = x_pos(:)';
    m_r1 = m_r1(:)'; s_r1 = s_r1(:)';
    m_r3 = m_r3(:)'; s_r3 = s_r3(:)';
    m_tt = m_tt(:)'; s_tt = s_tt(:)';

    axes(ax); hold(ax, 'on');

    % Noise-floor band (drawn first so it sits behind the markers) plus explicit
    % noise pass-rate data points at every beam slot. The noise-only null is a
    % single session-level value, so the points sit at a constant level across beams.
    h_nf = [];
    h_np = [];
    if ~isempty(noise_floor)
        nf_m = noise_floor.mean_pass_rate;
        nf_s = noise_floor.std_pass_rate;
        xl   = [min(x_pos) - 0.5, max(x_pos) + 0.5];
        h_nf = patch(ax, [xl(1), xl(2), xl(2), xl(1)], ...
            [nf_m - nf_s, nf_m - nf_s, nf_m + nf_s, nf_m + nf_s], grey, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot(ax, xl, [nf_m, nf_m], ':', 'Color', grey, 'LineWidth', 1.3);
        h_np = plot(ax, x_pos, nf_m * ones(size(x_pos)), 'd', 'Color', grey, ...
            'MarkerFaceColor', grey, 'MarkerSize', 6, 'LineStyle', 'none');
    end

    % Per-beam connector between the two recon means (colour = greater one).
    for n = 1:numel(x_pos)
        a = m_r1(n); b = m_r3(n);
        if isfinite(a) && isfinite(b)
            if a >= b, lc = green; else, lc = red; end
            plot(ax, [x_pos(n), x_pos(n)], [a, b], '-', 'Color', lc, 'LineWidth', 1.5);
        end
    end

    h_r1 = errorbar(ax, x_pos, m_r1, s_r1, 'o', 'Color', green, 'MarkerFaceColor', green, ...
        'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);
    h_r3 = errorbar(ax, x_pos, m_r3, s_r3, 'o', 'Color', red, 'MarkerFaceColor', red, ...
        'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);
    h_tt = errorbar(ax, x_pos, m_tt, s_tt, 'o', 'Color', blue, 'MarkerFaceColor', blue, ...
        'MarkerSize', 7, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);

    yline(90, 'k--', '90%', 'LineWidth', 1.0, 'FontSize', 8, ...
        'LabelHorizontalAlignment', 'left');
    apply_beam_axis(x_pos, x_labels);

    hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on');
    ylim(ax, [0, 105]);
    xlabel(ax, 'Beam number');
    ylabel(ax, metric.axis_label);
    title(ax, sprintf('Beam %s Summary (mean \\pm std over segments)   |   %s / %s', ...
        metric.title, strrep(patient_id, '_', '\_'), strrep(session, '_', '\_')), ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'tex');

    leg_h = [h_r1, h_r3, h_tt];
    leg_s = {'Recon CT\_1 vs Truth CT\_1', 'Recon CT\_3 vs Truth CT\_1', ...
             'Truth CT\_1 vs Truth CT\_3'};
    if ~isempty(h_np)
        leg_h(end + 1) = h_np;
        if isfield(noise_floor, 'is_fallback') && noise_floor.is_fallback
            leg_s{end + 1} = sprintf('Noise pass rate %.1f%% (assumed)', ...
                noise_floor.mean_pass_rate);
        else
            leg_s{end + 1} = sprintf('Noise pass rate %.1f \\pm %.1f%%', ...
                noise_floor.mean_pass_rate, noise_floor.std_pass_rate);
        end
    end
    legend(leg_h, leg_s, 'Location', 'best', 'FontSize', 9);
end

function render_differentials(ax, x_pos, x_labels, rd_m, rd_s, td_m, td_s, ...
        m_r1, noise_floor, metric, patient_id, session)
%RENDER_DIFFERENTIALS Per-beam pass-rate differentials into axes ax.
%  Recon differential (recon1 - recon3): green where >= 0 (change detected the
%  correct way), red where < 0. Truth differential 100 - (truth1 vs truth3), the
%  true change magnitude referenced to truth CT_1, in blue. Both carry std error
%  bars; the trailing slot is the pooled "All" entry. When noise_floor is supplied,
%  a noise-reconstruction differential (m_r1 CT_1 pass rate minus the session noise
%  pass rate) is drawn as one grey data point per slot.
    green = [0.15, 0.60, 0.20];
    red   = [0.80, 0.15, 0.15];
    blue  = [0.20, 0.40, 0.80];
    grey  = [0.35, 0.35, 0.35];

    x_pos = x_pos(:)';
    rd_m = rd_m(:)'; rd_s = rd_s(:)';
    td_m = td_m(:)'; td_s = td_s(:)';
    m_r1 = m_r1(:)';

    axes(ax); hold(ax, 'on');
    yline(0, 'k-', 'LineWidth', 1.0);

    % Recon differential: colour each point by the sign of its mean.
    for i = 1:numel(x_pos)
        if ~isfinite(rd_m(i)), continue; end
        if rd_m(i) >= 0, c = green; else, c = red; end
        errorbar(ax, x_pos(i), rd_m(i), rd_s(i), 'o', 'Color', c, 'MarkerFaceColor', c, ...
            'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);
    end

    % Truth differential (blue).
    h_td = errorbar(ax, x_pos, td_m, td_s, 's', 'Color', blue, 'MarkerFaceColor', blue, ...
        'MarkerSize', 7, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);

    % Noise-reconstruction differential: how far the real CT_1 recon sits above the
    % noise-only null (CT_1 pass rate minus the session noise pass rate). One grey
    % point per beam slot (plus the pooled "All"); the noise level is a single value.
    h_nd = [];
    if ~isempty(noise_floor)
        nd   = m_r1 - noise_floor.mean_pass_rate;
        h_nd = plot(ax, x_pos, nd, 'd', 'Color', grey, 'MarkerFaceColor', grey, ...
            'MarkerSize', 7, 'LineStyle', 'none');
    end

    % Legend proxies for the two-colour recon differential.
    h_pos = plot(ax, nan, nan, 'o', 'Color', green, 'MarkerFaceColor', green, ...
        'LineStyle', 'none', 'MarkerSize', 8);
    h_neg = plot(ax, nan, nan, 'o', 'Color', red, 'MarkerFaceColor', red, ...
        'LineStyle', 'none', 'MarkerSize', 8);

    apply_beam_axis(x_pos, x_labels);
    hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Beam number');
    ylabel(ax, sprintf('%s differential (%%)', metric.title));
    title(ax, sprintf(['Change-Detection Differentials (mean \\pm std over segments)' ...
        '   |   %s / %s'], strrep(patient_id, '_', '\_'), strrep(session, '_', '\_')), ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'tex');
    leg_h = [h_pos, h_neg, h_td];
    leg_s = {'Recon CT\_1 - CT\_3 > 0 (correct)', 'Recon CT\_1 - CT\_3 < 0', ...
             '100 - (Truth CT\_1 vs CT\_3)'};
    if ~isempty(h_nd)
        leg_h(end + 1) = h_nd;
        leg_s{end + 1} = 'Recon CT\_1 - Noise';
    end
    legend(leg_h, leg_s, 'Location', 'best', 'FontSize', 9);
end

function [rd_m, rd_s, td_m, td_s] = beam_differentials(seg_pass_by_beam, ...
        all_seg_pass, d_r1, d_r3, d_tt)
%BEAM_DIFFERENTIALS Per-beam mean/std of the paired differentials, plus a trailing
%  pooled "All" entry. Differentials are formed per SEGMENT and then reduced, so
%  the std is the spread of the per-segment difference (a paired statistic), not a
%  quadrature sum of two series' stds.
    nB   = numel(seg_pass_by_beam);
    rd_m = nan(1, nB + 1); rd_s = nan(1, nB + 1);
    td_m = nan(1, nB + 1); td_s = nan(1, nB + 1);
    for n = 1:nB
        [rd_m(n), rd_s(n), td_m(n), td_s(n)] = diff_stats(seg_pass_by_beam{n}, d_r1, d_r3, d_tt);
    end
    [rd_m(end), rd_s(end), td_m(end), td_s(end)] = diff_stats(all_seg_pass, d_r1, d_r3, d_tt);
end

function [rd_m, rd_s, td_m, td_s] = diff_stats(M, d_r1, d_r3, d_tt)
%DIFF_STATS Mean/std of the recon and truth differentials over a [nSeg x nComp] set.
%  recon diff = (recon1 vs truth1) - (recon3 vs truth1); truth diff = 100 - (truth1
%  vs truth3). Empty input yields NaNs.
    if isempty(M)
        rd_m = NaN; rd_s = NaN; td_m = NaN; td_s = NaN;
        return;
    end
    recon_d = M(:, d_r1) - M(:, d_r3);   % paired, per segment
    truth_d = 100 - M(:, d_tt);
    rd_m = mean(recon_d, 'omitnan'); rd_s = std(recon_d, 0, 'omitnan');
    td_m = mean(truth_d, 'omitnan'); td_s = std(truth_d, 0, 'omitnan');
end

function plot_recon_fidelity_summary(x_pos, x_labels, m_11, s_11, m_33, s_33, ...
        metric, patient_id, session)
%PLOT_RECON_FIDELITY_SUMMARY Own-CT reconstruction fidelity per beam.
%  Each recon is gamma-compared against the truth of its OWN CT - truth_CT1 vs
%  recon_CT1 (green) and truth_CT3 vs recon_CT3 (red) - so a series sitting
%  consistently higher means that CT's reconstruction is consistently better.
%  Points within a series are joined by straight lines; the two series are NOT
%  connected to each other. Same axis convention as the summary figure: the final
%  slot is the pooled "All" entry, fenced off by a dashed vertical rule.
    green = [0.15, 0.60, 0.20];
    red   = [0.80, 0.15, 0.15];

    x_pos = x_pos(:)';
    m_11 = m_11(:)'; s_11 = s_11(:)';
    m_33 = m_33(:)'; s_33 = s_33(:)';

    figure('Name', 'Reconstruction Fidelity Summary (all segments)', 'Color', 'w', ...
        'NumberTitle', 'off', ...
        'Position', [140, 140, max(720, 55 * numel(x_pos) + 220), 480]);
    hold on;

    % A NaN gap before the final slot keeps the connecting line within the beam
    % sequence instead of running it into the pooled "All" entry.
    [xg, m11g, s11g] = gap_before_last(x_pos, m_11, s_11);
    [~,  m33g, s33g] = gap_before_last(x_pos, m_33, s_33);

    h_11 = errorbar(xg, m11g, s11g, '-o', 'Color', green, 'MarkerFaceColor', green, ...
        'MarkerSize', 8, 'CapSize', 7, 'LineWidth', 1.2);
    h_33 = errorbar(xg, m33g, s33g, '-o', 'Color', red, 'MarkerFaceColor', red, ...
        'MarkerSize', 8, 'CapSize', 7, 'LineWidth', 1.2);

    yline(90, 'k--', '90%', 'LineWidth', 1.0, 'FontSize', 8, ...
        'LabelHorizontalAlignment', 'left');
    apply_beam_axis(x_pos, x_labels);

    hold off; grid on; box on;
    ylim([0, 105]);
    xlabel('Beam number');
    ylabel(metric.axis_label);
    title(sprintf(['Reconstruction Fidelity vs Own-CT Truth (mean \\pm std over ' ...
        'segments)   |   %s / %s'], ...
        strrep(patient_id, '_', '\_'), strrep(session, '_', '\_')), ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'tex');
    legend([h_11, h_33], ...
        {'Truth CT\_1 vs Recon CT\_1', 'Truth CT\_3 vs Recon CT\_3'}, ...
        'Location', 'best', 'FontSize', 9);
    drawnow;
end

function [xg, yg, eg] = gap_before_last(x, y, e)
%GAP_BEFORE_LAST Insert a NaN sample before the last point to break the line.
%  The NaN carries no marker and no error bar, so the last point still plots but
%  is not joined to the preceding series.
    if numel(x) < 2
        xg = x; yg = y; eg = e;
        return;
    end
    xg = [x(1:end-1), NaN, x(end)];
    yg = [y(1:end-1), NaN, y(end)];
    eg = [e(1:end-1), NaN, e(end)];
end

function apply_beam_axis(x_pos, x_labels)
%APPLY_BEAM_AXIS Tick labels for the beam slots + a rule before the "All" slot.
    if isempty(x_pos), return; end
    xlim([min(x_pos) - 0.5, max(x_pos) + 0.5]);
    xticks(x_pos);
    xticklabels(x_labels);
    if numel(x_pos) > 1
        xline(x_pos(end) - 0.5, ':', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.0);
    end
end

%% =========================================================================
%  RANDOM PER-SEGMENT PANEL VISUALIZATION (own cache; re-expands folded maps)
%% =========================================================================

function show_random_segment_panels(CONFIG, crit, ct1_str, ct3_str, metric)
%SHOW_RANDOM_SEGMENT_PANELS N random segments as tabs; one row of truth/recon/metric
%  panels per requested comparison. The metric map is the Step-2.5 folded gamma-index
%  ('gamma_index') or local-SSIM ('ssim') map re-expanded from the segment's CT_1
%  recon file. Replots from CONFIG.random_cache_file when it matches; otherwise
%  (re)loads recon volumes from disk, re-expands the folded maps, caches the
%  selected slices, and plots.
    fprintf('\n----------- RANDOM SEGMENT PANELS -----------\n');

    % Resolve the viz-cache path (own file, beside the recon doses).
    cache_path = CONFIG.random_cache_file;
    if isempty(fileparts(cache_path))
        cache_path = fullfile(CONFIG.output_dir, cache_path);
    end

    picks = [];
    if CONFIG.use_cache && exist(cache_path, 'file') == 2
        [ok, cpicks, why] = try_load_random_cache(cache_path, CONFIG, crit);
        if ok
            picks = cpicks;
            fprintf('[CACHE] Reusing %d random panel(s) from: %s\n', ...
                numel(picks), cache_path);
        else
            fprintf('[CACHE] Random-viz cache not reusable (%s); rebuilding.\n', why);
        end
    end

    if isempty(picks)
        picks = build_random_picks(CONFIG, ct1_str, ct3_str);
        if isempty(picks)
            warning('study_pass_rates_allsegments:NoRandomPicks', ...
                'No segments available for the random panel view; nothing to plot.');
            return;
        end
        if CONFIG.save_results
            meta = struct('config_hash', CONFIG.config_hash, ...
                'comparisons', {CONFIG.random_comparisons(:)'}, ...
                'n_random_requested', CONFIG.n_random, ...
                'seed', CONFIG.random_seed, 'normalize', logical(CONFIG.normalize), ...
                'crit', crit, 'eval_method', CONFIG.eval_method, ...
                'plan_type', CONFIG.plan_type, ...
                'patient_id', CONFIG.patient_id, 'session', CONFIG.session); %#ok<NASGU>
            try
                save(cache_path, 'picks', 'meta', '-v7.3');
                fprintf('[CACHE] Saved %d random panel(s) to: %s\n', ...
                    numel(picks), cache_path);
            catch ME
                warning('study_pass_rates_allsegments:RandomCacheSave', ...
                    'Could not save random-viz cache: %s', ME.message);
            end
        end
    end

    plot_random_picks(picks, CONFIG, crit, metric);
end

function [ok, picks, why] = try_load_random_cache(cache_path, CONFIG, crit)
%TRY_LOAD_RANDOM_CACHE Load cached random panels; validate against current CONFIG.
    ok = false; picks = []; why = '';
    try
        C = load(cache_path);
    catch ME
        why = sprintf('load error: %s', ME.message); return;
    end
    if ~isfield(C, 'picks') || ~isfield(C, 'meta')
        why = 'missing picks/meta'; return;
    end
    m = C.meta;
    if ~isequal(m.config_hash, CONFIG.config_hash), why = 'config_hash differs'; return; end
    if ~isequal(m.comparisons(:)', CONFIG.random_comparisons(:)')
        why = 'comparison set differs'; return;
    end
    if ~isequal(m.n_random_requested, CONFIG.n_random), why = 'N differs'; return; end
    if ~isequal(m.seed, CONFIG.random_seed), why = 'seed differs'; return; end
    if ~isequal(logical(m.normalize), logical(CONFIG.normalize)), why = 'normalize differs'; return; end
    if ~isequal(m.crit, crit), why = 'criterion differs'; return; end
    cached_method = 'gamma_index';
    if isfield(m, 'eval_method') && ~isempty(m.eval_method)
        cached_method = lower(char(m.eval_method));
    end
    if ~strcmp(cached_method, CONFIG.eval_method), why = 'eval_method differs'; return; end
    if isfield(m, 'plan_type') && ~strcmpi(char(m.plan_type), char(CONFIG.plan_type))
        why = 'plan_type differs'; return;
    end
    picks = C.picks;
    ok = true;
end

function picks = build_random_picks(CONFIG, ct1_str, ct3_str)
%BUILD_RANDOM_PICKS Randomly select CONFIG.n_random segments and build their
%  truth/recon/metric display slices for each CONFIG.random_comparisons entry.
%  Loads beams lazily (seeded shuffle) until enough segments are collected; the
%  metric map for each panel is the Step-2.5 folded map, not a recomputation.
    picks = struct('beam', {}, 'seg', {}, 'spacing', {}, 'panels', {});

    specs = resolve_comparisons(CONFIG.comparisons, CONFIG.random_comparisons);
    if isempty(specs)
        warning('study_pass_rates_allsegments:NoRandomComparisons', ...
            'None of CONFIG.random_comparisons matched CONFIG.comparisons.');
        return;
    end

    rng(CONFIG.random_seed);                 % reproducible selection
    beams = CONFIG.beams(:)';
    beams = beams(randperm(numel(beams)));   % shuffle beam visit order
    need  = CONFIG.n_random;

    for b = beams
        if numel(picks) >= need, break; end
        try
            out = load_beam_set(CONFIG, b);
        catch ME
            warning('study_pass_rates_allsegments:SkipBeam', ...
                'Random view: skipping beam #%d (load failed): %s', b, ME.message);
            continue;
        end
        spacing  = out.metadata.spacing(:)';
        seg_data = group_beam_segments(out.fields, spacing, ct1_str, ct3_str, ...
            CONFIG.normalize, b);
        if isempty(seg_data), continue; end

        seg_data = seg_data(randperm(numel(seg_data)));   % random segment order
        for k = 1:numel(seg_data)
            if numel(picks) >= need, break; end
            S = seg_data{k};
            pk = struct('beam', b, 'seg', S.seg, 'spacing', spacing, ...
                'panels', build_panels(S, specs, CONFIG.eval_method, CONFIG.config_hash));
            picks(end+1) = pk; %#ok<AGROW>
        end
    end

    fprintf('[RANDOM] Selected %d segment(s) (requested %d) across beams.\n', ...
        numel(picks), need);
end

function panels = build_panels(S, specs, eval_method, hash)
%BUILD_PANELS One panel entry per comparison spec: truth/recon slices, the metric-map
%  slice (at the reference max-dose slice), the stored score %, and labels. The
%  metric map is re-expanded from the segment's CT_1 recon file (Step-2.5 folded
%  gamma-index or local-SSIM values) rather than recomputed. .gamma_img holds
%  whichever map and .passpct the matching stored score.
    if nargin < 3 || isempty(eval_method), eval_method = 'gamma_index'; end
    panels = struct('name', {}, 'ref_label', {}, 'tgt_label', {}, ...
        'ref_img', {}, 'tgt_img', {}, 'gamma_img', {}, 'passpct', {}, ...
        'dose_max', {}, 'iz', {});
    for c = 1:size(specs, 1)
        ref = get_dose_field(S, specs{c, 2});
        tgt = get_dose_field(S, specs{c, 3});

        % Metric map + stored score from the Step-2.5 fold (no recompute). When the
        % fold is absent/mismatched, show the images with a blank (NaN) metric map.
        [gmap, passpct] = expand_folded_map(S.ct1_recon_file, specs{c, 1}, ...
            eval_method, hash);
        if isempty(gmap) || ~isequal(size(gmap), size(ref))
            warning('study_pass_rates_allsegments:NoFoldedMetric', ...
                ['Beam #%d seg %d: no Step-2.5 folded %s map for "%s"; ' ...
                 'showing a blank metric panel (run/rerun step25_segment_metrics).'], ...
                S.beam, S.seg, eval_method, specs{c, 1});
            gmap    = nan(size(ref));
            passpct = NaN;
        end

        % Display slice: axial slice (dim 3) of the reference's peak dose.
        [~, iz] = max_slice(ref);

        P = struct();
        P.name      = specs{c, 1};
        P.ref_label = specs{c, 4};
        P.tgt_label = specs{c, 5};
        P.ref_img   = squeeze(ref(:, :, iz));
        P.tgt_img   = squeeze(tgt(:, :, iz));
        P.gamma_img = squeeze(gmap(:, :, iz));
        P.passpct   = passpct;
        P.dose_max  = max(ref(:));
        P.iz        = iz;
        panels(c) = P; %#ok<AGROW>
    end
end

function [M, passpct] = expand_folded_map(ct1_file, comp_name, eval_method, hash)
%EXPAND_FOLDED_MAP Re-expand a Step-2.5 folded metric map to a full volume.
%  Reads the 'segment_metrics' variable folded into a segment's CT_1 recon .mat,
%  finds the named comparison, and rebuilds the dense gamma-index / local-SSIM
%  volume (NaN outside the stored 10% eval mask), returning the stored scalar
%  score alongside it. Returns [] / NaN when the file, the variable, the matching
%  config hash, or the named comparison is not present (caller falls back to a
%  blank panel -- nothing is recomputed here).
    M = []; passpct = NaN;
    if isempty(ct1_file) || exist(ct1_file, 'file') ~= 2
        return;
    end
    try
        L = load(ct1_file, 'segment_metrics');
    catch
        return;
    end
    if ~isfield(L, 'segment_metrics'), return; end
    sm = L.segment_metrics;

    % Only use a fold produced for this config hash.
    if ~isempty(hash) && isfield(sm, 'config_hash') && ~isempty(sm.config_hash) ...
            && ~strcmpi(char(sm.config_hash), char(hash))
        return;
    end
    if ~isfield(sm, 'comparison') || ~isfield(sm, 'vol_size'), return; end

    c = [];
    for k = 1:numel(sm.comparison)
        if strcmpi(sm.comparison(k).name, comp_name)
            c = sm.comparison(k);
            break;
        end
    end
    if isempty(c), return; end

    M = nan(sm.vol_size);
    if strcmp(eval_method, 'ssim')
        if ~isempty(c.ssim_vals), M(c.mask_idx) = double(c.ssim_vals); end
        passpct = c.ssim_mean;
    else
        if ~isempty(c.gamma_vals), M(c.mask_idx) = double(c.gamma_vals); end
        passpct = c.gamma_pass_rate;
    end
end

function [pk, iz] = max_slice(vol)
%MAX_SLICE Index (dim 3) of the slice holding the volume's maximum value.
    [pk, lin] = max(vol(:));
    [~, ~, iz] = ind2sub(size(vol), lin);
    if isempty(iz) || ~isfinite(pk), iz = max(1, round(size(vol, 3) / 2)); end
end

function specs = resolve_comparisons(all_comps, names)
%RESOLVE_COMPARISONS Rows of all_comps (name,ref,tgt,ref_label,tgt_label) for names.
    specs = cell(0, 5);
    for i = 1:numel(names)
        r = comp_col(all_comps, names{i});
        if isnan(r)
            warning('study_pass_rates_allsegments:UnknownRandomComparison', ...
                'Random comparison "%s" not found in CONFIG.comparisons; skipping.', names{i});
            continue;
        end
        specs(end+1, :) = all_comps(r, 1:5); %#ok<AGROW>
    end
end

function plot_random_picks(picks, CONFIG, crit, metric)
%PLOT_RANDOM_PICKS One window, one tab per pick; each tab has one ROW of
%  truth | recon | metric panels per comparison (first comparison on top). The
%  metric column is a gamma-index map ('gamma_index') or a local-SSIM map ('ssim').
%  With the default CONFIG.random_comparisons that puts recon_CT1 (Dose 1 on recon
%  1) on the top row and recon_CT3 (Dose 2 on recon 1) on the bottom row.
    np = numel(picks);
    if np == 0, return; end
    if nargin < 4 || isempty(metric)
        metric = build_metric_descriptor(CONFIG.eval_method, crit);
    end

    % Metric-column appearance (colormap / colour limits / title + tab suffix).
    if metric.is_ssim
        metric_cmap = parula(256);
        metric_clim = [0, 1];               % local SSIM, 1 = identical
        tab_suffix  = 'local SSIM';
    else
        metric_cmap = gamma_colormap(256);
        metric_clim = [0, 2];               % gamma, <=1 passes
        tab_suffix  = sprintf('gamma %g%%/%g mm', crit, crit);
    end

    nc_max = max(arrayfun(@(k) numel(picks(k).panels), 1:np));
    fig = figure('Name', sprintf('Random segment panels  |  %s / %s', ...
        CONFIG.patient_id, CONFIG.session), 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 900, max(300, 300 * nc_max)]);
    tg = uitabgroup(fig);

    for i = 1:np
        pk = picks(i);
        P  = pk.panels;
        nc = numel(P);

        tab = uitab(tg, 'Title', sprintf('B%d S%d', pk.beam, pk.seg));
        % nc rows (one comparison each) x 3 columns (truth | recon | metric). Tiles
        % fill row-major, so comparison c occupies row c.
        tl  = tiledlayout(tab, nc, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl, sprintf('Beam %d, Segment %d   |   %s', ...
            pk.beam, pk.seg, tab_suffix), 'FontWeight', 'bold');

        for c = 1:nc
            p = P(c);
            cmax = p.dose_max; if ~(cmax > 0), cmax = 1; end

            % Truth dose (reference).
            ax = nexttile(tl);
            imagesc(ax, p.ref_img); axis(ax, 'image', 'off');
            colormap(ax, hot); clim(ax, [0, cmax]); colorbar(ax);
            title(ax, {p.ref_label, 'truth'}, 'Interpreter', 'tex');

            % Reconstruction (target).
            ax = nexttile(tl);
            imagesc(ax, p.tgt_img); axis(ax, 'image', 'off');
            colormap(ax, hot); clim(ax, [0, cmax]); colorbar(ax);
            title(ax, {p.tgt_label, 'recon'}, 'Interpreter', 'tex');

            % Metric map (gamma index <=1 passes, or local SSIM 0..1).
            ax = nexttile(tl);
            imagesc(ax, p.gamma_img); axis(ax, 'image', 'off');
            colormap(ax, metric_cmap); clim(ax, metric_clim); colorbar(ax);
            if metric.is_ssim
                title(ax, sprintf('SSIM map  (mean %.1f%%)', p.passpct), ...
                    'Interpreter', 'tex');
            else
                title(ax, sprintf('\\gamma index  (%.1f%% \\leq 1)', p.passpct), ...
                    'Interpreter', 'tex');
            end
        end
    end
    drawnow;
end

function cmap = gamma_colormap(n)
%GAMMA_COLORMAP Green for gamma<=1 (pass), ramping green->red for 1<gamma<=2 (fail).
    if nargin < 1 || isempty(n), n = 256; end
    x    = linspace(0, 2, n)';
    cmap = zeros(n, 3);
    green = [0.15, 0.60, 0.20];
    red   = [0.80, 0.15, 0.15];
    for i = 1:n
        if x(i) <= 1
            cmap(i, :) = green;                 % solid pass colour
        else
            f = min(x(i) - 1, 1);               % 0..1 across gamma 1..2
            cmap(i, :) = (1 - f) * green + f * red;
        end
    end
end

%% =========================================================================
%  SEGMENT LOADING / GROUPING (recon volumes for the random panels only)
%% =========================================================================

function out = load_beam_set(CONFIG, beam)
%LOAD_BEAM_SET Load every segment/CT of a beam via load_recon_dose_data (set).
    args = {'Mode', 'set', 'Beam', beam, 'IncludeEthos', false, 'IncludeCBCT', false};
    if isfield(CONFIG, 'plan_type') && ~isempty(CONFIG.plan_type) ...
            && ~strcmpi(CONFIG.plan_type, 'any')
        args = [args, {'PlanType', CONFIG.plan_type}];
    end
    if isfield(CONFIG, 'config_hash') && ~isempty(CONFIG.config_hash)
        args = [args, {'Hash', CONFIG.config_hash}];
    end
    out = load_recon_dose_data(CONFIG.patient_id, CONFIG.session, CONFIG, args{:});
end

function seg_data = group_beam_segments(fields, spacing, ct1_str, ct3_str, normalize, beam)
%GROUP_BEAM_SEGMENTS Group a beam's loaded fields into per-segment volume sets.
%  Returns a cell array; each element S has .recon_CT1/.recon_CT3/.rs_CT1/.rs_CT3
%  (double), .spacing, .seg, .beam and .ct1_recon_file (the file carrying the
%  Step-2.5 folded metrics). Segments missing the CT_1/CT_3 pair, with mismatched
%  grids, or without a CT_1 recon file path are skipped (warning). When `normalize`
%  is true each recon is scaled by the least-squares gain to its OWN-CT truth over
%  that truth's 10% region -- the same scaling Step 2.5 applied.
    seg_ids = arrayfun(@(k) seg_of_field(fields(k)), 1:numel(fields));
    uniq    = unique(seg_ids(~isnan(seg_ids)));
    uniq    = uniq(:).';   % ensure row so the loop iterates one segment at a time

    seg_data = {};
    for s = uniq
        sub = fields(seg_ids == s);
        iA  = find_field_by_ct(sub, ct1_str);
        iB  = find_field_by_ct(sub, ct3_str);
        if isempty(iA) || isempty(iB)
            continue;   % missing CT_1/CT_3 pair for this segment -> skip entry
        end
        fA = sub(iA);
        fB = sub(iB);

        S = struct();
        S.recon_CT1      = double(fA.recon_dose);
        S.recon_CT3      = double(fB.recon_dose);
        S.rs_CT1         = double(fA.rs_dose);
        S.rs_CT3         = double(fB.rs_dose);
        S.spacing        = spacing;
        S.seg            = s;
        S.beam           = beam;
        S.ct1_recon_file = fA.recon_file;   % holds the Step-2.5 folded metrics

        % Skip entries whose recon/truth grids disagree (bad / missing data).
        if ~isequal(size(S.recon_CT1), size(S.rs_CT1)) || ...
           ~isequal(size(S.recon_CT3), size(S.rs_CT3)) || ...
           ~isequal(size(S.rs_CT1),    size(S.rs_CT3))
            warning('study_pass_rates_allsegments:GridMismatch', ...
                'Beam #%d seg %d: grid mismatch; skipping segment.', beam, s);
            continue;
        end

        % Least-squares relative normalization (recon -> own-CT truth).
        if normalize
            S.recon_CT1 = S.recon_CT1 * least_squares_gain(S.rs_CT1, S.recon_CT1);
            S.recon_CT3 = S.recon_CT3 * least_squares_gain(S.rs_CT3, S.recon_CT3);
        end

        seg_data{end+1} = S; %#ok<AGROW>
    end
end

function s = seg_of_field(fld)
%SEG_OF_FIELD Segment number of a loaded field (rtplan, else filename token).
    s = NaN;
    if isfield(fld, 'rtplan') && isfield(fld.rtplan, 'seg_num') ...
            && ~isempty(fld.rtplan.seg_num) && isnumeric(fld.rtplan.seg_num) ...
            && ~isnan(fld.rtplan.seg_num)
        s = double(fld.rtplan.seg_num);
        return;
    end
    if isfield(fld, 'source_mat_filename')
        tok = regexp(char(fld.source_mat_filename), '_B\d+_(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok), s = str2double(tok{1}); end
    end
end

function idx = find_field_by_ct(fields, want_ct)
%FIND_FIELD_BY_CT Index of the first field whose CT label is want_ct.
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
        tok = regexp(char(fld.source_mat_filename), 'CT[_-]?(\d+)', 'tokens', 'once');
        if ~isempty(tok), lbl = sprintf('CT_%s', tok{1}); end
    end
    if isempty(lbl), lbl = 'unknown'; end
end

function g = least_squares_gain(rs_truth, recon)
%LEAST_SQUARES_GAIN Scalar gain aligning a recon to its RS truth (relative norm).
%  g = sum(rs.*recon)/sum(recon.^2) over the truth's 10% low-dose region; the g
%  minimizing ||rs - g*recon||^2 there. Falls back to 1 for empty/zero inputs.
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

function v = get_dose_field(S, fieldname)
%GET_DOSE_FIELD Fetch a named per-segment volume for a comparison spec.
    switch fieldname
        case 'rs_CT1',    v = S.rs_CT1;
        case 'rs_CT3',    v = S.rs_CT3;
        case 'recon_CT1', v = S.recon_CT1;
        case 'recon_CT3', v = S.recon_CT3;
        otherwise
            error('study_pass_rates_allsegments:BadField', ...
                'Unknown comparison volume "%s".', fieldname);
    end
end
