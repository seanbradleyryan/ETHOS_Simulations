%% =========================================================================
%  STUDY_PASS_RATES_ALLSEGMENTS.m
%  Per-BEAM photoacoustic gamma summary, averaged over ALL segments in the beam.
%
%  A companion to study_pass_rates_individual.m. Instead of one field per beam,
%  this loads EVERY segment of each beam (both CT images) from the already-
%  reconstructed doses on disk, computes the gamma pass rate of each configured
%  comparison for every segment, then reduces to a per-beam MEAN +/- STD across
%  segments. Outputs are two per-beam summary figures (change detection and
%  reconstruction fidelity) plus an OPTIONAL window of N random per-segment panels
%  (truth | recon | gamma index). Both the pass-rate results and the random panels
%  are cached so a matching re-run replots without recomputing gamma.
%
%  Comparisons (CONFIG.comparisons), each per segment:
%    truth1_vs_truth3  RayStation truth CT_1 vs truth CT_3   (green connector ref)
%    truth1_vs_recon1  RayStation truth CT_1 vs recon CT_1   (green)
%    truth1_vs_recon3  RayStation truth CT_1 vs recon CT_3   (red)
%    truth3_vs_recon3  RayStation truth CT_3 vs recon CT_3   (own-CT fidelity)
%
%  When CONFIG.normalize is true, each recon is first rescaled by the least-
%  squares gain that best matches it to its OWN-CT truth over that truth's 10%
%  region (recon_CT1->rs_CT1, recon_CT3->rs_CT3), cancelling the stored absolute
%  scale (CONFIG.correction_factor); the RS truths are never rescaled.
%
%  FIGURE 1 - change detection (x = beam number, y = pass rate %, single crit):
%    - recon_CT1 vs truth_CT1  (green markers, std error bars),
%    - recon_CT3 vs truth_CT1  (red markers, std error bars),
%    - truth_CT1 vs truth_CT3  (blue markers, std error bars),
%    - a per-beam vertical connector between the two recon means, coloured like
%      whichever mean is greater.
%
%  FIGURE 2 - reconstruction fidelity (each recon vs its OWN-CT truth):
%    - truth_CT1 vs recon_CT1 (green), truth_CT3 vs recon_CT3 (red), each series
%      joined by straight lines across beams. No vertical connectors.
%
%  Both figures carry a trailing "All" x-axis entry: the mean +/- std pooled over
%  EVERY segment of EVERY processed beam, separated from the per-beam entries by
%  a dashed vertical rule.
%
%  OPERATIONAL:
%    - Progress notifications + a total-runtime record are printed to console.
%    - The whole console session is mirrored to CONFIG.log_file via diary, and
%      results + log are written beside the recon doses in
%      SimulationResults/[PatientID]/[Session]/[method]/.
%    - ALL CalcGamma console output is suppressed (evalc wrapper).
%    - Missing entries (absent beam, absent segment, missing CT_1/CT_3 pair or
%      recon file, grid mismatch) are skipped rather than raised.
%    - Per-beam segment gamma runs in parallel (parfor) or serially per
%      CONFIG.use_parallel.
%
%  NOTE: HIPAA / remote-execution - this file is WRITTEN here but must be RUN on
%  the remote device. Do not execute locally.
%  =========================================================================

clear; clc; close all;
run_timer = tic;   % program runtime record

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';
CONFIG.treatment_site = 'Pancreas';

% Beams to summarize (all segments of each are loaded and averaged).
CONFIG.beams = 1:17;

% Restrict to a single plan type so segment grouping is unambiguous.
CONFIG.plan_type = 'reference';   % 'reference' | 'adapted' | 'any'

% The two CT image indices (lower -> *_CT1 volumes, higher -> *_CT3 volumes).
CONFIG.ct_pair = [1, 3];

% Explicit recon config-hash override ('' => auto-discover on disk via loader).
CONFIG.config_hash = '505ae853';

% Needed by load_recon_dose_data to resolve the method folder / hash on disk.
CONFIG.gruneisen_method = 'threshold_2';

% --- Comparisons to run (configurable list) ---------------------------------
%  {name, ref_field, tgt_field, ref_label, tgt_label}. ref/tgt name a per-segment
%  volume: 'rs_CT1' | 'rs_CT3' | 'recon_CT1' | 'recon_CT3'. Gamma uses ref_field
%  as the reference and builds the 10% eval mask from it.
CONFIG.comparisons = { ...
    'truth1_vs_truth3', 'rs_CT1', 'rs_CT3',    'Truth CT\_1', 'Truth CT\_3'; ...
    'truth1_vs_recon1', 'rs_CT1', 'recon_CT1', 'Truth CT\_1', 'Recon CT\_1'; ...
    'truth1_vs_recon3', 'rs_CT1', 'recon_CT3', 'Truth CT\_1', 'Recon CT\_3'; ...
    'truth3_vs_recon3', 'rs_CT3', 'recon_CT3', 'Truth CT\_3', 'Recon CT\_3'  ...
};

% Gamma criterion n (evaluated as n%/n mm). Single value for the summary.
CONFIG.gamma_n = 3;   % 3%/3 mm

% Least-squares relative normalization of each recon to its own-CT truth.
CONFIG.normalize = true;

% Parallelization: true -> parfor over a beam's segments; false -> serial.
% CPU gamma (see quiet_gamma_pass), so workers scale with physical cores.
CONFIG.use_parallel = true;
CONFIG.num_workers  = 64;   % default local pool size

% --- Output ---
% Results and the console log are cached beside the recon doses, i.e. in
% SimulationResults/[PatientID]/[Session]/[gruneisen_method]/. Set output_dir to
% a path to override that location.
CONFIG.save_results = true;
CONFIG.output_dir   = '';   % '' => the recon-dose directory
CONFIG.output_file  = 'pass_rates_allsegments_results.mat';
CONFIG.log_file     = 'pass_rates_allsegments_log.txt';

% --- Cache ---
% When true, if CONFIG.output_file already exists AND its saved run matches the
% current CONFIG (hash, criterion, comparisons, beams, normalize), the per-beam
% gamma is SKIPPED and the summary plot / console table are regenerated straight
% from the cached results. Set false to always recompute.
CONFIG.use_cache = true;

% --- Random per-segment panel visualization ---
% When true, ALSO show N randomly chosen segments, each in its own tab of a single
% window, with one row of panels per tab: truth dose | reconstruction | gamma index
% for every comparison in CONFIG.random_comparisons. The selected panels are cached
% to CONFIG.random_cache_file (keyed by hash/criterion/comparisons/N/seed/normalize)
% so a re-run with the same settings replots WITHOUT recomputing gamma; on a miss
% the recon volumes are (re)loaded from disk and the gamma maps recomputed. This
% cache is INDEPENDENT of the pass-rate results cache above (which stores no volumes).
CONFIG.plot_random_segments = true;
CONFIG.n_random             = 5;     % number of random segments to display
CONFIG.random_seed          = 42;    % reproducible selection (used in the cache key)
% Which comparisons to draw per segment (names from CONFIG.comparisons). Default:
% the two recon comparisons -> "Dose 1 on recon 1" and "Dose 2 on recon 1".
CONFIG.random_comparisons   = {'truth1_vs_recon1', 'truth1_vs_recon3'};
CONFIG.random_cache_file    = 'pass_rates_allsegments_random_viz.mat';

%% ===================== SETUP ============================================

% Cache directory: the recon-dose folder unless explicitly overridden.
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

beams   = CONFIG.beams(:)';
nBeams  = numel(beams);
nComp   = size(CONFIG.comparisons, 1);
crit    = CONFIG.gamma_n(1);

has_gamma = (exist('CalcGamma', 'file') == 2);
if ~has_gamma
    error('study_pass_rates_allsegments:NoCalcGamma', ...
        'CalcGamma not found on the path; cannot compute pass rates.');
end

% Column indices (into CONFIG.comparisons / the pass matrices) for the plots.
d_tt = comp_col(CONFIG.comparisons, 'truth1_vs_truth3');   % fig 1, blue
d_r1 = comp_col(CONFIG.comparisons, 'truth1_vs_recon1');   % fig 1 green, fig 2 green
d_r3 = comp_col(CONFIG.comparisons, 'truth1_vs_recon3');   % fig 1, red
d_33 = comp_col(CONFIG.comparisons, 'truth3_vs_recon3');   % fig 2, red

if any(isnan([d_tt, d_r1, d_r3, d_33]))
    error('study_pass_rates_allsegments:MissingComparison', ...
        ['CONFIG.comparisons must define truth1_vs_truth3, truth1_vs_recon1, ' ...
         'truth1_vs_recon3 and truth3_vs_recon3 for the summary figures.']);
end

fprintf('============================================================\n');
fprintf(' STUDY_PASS_RATES_ALLSEGMENTS\n');
fprintf(' Patient %s | %s | plan=%s | hash=%s\n', ...
    CONFIG.patient_id, CONFIG.session, CONFIG.plan_type, CONFIG.config_hash);
fprintf(' Beams: %s  (%d)\n', mat2str(beams), nBeams);
fprintf(' Criterion: %g%%/%gmm | normalize=%d | parallel=%d\n', ...
    crit, crit, CONFIG.normalize, CONFIG.use_parallel);
fprintf(' Cache dir: %s\n', CONFIG.output_dir);
fprintf(' Log file : %s\n', log_path);
fprintf('============================================================\n');

% Resolve the results/cache path once (used for both cache-read and save). Results
% live beside the recon doses (CONFIG.output_dir) unless an absolute name is given.
out_path = CONFIG.output_file;
if isempty(fileparts(out_path))
    out_path = fullfile(CONFIG.output_dir, out_path);
end

%% ===================== CACHE CHECK =====================================
% If a matching cached run exists, load it and skip straight to plotting.

loaded_from_cache = false;
beam_list = [];
mean_pass = [];
std_pass  = [];
nseg_used = [];
nproc     = 0;
total_seg_evals = 0;

if CONFIG.use_cache && exist(out_path, 'file') == 2
    [ok, cached, why] = try_load_cache(out_path, CONFIG, crit);
    if ok
        beam_list = cached.beams;
        mean_pass = cached.mean_pass;
        std_pass  = cached.std_pass;
        nseg_used = cached.n_segments;
        seg_pass_by_beam = cached.seg_pass;   % {1 x nBeams}, each [nSeg x nComp]
        nproc     = numel(beam_list);
        total_seg_evals = sum(nseg_used) * nComp;
        loaded_from_cache = true;
        fprintf('\n[CACHE] Reusing results from: %s\n', out_path);
        fprintf('[CACHE] %d beam(s), %d segment(s); skipping gamma recomputation.\n', ...
            nproc, sum(nseg_used));
    else
        fprintf('\n[CACHE] Existing file not reusable (%s); recomputing.\n', why);
    end
end

if ~loaded_from_cache

if CONFIG.use_parallel
    ensure_pool(CONFIG.num_workers);
end

% Per-beam aggregates (trimmed to processed beams afterwards).
beam_list = nan(1, nBeams);
mean_pass = nan(nBeams, nComp);
std_pass  = nan(nBeams, nComp);
nseg_used = zeros(1, nBeams);
nproc     = 0;
total_seg_evals = 0;

% Per-segment pass rates kept per beam ([nSeg x nComp] each) so the pooled "All"
% entry is a true mean/std over every segment rather than a mean of beam means.
seg_pass_by_beam = cell(1, nBeams);

%% ===================== PER-BEAM PROCESSING ==============================

for bi = 1:nBeams
    b = beams(bi);
    fprintf('\n[Beam %d/%d] beam #%d: loading all segments...\n', bi, nBeams, b);

    % ---- Load every segment of this beam (both CTs). Skip on any load error. ----
    try
        out = load_beam_set(CONFIG, b);
    catch ME
        warning('study_pass_rates_allsegments:SkipBeam', ...
            'Skipping beam #%d (load failed): %s', b, ME.message);
        continue;
    end

    fields  = out.fields;
    spacing = out.metadata.spacing(:)';

    % ---- Group fields by segment; build a normalized volume set per segment. ----
    seg_data = group_beam_segments(fields, spacing, ct1_str, ct3_str, ...
        CONFIG.normalize, b);

    ns = numel(seg_data);
    if ns == 0
        warning('study_pass_rates_allsegments:NoSegments', ...
            'Beam #%d: no usable segments found; skipping beam.', b);
        continue;
    end

    % ---- Gamma pass rate per (segment x comparison). CalcGamma is silenced. ----
    comparisons = CONFIG.comparisons;
    pass_seg    = nan(ns, nComp);
    if CONFIG.use_parallel
        parfor si = 1:ns
            pass_seg(si, :) = seg_pass_rates(seg_data{si}, comparisons, crit);
        end
    else
        for si = 1:ns
            pass_seg(si, :) = seg_pass_rates(seg_data{si}, comparisons, crit);
        end
    end

    % ---- Reduce to per-beam mean +/- std across segments. ----
    nproc = nproc + 1;
    beam_list(nproc)   = b;
    mean_pass(nproc,:) = mean(pass_seg, 1, 'omitnan');
    std_pass(nproc,:)  = std(pass_seg, 0, 1, 'omitnan');
    nseg_used(nproc)   = ns;
    seg_pass_by_beam{nproc} = pass_seg;
    total_seg_evals    = total_seg_evals + ns * nComp;

    elapsed = toc(run_timer);
    eta     = elapsed / bi * (nBeams - bi);
    fprintf(['  beam #%d done: %d segment(s) used | mean pass ' ...
             'recon1=%.1f%% recon3=%.1f%% truth3=%.1f%% | elapsed %.1fs, ETA %.1fs\n'], ...
        b, ns, mean_pass(nproc, d_r1), mean_pass(nproc, d_r3), ...
        mean_pass(nproc, d_tt), elapsed, eta);
end

% Trim aggregates to processed beams.
beam_list = beam_list(1:nproc);
mean_pass = mean_pass(1:nproc, :);
std_pass  = std_pass(1:nproc, :);
nseg_used = nseg_used(1:nproc);
seg_pass_by_beam = seg_pass_by_beam(1:nproc);

if nproc == 0
    error('study_pass_rates_allsegments:NoBeams', ...
        'No beams could be processed (all skipped).');
end

end   % if ~loaded_from_cache

% Sort so the axis is well-defined regardless of how CONFIG.beams was written.
% Applies to both the freshly-computed and the cache-loaded paths.
[beam_list, ord]  = sort(beam_list);
mean_pass         = mean_pass(ord, :);
std_pass          = std_pass(ord, :);
nseg_used         = nseg_used(ord);
seg_pass_by_beam  = seg_pass_by_beam(ord);

%% ===================== POOLED "ALL BEAMS" AGGREGATE =====================
%  Mean +/- std over every segment of every processed beam (each segment weighted
%  equally), plotted as the trailing "All" x-axis entry on both figures. Computed
%  for both the freshly-computed and the cache-loaded paths (seg_pass_by_beam is
%  restored from the cache in the latter).

all_seg_pass = vertcat(seg_pass_by_beam{:});   % [sum(nseg_used) x nComp]
all_mean     = mean(all_seg_pass, 1, 'omitnan');
all_std      = std(all_seg_pass, 0, 1, 'omitnan');
all_nseg     = size(all_seg_pass, 1);

%% ===================== CONSOLE SUMMARY (BY BEAM) =======================

fprintf('\n==================== BEAM MEAN PASS RATES ====================\n');
fprintf('(global gamma %g%%/%gmm, mean +/- std over segments)\n', crit, crit);
for n = 1:nproc
    fprintf('\n----- [beam #%d]  (%d segments) -----\n', beam_list(n), nseg_used(n));
    for d = 1:nComp
        fprintf('  %-18s (%s vs %s)   %6.2f%% +/- %5.2f%%\n', ...
            CONFIG.comparisons{d,1}, CONFIG.comparisons{d,2}, CONFIG.comparisons{d,3}, ...
            mean_pass(n, d), std_pass(n, d));
    end
end

fprintf('\n----- [ALL BEAMS]  (%d segments over %d beams) -----\n', all_nseg, nproc);
for d = 1:nComp
    fprintf('  %-18s (%s vs %s)   %6.2f%% +/- %5.2f%%\n', ...
        CONFIG.comparisons{d,1}, CONFIG.comparisons{d,2}, CONFIG.comparisons{d,3}, ...
        all_mean(d), all_std(d));
end
fprintf('\n=============================================================\n');

%% ===================== SUMMARY PLOTS ====================================
%  x positions are 1..nproc for the beams plus one trailing slot for "All", so
%  the pooled entry sits on the same axis without colliding with a beam number.

x_pos    = 1:(nproc + 1);
x_labels = [arrayfun(@(b) sprintf('%d', b), beam_list, 'UniformOutput', false), {'All'}];

% Figure 1: change detection (both recons vs the CT_1 truth, plus truth-vs-truth).
plot_beam_pass_rate_summary(x_pos, x_labels, ...
    [mean_pass(:, d_r1); all_mean(d_r1)], [std_pass(:, d_r1); all_std(d_r1)], ...
    [mean_pass(:, d_r3); all_mean(d_r3)], [std_pass(:, d_r3); all_std(d_r3)], ...
    [mean_pass(:, d_tt); all_mean(d_tt)], [std_pass(:, d_tt); all_std(d_tt)], ...
    crit, CONFIG.patient_id, CONFIG.session);

% Figure 2: reconstruction fidelity (each recon against its OWN-CT truth).
plot_recon_fidelity_summary(x_pos, x_labels, ...
    [mean_pass(:, d_r1); all_mean(d_r1)], [std_pass(:, d_r1); all_std(d_r1)], ...
    [mean_pass(:, d_33); all_mean(d_33)], [std_pass(:, d_33); all_std(d_33)], ...
    crit, CONFIG.patient_id, CONFIG.session);

%% ============= RANDOM PER-SEGMENT PANEL VISUALIZATION ==================
% N random segments, each in its own tab: one row of truth | recon | gamma per
% requested comparison. Uses its own cache (recomputes gamma only on a miss).

if CONFIG.plot_random_segments
    show_random_segment_panels(CONFIG, crit, ct1_str, ct3_str);
end

%% ========================= SAVE RESULTS ================================

total_runtime = toc(run_timer);

if CONFIG.save_results && ~loaded_from_cache
    RESULTS = struct();
    RESULTS.config        = CONFIG;
    RESULTS.gamma_crit    = crit;
    RESULTS.beams         = beam_list;
    RESULTS.comparisons   = CONFIG.comparisons(:,1)';
    RESULTS.mean_pass     = mean_pass;   % [nBeams x nComp]
    RESULTS.std_pass      = std_pass;    % [nBeams x nComp]
    RESULTS.n_segments    = nseg_used;
    RESULTS.seg_pass      = seg_pass_by_beam;   % {1 x nBeams}, each [nSeg x nComp]
    RESULTS.all_mean      = all_mean;    % [1 x nComp], pooled over every segment
    RESULTS.all_std       = all_std;     % [1 x nComp]
    RESULTS.all_n_segments = all_nseg;
    RESULTS.log_file      = log_path;
    RESULTS.total_runtime_s = total_runtime;

    save(out_path, '-struct', 'RESULTS', '-v7.3');   % out_path resolved above
    fprintf('\nResults saved to: %s\n', out_path);
end

fprintf('\nTotal runtime: %.1f s (%.2f min) | %d beam(s), %d segment-evaluation(s).\n', ...
    total_runtime, total_runtime/60, nproc, total_seg_evals);
fprintf('Console log written to: %s\n', log_path);
diary off;


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function ensure_pool(desired_workers)
%ENSURE_POOL Start a local CPU pool of desired_workers if none exists.
%  Clamps the request to the cluster's NumWorkers limit; non-fatal if the
%  Parallel Computing Toolbox is absent (loops then run serially).
    if nargin < 1 || isempty(desired_workers), desired_workers = 64; end
    if exist('parpool', 'file') ~= 2
        fprintf('  [WARN] Parallel Computing Toolbox not found; loops run serially.\n');
        return;
    end
    try
        if isempty(gcp('nocreate'))
            nw = desired_workers;
            try
                c  = parcluster('local');
                nw = min(desired_workers, c.NumWorkers);
                if nw < desired_workers
                    fprintf('  [Pool] Requested %d workers; cluster caps at %d.\n', ...
                        desired_workers, nw);
                end
            catch
            end
            parpool('local', nw);
            fprintf('  [Pool] Started local pool with %d worker(s).\n', nw);
        end
    catch ME
        fprintf('  [WARN] Could not start parallel pool (%s). Running serially.\n', ME.message);
    end
end

function [ok, cached, why] = try_load_cache(out_path, CONFIG, crit)
%TRY_LOAD_CACHE Load a saved results .mat and check it matches the current run.
%  Returns ok=true and a struct with .beams/.mean_pass/.std_pass/.n_segments/
%  .seg_pass (restricted to the requested beams, in order) when the cached run was
%  produced with the same gamma criterion, comparison set, requested beams, config
%  hash and normalization; otherwise ok=false and a short reason in `why`. Any
%  load/format error is reported (non-fatal) rather than raised.
    ok = false; cached = struct(); why = '';
    try
        cached = load(out_path);
    catch ME
        why = sprintf('load error: %s', ME.message);
        return;
    end

    req = {'beams', 'mean_pass', 'std_pass', 'n_segments', 'seg_pass', ...
           'gamma_crit', 'comparisons', 'config'};
    miss = req(~isfield(cached, req));
    if ~isempty(miss)
        why = sprintf('missing field(s): %s', strjoin(miss, ', '));
        return;
    end

    % Gamma criterion must match.
    if ~isequal(cached.gamma_crit, crit)
        why = sprintf('criterion %g%% ~= cached %g%%', crit, cached.gamma_crit);
        return;
    end

    % Comparison set (names, in order) must match.
    if ~isequal(cached.comparisons(:)', CONFIG.comparisons(:,1)')
        why = 'comparison set differs';
        return;
    end

    % Every requested beam must be present in the cache.
    if ~all(ismember(CONFIG.beams(:)', cached.beams(:)'))
        why = 'requested beams not all present in cache';
        return;
    end

    % Key inputs that change the numbers must match the saved CONFIG.
    cc = cached.config;
    if isfield(cc, 'config_hash') && ~isequal(cc.config_hash, CONFIG.config_hash)
        why = 'config_hash differs';
        return;
    end
    if isfield(cc, 'normalize') && ~isequal(logical(cc.normalize), logical(CONFIG.normalize))
        why = 'normalize flag differs';
        return;
    end
    if isfield(cc, 'plan_type') && ~strcmpi(char(cc.plan_type), char(CONFIG.plan_type))
        why = 'plan_type differs';
        return;
    end

    % Restrict the cached aggregates to exactly the requested beams, in order.
    [tf, loc] = ismember(CONFIG.beams(:)', cached.beams(:)');
    loc = loc(tf);
    cached.beams      = cached.beams(loc);
    cached.mean_pass  = cached.mean_pass(loc, :);
    cached.std_pass   = cached.std_pass(loc, :);
    cached.n_segments = cached.n_segments(loc);
    cached.seg_pass   = cached.seg_pass(loc);   % {1 x nBeams} per-segment matrices

    ok = true;
end

function seg_data = group_beam_segments(fields, spacing, ct1_str, ct3_str, normalize, beam)
%GROUP_BEAM_SEGMENTS Group a beam's loaded fields into per-segment volume sets.
%  Returns a cell array; each element S has .recon_CT1/.recon_CT3/.rs_CT1/.rs_CT3
%  (double), .spacing and .seg. Segments missing the CT_1/CT_3 pair or with
%  mismatched grids are skipped (warning). When `normalize` is true each recon is
%  scaled by the least-squares gain to its OWN-CT truth over that truth's 10% region.
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
        S.recon_CT1 = double(fA.recon_dose);
        S.recon_CT3 = double(fB.recon_dose);
        S.rs_CT1    = double(fA.rs_dose);
        S.rs_CT3    = double(fB.rs_dose);
        S.spacing   = spacing;
        S.seg       = s;

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

function p = seg_pass_rates(S, comparisons, crit)
%SEG_PASS_RATES Pass rate (%) of every comparison for one segment.
%  Reference = comparison ref volume; 10% reference cutoff eval mask; global
%  gamma at crit%/crit mm. Returns a 1 x nComp row (NaN where gamma failed).
    np = size(comparisons, 1);
    p  = nan(1, np);
    for d = 1:np
        ref = get_dose_field(S, comparisons{d, 2});
        tgt = get_dose_field(S, comparisons{d, 3});
        if max(ref(:)) > 0
            mask = ref >= 0.10 * max(ref(:));
        else
            mask = tgt >= 0.10 * max(tgt(:));
        end
        p(d) = quiet_gamma_pass(ref, tgt, mask, crit, S.spacing);
    end
end

function p = quiet_gamma_pass(ref, tgt, mask, crit, spacing)
%QUIET_GAMMA_PASS Global gamma pass rate (%), with CalcGamma output suppressed.
%  evalc captures/discards all console output CalcGamma would otherwise print.
%  Forces CPU computation ('cpu', 1) so many parfor workers do not contend over
%  the GPU.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(ref));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(tgt));
    gmap = [];
    try
        evalc(['gmap = CalcGamma(ref_struct, tgt_struct, crit, crit, ', ...
               '''local'', 0, ''limit'', crit*2, ''restrict'', 1, ''cpu'', 1);']);
        p = 100 * mean(gmap(mask) <= 1);
    catch
        p = NaN;
    end
end

function d = comp_col(comparisons, name)
%COMP_COL Column index (in comparisons / the pass matrices) for a comparison name.
    d = NaN;
    for r = 1:size(comparisons, 1)
        if strcmpi(comparisons{r, 1}, name)
            d = r;
            return;
        end
    end
end

function plot_beam_pass_rate_summary(x_pos, x_labels, m_r1, s_r1, m_r3, s_r3, ...
        m_tt, s_tt, crit, patient_id, session)
%PLOT_BEAM_PASS_RATE_SUMMARY Mean pass rate (%) per beam, std error bars.
%  Three series at the given n%/n mm criterion:
%    recon_CT1 vs truth_CT1 (green), recon_CT3 vs truth_CT1 (red),
%    truth_CT1 vs truth_CT3 (blue), each with std error bars. A per-beam vertical
%    connector joins the two recon means, coloured like whichever mean is greater.
%  x_pos are consecutive slots and x_labels their tick labels; the final slot is
%  the pooled "All" entry and is fenced off by a dashed vertical rule.
    green = [0.15, 0.60, 0.20];
    red   = [0.80, 0.15, 0.15];
    blue  = [0.20, 0.40, 0.80];

    x_pos = x_pos(:)';
    m_r1 = m_r1(:)'; s_r1 = s_r1(:)';
    m_r3 = m_r3(:)'; s_r3 = s_r3(:)';
    m_tt = m_tt(:)'; s_tt = s_tt(:)';

    figure('Name', 'Beam Pass-Rate Summary (all segments)', 'Color', 'w', ...
        'NumberTitle', 'off', ...
        'Position', [100, 100, max(720, 55 * numel(x_pos) + 220), 480]);
    hold on;

    % Per-beam connector between the two recon means (colour = greater one).
    for n = 1:numel(x_pos)
        a = m_r1(n);
        b = m_r3(n);
        if isfinite(a) && isfinite(b)
            if a >= b, lc = green; else, lc = red; end
            plot([x_pos(n), x_pos(n)], [a, b], '-', 'Color', lc, 'LineWidth', 1.5);
        end
    end

    h_r1 = errorbar(x_pos, m_r1, s_r1, 'o', 'Color', green, 'MarkerFaceColor', green, ...
        'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);
    h_r3 = errorbar(x_pos, m_r3, s_r3, 'o', 'Color', red, 'MarkerFaceColor', red, ...
        'MarkerSize', 8, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);
    h_tt = errorbar(x_pos, m_tt, s_tt, 'o', 'Color', blue, 'MarkerFaceColor', blue, ...
        'MarkerSize', 7, 'LineStyle', 'none', 'CapSize', 7, 'LineWidth', 1.2);

    yline(90, 'k--', '90%', 'LineWidth', 1.0, 'FontSize', 8, ...
        'LabelHorizontalAlignment', 'left');
    apply_beam_axis(x_pos, x_labels);

    hold off; grid on; box on;
    ylim([0, 105]);
    xlabel('Beam number');
    ylabel(sprintf('Gamma pass rate (%%)  @ %g%%/%g mm', crit, crit));
    title(sprintf('Beam Pass-Rate Summary (mean \\pm std over segments)   |   %s / %s', ...
        strrep(patient_id, '_', '\_'), strrep(session, '_', '\_')), ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'tex');
    legend([h_r1, h_r3, h_tt], ...
        {'Recon CT\_1 vs Truth CT\_1', 'Recon CT\_3 vs Truth CT\_1', ...
         'Truth CT\_1 vs Truth CT\_3'}, 'Location', 'best', 'FontSize', 9);
    drawnow;
end

%% =========================================================================
%  RANDOM PER-SEGMENT PANEL VISUALIZATION (own cache)
%% =========================================================================

function show_random_segment_panels(CONFIG, crit, ct1_str, ct3_str)
%SHOW_RANDOM_SEGMENT_PANELS N random segments as tabs; one row of truth/recon/gamma
%  panels per requested comparison. Replots from CONFIG.random_cache_file when it
%  matches (no gamma recompute); otherwise (re)loads recon volumes from disk,
%  recomputes the gamma maps, caches the selected slices, and plots.
    fprintf('\n----------- RANDOM SEGMENT PANELS -----------\n');

    % Resolve the viz-cache path (own file, independent of the results cache).
    % Cached beside the recon doses (CONFIG.output_dir), like the results cache.
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
            fprintf('[CACHE] Random-viz cache not reusable (%s); recomputing.\n', why);
        end
    end

    if isempty(picks)
        picks = build_random_picks(CONFIG, crit, ct1_str, ct3_str);
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
                'crit', crit, 'plan_type', CONFIG.plan_type, ...
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

    plot_random_picks(picks, CONFIG, crit);
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
    if isfield(m, 'plan_type') && ~strcmpi(char(m.plan_type), char(CONFIG.plan_type))
        why = 'plan_type differs'; return;
    end
    picks = C.picks;
    ok = true;
end

function picks = build_random_picks(CONFIG, crit, ct1_str, ct3_str)
%BUILD_RANDOM_PICKS Randomly select CONFIG.n_random segments and compute their
%  truth/recon/gamma display slices for each CONFIG.random_comparisons entry.
%  Loads beams lazily (seeded shuffle) until enough segments are collected.
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
                'panels', build_panels(S, specs, crit));
            picks(end+1) = pk; %#ok<AGROW>
        end
    end

    fprintf('[RANDOM] Selected %d segment(s) (requested %d) across beams.\n', ...
        numel(picks), need);
end

function panels = build_panels(S, specs, crit)
%BUILD_PANELS One panel entry per comparison spec: truth/recon slices, gamma slice
%  (at the reference max-dose slice), pass %, and labels.
    panels = struct('name', {}, 'ref_label', {}, 'tgt_label', {}, ...
        'ref_img', {}, 'tgt_img', {}, 'gamma_img', {}, 'passpct', {}, ...
        'dose_max', {}, 'iz', {});
    for c = 1:size(specs, 1)
        ref = get_dose_field(S, specs{c, 2});
        tgt = get_dose_field(S, specs{c, 3});

        if max(ref(:)) > 0
            mask = ref >= 0.10 * max(ref(:));
        else
            mask = tgt >= 0.10 * max(tgt(:));
        end
        gmap = quiet_gamma_map(ref, tgt, crit, S.spacing);
        if isempty(gmap), gmap = nan(size(ref)); end
        if any(mask(:))
            passpct = 100 * mean(gmap(mask) <= 1, 'omitnan');
        else
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

function gmap = quiet_gamma_map(ref, tgt, crit, spacing)
%QUIET_GAMMA_MAP Full global gamma-index volume (CalcGamma output suppressed).
%  Same CalcGamma call as quiet_gamma_pass but returns the map instead of a scalar.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(ref));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(tgt));
    gmap = [];
    try
        evalc(['gmap = CalcGamma(ref_struct, tgt_struct, crit, crit, ', ...
               '''local'', 0, ''limit'', crit*2, ''restrict'', 1, ''cpu'', 1);']);
    catch
        gmap = [];
    end
end

function plot_random_picks(picks, CONFIG, crit)
%PLOT_RANDOM_PICKS One window, one tab per pick; each tab is a single row of
%  truth | recon | gamma panels for every comparison.
    np = numel(picks);
    if np == 0, return; end

    fig = figure('Name', sprintf('Random segment panels  |  %s / %s', ...
        CONFIG.patient_id, CONFIG.session), 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 1500, 380]);
    tg = uitabgroup(fig);
    gmap_cmap = gamma_colormap(256);

    for i = 1:np
        pk = picks(i);
        P  = pk.panels;
        nc = numel(P);

        tab = uitab(tg, 'Title', sprintf('B%d S%d', pk.beam, pk.seg));
        tl  = tiledlayout(tab, 1, 3 * nc, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl, sprintf('Beam %d, Segment %d   |   gamma %g%%/%g mm', ...
            pk.beam, pk.seg, crit, crit), 'FontWeight', 'bold');

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

            % Gamma index (<=1 passes).
            ax = nexttile(tl);
            imagesc(ax, p.gamma_img); axis(ax, 'image', 'off');
            colormap(ax, gmap_cmap); clim(ax, [0, 2]); colorbar(ax);
            title(ax, sprintf('\\gamma index  (%.1f%% \\leq 1)', p.passpct), ...
                'Interpreter', 'tex');
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

function plot_recon_fidelity_summary(x_pos, x_labels, m_11, s_11, m_33, s_33, ...
        crit, patient_id, session)
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
    ylabel(sprintf('Gamma pass rate (%%)  @ %g%%/%g mm', crit, crit));
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
