%% =========================================================================
%  STUDY_PASS_RATES_ALLSEGMENTS.m
%  Per-BEAM photoacoustic gamma summary, averaged over ALL segments in the beam.
%
%  A companion to study_pass_rates_individual.m. Instead of one field per beam,
%  this loads EVERY segment of each beam (both CT images) from the already-
%  reconstructed doses on disk, computes the gamma pass rate of each configured
%  comparison for every segment, then reduces to a per-beam MEAN +/- STD across
%  segments. The only output is the beam pass-rate summary figure (mean markers
%  with std error bars); all dose/sensor/gamma-index/curve plotting is omitted.
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
                'Beam #%d seg %d: grid mismatch; skipping segment.', b, s);
            continue;
        end

        % Least-squares relative normalization (recon -> own-CT truth).
        if CONFIG.normalize
            S.recon_CT1 = S.recon_CT1 * least_squares_gain(S.rs_CT1, S.recon_CT1);
            S.recon_CT3 = S.recon_CT3 * least_squares_gain(S.rs_CT3, S.recon_CT3);
        end

        seg_data{end+1} = S; %#ok<SAGROW>
    end

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

% Beams are processed in ascending order, but sort so the axis is well-defined
% regardless of how CONFIG.beams was written.
[beam_list, ord]  = sort(beam_list);
mean_pass         = mean_pass(ord, :);
std_pass          = std_pass(ord, :);
nseg_used         = nseg_used(ord);
seg_pass_by_beam  = seg_pass_by_beam(ord);

%% ===================== POOLED "ALL BEAMS" AGGREGATE =====================
%  Mean +/- std over every segment of every processed beam (each segment weighted
%  equally), plotted as the trailing "All" x-axis entry on both figures.

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

%% ========================= SAVE RESULTS ================================

total_runtime = toc(run_timer);

if CONFIG.save_results
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

    out_path = CONFIG.output_file;
    if isempty(fileparts(out_path))
        out_path = fullfile(CONFIG.output_dir, out_path);
    end
    save(out_path, '-struct', 'RESULTS', '-v7.3');
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
