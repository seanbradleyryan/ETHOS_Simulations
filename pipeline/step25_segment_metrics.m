function results = step25_segment_metrics(patient_id, session, config)
%STEP25_SEGMENT_METRICS Per-segment gamma + SSIM metrics folded into recon files.
%
%   results = step25_segment_metrics(patient_id, session, config)
%
%   PURPOSE:
%   Adapt the batch analysis of study_pass_rates_allsegments.m into a pipeline
%   step that runs after Step 2 (k-Wave reconstruction). For every beam/segment
%   whose CT_1 and CT_3 reconstructions exist on disk, it computes the SAME four
%   comparisons the study plots -- but saves the RAW per-segment data instead of
%   figures, keyed to the simulation config hash. Both the global gamma index
%   (at CONFIG.gamma_dose_pct/gamma_dist_mm) and the local SSIM map are computed
%   for each comparison.
%
%   Why a separate post-Step-2 step (not inside the sim parfor): two of the four
%   comparisons are CROSS-CT (recon_CT3 vs truth_CT1, truth_CT1 vs truth_CT3),
%   so they need BOTH the CT_1 and CT_3 reconstructions of the same segment. The
%   sim parfor produces one field (one CT) at a time, so the pair is not yet
%   available there. This step pairs them and parallelizes across CPUs (parfor
%   over a beam's segments) while the sim's GPU work is already done.
%
%   COMPARISONS (per segment, reference builds the 10% eval mask):
%       truth1_vs_truth3   rs_CT1    vs rs_CT3
%       truth1_vs_recon1   rs_CT1    vs recon_CT1
%       truth1_vs_recon3   rs_CT1    vs recon_CT3
%       truth3_vs_recon3   rs_CT3    vs recon_CT3
%
%   NORMALIZATION (CONFIG.metrics_normalize, default true): each recon is scaled
%   by the least-squares gain that best matches it to its OWN-CT truth over that
%   truth's 10% region (recon_CT1->rs_CT1, recon_CT3->rs_CT3), exactly as the
%   study does. The RS truths are never rescaled. Gains are stored per segment.
%
%   STORAGE (masked-region-only, folded into the CT_1 recon file):
%   The per-segment result is appended as a variable 'segment_metrics' into that
%   segment's CT_1 reconstruction .mat (<base_CT1>_recon_<hash>.mat). Only the
%   voxels inside each comparison's 10% eval mask are stored (linear indices +
%   single-precision values + volume size), so the full gamma-index / SSIM map
%   can be re-expanded on demand without keeping hundreds of GB of dense volumes.
%   Nothing is written to the CT_3 recon file. See the context .md files.
%
%   RESUMABLE: a segment whose CT_1 recon file already carries 'segment_metrics'
%   is skipped (its stored scalars are still folded into the summary), so the
%   step can be re-run piecewise as more fields finish, just like the sim.
%   Set CONFIG.metrics_overwrite = true to force recomputation.
%
%   INPUTS:
%       patient_id - Char/string patient identifier (e.g. '1194203').
%       session    - Char/string session name (e.g. 'Session_1').
%       config     - Struct with (defaults filled if absent):
%           .working_dir           REQUIRED. Base directory path.
%           .gruneisen_method      Simulation-results subfolder (recon location).
%           .gamma_dose_pct        [3.0]  Gamma dose-difference criterion (%).
%           .gamma_dist_mm         [3.0]  Gamma distance-to-agreement (mm).
%           .gamma_dose_cutoff_pct [10.0] Low-dose eval cutoff (% of ref max).
%           .metrics_plan_type     ['reference'] 'reference'|'adapted'|'any'.
%           .metrics_ct_pair       [1 3]  The two CT indices (lo->CT1, hi->CT3).
%           .metrics_normalize     [true] Least-squares recon->own-truth gain.
%           .metrics_beams         []     Beams to process ([] => all on disk).
%           .metrics_overwrite     [false] Recompute even if already folded in.
%           .metrics_write_summary [true] Write the per-beam rollup .mat.
%           .use_parallel          [true] parfor over a beam's segments (CPU).
%           .num_parallel_workers  [8]    Pool size when none is running.
%
%   OUTPUTS:
%       results - Summary struct (also saved to segment_metrics_summary_<hash>.mat
%                 when metrics_write_summary): per-beam and pooled mean/std of the
%                 gamma pass rate (%) and mean local SSIM (%), per comparison,
%                 plus the per-segment matrices. This is the data the study's
%                 plots are built from -- no figures are produced here.
%
%   DEPENDENCIES:
%       - CalcGamma.m, compute_local_ssim.m (utils), load_recon_dose_data.m
%       - Pipeline outputs from pipeline_simulate (per-field recon files).
%
%   See also: study_pass_rates_allsegments, load_recon_dose_data,
%             compute_local_ssim, CalcGamma, step3_analysis, pipeline_simulate

    %% ======================== INPUT VALIDATION ========================
    if ~ischar(patient_id) && ~isstring(patient_id)
        error('step25_segment_metrics:InvalidInput', ...
            'patient_id must be a string or character array.');
    end
    patient_id = char(patient_id);
    if ~ischar(session) && ~isstring(session)
        error('step25_segment_metrics:InvalidInput', ...
            'session must be a string or character array.');
    end
    session = char(session);
    if ~isstruct(config) || ~isfield(config, 'working_dir')
        error('step25_segment_metrics:MissingConfig', ...
            'config must be a struct containing working_dir.');
    end

    %% ======================== CONFIG DEFAULTS ========================
    config = default_field(config, 'gamma_dose_pct',        3.0);
    config = default_field(config, 'gamma_dist_mm',         3.0);
    config = default_field(config, 'gamma_dose_cutoff_pct', 10.0);
    config = default_field(config, 'metrics_plan_type',     'reference');
    config = default_field(config, 'metrics_ct_pair',       [1, 3]);
    config = default_field(config, 'metrics_normalize',     true);
    config = default_field(config, 'metrics_beams',         []);
    config = default_field(config, 'metrics_overwrite',     false);
    config = default_field(config, 'metrics_write_summary', true);
    config = default_field(config, 'use_parallel',          true);
    config = default_field(config, 'num_parallel_workers',  8);

    if exist('CalcGamma', 'file') ~= 2
        error('step25_segment_metrics:NoCalcGamma', ...
            'CalcGamma not found on the path; cannot compute gamma pass rates.');
    end

    dose_pct   = config.gamma_dose_pct;
    dist_mm    = config.gamma_dist_mm;
    cutoff_fr  = config.gamma_dose_cutoff_pct / 100;
    normalize  = logical(config.metrics_normalize);
    ct_lo      = min(config.metrics_ct_pair);
    ct_hi      = max(config.metrics_ct_pair);
    ct1_str    = sprintf('CT_%d', ct_lo);
    ct3_str    = sprintf('CT_%d', ct_hi);

    % Fixed comparison set (name, reference volume, target volume). The
    % reference builds the 10% eval mask and is CalcGamma's reference.
    comparisons = { ...
        'truth1_vs_truth3', 'rs_CT1', 'rs_CT3'; ...
        'truth1_vs_recon1', 'rs_CT1', 'recon_CT1'; ...
        'truth1_vs_recon3', 'rs_CT1', 'recon_CT3'; ...
        'truth3_vs_recon3', 'rs_CT3', 'recon_CT3'  ...
    };
    nComp = size(comparisons, 1);

    fprintf('\n=========================================================\n');
    fprintf('  [STEP 2.5] Per-Segment Gamma + SSIM Metrics\n');
    fprintf('  Patient %s | %s | plan=%s\n', patient_id, session, config.metrics_plan_type);
    fprintf('  Gamma: %.1f%%/%.1f mm (cutoff %.1f%% of ref max) | normalize=%d\n', ...
        dose_pct, dist_mm, config.gamma_dose_cutoff_pct, normalize);
    fprintf('=========================================================\n');

    %% ======================== RESOLVE BEAM LIST ========================
    beams = config.metrics_beams(:)';
    if isempty(beams)
        [field_index, ~] = list_processed_field_doses(patient_id, session, config);
        beams = unique([field_index.beam_index]);
        beams = beams(~isnan(beams));
    end
    if isempty(beams)
        error('step25_segment_metrics:NoBeams', ...
            'No beams found to analyze for %s / %s.', patient_id, session);
    end
    nBeams = numel(beams);

    % CPU pool for the per-segment gamma/SSIM. Reuse any running pool; the sim's
    % pool is fine (just fewer workers than physical cores).
    if config.use_parallel
        ensure_cpu_pool(config.num_parallel_workers);
    end

    %% ======================== PER-BEAM PROCESSING ========================
    run_timer = tic;
    beam_list       = nan(1, nBeams);
    n_segments      = zeros(1, nBeams);
    seg_gamma_pass  = cell(1, nBeams);   % {beam} = [nSeg x nComp] gamma pass %
    seg_ssim_mean   = cell(1, nBeams);   % {beam} = [nSeg x nComp] mean SSIM %
    nproc           = 0;
    sim_dir         = '';
    hash8           = '';

    for bi = 1:nBeams
        b = beams(bi);
        fprintf('\n[Beam %d/%d] beam #%d: loading segments...\n', bi, nBeams, b);

        try
            beam_set = load_beam_set(patient_id, session, config, b);
        catch ME
            warning('step25_segment_metrics:SkipBeam', ...
                'Skipping beam #%d (load failed): %s', b, ME.message);
            continue;
        end

        % Capture the recon location / hash from the first successful load.
        if isempty(hash8)
            hash8   = beam_set.config_hash;
            sim_dir = fullfile(config.working_dir, 'SimulationResults', ...
                patient_id, session, beam_set.gruneisen_method);
        end

        spacing  = beam_set.metadata.spacing(:)';
        seg_data = group_beam_segments(beam_set.fields, spacing, ct1_str, ct3_str, ...
            normalize, b);
        ns = numel(seg_data);
        if ns == 0
            warning('step25_segment_metrics:NoSegments', ...
                'Beam #%d: no usable CT_%d/CT_%d segment pairs; skipping.', ...
                b, ct_lo, ct_hi);
            continue;
        end

        gamma_pass = nan(ns, nComp);
        ssim_mean  = nan(ns, nComp);
        overwrite  = config.metrics_overwrite;

        if config.use_parallel
            parfor si = 1:ns
                [gamma_pass(si, :), ssim_mean(si, :)] = process_segment( ...
                    seg_data{si}, comparisons, dose_pct, dist_mm, cutoff_fr, ...
                    normalize, hash8, overwrite);
            end
        else
            for si = 1:ns
                [gamma_pass(si, :), ssim_mean(si, :)] = process_segment( ...
                    seg_data{si}, comparisons, dose_pct, dist_mm, cutoff_fr, ...
                    normalize, hash8, overwrite);
            end
        end

        nproc = nproc + 1;
        beam_list(nproc)      = b;
        n_segments(nproc)     = ns;
        seg_gamma_pass{nproc} = gamma_pass;
        seg_ssim_mean{nproc}  = ssim_mean;

        elapsed = toc(run_timer);
        eta     = elapsed / bi * (nBeams - bi);
        fprintf(['  beam #%d done: %d segment(s) | mean gamma recon1=%.1f%% ' ...
                 'recon3=%.1f%% | mean SSIM recon1=%.1f%% recon3=%.1f%% | ' ...
                 'elapsed %.1fs, ETA %.1fs\n'], b, ns, ...
            mean(gamma_pass(:, 2), 'omitnan'), mean(gamma_pass(:, 3), 'omitnan'), ...
            mean(ssim_mean(:, 2),  'omitnan'), mean(ssim_mean(:, 3),  'omitnan'), ...
            elapsed, eta);
    end

    if nproc == 0
        error('step25_segment_metrics:NoBeamsProcessed', ...
            'No beams could be processed (all skipped).');
    end

    beam_list      = beam_list(1:nproc);
    n_segments     = n_segments(1:nproc);
    seg_gamma_pass = seg_gamma_pass(1:nproc);
    seg_ssim_mean  = seg_ssim_mean(1:nproc);

    %% ======================== SUMMARY ROLLUP ========================
    % Per-beam mean/std over segments, plus a pooled "all segments" aggregate.
    % This is the data the study's plots are drawn from -- no figures here.
    results = build_summary(patient_id, session, hash8, config, comparisons, ...
        beam_list, n_segments, seg_gamma_pass, seg_ssim_mean);

    if config.metrics_write_summary && ~isempty(sim_dir)
        summary_path = fullfile(sim_dir, ...
            sprintf('segment_metrics_summary_%s.mat', hash8));
        save(summary_path, '-struct', 'results', '-v7.3');
        fprintf('\n[STEP 2.5] Summary saved: %s\n', summary_path);
    end

    fprintf('[STEP 2.5] Complete: %d beam(s), %d segment(s) in %.1f s.\n', ...
        nproc, sum(n_segments), toc(run_timer));
end


%% =========================================================================
%  PER-SEGMENT WORK (runs on a parfor worker)
%% =========================================================================

function [gamma_pass, ssim_mean] = process_segment(S, comparisons, ...
        dose_pct, dist_mm, cutoff_fr, normalize, hash8, overwrite)
%PROCESS_SEGMENT Compute the four comparisons for one segment and fold the raw
%  masked result into the segment's CT_1 recon file. Returns the per-comparison
%  gamma pass rate (%) and mean local SSIM (%) for the summary rollup. A segment
%  already folded in (and not overwriting) is not recomputed: its stored scalars
%  are read back so the summary stays complete on a resumed run.
    nComp      = size(comparisons, 1);
    gamma_pass = nan(1, nComp);
    ssim_mean  = nan(1, nComp);

    ct1_file = S.ct1_recon_file;

    if ~overwrite && has_folded_metrics(ct1_file)
        [gamma_pass, ssim_mean] = read_folded_scalars(ct1_file, nComp);
        if all(~isnan(gamma_pass)) || all(~isnan(ssim_mean))
            return;
        end
    end

    % Cross-instance safety: several pipeline_simulate copies can reach Step 2.5
    % at once and would otherwise append to the SAME CT_1 recon .mat concurrently
    % (corrupting it, recon_dose included). Claim an atomic directory lock next to
    % the file; if a sibling holds it, leave the segment to them (NaN row, picked
    % up on a later resumed run). Within one run each segment is a distinct file,
    % so parfor workers never contend -- this guards only cross-instance overlap.
    lock_dir = [ct1_file, '.metrics_lock'];
    [ok, ~, msgid] = mkdir(lock_dir);
    if ~ok || strcmpi(msgid, 'MATLAB:MKDIR:DirectoryExists')
        return;   % a sibling instance is computing this segment
    end

    vol_size = size(S.rs_CT1);

    comp = struct('name', {}, 'ref_field', {}, 'tgt_field', {}, ...
        'mask_idx', {}, 'gamma_vals', {}, 'gamma_pass_rate', {}, ...
        'ssim_vals', {}, 'ssim_mean', {});

    for d = 1:nComp
        ref = get_dose_field(S, comparisons{d, 2});
        tgt = get_dose_field(S, comparisons{d, 3});

        % 10%-of-reference eval mask (fall back to the target if the reference
        % is all-zero), matching the study.
        if max(ref(:)) > 0
            mask = ref >= cutoff_fr * max(ref(:));
        else
            mask = tgt >= cutoff_fr * max(tgt(:));
        end
        mask_idx = uint32(find(mask));

        % Global gamma index (CPU, CalcGamma output suppressed).
        gmap = quiet_gamma_map(ref, tgt, dose_pct, dist_mm, S.spacing);
        if isempty(gmap)
            g_vals = single([]);
            g_pass = NaN;
        else
            g_vals = single(gmap(mask));
            g_pass = 100 * mean(double(g_vals) <= 1);
        end

        % Local SSIM map + mean over the same eval mask.
        [smap, s_frac] = compute_local_ssim(ref, tgt, mask);
        if isempty(smap)
            s_vals = single([]);
            s_pct  = NaN;
        else
            s_vals = single(smap(mask));
            s_pct  = 100 * s_frac;
        end

        gamma_pass(d) = g_pass;
        ssim_mean(d)  = s_pct;

        comp(d) = struct('name', comparisons{d, 1}, ...
            'ref_field', comparisons{d, 2}, 'tgt_field', comparisons{d, 3}, ...
            'mask_idx', mask_idx, 'gamma_vals', g_vals, ...
            'gamma_pass_rate', g_pass, 'ssim_vals', s_vals, 'ssim_mean', s_pct);
    end

    % Assemble the folded record. Only masked voxels are stored; expand a map on
    % demand with:  M = nan(segment_metrics.vol_size); M(c.mask_idx) = c.gamma_vals;
    segment_metrics = struct();
    segment_metrics.config_hash           = hash8;
    segment_metrics.beam_num              = S.beam;
    segment_metrics.seg_num               = S.seg;
    segment_metrics.spacing               = S.spacing;
    segment_metrics.vol_size              = vol_size;
    segment_metrics.normalize             = normalize;
    segment_metrics.recon_ct1_gain        = S.gain_CT1;
    segment_metrics.recon_ct3_gain        = S.gain_CT3;
    segment_metrics.gamma_dose_pct        = dose_pct;
    segment_metrics.gamma_dist_mm         = dist_mm;
    segment_metrics.gamma_dose_cutoff_pct = cutoff_fr * 100;
    segment_metrics.comparison            = comp;   %#ok<STRNU> saved below

    % Fold into the CT_1 recon file (append; only this file, per convention),
    % then release the lock whether or not the save succeeded.
    try
        save(ct1_file, 'segment_metrics', '-append');
    catch ME
        warning('step25_segment_metrics:SaveFailed', ...
            'Could not fold metrics into %s: %s', ct1_file, ME.message);
    end
    if isfolder(lock_dir), rmdir(lock_dir, 's'); end
end

function [gamma_pass, ssim_mean] = read_folded_scalars(ct1_file, nComp)
%READ_FOLDED_SCALARS Per-comparison gamma pass % and mean SSIM % already folded
%  into a CT_1 recon file, so a resumed run keeps the summary complete without
%  recomputing. Returns NaN rows if the stored form is absent/unreadable.
    gamma_pass = nan(1, nComp);
    ssim_mean  = nan(1, nComp);
    try
        L  = load(ct1_file, 'segment_metrics');
        gp = [L.segment_metrics.comparison.gamma_pass_rate];
        sm = [L.segment_metrics.comparison.ssim_mean];
        if numel(gp) == nComp, gamma_pass = gp; end
        if numel(sm) == nComp, ssim_mean  = sm; end
    catch
        % leave NaN rows -> the caller recomputes
    end
end


%% =========================================================================
%  SEGMENT GROUPING (adapted from study_pass_rates_allsegments)
%% =========================================================================

function seg_data = group_beam_segments(fields, spacing, ct1_str, ct3_str, normalize, beam)
%GROUP_BEAM_SEGMENTS Group a beam's loaded fields into per-segment volume sets.
%  Each element S carries recon_CT1/recon_CT3/rs_CT1/rs_CT3 (double), spacing,
%  seg, beam, the CT_1 recon file path (fold target) and the applied gains.
%  Segments missing the CT_1/CT_3 pair or with mismatched grids are skipped.
    seg_ids = arrayfun(@(k) seg_of_field(fields(k)), 1:numel(fields));
    uniq    = unique(seg_ids(~isnan(seg_ids)));
    uniq    = uniq(:).';

    seg_data = {};
    for s = uniq
        sub = fields(seg_ids == s);
        iA  = find_field_by_ct(sub, ct1_str);
        iB  = find_field_by_ct(sub, ct3_str);
        if isempty(iA) || isempty(iB)
            continue;   % need both CTs for the cross-CT comparisons
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
        S.ct1_recon_file = fA.recon_file;   % fold destination
        S.gain_CT1       = 1;
        S.gain_CT3       = 1;

        if ~isequal(size(S.recon_CT1), size(S.rs_CT1)) || ...
           ~isequal(size(S.recon_CT3), size(S.rs_CT3)) || ...
           ~isequal(size(S.rs_CT1),    size(S.rs_CT3))
            warning('step25_segment_metrics:GridMismatch', ...
                'Beam #%d seg %d: grid mismatch; skipping segment.', beam, s);
            continue;
        end
        if isempty(S.ct1_recon_file)
            warning('step25_segment_metrics:NoCT1File', ...
                'Beam #%d seg %d: CT_1 recon file path missing; skipping.', beam, s);
            continue;
        end

        if normalize
            S.gain_CT1  = least_squares_gain(S.rs_CT1, S.recon_CT1);
            S.gain_CT3  = least_squares_gain(S.rs_CT3, S.recon_CT3);
            S.recon_CT1 = S.recon_CT1 * S.gain_CT1;
            S.recon_CT3 = S.recon_CT3 * S.gain_CT3;
        end

        seg_data{end+1} = S; %#ok<AGROW>
    end
end

function beam_set = load_beam_set(patient_id, session, config, beam)
%LOAD_BEAM_SET Load every segment/CT of a beam via load_recon_dose_data (set).
%  RS truth is needed; ETHOS/CBCT are not (skipped for speed).
    args = {'Mode', 'set', 'Beam', beam, 'IncludeEthos', false, ...
        'IncludeCBCT', false, 'IncludeRS', true};
    if ~strcmpi(config.metrics_plan_type, 'any')
        args = [args, {'PlanType', config.metrics_plan_type}];
    end
    if isfield(config, 'config_hash') && ~isempty(config.config_hash)
        args = [args, {'Hash', config.config_hash}];
    end
    beam_set = load_recon_dose_data(patient_id, session, config, args{:});
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
%  g minimizes ||rs - g*recon||^2 over the truth's 10% region. Falls back to 1.
    rs_truth = double(rs_truth);
    recon    = double(recon);
    if max(rs_truth(:)) > 0
        mask = rs_truth >= 0.10 * max(rs_truth(:));
    else
        mask = true(size(rs_truth));
    end
    r     = recon(mask);
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
            error('step25_segment_metrics:BadField', ...
                'Unknown comparison volume "%s".', fieldname);
    end
end


%% =========================================================================
%  GAMMA / SSIM PRIMITIVES
%% =========================================================================

function gmap = quiet_gamma_map(ref, tgt, dose_pct, dist_mm, spacing)
%QUIET_GAMMA_MAP Full global gamma-index volume (CalcGamma output suppressed).
%  Global gamma ('local',0), DTA search limited to 2x the criterion, restricted
%  search for speed, forced onto the CPU so many parfor workers do not contend
%  over the GPU. Returns [] on failure. Matches study_pass_rates_allsegments.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(ref));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(tgt));
    limit      = 2 * dist_mm;
    gmap = [];
    try
        evalc(['gmap = CalcGamma(ref_struct, tgt_struct, dose_pct, dist_mm, ', ...
               '''local'', 0, ''limit'', limit, ''restrict'', 1, ''cpu'', 1);']);
    catch
        gmap = [];
    end
end


%% =========================================================================
%  SUMMARY ROLLUP
%% =========================================================================

function S = build_summary(patient_id, session, hash8, config, comparisons, ...
        beam_list, n_segments, seg_gamma_pass, seg_ssim_mean)
%BUILD_SUMMARY Per-beam mean/std over segments + pooled aggregate, both metrics.
    nComp = size(comparisons, 1);
    nB    = numel(beam_list);

    gamma_mean = nan(nB, nComp); gamma_std = nan(nB, nComp);
    ssim_mean  = nan(nB, nComp); ssim_std  = nan(nB, nComp);
    for n = 1:nB
        gamma_mean(n, :) = mean(seg_gamma_pass{n}, 1, 'omitnan');
        gamma_std(n, :)  = std(seg_gamma_pass{n}, 0, 1, 'omitnan');
        ssim_mean(n, :)  = mean(seg_ssim_mean{n}, 1, 'omitnan');
        ssim_std(n, :)   = std(seg_ssim_mean{n}, 0, 1, 'omitnan');
    end

    all_gamma = vertcat(seg_gamma_pass{:});   % [sum(nSeg) x nComp]
    all_ssim  = vertcat(seg_ssim_mean{:});

    S = struct();
    S.patient_id    = patient_id;
    S.session       = session;
    S.config_hash   = hash8;
    S.plan_type     = config.metrics_plan_type;
    S.normalize     = logical(config.metrics_normalize);
    S.comparisons   = comparisons(:, 1)';
    S.beams         = beam_list;
    S.n_segments    = n_segments;
    S.all_n_segments = size(all_gamma, 1);

    S.gamma_dose_pct        = config.gamma_dose_pct;
    S.gamma_dist_mm         = config.gamma_dist_mm;
    S.gamma_dose_cutoff_pct = config.gamma_dose_cutoff_pct;

    S.gamma.mean_pass = gamma_mean;   % [nBeam x nComp] gamma pass rate %
    S.gamma.std_pass  = gamma_std;
    S.gamma.all_mean  = mean(all_gamma, 1, 'omitnan');
    S.gamma.all_std   = std(all_gamma, 0, 1, 'omitnan');
    S.gamma.seg_pass  = seg_gamma_pass;   % {1 x nBeam}, each [nSeg x nComp]

    S.ssim.mean       = ssim_mean;    % [nBeam x nComp] mean local SSIM %
    S.ssim.std        = ssim_std;
    S.ssim.all_mean   = mean(all_ssim, 1, 'omitnan');
    S.ssim.all_std    = std(all_ssim, 0, 1, 'omitnan');
    S.ssim.seg_mean   = seg_ssim_mean;    % {1 x nBeam}, each [nSeg x nComp]
end


%% =========================================================================
%  SMALL HELPERS
%% =========================================================================

function config = default_field(config, name, value)
%DEFAULT_FIELD Set config.(name)=value only if absent/empty.
    if ~isfield(config, name) || isempty(config.(name))
        config.(name) = value;
    end
end

function tf = has_folded_metrics(mat_file)
%HAS_FOLDED_METRICS True if the .mat already carries a 'segment_metrics' var.
    tf = false;
    if isempty(mat_file) || ~isfile(mat_file), return; end
    try
        tf = ismember('segment_metrics', who('-file', mat_file));
    catch
        tf = false;
    end
end

function ensure_cpu_pool(desired_workers)
%ENSURE_CPU_POOL Reuse a running pool; otherwise start a local one. Non-fatal if
%  the Parallel Computing Toolbox is absent (the caller's loop runs serially).
    if exist('parpool', 'file') ~= 2
        fprintf('  [WARN] Parallel Computing Toolbox not found; running serially.\n');
        return;
    end
    if ~isempty(gcp('nocreate'))
        return;   % reuse the sim's pool as-is
    end
    try
        parpool('local', desired_workers);
        fprintf('  [Pool] Started local pool with %d worker(s).\n', desired_workers);
    catch ME
        fprintf('  [WARN] Could not start parallel pool (%s). Running serially.\n', ...
            ME.message);
    end
end
