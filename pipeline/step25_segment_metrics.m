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
%   over a beam's segments; beam mode: over the beams) while the sim's GPU work
%   is already done.
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
%   BEAM MODE (config.reconstruction_mode = 'beam'): Step 2 simulated each beam's
%   summed dose instead of each segment, so each beam is treated as ONE
%   "segment" (seg_num 0 = the whole beam): truths are the beam totals
%   total_dose_<CT>_B<n>.mat, recons are total_dose_<CT>_B<n>_recon_<hash>.mat,
%   and the metrics fold into the CT_1 beam recon. All beams are loaded first,
%   then ONE parfor runs over the beams (a parfor over a beam's single segment
%   would run the beams one at a time). metrics_plan_type is not
%   applied (the beam totals sum every segment on that CT). When
%   config.config_hash is set, the mode is read from that hash's registered
%   config (hashtoconfig), not from config.reconstruction_mode, so a config
%   left at the 'segment' default cannot load an older segment run instead.
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
%           .config_hash           ['']   8-char sim hash of the recons to use
%                                  ('' = from the config). Also sets the mode.
%           .reconstruction_mode   ['segment'] 'segment' | 'beam'. Overridden
%                                  by the hash's registered mode when
%                                  .config_hash is set. Beam mode also needs
%                                  .gruneisen_method and either .config_hash
%                                  or the full sim CONFIG.
%           .metrics_overwrite     [false] Recompute even if already folded in.
%           .metrics_write_summary [true] Write the per-beam rollup .mat.
%           .metrics_noise_floor   [true] Compute the noise-only null pass rate
%                                  (once per session; cached by the sim hash).
%           .metrics_noise_minutes [30]   Ensemble time budget (min) handed to
%                                  noise_ensemble_error_bars (TimeBudgetMin).
%           .use_parallel          [true] parfor over a beam's segments, or
%                                  over the beams in beam mode (CPU).
%           .metrics_auto_workers  [true] size a NEW pool to the physical CPU
%                                  core count (best for this CPU-bound work);
%                                  set false to use num_parallel_workers exactly.
%           .num_parallel_workers  [8]    pool size when auto-sizing is off or
%                                  core detection fails (fallback). Ignored when
%                                  a pool is already running (that pool is reused).
%
%   OUTPUTS:
%       results - Summary struct (also saved to segment_metrics_summary_<hash>.mat
%                 when metrics_write_summary): per-beam and pooled mean/std of the
%                 gamma pass rate (%) and mean local SSIM (%), per comparison,
%                 plus the per-segment matrices and `.noise_floor` (the noise-only
%                 null ensemble: mean/std pass rate for the session, or [] when
%                 disabled/unavailable). This is the data the study's plots are
%                 built from -- no figures are produced here. A plotting script
%                 draws `.noise_floor.mean_pass_rate` as the horizontal null line
%                 and the differential (truth1_vs_recon1 pass) - (noise mean).
%                 `.full_vs_sum` (per CT): gamma %, SSIM %, sum ratio of the
%                 summed field recons vs the Step 2.05 full-field recon.
%                 `.timing` holds the gamma and SSIM seconds of every segment
%                 computed this run (worker time), plus their averages/totals.
%
%   DEPENDENCIES:
%       - CalcGamma.m, compute_local_ssim.m (utils), load_recon_dose_data.m,
%         load_beam_total_set.m (beam mode), hashtoconfig.m (mode of a hash)
%       - noise_ensemble_error_bars.m (utils) for the noise-only null floor.
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
    config = default_field(config, 'config_hash',           '');
    config = default_field(config, 'reconstruction_mode',   'segment');
    config = default_field(config, 'metrics_overwrite',     false);
    config = default_field(config, 'metrics_write_summary', true);
    config = default_field(config, 'metrics_noise_floor',   true);
    config = default_field(config, 'metrics_noise_minutes', 30);
    config = default_field(config, 'use_parallel',          true);
    config = default_field(config, 'metrics_auto_workers',  true);
    config = default_field(config, 'num_parallel_workers',  8);

    if exist('CalcGamma', 'file') ~= 2
        error('step25_segment_metrics:NoCalcGamma', ...
            'CalcGamma not found on the path; cannot compute gamma pass rates.');
    end

    % A given hash decides the mode: its registered config says whether Step 2
    % ran per segment or per beam (the hash leaves reconstruction_mode out when
    % it is 'segment'). Trusting config.reconstruction_mode instead lets a config
    % left at the 'segment' default load an older segment run's recons.
    if ~isempty(config.config_hash)
        run_cfg  = hashtoconfig(config.config_hash, config);
        run_mode = 'segment';
        if isfield(run_cfg, 'reconstruction_mode')
            run_mode = run_cfg.reconstruction_mode;
        end
        if ~strcmp(run_mode, config.reconstruction_mode)
            warning('step25_segment_metrics:ModeFromHash', ...
                'config.reconstruction_mode = ''%s'' but hash %s is a ''%s''-mode run; using ''%s''.', ...
                config.reconstruction_mode, config.config_hash, run_mode, run_mode);
        end
        config.reconstruction_mode = run_mode;
    end
    beam_mode = strcmp(config.reconstruction_mode, 'beam');

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

    fprintf('\n=========================================================\n');
    fprintf('  [STEP 2.5] Per-Segment Gamma + SSIM Metrics\n');
    fprintf('  Patient %s | %s | plan=%s | mode=%s\n', patient_id, session, ...
        config.metrics_plan_type, config.reconstruction_mode);
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
    have_pool = false;
    n_workers = 1;   % used by the time-remaining estimate
    if config.use_parallel
        if config.metrics_auto_workers
            nworkers = detect_num_cores(config.num_parallel_workers);
        else
            nworkers = config.num_parallel_workers;
        end
        ensure_cpu_pool(nworkers);
        if exist('gcp', 'file') == 2 && ~isempty(gcp('nocreate'))
            have_pool = true;
            pool      = gcp('nocreate');
            n_workers = pool.NumWorkers;
        end
    end

    %% ======================== PER-BEAM PROCESSING ========================
    run_timer = tic;
    % TIMING: gamma / SSIM seconds of every segment computed this run (same
    % running-average ETA as Step 2's per-field sims in pipeline_simulate).
    all_gamma_sec   = zeros(1, 0);
    all_ssim_sec    = zeros(1, 0);
    beam_list       = nan(1, nBeams);
    n_segments      = zeros(1, nBeams);
    seg_gamma_pass  = cell(1, nBeams);   % {beam} = [nSeg x nComp] gamma pass %
    seg_ssim_mean   = cell(1, nBeams);   % {beam} = [nSeg x nComp] mean SSIM %
    nproc           = 0;
    sim_dir         = '';
    hash8           = '';
    recon_sum       = struct();          % .CT_1 / .CT_3 summed field recons
    spacing         = [];
    queued_segs     = {};                % beam mode: each beam's one segment

    for bi = 1:nBeams
        b = beams(bi);
        try
            if beam_mode
                fprintf('\n[Beam %d/%d] beam #%d: loading beam totals + beam recons...\n', bi, nBeams, b);
                beam_set = load_beam_total_set(patient_id, session, config, b, ct1_str, ct3_str);
            else
                fprintf('\n[Beam %d/%d] beam #%d: loading segments...\n', bi, nBeams, b);
                beam_set = load_beam_set(patient_id, session, config, b);
            end
        catch ME
            warning('step25_segment_metrics:SkipBeam', ...
                'Skipping beam #%d (load failed): %s', b, ME.message);
            continue;
        end

        % Capture the recon location / hash from the first successful load.
        if isempty(hash8)
            hash8   = beam_set.config_hash;
            sim_dir = get_simulation_dir(config.working_dir, patient_id, session, hash8);
        end

        % Running per-CT sum of the raw (un-normalized) field recons, compared
        % below with the Step 2.05 full-field (total RS dose) recon.
        for k = 1:numel(beam_set.fields)
            lbl = field_ct_label(beam_set.fields(k));
            if ~strcmpi(lbl, ct1_str) && ~strcmpi(lbl, ct3_str)
                continue;
            end
            if ~isfield(recon_sum, lbl)
                recon_sum.(lbl) = zeros(size(beam_set.fields(k).recon_dose));
            end
            recon_sum.(lbl) = recon_sum.(lbl) + double(beam_set.fields(k).recon_dose);
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

        % Beam mode: queue the beam's one segment (seg 0); every beam runs in
        % one parfor after this loop.
        if beam_mode
            queued_segs{end+1} = seg_data{1}; %#ok<AGROW>
            continue;
        end

        % TIMING: seed the live estimate with this beam's segment count, a guess
        % for the beams not yet loaded (mean segments/beam so far), the worker
        % count, and the gamma/SSIM times recorded so far.
        later_est = round(mean([n_segments(1:nproc), ns]) * (nBeams - bi));
        report_segment_progress(struct('kind', 'init', 'n_segments', ns, ...
            'later_est', later_est, 'workers', n_workers, ...
            'gamma_times', all_gamma_sec, 'ssim_times', all_ssim_sec));

        [gamma_pass, ssim_mean, gamma_sec, ssim_sec] = compute_segments(seg_data, ...
            have_pool, comparisons, dose_pct, dist_mm, cutoff_fr, normalize, hash8, ...
            config.metrics_overwrite);

        % TIMING: keep the times of the segments actually computed.
        computed      = ~isnan(gamma_sec);
        all_gamma_sec = [all_gamma_sec, gamma_sec(computed)']; %#ok<AGROW>
        all_ssim_sec  = [all_ssim_sec,  ssim_sec(computed)'];  %#ok<AGROW>

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

    %% ============ BEAM MODE: ONE PARFOR OVER ALL BEAMS ============
    % Each beam is a single "segment", so a parfor over one beam's segments
    % would run the beams one at a time. The loop above only loaded and queued
    % them; here all beams run in parallel (one worker per beam).
    if beam_mode && ~isempty(queued_segs)
        nq = numel(queued_segs);
        fprintf('\n[STEP 2.5] Beam mode: gamma + SSIM for %d beam(s) in parallel...\n', nq);
        report_segment_progress(struct('kind', 'init', 'n_segments', nq, ...
            'later_est', 0, 'workers', n_workers, ...
            'gamma_times', all_gamma_sec, 'ssim_times', all_ssim_sec));

        [gamma_pass, ssim_mean, gamma_sec, ssim_sec] = compute_segments(queued_segs, ...
            have_pool, comparisons, dose_pct, dist_mm, cutoff_fr, normalize, hash8, ...
            config.metrics_overwrite);

        computed      = ~isnan(gamma_sec);
        all_gamma_sec = gamma_sec(computed)';
        all_ssim_sec  = ssim_sec(computed)';

        for q = 1:nq
            nproc = nproc + 1;
            beam_list(nproc)      = queued_segs{q}.beam;
            n_segments(nproc)     = 1;
            seg_gamma_pass{nproc} = gamma_pass(q, :);
            seg_ssim_mean{nproc}  = ssim_mean(q, :);
            fprintf(['  beam #%d done: gamma recon1=%.1f%% recon3=%.1f%% | ' ...
                     'SSIM recon1=%.1f%% recon3=%.1f%%\n'], beam_list(nproc), ...
                gamma_pass(q, 2), gamma_pass(q, 3), ssim_mean(q, 2), ssim_mean(q, 3));
        end
        queued_segs = {};   % free the beam volumes before the noise floor
    end

    if nproc == 0
        error('step25_segment_metrics:NoBeamsProcessed', ...
            'No beams could be processed (all skipped).');
    end

    beam_list      = beam_list(1:nproc);
    n_segments     = n_segments(1:nproc);
    seg_gamma_pass = seg_gamma_pass(1:nproc);
    seg_ssim_mean  = seg_ssim_mean(1:nproc);

    %% ======================== NOISE-ONLY NULL FLOOR ========================
    % Null hypothesis: the gamma pass rate of the CT_1 truth against a NOISE-ONLY
    % reconstruction (true acoustic signal nulled). noise_ensemble_error_bars runs
    % an ensemble of such reconstructions and returns the mean +/- std pass rate.
    % Computed ONCE for the whole session (the ensemble is cached by the sim config
    % hash -- the same hash the recons use -- so every beam/segment shares it), and
    % stored on the summary as `noise_floor` for later comparison with the per-beam
    % pass rates. Only meaningful for gamma, so it uses gamma_dose_pct/gamma_dist_mm.
    noise_floor = [];
    if config.metrics_noise_floor
        fprintf('\n[STEP 2.5] Computing noise-only null floor (once per session)...\n');
        noise_floor = compute_session_noise_floor(patient_id, session, config, ...
            dose_pct, dist_mm, hash8);
        if ~isempty(noise_floor)
            fprintf('[STEP 2.5] Noise floor: %.2f +/- %.2f %% over %d sample(s).\n', ...
                noise_floor.mean_pass_rate, noise_floor.std_pass_rate, ...
                noise_floor.num_samples);
        end
    end

    %% ============ FULL-FIELD RECON vs SUM OF FIELD RECONS ============
    % Step 2.05 simulates each CT's total RS dose in one k-Wave run. Linearity
    % says that should equal the sum of the per-field recons. Reference = the
    % full-field recon; target = the summed recon (both raw, no gain).
    fprintf('\n[STEP 2.5] Comparing full-field recon with summed field recons...\n');
    full_vs_sum = struct('ct_label', {}, 'gamma_pass_rate', {}, 'ssim_mean', {}, ...
        'sum_ratio', {}, 'max_full_Gy', {}, 'max_sum_Gy', {});
    ct_labels = {ct1_str, ct3_str};
    for k = 1:2
        lbl       = ct_labels{k};
        full_file = fullfile(sim_dir, sprintf('full_field_dose_%s_%s.mat', lbl, hash8));
        if ~isfile(full_file) || ~isfield(recon_sum, lbl)
            fprintf('  %s: skipped (no full-field recon or no field recons).\n', lbl);
            continue;
        end
        L    = load(full_file, 'recon_dose');
        ref  = double(L.recon_dose);
        tgt  = recon_sum.(lbl);
        if ~isequal(size(ref), size(tgt))
            warning('step25_segment_metrics:FullFieldGrid', ...
                '%s: full-field and summed recon grids differ; skipping.', lbl);
            continue;
        end
        mask = ref >= cutoff_fr * max(ref(:));

        gmap = quiet_gamma_map(ref, tgt, dose_pct, dist_mm, spacing);
        g_pass = NaN;
        if ~isempty(gmap)
            g_pass = 100 * mean(gmap(mask) <= 1);
        end
        [~, s_frac] = compute_local_ssim(ref, tgt, mask);

        full_vs_sum(end+1) = struct('ct_label', lbl, 'gamma_pass_rate', g_pass, ...
            'ssim_mean', 100 * s_frac, 'sum_ratio', sum(tgt(:)) / sum(ref(:)), ...
            'max_full_Gy', max(ref(:)), 'max_sum_Gy', max(tgt(:))); %#ok<AGROW>
        fprintf('  %s: gamma %.1f%% | SSIM %.1f%% | sum(fields)/full = %.4f\n', ...
            lbl, g_pass, 100 * s_frac, full_vs_sum(end).sum_ratio);
    end

    %% ======================== SUMMARY ROLLUP ========================
    % Per-beam mean/std over segments, plus a pooled "all segments" aggregate.
    % This is the data the study's plots are drawn from -- no figures here.
    results = build_summary(patient_id, session, hash8, config, comparisons, ...
        beam_list, n_segments, seg_gamma_pass, seg_ssim_mean);
    results.noise_floor = noise_floor;   % [] when disabled/unavailable
    results.full_vs_sum = full_vs_sum;   % full-field recon vs summed field recons

    % TIMING: per-segment gamma / SSIM worker times (averages NaN if none computed).
    results.timing = struct();
    results.timing.gamma_sec_per_segment = all_gamma_sec;
    results.timing.ssim_sec_per_segment  = all_ssim_sec;
    results.timing.avg_gamma_sec         = mean(all_gamma_sec);
    results.timing.avg_ssim_sec          = mean(all_ssim_sec);
    results.timing.total_gamma_sec       = sum(all_gamma_sec);
    results.timing.total_ssim_sec        = sum(all_ssim_sec);

    if config.metrics_write_summary && ~isempty(sim_dir)
        summary_path = fullfile(sim_dir, ...
            sprintf('segment_metrics_summary_%s.mat', hash8));
        save(summary_path, '-struct', 'results', '-v7.3');
        fprintf('\n[STEP 2.5] Summary saved: %s\n', summary_path);
    end

    fprintf('[STEP 2.5] Complete: %d beam(s), %d segment(s) in %.1f s.\n', ...
        nproc, sum(n_segments), toc(run_timer));
    fprintf(['[STEP 2.5] Timing: %d segment(s) computed | gamma %.1f s total ' ...
             '(%.1f s/segment) | SSIM %.1f s total (%.1f s/segment)\n'], ...
        numel(all_gamma_sec), results.timing.total_gamma_sec, ...
        results.timing.avg_gamma_sec, results.timing.total_ssim_sec, ...
        results.timing.avg_ssim_sec);
end


function report_segment_progress(msg)
%REPORT_SEGMENT_PROGRESS Live Step 2.5 time-remaining estimate (runs on the client).
%  Same scheme as Step 2's report_field_progress in pipeline_simulate. An 'init'
%  message (sent before each beam's segment loop) seeds the countdown, the
%  estimated segments in later beams, the worker count and the gamma/SSIM times
%  recorded so far. Each segment then sends a 'done' message with its gamma and
%  SSIM seconds; NaN means it was skipped (already folded in or held by a
%  sibling instance), which only advances the countdown.
%  ETA = segments left * (avg gamma + avg SSIM s/segment) / workers.
    persistent remaining later_est workers gamma_times ssim_times
    if strcmp(msg.kind, 'init')
        remaining   = msg.n_segments;
        later_est   = msg.later_est;
        workers     = max(msg.workers, 1);
        gamma_times = msg.gamma_times;
        ssim_times  = msg.ssim_times;
        return;
    end
    remaining = max(remaining - 1, 0);
    if isnan(msg.gamma_sec)
        return;
    end
    gamma_times(end+1) = msg.gamma_sec;
    ssim_times(end+1)  = msg.ssim_sec;
    avg_gamma = mean(gamma_times);
    avg_ssim  = mean(ssim_times);
    eta_min   = (remaining + later_est) * (avg_gamma + avg_ssim) / workers / 60;
    fprintf(['    [TIMING] beam #%d seg %d: gamma %.1f s, SSIM %.1f s  --  ' ...
             '%d left in beam (+~%d in later beams), est. %.1f min remaining ' ...
             '(@ %.1f s gamma + %.1f s SSIM per segment, %d worker(s))\n'], ...
        msg.beam, msg.seg, msg.gamma_sec, msg.ssim_sec, remaining, later_est, ...
        eta_min, avg_gamma, avg_ssim, workers);
end


%% =========================================================================
%  PER-SEGMENT WORK (runs on a parfor worker)
%% =========================================================================

function [gamma_pass, ssim_mean, gamma_sec, ssim_sec] = compute_segments(seg_data, ...
        use_pool, comparisons, dose_pct, dist_mm, cutoff_fr, normalize, hash8, overwrite)
%COMPUTE_SEGMENTS Run process_segment on every entry of seg_data: parfor when a
%  pool is up, else a plain loop. Segment mode passes one beam's segments; beam
%  mode passes every beam (one entry each). Output rows follow seg_data order.
    ns         = numel(seg_data);
    nComp      = size(comparisons, 1);
    gamma_pass = nan(ns, nComp);
    ssim_mean  = nan(ns, nComp);
    gamma_sec  = nan(ns, 1);   % NaN = segment not computed here (skipped)
    ssim_sec   = nan(ns, 1);

    if use_pool
        % Workers send one 'done' message per segment; the client prints it.
        progress_queue = parallel.pool.DataQueue;
        afterEach(progress_queue, @report_segment_progress);
        parfor si = 1:ns
            S = seg_data{si};
            [gamma_pass(si, :), ssim_mean(si, :), gamma_sec(si), ssim_sec(si)] = ...
                process_segment(S, comparisons, dose_pct, dist_mm, cutoff_fr, ...
                normalize, hash8, overwrite);
            send(progress_queue, struct('kind', 'done', 'beam', S.beam, 'seg', S.seg, ...
                'gamma_sec', gamma_sec(si), 'ssim_sec', ssim_sec(si)));
        end
    else
        for si = 1:ns
            S = seg_data{si};
            [gamma_pass(si, :), ssim_mean(si, :), gamma_sec(si), ssim_sec(si)] = ...
                process_segment(S, comparisons, dose_pct, dist_mm, cutoff_fr, ...
                normalize, hash8, overwrite);
            report_segment_progress(struct('kind', 'done', 'beam', S.beam, 'seg', S.seg, ...
                'gamma_sec', gamma_sec(si), 'ssim_sec', ssim_sec(si)));
        end
    end
end

function [gamma_pass, ssim_mean, gamma_sec, ssim_sec] = process_segment(S, comparisons, ...
        dose_pct, dist_mm, cutoff_fr, normalize, hash8, overwrite)
%PROCESS_SEGMENT Compute the four comparisons for one segment and fold the raw
%  masked result into the segment's CT_1 recon file. Returns the per-comparison
%  gamma pass rate (%) and mean local SSIM (%) for the summary rollup, plus the
%  seconds spent on gamma and on SSIM (summed over the comparisons; NaN when the
%  segment is not computed here). A segment already folded in (and not
%  overwriting) is not recomputed: its stored scalars are read back so the
%  summary stays complete on a resumed run.
    nComp      = size(comparisons, 1);
    gamma_pass = nan(1, nComp);
    ssim_mean  = nan(1, nComp);
    gamma_sec  = NaN;
    ssim_sec   = NaN;

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

    vol_size  = size(S.rs_CT1);
    gamma_sec = 0;
    ssim_sec  = 0;

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
        t_metric  = tic;
        gmap      = quiet_gamma_map(ref, tgt, dose_pct, dist_mm, S.spacing);
        gamma_sec = gamma_sec + toc(t_metric);
        if isempty(gmap)
            g_vals = single([]);
            g_pass = NaN;
        else
            g_vals = single(gmap(mask));
            g_pass = 100 * mean(double(g_vals) <= 1);
        end

        % Local SSIM map + mean over the same eval mask.
        t_metric       = tic;
        [smap, s_frac] = compute_local_ssim(ref, tgt, mask);
        ssim_sec       = ssim_sec + toc(t_metric);
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
%  Global gamma ('local',0), CalcGamma's default search limit (2x the DTA,
%  gamma capped at 2), restricted search for speed, forced onto the CPU so many
%  parfor workers do not contend over the GPU. Returns [] on failure.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(ref));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(tgt));
    gmap = [];
    try
        evalc(['gmap = CalcGamma(ref_struct, tgt_struct, dose_pct, dist_mm, ', ...
               '''local'', 0, ''restrict'', 1, ''cpu'', 1);']);
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
%  NOISE-ONLY NULL FLOOR
%% =========================================================================

function ens = compute_session_noise_floor(patient_id, session, config, ...
        dose_pct, dist_mm, hash8)
%COMPUTE_SESSION_NOISE_FLOOR Noise-only null ensemble for the whole session.
%  Returns the noise_ensemble_error_bars struct (mean/std pass rate etc.) or []
%  on failure. The ensemble is cached by the sim config hash, so it is computed
%  once and reused across every beam/segment. Because several pipeline_simulate
%  instances can reach Step 2.5 together, the FIRST (uncached) ~30 min compute is
%  guarded by a lock so only one instance pays it; siblings skip and pick up the
%  cached value on a later run. A present cache is loaded without locking.
    ens = [];

    % noise_ensemble_error_bars caches under config.patient_id/session.
    config.patient_id = patient_id;
    config.session    = session;

    % Same path noise_ensemble_error_bars writes to (its default CacheDir).
    nf_hash    = compute_sim_config_hash(config);   % equals hash8 for the sim CONFIG
    cache_file = fullfile(get_simulation_dir(config.working_dir, patient_id, session, nf_hash), ...
        sprintf('noise_ensemble_%s_%s_%s.mat', patient_id, session, nf_hash));

    have_lock = false;
    lock_dir  = [cache_file, '.lock'];
    if ~isfile(cache_file)
        [ok, ~, msgid] = mkdir(lock_dir);
        have_lock = ok && ~strcmpi(msgid, 'MATLAB:MKDIR:DirectoryExists');
        if ~have_lock && ~isfile(cache_file)
            fprintf('[STEP 2.5] Noise floor: a sibling instance is computing it; skipping this run.\n');
            return;
        end
    end

    try
        ens = run_noise_ensemble(patient_id, session, config, dose_pct, dist_mm, hash8);
    catch ME
        warning('step25_segment_metrics:NoiseFloorFailed', ...
            'Noise floor failed (%s); the summary will omit it.', ME.message);
        ens = [];
    end

    if have_lock && isfolder(lock_dir)
        rmdir(lock_dir, 's');
    end
end

function ens = run_noise_ensemble(patient_id, session, config, dose_pct, dist_mm, hash8)
%RUN_NOISE_ENSEMBLE Load the reference-CT geometry + summed RS truth and hand them
%  to noise_ensemble_error_bars. Mirrors study_pass_rates_allsegments' noise-floor
%  setup, but uses the real sim CONFIG (so the ensemble hash already matches the
%  recon hash -- no build_noise_sim_config indirection needed).
    tot = load_recon_dose_data(patient_id, session, config, ...
        'Mode', 'total', 'IncludeCBCT', true, 'IncludeEthos', false, ...
        'IncludeRS', true, 'Hash', hash8);

    ct1_field = sprintf('CT_%d', min(config.metrics_ct_pair));
    if ~isfield(tot, 'cbct') || ~isfield(tot.cbct, ct1_field)
        error('step25_segment_metrics:NoiseFloorNoCBCT', ...
            'Reference CBCT geometry %s not available for the noise floor.', ct1_field);
    end

    sct        = tot.cbct.(ct1_field);
    ref_truth  = double(tot.rs_dose);
    spacing_mm = tot.metadata.spacing(:)';
    if ~isfield(sct, 'spacing') || isempty(sct.spacing)
        sct.spacing = spacing_mm;
    end
    if ~isfield(sct, 'origin') || isempty(sct.origin)
        if isfield(tot.metadata, 'origin') && ~isempty(tot.metadata.origin)
            sct.origin = tot.metadata.origin;
        else
            sct.origin = [0, 0, 0];
        end
    end

    beam_meta = [];
    if isfield(tot.metadata, 'beam_metadata')
        beam_meta = tot.metadata.beam_metadata;
    end

    % Plan-level floor: placement uses the summed dose + all beams, so a single
    % representative gantry is fine (the reference beam's, else 0).
    gantry_angle = 0;
    if ~isempty(beam_meta) && isfield(beam_meta(1), 'gantry_angle') ...
            && ~isempty(beam_meta(1).gantry_angle)
        gantry_angle = beam_meta(1).gantry_angle;
    end

    crit_label = sprintf('%g%%/%g mm', dose_pct, dist_mm);
    ens = noise_ensemble_error_bars(config, ref_truth, spacing_mm, sct, ...
        gantry_angle, beam_meta, ...
        'TimeBudgetMin', config.metrics_noise_minutes, ...
        'GammaCriteria', {dose_pct, dist_mm, crit_label}, ...
        'Normalize', logical(config.metrics_normalize));
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
%ENSURE_CPU_POOL Reuse a running pool; otherwise start a local one sized to
%  desired_workers (clamped to the local cluster's NumWorkers cap). Non-fatal if
%  the Parallel Computing Toolbox is absent (the caller's loop runs serially).
    if exist('parpool', 'file') ~= 2
        fprintf('  [WARN] Parallel Computing Toolbox not found; running serially.\n');
        return;
    end
    if ~isempty(gcp('nocreate'))
        return;   % reuse a running pool (e.g. the sim's) as-is
    end
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
    try
        parpool('local', nw);
        fprintf('  [Pool] Started local pool with %d worker(s).\n', nw);
    catch ME
        fprintf('  [WARN] Could not start parallel pool (%s). Running serially.\n', ...
            ME.message);
    end
end

function n = detect_num_cores(default_workers)
%DETECT_NUM_CORES Physical CPU core count for sizing the worker pool.
%  Prefers the local cluster profile's NumWorkers (what parpool uses by default,
%  = physical cores under the default profile), then feature('numcores'); falls
%  back to default_workers if neither is available or returns something odd.
%  Physical (not logical/hyperthreaded) cores is the right target: the per-segment
%  gamma/SSIM is CPU- and memory-bandwidth-bound, so hyperthreads add contention
%  rather than throughput, and each MATLAB worker already runs single-threaded.
    if nargin < 1 || isempty(default_workers), default_workers = 8; end
    n = [];
    try
        c = parcluster('local');
        n = c.NumWorkers;
    catch
        n = [];
    end
    if isempty(n) || ~isfinite(n) || n < 1
        try
            n = feature('numcores');   % undocumented but widely available
        catch
            n = [];
        end
    end
    if isempty(n) || ~isfinite(n) || n < 1
        n = default_workers;
    end
    n = double(round(n));
end
