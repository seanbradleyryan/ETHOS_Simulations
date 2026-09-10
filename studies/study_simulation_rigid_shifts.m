%% =========================================================================
%  STUDY_SIMULATION_RIGID_SHIFTS.m
%  Rigid lateral-shift sensitivity of the IRAI gamma / SSIM metrics.
%
%  PURPOSE
%    Characterize how quickly reconstruction agreement collapses when the dose
%    is rigidly displaced sideways in the patient, with the anatomy AND the
%    ultrasound array held fixed. For one random segment per beam (CT_1 only)
%    this script:
%      1. Reconstructs the un-shifted field with a simple 10-iteration time
%         reversal recon and scores it (gamma + local SSIM) against its own
%         RayStation truth.
%      2. Re-runs the reconstruction with the SOURCE dose translated laterally
%         left, then right, by 2^n mm (n = 0,1,2,...), stopping each direction
%         once the dose centroid leaves the computational grid. Every shifted
%         recon is scored against the SAME original (un-shifted) truth.
%      3. Plots gamma pass rate and mean local SSIM as functions of the signed
%         lateral shift, one curve per beam, in a single figure (FIGURE 1).
%      4. Plots, per beam, a coronal max-projection overlay of the true dose
%         (green) over one representative translated dose (red), all beams in a
%         single, SEPARATE figure (FIGURE 2).
%
%    The per-shift reconstructions are the expensive part; they are run in a
%    parfor loop with one worker pinned per GPU, so a multi-GPU node reconstructs
%    several shifts at once. All results (metrics + overlay projections + CONFIG)
%    are cached to a .mat so both figures can be regenerated without re-running
%    k-Wave.
%
%  This is a fresh k-Wave study: it GENERATES its own recons via
%  run_single_field_simulation. It does not read pre-computed recons.
%
%  FIXED-GEOMETRY DESIGN (why the numbers mean what they mean)
%    - The patient anatomy (CT_1) is fixed, so the acoustic medium is built ONCE
%      and reused for every beam and every shift.
%    - The ultrasound array (determine_sensor_mask_lateral) is placed ONCE at the
%      plan level and FROZEN (passed as precomputed_sensor), so it never chases
%      the translated dose. Only the dose position changes, which is exactly the
%      quantity this study probes.
%    - Each shifted recon is compared to the fixed original truth, so as the dose
%      slides away the eval-mask overlap shrinks and the pass rate falls.
%
%  LATERAL AXIS (flagged assumption)
%    Per CLAUDE.md the MATLAB array is (row,col,slice) = (Y,X,Z), so the lateral
%    (left-right) axis is dim 2 (columns): higher column index = patient LEFT,
%    lower = patient RIGHT (HFS). A "+" signed shift below is toward patient
%    LEFT, "-" toward patient RIGHT. The script prints this each run so the
%    anatomical direction can be sanity-checked.
%
%  NOTE: HIPAA / remote-execution - this file is WRITTEN here but must be RUN on
%  the remote device. Do not execute locally.
%
%  RUNTIME: each shift is a full forward + 10-iteration time-reversal recon, so
%  the full sweep (all beams, both directions) is a multi-hour GPU job even with
%  multiple GPUs. Trim CONFIG.beams for a smaller first pass.
%  =========================================================================

clear; clc; close all;

% Script lives in studies/, one level below the repo root. run_single_field_simulation
% lives at the repo root; the helpers live in utils/ and pipeline/.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'utils')));
addpath(genpath(fullfile(repoRoot, 'pipeline')));
run_timer = tic;

%% ========================= CONFIGURATION ================================

% --- Machine / data selection ---
CONFIG.working_dir      = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id       = '1194203';
CONFIG.session          = 'Session_1';
CONFIG.treatment_site   = 'Pancreas';
CONFIG.gruneisen_method = 'threshold_2';

% --- Segment selection (one random segment per beam, CT_1 only) ---
CONFIG.beams        = 1:17;         % beams to include (trim for a faster run)
CONFIG.plan_type    = 'reference';  % 'reference' | 'adapted' | 'any'
CONFIG.ct_label     = 'CT_1';       % use CT 1 for all calculations
CONFIG.random_seed  = 42;           % reproducible per-beam segment pick

% --- Reconstruction (simple full-access time reversal) ---
CONFIG.num_time_reversal_iter = 10;
CONFIG.reconstruction_method  = 'tr';

% --- Sensor (realistic tilted array, placed once and frozen) ---
% Kept at the pipeline default so the study uses the same array the real recons
% use. If you switch to the ANTERIOR 'determine_sensor_mask', also set
% engine_config.force_turn_angle = 290 (the known anterior-rotation workaround);
% the lateral method used here needs no such override.
CONFIG.sensor_placement_method = 'determine_sensor_mask_lateral';

% --- Resolution ---
CONFIG.downscale_factor = 1;        % 1 = full resolution

% --- Shift sweep ---
CONFIG.base_shift_mm     = 2;       % geometric base: shift_n = base^n mm (n = 0,1,2,...)
CONFIG.body_mask_shifted = true;    % re-intersect each translated dose with the fixed body

% --- Metrics ---
CONFIG.gamma_n   = 3;               % gamma criterion (n %/n mm)
CONFIG.normalize = true;            % lsq-scale each recon to its own source dose

% --- Parallelization (multi-GPU) ---
% One worker is pinned per GPU and the per-shift reconstructions run in a parfor.
% num_workers = [] sizes the pool to gpuDeviceCount. With <=1 GPU the loop runs
% serially on the client (no pool) to avoid device contention.
CONFIG.use_parallel = true;
CONFIG.num_workers  = [];

% --- Overlay figure (FIGURE 2) ---
% Signed shift (mm) whose translated dose is drawn (red) over the true dose
% (green). [] = each beam's largest available leftward (+) shift. A number picks
% the nearest available shift for every beam.
CONFIG.overlay_shift_mm = [];

% --- Output / cache ---
CONFIG.save_results = true;
% When true and a matching results .mat already exists (same sim-affecting
% CONFIG), the whole k-Wave sweep is skipped and both figures are regenerated
% straight from the cache.
CONFIG.use_cache = true;

%% ========================= SETUP ========================================

rng(CONFIG.random_seed);

fprintf('============================================================\n');
fprintf(' STUDY_SIMULATION_RIGID_SHIFTS\n');
fprintf(' Patient %s | %s | plan=%s | %s\n', ...
    CONFIG.patient_id, CONFIG.session, CONFIG.plan_type, CONFIG.ct_label);
fprintf(' Beams: %s\n', mat2str(CONFIG.beams));
fprintf(' Recon: %s, %d iterations | sensor: %s | downscale: %g\n', ...
    CONFIG.reconstruction_method, CONFIG.num_time_reversal_iter, ...
    CONFIG.sensor_placement_method, CONFIG.downscale_factor);
fprintf(' Lateral axis = array dim 2 (X). +mm = patient LEFT, -mm = patient RIGHT (HFS).\n');
fprintf('============================================================\n');

% Engine config: start from the shared defaults, override only the study knobs.
% normalize/plot are forced off: the engine's 'normalize' means max-normalize,
% but this study does its own least-squares normalization at the metric stage.
engine_config = get_default_config();
engine_config.working_dir             = CONFIG.working_dir;
engine_config.patient_id              = CONFIG.patient_id;
engine_config.session                 = CONFIG.session;
engine_config.gruneisen_method        = CONFIG.gruneisen_method;
engine_config.sensor_placement_method = CONFIG.sensor_placement_method;
engine_config.reconstruction_method   = CONFIG.reconstruction_method;
engine_config.num_time_reversal_iter  = CONFIG.num_time_reversal_iter;
engine_config.blind_recon             = false;
engine_config.downscale_factor        = CONFIG.downscale_factor;
engine_config.normalize               = false;
engine_config.plot_results            = false;

processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
    CONFIG.patient_id, CONFIG.session, 'processed');
out_dir = fullfile(CONFIG.working_dir, 'SimulationResults', ...
    CONFIG.patient_id, CONFIG.session, CONFIG.gruneisen_method);
if ~isfolder(out_dir)
    mkdir(out_dir);
end
out_stem = fullfile(out_dir, sprintf('rigid_shift_study_%s_%s', ...
    CONFIG.patient_id, CONFIG.session));
out_mat  = [out_stem, '.mat'];

%% ========================= CACHE CHECK ==================================

loaded_from_cache = false;
if CONFIG.use_cache && isfile(out_mat)
    S = load(out_mat);
    if isfield(S, 'config') && isfield(S, 'results') && cache_is_valid(S.config, CONFIG)
        results    = S.results;
        spacing_mm = S.spacing_mm;
        lat_dim    = S.lat_dim;
        loaded_from_cache = true;
        fprintf('\n[CACHE] Reusing results from %s (%d beams); skipping k-Wave.\n', ...
            out_mat, numel(results));
    else
        fprintf('\n[CACHE] Existing %s does not match this CONFIG; recomputing.\n', out_mat);
    end
end

%% ===================== COMPUTE (skipped on a cache hit) =================

if ~loaded_from_cache

% ---- Segment selection: one random segment per beam. ----
[field_index, ~] = list_processed_field_doses(CONFIG.patient_id, CONFIG.session, CONFIG);
chosen = pick_random_segment_per_beam(field_index, CONFIG.beams, ...
    CONFIG.plan_type, CONFIG.ct_label);
if isempty(chosen)
    error('study_simulation_rigid_shifts:NoSegments', ...
        'No segments matched plan_type=%s / ct_label=%s for the requested beams.', ...
        CONFIG.plan_type, CONFIG.ct_label);
end
fprintf('\nSelected %d beam(s):\n', numel(chosen));
for k = 1:numel(chosen)
    fprintf('  Beam %2d -> segment %d (%s)\n', chosen(k).beam, chosen(k).segment, ...
        chosen(k).source_mat_filename);
end

% ---- Shared geometry: CT_1 CBCT, medium (once), frozen plan sensor (once). ----
sct = load_cbct1(processed_dir);
spacing_mm     = sct.spacing(:)';
lat_dim        = 2;                  % lateral = dim 2 (X, columns)
lat_spacing_mm = spacing_mm(lat_dim);
n_lat          = size(sct.bodyMask, lat_dim);
if ~isfield(sct, 'origin') || isempty(sct.origin), sct.origin = [0, 0, 0]; end

beam_metadata = load_beam_metadata(processed_dir);

fprintf('\nBuilding acoustic medium (once, CT_1)...\n');
medium = build_medium_with_bath(sct, engine_config);

fprintf('Placing sensor array (once, frozen for all beams/shifts)...\n');
placement_dose = load_summed_rs_dose(processed_dir, chosen(1).file, size(sct.bodyMask));
precomp = precompute_plan_sensor(sct, placement_dose, beam_metadata, engine_config);
if isempty(precomp)
    fprintf('  Sensor method is deterministic; the engine builds a fixed mask inline.\n');
else
    fprintf('  Frozen sensor mask: %d active points.\n', sum(precomp.sensor_mask(:)));
end

% ---- Parallel pool: one worker pinned per GPU (multi-GPU only). ----
num_gpus = 0;
if engine_config.use_gpu
    try, num_gpus = gpuDeviceCount; catch, num_gpus = 0; end
end
if CONFIG.use_parallel && num_gpus > 1
    nw = CONFIG.num_workers;
    if isempty(nw), nw = num_gpus; end
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= nw
        if ~isempty(pool), delete(pool); end
        pool = parpool('local', nw);
    end
    spmd
        gpuDevice(mod(labindex - 1, num_gpus) + 1);   % pin each worker to a GPU
    end
    par_M = pool.NumWorkers;
    fprintf('[PARALLEL] Pool of %d worker(s) across %d GPU(s).\n', pool.NumWorkers, num_gpus);
else
    par_M = 0;   % serial on the client (single GPU or parallel disabled)
    fprintf('[SERIAL] Running reconstructions serially (%d GPU detected).\n', num_gpus);
end

% Scalars pulled out of CONFIG so the parfor body classifies them cleanly.
gamma_n          = CONFIG.gamma_n;
normalize        = CONFIG.normalize;
body_mask_shifted = CONFIG.body_mask_shifted;

%% ===================== PER-BEAM SHIFT SWEEP =============================

results = struct('beam', {}, 'segment', {}, 'shift_mm', {}, 'shift_vox', {}, ...
    'gamma', {}, 'ssim', {}, 'overlay_truth_mip', {}, 'overlay_shift_mip', {}, ...
    'overlay_shift_mm', {});

for bi = 1:numel(chosen)
    b = chosen(bi).beam;
    fprintf('\n[Beam %d/%d] beam #%d, segment %d\n', bi, numel(chosen), b, chosen(bi).segment);

    % ---- Load the field dose; body-mask it to get the fixed truth. ----
    fd = load_field_dose_file(chosen(bi).file);
    base_dose = double(fd.dose_Gy);                        % rigid source (unmasked)
    if ~isequal(size(base_dose), size(sct.bodyMask))
        warning('study_simulation_rigid_shifts:GridMismatch', ...
            'Beam #%d dose grid %s ~= CBCT_1 grid %s; skipping.', ...
            b, mat2str(size(base_dose)), mat2str(size(sct.bodyMask)));
        continue;
    end
    truth0 = base_dose .* double(sct.bodyMask);            % fixed comparison truth
    if max(truth0(:)) <= 0
        warning('study_simulation_rigid_shifts:EmptyTruth', ...
            'Beam #%d segment %d has no in-body dose; skipping.', b, chosen(bi).segment);
        continue;
    end
    eval_mask = truth0 >= 0.10 * max(truth0(:));           % fixed 10% eval region
    field_template = make_field_template(fd, sct, spacing_mm);

    % ---- Signed shift list: 0, then +/- base^n mm until the centroid exits. ----
    centroid0 = lateral_centroid(base_dose, lat_dim);
    [shift_mm_list, shift_vox_list] = build_shift_list( ...
        centroid0, n_lat, lat_spacing_mm, CONFIG.base_shift_mm);
    nShift = numel(shift_mm_list);
    fprintf('  Centroid at voxel %.1f/%d -> %d shift(s): %s mm\n', ...
        centroid0, n_lat, nShift, mat2str(shift_mm_list));

    gamma_vals = nan(1, nShift);
    ssim_vals  = nan(1, nShift);

    % ---- One forward + 10-iter TR reconstruction per shift, across the GPUs. ----
    parfor (si = 1:nShift, par_M)
        shift_vox = shift_vox_list(si);

        shifted_src = translate_volume(base_dose, shift_vox, lat_dim);
        if body_mask_shifted
            shifted_src = shifted_src .* double(sct.bodyMask);
        end
        field_shift = field_template;
        field_shift.dose_Gy = shifted_src;

        try
            recon = run_recon_quietly(field_shift, sct, medium, ...
                beam_metadata, engine_config, precomp);
            if normalize
                recon = recon * least_squares_gain(shifted_src, recon);
            end
            gamma_vals(si) = gamma_pass_rate(truth0, recon, eval_mask, gamma_n, spacing_mm);
            [~, ssim_mean] = compute_local_ssim(truth0, recon, eval_mask);
            ssim_vals(si)  = 100 * ssim_mean;
        catch ME
            warning('study_simulation_rigid_shifts:ShiftFailed', ...
                'Beam #%d shift %+g mm failed (%s); recording NaN.', ...
                b, shift_mm_list(si), ME.message);
        end
    end

    for si = 1:nShift
        fprintf('    shift %+7.1f mm (%+d vox): gamma %5.1f%% | SSIM %5.1f%%\n', ...
            shift_mm_list(si), shift_vox_list(si), gamma_vals(si), ssim_vals(si));
    end

    % ---- Overlay projections for FIGURE 2 (cheap; no reconstruction). ----
    ov_mm  = choose_overlay_shift(shift_mm_list, CONFIG.overlay_shift_mm);
    ov_vox = shift_vox_list(find(shift_mm_list == ov_mm, 1));
    ov_src = translate_volume(base_dose, ov_vox, lat_dim);
    if CONFIG.body_mask_shifted
        ov_src = ov_src .* double(sct.bodyMask);
    end

    results(end+1) = struct('beam', b, 'segment', chosen(bi).segment, ...
        'shift_mm', shift_mm_list, 'shift_vox', shift_vox_list, ...
        'gamma', gamma_vals, 'ssim', ssim_vals, ...
        'overlay_truth_mip', coronal_mip(truth0), ...
        'overlay_shift_mip', coronal_mip(ov_src), ...
        'overlay_shift_mm', ov_mm); %#ok<SAGROW>

    fprintf('  Beam #%d done | elapsed %.0f s\n', b, toc(run_timer));
end

if isempty(results)
    error('study_simulation_rigid_shifts:NoResults', ...
        'No beams produced results (all skipped).');
end

end   % if ~loaded_from_cache

%% ========================= PLOTS (two figures) =========================

fig_metrics = plot_shift_sensitivity(results, CONFIG);   % FIGURE 1
fig_overlay = plot_overlay_grid(results, CONFIG);        % FIGURE 2 (separate)

%% ========================= SAVE RESULTS ================================

if CONFIG.save_results && ~loaded_from_cache
    RESULTS = struct();
    RESULTS.config          = CONFIG;
    RESULTS.results         = results;
    RESULTS.spacing_mm      = spacing_mm;
    RESULTS.lat_dim         = lat_dim;
    RESULTS.total_runtime_s = toc(run_timer);

    save(out_mat, '-struct', 'RESULTS', '-v7.3');
    fprintf('\nResults saved to: %s\n', out_mat);

    savefig(fig_metrics, [out_stem, '_metrics.fig']);
    saveas(fig_metrics,  [out_stem, '_metrics.png']);
    savefig(fig_overlay, [out_stem, '_overlay.fig']);
    saveas(fig_overlay,  [out_stem, '_overlay.png']);
    fprintf('Figures saved to: %s_{metrics,overlay}.{fig,png}\n', out_stem);
end

fprintf('\nTotal runtime: %.1f s (%.2f min) | %d beam(s).\n', ...
    toc(run_timer), toc(run_timer)/60, numel(results));


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function recon = run_recon_quietly(field_shift, sct, medium, beam_metadata, engine_config, precomp)
%RUN_RECON_QUIETLY One k-Wave forward + time-reversal recon, banners silenced.
%  Wrapping the evalc in a named function (rather than a bare string in the
%  parfor body) keeps every input a normal argument, so parfor still ships them
%  to the workers (mirrors pipeline_simulate's run_field_quietly). Gathers the
%  result off the GPU so the caller gets a plain CPU double.
    [~, recon] = evalc(['run_single_field_simulation(field_shift, sct, ', ...
        'medium, beam_metadata, engine_config, precomp)']);
    recon = double(gather(recon));
end

function chosen = pick_random_segment_per_beam(field_index, beams, plan_type, ct_label)
%PICK_RANDOM_SEGMENT_PER_BEAM One random matching segment per requested beam.
%  Filters the dose-file index by plan type and CT label (parsed from the
%  filename, same token grammar as load_recon_dose_data) and picks one segment
%  per beam using the already-seeded RNG. Beams with no match are skipped.
    chosen = struct('beam', {}, 'segment', {}, 'file', {}, 'source_mat_filename', {});
    for b = beams(:)'
        cand = [];
        for i = 1:numel(field_index)
            if field_index(i).beam_index ~= b, continue; end
            [pt, ct] = parse_dose_tokens(field_index(i).source_mat_filename);
            if ~strcmpi(plan_type, 'any') && ~strcmpi(pt, plan_type), continue; end
            if ~strcmpi(ct_label, 'any') && ~strcmpi(ct, ct_label), continue; end
            cand(end+1) = i; %#ok<AGROW>
        end
        if isempty(cand)
            warning('study_simulation_rigid_shifts:NoBeamMatch', ...
                'No segment matched beam #%d; skipping.', b);
            continue;
        end
        pick = cand(randi(numel(cand)));
        chosen(end+1) = struct('beam', b, ...
            'segment', field_index(pick).segment, ...
            'file', field_index(pick).file, ...
            'source_mat_filename', field_index(pick).source_mat_filename); %#ok<AGROW>
    end
end

function [plan_type, ct_label] = parse_dose_tokens(name)
%PARSE_DOSE_TOKENS Plan type + CT label from a dose_*.mat filename ('' if absent).
    plan_type = '';
    ct_label  = '';
    tok = regexp(char(name), ...
        '_(adapted|reference)(?:_CT_(\d+))?_B\d+_\d+\.mat$', 'tokens', 'once');
    if ~isempty(tok)
        plan_type = tok{1};
        if numel(tok) >= 2 && ~isempty(tok{2})
            ct_label = sprintf('CT_%s', tok{2});
        end
    end
end

function sct = load_cbct1(processed_dir)
%LOAD_CBCT1 Load the resampled CT_1 geometry (CBCT1_resampled) from the processed dir.
    f = fullfile(processed_dir, 'CBCT1_resampled.mat');
    if ~isfile(f)
        error('study_simulation_rigid_shifts:NoCBCT1', 'CBCT1_resampled.mat not found: %s', f);
    end
    L = load(f, 'CBCT1_resampled');
    if ~isfield(L, 'CBCT1_resampled')
        error('study_simulation_rigid_shifts:NoCBCT1Struct', ...
            'CBCT1_resampled not present in %s', f);
    end
    sct = L.CBCT1_resampled;
    if ~isfield(sct, 'cubeHU') || ~isfield(sct, 'bodyMask')
        error('study_simulation_rigid_shifts:BadCBCT1', ...
            'CBCT1_resampled missing cubeHU/bodyMask.');
    end
    if ~isfield(sct, 'couchMask') || isempty(sct.couchMask)
        sct.couchMask = false(size(sct.bodyMask));
    end
end

function beam_metadata = load_beam_metadata(processed_dir)
%LOAD_BEAM_METADATA Plan beam metadata from metadata.mat ([] if unavailable).
    beam_metadata = [];
    f = fullfile(processed_dir, 'metadata.mat');
    if ~isfile(f), return; end
    try
        md = load(f, 'metadata');
        if isfield(md, 'metadata') && isfield(md.metadata, 'beam_metadata')
            beam_metadata = md.metadata.beam_metadata;
        end
    catch ME
        warning('study_simulation_rigid_shifts:MetadataLoad', ...
            'Failed to load metadata.mat: %s', ME.message);
    end
end

function d = load_summed_rs_dose(processed_dir, fallback_file, target_dims)
%LOAD_SUMMED_RS_DOSE Summed CT_1 RayStation dose for the sensor placement.
%  Uses total_rs_dose.mat (sparse-aware) when present; otherwise falls back to
%  the first selected segment's dose so the frozen sensor still has a placement
%  dose. Falls back too if the summed grid does not match the CBCT grid.
    d = [];
    f = fullfile(processed_dir, 'total_rs_dose.mat');
    if isfile(f)
        L = load(f);
        if isfield(L, 'total_rs_dose_sparse') && isfield(L, 'total_rs_dose_dims')
            d = reshape(full(L.total_rs_dose_sparse), L.total_rs_dose_dims);
        elseif isfield(L, 'total_rs_dose')
            d = L.total_rs_dose;
        end
    end
    if ~isempty(d) && ~isequal(size(d), target_dims)
        warning('study_simulation_rigid_shifts:SummedGridMismatch', ...
            'total_rs_dose grid %s ~= CBCT grid %s; using a segment dose for placement.', ...
            mat2str(size(d)), mat2str(target_dims));
        d = [];
    end
    if isempty(d)
        if ~isfile(f)
            warning('study_simulation_rigid_shifts:NoSummedDose', ...
                'total_rs_dose.mat not found; using a segment dose for sensor placement.');
        end
        fd = load_field_dose_file(fallback_file);
        d  = fd.dose_Gy;
    end
    d = double(d);
end

function precomp = precompute_plan_sensor(sct, placement_dose, beam_metadata, config)
%PRECOMPUTE_PLAN_SENSOR Build ONE frozen sensor mask (mirrors pipeline_simulate).
%  Only the determine_sensor_mask[_lateral] methods are precomputed: their
%  placement is deterministic from the placement dose, so it is computed once and
%  reused for every beam/shift (passed to the engine as precomputed_sensor).
%  Returns [] for deterministic grid-plane methods, which the engine builds inline.
    precomp = [];
    method = config.sensor_placement_method;
    if ~ismember(method, {'determine_sensor_mask', 'determine_sensor_mask_lateral'})
        return;
    end

    % SCT-like + field-dose structs the placement routines expect (mirrors the
    % inline construction in run_single_field_simulation / compute_plan_sensor_mask).
    sct_s = struct('cubeHU', sct.cubeHU, 'bodyMask', sct.bodyMask, ...
        'couchMask', sct.couchMask, 'origin', sct.origin, 'spacing', sct.spacing(:)');

    fd = struct();
    fd.dose_Gy       = placement_dose;
    fd.total_dose_Gy = placement_dose;
    fd.gantry_angle  = 0;
    fd.origin        = sct.origin;
    fd.spacing       = sct.spacing(:)';
    fd.dimensions    = size(sct.bodyMask);

    if strcmp(method, 'determine_sensor_mask_lateral')
        [sensor_mask, sensor_info] = determine_sensor_mask_lateral(sct_s, fd, beam_metadata, config);
    else
        [sensor_mask, sensor_info] = determine_sensor_mask(sct_s, fd, beam_metadata, config);
    end

    precomp = struct();
    precomp.sensor_mask     = sensor_mask;
    precomp.sensor_info     = sensor_info;
    precomp.ct_label        = 'CT_1';
    precomp.base_dimensions = size(sct.bodyMask);
    precomp.config_hash     = '';
end

function tmpl = make_field_template(fd, sct, spacing_mm)
%MAKE_FIELD_TEMPLATE field_dose struct for the engine, minus dose_Gy.
%  Carries the meterset / gantry / jaw / isocenter / geometry the engine and the
%  sensor exclusion zones need; the per-shift dose is assigned by the caller.
    tmpl = struct();
    tmpl.dose_Gy      = [];   % filled per shift
    tmpl.spacing      = spacing_mm;
    tmpl.dimensions   = size(sct.bodyMask);
    tmpl.origin       = pick_field(fd, 'origin', sct.origin);
    tmpl.gantry_angle = pick_field(fd, 'gantry_angle', 0);
    ms = pick_field(fd, 'meterset', 100);
    if isempty(ms) || ms <= 0, ms = 100; end
    tmpl.meterset = ms;
    for f = {'isocenter', 'jaw_x', 'jaw_y'}
        if isfield(fd, f{1}) && ~isempty(fd.(f{1}))
            tmpl.(f{1}) = fd.(f{1});
        end
    end
end

function [shift_mm_list, shift_vox_list] = build_shift_list(centroid0, n_lat, lat_spacing_mm, base_mm)
%BUILD_SHIFT_LIST Signed shifts: 0, then +/- base^n mm (n=0,1,2,...).
%  Each direction stops when the rigid dose centroid (centroid0 + shift, in
%  voxels) would leave [1, n_lat]. + = toward higher index (patient left).
    shift_mm_list  = 0;
    shift_vox_list = 0;
    for s = [1, -1]                       % +1 = left, -1 = right
        n = 0;
        while n <= 20                     % safety cap (base^20 mm is enormous)
            mm        = base_mm ^ n;
            shift_vox = s * round(mm / lat_spacing_mm);
            n = n + 1;
            if shift_vox == 0, continue; end          % sub-voxel: skip to next n
            if (centroid0 + shift_vox) < 1 || (centroid0 + shift_vox) > n_lat
                break;                                % centroid off grid: done
            end
            shift_mm_list(end+1)  = s * mm;           %#ok<AGROW>
            shift_vox_list(end+1) = shift_vox;        %#ok<AGROW>
        end
    end
end

function out = translate_volume(vol, shift, dim)
%TRANSLATE_VOLUME Integer shift of vol along dim with ZERO fill (no wraparound).
%  Positive shift moves content toward higher indices. Content pushed off the
%  grid is discarded (dose that leaves the patient is lost, as intended).
    out = zeros(size(vol), 'like', vol);
    n = size(vol, dim);
    if abs(shift) >= n
        return;                           % everything shifted off the grid
    end
    idx_src = max(1, 1 - shift) : min(n, n - shift);
    idx_dst = idx_src + shift;
    subs_src = repmat({':'}, 1, ndims(vol));
    subs_dst = repmat({':'}, 1, ndims(vol));
    subs_src{dim} = idx_src;
    subs_dst{dim} = idx_dst;
    out(subs_dst{:}) = vol(subs_src{:});
end

function c = lateral_centroid(dose, dim)
%LATERAL_CENTROID Dose-weighted centroid (voxel index) along dim. NaN if no dose.
    w = double(dose);
    w(w < 0) = 0;
    otherdims = setdiff(1:ndims(dose), dim);
    profile = w;
    for d = otherdims
        profile = sum(profile, d);
    end
    profile = profile(:);
    total = sum(profile);
    if total <= 0
        c = NaN;
        return;
    end
    c = sum((1:numel(profile))' .* profile) / total;
end

function mm = choose_overlay_shift(shift_mm_list, requested)
%CHOOSE_OVERLAY_SHIFT Signed shift (mm) to draw in FIGURE 2 for this beam.
%  requested = [] -> the largest available leftward (+) shift, else the largest
%  magnitude; a number -> the nearest available shift.
    if isempty(requested)
        pos = shift_mm_list(shift_mm_list > 0);
        if ~isempty(pos)
            mm = max(pos);
        else
            [~, ix] = max(abs(shift_mm_list));
            mm = shift_mm_list(ix);
        end
    else
        [~, ix] = min(abs(shift_mm_list - requested));
        mm = shift_mm_list(ix);
    end
end

function mip = coronal_mip(vol)
%CORONAL_MIP Coronal max-intensity projection (project along dim 1 = Y).
%  Returns a single-precision [Z x X] image with superior at the top and the
%  lateral (X) axis horizontal, so a lateral shift shows as horizontal motion.
    m   = squeeze(max(double(vol), [], 1));   % [X x Z]
    mip = single(flipud(m'));                 % [Z x X], superior up
end

function p = gamma_pass_rate(ref, tgt, mask, crit, spacing)
%GAMMA_PASS_RATE Global gamma pass rate (%) over mask, CalcGamma output silenced.
%  Same call as study_pass_rates_allsegments: global gamma (local,0) at crit%/crit
%  mm, DTA search capped at 2*crit, restrict on, forced onto the CPU.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(ref));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', double(tgt));
    gmap = []; %#ok<NASGU>
    try
        evalc(['gmap = CalcGamma(ref_struct, tgt_struct, crit, crit, ', ...
               '''local'', 0, ''limit'', crit*2, ''restrict'', 1, ''cpu'', 1);']);
        if any(mask(:))
            p = 100 * mean(gmap(mask) <= 1);
        else
            p = NaN;
        end
    catch
        p = NaN;
    end
end

function tf = cache_is_valid(saved, CONFIG)
%CACHE_IS_VALID True when a saved run's CONFIG matches the sim-affecting fields.
%  Overlay-only knobs (overlay_shift_mm) are intentionally excluded: they change
%  the picture, not the reconstructions, but the overlay projections are baked
%  into the cache at their chosen shift -- clear the .mat to redraw a new shift.
    tf = false;
    keys = {'patient_id', 'session', 'beams', 'plan_type', 'ct_label', ...
        'random_seed', 'num_time_reversal_iter', 'reconstruction_method', ...
        'sensor_placement_method', 'downscale_factor', 'base_shift_mm', ...
        'body_mask_shifted', 'gamma_n', 'normalize'};
    for i = 1:numel(keys)
        k = keys{i};
        if ~isfield(saved, k) || ~isfield(CONFIG, k) || ~isequal(saved.(k), CONFIG.(k))
            return;
        end
    end
    tf = true;
end

function v = pick_field(s, f, d)
%PICK_FIELD s.(f) when present and non-empty, else default d.
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function fig = plot_shift_sensitivity(results, CONFIG)
%PLOT_SHIFT_SENSITIVITY FIGURE 1: gamma and SSIM vs signed shift, one line/beam.
%  Two panels; x = 0 baseline; a 90% reference line on the gamma panel.
%  + shift = patient left, - shift = patient right.
    fig = figure('Name', 'Rigid lateral-shift sensitivity', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1200, 520]);

    ax1 = subplot(1, 2, 1); hold(ax1, 'on'); grid(ax1, 'on');
    ax2 = subplot(1, 2, 2); hold(ax2, 'on'); grid(ax2, 'on');

    colors = lines(numel(results));
    for bi = 1:numel(results)
        R = results(bi);
        [xs, ord] = sort(R.shift_mm);
        lbl = sprintf('Beam %d (seg %d)', R.beam, R.segment);
        plot(ax1, xs, R.gamma(ord), '-o', 'Color', colors(bi, :), ...
            'MarkerFaceColor', colors(bi, :), 'DisplayName', lbl);
        plot(ax2, xs, R.ssim(ord),  '-o', 'Color', colors(bi, :), ...
            'MarkerFaceColor', colors(bi, :), 'DisplayName', lbl);
    end

    yline(ax1, 90, '--r', '90%');
    xline(ax1, 0, ':', 'Color', [0.5 0.5 0.5]);
    xline(ax2, 0, ':', 'Color', [0.5 0.5 0.5]);

    xlabel(ax1, 'Lateral shift (mm)   [+ = patient LEFT]');
    ylabel(ax1, sprintf('Gamma pass rate (%%)  @ %g%%/%g mm', CONFIG.gamma_n, CONFIG.gamma_n));
    title(ax1, 'Gamma vs rigid lateral shift');
    ylim(ax1, [0, 100]);

    xlabel(ax2, 'Lateral shift (mm)   [+ = patient LEFT]');
    ylabel(ax2, 'Mean local SSIM (%)');
    title(ax2, 'SSIM vs rigid lateral shift');

    legend(ax1, 'show', 'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 8);
    sgtitle(sprintf('Rigid lateral-shift sensitivity  |  %s / %s  |  %s, %d-iter TR', ...
        CONFIG.patient_id, CONFIG.session, CONFIG.reconstruction_method, ...
        CONFIG.num_time_reversal_iter), 'FontWeight', 'bold');
    drawnow;
end

function fig = plot_overlay_grid(results, CONFIG)
%PLOT_OVERLAY_GRID FIGURE 2 (separate): true dose (green) over translated dose
%  (red), coronal max-projection, one panel per beam in a single figure. Overlap
%  reads yellow, so the green/red separation is the rigid lateral shift.
    nB = numel(results);
    ncols = min(6, nB);
    nrows = ceil(nB / ncols);

    fig = figure('Name', 'Truth (green) vs translated dose (red)', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [120, 120, 240 * ncols + 80, 240 * nrows + 80]);

    for bi = 1:nB
        R = results(bi);
        ax = subplot(nrows, ncols, bi);
        image(ax, overlay_rgb(R.overlay_truth_mip, R.overlay_shift_mip));
        axis(ax, 'image'); set(ax, 'XTick', [], 'YTick', []);
        title(ax, sprintf('Beam %d  (%+g mm)', R.beam, R.overlay_shift_mm), 'FontSize', 8);
    end

    sgtitle(sprintf(['Truth (green) vs translated dose (red), coronal MIP  |  %s / %s' ...
        '   [X horizontal: + = patient LEFT]'], CONFIG.patient_id, CONFIG.session), ...
        'FontWeight', 'bold');
    drawnow;
end

function rgb = overlay_rgb(truth_mip, shift_mip)
%OVERLAY_RGB RGB composite: red = translated dose, green = truth (each max-normed).
    t = double(truth_mip);
    s = double(shift_mip);
    if max(t(:)) > 0, t = t / max(t(:)); end
    if max(s(:)) > 0, s = s / max(s(:)); end
    rgb = zeros([size(t), 3]);
    rgb(:, :, 1) = s;   % red   = translated dose
    rgb(:, :, 2) = t;   % green = truth
end

function medium = build_medium_with_bath(sct, config)
%BUILD_MEDIUM_WITH_BATH create_acoustic_medium + force the coupling bath (outside
%  body / couch) to the uniform medium. Inlined copy of the sanctioned helper in
%  run_standalone_field.m / noise_ensemble_error_bars.m (the single canonical
%  medium builder). Whole-grid force-uniform toggles are left to the engine.
    medium = create_acoustic_medium(sct, config);
    ud = pick_field(config, 'uniform_density',     1000);
    uc = pick_field(config, 'uniform_sound_speed', 1540);
    ua = pick_field(config, 'uniform_alpha_coeff', 0);
    ug = pick_field(config, 'uniform_gruneisen',   1.0);
    if isfield(sct, 'bodyMask')
        outside = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask') && ~isempty(sct.couchMask)
            outside = outside | logical(sct.couchMask);
        end
        medium.density(outside)     = ud;
        medium.sound_speed(outside) = uc;
        if numel(medium.alpha_coeff) > 1
            medium.alpha_coeff(outside) = ua;
        end
        medium.gruneisen(outside)   = ug;
    end
end
