function results = step24_medium_null(patient_id, session, config)
%STEP24_MEDIUM_NULL Constant-medium null: each CT_3 dose simulated and reconstructed in the CT_1 medium.
%
%   results = step24_medium_null(patient_id, session, config)
%   results = step24_medium_null()   % standalone: 1194203, Session_1
%
%   PURPOSE:
%   Recon CT_1 and recon CT_3 differ for two reasons: the dose changed, and the
%   CT_3 wave travelled through a different body (the CT_3 medium) before it was
%   reconstructed on CT_1. This step removes the second reason. For every CT_3
%   field it:
%       1. computes the initial pressure of the CT_3 dose on CT_3 tissue,
%          p0 = D3 * Gamma3 * rho3 (the same p0 Step 2 uses for that field),
%       2. forward-simulates that p0 through the CT_1 medium,
%       3. reconstructs on the CT_1 medium (pressure -> dose with Gamma1 * rho1,
%          CT_1 body mask), as the blind CT_3 recon of Step 2 does,
%       4. scores the result against Truth CT_1: global gamma + local SSIM.
%   The "medium-null recon" then differs from recon CT_1 only by the dose (incl.
%   the CT_3 tissue Gamma*rho in p0), and from recon CT_3 only by the
%   propagation medium and body outline. Each also has its own noise draw.
%   Step 3 draws its medium-null figures from these files (figures/medium_null/).
%   This is a grader-side control: it uses CT_3, so it is not a detector test.
%
%   HOW p0 STAYS ON CT_3: run_single_field_simulation builds p0 from the forward
%   medium as D * gruneisen * density. Here the forward medium is the CT_1 medium
%   with its gruneisen replaced by Gamma3 * rho3 / rho1. Then gruneisen * density
%   = Gamma3 * rho3, while density, sound speed and attenuation (all the wave
%   sees) stay CT_1. Blind mode (config.blind_recon) reconstructs and converts
%   pressure to dose with the unmodified CT_1 medium.
%
%   CACHE: one output per CT_3 field, skipped when it is on disk (resumable).
%   Each field is claimed with an atomic <output>.lock folder, so several
%   pipeline instances can share the work. A crashed run leaves its lock folder
%   behind: delete it to retry that field.
%
%   INPUTS:
%       patient_id - Char/string patient ID (e.g. '1194203').
%       session    - Char/string session name (e.g. 'Session_1').
%       config     - Simulation CONFIG, as built in pipeline_simulate
%                    (defaults in []):
%           .working_dir            REQUIRED. Root holding RayStationFiles/.
%           .gruneisen_method       REQUIRED. Tissue model (create_acoustic_medium, hash).
%           .config_hash            [computed] 8-char sim hash (compute_sim_config_hash).
%           .reconstruction_mode    ['segment'] 'segment' | 'beam' (beam totals).
%           .metrics_plan_type      ['reference'] Segment mode: CT_3 fields of
%                                   this plan only ('any' = all), as Step 2.5.
%           .metrics_normalize      [true]  Scale the recon by its least-squares
%                                   gain to the CT_3 truth (same rule as recon CT_3).
%           .gamma_dose_pct         [3.0]   Gamma dose criterion (%).
%           .gamma_dist_mm          [3.0]   Gamma distance criterion (mm).
%           .gamma_dose_cutoff_pct  [10.0]  Eval mask: Truth CT_1 >= this % of its max.
%           .use_parallel           [false] parfor over fields (one GPU per worker).
%           .num_parallel_workers   Pool size when no pool is running.
%           .step24_overwrite       [false] Redo fields already on disk.
%           Plus every field run_single_field_simulation and
%           create_acoustic_medium read. The force_uniform_* overrides are not
%           supported (they would change p0).
%
%   OUTPUTS:
%       results - Struct: .patient_id, .session, .config_hash, .n_fields,
%                 .n_cached, .n_done, .n_locked, .n_failed
%
%       Per CT_3 field, in get_simulation_dir(..., hash) next to the Step 2 recons:
%           <base>_medium_null_<hash>.mat   (<base> = the CT_3 dose file name,
%                                            e.g. dose_..._reference_CT_3_B5_12
%                                            or total_dose_CT_3_B5 in beam mode)
%               recon_dose   medium-null recon (Gy, no gain applied)
%               config_hash, noise_stats (as Step 2)
%               null_scores  .beam_num, .seg_num (0 in beam mode), .plan_type,
%                            .gamma_pass_rate, .ssim_mean (%), .recon_gain
%               null_maps    .mask_idx, .gamma_vals, .ssim_vals (eval mask only),
%                            .vol_size. Re-expand: M = nan(m.vol_size);
%                            M(m.mask_idx) = m.gamma_vals;
%
%   ALGORITHM:
%       1. List the CT_3 fields (segment mode) or CT_3 beam totals (beam mode).
%       2. Load both CBCTs; build the CT_1 medium and the hybrid forward medium.
%       3. Load the Step 2 plan sensor mask.
%       4. For each field not on disk: claim, simulate (blind, recon on CT_1),
%          LS-scale to Truth CT_3, gamma + SSIM vs Truth CT_1, save, release.
%          The first field runs on this client (fills the CT_1 DAS / adjoint
%          caches if missing); the rest run in a parfor when use_parallel.
%
%   EXAMPLE:
%       addpath(genpath('utils')); addpath('pipeline');
%       CONFIG = get_default_config();   % same sim fields (hash) as pipeline_simulate
%       r = step24_medium_null('1194203', 'Session_1', CONFIG);
%
%   DEPENDENCIES:
%       run_single_field_simulation, create_acoustic_medium, k-Wave,
%       list_processed_field_doses, load_field_dose_file, compute_sim_config_hash,
%       get_simulation_dir, least_squares_gain, CalcGamma, compute_local_ssim.
%       Step 2 sensor_mask_<hash>.mat (determine_sensor_mask methods).
%
%   See also: pipeline_simulate, step205_run_full_field_simulation,
%             step25_segment_metrics, step3_analysis

    %% ======================== STANDALONE DEFAULTS ========================
    if nargin == 0
        repo_root = fileparts(fileparts(mfilename('fullpath')));
        addpath(genpath(fullfile(repo_root, 'utils')));
        patient_id = '1194203';
        session    = 'Session_1';
        config     = get_default_config();
        config.plot_results = false;   % no engine figures, as in pipeline_simulate
    end

    %% ======================== INPUT VALIDATION ========================
    if ~ischar(patient_id) && ~isstring(patient_id)
        error('step24_medium_null:InvalidInput', ...
            'patient_id must be a string or character array.');
    end
    patient_id = char(patient_id);
    if ~ischar(session) && ~isstring(session)
        error('step24_medium_null:InvalidInput', ...
            'session must be a string or character array.');
    end
    session = char(session);
    if ~isstruct(config) || ~isfield(config, 'working_dir') || ~isfield(config, 'gruneisen_method')
        error('step24_medium_null:MissingConfig', ...
            'config must be a struct containing working_dir and gruneisen_method.');
    end

    %% ======================== CONFIG DEFAULTS ========================
    if ~isfield(config, 'config_hash') || isempty(config.config_hash)
        config.config_hash = compute_sim_config_hash(config);
    end
    config = default_field(config, 'reconstruction_mode',   'segment');
    config = default_field(config, 'metrics_plan_type',     'reference');
    config = default_field(config, 'metrics_normalize',     true);
    config = default_field(config, 'gamma_dose_pct',        3.0);
    config = default_field(config, 'gamma_dist_mm',         3.0);
    config = default_field(config, 'gamma_dose_cutoff_pct', 10.0);
    config = default_field(config, 'use_parallel',          false);
    config = default_field(config, 'step24_overwrite',      false);
    overrides = {'force_uniform_density', 'force_uniform_sound_speed', ...
                 'force_uniform_attenuation', 'force_uniform_gruneisen'};
    for k = 1:numel(overrides)
        if isfield(config, overrides{k}) && config.(overrides{k})
            error('step24_medium_null:Unsupported', ...
                '%s = true would change p0; the medium null does not support it.', overrides{k});
        end
    end

    hash8         = config.config_hash;
    beam_mode     = strcmp(config.reconstruction_mode, 'beam');
    sim_dir       = get_simulation_dir(config.working_dir, patient_id, session, hash8);
    processed_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session, 'processed');

    fprintf('\n=========================================================\n');
    fprintf('  [STEP 2.4] Medium null: CT_3 dose, forward + recon in the CT_1 medium\n');
    fprintf('  Patient %s | %s | hash %s | mode %s\n', patient_id, session, hash8, ...
        config.reconstruction_mode);
    fprintf('=========================================================\n');

    %% ======================== CT_3 FIELDS ========================
    % One job per CT_3 segment (segment mode) or per CT_3 beam total (beam mode).
    % Truth CT_1 = the same file name with the CT number swapped.
    [field_index, dose_metadata] = list_processed_field_doses(patient_id, session, config);
    if beam_mode || strcmpi(config.metrics_plan_type, 'any')
        ct3_tag = '_CT_3_B';
    else
        ct3_tag = sprintf('_%s_CT_3_B', lower(config.metrics_plan_type));
    end

    jobs = struct('name', {}, 'dose_file', {}, 'meta_file', {}, 'truth1_file', {}, ...
        'out_file', {}, 'beam', {}, 'seg', {}, 'plan_type', {});
    for i = 1:numel(field_index)
        f = field_index(i);
        if ~contains(f.source_mat_filename, ct3_tag), continue; end
        if beam_mode
            if any([jobs.beam] == f.beam_index), continue; end   % one job per beam
            name = sprintf('total_dose_CT_3_B%d.mat', f.beam_index);
            seg  = 0;
            plan = '';
        else
            name = f.source_mat_filename;
            seg  = f.segment;
            tok  = regexp(name, '_(adapted|reference)_CT_', 'tokens', 'once');
            plan = '';
            if ~isempty(tok), plan = tok{1}; end
        end
        [~, stem] = fileparts(name);
        jobs(end + 1) = struct('name', name, ...
            'dose_file',   fullfile(processed_dir, name), ...
            'meta_file',   f.file, ...
            'truth1_file', fullfile(processed_dir, strrep(name, '_CT_3_B', '_CT_1_B')), ...
            'out_file',    fullfile(sim_dir, sprintf('%s_medium_null_%s.mat', stem, hash8)), ...
            'beam', f.beam_index, 'seg', seg, 'plan_type', plan); %#ok<AGROW>
    end

    results             = struct();
    results.patient_id  = patient_id;
    results.session     = session;
    results.config_hash = hash8;
    results.n_fields    = numel(jobs);

    if config.step24_overwrite
        todo = jobs;
    else
        todo = jobs(~isfile({jobs.out_file}));
    end
    n_todo = numel(todo);
    results.n_cached = numel(jobs) - n_todo;
    fprintf('  %d CT_3 field(s): %d on disk, %d to simulate.\n', numel(jobs), results.n_cached, n_todo);
    if n_todo == 0
        [results.n_done, results.n_locked, results.n_failed] = deal(0);
        return;
    end

    %% ======================== MEDIA + SENSOR ========================
    L = load(fullfile(processed_dir, 'CBCT1_resampled.mat'), 'CBCT1_resampled');
    shared.cbct1 = L.CBCT1_resampled;   % recon geometry: body mask, DAS cache name (CT_1)
    L = load(fullfile(processed_dir, 'CBCT3_resampled.mat'), 'CBCT3_resampled');
    medium3 = create_acoustic_medium(L.CBCT3_resampled, config);
    shared.medium1 = create_acoustic_medium(shared.cbct1, config);

    % Forward medium: CT_1 everywhere except gruneisen, chosen so that
    % gruneisen * density (CT_1) = Gamma3 * rho3, i.e. p0 is the CT_3 p0.
    shared.medium_fwd = shared.medium1;
    shared.medium_fwd.gruneisen = medium3.gruneisen .* medium3.density ./ shared.medium1.density;
    clear medium3 L

    sensor_file   = fullfile(sim_dir, sprintf('sensor_mask_%s.mat', hash8));
    shared.sensor = [];
    if isfile(sensor_file)
        L = load(sensor_file, 'precomputed_sensor');
        shared.sensor = L.precomputed_sensor;
        fprintf('  Using precomputed plan sensor: %s\n', sensor_file);
    elseif ismember(config.sensor_placement_method, ...
            {'determine_sensor_mask', 'determine_sensor_mask_lateral'})
        error('step24_medium_null:NoSensorMask', ...
            'Sensor mask not found: %s\nRun Step 2 (pipeline_simulate) with this CONFIG first.', ...
            sensor_file);
    end

    shared.beam_metadata = [];
    if isfield(dose_metadata, 'beam_metadata')
        shared.beam_metadata = dose_metadata.beam_metadata;
    end
    shared.spacing = dose_metadata.spacing(:)';
    shared.hash8   = hash8;

    %% ======================== SIMULATE + SCORE ========================
    % The first field runs here, so a missing CT_1 DAS matrix / adjoint power
    % cache is built once before parallel workers start.
    use_parfor = config.use_parallel && n_todo > 1;
    status     = cell(1, n_todo);
    status{1}  = run_null_field(todo(1), shared, config, false);
    if use_parfor
        if isempty(gcp('nocreate'))
            parpool(config.num_parallel_workers);
        end
        if config.use_gpu && gpuDeviceCount > 0
            num_gpus = gpuDeviceCount;
            spmd
                gpuDevice(mod(labindex - 1, num_gpus) + 1);   % one GPU per worker, round-robin
            end
        end
        rest        = todo(2:end);
        rest_status = cell(1, numel(rest));
        parfor k = 1:numel(rest)
            rest_status{k} = run_null_field(rest(k), shared, config, true);
        end
        status(2:end) = rest_status;
    else
        for k = 2:n_todo
            status{k} = run_null_field(todo(k), shared, config, false);
        end
    end

    results.n_done   = nnz(strcmp(status, 'done'));
    results.n_locked = nnz(strcmp(status, 'locked'));
    results.n_failed = n_todo - results.n_done - results.n_locked;
    fprintf('[STEP 2.4] Done: %d simulated, %d held by another instance, %d failed.\n', ...
        results.n_done, results.n_locked, results.n_failed);
    if results.n_locked > 0
        fprintf('           If no other run is active, delete the *_medium_null_*.mat.lock folders in %s\n', sim_dir);
    end
end


function status = run_null_field(job, shared, config, quiet)
% One CT_3 field: claim it, simulate in the CT_1 medium, score vs Truth CT_1, save.
% Returns 'done', 'locked' or 'failed: <message>'. quiet = true (parfor) hides the
% engine's console output.
    lock_dir = [job.out_file, '.lock'];
    [ok, ~, msgid] = mkdir(lock_dir);
    if ~ok || strcmpi(msgid, 'MATLAB:MKDIR:DirectoryExists')
        status = 'locked';
        return;
    end

    try
        % Field metadata (meterset, gantry) from the segment file. Beam mode
        % swaps in the beam's summed dose.
        fd = load_field_dose_file(job.meta_file);
        if ~strcmp(job.dose_file, job.meta_file)
            fd.dose_Gy = load_dose_volume(job.dose_file);
        end
        truth1 = load_dose_volume(job.truth1_file);
        truth3 = double(fd.dose_Gy);

        sim_config             = config;
        sim_config.blind_recon = true;   % reconstruct + convert on medium1 (true CT_1)
        if quiet
            [~, recon_dose, sim_results] = evalc(['run_single_field_simulation(fd, ' ...
                'shared.cbct1, shared.medium_fwd, shared.beam_metadata, sim_config, ' ...
                'shared.sensor, shared.medium1)']);
        else
            [recon_dose, sim_results] = run_single_field_simulation(fd, shared.cbct1, ...
                shared.medium_fwd, shared.beam_metadata, sim_config, shared.sensor, shared.medium1);
        end
        if isfield(sim_results, 'failure_report')
            error('step24_medium_null:SimFailed', '%s', sim_results.failure_report);
        end
        recon_dose = double(gather(recon_dose));

        % Same gain rule as recon CT_3 in Step 2.5: fit to the dose that made it.
        gain = 1;
        if config.metrics_normalize
            gain = least_squares_gain(truth3, recon_dose);
        end
        target = recon_dose * gain;

        % Eval mask from the reference (Truth CT_1), as Step 2.5.
        cutoff_fr = config.gamma_dose_cutoff_pct / 100;
        if max(truth1(:)) > 0
            mask = truth1 >= cutoff_fr * max(truth1(:));
        else
            mask = target >= cutoff_fr * max(target(:));
        end

        % Global gamma on the CPU, same CalcGamma call as Step 2.5.
        ref_struct = struct('start', [0, 0, 0], 'width', shared.spacing, 'data', truth1);
        tgt_struct = struct('start', [0, 0, 0], 'width', shared.spacing, 'data', target);
        dose_pct   = config.gamma_dose_pct;
        dist_mm    = config.gamma_dist_mm;
        [~, gmap]  = evalc(['CalcGamma(ref_struct, tgt_struct, dose_pct, dist_mm, ' ...
            '''local'', 0, ''restrict'', 1, ''cpu'', 1)']);
        g_vals = single(gmap(mask));

        [smap, s_frac] = compute_local_ssim(truth1, target, mask);
        s_vals = single([]);
        if ~isempty(smap)
            s_vals = single(smap(mask));
        end

        null_scores = struct('beam_num', job.beam, 'seg_num', job.seg, ...
            'plan_type', job.plan_type, ...
            'gamma_pass_rate', 100 * mean(double(g_vals) <= 1), ...
            'ssim_mean', 100 * s_frac, 'recon_gain', gain);
        null_maps = struct('mask_idx', uint32(find(mask)), 'gamma_vals', g_vals, ...
            'ssim_vals', s_vals, 'vol_size', size(truth1)); %#ok<NASGU> saved below

        noise_stats = struct('sensor_signal_peak', NaN, 'noise_amp', NaN, 'snr', NaN);
        names = fieldnames(noise_stats);
        for n = 1:numel(names)
            if isfield(sim_results, names{n}) && ~isempty(sim_results.(names{n}))
                noise_stats.(names{n}) = double(gather(sim_results.(names{n})));
            end
        end

        % Temp file + rename, so a half-written file is never seen as done.
        config_hash = shared.hash8; %#ok<NASGU> saved variable
        tmp_file    = [job.out_file(1:end - 4), '.tmp.mat'];
        save(tmp_file, 'recon_dose', 'config_hash', 'noise_stats', 'null_scores', ...
            'null_maps', '-v7.3');
        movefile(tmp_file, job.out_file, 'f');

        fprintf('  %s: gamma %.1f%% | SSIM %.1f%% vs Truth CT_1\n', job.name, ...
            null_scores.gamma_pass_rate, null_scores.ssim_mean);
        status = 'done';
    catch ME
        fprintf('  %s: FAILED: %s\n', job.name, ME.message);
        status = ['failed: ', ME.message];
    end
    if isfolder(lock_dir), rmdir(lock_dir, 's'); end
end


function dose = load_dose_volume(fpath)
% Dense dose (Gy) from a Step 1.5 file: a field dose_*.mat or a beam total_dose_*.mat.
    [~, name] = fileparts(fpath);
    if startsWith(name, 'total_dose_')
        S = load(fpath);
        if isfield(S, 'beam_total_sparse')
            dose = reshape(full(S.beam_total_sparse), S.beam_total_dims);
        else
            dose = S.beam_total;
        end
    else
        fd   = load_field_dose_file(fpath);
        dose = fd.dose_Gy;
    end
    dose = double(dose);
end


function config = default_field(config, name, value)
    if ~isfield(config, name) || isempty(config.(name))
        config.(name) = value;
    end
end
