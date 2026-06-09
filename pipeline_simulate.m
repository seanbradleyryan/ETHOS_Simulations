%% =========================================================================
%  PIPELINE_SIMULATE.m
%  ETHOS Photoacoustic Pipeline  Simulation and Analysis Steps
%  =========================================================================
%
%  PURPOSE:
%  Runs all steps that depend on RayStation field dose exports.  Must be
%  run AFTER:
%    1. pipeline_setup.m has completed (Steps 0 / 0.5 / 0.6)
%    2. Exploded-segment RTPLAN files have been imported into RayStation
%    3. Field doses have been exported as Plan_Field*_Beam*_B*_S*.dcm into
%       RayStationFiles/<PatientID>/<Session>/
%
%  STEPS EXECUTED:
%    Step 2    k-Wave photoacoustic simulation (forward + time-reversal)
%    Step 3    Gamma analysis, SSIM, and visualization
%
%  PREREQUISITES:
%    - MATLAB R2022a or later
%    - k-Wave Toolbox (http://www.k-wave.org)
%    - Image Processing Toolbox
%    - Parallel Computing Toolbox (optional, for CONFIG.use_parallel)
%    - pipeline_setup.m must have been run successfully
%    - pipeline_compress.m must have been run on the Windows work laptop
%      (processes field doses via Step 1.5 and uploads to cluster)
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: March 2026
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION =================================

% --- Patient and Session Selection ---
CONFIG.patients        = {'1194203'};
CONFIG.sessions        = {'Session_1'};
CONFIG.treatment_site  = 'Pancreas';

% --- Directory Paths ---
CONFIG.working_dir  = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.matrad_path  = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

% --- Gruneisen Parameter Method ---
% Options: 'uniform' | 'threshold_1' (9-tissue) | 'threshold_2' (4-tissue)
CONFIG.gruneisen_method = 'threshold_2';

% --- Acoustic Simulation Parameters ---
CONFIG.dose_per_pulse_cGy       = 0.16;   % cGy per LINAC pulse
CONFIG.pml_size                 = 10;     % PML thickness (voxels)
CONFIG.cfl_number               = 0.3;    % CFL stability criterion
CONFIG.use_gpu                  = true;   % GPU acceleration
CONFIG.num_time_reversal_iter   = 1;      % Time-reversal iterations per field

% --- Sensor Placement ---
% Controls the k-Wave sensor mask geometry used in every per-field simulation.
%   'full_plane_anterior' : Full YZ plane at x = sensor_x_index.
%   'full_plane_lateral'  : Full XZ plane at y = sensor_y_index.
%   'spherical'           : Spherical shell.
CONFIG.sensor_placement_method = 'determine_sensor_mask';
CONFIG.sensor_x_index = 20;   % Used by 'full_plane_anterior'
CONFIG.sensor_y_index = 20;   % Used by 'full_plane_lateral'

% --- Pulse Convolution / Noise / Deconvolution ---
% Mimics a finite transducer impulse response applied to forward sensor data.
% Set convolution_kernel to 0 to disable the entire block.
CONFIG.convolution_kernel  = 0;   % Gaussian sigma in seconds (4 us)
CONFIG.conv_noise_level    = 0.01;   % Noise amplitude as fraction of peak sensor signal
CONFIG.conv_deconv_lambda  = 1e-4;   % Wiener regularization for deconvolution

% --- Gamma Analysis Parameters ---
CONFIG.gamma_dose_pct        = 3.0;    % Dose difference threshold (%)
CONFIG.gamma_dist_mm         = 3.0;    % Distance-to-agreement (mm)
CONFIG.gamma_dose_cutoff_pct = 10.0;   % Ignore voxels below this % of max

% --- SSIM & Visualization Parameters ---
CONFIG.analysis_compute_ssim = true;   % Compute SSIM alongside gamma
CONFIG.analysis_plot_results = true;   % Generate PNG figures
CONFIG.analysis_plot_slices  = 'auto'; % Slice indices or 'auto' (25/50/75th %ile)
CONFIG.sensor_mode           = CONFIG.sensor_placement_method;  % passed to step3

% --- Parallel Processing ---
CONFIG.use_parallel          = true;
CONFIG.num_parallel_workers  = 8;

% --- Multi-Instance Coordination ---
% Multiple copies of this script (e.g. one per remote-desktop session, all
% sharing working_dir) coordinate per-field through lock + status files under
% SimulationResults/<patient>/<session>/<method>/field_status/. Before a field
% is processed its lock is claimed atomically and an 'in_progress' status is
% written so sibling instances skip it; on completion the status is overwritten
% with 'complete'. A lock whose status file is older than stale_claim_minutes
% is assumed orphaned (crashed sibling) and may be overtaken.
CONFIG.stale_claim_minutes   = 120;

% --- Pipeline Control Flags ---
CONFIG.run_step2   = true;    % Step 2  : k-Wave simulation
CONFIG.run_step3   = true;    % Step 3  : Gamma analysis

% --- Define Tissue Property Tables ---
CONFIG.tissue_tables = define_tissue_tables();

%% ========================= INITIALIZATION ================================

fprintf('=========================================================\n');
fprintf('  ETHOS Pipeline  Simulate & Analyse (Steps 2 / 3)\n');
fprintf('=========================================================\n');
fprintf('  Started: %s\n', datetime('now'));
fprintf('  Working directory: %s\n', CONFIG.working_dir);
fprintf('  Gruneisen method: %s\n', CONFIG.gruneisen_method);
fprintf('=========================================================\n\n');

% Add required paths
addpath(genpath(CONFIG.matrad_path));
addpath(genpath(fullfile(CONFIG.working_dir, 'PipelineScripts')));

% Verify k-Wave installation
if ~exist('kWaveGrid', 'file')
    error('k-Wave toolbox not found. Please add k-Wave to the MATLAB path.');
end

% Compute a stable hash of the sim-affecting CONFIG fields.  Every cached
% artifact for this run is keyed on this hash, so flipping any sim-relevant
% knob between runs yields a parallel set of files rather than silently
% reusing stale ones.
[CONFIG_HASH, CONFIG_CANONICAL] = compute_sim_config_hash(CONFIG);
fprintf('  Config hash: %s\n', CONFIG_HASH);
fprintf('=========================================================\n\n');

% Initialize results structure
RESULTS                 = struct();
RESULTS.timestamp       = datetime('now');
RESULTS.config          = CONFIG;
RESULTS.config_hash     = CONFIG_HASH;
RESULTS.config_canonical = CONFIG_CANONICAL;
RESULTS.patients        = struct();

%% ========================= MAIN PROCESSING LOOP ==========================

for p_idx = 1:length(CONFIG.patients)
    patient_id = CONFIG.patients{p_idx};

    for s_idx = 1:length(CONFIG.sessions)
        session = CONFIG.sessions{s_idx};

        fprintf('\n=== Processing: Patient %s, %s ===\n', patient_id, session);

        result_key = make_result_key(patient_id, session);
        RESULTS.patients.(result_key) = init_patient_result(patient_id, session);

        % Open per-session log file (appends across runs)
        ensure_config_registry(CONFIG_HASH, CONFIG_CANONICAL, ...
            patient_id, session, CONFIG);
        log_fid = open_simulation_log(patient_id, session, CONFIG, ...
            CONFIG_HASH, CONFIG_CANONICAL);
        log_msg(log_fid, 'Run started  patient %s, session %s, config %s', ...
            patient_id, session, CONFIG_HASH);

        try

            %% ============================================================
            %  Check RayStation files are present before doing anything
            %% ============================================================
            rs_dir = get_raystation_directory(patient_id, session, CONFIG);
            [rs_files_exist, num_dose_files] = check_raystation_files(rs_dir);

            if ~rs_files_exist
                fprintf('[STEP 1] WAITING: No RayStation field dose files found.\n');
                fprintf('         Expected location: %s\n', rs_dir);
                fprintf('         Please complete the RayStation export, then re-run.\n');
                RESULTS.patients.(result_key).status = 'awaiting_raystation';
                continue;
            end
            fprintf('[STEP 1] Found %d dose file(s) in %s\n', num_dose_files, rs_dir);

            %% ============================================================
            %  STEP 1.5: Load Processed Field Doses (run pipeline_compress.m
            %            on the Windows work laptop before this step)
            %% ============================================================
            fprintf('\n[STEP 1.5] Indexing previously processed field doses...\n');
            % Index the dose files WITHOUT loading the heavy dose arrays. Each
            % field's dose is loaded on demand inside the simulation loop (one
            % dose per worker via load_field_dose_file) so the full set is never
            % broadcast into every parfor worker.
            [field_index, dose_metadata] = ...
                list_processed_field_doses(patient_id, session, CONFIG);

            % Load the per-dose CBCT geometries (CT_1 / CT_3). Each field's
            % simulation uses the CBCT the dose was actually computed on
            % (field_dose.ct_label) -- the SCT is not used anywhere.
            cbct_by_label = load_cbct_resampled(patient_id, session, CONFIG);

            % Extract beam metadata for sensor placement
            if isfield(dose_metadata, 'beam_metadata') && ~isempty(dose_metadata.beam_metadata)
                beam_metadata = dose_metadata.beam_metadata;
                fprintf('           Beam metadata: %d beams loaded.\n', length(beam_metadata));
            else
                beam_metadata = [];
                fprintf('           [WARNING] No beam metadata  legacy sensor placement.\n');
            end

            %% ============================================================
            %  STEP 2: k-Wave Photoacoustic Simulation
            %% ============================================================
            if CONFIG.run_step2
                fprintf('\n[STEP 2] Running k-Wave simulations...\n');
                sim_start = tic;

                % Build one acoustic medium per CBCT geometry (CT_1 / CT_3).
                % Each field picks the medium matching field_dose.ct_label.
                medium_by_label = struct();
                cbct_labels = fieldnames(cbct_by_label);
                for li = 1:numel(cbct_labels)
                    lbl = cbct_labels{li};
                    fprintf('         Building acoustic medium for %s...\n', lbl);
                    medium_by_label.(lbl) = create_acoustic_medium(cbct_by_label.(lbl), CONFIG);
                end

                valid_field_indices = 1:numel(field_index);
                num_fields          = length(valid_field_indices);
                grid_dims           = dose_metadata.dimensions;

                fprintf('         Processing %d field(s)...\n', num_fields);

                % --- Load per-hash cache manifest and cross-check against
                %     reconstruction files on disk for the ACTIVE hash only.
                %     Files belonging to a different config are invisible here.
                cache_manifest = load_cache_manifest(patient_id, session, CONFIG, CONFIG_HASH);
                completed_idxs = get_completed_from_manifest(cache_manifest);

                file_idxs = [];
                for fi = valid_field_indices(:).'
                    expected = expected_recon_path(field_index(fi), ...
                        patient_id, session, CONFIG, CONFIG_HASH);
                    if isfile(expected)
                        file_idxs(end+1) = fi; %#ok<AGROW>
                    end
                end
                completed_idxs = union(completed_idxs, file_idxs);

                pending_idxs = setdiff(valid_field_indices, completed_idxs);
                cached_idxs  = intersect(valid_field_indices, completed_idxs);
                num_pending  = numel(pending_idxs);

                print_cache_status(cache_manifest, cached_idxs, pending_idxs, ...
                    field_index, log_fid, CONFIG_HASH);

                % Ensure the cross-instance status directory exists and read the
                % stale-claim threshold (minutes) used to overtake locks left
                % behind by crashed sibling instances.
                status_dir = get_status_directory(patient_id, session, CONFIG);
                if ~exist(status_dir, 'dir'), mkdir(status_dir); end
                stale_minutes = CONFIG.stale_claim_minutes;

                if CONFIG.use_parallel && num_pending > 1
                    fprintf('         Using parallel processing (parfor)...\n');
                    log_msg(log_fid, 'Parallel mode: %d fields pending', num_pending);
                    pool = gcp('nocreate');
                    if isempty(pool) || pool.NumWorkers ~= CONFIG.num_parallel_workers
                        if ~isempty(pool), delete(pool); end
                        parpool(CONFIG.num_parallel_workers);
                    end
                    % Pre-slice the lightweight per-field metadata so parfor
                    % ships one small struct per iteration. The heavy dose_Gy is
                    % NOT held here -- each worker loads its own dose from disk
                    % (load_field_dose_file) rather than broadcasting all doses.
                    pending_meta = field_index(pending_idxs);

                    field_computed = false(num_pending, 1);
                    field_elapsed  = zeros(num_pending, 1);
                    field_gantry   = nan(num_pending, 1);
                    t_parfor       = tic;
                    parfor f = 1:num_pending
                        fi   = pending_idxs(f);
                        meta = pending_meta(f);
                        recon_path  = expected_recon_path(meta, patient_id, session, CONFIG, CONFIG_HASH);
                        lock_path   = field_lock_path(meta, patient_id, session, CONFIG, CONFIG_HASH);
                        status_path = field_status_path(meta, patient_id, session, CONFIG, CONFIG_HASH);

                        % A sibling instance may have finished this field after
                        % the pending list was built -- skip if so.
                        if isfile(recon_path)
                            fprintf('           Field %d already on disk; skipping.\n', fi);
                            continue;
                        end

                        % Atomically claim the field. If another instance holds a
                        % live lock, leave it to them.
                        [claimed, why] = claim_field(lock_path, status_path, stale_minutes);
                        if ~claimed
                            fprintf('           Field %d held by another instance; skipping.\n', fi);
                            continue;
                        end

                        % Load THIS field's dose from disk (one dose in this
                        % worker, not the whole set). On a load failure release
                        % the claim so a corrupt file neither orphans the lock nor
                        % aborts the patient.
                        try
                            fd = load_field_dose_file(meta.file);
                        catch loadME
                            fprintf('           Field %d load failed (%s); releasing claim.\n', ...
                                fi, loadME.message);
                            if isfile(status_path), delete(status_path);   end
                            if isfolder(lock_path), rmdir(lock_path, 's'); end
                            continue;
                        end
                        field_gantry(f) = fd.gantry_angle;

                        % First thing after load: announce IN PROGRESS so sibling
                        % instances see this field is taken.
                        write_field_status(status_path, 'in_progress', fi, fd.gantry_angle, CONFIG_HASH);
                        fprintf('           Field %d (gantry: %.1f deg) [claim: %s]...\n', ...
                            fi, fd.gantry_angle, why);

                        t_field = tic;
                        [cbct_fd, medium_fd] = select_cbct_for_field(fd, ...
                            cbct_by_label, medium_by_label);
                        [rd, ~] = run_single_field_simulation(fd, cbct_fd, ...
                            medium_fd, beam_metadata, CONFIG);
                        rd = gather(rd);
                        save_field_reconstruction(rd, fd, patient_id, session, CONFIG, CONFIG_HASH);
                        save_simulated_dose(rd, fd, patient_id, session, CONFIG, CONFIG_HASH);

                        % Overwrite IN PROGRESS with COMPLETE.
                        write_field_status(status_path, 'complete', fi, fd.gantry_angle, CONFIG_HASH);
                        field_computed(f) = true;
                        field_elapsed(f)  = toc(t_field);
                    end
                    elapsed_parfor = toc(t_parfor);
                    for f = 1:num_pending
                        if ~field_computed(f), continue; end
                        fi = pending_idxs(f);
                        % Merge the real gantry (read in the worker) into the
                        % lightweight metadata for the manifest entry.
                        mi = field_index(fi);
                        mi.gantry_angle = field_gantry(f);
                        cache_manifest = update_manifest_field(cache_manifest, fi, ...
                            mi, patient_id, session, field_elapsed(f));
                        log_msg(log_fid, 'Field %d complete [%s] (parallel, %.1f s)', ...
                            fi, CONFIG_HASH, field_elapsed(f));
                    end
                    save_cache_manifest(cache_manifest, patient_id, session, CONFIG, CONFIG_HASH);
                    fprintf('         Parallel block done in %.1f s (%d field(s) computed here).\n', ...
                        elapsed_parfor, nnz(field_computed));
                else
                    fprintf('         Using serial processing...\n');
                    log_msg(log_fid, 'Serial mode: %d fields pending', num_pending);
                    for f = 1:num_pending
                        fi          = pending_idxs(f);
                        meta        = field_index(fi);
                        recon_path  = expected_recon_path(meta, patient_id, session, CONFIG, CONFIG_HASH);
                        lock_path   = field_lock_path(meta, patient_id, session, CONFIG, CONFIG_HASH);
                        status_path = field_status_path(meta, patient_id, session, CONFIG, CONFIG_HASH);

                        if isfile(recon_path)
                            fprintf('           Field %d/%d already on disk; skipping.\n', f, num_pending);
                            continue;
                        end

                        [claimed, why] = claim_field(lock_path, status_path, stale_minutes);
                        if ~claimed
                            fprintf('           Field %d/%d held by another instance; skipping.\n', f, num_pending);
                            log_msg(log_fid, 'Field %d held by another instance; skipping.', fi);
                            continue;
                        end

                        % Load this field's dose from disk (one dose at a time).
                        % Release the claim on a load failure so a corrupt file
                        % neither orphans the lock nor aborts the patient.
                        try
                            fd = load_field_dose_file(meta.file);
                        catch loadME
                            fprintf('           Field %d/%d load failed (%s); releasing claim.\n', ...
                                f, num_pending, loadME.message);
                            log_msg(log_fid, 'Field %d load failed: %s', fi, loadME.message);
                            if isfile(status_path), delete(status_path);   end
                            if isfolder(lock_path), rmdir(lock_path, 's'); end
                            continue;
                        end

                        % First thing: announce IN PROGRESS to sibling instances.
                        write_field_status(status_path, 'in_progress', fi, fd.gantry_angle, CONFIG_HASH);
                        log_msg(log_fid, 'Field %d IN PROGRESS (gantry %.1f deg, host %s)', ...
                            fi, fd.gantry_angle, get_hostname());
                        fprintf('           Field %d/%d (gantry: %.1f deg) [claim: %s]...\n', ...
                            f, num_pending, fd.gantry_angle, why);

                        t_field = tic;
                        [cbct_fd, medium_fd] = select_cbct_for_field(fd, ...
                            cbct_by_label, medium_by_label);
                        [recon_dose, ~] = run_single_field_simulation(fd, cbct_fd, ...
                            medium_fd, beam_metadata, CONFIG);
                        save_field_reconstruction(recon_dose, fd, patient_id, session, CONFIG, CONFIG_HASH);
                        save_simulated_dose(recon_dose, fd, patient_id, session, CONFIG, CONFIG_HASH);

                        % Overwrite IN PROGRESS with COMPLETE.
                        write_field_status(status_path, 'complete', fi, fd.gantry_angle, CONFIG_HASH);
                        elapsed_f = toc(t_field);
                        cache_manifest = update_manifest_field(cache_manifest, fi, ...
                            fd, patient_id, session, elapsed_f);
                        save_cache_manifest(cache_manifest, patient_id, session, CONFIG, CONFIG_HASH);
                        log_msg(log_fid, 'Field %d complete [%s] (%.1f s)', ...
                            fi, CONFIG_HASH, elapsed_f);
                        fprintf('             Done in %.1f s\n', elapsed_f);
                    end
                end

                % --- Assemble the total reconstruction from every per-field
                %     recon on disk for the active config hash. This picks up
                %     fields completed by sibling instances, not just this one.
                total_recon  = zeros(grid_dims);
                missing_idxs = [];
                for fi = valid_field_indices(:).'
                    rp = expected_recon_path(field_index(fi), patient_id, session, CONFIG, CONFIG_HASH);
                    if isfile(rp)
                        cd_ = load(rp, 'recon_dose');
                        total_recon = total_recon + cd_.recon_dose;
                    else
                        missing_idxs(end+1) = fi; %#ok<AGROW>
                    end
                end

                % If any field is still unprocessed (a sibling instance owns it
                % and hasn't finished), defer total assembly and analysis. Re-run
                % once the siblings complete to pick the remaining fields up.
                if ~isempty(missing_idxs)
                    fprintf('\n[STEP 2] DEFERRED: %d field(s) still pending: %s\n', ...
                        numel(missing_idxs), mat2str(missing_idxs));
                    fprintf('         These are likely in progress on sibling instances.\n');
                    fprintf('         Re-run this script once they finish to assemble the\n');
                    fprintf('         total reconstruction and run the analysis steps.\n');
                    log_msg(log_fid, ['Step 2 deferred: %d field(s) pending on sibling ' ...
                        'instances %s -- total not assembled'], ...
                        numel(missing_idxs), mat2str(missing_idxs));
                    RESULTS.patients.(result_key).status = 'awaiting_siblings';
                    close_simulation_log(log_fid);
                    continue;
                end

                sim_time = toc(sim_start);
                save_total_reconstruction(total_recon, dose_metadata, ...
                    patient_id, session, CONFIG, CONFIG_HASH);

                RESULTS.patients.(result_key).simulation_time_sec  = sim_time;
                RESULTS.patients.(result_key).total_recon_max_Gy   = max(total_recon(:));

                log_msg(log_fid, 'Step 2 complete  %.1f s total, max recon dose %.4f Gy', ...
                    sim_time, max(total_recon(:)));
                fprintf('[STEP 2] Complete (%.1f s). Max recon dose: %.4f Gy\n', ...
                    sim_time, max(total_recon(:)));
            else
                fprintf('\n[STEP 2] Loading previously computed reconstruction...\n');
                [total_recon, ~] = load_total_reconstruction(patient_id, session, CONFIG, CONFIG_HASH);
            end

            %% ============================================================
            %  STEP 3: Gamma Analysis, SSIM & Visualization
            %% ============================================================
            if CONFIG.run_step3
                total_recon = double(gather(total_recon));
                fprintf('\n[STEP 3] Running analysis (gamma + SSIM)...\n');

                step3_results = step3_analysis(patient_id, session, CONFIG);

                RESULTS.patients.(result_key).gamma_ethos_vs_rs_pass_pct = ...
                    step3_results.ethos_vs_rs.gamma.pass_rate;
                RESULTS.patients.(result_key).gamma_rs_vs_recon_pass_pct = ...
                    step3_results.rs_vs_recon.gamma.pass_rate;

                if CONFIG.analysis_compute_ssim && ~isempty(step3_results.ethos_vs_rs.ssim)
                    RESULTS.patients.(result_key).ssim_ethos_vs_rs = ...
                        step3_results.ethos_vs_rs.ssim.ssim_3d;
                    RESULTS.patients.(result_key).ssim_rs_vs_recon = ...
                        step3_results.rs_vs_recon.ssim.ssim_3d;
                end

                fprintf('[STEP 3] Complete.\n');
                fprintf('         ETHOS vs RS  gamma pass: %.1f%%\n', ...
                    step3_results.ethos_vs_rs.gamma.pass_rate);
                fprintf('         RS vs Recon  gamma pass: %.1f%%\n', ...
                    step3_results.rs_vs_recon.gamma.pass_rate);
                if CONFIG.analysis_compute_ssim && ~isempty(step3_results.ethos_vs_rs.ssim)
                    fprintf('         ETHOS vs RS  SSIM:       %.4f\n', ...
                        step3_results.ethos_vs_rs.ssim.ssim_3d);
                    fprintf('         RS vs Recon  SSIM:       %.4f\n', ...
                        step3_results.rs_vs_recon.ssim.ssim_3d);
                end
            end

            RESULTS.patients.(result_key).status = 'complete';
            log_msg(log_fid, 'Patient %s / %s: COMPLETE', patient_id, session);
            close_simulation_log(log_fid);
            fprintf('\n=== %s/%s: SIMULATION COMPLETE ===\n', patient_id, session);

        catch ME
            fprintf('\n[ERROR] %s\n', ME.message);
            for k = 1:length(ME.stack)
                fprintf('        In %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            log_msg(log_fid, 'ERROR: %s', ME.message);
            for k = 1:length(ME.stack)
                log_msg(log_fid, '  at %s (line %d)', ME.stack(k).name, ME.stack(k).line);
            end
            close_simulation_log(log_fid);
            RESULTS.patients.(result_key).status        = 'error';
            RESULTS.patients.(result_key).error.message = ME.message;
            RESULTS.patients.(result_key).error.stack   = ME.stack;
            fprintf('\n=== %s/%s: SIMULATION FAILED ===\n', patient_id, session);
            continue;
        end

    end  % session loop
end  % patient loop

%% ========================= FINALIZE ======================================

fprintf('\n=========================================================\n');
fprintf('  Simulation Pipeline Complete\n');
fprintf('=========================================================\n');

results_filename = sprintf('pipeline_results_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
results_path     = fullfile(CONFIG.working_dir, 'AnalysisResults', results_filename);
if ~exist(fullfile(CONFIG.working_dir, 'AnalysisResults'), 'dir')
    mkdir(fullfile(CONFIG.working_dir, 'AnalysisResults'));
end
save(results_path, 'RESULTS', '-v7.3');
fprintf('  Results saved: %s\n', results_path);

generate_simulation_summary(RESULTS);
fprintf('\n=========================================================\n\n');


%% =========================================================================
%  TISSUE PROPERTY TABLES
%% =========================================================================

function tables = define_tissue_tables()
    % THRESHOLD_1: Detailed 9-tissue segmentation
    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air','Lung','Fat','Water','Blood','Muscle','SoftTissue','Bone','Metal'};
    tables.threshold_1.density       = [1.2,  400,  920,  1000, 1060, 1050, 1040, 1900, 7800];
    tables.threshold_1.sound_speed   = [343,  600,  1450, 1480, 1575, 1580, 1540, 3200, 5900];
    tables.threshold_1.alpha_coeff   = [0,    0.5,  0.48, 0.0022, 0.2, 0.5, 0.5,  10,   0];
    tables.threshold_1.alpha_power   = [1.0,  1.5,  1.5,  2.0,  1.3,  1.0, 1.1,  1.0,  1.0];
    tables.threshold_1.gruneisen     = [0,    0.5,  0.7,  0.11, 0.15, 0.2, 1.0,  0,    0];

    % THRESHOLD_2: Simplified 4-tissue model
    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water','Fat','SoftTissue','Bone'};
    tables.threshold_2.density       = [1000,   920,  1040,        1900];
    tables.threshold_2.sound_speed   = [1480,   1450, 1540,        3200];
    tables.threshold_2.alpha_coeff   = [0.0022, 0.48, 0.5,         10];
    tables.threshold_2.alpha_power   = [2.0,    1.5,  1.1,         1.0];
    tables.threshold_2.gruneisen     = [0.11,   0.7,  1.0,         1.0];

    % UNIFORM: single water-like medium everywhere
    tables.uniform = struct();
    tables.uniform.density      = 1000;
    tables.uniform.sound_speed  = 1540;
    tables.uniform.alpha_coeff  = 0;
    tables.uniform.alpha_power  = 1.1;
    tables.uniform.gruneisen    = 1.0;
end


%% =========================================================================
%  STEP 2 HELPER FUNCTIONS
%% =========================================================================

function cbct_by_label = load_cbct_resampled(patient_id, session, config)
    % Load the per-dose CBCT geometries (CBCT1/CBCT3_resampled) produced by
    % step15_process_doses, keyed by CT label ('CT_1' / 'CT_3'). These are
    % already resampled to the dose grid and carry cubeHU, cubeDensity,
    % bodyMask, couchMask, origin, spacing and dimensions.
    processed_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        patient_id, session, 'processed');

    specs = { 'CBCT1_resampled.mat', 'CBCT1_resampled'; ...
              'CBCT3_resampled.mat', 'CBCT3_resampled' };

    cbct_by_label = struct();
    for k = 1:size(specs, 1)
        fpath = fullfile(processed_dir, specs{k, 1});
        if ~isfile(fpath)
            error('pipeline_simulate:CBCTNotFound', ...
                'Required CBCT geometry file not found: %s', fpath);
        end
        loaded = load(fpath, specs{k, 2});
        cbct   = loaded.(specs{k, 2});
        if isfield(cbct, 'ct_label') && ~isempty(cbct.ct_label)
            lbl = strrep(cbct.ct_label, '-', '_');
        else
            error('pipeline_simulate:CBCTNoLabel', ...
                '%s is missing a ct_label field.', specs{k, 1});
        end
        cbct_by_label.(lbl) = cbct;
        fprintf('           Loaded CBCT geometry %s from %s\n', lbl, specs{k, 1});
    end
end

function [cbct, medium] = select_cbct_for_field(field_dose, cbct_by_label, medium_by_label)
    % Pick the CBCT geometry + acoustic medium matching the CBCT the dose was
    % computed on (field_dose.ct_label).
    if ~isfield(field_dose, 'ct_label') || isempty(field_dose.ct_label)
        error('pipeline_simulate:NoCtLabel', ...
            'Field dose (gantry %.1f deg) has no ct_label; cannot select CBCT.', ...
            field_dose.gantry_angle);
    end
    lbl = strrep(field_dose.ct_label, '-', '_');
    if ~isfield(cbct_by_label, lbl) || ~isfield(medium_by_label, lbl)
        error('pipeline_simulate:CtLabelUnmatched', ...
            'No CBCT/medium loaded for ct_label "%s".', field_dose.ct_label);
    end
    cbct   = cbct_by_label.(lbl);
    medium = medium_by_label.(lbl);
end

function save_field_reconstruction(recon_dose, field_dose, patient_id, session, config, hash8)
    % Save reconstruction as <input_basename>_recon_<hash8>.mat next to other
    % per-config outputs.  Only the hash is embedded; the full config lives
    % in <sim_dir>/config_registry.json.
    out_path    = expected_recon_path(field_dose, patient_id, session, config, hash8);
    config_hash = hash8;  %#ok<NASGU>  saved variable name
    save(out_path, 'recon_dose', 'config_hash', '-v7.3');
end

function save_simulated_dose(recon_dose, field_dose, patient_id, session, config, hash8)
    %SAVE_SIMULATED_DOSE Save per-field reconstruction next to the processed
    %  input dose files, named <input_basename>_sim_<hash8>.mat. Variable inside
    %  is sim_dose. Falls back to a synthetic name if the source filename was
    %  not propagated by load_processed_data.  Only the hash is embedded.
    sim_dose_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        patient_id, session, 'processed', 'simulated_doses');
    if ~exist(sim_dose_dir, 'dir')
        mkdir(sim_dose_dir);
    end

    out_path    = expected_sim_path(field_dose, patient_id, session, config, hash8);
    sim_dose    = recon_dose;  %#ok<NASGU>  saved variable name
    config_hash = hash8;       %#ok<NASGU>  saved variable name
    save(out_path, 'sim_dose', 'config_hash', '-v7.3');
end

function save_total_reconstruction(total_recon, metadata, patient_id, session, config, hash8)
    out_path    = expected_total_recon_path(patient_id, session, config, hash8);
    config_hash = hash8;  %#ok<NASGU>  saved variable name
    save(out_path, 'total_recon', 'metadata', 'config_hash', '-v7.3');
end

function [total_recon, metadata] = load_total_reconstruction(patient_id, session, config, hash8)
    data        = load(expected_total_recon_path(patient_id, session, config, hash8));
    total_recon = data.total_recon;
    metadata    = data.metadata;
end


%% =========================================================================
%  DIRECTORY / PATH HELPERS
%% =========================================================================

function rs_dir = get_raystation_directory(patient_id, session, config)
    rs_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session);
end

function sim_dir = get_simulation_directory(patient_id, session, config)
    sim_dir = fullfile(config.working_dir, 'SimulationResults', patient_id, ...
        session, config.gruneisen_method);
    if ~exist(sim_dir, 'dir')
        mkdir(sim_dir);
    end
end

function [exists, num_files] = check_raystation_files(rs_dir)
    exists    = false;
    num_files = 0;
    if ~exist(rs_dir, 'dir'), return; end
    % Primary: processed .mat files with new naming: dose_[id]_[session]_[adapted|reference]_B[n]_S[m].mat
    processed_dir = fullfile(rs_dir, 'processed');
    if exist(processed_dir, 'dir')
        rd_files = dir(fullfile(processed_dir, 'dose_*.mat'));
        if ~isempty(rd_files)
            num_files = numel(rd_files);
            exists    = true;
            return;
        end
    end
    % Fallback: DICOM files in rs_dir
    rd_files = dir(fullfile(rs_dir, 'dose_*.dcm'));
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'Plan_Field*_Beam*_B*_S*.dcm'));
    end
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'RD.*.dcm'));
    end
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'RD*.dcm'));
    end
    num_files = numel(rd_files);
    exists    = (num_files > 0);
end

function key = make_result_key(patient_id, session)
    key = sprintf('P%s_%s', patient_id, strrep(session, '_', ''));
end

function result = init_patient_result(patient_id, session)
    result            = struct();
    result.patient_id = patient_id;
    result.session    = session;
    result.status     = 'started';
    result.start_time = datetime('now');
end

%% =========================================================================
%  LOGGING HELPERS
%% =========================================================================

function log_fid = open_simulation_log(patient_id, session, config, hash8, canonical)
    % Opens (or creates) a per-session append-mode log file and writes a
    % header that includes the active config hash and a dump of every
    % sim-affecting CONFIG field.
    log_dir  = get_simulation_directory(patient_id, session, config);
    log_path = fullfile(log_dir, 'simulation_log.txt');
    log_fid  = fopen(log_path, 'a');
    if log_fid == -1
        warning('pipeline_simulate:logOpen', 'Could not open log file: %s', log_path);
        log_fid = -1;
        return;
    end
    fprintf(log_fid, '\n========================================\n');
    fprintf(log_fid, 'Run started : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(log_fid, 'Patient     : %s\n', patient_id);
    fprintf(log_fid, 'Session     : %s\n', session);
    if nargin >= 4 && ~isempty(hash8)
        fprintf(log_fid, 'Config hash : %s\n', hash8);
    end
    if nargin >= 5 && isstruct(canonical)
        fprintf(log_fid, 'Config used :\n');
        names = fieldnames(canonical);
        for i = 1:numel(names)
            fprintf(log_fid, '  %-32s = %s\n', names{i}, ...
                format_config_value(canonical.(names{i})));
        end
    end
    fprintf(log_fid, '========================================\n');
end

function s = format_config_value(v)
    % Stringify a canonical config value for the log header.
    if ischar(v) || isstring(v)
        s = char(v);
    elseif islogical(v)
        if v, s = 'true'; else, s = 'false'; end
    elseif isnumeric(v) && isscalar(v)
        if v == fix(v) && abs(v) < 1e9
            s = sprintf('%d', v);
        else
            s = sprintf('%.6g', v);
        end
    elseif isnumeric(v)
        s = mat2str(v);
    else
        try
            s = jsonencode(v);
        catch
            s = '<unprintable>';
        end
    end
end

function log_msg(log_fid, fmt, varargin)
    % Writes a timestamped line to the log file (and only the log file).
    if log_fid < 1, return; end
    msg = sprintf(fmt, varargin{:});
    fprintf(log_fid, '[%s] %s\n', datestr(now, 'HH:MM:SS'), msg);
end

function close_simulation_log(log_fid)
    if log_fid < 1, return; end
    fprintf(log_fid, 'Run ended   : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(log_fid, '========================================\n');
    fclose(log_fid);
end


%% =========================================================================
%  CROSS-INSTANCE CLAIM / STATUS HELPERS
%% =========================================================================
%  Let several copies of this script (one per remote-desktop session, sharing
%  working_dir) divide the per-field work without stepping on each other. Each
%  field is gated by an atomically-created lock directory and described by a
%  human-readable .status file that any instance can read.

function [claimed, reason] = claim_field(lock_path, status_path, stale_minutes)
    % Atomically claim a field using directory creation, which is atomic on
    % POSIX/NFS. Returns claimed=true only if THIS process created the lock.
    % If a sibling holds the lock but its status file is older than
    % stale_minutes (e.g. the sibling crashed mid-field), the lock is overtaken
    % so the field is not orphaned forever.
    [ok, ~, msgid] = mkdir(lock_path);
    if ok && ~strcmpi(msgid, 'MATLAB:MKDIR:DirectoryExists')
        claimed = true; reason = 'new';
        return;
    end
    claimed = false; reason = 'held';
    if stale_minutes > 0 && isfile(status_path)
        info = dir(status_path);
        if ~isempty(info)
            age_min = (now - info(1).datenum) * 24 * 60;
            if age_min > stale_minutes
                claimed = true; reason = 'stale_overtake';
            end
        end
    end
end

function write_field_status(status_path, state, field_idx, gantry, hash8)
    % Overwrite the per-field status file with the current state
    % ('in_progress' | 'complete'). This file is the cross-instance signal
    % board sibling processes read to see what is taken and what is done.
    fid = fopen(status_path, 'w');
    if fid < 0, return; end
    fprintf(fid, 'status  = %s\n', state);
    fprintf(fid, 'field   = %d\n', field_idx);
    fprintf(fid, 'gantry  = %.1f\n', gantry);
    fprintf(fid, 'config  = %s\n', hash8);
    fprintf(fid, 'host    = %s\n', get_hostname());
    fprintf(fid, 'pid     = %d\n', feature('getpid'));
    fprintf(fid, 'updated = %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fclose(fid);
end

function h = get_hostname()
    % Best-effort machine name so a status file says which session owns a field.
    h = getenv('HOSTNAME');
    if isempty(h), h = getenv('COMPUTERNAME'); end
    if isempty(h)
        try
            [st, out] = system('hostname');
            if st == 0, h = strtrim(out); end
        catch
            h = '';
        end
    end
    if isempty(h), h = 'unknown'; end
end


%% =========================================================================
%  CACHE MANIFEST HELPERS
%% =========================================================================

function manifest = load_cache_manifest(patient_id, session, config, hash8)
    % Loads the per-(session,config) cache manifest (or creates a fresh one).
    manifest_path = expected_manifest_path(patient_id, session, config, hash8);
    if exist(manifest_path, 'file')
        data     = load(manifest_path, 'manifest');
        manifest = data.manifest;
        if ~isfield(manifest, 'config_hash') || ~strcmp(manifest.config_hash, hash8)
            % Stale or migrated manifest -- treat as fresh.
            manifest = fresh_manifest(hash8);
        end
        if ~isfield(manifest, 'fields')
            manifest.fields = struct();
        end
    else
        manifest = fresh_manifest(hash8);
    end
end

function manifest = fresh_manifest(hash8)
    manifest              = struct();
    manifest.config_hash  = hash8;
    manifest.fields       = struct();   % sub-struct keyed by 'f_NNN'
    manifest.created      = datetime('now');
    manifest.last_updated = datetime('now');
end

function save_cache_manifest(manifest, patient_id, session, config, hash8)
    manifest.last_updated = datetime('now');
    manifest.config_hash  = hash8;
    save(expected_manifest_path(patient_id, session, config, hash8), 'manifest');
end

function manifest = update_manifest_field(manifest, field_idx, field_dose, ...
        patient_id, session, elapsed_sec)
    % Records a completed field in the manifest.  Stores gantry angle and
    % source basename so the manifest is self-describing.
    key  = sprintf('f_%03d', field_idx);
    base = get_field_basename(field_dose, patient_id, session);
    manifest.fields.(key).field_idx       = field_idx;
    manifest.fields.(key).gantry_angle    = field_dose.gantry_angle;
    manifest.fields.(key).source_basename = base;
    manifest.fields.(key).completed_at    = datetime('now');
    manifest.fields.(key).elapsed_sec     = elapsed_sec;
end

function completed = get_completed_from_manifest(manifest)
    % Returns a numeric array of field indices recorded in the manifest.
    completed = [];
    if ~isfield(manifest, 'fields'), return; end
    keys = fieldnames(manifest.fields);
    for i = 1:numel(keys)
        completed(end+1) = manifest.fields.(keys{i}).field_idx; %#ok<AGROW>
    end
end

function print_cache_status(manifest, cached_idxs, pending_idxs, field_index, log_fid, hash8)
    % Prints a formatted table of cached vs pending fields for the active hash.
    % field_index is the lightweight struct array from list_processed_field_doses;
    % its gantry_angle may be NaN (the real value is read when the dose loads).
    fprintf('\n         --- Cache / Resume Status ---\n');
    fprintf('         Active config hash: %s\n', hash8);
    if ~isempty(cached_idxs)
        fprintf('         Previously completed (this config):\n');
        for i = 1:numel(cached_idxs)
            fi  = cached_idxs(i);
            key = sprintf('f_%03d', fi);
            if isfield(manifest.fields, key)
                e = manifest.fields.(key);
                fprintf('           [DONE] Field %3d  gantry %6.1f deg  completed %s  (%.0f s)\n', ...
                    fi, e.gantry_angle, ...
                    datestr(e.completed_at, 'yyyy-mm-dd HH:MM:SS'), e.elapsed_sec);
                log_msg(log_fid, 'CACHE HIT [%s] field %d (gantry %.1f deg), completed %s', ...
                    hash8, fi, e.gantry_angle, ...
                    datestr(e.completed_at, 'yyyy-mm-dd HH:MM:SS'));
            else
                fprintf('           [DONE] Field %3d  (file on disk, no manifest entry)\n', fi);
                log_msg(log_fid, 'CACHE HIT [%s] field %d (file present, no manifest entry)', ...
                    hash8, fi);
            end
        end
    end
    if ~isempty(pending_idxs)
        fprintf('         Pending:\n');
        for i = 1:numel(pending_idxs)
            fi = pending_idxs(i);
            if ~isempty(field_index) && numel(field_index) >= fi && ...
                    ~isnan(field_index(fi).gantry_angle)
                fprintf('           [TODO] Field %3d  gantry %6.1f deg\n', ...
                    fi, field_index(fi).gantry_angle);
            else
                fprintf('           [TODO] Field %3d\n', fi);
            end
        end
    end
    fprintf('         Total: %d done, %d pending\n', numel(cached_idxs), numel(pending_idxs));
    fprintf('         ----------------------------\n\n');
end


%% =========================================================================
%  SUMMARY
%% =========================================================================

function generate_simulation_summary(results)
    fprintf('\n--- Simulation Summary ---\n');
    pf = fieldnames(results.patients);
    for i = 1:numel(pf)
        p = results.patients.(pf{i});
        fprintf('\n%s / %s:\n', p.patient_id, p.session);
        fprintf('  Status: %s\n', p.status);
        if strcmp(p.status, 'complete')
            if isfield(p, 'num_fields')
                fprintf('  Fields processed: %d\n', p.num_fields);
            end
            if isfield(p, 'simulation_time_sec')
                fprintf('  Simulation time: %.1f s\n', p.simulation_time_sec);
            end
            if isfield(p, 'gamma_ethos_vs_rs_pass_pct')
                fprintf('  Gamma ETHOS vs RS:  %.1f%%\n', p.gamma_ethos_vs_rs_pass_pct);
            end
            if isfield(p, 'gamma_rs_vs_recon_pass_pct')
                fprintf('  Gamma RS vs Recon:  %.1f%%\n', p.gamma_rs_vs_recon_pass_pct);
            end
            if isfield(p, 'ssim_ethos_vs_rs')
                fprintf('  SSIM  ETHOS vs RS:  %.4f\n', p.ssim_ethos_vs_rs);
            end
            if isfield(p, 'ssim_rs_vs_recon')
                fprintf('  SSIM  RS vs Recon:  %.4f\n', p.ssim_rs_vs_recon);
            end
        elseif strcmp(p.status, 'error')
            fprintf('  Error: %s\n', p.error.message);
        end
    end
end


%% =========================================================================
%  HASH->CONFIG REGISTRY
%% =========================================================================
%  The compute_sim_config_hash helper (and its sim_config_fields allow-list)
%  lives in a standalone file so step3_analysis.m and any future consumer
%  can compute the same hash without duplicating the field list.

function key = registry_key(hash_hex)
    % Struct field names cannot start with a digit, so prefix the hash with
    % 'h_' when used as a struct/JSON key.  The bare hash remains in filenames.
    key = ['h_' hash_hex];
end

function p = config_registry_path(patient_id, session, config)
    sim_dir = get_simulation_directory(patient_id, session, config);
    p       = fullfile(sim_dir, 'config_registry.json');
end

function ensure_config_registry(hash_hex, canonical, patient_id, session, config)
    % Append (or refresh) one hash->config entry in
    % <sim_dir>/config_registry.json.  Existing entries are left intact
    % except for last_updated, which is bumped to the current time.
    registry = load_config_registry(patient_id, session, config);
    key      = registry_key(hash_hex);
    now_str  = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    if isfield(registry, key)
        registry.(key).last_updated = now_str;
    else
        entry.first_seen   = now_str;
        entry.last_updated = now_str;
        entry.config       = canonical;
        registry.(key)     = entry;
    end

    json          = jsonencode(registry, 'PrettyPrint', true);
    registry_path = config_registry_path(patient_id, session, config);
    fid           = fopen(registry_path, 'w');
    if fid < 0
        warning('pipeline_simulate:registry', ...
            'Could not write config registry: %s', registry_path);
        return;
    end
    fprintf(fid, '%s', json);
    fclose(fid);
end

function registry = load_config_registry(patient_id, session, config)
    % Read <sim_dir>/config_registry.json or return an empty struct.
    registry      = struct();
    registry_path = config_registry_path(patient_id, session, config);
    if ~isfile(registry_path), return; end
    fid = fopen(registry_path, 'r');
    if fid < 0, return; end
    raw = fread(fid, '*char')';
    fclose(fid);
    if isempty(strtrim(raw)), return; end
    try
        decoded = jsondecode(raw);
        if isstruct(decoded)
            registry = decoded;
        end
    catch ME
        warning('pipeline_simulate:registry', ...
            'Failed to parse %s: %s', registry_path, ME.message);
    end
end

function canonical = decode_config_hash(hash_hex, patient_id, session, config)
    % Look up the canonical config that produced the given hash by reading
    % the per-session config_registry.json.  Returns [] when the hash is
    % unknown.
    canonical = [];
    registry  = load_config_registry(patient_id, session, config);
    key       = registry_key(hash_hex);
    if isfield(registry, key) && isfield(registry.(key), 'config')
        canonical = registry.(key).config;
    end
end


%% =========================================================================
%  PER-FIELD FILE NAMING (input/output paired by basename)
%% =========================================================================

function base = get_field_basename(field_dose, patient_id, session)
    % Derive the stem used for <stem>_recon_<hash>.mat and <stem>_sim_<hash>.mat
    % so reconstructions sit alongside their input dose by name.  Mirrors
    % step15_process_doses naming via field_dose.source_mat_filename.
    if isfield(field_dose, 'source_mat_filename') && ~isempty(field_dose.source_mat_filename)
        [~, base, ~] = fileparts(field_dose.source_mat_filename);
        return;
    end
    if isfield(field_dose, 'beam_index') && ~isempty(field_dose.beam_index)
        base = sprintf('dose_%s_%s_field_%03d', patient_id, session, ...
            field_dose.beam_index);
    else
        base = sprintf('dose_%s_%s_field', patient_id, session);
    end
end

function p = expected_recon_path(field_dose, patient_id, session, config, hash8)
    sim_dir = get_simulation_directory(patient_id, session, config);
    base    = get_field_basename(field_dose, patient_id, session);
    p       = fullfile(sim_dir, sprintf('%s_recon_%s.mat', base, hash8));
end

function p = expected_sim_path(field_dose, patient_id, session, config, hash8)
    sim_dose_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        patient_id, session, 'processed', 'simulated_doses');
    base = get_field_basename(field_dose, patient_id, session);
    p    = fullfile(sim_dose_dir, sprintf('%s_sim_%s.mat', base, hash8));
end

function p = expected_total_recon_path(patient_id, session, config, hash8)
    sim_dir = get_simulation_directory(patient_id, session, config);
    p       = fullfile(sim_dir, sprintf('total_recon_dose_%s.mat', hash8));
end

function p = expected_manifest_path(patient_id, session, config, hash8)
    sim_dir = get_simulation_directory(patient_id, session, config);
    p       = fullfile(sim_dir, sprintf('cache_manifest_%s.mat', hash8));
end

function status_dir = get_status_directory(patient_id, session, config)
    % Per-(session, config) directory holding the cross-instance lock dirs and
    % .status files. The caller is responsible for creating it once.
    status_dir = fullfile(get_simulation_directory(patient_id, session, config), ...
        'field_status');
end

function p = field_status_path(field_dose, patient_id, session, config, hash8)
    base = get_field_basename(field_dose, patient_id, session);
    p    = fullfile(get_status_directory(patient_id, session, config), ...
        sprintf('%s_%s.status', base, hash8));
end

function p = field_lock_path(field_dose, patient_id, session, config, hash8)
    base = get_field_basename(field_dose, patient_id, session);
    p    = fullfile(get_status_directory(patient_id, session, config), ...
        sprintf('%s_%s.lock', base, hash8));
end
