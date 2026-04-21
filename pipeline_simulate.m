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
CONFIG.enable_spherical_correction = false;    % Compute PSF filter via get_psf
CONFIG.regularization_lambda       = 0.01;   % Wiener regularization for PSF

% --- Sensor Placement ---
% Controls the k-Wave sensor mask geometry used in every per-field simulation.
%   'full_plane_anterior' : Full YZ plane at x = sensor_x_index.
%   'full_plane_lateral'  : Full XZ plane at y = sensor_y_index.
%   'spherical'           : Spherical shell (no PSF correction applied).
CONFIG.sensor_placement_method = 'full_plane_lateral';
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

% Initialize results structure
RESULTS           = struct();
RESULTS.timestamp = datetime('now');
RESULTS.config    = CONFIG;
RESULTS.patients  = struct();

%% ========================= MAIN PROCESSING LOOP ==========================

for p_idx = 1:length(CONFIG.patients)
    patient_id = CONFIG.patients{p_idx};

    for s_idx = 1:length(CONFIG.sessions)
        session = CONFIG.sessions{s_idx};

        fprintf('\n=== Processing: Patient %s, %s ===\n', patient_id, session);

        result_key = make_result_key(patient_id, session);
        RESULTS.patients.(result_key) = init_patient_result(patient_id, session);

        % Open per-session log file (appends across runs)
        log_fid = open_simulation_log(patient_id, session, CONFIG);
        log_msg(log_fid, 'Run started  patient %s, session %s', patient_id, session);

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
            fprintf('\n[STEP 1.5] Loading previously processed data...\n');
            [field_doses, sct_resampled, total_rs_dose, dose_metadata] = ...
                load_processed_data(patient_id, session, CONFIG);

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

                % Build acoustic medium from CT
                medium = create_acoustic_medium(sct_resampled, CONFIG);

                % Pre-compute PSF correction filter (once for all fields)
                if CONFIG.enable_spherical_correction
                    fprintf('         Computing PSF correction filter...\n');
                    psf_filter = get_psf(total_rs_dose, sct_resampled, medium, CONFIG);
                    fprintf('         PSF filter ready (%.1f s).\n', ...
                        psf_filter.computation_time_s);
                else
                    psf_filter = [];
                    fprintf('         PSF correction: disabled.\n');
                end

                valid_field_indices = find(~cellfun(@isempty, field_doses));
                num_fields          = length(valid_field_indices);
                grid_dims           = dose_metadata.dimensions;

                fprintf('         Processing %d field(s)...\n', num_fields);

                % --- Load cache manifest and cross-check against files on disk ---
                sim_dir_cache  = get_simulation_directory(patient_id, session, CONFIG);
                cache_manifest = load_cache_manifest(patient_id, session, CONFIG);
                completed_idxs = get_completed_from_manifest(cache_manifest);

                % Cross-check: also scan files in case manifest is missing entries
                existing_recon = dir(fullfile(sim_dir_cache, 'field_recon_*.mat'));
                file_idxs      = [];
                for ef = 1:numel(existing_recon)
                    tok = regexp(existing_recon(ef).name, 'field_recon_(\d+)\.mat', 'tokens');
                    if ~isempty(tok) && ~isempty(tok{1})
                        file_idxs(end+1) = str2double(tok{1}{1}); %#ok<AGROW>
                    end
                end
                completed_idxs = union(completed_idxs, file_idxs);

                pending_idxs = setdiff(valid_field_indices, completed_idxs);
                cached_idxs  = intersect(valid_field_indices, completed_idxs);
                num_pending  = numel(pending_idxs);

                print_cache_status(cache_manifest, cached_idxs, pending_idxs, field_doses, log_fid);

                total_recon = zeros(grid_dims);

                if CONFIG.use_parallel && num_pending > 1
                    fprintf('         Using parallel processing (parfor)...\n');
                    log_msg(log_fid, 'Parallel mode: %d fields pending', num_pending);
                    pool = gcp('nocreate');
                    if isempty(pool) || pool.NumWorkers ~= CONFIG.num_parallel_workers
                        if ~isempty(pool), delete(pool); end
                        parpool(CONFIG.num_parallel_workers);
                    end
                    recon_doses    = cell(num_pending, 1);
                    t_parfor       = tic;
                    parfor f = 1:num_pending
                        fi = pending_idxs(f);
                        fprintf('           Field %d (gantry: %.1f deg)...\n', ...
                            fi, field_doses{fi}.gantry_angle);
                        [recon_doses{f}, ~] = run_single_field_simulation(...
                            field_doses{fi}, sct_resampled, medium, ...
                            beam_metadata, CONFIG, psf_filter);
                        recon_doses{f} = gather(recon_doses{f});
                        save_field_reconstruction(recon_doses{f}, fi, ...
                            patient_id, session, CONFIG);
                    end
                    elapsed_parfor   = toc(t_parfor);
                    elapsed_per_field = elapsed_parfor / max(num_pending, 1);
                    for f = 1:num_pending
                        fi = pending_idxs(f);
                        total_recon    = total_recon + recon_doses{f};
                        cache_manifest = update_manifest_field(cache_manifest, fi, ...
                            field_doses{fi}.gantry_angle, elapsed_per_field);
                        log_msg(log_fid, 'Field %d complete (parallel, avg %.1f s/field)', ...
                            fi, elapsed_per_field);
                    end
                    save_cache_manifest(cache_manifest, patient_id, session, CONFIG);
                    fprintf('         Parallel block done in %.1f s (avg %.1f s/field).\n', ...
                        elapsed_parfor, elapsed_per_field);
                else
                    fprintf('         Using serial processing...\n');
                    log_msg(log_fid, 'Serial mode: %d fields pending', num_pending);
                    for f = 1:num_pending
                        fi       = pending_idxs(f);
                        t_field  = tic;
                        fprintf('           Field %d/%d (gantry: %.1f deg)...\n', ...
                            f, num_pending, field_doses{fi}.gantry_angle);
                        log_msg(log_fid, 'Field %d start (gantry %.1f deg)', ...
                            fi, field_doses{fi}.gantry_angle);
                        [recon_dose, ~] = run_single_field_simulation(...
                            field_doses{fi}, sct_resampled, medium, ...
                            beam_metadata, CONFIG, psf_filter);
                        total_recon = total_recon + recon_dose;
                        save_field_reconstruction(recon_dose, fi, patient_id, session, CONFIG);
                        elapsed_f  = toc(t_field);
                        cache_manifest = update_manifest_field(cache_manifest, fi, ...
                            field_doses{fi}.gantry_angle, elapsed_f);
                        save_cache_manifest(cache_manifest, patient_id, session, CONFIG);
                        log_msg(log_fid, 'Field %d complete (%.1f s)', fi, elapsed_f);
                        fprintf('             Done in %.1f s\n', elapsed_f);
                    end
                end

                % Accumulate cached field reconstructions
                if ~isempty(cached_idxs)
                    fprintf('         Loading %d cached reconstruction(s)...\n', numel(cached_idxs));
                    for c = 1:numel(cached_idxs)
                        cd_ = load(fullfile(sim_dir_cache, ...
                            sprintf('field_recon_%03d.mat', cached_idxs(c))));
                        total_recon = total_recon + cd_.recon_dose;
                    end
                end

                sim_time = toc(sim_start);
                save_total_reconstruction(total_recon, dose_metadata, patient_id, session, CONFIG);

                RESULTS.patients.(result_key).simulation_time_sec  = sim_time;
                RESULTS.patients.(result_key).total_recon_max_Gy   = max(total_recon(:));

                log_msg(log_fid, 'Step 2 complete  %.1f s total, max recon dose %.4f Gy', ...
                    sim_time, max(total_recon(:)));
                fprintf('[STEP 2] Complete (%.1f s). Max recon dose: %.4f Gy\n', ...
                    sim_time, max(total_recon(:)));
            else
                fprintf('\n[STEP 2] Loading previously computed reconstruction...\n');
                [total_recon, ~] = load_total_reconstruction(patient_id, session, CONFIG);
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

function save_field_reconstruction(recon_dose, field_idx, patient_id, session, config)
    sim_dir  = get_simulation_directory(patient_id, session, config);
    filename = sprintf('field_recon_%03d.mat', field_idx);
    save(fullfile(sim_dir, filename), 'recon_dose', '-v7.3');
end

function save_total_reconstruction(total_recon, metadata, patient_id, session, config)
    sim_dir = get_simulation_directory(patient_id, session, config);
    save(fullfile(sim_dir, 'total_recon_dose.mat'), 'total_recon', 'metadata', '-v7.3');
end

function [total_recon, metadata] = load_total_reconstruction(patient_id, session, config)
    sim_dir = get_simulation_directory(patient_id, session, config);
    data    = load(fullfile(sim_dir, 'total_recon_dose.mat'));
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

function log_fid = open_simulation_log(patient_id, session, config)
    % Opens (or creates) a per-session append-mode log file.
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
    fprintf(log_fid, '========================================\n');
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
%  CACHE MANIFEST HELPERS
%% =========================================================================

function manifest = load_cache_manifest(patient_id, session, config)
    % Loads the per-session cache manifest (or creates a fresh one).
    sim_dir       = get_simulation_directory(patient_id, session, config);
    manifest_path = fullfile(sim_dir, 'cache_manifest.mat');
    if exist(manifest_path, 'file')
        data     = load(manifest_path, 'manifest');
        manifest = data.manifest;
    else
        manifest              = struct();
        manifest.fields       = struct();   % sub-struct keyed by 'f_NNN'
        manifest.created      = datetime('now');
        manifest.last_updated = datetime('now');
    end
end

function save_cache_manifest(manifest, patient_id, session, config)
    sim_dir            = get_simulation_directory(patient_id, session, config);
    manifest.last_updated = datetime('now');
    save(fullfile(sim_dir, 'cache_manifest.mat'), 'manifest');
end

function manifest = update_manifest_field(manifest, field_idx, gantry_angle, elapsed_sec)
    % Records a completed field in the manifest.
    key                            = sprintf('f_%03d', field_idx);
    manifest.fields.(key).field_idx    = field_idx;
    manifest.fields.(key).gantry_angle = gantry_angle;
    manifest.fields.(key).completed_at = datetime('now');
    manifest.fields.(key).elapsed_sec  = elapsed_sec;
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

function print_cache_status(manifest, cached_idxs, pending_idxs, field_doses, log_fid)
    % Prints a formatted table of cached vs pending fields.
    fprintf('\n         --- Cache / Resume Status ---\n');
    if ~isempty(cached_idxs)
        fprintf('         Previously completed:\n');
        for i = 1:numel(cached_idxs)
            fi  = cached_idxs(i);
            key = sprintf('f_%03d', fi);
            if isfield(manifest.fields, key)
                e = manifest.fields.(key);
                fprintf('           [DONE] Field %3d  gantry %6.1f deg  completed %s  (%.0f s)\n', ...
                    fi, e.gantry_angle, ...
                    datestr(e.completed_at, 'yyyy-mm-dd HH:MM:SS'), e.elapsed_sec);
                log_msg(log_fid, 'CACHE HIT  field %d (gantry %.1f deg), completed %s', ...
                    fi, e.gantry_angle, datestr(e.completed_at, 'yyyy-mm-dd HH:MM:SS'));
            else
                fprintf('           [DONE] Field %3d  (file on disk, no manifest entry)\n', fi);
                log_msg(log_fid, 'CACHE HIT  field %d (file present, no manifest entry)', fi);
            end
        end
    end
    if ~isempty(pending_idxs)
        fprintf('         Pending:\n');
        for i = 1:numel(pending_idxs)
            fi = pending_idxs(i);
            if ~isempty(field_doses) && numel(field_doses) >= fi && ~isempty(field_doses{fi})
                fprintf('           [TODO] Field %3d  gantry %6.1f deg\n', ...
                    fi, field_doses{fi}.gantry_angle);
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
