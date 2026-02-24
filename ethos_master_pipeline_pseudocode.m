%% =========================================================================
%  ETHOS_MASTER_PIPELINE.m
%  Photoacoustic Dose Reconstruction Pipeline for Varian ETHOS
%  =========================================================================
%
%  USAGE:
%    1. Configure parameters in the CONFIGURATION section
%    2. Run this script
%    3. After Step 0.5 completes, manually export field doses from Raystation
%    4. Re-run script to continue from Step 1.5
%
%  PREREQUISITES:
%    - MATLAB R2022a or later
%    - k-Wave Toolbox (http://www.k-wave.org)
%    - Image Processing Toolbox
%    - Parallel Computing Toolbox (optional)
%
%  AUTHOR: [Your Name]
%  DATE: February 2026
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION =================================

% --- Patient and Session Selection ---
CONFIG.patients = {'1194203'};
CONFIG.sessions = {'Session_1'};
CONFIG.treatment_site = 'Pancreas';

% --- Directory Paths ---
CONFIG.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.matrad_path = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

% --- Gruneisen Parameter Method ---
% Options: 'uniform', 'threshold_1', 'threshold_2'
CONFIG.gruneisen_method = 'threshold_2';

% --- Acoustic Simulation Parameters ---
CONFIG.dose_per_pulse_cGy = 0.16;       % cGy per LINAC pulse
CONFIG.pml_size = 10;                    % Perfectly matched layer (voxels)
CONFIG.cfl_number = 0.3;                 % CFL stability criterion
CONFIG.use_gpu = true;                   % GPU acceleration
CONFIG.num_time_reversal_iter = 1;       % Time reversal iterations
CONFIG.enable_spherical_correction = true;   % Compute PSF correction via get_psf
CONFIG.regularization_lambda = 0.01;         % Wiener regularization for PSF filter

% --- Sensor Placement Parameters ---
CONFIG.sensor_size_cm       = [10, 10];  % Physical sensor [X, Z] in cm
CONFIG.sensor_standoff_mm   = 5;         % Gap between body surface and sensor (mm)
CONFIG.element_size_mm      = [];        % Element patch size for averaging (mm)
                                         % Empty = no averaging; set to sweep parametrically
CONFIG.jaw_margin_mm        = 10;        % Extra margin around jaw projection (mm)
CONFIG.sensor_placement     = 'anterior'; % Placement side
CONFIG.sensor_mode          = 'full_lateral_plane';  % Sensor geometry for dose visualization
                                                      %   'full_lateral_plane': sensor.mask(1,:,:) = 1
                                                      %   (entire first lateral face of the grid)

% --- Gamma Analysis Parameters ---
CONFIG.gamma_dose_pct = 3.0;             % Dose difference threshold (%)
CONFIG.gamma_dist_mm = 3.0;              % Distance-to-agreement (mm)
CONFIG.gamma_dose_cutoff_pct = 10.0;     % Ignore voxels below this % of max

% --- SSIM & Visualization Parameters ---
CONFIG.analysis_compute_ssim = true;     % Compute SSIM alongside gamma
CONFIG.analysis_plot_results = false;    % Generate PNG figures
CONFIG.analysis_plot_slices  = 'auto';   % Slice indices or 'auto' (25/50/75th %ile)

% --- MLC Gap Correction Parameters ---
CONFIG.mlc_min_gap_mm = 0.5;             % Minimum allowed gap
CONFIG.mlc_expansion_mm = 0.4;           % Expansion per side
CONFIG.mlc_position_range = [-140, 140]; % Valid leaf position range (mm)

% --- Pipeline Control Flags ---
CONFIG.run_step0  = false;   % Sort DICOM files
CONFIG.run_step05 = false;   % Fix MLC gaps
CONFIG.run_step15 = false;   % Process doses and resample CT
CONFIG.run_step2  = false;   % k-Wave simulation
CONFIG.run_step3  = false;   % Gamma analysis
CONFIG.use_parallel = false; % Use parfor for simulations

% --- Define Tissue Property Tables ---
CONFIG.tissue_tables = define_tissue_tables();

%% ========================= INITIALIZATION ================================

fprintf('=========================================================\n');
fprintf('  ETHOS Photoacoustic Dose Reconstruction Pipeline\n');
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
    error('k-Wave toolbox not found. Please add k-Wave to MATLAB path.');
end

% Initialize results structure
RESULTS = struct();
RESULTS.timestamp = datetime('now');
RESULTS.config = CONFIG;
RESULTS.patients = struct();

%% ========================= MAIN PROCESSING LOOP ==========================

for p_idx = 1:length(CONFIG.patients)
    patient_id = CONFIG.patients{p_idx};
    
    for s_idx = 1:length(CONFIG.sessions)
        session = CONFIG.sessions{s_idx};
        
        fprintf('\n=== Processing: Patient %s, %s ===\n', patient_id, session);
        
        % Create unique key for results storage
        result_key = make_result_key(patient_id, session);
        RESULTS.patients.(result_key) = init_patient_result(patient_id, session);
        
        try
            % ============================================================
            % STEP 0: Sort DICOM Files
            % ============================================================
            if CONFIG.run_step0
                fprintf('\n[STEP 0] Sorting DICOM files...\n');
                
                sct_dir = step0_sort_dicom(patient_id, session, CONFIG);
                
                RESULTS.patients.(result_key).sct_dir = sct_dir;
                fprintf('[STEP 0] Complete. Output: %s\n', sct_dir);
            else
                sct_dir = get_sct_directory(patient_id, session, CONFIG);
                fprintf('[STEP 0] Skipped. Using existing: %s\n', sct_dir);
            end
            
            % ============================================================
            % STEP 0.5: Fix MLC Gaps in RTPLAN
            % ============================================================
            if CONFIG.run_step05
                fprintf('\n[STEP 0.5] Fixing MLC gaps...\n');
                
                [adjusted_rtplan_path, num_corrections] = ...
                    step05_fix_mlc_gaps(patient_id, session, CONFIG);
                
                RESULTS.patients.(result_key).adjusted_rtplan = adjusted_rtplan_path;
                RESULTS.patients.(result_key).mlc_corrections = num_corrections;
                fprintf('[STEP 0.5] Complete. Corrections: %d\n', num_corrections);
            end
            
            % ============================================================
            % STEP 1: Manual Raystation Export (External)
            % ============================================================
            fprintf('\n[STEP 1] Checking for Raystation files...\n');
            
            rs_dir = get_raystation_directory(patient_id, session, CONFIG);
            [rs_files_exist, num_dose_files] = check_raystation_files(rs_dir);
            
            if ~rs_files_exist
                fprintf('[STEP 1] WAITING: No Raystation files found.\n');
                fprintf('         Expected location: %s\n', rs_dir);
                fprintf('         Please:\n');
                fprintf('           1. Import adjusted RTPLAN into Raystation\n');
                fprintf('           2. Recalculate dose for each field\n');
                fprintf('           3. Export field doses as RD*.dcm\n');
                fprintf('           4. Re-run this pipeline\n');
                
                RESULTS.patients.(result_key).status = 'awaiting_raystation';
                continue;  % Skip to next patient/session
            end
            
            fprintf('[STEP 1] Found %d dose files in %s\n', num_dose_files, rs_dir);
            
            % ============================================================
            % STEP 1.5: Process Field Doses and Resample CT
            % ============================================================
            if CONFIG.run_step15
                fprintf('\n[STEP 1.5] Processing field doses and resampling CT...\n');
                
                % Process doses - saves individual field_dose_XXX.mat files (memory constraint)
                % Returns cell array loaded from individual files, plus sct and total dose
                [field_doses, sct_resampled, total_rs_dose, dose_metadata] = ...
                    step15_process_doses(patient_id, session, CONFIG);
                
                num_valid_fields = sum(~cellfun(@isempty, field_doses));
                RESULTS.patients.(result_key).num_fields = num_valid_fields;
                RESULTS.patients.(result_key).dose_grid_size = dose_metadata.dimensions;
                RESULTS.patients.(result_key).total_rs_dose_max_Gy = max(total_rs_dose(:));
                
                fprintf('[STEP 1.5] Complete. Processed %d fields (saved as individual files).\n', num_valid_fields);
                fprintf('           Dose grid: [%d x %d x %d]\n', dose_metadata.dimensions);
                fprintf('           Max RS dose: %.4f Gy\n', max(total_rs_dose(:)));
            else
                fprintf('\n[STEP 1.5] Loading previously processed data...\n');
                % Load from individual files
                [field_doses, sct_resampled, total_rs_dose, dose_metadata] = ...
                    load_processed_data(patient_id, session, CONFIG);
            end
            
            % Extract beam metadata (isocenter + jaw positions) for sensor placement.
            % beam_metadata is stored in dose_metadata by step15_process_doses.
            if isfield(dose_metadata, 'beam_metadata') && ~isempty(dose_metadata.beam_metadata)
                beam_metadata = dose_metadata.beam_metadata;
                fprintf('         Beam metadata: %d beams (with isocenter/jaw data for sensor placement)\n', ...
                    length(beam_metadata));
            else
                beam_metadata = [];
                fprintf('         [WARNING] No beam metadata in dose_metadata — sensor placement will use legacy mode\n');
            end
            
            % ============================================================
            % STEP 2: k-Wave Photoacoustic Simulation
            % ============================================================
            if CONFIG.run_step2
                fprintf('\n[STEP 2] Running k-Wave simulations...\n');
                
                sim_start = tic;
                
                % Create acoustic medium from CT
                medium = create_acoustic_medium(sct_resampled, CONFIG);

                % Pre-compute PSF correction filter (once for all fields)
                if CONFIG.enable_spherical_correction
                    fprintf('         Computing PSF correction filter...\n');
                    psf_filter = get_psf(total_rs_dose, sct_resampled, medium, CONFIG);
                    fprintf('         PSF filter computed (%.1f s)\n', psf_filter.computation_time_s);
                else
                    psf_filter = [];
                    fprintf('         PSF correction: disabled\n');
                end

                % Get number of valid fields
                valid_field_indices = find(~cellfun(@isempty, field_doses));
                num_fields = length(valid_field_indices);
                
                fprintf('         Processing %d fields...\n', num_fields);
                
                % Initialize total reconstructed dose
                grid_dims = dose_metadata.dimensions;
                
                if CONFIG.use_parallel && num_fields > 1
                    % Parallel processing
                    fprintf('         Using parallel processing (parfor)...\n');
                    
                    recon_doses = cell(num_fields, 1);
                    
                    parfor f_idx = 1:num_fields
                        field_idx = valid_field_indices(f_idx);
                        
                        fprintf('           Field %d (gantry: %.1f°)...\n', ...
                            field_idx, field_doses{field_idx}.gantry_angle);
                        
                        [recon_doses{f_idx}, ~] = run_single_field_simulation(...
                            field_doses{field_idx}, sct_resampled, medium, ...
                            beam_metadata, CONFIG, psf_filter);
                        
                        % Save individual field result
                        save_field_reconstruction(recon_doses{f_idx}, field_idx, ...
                            patient_id, session, CONFIG);
                    end
                    
                    % Sum all reconstructed doses
                    total_recon = zeros(grid_dims);
                    for f_idx = 1:num_fields
                        total_recon = total_recon + recon_doses{f_idx};
                    end
                    
                else
                    % Serial processing
                    fprintf('         Using serial processing...\n');
                    
                    total_recon = zeros(grid_dims);
                    sim_results_all = cell(num_fields, 1);
                    
                    for f_idx = 1:num_fields
                        field_idx = valid_field_indices(f_idx);
                        
                        fprintf('           Field %d/%d (gantry: %.1f°)...\n', ...
                            f_idx, num_fields, field_doses{field_idx}.gantry_angle);
                        
                        [recon_dose, sim_results] = run_single_field_simulation(...
                            field_doses{field_idx}, sct_resampled, medium, ...
                            beam_metadata, CONFIG, psf_filter);
                        
                        total_recon = total_recon + recon_dose;
                        sim_results_all{f_idx} = sim_results;
                        
                        % Save individual field result (with sensor geometry)
                        save_field_reconstruction(recon_dose, field_idx, ...
                            patient_id, session, CONFIG);
                    end
                end
                
                sim_time = toc(sim_start);
                
                % Save total reconstructed dose
                save_total_reconstruction(total_recon, dose_metadata, ...
                    patient_id, session, CONFIG);
                
                RESULTS.patients.(result_key).simulation_time_sec = sim_time;
                RESULTS.patients.(result_key).total_recon_max_Gy = max(total_recon(:));
                
                fprintf('[STEP 2] Complete. Time: %.1f sec\n', sim_time);
                fprintf('         Max reconstructed dose: %.4f Gy\n', max(total_recon(:)));
            else
                fprintf('\n[STEP 2] Loading previously computed reconstruction...\n');
                [total_recon, ~] = load_total_reconstruction(patient_id, session, CONFIG);
            end
            
            % ============================================================
            % STEP 3: Gamma Analysis, SSIM & Visualization
            % ============================================================
            if CONFIG.run_step3
                total_recon = double(gather(total_recon));
                fprintf('\n[STEP 3] Running analysis (gamma + SSIM)...\n');
                
                % Call step3_analysis â€” loads all data, computes gamma and
                % SSIM, optionally generates figures, saves results to disk.
                % Signature matches: step3_analysis(patient_id, session, config)
                step3_results = step3_analysis(patient_id, session, CONFIG);
                
                % Store gamma pass rates in pipeline RESULTS
                RESULTS.patients.(result_key).gamma_ethos_vs_rs_pass_pct = ...
                    step3_results.ethos_vs_rs.gamma.pass_rate;
                RESULTS.patients.(result_key).gamma_rs_vs_recon_pass_pct = ...
                    step3_results.rs_vs_recon.gamma.pass_rate;
                
                % Store SSIM if computed
                if CONFIG.analysis_compute_ssim && ~isempty(step3_results.ethos_vs_rs.ssim)
                    RESULTS.patients.(result_key).ssim_ethos_vs_rs = ...
                        step3_results.ethos_vs_rs.ssim.ssim_3d;
                    RESULTS.patients.(result_key).ssim_rs_vs_recon = ...
                        step3_results.rs_vs_recon.ssim.ssim_3d;
                end
                
                fprintf('[STEP 3] Complete.\n');
                fprintf('         ETHOS vs RS gamma pass rate: %.1f%%\n', ...
                    step3_results.ethos_vs_rs.gamma.pass_rate);
                fprintf('         RS vs Recon gamma pass rate: %.1f%%\n', ...
                    step3_results.rs_vs_recon.gamma.pass_rate);
                if CONFIG.analysis_compute_ssim && ~isempty(step3_results.ethos_vs_rs.ssim)
                    fprintf('         ETHOS vs RS SSIM: %.4f\n', ...
                        step3_results.ethos_vs_rs.ssim.ssim_3d);
                    fprintf('         RS vs Recon SSIM: %.4f\n', ...
                        step3_results.rs_vs_recon.ssim.ssim_3d);
                end
            end
            
            % Mark patient/session as complete
            RESULTS.patients.(result_key).status = 'complete';
            fprintf('\n=== %s/%s: COMPLETE ===\n', patient_id, session);
            
        catch ME
            % Handle errors gracefully
            fprintf('\n[ERROR] %s\n', ME.message);
            for k = 1:length(ME.stack)
                fprintf('        In %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            
            RESULTS.patients.(result_key).status = 'error';
            RESULTS.patients.(result_key).error.message = ME.message;
            RESULTS.patients.(result_key).error.stack = ME.stack;
            
            fprintf('\n=== %s/%s: FAILED ===\n', patient_id, session);
            continue;  % Continue with next patient/session
        end
        
    end  % session loop
end  % patient loop

%% ========================= FINALIZE ======================================

fprintf('\n=========================================================\n');
fprintf('  Pipeline Complete\n');
fprintf('=========================================================\n');

% Save results
results_filename = sprintf('pipeline_results_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
results_path = fullfile(CONFIG.working_dir, 'AnalysisResults', results_filename);

if ~exist(fullfile(CONFIG.working_dir, 'AnalysisResults'), 'dir')
    mkdir(fullfile(CONFIG.working_dir, 'AnalysisResults'));
end

save(results_path, 'RESULTS', '-v7.3');
fprintf('  Results saved: %s\n', results_path);

% Generate summary
generate_pipeline_summary(RESULTS);

fprintf('\n=========================================================\n\n');


%% ========================= HELPER FUNCTIONS ==============================

function tables = define_tissue_tables()
    % Define tissue property lookup tables for different thresholding methods
    
    % THRESHOLD_1: Detailed tissue segmentation
    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air', 'Lung', 'Fat', 'Water', 'Blood', 'Muscle', 'SoftTissue', 'Bone', 'Metal'};
    tables.threshold_1.density       = [1.2,   400,   920,   1000,  1060,   1050,    1040,        1900,  7800];
    tables.threshold_1.sound_speed   = [343,   600,   1450,  1480,  1575,   1580,    1540,        3200,  5900];
    tables.threshold_1.alpha_coeff   = [0,     0.5,   0.48,  0.0022, 0.2,   0.5,     0.5,         10,    0];
    tables.threshold_1.alpha_power   = [1.0,   1.5,   1.5,   2.0,   1.3,    1.0,     1.1,         1.0,   1.0];
    tables.threshold_1.gruneisen     = [0,     0.5,   0.7,   0.11,  0.15,   0.2,     1.0,         0,     0];
    
    % THRESHOLD_2: Simplified 4-tissue model
    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water', 'Fat', 'SoftTissue', 'Bone'};
    tables.threshold_2.density       = [1000,    920,  1040,          1900];
    tables.threshold_2.sound_speed   = [1480,    1450, 1540,          3200];
    tables.threshold_2.alpha_coeff   = [0.0022,  0.48, 0.5,           10];
    tables.threshold_2.alpha_power   = [2.0,     1.5,  1.1,           1.0];
    tables.threshold_2.gruneisen     = [0.11,    0.7,  1.0,           0];
    
    % UNIFORM: Single property values
    tables.uniform = struct();
    tables.uniform.density      = 1000;
    tables.uniform.sound_speed  = 1540;
    tables.uniform.alpha_coeff  = 0.5;
    tables.uniform.alpha_power  = 1.1;
    tables.uniform.gruneisen    = 1.0;
end

function key = make_result_key(patient_id, session)
    % Create valid MATLAB struct field name from patient_id and session
    key = sprintf('P%s_%s', patient_id, strrep(session, '_', ''));
end

function result = init_patient_result(patient_id, session)
    % Initialize result structure for a patient/session
    result = struct();
    result.patient_id = patient_id;
    result.session = session;
    result.status = 'started';
    result.start_time = datetime('now');
end

function sct_dir = get_sct_directory(patient_id, session, config)
    % Return path to SCT directory
    sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
        config.treatment_site, session, 'sct');
end

function rs_dir = get_raystation_directory(patient_id, session, config)
    % Return path to Raystation output directory
    rs_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session);
end

function sim_dir = get_simulation_directory(patient_id, session, config)
    % Return path to simulation output directory
    sim_dir = fullfile(config.working_dir, 'SimulationResults', patient_id, ...
        session, config.gruneisen_method);
    if ~exist(sim_dir, 'dir')
        mkdir(sim_dir);
    end
end

function analysis_dir = get_analysis_directory(patient_id, session, config)
    % Return path to analysis output directory
    analysis_dir = fullfile(config.working_dir, 'AnalysisResults', patient_id, session);
    if ~exist(analysis_dir, 'dir')
        mkdir(analysis_dir);
    end
end

function [exists, num_files] = check_raystation_files(rs_dir)
    % Check if Raystation dose files exist
    exists = false; 
    num_files = 0; 
    if  exist(rs_dir, 'dir') > 0
        exists = true;
    
    
    rd_files = dir(fullfile(rs_dir, 'RD.*.dcm'));
    num_files = length(rd_files)
    end
end

function generate_pipeline_summary(results)
    % Print summary of pipeline results
    
    fprintf('\n--- Pipeline Summary ---\n');
    
    patient_fields = fieldnames(results.patients);
    
    for i = 1:length(patient_fields)
        p = results.patients.(patient_fields{i});
        fprintf('\n%s / %s:\n', p.patient_id, p.session);
        fprintf('  Status: %s\n', p.status);
        
        if strcmp(p.status, 'complete')
            if isfield(p, 'num_fields')
                fprintf('  Fields processed: %d\n', p.num_fields);
            end
            if isfield(p, 'simulation_time_sec')
                fprintf('  Simulation time: %.1f sec\n', p.simulation_time_sec);
            end
            if isfield(p, 'gamma_ethos_vs_rs_pass_pct')
                fprintf('  Gamma (ETHOS vs RS): %.1f%%\n', p.gamma_ethos_vs_rs_pass_pct);
            end
            if isfield(p, 'gamma_rs_vs_recon_pass_pct')
                fprintf('  Gamma (RS vs Recon): %.1f%%\n', p.gamma_rs_vs_recon_pass_pct);
            end
            if isfield(p, 'ssim_ethos_vs_rs')
                fprintf('  SSIM  (ETHOS vs RS): %.4f\n', p.ssim_ethos_vs_rs);
            end
            if isfield(p, 'ssim_rs_vs_recon')
                fprintf('  SSIM  (RS vs Recon): %.4f\n', p.ssim_rs_vs_recon);
            end
        elseif strcmp(p.status, 'error')
            fprintf('  Error: %s\n', p.error.message);
        end
    end
end


%% ========================= STEP 2 HELPER FUNCTIONS =======================
% These functions support Step 2 (k-Wave simulation) and are called
% directly from the pipeline script. Step 3 functions are implemented
% in step3_analysis.m and called via step3_analysis(patient_id, session, CONFIG).

function save_field_reconstruction(recon_dose, field_idx, patient_id, session, config)
    % Save individual field reconstruction
    sim_dir = get_simulation_directory(patient_id, session, config);
    filename = sprintf('field_recon_%03d.mat', field_idx);
    save(fullfile(sim_dir, filename), 'recon_dose', '-v7.3');
end

function save_total_reconstruction(total_recon, metadata, patient_id, session, config)
    % Save total reconstructed dose
    sim_dir = get_simulation_directory(patient_id, session, config);
    save(fullfile(sim_dir, 'total_recon_dose.mat'), 'total_recon', 'metadata', '-v7.3');
end

function [total_recon, metadata] = load_total_reconstruction(patient_id, session, config)
    % Load previously computed total reconstruction
    sim_dir = get_simulation_directory(patient_id, session, config);
    data = load(fullfile(sim_dir, 'total_recon_dose.mat'));
    total_recon = data.total_recon;
    metadata = data.metadata;
end