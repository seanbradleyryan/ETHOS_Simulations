%% =========================================================================
%  TEST_STEP2_KWAVE_SIMULATION.m
%  Test script for Step 2 k-Wave simulation functions
%  =========================================================================
%
%  This script demonstrates how to use the Step 2 k-Wave simulation
%  functions either standalone or as part of the master pipeline.
%
%  USAGE:
%    1. Set patient_id and session below
%    2. Run this script
%
%  PREREQUISITES:
%    - Step 1.5 must have been completed (processed data exists)
%    - k-Wave toolbox must be on MATLAB path
%
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION =================================

% Patient and session to process
patient_id = '1194203';
session = 'Session_1';

% Configuration struct
CONFIG = struct();
CONFIG.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.treatment_site = 'Pancreas';

% Gruneisen method: 'uniform', 'threshold_1', 'threshold_2'
CONFIG.gruneisen_method = 'threshold_2';

% Simulation parameters
CONFIG.dose_per_pulse_cGy = 0.16;   % cGy per LINAC pulse
CONFIG.pml_size = 10;                % PML thickness (voxels)
CONFIG.cfl_number = 0.3;             % CFL stability criterion
CONFIG.use_gpu = true;               % GPU acceleration
CONFIG.use_parallel = false;         % Use parfor for multiple fields
CONFIG.plot_sim = false;             % Show k-Wave plots during simulation

% Define tissue tables (or use defaults)
CONFIG.tissue_tables = define_tissue_tables();

%% ========================= VERIFY PREREQUISITES ==========================

fprintf('=========================================================\n');
fprintf('  Step 2 Test Script - k-Wave Simulation\n');
fprintf('=========================================================\n\n');

% Check for k-Wave
if ~exist('kWaveGrid', 'file')
    error('k-Wave toolbox not found. Please add k-Wave to MATLAB path.');
end
fprintf('[OK] k-Wave toolbox found\n');

% Check for processed data
processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
    patient_id, session, 'processed');

if ~exist(processed_dir, 'dir')
    error('Processed directory not found: %s\nRun step15_process_doses first.', processed_dir);
end
fprintf('[OK] Processed directory exists: %s\n', processed_dir);

% Check for field dose files
field_files = dir(fullfile(processed_dir, 'field_dose_*.mat'));
fprintf('[OK] Found %d field dose files\n', length(field_files));

% Check for SCT resampled
if ~exist(fullfile(processed_dir, 'sct_resampled.mat'), 'file')
    error('SCT resampled file not found');
end
fprintf('[OK] SCT resampled file exists\n');

% Check GPU availability
if CONFIG.use_gpu
    try
        g = gpuDevice;
        fprintf('[OK] GPU available: %s (%.1f GB)\n', g.Name, g.TotalMemory/1e9);
    catch
        fprintf('[WARN] GPU not available, will use CPU\n');
        CONFIG.use_gpu = false;
    end
end

fprintf('\n');

%% ========================= RUN SIMULATION ================================

fprintf('Starting Step 2 simulation...\n\n');

% Run the full simulation
[total_recon, field_recons, sim_log] = step2_run_all_simulations(patient_id, session, CONFIG);

%% ========================= DISPLAY RESULTS ===============================

fprintf('\n=========================================================\n');
fprintf('  Simulation Results Summary\n');
fprintf('=========================================================\n');
fprintf('  Total time: %.1f seconds\n', sim_log.total_time_sec);
fprintf('  Fields processed: %d\n', sim_log.num_fields_processed);
fprintf('  Max reconstructed dose: %.4f Gy\n', sim_log.total_recon_max_Gy);
fprintf('=========================================================\n\n');

%% ========================= VISUALIZATION =================================

% Load original doses for comparison
fprintf('Loading original doses for comparison...\n');

% Load total Raystation dose
total_rs_file = fullfile(processed_dir, 'total_rs_dose.mat');
if exist(total_rs_file, 'file')
    rs_data = load(total_rs_file);
    if isfield(rs_data, 'total_rs_dose')
        total_rs_dose = rs_data.total_rs_dose;
    elseif isfield(rs_data, 'total_dose')
        total_rs_dose = rs_data.total_dose;
    else
        total_rs_dose = [];
    end
else
    total_rs_dose = [];
end

% Find slice with maximum dose
[~, max_idx] = max(total_recon(:));
[max_x, max_y, max_z] = ind2sub(size(total_recon), max_idx);

%% Create comparison figure
figure('Name', 'Step 2 Results', 'Position', [100, 100, 1400, 900], 'Color', 'w');

% Axial slice through max dose
subplot(2, 3, 1);
imagesc(squeeze(total_recon(:, :, max_z)));
colorbar;
title(sprintf('Reconstructed - Axial (z=%d)', max_z));
xlabel('Y'); ylabel('X');
axis image;

if ~isempty(total_rs_dose)
    subplot(2, 3, 2);
    imagesc(squeeze(total_rs_dose(:, :, max_z)));
    colorbar;
    title(sprintf('Original (RS) - Axial (z=%d)', max_z));
    xlabel('Y'); ylabel('X');
    axis image;
    
    subplot(2, 3, 3);
    diff_dose = total_recon - total_rs_dose;
    imagesc(squeeze(diff_dose(:, :, max_z)));
    colorbar;
    title('Difference (Recon - Orig)');
    xlabel('Y'); ylabel('X');
    axis image;
end

% Coronal slice through max dose
subplot(2, 3, 4);
imagesc(squeeze(total_recon(:, max_y, :)));
colorbar;
title(sprintf('Reconstructed - Coronal (y=%d)', max_y));
xlabel('Z'); ylabel('X');
axis image;

if ~isempty(total_rs_dose)
    subplot(2, 3, 5);
    imagesc(squeeze(total_rs_dose(:, max_y, :)));
    colorbar;
    title(sprintf('Original (RS) - Coronal (y=%d)', max_y));
    xlabel('Z'); ylabel('X');
    axis image;
end

% Dose profile comparison
subplot(2, 3, 6);
recon_profile = squeeze(total_recon(max_x, :, max_z));
plot(recon_profile, 'b-', 'LineWidth', 2, 'DisplayName', 'Reconstructed');
hold on;
if ~isempty(total_rs_dose)
    orig_profile = squeeze(total_rs_dose(max_x, :, max_z));
    plot(orig_profile, 'r--', 'LineWidth', 2, 'DisplayName', 'Original');
end
xlabel('Y Index');
ylabel('Dose (Gy)');
title('Dose Profile Comparison');
legend('Location', 'best');
grid on;

sgtitle(sprintf('Patient: %s, Session: %s, Method: %s', ...
    patient_id, session, CONFIG.gruneisen_method), 'FontSize', 14);

%% Save figure
sim_dir = fullfile(CONFIG.working_dir, 'SimulationResults', ...
    patient_id, session, CONFIG.gruneisen_method);
saveas(gcf, fullfile(sim_dir, 'reconstruction_comparison.png'));
saveas(gcf, fullfile(sim_dir, 'reconstruction_comparison.fig'));
fprintf('Visualization saved to: %s\n', sim_dir);

%% ========================= METRICS =======================================

if ~isempty(total_rs_dose)
    fprintf('\n--- Reconstruction Metrics ---\n');
    
    % Only evaluate where dose is significant
    dose_threshold = 0.01 * max(total_rs_dose(:));
    valid_mask = total_rs_dose > dose_threshold;
    
    if any(valid_mask(:))
        orig_valid = total_rs_dose(valid_mask);
        recon_valid = total_recon(valid_mask);
        diff_valid = recon_valid - orig_valid;
        
        % Mean Absolute Error
        mae = mean(abs(diff_valid));
        fprintf('  Mean Absolute Error: %.4f Gy\n', mae);
        
        % Root Mean Square Error
        rmse = sqrt(mean(diff_valid.^2));
        fprintf('  Root Mean Square Error: %.4f Gy\n', rmse);
        
        % Mean Relative Error
        mre = mean(abs(diff_valid) ./ orig_valid) * 100;
        fprintf('  Mean Relative Error: %.2f%%\n', mre);
        
        % Correlation coefficient
        r = corrcoef(orig_valid, recon_valid);
        fprintf('  Correlation Coefficient: %.4f\n', r(1, 2));
        
        % Max dose comparison
        fprintf('  Original max dose: %.4f Gy\n', max(total_rs_dose(:)));
        fprintf('  Reconstructed max dose: %.4f Gy\n', max(total_recon(:)));
        fprintf('  Max dose ratio: %.2f%%\n', 100 * max(total_recon(:)) / max(total_rs_dose(:)));
    end
end

fprintf('\n=========================================================\n');
fprintf('  Test Complete\n');
fprintf('=========================================================\n');


%% ========================= LOCAL FUNCTION ================================

function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Create tissue property lookup tables
    tables = struct();
    
    % UNIFORM
    tables.uniform = struct();
    tables.uniform.density = 1000;
    tables.uniform.sound_speed = 1540;
    tables.uniform.alpha_coeff = 0.5;
    tables.uniform.alpha_power = 1.1;
    tables.uniform.gruneisen = 1.0;
    
    % THRESHOLD_1: Detailed (9 tissues)
    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1024, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air', 'Lung', 'Fat', 'Water', 'Blood', 'Muscle', 'SoftTissue', 'Bone', 'Metal'};
    tables.threshold_1.density       = [1.2,   400,   920,   1000,  1060,   1050,    1040,        1900,  7800];
    tables.threshold_1.sound_speed   = [343,   600,   1450,  1480,  1575,   1580,    1540,        3200,  5900];
    tables.threshold_1.alpha_coeff   = [0,     0.5,   0.48,  0.0022, 0.2,   0.5,     0.5,         10,    0];
    tables.threshold_1.alpha_power   = [1.0,   1.5,   1.5,   2.0,   1.3,    1.0,     1.1,         1.0,   1.0];
    tables.threshold_1.gruneisen     = [0,     0.5,   0.7,   0.11,  0.15,   0.2,     1.0,         0,     0];
    
    % THRESHOLD_2: Simplified (4 tissues)
    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1024, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water', 'Fat', 'SoftTissue', 'Bone'};
    tables.threshold_2.density       = [1000,    920,  1040,          1900];
    tables.threshold_2.sound_speed   = [1480,    1450, 1540,          3200];
    tables.threshold_2.alpha_coeff   = [0.0022,  0.48, 0.5,           10];
    tables.threshold_2.alpha_power   = [2.0,     1.5,  1.1,           1.0];
    tables.threshold_2.gruneisen     = [0.11,    0.7,  1.0,           0];
end
