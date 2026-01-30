%% visualize_rtdose_overlay.m
% Script to visualize exported RT dose files
% For analyzing individual beam/segment doses from RayStation exports
%
% This script:
% 1. Loads selected RT dose files
% 2. Finds the slice with maximum dose for each file
% 3. Creates visualizations showing only the dose distribution
%
% Author: Generated for ETHOS Simulations project
% Date: January 2025

clear; clc; close all;

%% Configuration
% Working directory base path
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';

% Patient and session
patient_id = '1194203';
session = 'Session_1';

% Path to exported dose files
dose_folder = fullfile(wd, 'RayStationFiles', patient_id, session);

% Number of dose files to randomly select and plot
num_doses_to_plot = 3;  % Pick a few dose files to visualize

% Visualization parameters
dose_colormap = 'jet';           % Colormap for dose display

%% Find and Load RT Dose Files
dose_data = {};
dose_labels = {};

fprintf('\nSearching for dose files in: %s\n', dose_folder);

% Find all RD*.dcm files
all_dose_files = dir(fullfile(dose_folder, 'RD*.dcm'));

if isempty(all_dose_files)
    error('No RD*.dcm files found in: %s', dose_folder);
end

fprintf('Found %d total dose files\n', length(all_dose_files));

% Randomly select a few to plot
num_to_select = min(num_doses_to_plot, length(all_dose_files));
if length(all_dose_files) > num_to_select
    % Random selection
    selected_indices = randperm(length(all_dose_files), num_to_select);
else
    % Use all available
    selected_indices = 1:length(all_dose_files);
end

fprintf('Randomly selected %d dose files to visualize\n\n', num_to_select);

% Load selected dose files
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    dose_file = all_dose_files(idx);
    dose_filepath = fullfile(dose_folder, dose_file.name);
    
    fprintf('Loading [%d/%d]: %s\n', i, num_to_select, dose_file.name);
    
    % Load RT dose
    [dose_volume, dose_info] = load_rtdose(dose_filepath);
    
    % Store dose data
    dose_data{end+1} = dose_volume;
    dose_labels{end+1} = strrep(dose_file.name, '.dcm', '');
    
    fprintf('  Dose volume size: %d x %d x %d\n', ...
        size(dose_volume, 1), size(dose_volume, 2), size(dose_volume, 3));
    fprintf('  Max dose: %.2f Gy\n', max(dose_volume(:)));
end

if isempty(dose_data)
    error('No dose files loaded.');
end

fprintf('\nLoaded %d dose distributions\n', length(dose_data));

%% Create Visualization

fprintf('\nCreating visualizations...\n');

% Create figure for each dose distribution
for dose_idx = 1:length(dose_data)
    
    figure('Name', dose_labels{dose_idx}, 'Position', [100, 100, 600, 500]);
    
    dose_vol = dose_data{dose_idx};
    
    % Find slice with maximum dose
    max_dose_per_slice = squeeze(max(max(dose_vol, [], 1), [], 2));
    [max_dose, max_slice_idx] = max(max_dose_per_slice);
    
    fprintf('Dose file %d (%s): Max dose = %.2f Gy at slice %d\n', ...
        dose_idx, dose_labels{dose_idx}, max_dose, max_slice_idx);
    
    % Get the slice with maximum dose
    dose_slice = dose_vol(:, :, max_slice_idx);
    
    % Display dose
    imagesc(dose_slice);
    colormap(dose_colormap);
    axis image;
    
    % Add colorbar
    cb = colorbar;
    ylabel(cb, 'Dose (Gy)', 'FontSize', 11);
    
    % Labels and title
    title(sprintf('%s\nSlice %d - Max Dose: %.2f Gy', ...
        dose_labels{dose_idx}, max_slice_idx, max_dose), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel('X (pixels)', 'FontSize', 10);
    ylabel('Y (pixels)', 'FontSize', 10);
    
    % Add grid for better readability
    grid on;
    set(gca, 'GridAlpha', 0.3, 'GridColor', 'w');
end

%% Create Comparison Figure (all doses, max dose slices)
if length(dose_data) > 1
    
    figure('Name', 'Dose Comparison - Max Dose Slices', ...
        'Position', [100, 100, 400*min(length(dose_data), 3), 400]);
    
    % Determine layout
    n_doses = length(dose_data);
    n_cols = min(n_doses, 3);
    n_rows = ceil(n_doses / n_cols);
    
    for dose_idx = 1:length(dose_data)
        subplot(n_rows, n_cols, dose_idx);
        
        dose_vol = dose_data{dose_idx};
        
        % Find slice with maximum dose
        max_dose_per_slice = squeeze(max(max(dose_vol, [], 1), [], 2));
        [max_dose, max_slice_idx] = max(max_dose_per_slice);
        
        % Get the slice with maximum dose
        dose_slice = dose_vol(:, :, max_slice_idx);
        
        % Display
        imagesc(dose_slice);
        colormap(dose_colormap);
        axis image;
        
        cb = colorbar;
        ylabel(cb, 'Dose (Gy)', 'FontSize', 9);
        
        title(sprintf('%s\nSlice %d: %.2f Gy', ...
            dose_labels{dose_idx}, max_slice_idx, max_dose), ...
            'FontSize', 10, 'Interpreter', 'none');
        
        grid on;
        set(gca, 'GridAlpha', 0.3, 'GridColor', 'w');
    end
    
    sgtitle('Dose Comparison - Maximum Dose Slices', ...
        'FontSize', 13, 'FontWeight', 'bold');
end

fprintf('\nVisualization complete!\n');

%% ========================================================================
%  Helper Functions
%% ========================================================================

function [dose_volume, dose_info] = load_rtdose(dose_filepath)
    % Load RT Dose DICOM file
    
    dose_info = dicominfo(dose_filepath);
    dose_raw = dicomread(dose_filepath);
    
    % Apply dose grid scaling
    if isfield(dose_info, 'DoseGridScaling')
        dose_volume = double(dose_raw) * dose_info.DoseGridScaling;
    else
        dose_volume = double(dose_raw);
    end
    
    % DICOM RT Dose may need permutation depending on storage order
    % Typically stored as [cols, rows, slices], need [rows, cols, slices]
    if ndims(dose_volume) == 3
        dose_volume = permute(dose_volume, [2, 1, 3]);
    end
end