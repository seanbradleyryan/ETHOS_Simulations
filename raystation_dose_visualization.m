%% visualize_rtdose_overlay.m
% Script to visualize exported RT dose files overlaid on planning CT
% For analyzing individual beam/segment doses from RayStation exports
%
% This script:
% 1. Loads planning CT (SCT) DICOM images
% 2. Loads selected RT dose files
% 3. Overlays dose distributions on CT slices
% 4. Creates multi-panel figures for comparison
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

% Path to exported dose files and CT
dose_folder = fullfile(wd, 'RayStationFiles', patient_id, session);
sct_folder = fullfile(wd, 'EthosExports', patient_id, 'Pancreas', session, 'sct');

% Number of dose files to randomly select and plot
num_doses_to_plot = 3;  % Pick a few dose files to visualize

% Visualization parameters
slice_indices = [50, 75, 100];  % CT slice indices to display (modify based on your data)
dose_colormap = 'jet';           % Colormap for dose overlay
dose_alpha = 0.6;                % Transparency of dose overlay (0-1)
dose_threshold = 0.1;            % Minimum dose to display (fraction of max dose)

%% Load Planning CT
fprintf('Loading planning CT from: %s\n', sct_folder);

% Find CT DICOM files
ct_files = dir(fullfile(sct_folder, '*.dcm'));
if isempty(ct_files)
    error('No DICOM files found in SCT folder');
end

% Read first file to get metadata
first_ct = dicominfo(fullfile(sct_folder, ct_files(1).name));

% Determine if this is CT or find CT files
ct_file_list = {};
for i = 1:length(ct_files)
    filepath = fullfile(sct_folder, ct_files(i).name);
    try
        info = dicominfo(filepath);
        % Check if it's a CT image
        if isfield(info, 'Modality') && strcmpi(info.Modality, 'CT')
            ct_file_list{end+1} = filepath;
        end
    catch
        continue;
    end
end

if isempty(ct_file_list)
    error('No CT images found in folder');
end

fprintf('Found %d CT slices\n', length(ct_file_list));

% Load CT volume
ct_volume = load_ct_volume(ct_file_list);
fprintf('CT volume size: %d x %d x %d\n', size(ct_volume, 1), size(ct_volume, 2), size(ct_volume, 3));

% Get CT metadata from first slice
ct_info = dicominfo(ct_file_list{1});
pixel_spacing = ct_info.PixelSpacing;  % [row, col] in mm
slice_thickness = ct_info.SliceThickness;  % in mm

fprintf('Pixel spacing: %.2f x %.2f mm\n', pixel_spacing(1), pixel_spacing(2));
fprintf('Slice thickness: %.2f mm\n', slice_thickness);

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
    error('No dose files loaded. Check file patterns and paths.');
end

fprintf('\nLoaded %d dose distributions\n', length(dose_data));

%% Create Visualization

% Adjust slice indices if needed
max_slices = size(ct_volume, 3);
slice_indices = slice_indices(slice_indices <= max_slices);

if isempty(slice_indices)
    slice_indices = round(linspace(1, max_slices, 3));
end

fprintf('\nCreating visualizations for slices: %s\n', num2str(slice_indices));

% Create figure for each dose distribution
for dose_idx = 1:length(dose_data)
    
    figure('Name', dose_labels{dose_idx}, 'Position', [100, 100, 1400, 400]);
    
    dose_vol = dose_data{dose_idx};
    max_dose = max(dose_vol(:));
    
    % Plot each slice
    for s = 1:length(slice_indices)
        subplot(1, length(slice_indices), s);
        
        slice_num = slice_indices(s);
        
        % Get CT slice
        ct_slice = ct_volume(:, :, slice_num);
        
        % Get dose slice (may need resampling if dimensions don't match)
        if size(dose_vol, 3) >= slice_num
            dose_slice = dose_vol(:, :, slice_num);
        else
            % Interpolate to match CT slice if dose has fewer slices
            z_ratio = size(dose_vol, 3) / size(ct_volume, 3);
            dose_slice_idx = round(slice_num * z_ratio);
            dose_slice_idx = max(1, min(dose_slice_idx, size(dose_vol, 3)));
            dose_slice = dose_vol(:, :, dose_slice_idx);
        end
        
        % Resample dose to CT grid if dimensions don't match
        if size(dose_slice, 1) ~= size(ct_slice, 1) || size(dose_slice, 2) ~= size(ct_slice, 2)
            [X_ct, Y_ct] = meshgrid(1:size(ct_slice, 2), 1:size(ct_slice, 1));
            [X_dose, Y_dose] = meshgrid(...
                linspace(1, size(ct_slice, 2), size(dose_slice, 2)), ...
                linspace(1, size(ct_slice, 1), size(dose_slice, 1)));
            dose_slice = interp2(X_dose, Y_dose, dose_slice, X_ct, Y_ct, 'linear', 0);
        end
        
        % Display CT
        imagesc(ct_slice);
        colormap(gca, 'gray');
        axis image;
        hold on;
        
        % Overlay dose with transparency
        dose_display = dose_slice;
        dose_display(dose_slice < dose_threshold * max_dose) = NaN;  % Threshold low doses
        
        h_dose = imagesc(dose_display, 'AlphaData', dose_alpha);
        colormap(h_dose.Parent, dose_colormap);
        
        % Add colorbar for dose
        cb = colorbar;
        ylabel(cb, 'Dose (Gy)', 'FontSize', 10);
        
        % Labels
        title(sprintf('Slice %d (z = %.1f mm)', slice_num, slice_num * slice_thickness), ...
            'FontSize', 11);
        xlabel('X (pixels)', 'FontSize', 9);
        ylabel('Y (pixels)', 'FontSize', 9);
        
        % Add dose statistics
        text(5, 15, sprintf('Max: %.2f Gy', max(dose_slice(:))), ...
            'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold', ...
            'BackgroundColor', [0 0 0 0.5]);
        
        hold off;
    end
    
    % Overall title
    sgtitle(sprintf('%s - Max dose: %.2f Gy', dose_labels{dose_idx}, max_dose), ...
        'FontSize', 12, 'FontWeight', 'bold');
end

%% Create Comparison Figure (all doses, single slice)
if length(dose_data) > 1
    
    comparison_slice = slice_indices(round(length(slice_indices)/2));  % Use middle slice
    
    figure('Name', 'Dose Comparison', 'Position', [100, 100, 400*length(dose_data), 400]);
    
    for dose_idx = 1:length(dose_data)
        subplot(1, length(dose_data), dose_idx);
        
        dose_vol = dose_data{dose_idx};
        max_dose = max(dose_vol(:));
        
        % Get slices
        ct_slice = ct_volume(:, :, comparison_slice);
        
        if size(dose_vol, 3) >= comparison_slice
            dose_slice = dose_vol(:, :, comparison_slice);
        else
            z_ratio = size(dose_vol, 3) / size(ct_volume, 3);
            dose_slice_idx = max(1, min(round(comparison_slice * z_ratio), size(dose_vol, 3)));
            dose_slice = dose_vol(:, :, dose_slice_idx);
        end
        
        % Resample if needed
        if size(dose_slice, 1) ~= size(ct_slice, 1) || size(dose_slice, 2) ~= size(ct_slice, 2)
            [X_ct, Y_ct] = meshgrid(1:size(ct_slice, 2), 1:size(ct_slice, 1));
            [X_dose, Y_dose] = meshgrid(...
                linspace(1, size(ct_slice, 2), size(dose_slice, 2)), ...
                linspace(1, size(ct_slice, 1), size(dose_slice, 1)));
            dose_slice = interp2(X_dose, Y_dose, dose_slice, X_ct, Y_ct, 'linear', 0);
        end
        
        % Display
        imagesc(ct_slice);
        colormap(gca, 'gray');
        axis image;
        hold on;
        
        dose_display = dose_slice;
        dose_display(dose_slice < dose_threshold * max_dose) = NaN;
        
        h_dose = imagesc(dose_display, 'AlphaData', dose_alpha);
        colormap(h_dose.Parent, dose_colormap);
        
        cb = colorbar;
        ylabel(cb, 'Dose (Gy)', 'FontSize', 9);
        
        title(dose_labels{dose_idx}, 'FontSize', 10, 'Interpreter', 'none');
        
        text(5, 15, sprintf('Max: %.2f Gy', max_dose), ...
            'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold', ...
            'BackgroundColor', [0 0 0 0.5]);
        
        hold off;
    end
    
    sgtitle(sprintf('Dose Comparison - Slice %d', comparison_slice), ...
        'FontSize', 12, 'FontWeight', 'bold');
end

fprintf('\nVisualization complete!\n');

%% ========================================================================
%  Helper Functions
%% ========================================================================

function ct_volume = load_ct_volume(ct_file_list)
    % Load CT DICOM files into 3D volume
    % Sorts by Image Position Patient (Z coordinate)
    
    num_slices = length(ct_file_list);
    
    % Read first slice to get dimensions
    first_slice = dicomread(ct_file_list{1});
    first_info = dicominfo(ct_file_list{1});
    
    [rows, cols] = size(first_slice);
    
    % Initialize volume
    ct_volume = zeros(rows, cols, num_slices);
    
    % Store slice positions for sorting
    slice_positions = zeros(num_slices, 1);
    
    % Load all slices
    for i = 1:num_slices
        info = dicominfo(ct_file_list{i});
        ct_volume(:, :, i) = dicomread(ct_file_list{i});
        
        % Get Z position
        if isfield(info, 'ImagePositionPatient')
            slice_positions(i) = info.ImagePositionPatient(3);
        else
            slice_positions(i) = i;  % Fallback to index
        end
    end
    
    % Sort slices by Z position
    [~, sort_idx] = sort(slice_positions);
    ct_volume = ct_volume(:, :, sort_idx);
end

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