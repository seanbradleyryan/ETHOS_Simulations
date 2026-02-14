%% VISUALIZE_PROCESSED_DATA
%
%   Script to visualize processed dose data and tissue masks at the 
%   maximum dose location. Shows sagittal, coronal, and transverse views
%   along with the body/region mask overlay.
%
%   USAGE:
%       1. Set the data_file path below to your processed data file
%       2. Run the script
%
%   OUTPUTS:
%       Figure 1: Three orthogonal views of dose at max dose location
%       Figure 2: Three orthogonal views of tissue/body mask at max dose location
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0

%% ======================== CONFIGURATION ========================

% Default file path - modify as needed
% Pattern: RayStationFiles/[PatientID]/[Session]/processed/
data_file = '/mnt/weka/home/80030361/ETHOS_Simulations/RayStationFiles/1194203/Session_1/processed/sct_resampled.mat';
dose_file = '/mnt/weka/home/80030361/ETHOS_Simulations/RayStationFiles/1194203/Session_1/processed/total_rs_dose.mat';

% Alternative: load a specific field dose
% dose_file = '/mnt/weka/home/80030361/ETHOS_Simulations/RayStationFiles/1194203/Session_1/processed/field_dose_001.mat';

%% ======================== LOAD DATA ========================

fprintf('Loading data files...\n');

% Load SCT resampled data (contains tissue masks)
if ~isfile(data_file)
    error('Data file not found: %s', data_file);
end
loaded = load(data_file);
sct_resampled = loaded.sct_resampled;
fprintf('  Loaded: %s\n', data_file);

% Load dose data
if ~isfile(dose_file)
    error('Dose file not found: %s', dose_file);
end
loaded = load(dose_file);
if isfield(loaded, 'total_rs_dose')
    dose_data = loaded.total_rs_dose;
    dose_label = 'Total Dose';
elseif isfield(loaded, 'field_dose')
    dose_data = loaded.field_dose.dose_Gy;
    dose_label = sprintf('Field %d Dose', loaded.field_dose.field_num);
else
    error('Could not find dose data in file');
end
fprintf('  Loaded: %s\n', dose_file);

%% ======================== FIND MAX DOSE LOCATION ========================

fprintf('\nFinding maximum dose location...\n');

[max_dose, max_idx] = max(dose_data(:));
[max_row, max_col, max_slice] = ind2sub(size(dose_data), max_idx);

fprintf('  Max dose: %.4f Gy\n', max_dose);
fprintf('  Location (row, col, slice): [%d, %d, %d]\n', max_row, max_col, max_slice);

% Get grid info
dims = size(dose_data);
spacing = sct_resampled.spacing;
origin = sct_resampled.origin;

fprintf('  Grid dimensions: [%d, %d, %d]\n', dims(1), dims(2), dims(3));
fprintf('  Spacing (mm): [%.2f, %.2f, %.2f]\n', spacing(1), spacing(2), spacing(3));

%% ======================== EXTRACT SLICES AT MAX DOSE ========================

% Transverse (axial) slice - fixed z
transverse_dose = squeeze(dose_data(:, :, max_slice));
transverse_hu = squeeze(sct_resampled.cubeHU(:, :, max_slice));

% Coronal slice - fixed y (row)
coronal_dose = squeeze(dose_data(max_row, :, :));
coronal_hu = squeeze(sct_resampled.cubeHU(max_row, :, :));

% Sagittal slice - fixed x (column)
sagittal_dose = squeeze(dose_data(:, max_col, :));
sagittal_hu = squeeze(sct_resampled.cubeHU(:, max_col, :));

%% ======================== FIGURE 1: DOSE VIEWS ========================

figure('Name', 'Dose Distribution at Max Dose Location', ...
    'Position', [50, 100, 1400, 500], 'Color', 'w');

% Transverse view
subplot(1, 3, 1);
imagesc(transverse_hu);
hold on;
contour(transverse_dose, 5, 'LineWidth', 1.5);
plot(max_col, max_row, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;
colormap(gca, gray);
colorbar;
axis equal tight;
title(sprintf('Transverse (Slice %d)', max_slice), 'FontSize', 12);
xlabel('Column (X)');
ylabel('Row (Y)');

% Coronal view
subplot(1, 3, 2);
imagesc(coronal_hu');
hold on;
contour(coronal_dose', 5, 'LineWidth', 1.5);
plot(max_col, max_slice, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;
colormap(gca, gray);
colorbar;
axis equal tight;
title(sprintf('Coronal (Row %d)', max_row), 'FontSize', 12);
xlabel('Column (X)');
ylabel('Slice (Z)');

% Sagittal view
subplot(1, 3, 3);
imagesc(sagittal_hu');
hold on;
contour(sagittal_dose', 5, 'LineWidth', 1.5);
plot(max_row, max_slice, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;
colormap(gca, gray);
colorbar;
axis equal tight;
title(sprintf('Sagittal (Col %d)', max_col), 'FontSize', 12);
xlabel('Row (Y)');
ylabel('Slice (Z)');

sgtitle(sprintf('%s - Max: %.4f Gy at [%d, %d, %d]', dose_label, max_dose, max_row, max_col, max_slice), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ======================== FIGURE 2: BODY/TISSUE MASK ========================

% Check for body mask
if isfield(sct_resampled, 'bodyMask')
    mask_data = sct_resampled.bodyMask;
    mask_label = 'Body Mask';
elseif isfield(sct_resampled, 'tissueMask')
    mask_data = sct_resampled.tissueMask > 0;  % Convert to logical
    mask_label = 'Tissue Mask';
else
    warning('No body or tissue mask found in data');
    mask_data = [];
end

if ~isempty(mask_data)
    % Extract mask slices
    transverse_mask = squeeze(mask_data(:, :, max_slice));
    coronal_mask = squeeze(mask_data(max_row, :, :));
    sagittal_mask = squeeze(mask_data(:, max_col, :));
    
    figure('Name', 'Body/Region Mask at Max Dose Location', ...
        'Position', [100, 150, 1400, 500], 'Color', 'w');
    
    % Transverse view
    subplot(1, 3, 1);
    % Overlay mask on HU
    rgb_img = createMaskOverlay(transverse_hu, transverse_mask);
    imshow(rgb_img);
    hold on;
    plot(max_col, max_row, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    title(sprintf('Transverse (Slice %d)', max_slice), 'FontSize', 12);
    xlabel('Column (X)');
    ylabel('Row (Y)');
    
    % Coronal view
    subplot(1, 3, 2);
    rgb_img = createMaskOverlay(coronal_hu', coronal_mask');
    imshow(rgb_img);
    hold on;
    plot(max_col, max_slice, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    title(sprintf('Coronal (Row %d)', max_row), 'FontSize', 12);
    xlabel('Column (X)');
    ylabel('Slice (Z)');
    
    % Sagittal view
    subplot(1, 3, 3);
    rgb_img = createMaskOverlay(sagittal_hu', sagittal_mask');
    imshow(rgb_img);
    hold on;
    plot(max_row, max_slice, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    title(sprintf('Sagittal (Col %d)', max_col), 'FontSize', 12);
    xlabel('Row (Y)');
    ylabel('Slice (Z)');
    
    sgtitle(sprintf('%s - Body voxels: %d / %d (%.1f%%)', mask_label, ...
        sum(mask_data(:)), numel(mask_data), 100*sum(mask_data(:))/numel(mask_data)), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

%% ======================== FIGURE 3: DOSE COLORMAP ========================

figure('Name', 'Dose with Colormap at Max Dose Location', ...
    'Position', [150, 200, 1400, 500], 'Color', 'w');

% Define dose colormap (hot-style)
dose_cmap = hot(256);

% Transverse view
subplot(1, 3, 1);
imagesc(transverse_dose);
hold on;
if ~isempty(mask_data)
    contour(transverse_mask, [0.5 0.5], 'g-', 'LineWidth', 1);
end
plot(max_col, max_row, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;
colormap(gca, dose_cmap);
cbar = colorbar;
ylabel(cbar, 'Dose (Gy)');
caxis([0, max_dose]);
axis equal tight;
title(sprintf('Transverse (Slice %d)', max_slice), 'FontSize', 12);
xlabel('Column (X)');
ylabel('Row (Y)');

% Coronal view
subplot(1, 3, 2);
imagesc(coronal_dose');
hold on;
if ~isempty(mask_data)
    contour(coronal_mask', [0.5 0.5], 'g-', 'LineWidth', 1);
end
plot(max_col, max_slice, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;
colormap(gca, dose_cmap);
cbar = colorbar;
ylabel(cbar, 'Dose (Gy)');
caxis([0, max_dose]);
axis equal tight;
title(sprintf('Coronal (Row %d)', max_row), 'FontSize', 12);
xlabel('Column (X)');
ylabel('Slice (Z)');

% Sagittal view
subplot(1, 3, 3);
imagesc(sagittal_dose');
hold on;
if ~isempty(mask_data)
    contour(sagittal_mask', [0.5 0.5], 'g-', 'LineWidth', 1);
end
plot(max_row, max_slice, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;
colormap(gca, dose_cmap);
cbar = colorbar;
ylabel(cbar, 'Dose (Gy)');
caxis([0, max_dose]);
axis equal tight;
title(sprintf('Sagittal (Col %d)', max_col), 'FontSize', 12);
xlabel('Row (Y)');
ylabel('Slice (Z)');

sgtitle(sprintf('%s Colormap - Max: %.4f Gy (green = body contour)', dose_label, max_dose), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nVisualization complete.\n');

%% ======================== HELPER FUNCTIONS ========================

function rgb_img = createMaskOverlay(hu_slice, mask_slice)
%CREATEMASKOVERLAY Create RGB image with mask overlay on HU background
%
%   Normalizes HU to grayscale and overlays mask in semi-transparent color

    % Normalize HU to [0, 1]
    hu_min = -1000;
    hu_max = 1000;
    hu_norm = (double(hu_slice) - hu_min) / (hu_max - hu_min);
    hu_norm = max(0, min(1, hu_norm));
    
    % Create grayscale RGB
    rgb_img = repmat(hu_norm, [1, 1, 3]);
    
    % Overlay mask in green with transparency
    mask_color = [0, 0.7, 0.3];  % Green tint
    alpha = 0.3;  % Transparency
    
    for c = 1:3
        channel = rgb_img(:, :, c);
        channel(mask_slice) = (1 - alpha) * channel(mask_slice) + alpha * mask_color(c);
        rgb_img(:, :, c) = channel;
    end
end