%% VISUALIZE_PROCESSED_DATA
%
%   Script to visualize processed dose data and tissue masks at the 
%   maximum dose location. Shows sagittal, coronal, and transverse views
%   along with the body/region mask overlay.
%
%   Also includes diagnostics for detecting streaking/gaps in body mask.
%
%   USAGE:
%       1. Set the data_file path below to your processed data file
%       2. Run the script
%
%   OUTPUTS:
%       Figure 1: Three orthogonal views of dose at max dose location
%       Figure 2: Three orthogonal views of tissue/body mask at max dose location
%       Figure 3: Dose colormap with body contour
%       Figure 4: Body mask diagnostics (slice coverage, gaps detection)
%       Figure 5: Side-by-side original field dose vs total RS dose
%       Figure 6: Gamma analysis (3%/3mm) - gamma map, pass/fail, histogram
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.2 (Added original dose comparison and gamma analysis)

%% ======================== CONFIGURATION ========================

% Default file path - modify as needed
% Pattern: RayStationFiles/[PatientID]/[Session]/processed/
data_file = '/mnt/weka/home/80030361/ETHOS_Simulations/RayStationFiles/1194203/Session_1/processed/sct_resampled.mat';
dose_file = '/mnt/weka/home/80030361/ETHOS_Simulations/RayStationFiles/1194203/Session_1/processed/total_rs_dose.mat';

% Alternative: load a specific field dose
% dose_file = '/mnt/weka/home/80030361/ETHOS_Simulations/RayStationFiles/1194203/Session_1/processed/field_dose_001.mat';

% Original field dose from ETHOS (DICOM) for gamma comparison
original_dose_dcm = '/mnt/weka/home/80030361/ETHOS_Simulations/RayStationFiles/1194203/Session_1/processed/original_field_dose.dcm';

% Gamma analysis parameters
gamma_percent = 3;   % Dose difference criterion (%)
gamma_dta_mm  = 3;   % Distance-to-agreement criterion (mm)
gamma_dose_threshold = 10;  % % of max dose below which to exclude from pass rate

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

% Load original field dose from DICOM
if isfile(original_dose_dcm)
    fprintf('  Loading original field dose DICOM...\n');
    orig_info = dicominfo(original_dose_dcm);
    orig_dose_raw = double(squeeze(dicomread(original_dose_dcm)));
    
    % Apply dose grid scaling
    if isfield(orig_info, 'DoseGridScaling')
        orig_dose_data = orig_dose_raw * orig_info.DoseGridScaling;
    else
        orig_dose_data = orig_dose_raw;
        warning('No DoseGridScaling found in original DICOM, using raw values');
    end
    
    % Extract geometry
    orig_origin = orig_info.ImagePositionPatient(:);  % [x, y, z] in mm
    orig_pixel_spacing = orig_info.PixelSpacing;      % [row, col] spacing
    if isfield(orig_info, 'GridFrameOffsetVector') && length(orig_info.GridFrameOffsetVector) >= 2
        orig_dz = abs(orig_info.GridFrameOffsetVector(2) - orig_info.GridFrameOffsetVector(1));
    elseif isfield(orig_info, 'SliceThickness')
        orig_dz = orig_info.SliceThickness;
    else
        orig_dz = spacing(3);  % Fallback to dose grid spacing
        warning('Cannot determine Z spacing for original dose, using %.3f mm', orig_dz);
    end
    orig_spacing = [orig_pixel_spacing(1); orig_pixel_spacing(2); orig_dz];
    
    fprintf('    Dimensions: [%d, %d, %d]\n', size(orig_dose_data, 1), size(orig_dose_data, 2), size(orig_dose_data, 3));
    fprintf('    Spacing (mm): [%.3f, %.3f, %.3f]\n', orig_spacing(1), orig_spacing(2), orig_spacing(3));
    fprintf('    Origin (mm): [%.3f, %.3f, %.3f]\n', orig_origin(1), orig_origin(2), orig_origin(3));
    fprintf('    Max dose: %.4f Gy\n', max(orig_dose_data(:)));
    
    has_original_dose = true;
else
    warning('Original field dose DICOM not found: %s', original_dose_dcm);
    has_original_dose = false;
end

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

%% ======================== FIGURE 4: BODY MASK DIAGNOSTICS ========================
% This figure helps diagnose streaking artifacts in the body mask

if ~isempty(mask_data)
    figure('Name', 'Body Mask Diagnostics - Streaking Analysis', ...
        'Position', [200, 50, 1400, 800], 'Color', 'w');
    
    % Calculate per-slice statistics
    num_slices = size(mask_data, 3);
    voxels_per_slice = zeros(num_slices, 1);
    for k = 1:num_slices
        voxels_per_slice(k) = sum(sum(mask_data(:, :, k)));
    end
    
    % Detect gaps (slices with significantly fewer voxels than neighbors)
    median_voxels = median(voxels_per_slice(voxels_per_slice > 0));
    gap_threshold = 0.5 * median_voxels;  % Less than 50% of median
    potential_gaps = find(voxels_per_slice > 0 & voxels_per_slice < gap_threshold);
    
    % Detect complete missing slices within body range
    first_slice = find(voxels_per_slice > 0, 1, 'first');
    last_slice = find(voxels_per_slice > 0, 1, 'last');
    missing_slices = find(voxels_per_slice(first_slice:last_slice) == 0) + first_slice - 1;
    
    % Subplot 1: Voxels per slice
    subplot(2, 3, 1);
    bar(voxels_per_slice);
    hold on;
    if ~isempty(potential_gaps)
        plot(potential_gaps, voxels_per_slice(potential_gaps), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    if ~isempty(missing_slices)
        plot(missing_slices, zeros(size(missing_slices)), 'rx', 'MarkerSize', 15, 'LineWidth', 2);
    end
    hold off;
    xlabel('Slice Index');
    ylabel('Voxels in Body Mask');
    title('Body Mask Voxels per Slice');
    legend_entries = {'Voxel count'};
    if ~isempty(potential_gaps)
        legend_entries{end+1} = 'Potential gaps';
    end
    if ~isempty(missing_slices)
        legend_entries{end+1} = 'Missing slices';
    end
    legend(legend_entries, 'Location', 'best');
    grid on;
    
    % Subplot 2: Slice-to-slice change (derivative)
    subplot(2, 3, 2);
    slice_diff = diff(voxels_per_slice);
    plot(2:num_slices, slice_diff, 'b-', 'LineWidth', 1);
    hold on;
    % Highlight large jumps
    large_jumps = find(abs(slice_diff) > 0.3 * median_voxels);
    if ~isempty(large_jumps)
        plot(large_jumps + 1, slice_diff(large_jumps), 'ro', 'MarkerSize', 8);
    end
    hold off;
    xlabel('Slice Index');
    ylabel('Change in Voxel Count');
    title('Slice-to-Slice Change (Detects Discontinuities)');
    grid on;
    
    % Subplot 3: Sagittal MIP of body mask
    subplot(2, 3, 3);
    sagittal_mip = squeeze(max(mask_data, [], 2));  % Max projection along columns
    imagesc(sagittal_mip');
    colormap(gca, gray);
    title('Sagittal MIP of Body Mask');
    xlabel('Row (Y)');
    ylabel('Slice (Z)');
    axis equal tight;
    
    % Subplot 4: Coronal MIP of body mask
    subplot(2, 3, 4);
    coronal_mip = squeeze(max(mask_data, [], 1));  % Max projection along rows
    imagesc(coronal_mip');
    colormap(gca, gray);
    title('Coronal MIP of Body Mask');
    xlabel('Column (X)');
    ylabel('Slice (Z)');
    axis equal tight;
    
    % Subplot 5: Example slices around a gap (if any)
    subplot(2, 3, 5);
    if ~isempty(potential_gaps)
        gap_slice = potential_gaps(1);
        % Show slice before, gap, and after
        slices_to_show = max(1, gap_slice-1):min(num_slices, gap_slice+1);
        montage_data = mask_data(:, :, slices_to_show);
        montage(montage_data, 'Size', [1, length(slices_to_show)]);
        title(sprintf('Slices around gap (slice %d)', gap_slice));
    elseif ~isempty(missing_slices)
        gap_slice = missing_slices(1);
        slices_to_show = max(1, gap_slice-1):min(num_slices, gap_slice+1);
        montage_data = mask_data(:, :, slices_to_show);
        montage(montage_data, 'Size', [1, length(slices_to_show)]);
        title(sprintf('Slices around missing (slice %d)', gap_slice));
    else
        imagesc(squeeze(mask_data(:, :, max_slice)));
        colormap(gca, gray);
        title(sprintf('Slice %d (no gaps detected)', max_slice));
    end
    axis equal tight;
    
    % Subplot 6: Diagnostic summary text
    subplot(2, 3, 6);
    axis off;
    
    % Build diagnostic text
    diag_text = {
        '=== BODY MASK DIAGNOSTICS ===', ...
        '', ...
        sprintf('Grid dimensions: [%d, %d, %d]', dims(1), dims(2), dims(3)), ...
        sprintf('Grid spacing (mm): [%.2f, %.2f, %.2f]', spacing(1), spacing(2), spacing(3)), ...
        '', ...
        sprintf('Total body voxels: %d', sum(mask_data(:))), ...
        sprintf('Body slice range: %d to %d', first_slice, last_slice), ...
        sprintf('Median voxels/slice: %.0f', median_voxels), ...
        '', ...
        '--- POTENTIAL ISSUES ---'
    };
    
    if isempty(potential_gaps) && isempty(missing_slices)
        diag_text{end+1} = 'No obvious gaps or missing slices detected.';
        diag_text{end+1} = '';
        diag_text{end+1} = 'If streaking is visible, possible causes:';
        diag_text{end+1} = '  1. RTSTRUCT contours have gaps in original data';
        diag_text{end+1} = '  2. Contour z-coords misaligned with dose grid';
        diag_text{end+1} = '  3. CT slice spacing differs from dose spacing';
        diag_text{end+1} = '';
        diag_text{end+1} = 'Check: Load RTSTRUCT directly and inspect';
        diag_text{end+1} = 'contour z-coordinates vs dose grid z-coords.';
    else
        if ~isempty(missing_slices)
            diag_text{end+1} = sprintf('Missing slices: %d', length(missing_slices));
            diag_text{end+1} = sprintf('  Indices: %s', mat2str(missing_slices(1:min(10,end))));
        end
        if ~isempty(potential_gaps)
            diag_text{end+1} = sprintf('Partial gaps: %d slices', length(potential_gaps));
            diag_text{end+1} = sprintf('  Indices: %s', mat2str(potential_gaps(1:min(10,end))));
        end
        diag_text{end+1} = '';
        diag_text{end+1} = 'LIKELY CAUSE: Contour z-coordinates do not';
        diag_text{end+1} = 'align with dose grid slice positions.';
        diag_text{end+1} = '';
        diag_text{end+1} = 'SOLUTION: Interpolate body mask in z-direction';
        diag_text{end+1} = 'or adjust slice matching tolerance.';
    end
    
    text(0.05, 0.95, diag_text, 'FontSize', 10, 'FontName', 'FixedWidth', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    sgtitle('Body Mask Diagnostics - Check for Streaking/Gaps', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Print diagnostic summary to console
    fprintf('\n=== BODY MASK DIAGNOSTICS ===\n');
    fprintf('Body slice range: %d to %d (of %d total)\n', first_slice, last_slice, num_slices);
    fprintf('Median voxels per slice: %.0f\n', median_voxels);
    if ~isempty(missing_slices)
        fprintf('WARNING: %d completely missing slices within body range\n', length(missing_slices));
        fprintf('  First few: %s\n', mat2str(missing_slices(1:min(5,end))));
    end
    if ~isempty(potential_gaps)
        fprintf('WARNING: %d slices with partial gaps (< 50%% of median)\n', length(potential_gaps));
        fprintf('  First few: %s\n', mat2str(potential_gaps(1:min(5,end))));
    end
    if isempty(missing_slices) && isempty(potential_gaps)
        fprintf('No obvious gaps detected in slice coverage.\n');
        fprintf('If streaking visible, the issue is likely in the original RTSTRUCT data.\n');
    end
end

fprintf('\nVisualization complete.\n');

%% ======================== FIGURE 5: ORIGINAL vs TOTAL DOSE COMPARISON ========================

if has_original_dose
    fprintf('\n--- Original Field Dose vs Total RS Dose Comparison ---\n');
    
    % Check if grids match; if not, resample original dose to the RS dose grid
    rs_dims = size(dose_data);
    rs_origin = sct_resampled.origin;
    rs_spacing = sct_resampled.spacing;
    
    if isequal(size(orig_dose_data), rs_dims)
        orig_dose_on_grid = orig_dose_data;
        fprintf('  Grid dimensions match - no resampling needed.\n');
    else
        fprintf('  Grid mismatch - resampling original dose to RS dose grid...\n');
        orig_dose_on_grid = resampleOrigDoseToGrid(orig_dose_data, orig_origin, orig_spacing, ...
            rs_origin, rs_spacing, rs_dims);
        fprintf('  Resampling complete. New dimensions: [%d, %d, %d]\n', size(orig_dose_on_grid));
    end
    
    % Find max dose location on original dose for slice selection
    [orig_max_dose, orig_max_idx] = max(orig_dose_on_grid(:));
    [orig_max_row, orig_max_col, orig_max_slice] = ind2sub(size(orig_dose_on_grid), orig_max_idx);
    
    % Use the slice of max dose from whichever has higher max
    if orig_max_dose >= max_dose
        comp_slice = orig_max_slice;
        comp_row = orig_max_row;
        comp_col = orig_max_col;
    else
        comp_slice = max_slice;
        comp_row = max_row;
        comp_col = max_col;
    end
    
    % Common color scale
    common_max = max(orig_max_dose, max_dose);
    
    % Extract transverse slices for comparison
    orig_trans = squeeze(orig_dose_on_grid(:, :, comp_slice));
    rs_trans = squeeze(dose_data(:, :, comp_slice));
    
    figure('Name', 'Original Field Dose vs Total RS Dose', ...
        'Position', [50, 50, 1200, 550], 'Color', 'w');
    
    dose_cmap_comp = hot(256);
    
    % Subplot 1: Original field dose
    subplot(1, 2, 1);
    imagesc(orig_trans);
    hold on;
    if ~isempty(mask_data)
        contour(squeeze(mask_data(:, :, comp_slice)), [0.5 0.5], 'g-', 'LineWidth', 1);
    end
    plot(comp_col, comp_row, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    colormap(gca, dose_cmap_comp);
    cbar = colorbar;
    ylabel(cbar, 'Dose (Gy)');
    caxis([0, common_max]);
    axis equal tight;
    title(sprintf('Original Field Dose (Max: %.4f Gy)', orig_max_dose), 'FontSize', 12);
    xlabel('Column (X)');
    ylabel('Row (Y)');
    
    % Subplot 2: Total RS dose
    subplot(1, 2, 2);
    imagesc(rs_trans);
    hold on;
    if ~isempty(mask_data)
        contour(squeeze(mask_data(:, :, comp_slice)), [0.5 0.5], 'g-', 'LineWidth', 1);
    end
    plot(comp_col, comp_row, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    colormap(gca, dose_cmap_comp);
    cbar = colorbar;
    ylabel(cbar, 'Dose (Gy)');
    caxis([0, common_max]);
    axis equal tight;
    title(sprintf('Total RS Dose (Max: %.4f Gy)', max_dose), 'FontSize', 12);
    xlabel('Column (X)');
    ylabel('Row (Y)');
    
    sgtitle(sprintf('Dose Comparison - Transverse Slice %d (green = body contour)', comp_slice), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    %% ======================== FIGURE 6: GAMMA ANALYSIS ========================
    
    fprintf('\n--- Gamma Analysis: %d%%/%dmm ---\n', gamma_percent, gamma_dta_mm);
    fprintf('  Reference: Original Field Dose\n');
    fprintf('  Evaluated: Total RS Dose\n');
    
    % Build CalcGamma input structures
    % Reference = original field dose
    reference.start = rs_origin(:)';   % [x, y, z] in mm
    reference.width = rs_spacing(:)';  % [dx, dy, dz] in mm
    reference.data  = orig_dose_on_grid;
    
    % Target = total RS dose
    target.start = rs_origin(:)';
    target.width = rs_spacing(:)';
    target.data  = dose_data;
    
    % Run 3D gamma (restrict=1 for speed, global gamma)
    fprintf('  Computing 3D gamma (this may take a while)...\n');
    tic;
    gamma_map = CalcGamma(reference, target, gamma_percent, gamma_dta_mm, ...
        'local', 0, 'restrict', 1);
    gamma_time = toc;
    fprintf('  Gamma computation completed in %.1f seconds.\n', gamma_time);
    
    % Calculate pass rate (exclude low-dose region)
    dose_threshold_Gy = (gamma_dose_threshold / 100) * max(orig_dose_on_grid(:));
    valid_mask = orig_dose_on_grid >= dose_threshold_Gy;
    gamma_valid = gamma_map(valid_mask);
    pass_rate = 100 * sum(gamma_valid <= 1) / numel(gamma_valid);
    
    fprintf('\n  === GAMMA RESULTS ===\n');
    fprintf('  Criteria: %d%%/%dmm (global)\n', gamma_percent, gamma_dta_mm);
    fprintf('  Dose threshold: %.4f Gy (%d%% of max)\n', dose_threshold_Gy, gamma_dose_threshold);
    fprintf('  Voxels evaluated: %d / %d\n', sum(valid_mask(:)), numel(valid_mask));
    fprintf('  Pass rate (gamma <= 1): %.2f%%\n', pass_rate);
    fprintf('  Mean gamma: %.4f\n', mean(gamma_valid));
    fprintf('  Max gamma: %.4f\n', max(gamma_valid));
    fprintf('  Median gamma: %.4f\n', median(gamma_valid));
    
    % Extract gamma at comparison slice
    gamma_trans = squeeze(gamma_map(:, :, comp_slice));
    
    figure('Name', 'Gamma Analysis: Original vs Total RS Dose', ...
        'Position', [100, 50, 1600, 550], 'Color', 'w');
    
    % Subplot 1: Gamma map (transverse)
    subplot(1, 3, 1);
    imagesc(gamma_trans);
    hold on;
    if ~isempty(mask_data)
        contour(squeeze(mask_data(:, :, comp_slice)), [0.5 0.5], 'k-', 'LineWidth', 1);
    end
    hold off;
    colormap(gca, jet(256));
    cbar = colorbar;
    ylabel(cbar, 'Gamma Index');
    caxis([0, 2]);
    axis equal tight;
    title(sprintf('Gamma Map (Slice %d)', comp_slice), 'FontSize', 12);
    xlabel('Column (X)');
    ylabel('Row (Y)');
    
    % Subplot 2: Pass/fail map (transverse)
    subplot(1, 3, 2);
    slice_mask = squeeze(valid_mask(:, :, comp_slice));
    pass_fail_img = nan(size(gamma_trans));
    pass_fail_img(slice_mask & gamma_trans <= 1) = 1;   % Pass = green
    pass_fail_img(slice_mask & gamma_trans > 1)  = 0;   % Fail = red
    
    imagesc(pass_fail_img);
    hold on;
    if ~isempty(mask_data)
        contour(squeeze(mask_data(:, :, comp_slice)), [0.5 0.5], 'k-', 'LineWidth', 1);
    end
    hold off;
    
    % Custom pass/fail colormap: red -> green
    pf_cmap = [1 0 0; 0 0.8 0];
    colormap(gca, pf_cmap);
    caxis([0, 1]);
    cbar = colorbar;
    set(cbar, 'Ticks', [0.25, 0.75], 'TickLabels', {'FAIL', 'PASS'});
    axis equal tight;
    
    % Slice pass rate
    gamma_slice_valid = gamma_trans(slice_mask);
    slice_pass_rate = 100 * sum(gamma_slice_valid <= 1) / numel(gamma_slice_valid);
    title(sprintf('Pass/Fail (Slice: %.1f%%)', slice_pass_rate), 'FontSize', 12);
    xlabel('Column (X)');
    ylabel('Row (Y)');
    
    % Subplot 3: Gamma histogram
    subplot(1, 3, 3);
    histogram(gamma_valid, 100, 'FaceColor', [0.3 0.5 0.8], 'EdgeAlpha', 0.3);
    hold on;
    xline(1, 'r--', 'LineWidth', 2);
    hold off;
    xlabel('Gamma Index');
    ylabel('Voxel Count');
    title('Gamma Histogram (Thresholded Volume)', 'FontSize', 12);
    
    % Add annotation with stats
    annotation_text = sprintf('Pass rate: %.2f%%\nMean: %.3f\nMedian: %.3f\nMax: %.3f', ...
        pass_rate, mean(gamma_valid), median(gamma_valid), max(gamma_valid));
    text(0.95, 0.95, annotation_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
        'FontSize', 10, 'FontName', 'FixedWidth', 'BackgroundColor', 'w', ...
        'EdgeColor', 'k');
    grid on;
    xlim([0, max(2, prctile(gamma_valid, 99))]);
    
    sgtitle(sprintf('Gamma Analysis: %d%%/%dmm - Overall Pass Rate: %.2f%%', ...
        gamma_percent, gamma_dta_mm, pass_rate), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
else
    fprintf('\nSkipping dose comparison and gamma analysis (no original dose loaded).\n');
end

fprintf('\nAll visualization and analysis complete.\n');

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

function resampled = resampleOrigDoseToGrid(orig_data, orig_origin, orig_spacing, ...
    target_origin, target_spacing, target_dims)
%RESAMPLEORIGDOSETOGRID Resample original dose to target grid via 3D interpolation
%
%   Uses trilinear interpolation to resample the original dose distribution
%   onto the target (RS dose) grid. Out-of-bounds voxels are set to zero.

    orig_dims = size(orig_data);
    
    % Build coordinate vectors for original data grid
    orig_x = orig_origin(1) + (0:orig_dims(2)-1) * orig_spacing(1);
    orig_y = orig_origin(2) + (0:orig_dims(1)-1) * orig_spacing(2);
    orig_z = orig_origin(3) + (0:orig_dims(3)-1) * orig_spacing(3);
    
    % Build coordinate vectors for target grid
    tar_x = target_origin(1) + (0:target_dims(2)-1) * target_spacing(1);
    tar_y = target_origin(2) + (0:target_dims(1)-1) * target_spacing(2);
    tar_z = target_origin(3) + (0:target_dims(3)-1) * target_spacing(3);
    
    % Create meshgrids for target
    [TarX, TarY, TarZ] = meshgrid(tar_x, tar_y, tar_z);
    
    % Interpolate
    resampled = interp3(orig_x, orig_y, orig_z, double(orig_data), ...
        TarX, TarY, TarZ, 'linear', 0);
end