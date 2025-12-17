%% ETHOS Field-by-Field Dose Visualization
% Purpose: Visualize and verify outputs from ETHOS IMRT Field-by-Field Dose Calculator
% Loads saved results and creates comprehensive plots for quality assurance
% Author: Generated for ETHOS dose analysis verification
% Date: 2025

clear; clc; close all;

%% Configuration
% Patient and session to visualize
patientID = '1194203';
sessionName = 'Session_1';

% Base directory (match original script)
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
dataPath = fullfile(wd, 'FieldDoses', patientID, sessionName);

% Verify data path exists
fprintf('=================================================\n');
fprintf('ETHOS Field Dose Visualization\n');
fprintf('=================================================\n');
fprintf('Patient ID: %s\n', patientID);
fprintf('Session: %s\n', sessionName);
fprintf('Data path: %s\n', dataPath);

if ~exist(dataPath, 'dir')
    error('Data directory does not exist: %s\nPlease run the dose calculation script first.', dataPath);
end
fprintf('Data directory found!\n\n');

%% Load Data
fprintf('[1/7] Loading saved data...\n');

% Load field doses
fieldDosesFile = fullfile(dataPath, 'fieldDoses.mat');
if ~exist(fieldDosesFile, 'file')
    error('fieldDoses.mat not found in: %s', dataPath);
end
load(fieldDosesFile, 'fieldDoses', 'stf', 'pln', 'ct', 'cst');
fprintf('  - Loaded: fieldDoses.mat\n');
fprintf('    Number of fields: %d\n', length(fieldDoses));

% Load reconstructed total dose
reconstructedFile = fullfile(dataPath, 'reconstructedDose.mat');
if exist(reconstructedFile, 'file')
    load(reconstructedFile, 'totalDose');
    fprintf('  - Loaded: reconstructedDose.mat\n');
    hasReconstructed = true;
else
    fprintf('  - Warning: reconstructedDose.mat not found\n');
    hasReconstructed = false;
end

% Load dose comparison (if available)
comparisonFile = fullfile(dataPath, 'doseComparison.mat');
if exist(comparisonFile, 'file')
    load(comparisonFile, 'comparison');
    fprintf('  - Loaded: doseComparison.mat\n');
    hasComparison = true;
else
    fprintf('  - Note: doseComparison.mat not found (comparison not available)\n');
    hasComparison = false;
end

fprintf('\n');

%% Extract Key Information
fprintf('[2/7] Extracting key information...\n');

% Get dose grid dimensions
if hasReconstructed
    doseGridSize = size(totalDose);
elseif ~isempty(fieldDoses{1})
    doseGridSize = size(fieldDoses{1}.physicalDose);
else
    error('No valid dose data found');
end

fprintf('  - Dose grid dimensions: %d x %d x %d\n', doseGridSize(1), doseGridSize(2), doseGridSize(3));

% Count valid fields
numValidFields = sum(~cellfun(@isempty, fieldDoses));
fprintf('  - Valid fields calculated: %d / %d\n', numValidFields, length(fieldDoses));

% Extract beam information
beamInfo = struct('beamIdx', {}, 'gantryAngle', {}, 'couchAngle', {}, 'maxDose', {});
for i = 1:length(fieldDoses)
    if ~isempty(fieldDoses{i})
        beamInfo(i).beamIdx = i;
        beamInfo(i).gantryAngle = fieldDoses{i}.gantryAngle;
        beamInfo(i).couchAngle = fieldDoses{i}.couchAngle;
        beamInfo(i).maxDose = max(fieldDoses{i}.physicalDose(:));
        fprintf('    Beam %d: Gantry=%.1f°, Couch=%.1f°, Max=%.2f Gy\n', ...
            i, beamInfo(i).gantryAngle, beamInfo(i).couchAngle, beamInfo(i).maxDose);
    end
end

fprintf('\n');

%% Plot 1: Individual Field Dose Distributions (Central Slices)
fprintf('[3/7] Creating individual field dose plots...\n');

% Find central slices
centralSlice.axial = round(doseGridSize(3) / 2);
centralSlice.sagittal = round(doseGridSize(1) / 2);
centralSlice.coronal = round(doseGridSize(2) / 2);

% Determine subplot layout
numFields = length(fieldDoses);
nRows = ceil(sqrt(numFields));
nCols = ceil(numFields / nRows);

% Axial view
figure('Name', 'Individual Field Doses - Axial View', 'Position', [100, 100, 1400, 900]);
for i = 1:numFields
    if ~isempty(fieldDoses{i})
        subplot(nRows, nCols, i);
        doseSlice = squeeze(fieldDoses{i}.physicalDose(:, :, centralSlice.axial));
        imagesc(doseSlice');
        axis equal tight;
        colorbar;
        colormap(jet);
        title(sprintf('Field %d (G:%.0f°, C:%.0f°)', ...
            i, fieldDoses{i}.gantryAngle, fieldDoses{i}.couchAngle));
        xlabel('X (pixels)');
        ylabel('Y (pixels)');
        caxis([0, max(doseSlice(:))]);
    end
end
sgtitle(sprintf('Individual Field Doses - Axial Slice %d/%d', centralSlice.axial, doseGridSize(3)), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('  - Created: Individual field doses (axial view)\n');

%% Plot 2: Field Dose Comparison (Side-by-Side Multi-View)
fprintf('[4/7] Creating multi-view field comparisons...\n');

% Select up to 4 fields for detailed comparison
fieldsToCompare = min(4, numFields);

figure('Name', 'Field Dose Multi-View Comparison', 'Position', [150, 150, 1600, 1000]);
for i = 1:fieldsToCompare
    if ~isempty(fieldDoses{i})
        % Axial view
        subplot(fieldsToCompare, 3, (i-1)*3 + 1);
        axialSlice = squeeze(fieldDoses{i}.physicalDose(:, :, centralSlice.axial));
        imagesc(axialSlice');
        axis equal tight;
        colorbar;
        colormap(jet);
        if i == 1
            title('Axial View');
        end
        ylabel(sprintf('Field %d\n(G:%.0f°)', i, fieldDoses{i}.gantryAngle));
        
        % Sagittal view
        subplot(fieldsToCompare, 3, (i-1)*3 + 2);
        sagittalSlice = squeeze(fieldDoses{i}.physicalDose(centralSlice.sagittal, :, :));
        imagesc(sagittalSlice');
        axis equal tight;
        colorbar;
        colormap(jet);
        if i == 1
            title('Sagittal View');
        end
        
        % Coronal view
        subplot(fieldsToCompare, 3, (i-1)*3 + 3);
        coronalSlice = squeeze(fieldDoses{i}.physicalDose(:, centralSlice.coronal, :));
        imagesc(coronalSlice');
        axis equal tight;
        colorbar;
        colormap(jet);
        if i == 1
            title('Coronal View');
        end
    end
end
sgtitle('Field Dose Multi-View Comparison', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('  - Created: Multi-view field comparison\n');

%% Plot 3: Reconstructed Total Dose
if hasReconstructed
    fprintf('[5/7] Creating reconstructed total dose plots...\n');
    
    % Find slice with maximum dose
    [~, maxIdx] = max(totalDose(:));
    [maxX, maxY, maxZ] = ind2sub(size(totalDose), maxIdx);
    
    figure('Name', 'Reconstructed Total Dose', 'Position', [200, 200, 1400, 500]);
    
    % Axial slice at max dose
    subplot(1, 3, 1);
    axialDose = squeeze(totalDose(:, :, maxZ));
    imagesc(axialDose');
    axis equal tight;
    colorbar;
    colormap(jet);
    title(sprintf('Axial (Z=%d, Max Dose)', maxZ));
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    hold on;
    plot(maxX, maxY, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
    
    % Sagittal slice at max dose
    subplot(1, 3, 2);
    sagittalDose = squeeze(totalDose(maxX, :, :));
    imagesc(sagittalDose');
    axis equal tight;
    colorbar;
    colormap(jet);
    title(sprintf('Sagittal (X=%d, Max Dose)', maxX));
    xlabel('Y (pixels)');
    ylabel('Z (pixels)');
    hold on;
    plot(maxY, maxZ, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
    
    % Coronal slice at max dose
    subplot(1, 3, 3);
    coronalDose = squeeze(totalDose(:, maxY, :));
    imagesc(coronalDose');
    axis equal tight;
    colorbar;
    colormap(jet);
    title(sprintf('Coronal (Y=%d, Max Dose)', maxY));
    xlabel('X (pixels)');
    ylabel('Z (pixels)');
    hold on;
    plot(maxX, maxZ, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
    
    sgtitle(sprintf('Reconstructed Total Dose (Max: %.2f Gy at [%d,%d,%d])', ...
        max(totalDose(:)), maxX, maxY, maxZ), 'FontSize', 14, 'FontWeight', 'bold');
    
    fprintf('  - Created: Reconstructed total dose visualization\n');
    fprintf('    Max dose: %.2f Gy at position [%d, %d, %d]\n', ...
        max(totalDose(:)), maxX, maxY, maxZ);
else
    fprintf('[5/7] Skipping reconstructed dose plots (not available)\n');
end

%% Plot 4: Dose Comparison (if available)
if hasComparison
    fprintf('[6/7] Creating dose comparison plots...\n');
    
    % Find slice with maximum difference
    [~, maxDiffIdx] = max(abs(comparison.difference(:)));
    [mdX, mdY, mdZ] = ind2sub(size(comparison.difference), maxDiffIdx);
    
    figure('Name', 'Dose Comparison: Calculated vs Reference', 'Position', [250, 250, 1600, 1000]);
    
    % Calculated dose
    subplot(2, 3, 1);
    calcSlice = squeeze(comparison.calculated(:, :, mdZ));
    imagesc(calcSlice');
    axis equal tight;
    colorbar;
    colormap(gca, jet);
    title('Calculated Dose (Axial)');
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    
    % Reference dose
    subplot(2, 3, 2);
    refSlice = squeeze(comparison.reference(:, :, mdZ));
    imagesc(refSlice');
    axis equal tight;
    colorbar;
    colormap(gca, jet);
    title('Reference Dose (Axial)');
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    
    % Difference
    subplot(2, 3, 3);
    diffSlice = squeeze(comparison.difference(:, :, mdZ));
    imagesc(diffSlice');
    axis equal tight;
    colorbar;
    colormap(gca, 'parula');
    title('Difference (Calc - Ref)');
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    
    % Profile comparison (through max difference point)
    subplot(2, 3, 4);
    calcProfile = squeeze(comparison.calculated(mdX, :, mdZ));
    refProfile = squeeze(comparison.reference(mdX, :, mdZ));
    plot(1:length(calcProfile), calcProfile, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
    hold on;
    plot(1:length(refProfile), refProfile, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
    grid on;
    xlabel('Y Position (pixels)');
    ylabel('Dose (Gy)');
    title(sprintf('Y-Profile at X=%d, Z=%d', mdX, mdZ));
    legend('Location', 'best');
    
    % Difference histogram
    subplot(2, 3, 5);
    histogram(comparison.difference(:), 50, 'EdgeColor', 'none');
    xlabel('Dose Difference (Gy)');
    ylabel('Frequency');
    title('Dose Difference Distribution');
    grid on;
    
    % Statistics table
    subplot(2, 3, 6);
    axis off;
    stats = {
        sprintf('Calculated Max: %.3f Gy', max(comparison.calculated(:)));
        sprintf('Reference Max: %.3f Gy', max(comparison.reference(:)));
        sprintf('Mean Absolute Diff: %.3f Gy', mean(abs(comparison.difference(:))));
        sprintf('Max Absolute Diff: %.3f Gy', max(abs(comparison.difference(:))));
        sprintf('RMS Difference: %.3f Gy', sqrt(mean(comparison.difference(:).^2)));
        sprintf('Mean Relative Diff: %.2f%%', 100*mean(abs(comparison.difference(:)))/max(comparison.reference(:)));
        sprintf('Max at: [%d, %d, %d]', mdX, mdY, mdZ);
    };
    text(0.1, 0.9, 'Comparison Statistics:', 'FontSize', 12, 'FontWeight', 'bold');
    for i = 1:length(stats)
        text(0.1, 0.85 - i*0.1, stats{i}, 'FontSize', 10);
    end
    
    sgtitle(sprintf('Dose Comparison - Slice Z=%d (Max Difference)', mdZ), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    fprintf('  - Created: Dose comparison visualization\n');
    fprintf('    Mean absolute difference: %.3f Gy\n', mean(abs(comparison.difference(:))));
    fprintf('    Max absolute difference: %.3f Gy\n', max(abs(comparison.difference(:))));
else
    fprintf('[6/7] Skipping dose comparison plots (not available)\n');
end

%% Plot 5: Field Contribution Analysis
fprintf('[7/7] Creating field contribution analysis...\n');

figure('Name', 'Field Contribution Analysis', 'Position', [300, 300, 1200, 800]);

% Bar chart of max doses per field
subplot(2, 2, 1);
maxDoses = cellfun(@(x) max(x.physicalDose(:)), fieldDoses(~cellfun(@isempty, fieldDoses)));
gantryAngles = [beamInfo.gantryAngle];
bar(1:length(maxDoses), maxDoses);
xlabel('Field Number');
ylabel('Max Dose (Gy)');
title('Maximum Dose per Field');
grid on;
set(gca, 'XTick', 1:length(maxDoses));

% Gantry angle distribution
subplot(2, 2, 2);
bar(1:length(gantryAngles), gantryAngles);
xlabel('Field Number');
ylabel('Gantry Angle (degrees)');
title('Gantry Angles');
grid on;
ylim([0, 360]);
set(gca, 'XTick', 1:length(gantryAngles));

% Dose contribution stacked (if reconstructed dose available)
if hasReconstructed
    subplot(2, 2, 3);
    % Sample along central axis
    centralProfile_z = round(doseGridSize(1)/2);
    centralProfile_y = round(doseGridSize(2)/2);
    
    profileData = zeros(doseGridSize(3), numValidFields);
    for i = 1:length(fieldDoses)
        if ~isempty(fieldDoses{i})
            profile = squeeze(fieldDoses{i}.physicalDose(centralProfile_z, centralProfile_y, :));
            profileData(:, i) = profile;
        end
    end
    
    area(1:doseGridSize(3), profileData);
    xlabel('Z Position (slice)');
    ylabel('Dose (Gy)');
    title('Stacked Field Contributions (Central Axis)');
    legend(arrayfun(@(x) sprintf('F%d', x), 1:numValidFields, 'UniformOutput', false), ...
        'Location', 'eastoutside');
    grid on;
    
    % Total dose profile
    subplot(2, 2, 4);
    totalProfile = squeeze(totalDose(centralProfile_z, centralProfile_y, :));
    plot(1:length(totalProfile), totalProfile, 'b-', 'LineWidth', 2);
    hold on;
    summedProfile = sum(profileData, 2);
    plot(1:length(summedProfile), summedProfile, 'r--', 'LineWidth', 2);
    xlabel('Z Position (slice)');
    ylabel('Dose (Gy)');
    title('Total Dose Profile Verification');
    legend('Reconstructed Total', 'Sum of Fields', 'Location', 'best');
    grid on;
else
    % Alternative plot: dose vs gantry angle
    subplot(2, 2, 3);
    scatter(gantryAngles, maxDoses, 100, 'filled');
    xlabel('Gantry Angle (degrees)');
    ylabel('Max Dose (Gy)');
    title('Max Dose vs Gantry Angle');
    grid on;
    
    % Field statistics
    subplot(2, 2, 4);
    axis off;
    text(0.1, 0.9, 'Field Statistics:', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('Number of fields: %d', numValidFields), 'FontSize', 10);
    text(0.1, 0.65, sprintf('Mean max dose: %.2f Gy', mean(maxDoses)), 'FontSize', 10);
    text(0.1, 0.55, sprintf('Std max dose: %.2f Gy', std(maxDoses)), 'FontSize', 10);
    text(0.1, 0.45, sprintf('Total max dose: %.2f Gy', sum(maxDoses)), 'FontSize', 10);
end

sgtitle('Field Contribution Analysis', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('  - Created: Field contribution analysis\n');

%% Summary Report
fprintf('\n=================================================\n');
fprintf('VISUALIZATION COMPLETE\n');
fprintf('=================================================\n');
fprintf('Summary:\n');
fprintf('  - Patient: %s, Session: %s\n', patientID, sessionName);
fprintf('  - Number of fields: %d\n', numValidFields);
fprintf('  - Dose grid: %d x %d x %d\n', doseGridSize(1), doseGridSize(2), doseGridSize(3));

if hasReconstructed
    fprintf('  - Reconstructed max dose: %.2f Gy\n', max(totalDose(:)));
end

if hasComparison
    fprintf('  - Mean dose difference: %.3f Gy\n', mean(abs(comparison.difference(:))));
    fprintf('  - Max dose difference: %.3f Gy\n', max(abs(comparison.difference(:))));
    fprintf('  - Relative difference: %.2f%%\n', ...
        100*mean(abs(comparison.difference(:)))/max(comparison.reference(:)));
end

fprintf('\nAll plots created successfully!\n');
fprintf('Total figures generated: %d\n', length(findall(0, 'Type', 'figure')));
fprintf('=================================================\n');
