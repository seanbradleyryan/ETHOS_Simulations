%% ETHOS IMRT Field-by-Field Dose Comprehensive Visualizer
% Purpose: Generate all 6 comprehensive figures from field dose calculations
% Requires: Output files from ETHOS_IMRT_Field_Calculator.m
% Author: Generated for ETHOS dose visualization
% Date: December 2025
%
% This script generates:
%   Figure 1: Field Dose Summary (statistics, geometry, contributions)
%   Figure 2: Individual Field Dose Distributions (spatial dose per beam)
%   Figure 3: Reconstructed Total Dose (3D analysis)
%   Figure 4: Dose Comparison (calculated vs reference)
%   Figure 5: Cumulative Dose Build-up (progressive accumulation)
%   Figure 6: DVH Analysis (dose-volume histograms)

clear; clc; close all;

%% Configuration
% Specify patient and session to visualize
patientID = '1194203';
sessionName = 'Session_1';

% Base directory (match your calculation script)
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
dataPath = fullfile(wd, 'FieldDoses', patientID, sessionName);

% Figure output directory
figPath = fullfile(dataPath, 'Figures');
if ~exist(figPath, 'dir')
    mkdir(figPath);
end

fprintf('========================================\n');
fprintf('ETHOS Field Dose Visualization\n');
fprintf('Patient: %s, Session: %s\n', patientID, sessionName);
fprintf('========================================\n\n');

%% Load Data
fprintf('[1/7] Loading data...\n');

% Load main field dose file
if ~exist(fullfile(dataPath, 'fieldDoses.mat'), 'file')
    error('Field dose file not found: %s', fullfile(dataPath, 'fieldDoses.mat'));
end

load(fullfile(dataPath, 'fieldDoses.mat'), 'fieldDoses', 'stf', 'pln', 'ct', 'cst');
fprintf('  - Loaded fieldDoses.mat\n');
fprintf('    Found %d fields\n', length(fieldDoses));

% Count successful fields
validFields = ~cellfun(@isempty, fieldDoses);
numValidFields = sum(validFields);
fprintf('    Valid fields: %d/%d\n', numValidFields, length(fieldDoses));

% Load reconstructed dose if available
hasReconstructed = false;
if exist(fullfile(dataPath, 'reconstructedDose.mat'), 'file')
    load(fullfile(dataPath, 'reconstructedDose.mat'), 'totalDose', 'calculatedGridSize');
    hasReconstructed = true;
    fprintf('  - Loaded reconstructedDose.mat\n');
    fprintf('    Grid size: %d x %d x %d\n', calculatedGridSize);
end

% Load comparison if available
hasComparison = false;
if exist(fullfile(dataPath, 'doseComparison.mat'), 'file')
    load(fullfile(dataPath, 'doseComparison.mat'), 'comparison');
    hasComparison = true;
    fprintf('  - Loaded doseComparison.mat\n');
end

fprintf('  ✓ Data loading complete\n\n');

%% Extract Field Information
fprintf('[2/7] Extracting field parameters...\n');

% Initialize arrays
numFields = length(fieldDoses);
maxDoses = zeros(numFields, 1);
nonzeroVoxels = zeros(numFields, 1);
gantryAngles = zeros(numFields, 1);
couchAngles = zeros(numFields, 1);
totalWeights = zeros(numFields, 1);

for i = 1:numFields
    if ~isempty(fieldDoses{i})
        maxDoses(i) = fieldDoses{i}.maxDose;
        nonzeroVoxels(i) = nnz(fieldDoses{i}.physicalDose);
        gantryAngles(i) = fieldDoses{i}.gantryAngle;
        couchAngles(i) = fieldDoses{i}.couchAngle;
        totalWeights(i) = sum(fieldDoses{i}.weights);
    end
end

% Filter valid fields
validIdx = validFields;
maxDoses_valid = maxDoses(validIdx);
nonzeroVoxels_valid = nonzeroVoxels(validIdx);
gantryAngles_valid = gantryAngles(validIdx);
totalWeights_valid = totalWeights(validIdx);

fprintf('  - Extracted parameters for %d fields\n', numValidFields);
fprintf('  - Max dose range: %.2f - %.2f Gy\n', min(maxDoses_valid), max(maxDoses_valid));
fprintf('  ✓ Parameter extraction complete\n\n');

%% FIGURE 1: Field Dose Summary
fprintf('[3/7] Creating Figure 1: Field Dose Summary...\n');

fig1 = figure('Position', [100, 100, 1600, 900], 'Name', 'Field Dose Summary');

% 1. Bar chart of maximum dose per field
subplot(2, 3, 1);
bar(find(validIdx), maxDoses_valid, 'FaceColor', [0.2 0.4 0.8]);
xlabel('Field Number');
ylabel('Maximum Dose (Gy)');
title('Maximum Dose Per Field');
grid on;
box on;

% 2. Polar plot of beam arrangement
subplot(2, 3, 2);
ax = polaraxes;
theta = deg2rad(gantryAngles_valid);
r = ones(size(theta));
polarplot(ax, theta, r, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.8 0.2 0.2], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
hold(ax, 'on');
% Add radial lines to origin
for i = 1:length(theta)
    polarplot(ax, [0 theta(i)], [0 1], 'k--', 'LineWidth', 1);
end
title(ax, 'Beam Arrangement (Gantry Angles)');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
rlim(ax, [0 1.2]);

% 3. Non-zero voxel coverage per field
subplot(2, 3, 3);
bar(find(validIdx), nonzeroVoxels_valid/1000, 'FaceColor', [0.2 0.8 0.4]);
xlabel('Field Number');
ylabel('Non-zero Voxels (×10³)');
title('Dose Coverage Per Field');
grid on;
box on;

% 4. Field weights comparison
subplot(2, 3, 4);
bar(find(validIdx), totalWeights_valid, 'FaceColor', [0.8 0.6 0.2]);
xlabel('Field Number');
ylabel('Total Weight (MU)');
title('Field Weights Comparison');
grid on;
box on;

% 5. Pie chart of relative dose contributions
subplot(2, 3, 5);
relativeContributions = maxDoses_valid / sum(maxDoses_valid) * 100;
validFieldNums = find(validIdx);
labels = cell(numValidFields, 1);
for i = 1:numValidFields
    labels{i} = sprintf('F%d (%.0f°)', validFieldNums(i), gantryAngles_valid(i));
end
pie(relativeContributions, labels);
title('Relative Dose Contributions');
colormap(jet(numValidFields));

% 6. Summary table
subplot(2, 3, 6);
axis off;
summaryText = {
    'FIELD SUMMARY',
    '═════════════════════════',
    sprintf('Total Fields: %d', numFields),
    sprintf('Valid Fields: %d', numValidFields),
    sprintf('Mean Max Dose: %.2f Gy', mean(maxDoses_valid)),
    sprintf('Total Max Dose: %.2f Gy', max(maxDoses_valid)),
    sprintf('Mean Coverage: %.1fk voxels', mean(nonzeroVoxels_valid)/1000),
    sprintf('Total MU: %.1f', sum(totalWeights_valid)),
    '',
    'BEAM GEOMETRY',
    '═════════════════════════',
    sprintf('Gantry Range: %.0f° - %.0f°', min(gantryAngles_valid), max(gantryAngles_valid)),
    sprintf('Mean Gantry: %.1f°', mean(gantryAngles_valid)),
    sprintf('Angular Spread: %.1f°', range(gantryAngles_valid)),
};
text(0.1, 0.9, summaryText, 'FontSize', 10, 'FontName', 'Courier', ...
    'VerticalAlignment', 'top', 'FontWeight', 'bold');

sgtitle(sprintf('Field Dose Summary - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(fig1, fullfile(figPath, 'Figure1_FieldSummary.png'));
saveas(fig1, fullfile(figPath, 'Figure1_FieldSummary.fig'));
fprintf('  ✓ Figure 1 saved\n\n');

%% FIGURE 2: Individual Field Dose Distributions
fprintf('[4/7] Creating Figure 2: Individual Field Dose Distributions...\n');

% Display first 4-6 fields
numFieldsToShow = min(6, numValidFields);
validFieldIndices = find(validIdx);
fieldsToPlot = validFieldIndices(1:numFieldsToShow);

fig2 = figure('Position', [100, 100, 1600, 900], 'Name', 'Individual Field Doses');

for i = 1:numFieldsToShow
    fieldIdx = fieldsToPlot(i);
    fieldDose = fieldDoses{fieldIdx}.physicalDose;
    
    % Get central axial slice
    centralSlice = round(size(fieldDose, 3) / 2);
    doseSlice = squeeze(fieldDose(:, :, centralSlice));
    
    subplot(2, 3, i);
    imagesc(doseSlice');
    axis equal tight;
    colorbar;
    colormap(jet);
    caxis([0 max(fieldDose(:))]);
    title(sprintf('Field %d - Gantry %.0f°\nMax: %.2f Gy', ...
        fieldIdx, fieldDoses{fieldIdx}.gantryAngle, fieldDoses{fieldIdx}.maxDose));
    xlabel('X (voxels)');
    ylabel('Y (voxels)');
end

sgtitle(sprintf('Individual Field Dose Distributions (Central Axial Slice) - Patient %s', patientID), ...
    'FontSize', 14, 'FontWeight', 'bold');

saveas(fig2, fullfile(figPath, 'Figure2_IndividualFields.png'));
saveas(fig2, fullfile(figPath, 'Figure2_IndividualFields.fig'));
fprintf('  ✓ Figure 2 saved\n\n');

%% FIGURE 3: Reconstructed Total Dose
fprintf('[5/7] Creating Figure 3: Reconstructed Total Dose...\n');

if hasReconstructed && max(totalDose(:)) > 0
    fig3 = figure('Position', [100, 100, 1800, 1000], 'Name', 'Reconstructed Total Dose');
    
    % Get dimensions
    [nx, ny, nz] = size(totalDose);
    cx = round(nx/2);
    cy = round(ny/2);
    cz = round(nz/2);
    maxDoseTotal = max(totalDose(:));
    
    % 1. Axial view
    subplot(2, 4, 1);
    imagesc(squeeze(totalDose(:, :, cz))');
    axis equal tight;
    colorbar;
    colormap(jet);
    title(sprintf('Axial View (z=%d)\nMax: %.2f Gy', cz, maxDoseTotal));
    xlabel('X'); ylabel('Y');
    
    % 2. Sagittal view
    subplot(2, 4, 2);
    imagesc(squeeze(totalDose(cx, :, :))');
    axis equal tight;
    colorbar;
    colormap(jet);
    title(sprintf('Sagittal View (x=%d)', cx));
    xlabel('Y'); ylabel('Z');
    
    % 3. Coronal view
    subplot(2, 4, 3);
    imagesc(squeeze(totalDose(:, cy, :))');
    axis equal tight;
    colorbar;
    colormap(jet);
    title(sprintf('Coronal View (y=%d)', cy));
    xlabel('X'); ylabel('Z');
    
    % 4. 3D isodose surfaces
    subplot(2, 4, 4);
    hold on;
    isodoseLevels = [0.95, 0.80, 0.50, 0.20];
    colors = [1 0 0; 1 0.5 0; 1 1 0; 0 1 0];
    legendEntries = {};
    
    for i = 1:length(isodoseLevels)
        level = isodoseLevels(i) * maxDoseTotal;
        try
            fv = isosurface(totalDose, level);
            if ~isempty(fv.vertices) && size(fv.vertices, 1) > 0
                p = patch(fv, 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
                legendEntries{end+1} = sprintf('%d%%', round(isodoseLevels(i)*100));
            end
        catch ME
            fprintf('    Warning: Could not create isosurface at %d%%: %s\n', ...
                round(isodoseLevels(i)*100), ME.message);
        end
    end
    
    view(3);
    axis vis3d equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Isodose Surfaces');
    if ~isempty(legendEntries)
        legend(legendEntries, 'Location', 'best');
    end
    grid on;
    lighting gouraud;
    camlight;
    
    % 5-7. Dose profiles through isocenter
    % X-axis profile
    subplot(2, 4, 5);
    profile_x = squeeze(totalDose(:, cy, cz));
    plot(1:nx, profile_x, 'b-', 'LineWidth', 2);
    xlabel('X (voxels)');
    ylabel('Dose (Gy)');
    title('X-axis Profile');
    grid on;
    ylim([0 maxDoseTotal*1.1]);
    
    % Y-axis profile
    subplot(2, 4, 6);
    profile_y = squeeze(totalDose(cx, :, cz));
    plot(1:ny, profile_y, 'r-', 'LineWidth', 2);
    xlabel('Y (voxels)');
    ylabel('Dose (Gy)');
    title('Y-axis Profile');
    grid on;
    ylim([0 maxDoseTotal*1.1]);
    
    % Z-axis profile
    subplot(2, 4, 7);
    profile_z = squeeze(totalDose(cx, cy, :));
    plot(1:nz, profile_z, 'g-', 'LineWidth', 2);
    xlabel('Z (voxels)');
    ylabel('Dose (Gy)');
    title('Z-axis Profile');
    grid on;
    ylim([0 maxDoseTotal*1.1]);
    
    % 8. Statistical summary
    subplot(2, 4, 8);
    axis off;
    
    % Calculate dose statistics
    doseAbove95 = sum(totalDose(:) >= 0.95*maxDoseTotal);
    doseAbove50 = sum(totalDose(:) >= 0.50*maxDoseTotal);
    doseAbove20 = sum(totalDose(:) >= 0.20*maxDoseTotal);
    totalVoxels = numel(totalDose);
    
    statsText = {
        'DOSE STATISTICS',
        '═════════════════════════',
        sprintf('Max Dose: %.2f Gy', maxDoseTotal),
        sprintf('Mean Dose: %.2f Gy', mean(totalDose(:))),
        sprintf('Median: %.2f Gy', median(totalDose(:))),
        sprintf('Std Dev: %.2f Gy', std(totalDose(:))),
        '',
        'COVERAGE VOLUMES',
        '═════════════════════════',
        sprintf('V95%%: %.1f%% of grid', doseAbove95/totalVoxels*100),
        sprintf('V50%%: %.1f%% of grid', doseAbove50/totalVoxels*100),
        sprintf('V20%%: %.1f%% of grid', doseAbove20/totalVoxels*100),
        '',
        'GRID INFO',
        '═════════════════════════',
        sprintf('Dimensions: %d×%d×%d', nx, ny, nz),
        sprintf('Total voxels: %d', totalVoxels),
    };
    text(0.1, 0.9, statsText, 'FontSize', 9, 'FontName', 'Courier', ...
        'VerticalAlignment', 'top', 'FontWeight', 'bold');
    
    sgtitle(sprintf('Reconstructed Total Dose Distribution - Patient %s', patientID), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    saveas(fig3, fullfile(figPath, 'Figure3_ReconstructedDose.png'));
    saveas(fig3, fullfile(figPath, 'Figure3_ReconstructedDose.fig'));
    fprintf('  ✓ Figure 3 saved\n\n');
else
    fprintf('  ⚠ Skipping Figure 3: No reconstructed dose available\n\n');
end

%% FIGURE 4: Dose Comparison
fprintf('[6/7] Creating Figure 4: Dose Comparison...\n');

if hasComparison && isfield(comparison, 'reference') && isfield(comparison, 'calculated')
    fig4 = figure('Position', [100, 100, 1800, 900], 'Name', 'Dose Comparison');
    
    calcDose = comparison.calculated;
    refDose = comparison.reference;
    
    % Check if we can do direct comparison
    if isequal(size(calcDose), size(refDose))
        diffDose = calcDose - refDose;
        
        % Get central slice
        cz = round(size(calcDose, 3) / 2);
        
        % 1. Calculated dose
        subplot(2, 3, 1);
        imagesc(squeeze(calcDose(:, :, cz))');
        axis equal tight;
        colorbar;
        colormap(jet);
        title(sprintf('Calculated Dose\nMax: %.2f Gy', max(calcDose(:))));
        xlabel('X'); ylabel('Y');
        
        % 2. Reference dose
        subplot(2, 3, 2);
        imagesc(squeeze(refDose(:, :, cz))');
        axis equal tight;
        colorbar;
        colormap(jet);
        title(sprintf('Reference Dose\nMax: %.2f Gy', max(refDose(:))));
        xlabel('X'); ylabel('Y');
        
        % 3. Difference map
        subplot(2, 3, 3);
        imagesc(squeeze(diffDose(:, :, cz))');
        axis equal tight;
        colorbar;
        colormap(redblue(256));
        maxAbsDiff = max(abs(diffDose(:)));
        caxis([-maxAbsDiff maxAbsDiff]);
        title(sprintf('Difference Map\nMax Diff: ±%.2f Gy', maxAbsDiff));
        xlabel('X'); ylabel('Y');
        
        % 4. Central axis profile comparison
        subplot(2, 3, 4);
        cx = round(size(calcDose, 1) / 2);
        cy = round(size(calcDose, 2) / 2);
        profile_calc = squeeze(calcDose(cx, cy, :));
        profile_ref = squeeze(refDose(cx, cy, :));
        plot(1:length(profile_calc), profile_calc, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
        hold on;
        plot(1:length(profile_ref), profile_ref, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
        xlabel('Z (voxels)');
        ylabel('Dose (Gy)');
        title('Central Axis Dose Profile');
        legend('Location', 'best');
        grid on;
        
        % 5. Scatter plot
        subplot(2, 3, 5);
        % Sample points for scatter (use every 10th point to avoid overplotting)
        sampIdx = 1:10:numel(calcDose);
        scatter(refDose(sampIdx), calcDose(sampIdx), 10, 'filled', 'MarkerFaceAlpha', 0.3);
        hold on;
        maxVal = max([refDose(:); calcDose(:)]);
        plot([0 maxVal], [0 maxVal], 'r--', 'LineWidth', 2);
        xlabel('Reference Dose (Gy)');
        ylabel('Calculated Dose (Gy)');
        title('Dose Correlation');
        axis equal;
        grid on;
        
        % Calculate correlation
        R = corrcoef(refDose(:), calcDose(:));
        text(0.1*maxVal, 0.9*maxVal, sprintf('R² = %.4f', R(1,2)^2), ...
            'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
        
        % 6. Quantitative metrics
        subplot(2, 3, 6);
        axis off;
        
        % Calculate gamma pass rate (optional, if you want to add)
        metricsText = {
            'COMPARISON METRICS',
            '═════════════════════════',
            sprintf('Mean Abs Diff: %.3f Gy', mean(abs(diffDose(:)))),
            sprintf('Max Difference: %.3f Gy', max(abs(diffDose(:)))),
            sprintf('RMS Difference: %.3f Gy', sqrt(mean(diffDose(:).^2))),
            sprintf('Correlation R²: %.4f', R(1,2)^2),
            '',
            'RELATIVE DIFFERENCE',
            '═════════════════════════',
            sprintf('Mean Rel Diff: %.2f%%', mean(abs(diffDose(:))./refDose(:)*100, 'omitnan')),
            sprintf('Max Rel Diff: %.2f%%', max(abs(diffDose(:))./refDose(:)*100, [], 'omitnan')),
            '',
            'DOSE STATISTICS',
            '═════════════════════════',
            'Calculated:',
            sprintf('  Max: %.2f Gy', max(calcDose(:))),
            sprintf('  Mean: %.2f Gy', mean(calcDose(:))),
            '',
            'Reference:',
            sprintf('  Max: %.2f Gy', max(refDose(:))),
            sprintf('  Mean: %.2f Gy', mean(refDose(:))),
        };
        text(0.1, 0.9, metricsText, 'FontSize', 9, 'FontName', 'Courier', ...
            'VerticalAlignment', 'top', 'FontWeight', 'bold');
        
    else
        % Size mismatch - show what we can
        subplot(2, 2, 1);
        imagesc(squeeze(max(calcDose, [], 3))');
        axis equal tight;
        title(sprintf('Calculated (Max Projection)\n%d×%d×%d', size(calcDose)));
        colorbar;
        colormap(jet);
        
        subplot(2, 2, 2);
        imagesc(squeeze(max(refDose, [], 3))');
        axis equal tight;
        title(sprintf('Reference (Max Projection)\n%d×%d×%d', size(refDose)));
        colorbar;
        colormap(jet);
        
        subplot(2, 2, [3 4]);
        axis off;
        text(0.5, 0.5, {'⚠ Grid Size Mismatch', '', ...
            sprintf('Calculated: %d × %d × %d', size(calcDose)), ...
            sprintf('Reference: %d × %d × %d', size(refDose)), ...
            '', 'Direct voxel-by-voxel comparison not possible', ...
            'Showing maximum intensity projections instead'}, ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    sgtitle(sprintf('Dose Comparison: Calculated vs Reference - Patient %s', patientID), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    saveas(fig4, fullfile(figPath, 'Figure4_DoseComparison.png'));
    saveas(fig4, fullfile(figPath, 'Figure4_DoseComparison.fig'));
    fprintf('  ✓ Figure 4 saved\n\n');
else
    fprintf('  ⚠ Skipping Figure 4: No comparison data available\n\n');
end

%% FIGURE 5: Cumulative Dose Build-up
fprintf('[7/7] Creating Figure 5: Cumulative Dose Build-up...\n');

if numValidFields > 0
    fig5 = figure('Position', [100, 100, 1800, 1000], 'Name', 'Cumulative Dose Build-up');
    
    % Show up to 6 cumulative stages
    numStages = min(6, numValidFields);
    validFieldIndices = find(validIdx);
    
    % Initialize cumulative dose
    cumulativeDose = zeros(size(fieldDoses{validFieldIndices(1)}.physicalDose));
    
    for stage = 1:numStages
        fieldIdx = validFieldIndices(stage);
        
        % Add current field to cumulative
        cumulativeDose = cumulativeDose + fieldDoses{fieldIdx}.physicalDose;
        
        % Get central slice
        cz = round(size(cumulativeDose, 3) / 2);
        doseSlice = squeeze(cumulativeDose(:, :, cz));
        
        subplot(2, 3, stage);
        imagesc(doseSlice');
        axis equal tight;
        colorbar;
        colormap(jet);
        if hasReconstructed
            caxis([0 max(totalDose(:))]);
        else
            caxis([0 max(cumulativeDose(:))]);
        end
        
        % Add percentage of total
        if hasReconstructed
            percentTotal = max(cumulativeDose(:)) / max(totalDose(:)) * 100;
            title(sprintf('After %d Fields (%.0f°)\nMax: %.2f Gy (%.0f%% of total)', ...
                stage, fieldDoses{fieldIdx}.gantryAngle, max(cumulativeDose(:)), percentTotal));
        else
            title(sprintf('After %d Fields (%.0f°)\nMax: %.2f Gy', ...
                stage, fieldDoses{fieldIdx}.gantryAngle, max(cumulativeDose(:))));
        end
        xlabel('X'); ylabel('Y');
    end
    
    sgtitle(sprintf('Cumulative Dose Build-up (Central Axial Slice) - Patient %s', patientID), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    saveas(fig5, fullfile(figPath, 'Figure5_CumulativeBuildUp.png'));
    saveas(fig5, fullfile(figPath, 'Figure5_CumulativeBuildUp.fig'));
    fprintf('  ✓ Figure 5 saved\n\n');
else
    fprintf('  ⚠ Skipping Figure 5: No valid fields available\n\n');
end

%% FIGURE 6: DVH Analysis
fprintf('[8/7] Creating Figure 6: DVH Analysis...\n');

if exist('cst', 'var') && hasReconstructed && max(totalDose(:)) > 0
    fig6 = figure('Position', [100, 100, 1400, 800], 'Name', 'DVH Analysis');
    
    % Calculate DVHs for all structures
    numStructures = size(cst, 1);
    colors = lines(numStructures);
    
    hold on;
    legendEntries = {};
    
    for i = 1:numStructures
        structName = cst{i, 2};
        
        % Get structure indices - handle different CST formats
        structIndices = [];
        if size(cst, 2) >= 4 && ~isempty(cst{i, 4})
            if iscell(cst{i, 4})
                if ~isempty(cst{i, 4}{1})
                    structIndices = cst{i, 4}{1};
                end
            elseif isnumeric(cst{i, 4})
                structIndices = cst{i, 4};
            end
        end
        
        if isempty(structIndices)
            fprintf('    Skipping structure %d (%s): no indices found\n', i, structName);
            continue;
        end
        
        % Validate indices
        if max(structIndices) > numel(totalDose) || min(structIndices) < 1
            fprintf('    Skipping structure %d (%s): invalid indices\n', i, structName);
            continue;
        end
        
        % Extract doses for this structure
        structDoses = totalDose(structIndices);
        
        if isempty(structDoses) || all(structDoses == 0)
            fprintf('    Skipping structure %d (%s): no dose data\n', i, structName);
            continue;
        end
        
        % Calculate cumulative DVH
        [counts, edges] = histcounts(structDoses, 100);
        cumCounts = flip(cumsum(flip(counts)));
        dvh_percent = [cumCounts / length(structDoses) * 100, 0];
        dose_bins = edges;
        
        % Plot DVH
        plot(dose_bins, dvh_percent, 'LineWidth', 2.5, 'Color', colors(i,:));
        legendEntries{end+1} = structName;
        
        fprintf('    Plotted DVH for: %s (%d voxels, max dose: %.2f Gy)\n', ...
            structName, length(structDoses), max(structDoses));
    end
    
    xlabel('Dose (Gy)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Volume (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Dose-Volume Histograms - Patient %s', patientID), ...
        'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    box on;
    
    if ~isempty(legendEntries)
        legend(legendEntries, 'Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 9);
    end
    
    xlim([0 max(totalDose(:))*1.05]);
    ylim([0 105]);
    
    % Add grid lines
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    
    saveas(fig6, fullfile(figPath, 'Figure6_DVH_Analysis.png'));
    saveas(fig6, fullfile(figPath, 'Figure6_DVH_Analysis.fig'));
    fprintf('  ✓ Figure 6 saved\n\n');
else
    if ~exist('cst', 'var')
        fprintf('  ⚠ Skipping Figure 6: CST data not available\n\n');
    elseif ~hasReconstructed
        fprintf('  ⚠ Skipping Figure 6: Reconstructed dose not available\n\n');
    else
        fprintf('  ⚠ Skipping Figure 6: Reconstructed dose is empty\n\n');
    end
end

%% Summary Report
fprintf('========================================\n');
fprintf('VISUALIZATION COMPLETE!\n');
fprintf('========================================\n\n');

fprintf('Output Directory: %s\n\n', figPath);

fprintf('Generated Figures:\n');
fprintf('  1. Field Dose Summary ✓\n');
fprintf('  2. Individual Field Dose Distributions ✓\n');
if hasReconstructed
    fprintf('  3. Reconstructed Total Dose ✓\n');
else
    fprintf('  3. Reconstructed Total Dose ✗ (no data)\n');
end
if hasComparison
    fprintf('  4. Dose Comparison ✓\n');
else
    fprintf('  4. Dose Comparison ✗ (no data)\n');
end
if numValidFields > 0
    fprintf('  5. Cumulative Dose Build-up ✓\n');
else
    fprintf('  5. Cumulative Dose Build-up ✗ (no valid fields)\n');
end
if exist('cst', 'var') && hasReconstructed
    fprintf('  6. DVH Analysis ✓\n');
else
    fprintf('  6. DVH Analysis ✗ (no data)\n');
end

fprintf('\n');
fprintf('Files saved as:\n');
fprintf('  - PNG format (for reports/presentations)\n');
fprintf('  - FIG format (for further editing in MATLAB)\n');
fprintf('\n========================================\n');

%% Helper function for red-blue diverging colormap
function cmap = redblue(m)
    % Red-blue diverging colormap centered at white
    if nargin < 1
        m = 256;
    end
    
    % Create colormap with smooth transition through white
    n = fix(0.5*m);
    
    % Blue to white
    r1 = linspace(0, 1, n)';
    g1 = linspace(0, 1, n)';
    b1 = ones(n, 1);
    
    % White to red
    r2 = ones(m-n, 1);
    g2 = linspace(1, 0, m-n)';
    b2 = linspace(1, 0, m-n)';
    
    cmap = [r1 g1 b1; r2 g2 b2];
end