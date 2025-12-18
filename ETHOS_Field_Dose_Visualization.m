%% ETHOS Field-by-Field Dose Visualizer
% Purpose: Load and visualize individual field doses calculated by ETHOS dose calculator
% Author: Generated for ETHOS dose analysis verification
% Date: 2025

clear; clc; close all;

%% Configuration
% Specify patient and session to analyze
patientID = '1194203';
sessionName = 'Session_1';

% Base paths
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
dataPath = fullfile(wd, 'FieldDoses', patientID, sessionName);

% Verify data directory exists
if ~exist(dataPath, 'dir')
    error('Data directory does not exist: %s', dataPath);
end

fprintf('Loading data from: %s\n', dataPath);

%% Load Data
fprintf('\n=== Loading Calculation Results ===\n');

% Load field doses
if exist(fullfile(dataPath, 'fieldDoses.mat'), 'file')
    load(fullfile(dataPath, 'fieldDoses.mat'), 'fieldDoses', 'stf', 'pln', 'ct', 'cst');
    fprintf('✓ Loaded fieldDoses.mat\n');
    fprintf('  - Number of fields: %d\n', length(fieldDoses));
    numSuccessful = sum(~cellfun(@isempty, fieldDoses));
    fprintf('  - Successfully calculated: %d/%d\n', numSuccessful, length(fieldDoses));
else
    error('fieldDoses.mat not found in %s', dataPath);
end

% Load reconstructed total dose
if exist(fullfile(dataPath, 'reconstructedDose.mat'), 'file')
    load(fullfile(dataPath, 'reconstructedDose.mat'), 'totalDose', 'calculatedGridSize');
    fprintf('✓ Loaded reconstructedDose.mat\n');
    fprintf('  - Total dose grid: %d x %d x %d\n', size(totalDose));
    fprintf('  - Max total dose: %.2f Gy\n', max(totalDose(:)));
    hasReconstructedDose = true;
else
    fprintf('⚠ reconstructedDose.mat not found\n');
    hasReconstructedDose = false;
end

% Load comparison data if available
if exist(fullfile(dataPath, 'doseComparison.mat'), 'file')
    load(fullfile(dataPath, 'doseComparison.mat'), 'comparison');
    fprintf('✓ Loaded doseComparison.mat\n');
    hasComparison = true;
    if isfield(comparison, 'reference')
        fprintf('  - Reference dose max: %.2f Gy\n', max(comparison.reference(:)));
    end
else
    fprintf('⚠ doseComparison.mat not found\n');
    hasComparison = false;
end

%% Figure 1: Field Dose Summary
fprintf('\n=== Creating Figure 1: Field Dose Summary ===\n');
figure('Name', 'Field Dose Summary', 'Position', [100, 100, 1400, 800]);

% Extract field information
numFields = length(fieldDoses);
fieldInfo = struct('idx', {}, 'gantry', {}, 'couch', {}, 'maxDose', {}, 'nonzeroVoxels', {});

validIdx = 1;
for i = 1:numFields
    if ~isempty(fieldDoses{i})
        fieldInfo(validIdx).idx = fieldDoses{i}.beamIdx;
        fieldInfo(validIdx).gantry = fieldDoses{i}.gantryAngle;
        fieldInfo(validIdx).couch = fieldDoses{i}.couchAngle;
        fieldInfo(validIdx).maxDose = fieldDoses{i}.maxDose;
        fieldInfo(validIdx).nonzeroVoxels = nnz(fieldDoses{i}.physicalDose);
        fieldInfo(validIdx).totalWeight = sum(fieldDoses{i}.weights);
        validIdx = validIdx + 1;
    end
end

% Subplot 1: Max dose per field
subplot(2,3,1);
bar([fieldInfo.maxDose], 'FaceColor', [0.2 0.6 0.8]);
xlabel('Field Number');
ylabel('Max Dose (Gy)');
title('Maximum Dose per Field');
grid on;
set(gca, 'XTick', 1:length(fieldInfo));

% Subplot 2: Gantry angles
subplot(2,3,2);
polarplot(deg2rad([fieldInfo.gantry]), [fieldInfo.maxDose], 'o-', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
title('Beam Arrangement (Gantry Angles)');
thetaticks(0:45:315);

% Subplot 3: Number of non-zero voxels per field
subplot(2,3,3);
bar([fieldInfo.nonzeroVoxels], 'FaceColor', [0.8 0.4 0.2]);
xlabel('Field Number');
ylabel('Non-zero Voxels');
title('Dose Coverage per Field');
grid on;
set(gca, 'XTick', 1:length(fieldInfo));

% Subplot 4: Field weights
subplot(2,3,4);
bar([fieldInfo.totalWeight], 'FaceColor', [0.4 0.7 0.4]);
xlabel('Field Number');
ylabel('Total Weight (MU)');
title('Field Weights');
grid on;
set(gca, 'XTick', 1:length(fieldInfo));

% Subplot 5: Dose contribution (percentage)
subplot(2,3,5);
if hasReconstructedDose
    doseContribution = ([fieldInfo.maxDose] / max(totalDose(:))) * 100;
else
    doseContribution = ([fieldInfo.maxDose] / sum([fieldInfo.maxDose])) * 100;
end
pie(doseContribution);
title('Relative Dose Contribution');
legend(arrayfun(@(x) sprintf('Field %d', x), [fieldInfo.idx], 'UniformOutput', false), ...
    'Location', 'eastoutside');

% Subplot 6: Field parameters table
subplot(2,3,6);
axis off;
tableData = cell(length(fieldInfo), 4);
for i = 1:length(fieldInfo)
    tableData{i,1} = sprintf('%d', fieldInfo(i).idx);
    tableData{i,2} = sprintf('%.1f°', fieldInfo(i).gantry);
    tableData{i,3} = sprintf('%.2f', fieldInfo(i).maxDose);
    tableData{i,4} = sprintf('%.1f', fieldInfo(i).totalWeight);
end
uitable('Data', tableData, ...
    'ColumnName', {'Field', 'Gantry', 'Max Dose (Gy)', 'Weight'}, ...
    'Units', 'normalized', ...
    'Position', [0.68, 0.05, 0.28, 0.4]);

sgtitle(sprintf('Field Dose Analysis: Patient %s - %s', patientID, sessionName), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% Figure 2: Individual Field Dose Distributions
fprintf('\n=== Creating Figure 2: Individual Field Dose Slices ===\n');

% Select central slice for visualization
if hasReconstructedDose
    centralSlice = round(size(totalDose, 3) / 2);
else
    centralSlice = round(size(fieldDoses{find(~cellfun(@isempty, fieldDoses), 1)}.physicalDose, 3) / 2);
end

% Calculate number of subplots needed
numValidFields = length(fieldInfo);
nRows = ceil(sqrt(numValidFields));
nCols = ceil(numValidFields / nRows);

figure('Name', 'Individual Field Dose Distributions', 'Position', [100, 100, 1600, 1000]);

for i = 1:numValidFields
    subplot(nRows, nCols, i);
    
    fieldIdx = fieldInfo(i).idx;
    doseSlice = fieldDoses{fieldIdx}.physicalDose(:, :, centralSlice);
    
    imagesc(doseSlice);
    colormap(jet);
    colorbar;
    axis equal tight;
    title(sprintf('Field %d (Gantry %.0f°)\nMax: %.2f Gy', ...
        fieldIdx, fieldInfo(i).gantry, fieldInfo(i).maxDose));
    xlabel('X (voxels)');
    ylabel('Y (voxels)');
    caxis([0, max(doseSlice(:))]);
end

sgtitle(sprintf('Field Dose Distributions - Central Slice (z=%d)', centralSlice), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% Figure 3: Reconstructed Total Dose
if hasReconstructedDose
    fprintf('\n=== Creating Figure 3: Reconstructed Total Dose ===\n');
    figure('Name', 'Reconstructed Total Dose', 'Position', [100, 100, 1400, 800]);
    
    % Three orthogonal views
    centralX = round(size(totalDose, 1) / 2);
    centralY = round(size(totalDose, 2) / 2);
    centralZ = round(size(totalDose, 3) / 2);
    
    % Axial view
    subplot(2,3,1);
    imagesc(totalDose(:, :, centralZ)');
    colormap(jet);
    colorbar;
    axis equal tight;
    title(sprintf('Axial (z=%d)', centralZ));
    xlabel('X'); ylabel('Y');
    
    % Sagittal view
    subplot(2,3,2);
    imagesc(squeeze(totalDose(centralX, :, :))');
    colormap(jet);
    colorbar;
    axis equal tight;
    title(sprintf('Sagittal (x=%d)', centralX));
    xlabel('Y'); ylabel('Z');
    
    % Coronal view
    subplot(2,3,3);
    imagesc(squeeze(totalDose(:, centralY, :))');
    colormap(jet);
    colorbar;
    axis equal tight;
    title(sprintf('Coronal (y=%d)', centralY));
    xlabel('X'); ylabel('Z');
    
    % 3D isodose view
    subplot(2,3,4);
    isodoseLevels = [0.2, 0.5, 0.8, 0.95] * max(totalDose(:));
    colors = [0.2 0.8 0.2; 0.8 0.8 0.2; 0.8 0.5 0.2; 0.8 0.2 0.2];
    
    hold on;
    for i = 1:length(isodoseLevels)
        isosurface(totalDose, isodoseLevels(i), 'FaceColor', colors(i,:), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
    view(3);
    axis equal tight;
    grid on;
    title('3D Isodose Surfaces');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend(arrayfun(@(x) sprintf('%.0f%%', x*100), [0.2, 0.5, 0.8, 0.95], ...
        'UniformOutput', false), 'Location', 'best');
    camlight; lighting gouraud;
    
    % Dose profile through isocenter
    subplot(2,3,5);
    profileX = totalDose(centralX, centralY, :);
    profileY = totalDose(centralX, :, centralZ);
    profileZ = totalDose(:, centralY, centralZ);
    
    plot(squeeze(profileX), 'b-', 'LineWidth', 2); hold on;
    plot(squeeze(profileY), 'r-', 'LineWidth', 2);
    plot(squeeze(profileZ), 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Position (voxels)');
    ylabel('Dose (Gy)');
    title('Dose Profiles Through Isocenter');
    legend('Z-axis', 'Y-axis', 'X-axis', 'Location', 'best');
    
    % Dose statistics
    subplot(2,3,6);
    axis off;
    statsText = {
        sprintf('Grid Size: %d × %d × %d', size(totalDose));
        sprintf('Max Dose: %.2f Gy', max(totalDose(:)));
        sprintf('Mean Dose: %.2f Gy', mean(totalDose(totalDose>0)));
        sprintf('Min Dose: %.2f Gy', min(totalDose(totalDose>0)));
        sprintf('Volume > 0 Gy: %.1f%%', 100*nnz(totalDose)/numel(totalDose));
        sprintf('Volume > 50%% max: %.1f%%', 100*nnz(totalDose>0.5*max(totalDose(:)))/numel(totalDose));
        ' ';
        sprintf('Number of Fields: %d', numValidFields);
        sprintf('Total Weight: %.1f MU', sum([fieldInfo.totalWeight]));
    };
    text(0.1, 0.9, statsText, 'VerticalAlignment', 'top', 'FontSize', 11, ...
        'FontName', 'FixedWidth');
    
    sgtitle('Reconstructed Total Dose Distribution', 'FontSize', 14, 'FontWeight', 'bold');
end

%% Figure 4: Dose Comparison (if reference available)
if hasComparison && isfield(comparison, 'difference')
    fprintf('\n=== Creating Figure 4: Dose Comparison with Reference ===\n');
    figure('Name', 'Dose Comparison', 'Position', [100, 100, 1600, 900]);
    
    centralSlice = round(size(comparison.calculated, 3) / 2);
    
    % Calculated dose
    subplot(2,3,1);
    imagesc(comparison.calculated(:, :, centralSlice)');
    colormap(jet);
    colorbar;
    axis equal tight;
    title(sprintf('Calculated Dose\nMax: %.2f Gy', max(comparison.calculated(:))));
    xlabel('X'); ylabel('Y');
    
    % Reference dose
    subplot(2,3,2);
    imagesc(comparison.reference(:, :, centralSlice)');
    colormap(jet);
    colorbar;
    axis equal tight;
    title(sprintf('Reference Dose\nMax: %.2f Gy', max(comparison.reference(:))));
    xlabel('X'); ylabel('Y');
    
    % Difference map
    subplot(2,3,3);
    diffSlice = comparison.difference(:, :, centralSlice)';
    imagesc(diffSlice);
    colormap(jet);
    colorbar;
    axis equal tight;
    title(sprintf('Difference (Calc - Ref)\nMax: %.2f Gy', max(abs(comparison.difference(:)))));
    xlabel('X'); ylabel('Y');
    caxis([-max(abs(diffSlice(:))), max(abs(diffSlice(:)))]);
    
    % Dose profile comparison
    subplot(2,3,4);
    profileCalc = squeeze(comparison.calculated(centralX, centralY, :));
    profileRef = squeeze(comparison.reference(centralX, centralY, :));
    plot(profileCalc, 'b-', 'LineWidth', 2); hold on;
    plot(profileRef, 'r--', 'LineWidth', 2);
    grid on;
    xlabel('Z Position (voxels)');
    ylabel('Dose (Gy)');
    title('Central Axis Dose Profile');
    legend('Calculated', 'Reference', 'Location', 'best');
    
    % Scatter plot
    subplot(2,3,5);
    calcVec = comparison.calculated(:);
    refVec = comparison.reference(:);
    % Sample points for faster plotting
    sampleIdx = randsample(length(calcVec), min(10000, length(calcVec)));
    scatter(refVec(sampleIdx), calcVec(sampleIdx), 10, 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    maxDose = max([max(calcVec), max(refVec)]);
    plot([0 maxDose], [0 maxDose], 'r--', 'LineWidth', 2);
    grid on;
    xlabel('Reference Dose (Gy)');
    ylabel('Calculated Dose (Gy)');
    title('Dose Correlation');
    axis equal tight;
    
    % Statistics table
    subplot(2,3,6);
    axis off;
    if isfield(comparison, 'metrics')
        statsText = {
            'Comparison Metrics:';
            ' ';
            sprintf('Mean Abs Diff: %.3f Gy', comparison.metrics.meanAbsDiff);
            sprintf('Max Diff: %.3f Gy', comparison.metrics.maxDiff);
            sprintf('RMS Diff: %.3f Gy', comparison.metrics.rmsDiff);
            ' ';
            sprintf('Calc Max: %.2f Gy', max(comparison.calculated(:)));
            sprintf('Ref Max: %.2f Gy', max(comparison.reference(:)));
            sprintf('Relative Diff: %.1f%%', 100*abs(max(comparison.calculated(:))-max(comparison.reference(:)))/max(comparison.reference(:)));
        };
    else
        statsText = {
            'Grid Size Mismatch';
            ' ';
            sprintf('Calculated: %d × %d × %d', comparison.calculatedGrid);
            sprintf('Reference: %d × %d × %d', comparison.referenceGrid);
            ' ';
            'Direct comparison not possible';
            'due to different dose grid';
            'resolutions.'
        };
    end
    text(0.1, 0.9, statsText, 'VerticalAlignment', 'top', 'FontSize', 11, ...
        'FontName', 'FixedWidth');
    
    sgtitle('Dose Comparison: Calculated vs Reference', 'FontSize', 14, 'FontWeight', 'bold');
end

%% Figure 5: Cumulative Dose Build-up
if hasReconstructedDose && numValidFields > 1
    fprintf('\n=== Creating Figure 5: Cumulative Dose Build-up ===\n');
    figure('Name', 'Cumulative Dose Build-up', 'Position', [100, 100, 1400, 800]);
    
    % Calculate cumulative dose
    cumulativeDose = zeros(size(fieldDoses{fieldInfo(1).idx}.physicalDose));
    maxDoses = zeros(1, numValidFields);
    
    centralSlice = round(size(cumulativeDose, 3) / 2);
    
    for i = 1:min(6, numValidFields)  % Show up to 6 stages
        subplot(2, 3, i);
        
        % Add this field to cumulative
        cumulativeDose = cumulativeDose + fieldDoses{fieldInfo(i).idx}.physicalDose;
        maxDoses(i) = max(cumulativeDose(:));
        
        imagesc(cumulativeDose(:, :, centralSlice)');
        colormap(jet);
        colorbar;
        axis equal tight;
        caxis([0, max(totalDose(:))]);
        
        if i < numValidFields
            title(sprintf('After Field %d (%.0f°)\nMax: %.2f Gy', ...
                i, fieldInfo(i).gantry, maxDoses(i)));
        else
            title(sprintf('All %d Fields\nMax: %.2f Gy', ...
                numValidFields, maxDoses(i)));
        end
        xlabel('X'); ylabel('Y');
    end
    
    sgtitle('Cumulative Dose Build-up', 'FontSize', 14, 'FontWeight', 'bold');
end

%% Figure 6: DVH Analysis (if structures available)
if exist('cst', 'var') && ~isempty(cst) && hasReconstructedDose
    fprintf('\n=== Creating Figure 6: DVH Analysis ===\n');
    
    try
        figure('Name', 'Dose-Volume Histogram', 'Position', [100, 100, 1200, 700]);
        
        % Get structure names
        numStructures = size(cst, 1);
        colors = lines(numStructures);
        
        hold on;
        legendEntries = {};
        
        for i = 1:numStructures
            structName = cst{i, 2};
            
            % Skip if no indices
            if isempty(cst{i, 4}) || ~isfield(cst{i, 4}, 'allIndices')
                continue;
            end
            
            % Get voxel indices for this structure
            voxelIdx = cst{i, 4}.allIndices;
            
            % Extract doses for these voxels
            if ~isempty(voxelIdx) && all(voxelIdx <= numel(totalDose))
                structDose = totalDose(voxelIdx);
                
                % Calculate DVH
                [counts, edges] = histcounts(structDose, 100);
                cumCounts = cumsum(counts(end:-1:1));
                cumCounts = cumCounts(end:-1:1);
                dvh = 100 * cumCounts / sum(counts);
                doseBins = edges(1:end-1);
                
                % Plot
                plot(doseBins, dvh, 'LineWidth', 2, 'Color', colors(i,:));
                legendEntries{end+1} = structName;
            end
        end
        
        hold off;
        grid on;
        xlabel('Dose (Gy)');
        ylabel('Volume (%)');
        title('Dose-Volume Histogram');
        legend(legendEntries, 'Location', 'best', 'Interpreter', 'none');
        xlim([0, max(totalDose(:))*1.1]);
        ylim([0, 105]);
        
    catch ME
        fprintf('Warning: Could not create DVH: %s\n', ME.message);
    end
end

%% Summary Report
fprintf('\n');
fprintf('========================================\n');
fprintf('VISUALIZATION SUMMARY\n');
fprintf('========================================\n');
fprintf('Patient ID: %s\n', patientID);
fprintf('Session: %s\n', sessionName);
fprintf('Number of fields: %d\n', numValidFields);
fprintf('Max dose: %.2f Gy\n', max(totalDose(:)));
fprintf('Data path: %s\n', dataPath);
fprintf('========================================\n');
fprintf('\nAll visualizations created successfully!\n');