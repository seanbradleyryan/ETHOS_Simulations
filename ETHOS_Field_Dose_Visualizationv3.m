%% ETHOS Field Dose Visualization and Analysis
% Purpose: Comprehensive visualization of field-by-field dose calculations
%          including comparison with reference RTDOSE
%
% Visualizations included:
%   1. Individual field dose distributions
%   2. Total calculated vs reference dose comparison
%   3. Dose difference maps
%   4. Dose profiles (axial, sagittal, coronal)
%   5. Dose-volume histograms (DVH)
%   6. Gamma analysis (2%/2mm, 3%/3mm)
%   7. Field contribution analysis
%   8. Statistical summary
%
% Author: Generated for ETHOS dose analysis
% Date: 2025

clear; clc; close all;

%% ======================== CONFIGURATION ========================
% Patient and session configuration
patientID = '1194203';
sessionName = 'Session_1';

% Base directories
baseDir = '/mnt/weka/home/80030361/ETHOS_Simulations';
fieldDoseDir = fullfile(baseDir, 'FieldDoses', patientID, sessionName);
outputDir = fullfile(baseDir, 'Visualizations', patientID, sessionName);

% Visualization settings
saveFigures = true;          % Save figures to disk
figureFormat = 'png';        % 'png', 'pdf', 'fig', 'eps'
figureResolution = 300;      % DPI for raster formats

% Gamma analysis parameters
gammaParams.doseThreshold = [2, 3];      % Dose difference criteria (%)
gammaParams.distThreshold = [2, 3];      % Distance-to-agreement criteria (mm)
gammaParams.doseCutoff = 10;             % Exclude voxels below this % of max dose
gammaParams.searchRadius = 5;            % Search radius in mm

% Colormap settings
doseCmap = 'jet';
diffCmap = 'coolwarm';  % Will create custom if not available

%% ======================== CREATE OUTPUT DIRECTORY ========================
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% ======================== LOAD DATA ========================
fprintf('=======================================================\n');
fprintf('  ETHOS Field Dose Visualization\n');
fprintf('  Patient: %s, Session: %s\n', patientID, sessionName);
fprintf('=======================================================\n\n');

fprintf('[1/8] Loading data...\n');

% Try to load resampled data first
fieldDoseFile = fullfile(fieldDoseDir, 'fieldDoses.mat');
comparisonFile = fullfile(fieldDoseDir, 'doseComparison.mat');
ctResampledFile = fullfile(fieldDoseDir, 'ctResampled.mat');

if ~exist(fieldDoseFile, 'file')
    error('Field dose file not found: %s\nRun ethos_field_dosev3.m first.', fieldDoseFile);
end

% Load field doses
loadedData = load(fieldDoseFile);

% Determine which format we have
if isfield(loadedData, 'fieldDosesResampled')
    fieldDoses = loadedData.fieldDosesResampled;
    fprintf('  - Loaded resampled field doses (RTDOSE grid)\n');
    isResampled = true;
else
    fieldDoses = loadedData.fieldDoses;
    fprintf('  - Loaded original field doses (CT grid)\n');
    isResampled = false;
end

% Load total doses
if isfield(loadedData, 'totalDoseResampled')
    totalCalculated = loadedData.totalDoseResampled;
elseif isfield(loadedData, 'totalDose')
    totalCalculated = loadedData.totalDose;
else
    % Reconstruct from field doses
    totalCalculated = [];
end

% Load reference dose
if isfield(loadedData, 'referenceDose')
    referenceDose = loadedData.referenceDose;
else
    referenceDose = [];
end

% Load dose grid info
if isfield(loadedData, 'doseGrid')
    doseGrid = loadedData.doseGrid;
    resolution = doseGrid.resolution;
else
    resolution = [2.5, 2.5, 2.5];  % Default
end

% Load CT if available
if isfield(loadedData, 'ctResampled_struct')
    ct = loadedData.ctResampled_struct;
    ctCube = ct.cubeHU{1};
elseif exist(ctResampledFile, 'file')
    ctData = load(ctResampledFile);
    ct = ctData.ctResampled_struct;
    ctCube = ct.cubeHU{1};
elseif isfield(loadedData, 'ct')
    ct = loadedData.ct;
    ctCube = ct.cubeHU{1};
else
    ctCube = [];
end

% Load comparison data if available
if exist(comparisonFile, 'file')
    compData = load(comparisonFile);
    comparison = compData.comparison;
    if isfield(comparison, 'difference')
        doseDifference = comparison.difference;
    end
end

% Count valid fields
numFields = length(fieldDoses);
validFields = find(~cellfun(@isempty, fieldDoses));
numValidFields = length(validFields);

fprintf('  - Valid fields: %d/%d\n', numValidFields, numFields);
fprintf('  - Grid resolution: [%.2f, %.2f, %.2f] mm\n', resolution(1), resolution(2), resolution(3));

% Reconstruct total if needed
if isempty(totalCalculated) && numValidFields > 0
    fprintf('  - Reconstructing total dose from fields...\n');
    sampleSize = size(fieldDoses{validFields(1)}.physicalDose);
    totalCalculated = zeros(sampleSize);
    for idx = validFields'
        totalCalculated = totalCalculated + fieldDoses{idx}.physicalDose;
    end
end

gridSize = size(totalCalculated);
fprintf('  - Grid dimensions: %d x %d x %d\n', gridSize(1), gridSize(2), gridSize(3));

%% ======================== FIGURE 1: INDIVIDUAL FIELD DOSES ========================
fprintf('\n[2/8] Generating individual field dose visualization...\n');

% Determine layout based on number of fields
nCols = min(6, numValidFields);
nRows = ceil(numValidFields / nCols);

fig1 = figure('Position', [50, 50, 250*nCols, 250*nRows], 'Color', 'w');
sgtitle(sprintf('Individual Field Doses - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

% Find slice with maximum total dose for consistent visualization
[~, maxIdx] = max(totalCalculated(:));
[~, ~, maxZ] = ind2sub(gridSize, maxIdx);

% Also find center of mass of dose distribution
[X, Y, Z] = ndgrid(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
totalDoseSum = sum(totalCalculated(:));
if totalDoseSum > 0
    comX = round(sum(X(:) .* totalCalculated(:)) / totalDoseSum);
    comY = round(sum(Y(:) .* totalCalculated(:)) / totalDoseSum);
    comZ = round(sum(Z(:) .* totalCalculated(:)) / totalDoseSum);
else
    comX = round(gridSize(1)/2);
    comY = round(gridSize(2)/2);
    comZ = round(gridSize(3)/2);
end

% Global color scale
allFieldMax = 0;
for idx = validFields'
    allFieldMax = max(allFieldMax, max(fieldDoses{idx}.physicalDose(:)));
end

for i = 1:numValidFields
    subplot(nRows, nCols, i);
    fieldIdx = validFields(i);
    fieldDose = fieldDoses{fieldIdx}.physicalDose;
    
    % Display axial slice at max dose location
    sliceData = squeeze(fieldDose(:, :, maxZ));
    imagesc(sliceData);
    colormap(gca, doseCmap);
    caxis([0, allFieldMax]);
    
    title(sprintf('Field %d (%.0f°)', fieldIdx, fieldDoses{fieldIdx}.gantryAngle), 'FontSize', 9);
    axis image off;
    
    if i == numValidFields
        cb = colorbar;
        cb.Label.String = 'Dose (Gy)';
    end
end

if saveFigures
    saveFigure(fig1, fullfile(outputDir, 'individual_field_doses'), figureFormat, figureResolution);
end

%% ======================== FIGURE 2: TOTAL DOSE COMPARISON ========================
fprintf('[3/8] Generating total dose comparison...\n');

fig2 = figure('Position', [100, 100, 1400, 500], 'Color', 'w');
sgtitle(sprintf('Total Dose Comparison - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

% Use center of mass slice
displayZ = comZ;

% Calculated dose
subplot(1, 3, 1);
calcSlice = squeeze(totalCalculated(:, :, displayZ));
imagesc(calcSlice);
colormap(gca, doseCmap);
colorbar;
title(sprintf('Calculated Total Dose (Z=%d)', displayZ));
xlabel('Y'); ylabel('X');
axis image;
caxis([0, max(totalCalculated(:))]);

% Reference dose
subplot(1, 3, 2);
if ~isempty(referenceDose)
    refSlice = squeeze(referenceDose(:, :, displayZ));
    imagesc(refSlice);
    colormap(gca, doseCmap);
    colorbar;
    title(sprintf('Reference RTDOSE (Z=%d)', displayZ));
    xlabel('Y'); ylabel('X');
    axis image;
    caxis([0, max(referenceDose(:))]);
else
    text(0.5, 0.5, 'Reference dose not available', 'HorizontalAlignment', 'center');
    axis off;
end

% Difference
subplot(1, 3, 3);
if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    diffSlice = squeeze(totalCalculated(:, :, displayZ) - referenceDose(:, :, displayZ));
    
    % Create symmetric colormap for difference
    maxDiff = max(abs(diffSlice(:)));
    imagesc(diffSlice);
    colormap(gca, createDivergingColormap());
    cb = colorbar;
    cb.Label.String = 'Difference (Gy)';
    caxis([-maxDiff, maxDiff]);
    title('Difference (Calc - Ref)');
    xlabel('Y'); ylabel('X');
    axis image;
else
    text(0.5, 0.5, 'Size mismatch or no reference', 'HorizontalAlignment', 'center');
    axis off;
end

if saveFigures
    saveFigure(fig2, fullfile(outputDir, 'total_dose_comparison'), figureFormat, figureResolution);
end

%% ======================== FIGURE 3: MULTI-PLANE VIEWS ========================
fprintf('[4/8] Generating multi-plane dose views...\n');

fig3 = figure('Position', [100, 100, 1600, 800], 'Color', 'w');
sgtitle(sprintf('Multi-Plane Dose Distribution - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

maxDoseVal = max(totalCalculated(:));

% Row 1: Calculated dose in 3 planes
% Axial
subplot(2, 4, 1);
imagesc(squeeze(totalCalculated(:, :, comZ)));
colormap(gca, doseCmap); colorbar;
title(sprintf('Calculated - Axial (Z=%d)', comZ));
xlabel('Y'); ylabel('X'); axis image;
caxis([0, maxDoseVal]);

% Coronal
subplot(2, 4, 2);
imagesc(squeeze(totalCalculated(:, comY, :)));
colormap(gca, doseCmap); colorbar;
title(sprintf('Calculated - Coronal (Y=%d)', comY));
xlabel('Z'); ylabel('X'); axis image;
caxis([0, maxDoseVal]);

% Sagittal
subplot(2, 4, 3);
imagesc(squeeze(totalCalculated(comX, :, :)));
colormap(gca, doseCmap); colorbar;
title(sprintf('Calculated - Sagittal (X=%d)', comX));
xlabel('Z'); ylabel('Y'); axis image;
caxis([0, maxDoseVal]);

% 3D isodose surfaces (if enough dose)
subplot(2, 4, 4);
if max(totalCalculated(:)) > 0
    isovalues = [0.2, 0.5, 0.8, 0.95] * maxDoseVal;
    try
        hold on;
        colors = [0 0 1; 0 1 0; 1 1 0; 1 0 0];
        alphas = [0.1, 0.2, 0.3, 0.5];
        for iv = 1:length(isovalues)
            if isovalues(iv) > 0
                p = patch(isosurface(totalCalculated, isovalues(iv)));
                set(p, 'FaceColor', colors(iv,:), 'EdgeColor', 'none', 'FaceAlpha', alphas(iv));
            end
        end
        view(3); axis equal tight;
        camlight; lighting gouraud;
        title('3D Isodose (20/50/80/95%)');
    catch
        text(0.5, 0.5, '3D rendering failed', 'HorizontalAlignment', 'center');
    end
else
    text(0.5, 0.5, 'No dose data', 'HorizontalAlignment', 'center');
end
axis off;

% Row 2: Reference dose (if available)
if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    maxRefVal = max(referenceDose(:));
    
    subplot(2, 4, 5);
    imagesc(squeeze(referenceDose(:, :, comZ)));
    colormap(gca, doseCmap); colorbar;
    title(sprintf('Reference - Axial (Z=%d)', comZ));
    xlabel('Y'); ylabel('X'); axis image;
    caxis([0, maxRefVal]);
    
    subplot(2, 4, 6);
    imagesc(squeeze(referenceDose(:, comY, :)));
    colormap(gca, doseCmap); colorbar;
    title(sprintf('Reference - Coronal (Y=%d)', comY));
    xlabel('Z'); ylabel('X'); axis image;
    caxis([0, maxRefVal]);
    
    subplot(2, 4, 7);
    imagesc(squeeze(referenceDose(comX, :, :)));
    colormap(gca, doseCmap); colorbar;
    title(sprintf('Reference - Sagittal (X=%d)', comX));
    xlabel('Z'); ylabel('Y'); axis image;
    caxis([0, maxRefVal]);
    
    % Percent difference
    subplot(2, 4, 8);
    percentDiff = zeros(gridSize);
    refMask = referenceDose > 0.1 * maxRefVal;
    percentDiff(refMask) = 100 * (totalCalculated(refMask) - referenceDose(refMask)) ./ referenceDose(refMask);
    
    imagesc(squeeze(percentDiff(:, :, comZ)));
    colormap(gca, createDivergingColormap());
    cb = colorbar;
    cb.Label.String = '% Difference';
    caxis([-50, 50]);
    title('Percent Difference - Axial');
    xlabel('Y'); ylabel('X'); axis image;
end

if saveFigures
    saveFigure(fig3, fullfile(outputDir, 'multiplane_views'), figureFormat, figureResolution);
end

%% ======================== FIGURE 4: DOSE PROFILES ========================
fprintf('[5/8] Generating dose profiles...\n');

fig4 = figure('Position', [100, 100, 1400, 400], 'Color', 'w');
sgtitle(sprintf('Dose Profiles Through Dose Center - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

% X-axis values in mm
xAxis = (1:gridSize(1)) * resolution(1);
yAxis = (1:gridSize(2)) * resolution(2);
zAxis = (1:gridSize(3)) * resolution(3);

% Profile along X (through comY, comZ)
subplot(1, 3, 1);
calcProfileX = squeeze(totalCalculated(:, comY, comZ));
plot(xAxis, calcProfileX, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    refProfileX = squeeze(referenceDose(:, comY, comZ));
    plot(xAxis, refProfileX, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('X Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('X Profile (Y=%d, Z=%d)', comY, comZ));
legend('Location', 'best');
grid on;

% Profile along Y (through comX, comZ)
subplot(1, 3, 2);
calcProfileY = squeeze(totalCalculated(comX, :, comZ));
plot(yAxis, calcProfileY, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    refProfileY = squeeze(referenceDose(comX, :, comZ));
    plot(yAxis, refProfileY, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('Y Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('Y Profile (X=%d, Z=%d)', comX, comZ));
legend('Location', 'best');
grid on;

% Profile along Z (through comX, comY) - Depth dose
subplot(1, 3, 3);
calcProfileZ = squeeze(totalCalculated(comX, comY, :));
plot(zAxis, calcProfileZ, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    refProfileZ = squeeze(referenceDose(comX, comY, :));
    plot(zAxis, refProfileZ, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('Z Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('Z Profile / PDD (X=%d, Y=%d)', comX, comY));
legend('Location', 'best');
grid on;

if saveFigures
    saveFigure(fig4, fullfile(outputDir, 'dose_profiles'), figureFormat, figureResolution);
end

%% ======================== FIGURE 5: DOSE VOLUME HISTOGRAM ========================
fprintf('[6/8] Generating dose-volume histograms...\n');

fig5 = figure('Position', [100, 100, 1000, 500], 'Color', 'w');
sgtitle(sprintf('Dose-Volume Histogram - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

% Calculate cumulative DVH
numBins = 200;
maxDoseForDVH = max([max(totalCalculated(:)), max(referenceDose(:))]) * 1.05;
doseBins = linspace(0, maxDoseForDVH, numBins);

% DVH for calculated dose
calcDVH = zeros(numBins, 1);
totalVoxels = numel(totalCalculated);
for i = 1:numBins
    calcDVH(i) = 100 * sum(totalCalculated(:) >= doseBins(i)) / totalVoxels;
end

subplot(1, 2, 1);
plot(doseBins, calcDVH, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;

if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    refDVH = zeros(numBins, 1);
    for i = 1:numBins
        refDVH(i) = 100 * sum(referenceDose(:) >= doseBins(i)) / totalVoxels;
    end
    plot(doseBins, refDVH, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end

xlabel('Dose (Gy)');
ylabel('Volume (%)');
title('Cumulative DVH - Entire Volume');
legend('Location', 'northeast');
grid on;
xlim([0, maxDoseForDVH]);
ylim([0, 100]);

% Differential DVH
subplot(1, 2, 2);
calcDiffDVH = -diff(calcDVH);
plot(doseBins(1:end-1), calcDiffDVH, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;

if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    refDiffDVH = -diff(refDVH);
    plot(doseBins(1:end-1), refDiffDVH, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end

xlabel('Dose (Gy)');
ylabel('Volume Fraction (%)');
title('Differential DVH');
legend('Location', 'northeast');
grid on;
xlim([0, maxDoseForDVH]);

if saveFigures
    saveFigure(fig5, fullfile(outputDir, 'dose_volume_histogram'), figureFormat, figureResolution);
end

%% ======================== FIGURE 6: GAMMA ANALYSIS ========================
fprintf('[7/8] Performing gamma analysis...\n');

if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    
    fig6 = figure('Position', [100, 100, 1200, 800], 'Color', 'w');
    sgtitle(sprintf('Gamma Analysis - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Perform gamma analysis for different criteria
    gammaResults = struct();
    
    for gIdx = 1:length(gammaParams.doseThreshold)
        dd = gammaParams.doseThreshold(gIdx);
        dta = gammaParams.distThreshold(gIdx);
        
        fprintf('  Computing gamma (%d%%/%dmm)...\n', dd, dta);
        
        % Simplified gamma calculation
        [gammaMap, passRate] = calculateGamma3D(totalCalculated, referenceDose, ...
            dd, dta, resolution, gammaParams.doseCutoff, gammaParams.searchRadius);
        
        gammaResults(gIdx).criteria = sprintf('%d%%/%dmm', dd, dta);
        gammaResults(gIdx).passRate = passRate;
        gammaResults(gIdx).gammaMap = gammaMap;
        
        fprintf('    Pass rate: %.1f%%\n', passRate);
    end
    
    % Display gamma maps
    for gIdx = 1:length(gammaParams.doseThreshold)
        % Axial slice
        subplot(2, 2, gIdx);
        gammaSlice = squeeze(gammaResults(gIdx).gammaMap(:, :, comZ));
        
        imagesc(gammaSlice);
        colormap(gca, 'jet');
        caxis([0, 2]);
        cb = colorbar;
        cb.Label.String = 'Gamma Index';
        
        title(sprintf('Gamma %s - Pass: %.1f%%', gammaResults(gIdx).criteria, gammaResults(gIdx).passRate));
        xlabel('Y'); ylabel('X');
        axis image;
        
        % Add pass/fail overlay
        hold on;
        failMask = gammaSlice > 1;
        [failY, failX] = find(failMask);
        if ~isempty(failX)
            plot(failX, failY, 'k.', 'MarkerSize', 1);
        end
    end
    
    % Gamma histogram
    subplot(2, 2, 3);
    for gIdx = 1:length(gammaParams.doseThreshold)
        gammaVals = gammaResults(gIdx).gammaMap(:);
        gammaVals = gammaVals(gammaVals > 0 & ~isnan(gammaVals));
        
        [counts, edges] = histcounts(gammaVals, 0:0.05:3);
        centers = (edges(1:end-1) + edges(2:end)) / 2;
        
        plot(centers, counts/sum(counts)*100, 'LineWidth', 2, ...
            'DisplayName', gammaResults(gIdx).criteria);
        hold on;
    end
    
    xline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Pass Threshold');
    xlabel('Gamma Index');
    ylabel('Frequency (%)');
    title('Gamma Index Distribution');
    legend('Location', 'northeast');
    grid on;
    xlim([0, 3]);
    
    % Pass rate summary
    subplot(2, 2, 4);
    criteriaLabels = {gammaResults.criteria};
    passRates = [gammaResults.passRate];
    
    bar(passRates);
    hold on;
    yline(95, 'r--', 'LineWidth', 2, 'DisplayName', '95% threshold');
    yline(90, 'y--', 'LineWidth', 2, 'DisplayName', '90% threshold');
    
    set(gca, 'XTickLabel', criteriaLabels);
    xlabel('Gamma Criteria');
    ylabel('Pass Rate (%)');
    title('Gamma Pass Rates');
    ylim([0, 100]);
    grid on;
    
    for i = 1:length(passRates)
        text(i, passRates(i)+2, sprintf('%.1f%%', passRates(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    if saveFigures
        saveFigure(fig6, fullfile(outputDir, 'gamma_analysis'), figureFormat, figureResolution);
    end
    
else
    fprintf('  Skipping gamma analysis - reference dose not available or size mismatch\n');
    gammaResults = [];
end

%% ======================== FIGURE 7: FIELD CONTRIBUTION ANALYSIS ========================
fprintf('[8/8] Generating field contribution analysis...\n');

fig7 = figure('Position', [100, 100, 1400, 600], 'Color', 'w');
sgtitle(sprintf('Field Contribution Analysis - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

% Calculate field contributions
fieldContributions = zeros(numValidFields, 1);
fieldMaxDoses = zeros(numValidFields, 1);
gantryAngles = zeros(numValidFields, 1);

for i = 1:numValidFields
    fieldIdx = validFields(i);
    fieldDose = fieldDoses{fieldIdx}.physicalDose;
    fieldContributions(i) = sum(fieldDose(:));
    fieldMaxDoses(i) = max(fieldDose(:));
    gantryAngles(i) = fieldDoses{fieldIdx}.gantryAngle;
end

totalContribution = sum(fieldContributions);
fieldPercentages = 100 * fieldContributions / totalContribution;

% Pie chart of contributions
subplot(2, 3, 1);
labels = arrayfun(@(x,y) sprintf('F%d (%.0f°)', x, y), validFields, gantryAngles, 'UniformOutput', false);
pie(fieldPercentages);
title('Field Dose Contributions');
legend(labels, 'Location', 'eastoutside', 'FontSize', 8);

% Bar chart of contributions
subplot(2, 3, 2);
bar(fieldPercentages);
xlabel('Field Index');
ylabel('Contribution (%)');
title('Field Contributions');
set(gca, 'XTick', 1:numValidFields, 'XTickLabel', arrayfun(@(x) sprintf('%d', x), validFields, 'UniformOutput', false));
grid on;

% Max dose per field
subplot(2, 3, 3);
bar(fieldMaxDoses);
xlabel('Field Index');
ylabel('Max Dose (Gy)');
title('Maximum Dose per Field');
set(gca, 'XTick', 1:numValidFields, 'XTickLabel', arrayfun(@(x) sprintf('%d', x), validFields, 'UniformOutput', false));
grid on;

% Polar plot of gantry angles with contribution
subplot(2, 3, 4);
theta = deg2rad(gantryAngles);
rho = fieldPercentages;
polarplot(theta, rho, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
hold on;
for i = 1:length(theta)
    polarplot([0, theta(i)], [0, rho(i)], 'b-', 'LineWidth', 1.5);
end
title('Gantry Angle vs Contribution');
rlim([0, max(rho)*1.2]);

% Cumulative dose buildup
subplot(2, 3, 5);
cumulativeDose = zeros(gridSize);
cumulativeMax = zeros(numValidFields, 1);
for i = 1:numValidFields
    fieldIdx = validFields(i);
    cumulativeDose = cumulativeDose + fieldDoses{fieldIdx}.physicalDose;
    cumulativeMax(i) = max(cumulativeDose(:));
end
plot(1:numValidFields, cumulativeMax, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Number of Fields');
ylabel('Cumulative Max Dose (Gy)');
title('Dose Buildup');
grid on;
xlim([1, numValidFields]);

% Field overlap analysis
subplot(2, 3, 6);
overlapCount = zeros(gridSize);
doseThresholdForOverlap = 0.01 * max(totalCalculated(:));
for i = 1:numValidFields
    fieldIdx = validFields(i);
    fieldDose = fieldDoses{fieldIdx}.physicalDose;
    overlapCount = overlapCount + (fieldDose > doseThresholdForOverlap);
end

imagesc(squeeze(overlapCount(:, :, comZ)));
colormap(gca, 'parula');
cb = colorbar;
cb.Label.String = 'Number of Fields';
title(sprintf('Field Overlap (Z=%d)', comZ));
xlabel('Y'); ylabel('X');
axis image;

if saveFigures
    saveFigure(fig7, fullfile(outputDir, 'field_contribution_analysis'), figureFormat, figureResolution);
end

%% ======================== STATISTICAL SUMMARY ========================
fprintf('\n=======================================================\n');
fprintf('  Statistical Summary\n');
fprintf('=======================================================\n');

% Basic statistics
fprintf('\nCalculated Dose Statistics:\n');
fprintf('  Max dose: %.4f Gy\n', max(totalCalculated(:)));
fprintf('  Mean dose: %.4f Gy\n', mean(totalCalculated(:)));
fprintf('  Median dose: %.4f Gy\n', median(totalCalculated(:)));
fprintf('  Std dev: %.4f Gy\n', std(totalCalculated(:)));

if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    fprintf('\nReference Dose Statistics:\n');
    fprintf('  Max dose: %.4f Gy\n', max(referenceDose(:)));
    fprintf('  Mean dose: %.4f Gy\n', mean(referenceDose(:)));
    fprintf('  Median dose: %.4f Gy\n', median(referenceDose(:)));
    fprintf('  Std dev: %.4f Gy\n', std(referenceDose(:)));
    
    % Comparison statistics
    doseDiff = totalCalculated - referenceDose;
    
    fprintf('\nComparison Statistics:\n');
    fprintf('  Mean difference: %.4f Gy\n', mean(doseDiff(:)));
    fprintf('  Mean absolute difference: %.4f Gy\n', mean(abs(doseDiff(:))));
    fprintf('  Max difference: %.4f Gy\n', max(abs(doseDiff(:))));
    fprintf('  RMS difference: %.4f Gy\n', sqrt(mean(doseDiff(:).^2)));
    
    % High dose region statistics
    highDoseMask = referenceDose > 0.5 * max(referenceDose(:));
    if any(highDoseMask(:))
        highDoseDiff = doseDiff(highDoseMask);
        highDoseRef = referenceDose(highDoseMask);
        relDiff = 100 * highDoseDiff ./ highDoseRef;
        
        fprintf('\nHigh Dose Region (>50%% max):\n');
        fprintf('  Voxels: %d\n', sum(highDoseMask(:)));
        fprintf('  Mean relative difference: %.2f%%\n', mean(relDiff));
        fprintf('  Max relative difference: %.2f%%\n', max(abs(relDiff)));
        
        % Correlation
        corrMatrix = corrcoef(totalCalculated(highDoseMask), referenceDose(highDoseMask));
        fprintf('  Correlation coefficient: %.4f\n', corrMatrix(1,2));
    end
    
    % Gamma results
    if exist('gammaResults', 'var') && ~isempty(gammaResults)
        fprintf('\nGamma Analysis Results:\n');
        for gIdx = 1:length(gammaResults)
            fprintf('  %s: %.1f%% pass rate\n', gammaResults(gIdx).criteria, gammaResults(gIdx).passRate);
        end
    end
end

fprintf('\nField Statistics:\n');
fprintf('  Number of fields: %d\n', numValidFields);
for i = 1:numValidFields
    fieldIdx = validFields(i);
    fprintf('  Field %d (%.0f°): Max=%.4f Gy, Contribution=%.1f%%\n', ...
        fieldIdx, gantryAngles(i), fieldMaxDoses(i), fieldPercentages(i));
end

%% ======================== SAVE SUMMARY REPORT ========================
fprintf('\nSaving summary report...\n');

% Create summary structure
summary = struct();
summary.patientID = patientID;
summary.sessionName = sessionName;
summary.timestamp = datetime('now');
summary.gridSize = gridSize;
summary.resolution = resolution;
summary.numFields = numValidFields;
summary.isResampled = isResampled;

summary.calculatedDose.max = max(totalCalculated(:));
summary.calculatedDose.mean = mean(totalCalculated(:));
summary.calculatedDose.std = std(totalCalculated(:));

if ~isempty(referenceDose) && isequal(size(totalCalculated), size(referenceDose))
    summary.referenceDose.max = max(referenceDose(:));
    summary.referenceDose.mean = mean(referenceDose(:));
    summary.referenceDose.std = std(referenceDose(:));
    
    summary.comparison.meanDiff = mean(doseDiff(:));
    summary.comparison.meanAbsDiff = mean(abs(doseDiff(:)));
    summary.comparison.maxAbsDiff = max(abs(doseDiff(:)));
    summary.comparison.rmsDiff = sqrt(mean(doseDiff(:).^2));
    
    if exist('gammaResults', 'var') && ~isempty(gammaResults)
        for gIdx = 1:length(gammaResults)
            fieldName = sprintf('gamma_%d_%d', gammaParams.doseThreshold(gIdx), gammaParams.distThreshold(gIdx));
            summary.gamma.(fieldName) = gammaResults(gIdx).passRate;
        end
    end
end

summary.fields = struct();
for i = 1:numValidFields
    fieldIdx = validFields(i);
    fieldName = sprintf('field_%d', fieldIdx);
    summary.fields.(fieldName).gantryAngle = gantryAngles(i);
    summary.fields.(fieldName).maxDose = fieldMaxDoses(i);
    summary.fields.(fieldName).contribution = fieldPercentages(i);
end

% Save summary
save(fullfile(outputDir, 'analysis_summary.mat'), 'summary');

% Write text report
reportFile = fullfile(outputDir, 'analysis_report.txt');
fid = fopen(reportFile, 'w');
fprintf(fid, 'ETHOS Field Dose Analysis Report\n');
fprintf(fid, '================================\n\n');
fprintf(fid, 'Patient ID: %s\n', patientID);
fprintf(fid, 'Session: %s\n', sessionName);
fprintf(fid, 'Generated: %s\n\n', char(summary.timestamp));
fprintf(fid, 'Grid Size: %d x %d x %d\n', gridSize(1), gridSize(2), gridSize(3));
fprintf(fid, 'Resolution: [%.2f, %.2f, %.2f] mm\n\n', resolution(1), resolution(2), resolution(3));

fprintf(fid, 'CALCULATED DOSE\n');
fprintf(fid, '  Max: %.4f Gy\n', summary.calculatedDose.max);
fprintf(fid, '  Mean: %.4f Gy\n', summary.calculatedDose.mean);
fprintf(fid, '  Std: %.4f Gy\n\n', summary.calculatedDose.std);

if isfield(summary, 'referenceDose')
    fprintf(fid, 'REFERENCE DOSE\n');
    fprintf(fid, '  Max: %.4f Gy\n', summary.referenceDose.max);
    fprintf(fid, '  Mean: %.4f Gy\n', summary.referenceDose.mean);
    fprintf(fid, '  Std: %.4f Gy\n\n', summary.referenceDose.std);
    
    fprintf(fid, 'COMPARISON\n');
    fprintf(fid, '  Mean Difference: %.4f Gy\n', summary.comparison.meanDiff);
    fprintf(fid, '  Mean Abs Difference: %.4f Gy\n', summary.comparison.meanAbsDiff);
    fprintf(fid, '  Max Abs Difference: %.4f Gy\n', summary.comparison.maxAbsDiff);
    fprintf(fid, '  RMS Difference: %.4f Gy\n\n', summary.comparison.rmsDiff);
    
    if isfield(summary, 'gamma')
        fprintf(fid, 'GAMMA ANALYSIS\n');
        gammaFields = fieldnames(summary.gamma);
        for gf = 1:length(gammaFields)
            fprintf(fid, '  %s: %.1f%% pass\n', gammaFields{gf}, summary.gamma.(gammaFields{gf}));
        end
        fprintf(fid, '\n');
    end
end

fprintf(fid, 'FIELD CONTRIBUTIONS\n');
for i = 1:numValidFields
    fieldIdx = validFields(i);
    fprintf(fid, '  Field %d (%.0f°): Max=%.4f Gy, Contribution=%.1f%%\n', ...
        fieldIdx, gantryAngles(i), fieldMaxDoses(i), fieldPercentages(i));
end

fclose(fid);

fprintf('\n=======================================================\n');
fprintf('  Visualization Complete\n');
fprintf('=======================================================\n');
fprintf('  Output directory: %s\n', outputDir);
fprintf('  Figures saved: %d\n', 7);
fprintf('  Summary saved: analysis_summary.mat, analysis_report.txt\n');
fprintf('=======================================================\n\n');

%% ======================== HELPER FUNCTIONS ========================

function saveFigure(fig, filename, format, resolution)
    % Save figure in specified format
    switch lower(format)
        case 'png'
            print(fig, filename, '-dpng', sprintf('-r%d', resolution));
        case 'pdf'
            print(fig, filename, '-dpdf', '-bestfit');
        case 'eps'
            print(fig, filename, '-depsc', '-tiff');
        case 'fig'
            savefig(fig, [filename, '.fig']);
        otherwise
            print(fig, filename, '-dpng', sprintf('-r%d', resolution));
    end
end

function cmap = createDivergingColormap(n)
    % Create blue-white-red diverging colormap
    if nargin < 1
        n = 256;
    end
    
    half = floor(n/2);
    
    % Blue to white
    r1 = linspace(0, 1, half)';
    g1 = linspace(0, 1, half)';
    b1 = ones(half, 1);
    
    % White to red
    r2 = ones(n-half, 1);
    g2 = linspace(1, 0, n-half)';
    b2 = linspace(1, 0, n-half)';
    
    cmap = [r1, g1, b1; r2, g2, b2];
end

function [gammaMap, passRate] = calculateGamma3D(calcDose, refDose, ddCrit, dtaCrit, resolution, doseCutoff, searchRadius)
    % Simplified 3D gamma calculation
    % ddCrit: dose difference criterion (%)
    % dtaCrit: distance-to-agreement criterion (mm)
    % resolution: [dx, dy, dz] in mm
    % doseCutoff: minimum dose threshold (% of max)
    % searchRadius: maximum search distance (mm)
    
    gridSize = size(refDose);
    gammaMap = nan(gridSize);
    
    % Normalize doses
    maxRefDose = max(refDose(:));
    doseCutoffAbs = doseCutoff / 100 * maxRefDose;
    
    % Create mask for evaluation
    evalMask = refDose > doseCutoffAbs;
    
    if ~any(evalMask(:))
        passRate = 100;
        return;
    end
    
    % Convert criteria
    ddCritAbs = ddCrit / 100 * maxRefDose;  % Absolute dose criterion
    
    % Search radius in voxels
    searchVoxels = ceil([searchRadius/resolution(1), searchRadius/resolution(2), searchRadius/resolution(3)]);
    
    % Get evaluation points
    [evalX, evalY, evalZ] = ind2sub(gridSize, find(evalMask));
    numEvalPoints = length(evalX);
    
    % Calculate gamma for each evaluation point
    gammaValues = zeros(numEvalPoints, 1);
    
    for i = 1:numEvalPoints
        x = evalX(i);
        y = evalY(i);
        z = evalZ(i);
        
        refVal = refDose(x, y, z);
        
        % Define search region
        xMin = max(1, x - searchVoxels(1));
        xMax = min(gridSize(1), x + searchVoxels(1));
        yMin = max(1, y - searchVoxels(2));
        yMax = min(gridSize(2), y + searchVoxels(2));
        zMin = max(1, z - searchVoxels(3));
        zMax = min(gridSize(3), z + searchVoxels(3));
        
        minGamma = inf;
        
        for xi = xMin:xMax
            for yi = yMin:yMax
                for zi = zMin:zMax
                    % Distance in mm
                    dist = sqrt(((xi-x)*resolution(1))^2 + ...
                               ((yi-y)*resolution(2))^2 + ...
                               ((zi-z)*resolution(3))^2);
                    
                    if dist <= searchRadius
                        calcVal = calcDose(xi, yi, zi);
                        
                        % Gamma calculation
                        doseTerm = ((calcVal - refVal) / ddCritAbs)^2;
                        distTerm = (dist / dtaCrit)^2;
                        gamma = sqrt(doseTerm + distTerm);
                        
                        minGamma = min(minGamma, gamma);
                    end
                end
            end
        end
        
        gammaValues(i) = minGamma;
        gammaMap(x, y, z) = minGamma;
    end
    
    % Calculate pass rate
    passRate = 100 * sum(gammaValues <= 1) / numEvalPoints;
end