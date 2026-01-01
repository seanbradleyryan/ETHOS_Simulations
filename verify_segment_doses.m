%% ETHOS Segment Dose Verification Visualizer
% Purpose: Visualize and verify output from calculate_field_dose.m
%
% This script provides comprehensive visualization tools to verify:
%   1. Segment dose distributions
%   2. Beam dose accumulations
%   3. Total dose vs reference comparison
%   4. MLC aperture shapes
%   5. Dose profiles (depth and lateral)
%   6. Gamma analysis
%   7. DVH comparisons
%   8. Segment weight analysis
%
% Usage:
%   1. Run calculate_field_dose.m first
%   2. Set patientID and sessionName below
%   3. Run this script
%
% Author: Generated for ETHOS dose analysis
% Date: 2025

clear; clc; close all;

%% Configuration
patientID = '1194203';
sessionName = 'Session_1';

% Base directory
baseDir = '/mnt/weka/home/80030361/ETHOS_Simulations';
dataDir = fullfile(baseDir, 'SegmentDoses', patientID, sessionName);
outputDir = fullfile(dataDir, 'Visualizations');

% Visualization options
saveAllFigures = true;
figureFormat = 'png';  % 'png', 'fig', or 'both'
dpi = 150;

% Gamma analysis parameters
gammaDoseThreshold = 3;    % % dose difference criterion
gammaDTAThreshold = 3;     % mm distance-to-agreement criterion
gammaDoseCutoff = 10;      % % of max dose (exclude low dose regions)

%% Create output directory
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load data
fprintf('=======================================================\n');
fprintf('  ETHOS Segment Dose Verification Visualizer\n');
fprintf('=======================================================\n\n');

fprintf('[1/9] Loading data...\n');

% Load segment doses (RTDOSE grid)
segmentDoseFile = fullfile(dataDir, 'segmentDoses.mat');
if exist(segmentDoseFile, 'file')
    fprintf('  Loading: segmentDoses.mat\n');
    data = load(segmentDoseFile);
    segmentDoses = data.segmentDosesResampled;
    beamDoses = data.beamDosesResampled;
    totalDose = data.totalDoseResampled;
    referenceDose = data.referenceDose;
    doseGrid = data.doseGrid;
    useResampled = true;
else
    % Fall back to CT grid data
    segmentDoseFile = fullfile(dataDir, 'segmentDoses_CTgrid.mat');
    if exist(segmentDoseFile, 'file')
        fprintf('  Loading: segmentDoses_CTgrid.mat (CT grid)\n');
        data = load(segmentDoseFile);
        segmentDoses = data.segmentDoses;
        beamDoses = data.beamDoses;
        totalDose = data.totalDose;
        referenceDose = [];
        doseGrid = [];
        useResampled = false;
    else
        error('No segment dose file found in %s', dataDir);
    end
end

% Load segment data (RTPLAN info)
segmentDataFile = fullfile(dataDir, 'segmentData.mat');
if exist(segmentDataFile, 'file')
    fprintf('  Loading: segmentData.mat\n');
    segData = load(segmentDataFile);
    segmentInfo = segData.segmentData;
else
    segmentInfo = [];
    fprintf('  WARNING: segmentData.mat not found\n');
end

% Load comparison if available
comparisonFile = fullfile(dataDir, 'doseComparison.mat');
if exist(comparisonFile, 'file')
    fprintf('  Loading: doseComparison.mat\n');
    compData = load(comparisonFile);
    comparison = compData.comparison;
else
    comparison = [];
end

% Get grid dimensions
gridSize = size(totalDose);
fprintf('  Dose grid size: %d x %d x %d\n', gridSize(1), gridSize(2), gridSize(3));

% Calculate grid coordinates (in mm)
if ~isempty(doseGrid)
    dx = doseGrid.resolution(1);
    dy = doseGrid.resolution(2);
    dz = doseGrid.resolution(3);
else
    dx = 2.5; dy = 2.5; dz = 2.5;  % Default
end
xVec = (0:gridSize(2)-1) * dx;
yVec = (0:gridSize(1)-1) * dy;
zVec = (0:gridSize(3)-1) * dz;

% Count valid data
numBeams = length(beamDoses);
numSegments = sum(~cellfun(@isempty, segmentDoses));
fprintf('  Beams: %d, Segments: %d\n\n', numBeams, numSegments);

%% [2/9] Overview: Total Dose vs Reference
fprintf('[2/9] Generating dose overview comparison...\n');

% Find slice with maximum dose
[maxDose, maxIdx] = max(totalDose(:));
[maxRow, maxCol, maxSlice] = ind2sub(gridSize, maxIdx);

fig1 = figure('Name', 'Dose Overview', 'Position', [50, 50, 1600, 900], 'Color', 'w');

% Calculated dose - Axial
subplot(2, 4, 1);
imagesc(xVec, yVec, squeeze(totalDose(:, :, maxSlice)));
colorbar;
colormap(gca, 'jet');
title(sprintf('Calculated - Axial (Z=%d)', maxSlice));
xlabel('X (mm)'); ylabel('Y (mm)');
axis image;
hold on;
plot(maxCol*dx, maxRow*dy, 'w+', 'MarkerSize', 15, 'LineWidth', 2);

% Calculated dose - Coronal
subplot(2, 4, 2);
imagesc(xVec, zVec, squeeze(totalDose(maxRow, :, :))');
colorbar;
colormap(gca, 'jet');
title(sprintf('Calculated - Coronal (Y=%d)', maxRow));
xlabel('X (mm)'); ylabel('Z (mm)');
axis image;

% Calculated dose - Sagittal
subplot(2, 4, 3);
imagesc(yVec, zVec, squeeze(totalDose(:, maxCol, :))');
colorbar;
colormap(gca, 'jet');
title(sprintf('Calculated - Sagittal (X=%d)', maxCol));
xlabel('Y (mm)'); ylabel('Z (mm)');
axis image;

% Dose histogram
subplot(2, 4, 4);
validDose = totalDose(totalDose > 0.01 * maxDose);
histogram(validDose, 50, 'FaceColor', [0.2 0.4 0.8]);
xlabel('Dose (Gy)');
ylabel('Voxels');
title('Dose Histogram (Calculated)');
grid on;

if ~isempty(referenceDose)
    % Reference dose - Axial
    subplot(2, 4, 5);
    imagesc(xVec, yVec, squeeze(referenceDose(:, :, maxSlice)));
    colorbar;
    colormap(gca, 'jet');
    title(sprintf('Reference - Axial (Z=%d)', maxSlice));
    xlabel('X (mm)'); ylabel('Y (mm)');
    axis image;
    
    % Reference dose - Coronal
    subplot(2, 4, 6);
    imagesc(xVec, zVec, squeeze(referenceDose(maxRow, :, :))');
    colorbar;
    colormap(gca, 'jet');
    title(sprintf('Reference - Coronal (Y=%d)', maxRow));
    xlabel('X (mm)'); ylabel('Z (mm)');
    axis image;
    
    % Difference
    subplot(2, 4, 7);
    doseDiff = totalDose - referenceDose;
    maxAbsDiff = max(abs(doseDiff(:)));
    imagesc(xVec, yVec, squeeze(doseDiff(:, :, maxSlice)));
    colorbar;
    colormap(gca, 'coolwarm');
    caxis([-maxAbsDiff, maxAbsDiff]);
    title('Difference (Calc - Ref)');
    xlabel('X (mm)'); ylabel('Y (mm)');
    axis image;
    
    % Scatter plot: Calculated vs Reference
    subplot(2, 4, 8);
    mask = referenceDose > 0.1 * max(referenceDose(:));
    scatter(referenceDose(mask), totalDose(mask), 1, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    maxVal = max([max(referenceDose(:)), max(totalDose(:))]);
    plot([0, maxVal], [0, maxVal], 'r-', 'LineWidth', 2);
    xlabel('Reference Dose (Gy)');
    ylabel('Calculated Dose (Gy)');
    title('Dose Correlation');
    axis equal;
    xlim([0, maxVal]);
    ylim([0, maxVal]);
    grid on;
    
    % Calculate correlation
    corrCoef = corrcoef(referenceDose(mask), totalDose(mask));
    text(0.1*maxVal, 0.9*maxVal, sprintf('R = %.4f', corrCoef(1,2)), 'FontSize', 12);
end

sgtitle(sprintf('Dose Overview: Patient %s - %s', patientID, sessionName), 'FontSize', 14);

if saveAllFigures
    saveFigure(fig1, fullfile(outputDir, '01_DoseOverview'), figureFormat, dpi);
end

%% [3/9] Beam-by-Beam Dose Contributions
fprintf('[3/9] Generating beam dose contributions...\n');

fig2 = figure('Name', 'Beam Contributions', 'Position', [50, 50, 1600, 900], 'Color', 'w');

numValidBeams = sum(~cellfun(@isempty, beamDoses));
nRows = ceil(sqrt(numValidBeams + 1));
nCols = ceil((numValidBeams + 1) / nRows);

beamIdx = 0;
beamMaxDoses = zeros(numValidBeams, 1);
beamNames = cell(numValidBeams, 1);

for i = 1:length(beamDoses)
    if ~isempty(beamDoses{i})
        beamIdx = beamIdx + 1;
        
        subplot(nRows, nCols, beamIdx);
        beamDose = beamDoses{i}.physicalDose;
        imagesc(xVec, yVec, squeeze(beamDose(:, :, maxSlice)));
        colorbar;
        colormap(gca, 'jet');
        
        titleStr = sprintf('Beam %d: %.1f°', i, beamDoses{i}.gantryAngle);
        if isfield(beamDoses{i}, 'beamName')
            titleStr = sprintf('%s\n%s', beamDoses{i}.beamName, titleStr);
        end
        title(titleStr, 'FontSize', 9);
        xlabel('X (mm)'); ylabel('Y (mm)');
        axis image;
        
        beamMaxDoses(beamIdx) = beamDoses{i}.maxDose;
        beamNames{beamIdx} = sprintf('B%d (%.0f°)', i, beamDoses{i}.gantryAngle);
    end
end

% Bar chart of beam contributions
subplot(nRows, nCols, beamIdx + 1);
bar(beamMaxDoses, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', beamNames, 'XTickLabelRotation', 45);
ylabel('Max Dose (Gy)');
title('Beam Max Doses');
grid on;

sgtitle(sprintf('Beam Dose Contributions (Axial Z=%d)', maxSlice), 'FontSize', 14);

if saveAllFigures
    saveFigure(fig2, fullfile(outputDir, '02_BeamContributions'), figureFormat, dpi);
end

%% [4/9] Segment Analysis
fprintf('[4/9] Generating segment analysis...\n');

if ~isempty(segmentInfo)
    fig3 = figure('Name', 'Segment Analysis', 'Position', [50, 50, 1400, 800], 'Color', 'w');
    
    % Collect segment data
    allSegmentWeights = [];
    allSegmentMetersets = [];
    segmentBeamLabels = [];
    
    for beamIdx = 1:length(segmentInfo.beams)
        beam = segmentInfo.beams{beamIdx};
        numSegs = length(beam.segments);
        
        for segIdx = 1:numSegs
            seg = beam.segments{segIdx};
            allSegmentWeights(end+1) = seg.segmentWeight;
            allSegmentMetersets(end+1) = seg.segmentMeterset;
            segmentBeamLabels(end+1) = beamIdx;
        end
    end
    
    % Segment weight distribution
    subplot(2, 3, 1);
    histogram(allSegmentWeights, 30, 'FaceColor', [0.4 0.7 0.4]);
    xlabel('Segment Weight');
    ylabel('Count');
    title('Segment Weight Distribution');
    grid on;
    
    % Segment meterset distribution
    subplot(2, 3, 2);
    histogram(allSegmentMetersets, 30, 'FaceColor', [0.7 0.4 0.4]);
    xlabel('Segment Meterset (MU)');
    ylabel('Count');
    title('Segment Meterset Distribution');
    grid on;
    
    % Segments per beam
    subplot(2, 3, 3);
    segsPerBeam = zeros(length(segmentInfo.beams), 1);
    for i = 1:length(segmentInfo.beams)
        segsPerBeam(i) = length(segmentInfo.beams{i}.segments);
    end
    bar(segsPerBeam, 'FaceColor', [0.5 0.5 0.8]);
    xlabel('Beam Index');
    ylabel('Number of Segments');
    title('Segments per Beam');
    grid on;
    
    % Cumulative segment weights by beam
    subplot(2, 3, 4);
    hold on;
    colors = lines(length(segmentInfo.beams));
    for beamIdx = 1:length(segmentInfo.beams)
        beam = segmentInfo.beams{beamIdx};
        weights = cellfun(@(s) s.segmentWeight, beam.segments);
        cumWeights = cumsum(weights);
        plot(1:length(weights), cumWeights, '-o', 'Color', colors(beamIdx,:), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('Beam %d', beamIdx));
    end
    xlabel('Segment Index');
    ylabel('Cumulative Weight');
    title('Cumulative Segment Weights');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    
    % Beam meterset pie chart
    subplot(2, 3, 5);
    beamMetersets = zeros(length(segmentInfo.beams), 1);
    beamLabels = cell(length(segmentInfo.beams), 1);
    for i = 1:length(segmentInfo.beams)
        beamMetersets(i) = segmentInfo.beams{i}.beamMeterset;
        beamLabels{i} = sprintf('B%d: %.0f MU', i, beamMetersets(i));
    end
    pie(beamMetersets, beamLabels);
    title('Beam Meterset Distribution');
    
    % Summary statistics text
    subplot(2, 3, 6);
    axis off;
    textStr = {
        sprintf('Total Beams: %d', length(segmentInfo.beams));
        sprintf('Total Segments: %d', segmentInfo.totalSegments);
        sprintf('Total Meterset: %.1f MU', sum(beamMetersets));
        '';
        sprintf('Avg Segments/Beam: %.1f', mean(segsPerBeam));
        sprintf('Min Segment Weight: %.4f', min(allSegmentWeights));
        sprintf('Max Segment Weight: %.4f', max(allSegmentWeights));
        sprintf('Mean Segment Weight: %.4f', mean(allSegmentWeights));
    };
    text(0.1, 0.9, textStr, 'FontSize', 11, 'VerticalAlignment', 'top', ...
        'FontName', 'FixedWidth');
    title('Summary Statistics');
    
    sgtitle(sprintf('Segment Analysis: Patient %s', patientID), 'FontSize', 14);
    
    if saveAllFigures
        saveFigure(fig3, fullfile(outputDir, '03_SegmentAnalysis'), figureFormat, dpi);
    end
end

%% [5/9] Individual Segment Doses (Sample)
fprintf('[5/9] Generating sample segment dose display...\n');

% Show first few segments from first beam
numSampleSegments = min(9, numSegments);
fig4 = figure('Name', 'Sample Segment Doses', 'Position', [50, 50, 1200, 1000], 'Color', 'w');

segCount = 0;
nRowsSeg = ceil(sqrt(numSampleSegments));
nColsSeg = ceil(numSampleSegments / nRowsSeg);

for i = 1:length(segmentDoses)
    if ~isempty(segmentDoses{i}) && segCount < numSampleSegments
        segCount = segCount + 1;
        
        subplot(nRowsSeg, nColsSeg, segCount);
        segDose = segmentDoses{i}.physicalDose;
        imagesc(xVec, yVec, squeeze(segDose(:, :, maxSlice)));
        colorbar;
        colormap(gca, 'jet');
        
        title(sprintf('B%d-S%d (W=%.3f)', ...
            segmentDoses{i}.beamIdx, segmentDoses{i}.segmentIdx, ...
            segmentDoses{i}.segmentWeight), 'FontSize', 9);
        xlabel('X'); ylabel('Y');
        axis image;
    end
end

sgtitle(sprintf('Sample Segment Doses (Axial Z=%d)', maxSlice), 'FontSize', 14);

if saveAllFigures
    saveFigure(fig4, fullfile(outputDir, '04_SampleSegmentDoses'), figureFormat, dpi);
end

%% [6/9] Dose Profiles
fprintf('[6/9] Generating dose profiles...\n');

fig5 = figure('Name', 'Dose Profiles', 'Position', [50, 50, 1400, 800], 'Color', 'w');

% Central axis depth dose (along Y through max)
subplot(2, 3, 1);
calcProfile = squeeze(totalDose(maxRow, maxCol, :));
plot(zVec, calcProfile, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if ~isempty(referenceDose)
    refProfile = squeeze(referenceDose(maxRow, maxCol, :));
    plot(zVec, refProfile, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('Z Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('Depth Profile (X=%d, Y=%d)', maxCol, maxRow));
legend('Location', 'best');
grid on;

% Lateral profile X (through max)
subplot(2, 3, 2);
calcProfileX = squeeze(totalDose(maxRow, :, maxSlice));
plot(xVec, calcProfileX, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if ~isempty(referenceDose)
    refProfileX = squeeze(referenceDose(maxRow, :, maxSlice));
    plot(xVec, refProfileX, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('X Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('Lateral Profile X (Y=%d, Z=%d)', maxRow, maxSlice));
legend('Location', 'best');
grid on;

% Lateral profile Y (through max)
subplot(2, 3, 3);
calcProfileY = squeeze(totalDose(:, maxCol, maxSlice));
plot(yVec, calcProfileY, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if ~isempty(referenceDose)
    refProfileY = squeeze(referenceDose(:, maxCol, maxSlice));
    plot(yVec, refProfileY, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('Y Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('Lateral Profile Y (X=%d, Z=%d)', maxCol, maxSlice));
legend('Location', 'best');
grid on;

% Profile differences
if ~isempty(referenceDose)
    subplot(2, 3, 4);
    plot(zVec, calcProfile - refProfile, 'k-', 'LineWidth', 1.5);
    xlabel('Z Position (mm)');
    ylabel('Dose Difference (Gy)');
    title('Depth Profile Difference');
    grid on;
    yline(0, 'r--');
    
    subplot(2, 3, 5);
    plot(xVec, calcProfileX - refProfileX, 'k-', 'LineWidth', 1.5);
    xlabel('X Position (mm)');
    ylabel('Dose Difference (Gy)');
    title('Lateral X Profile Difference');
    grid on;
    yline(0, 'r--');
    
    subplot(2, 3, 6);
    plot(yVec, calcProfileY - refProfileY, 'k-', 'LineWidth', 1.5);
    xlabel('Y Position (mm)');
    ylabel('Dose Difference (Gy)');
    title('Lateral Y Profile Difference');
    grid on;
    yline(0, 'r--');
end

sgtitle('Dose Profiles Through Maximum Dose Point', 'FontSize', 14);

if saveAllFigures
    saveFigure(fig5, fullfile(outputDir, '05_DoseProfiles'), figureFormat, dpi);
end

%% [7/9] Gamma Analysis
fprintf('[7/9] Performing gamma analysis...\n');

if ~isempty(referenceDose)
    fig6 = figure('Name', 'Gamma Analysis', 'Position', [50, 50, 1400, 600], 'Color', 'w');
    
    % Calculate gamma index (simplified 2D for display slice)
    refSlice = squeeze(referenceDose(:, :, maxSlice));
    calcSlice = squeeze(totalDose(:, :, maxSlice));
    
    % Normalize doses
    maxRefDose = max(referenceDose(:));
    
    % Dose difference criterion (as fraction)
    doseCrit = gammaDoseThreshold / 100;
    
    % DTA criterion (in voxels)
    dtaCritVoxels = gammaDTAThreshold / dx;
    
    % Calculate gamma (simplified version)
    [gammaMap, passRate] = calculateGamma2D(refSlice, calcSlice, maxRefDose, ...
        doseCrit, dtaCritVoxels, gammaDoseCutoff/100);
    
    % Gamma map
    subplot(1, 3, 1);
    imagesc(xVec, yVec, gammaMap);
    colorbar;
    caxis([0, 2]);
    colormap(gca, 'jet');
    title(sprintf('Gamma Map (%.0f%%/%.0fmm)', gammaDoseThreshold, gammaDTAThreshold));
    xlabel('X (mm)'); ylabel('Y (mm)');
    axis image;
    
    % Gamma pass/fail map
    subplot(1, 3, 2);
    passMap = gammaMap <= 1;
    imagesc(xVec, yVec, passMap);
    colorbar;
    colormap(gca, [1 0 0; 0 1 0]);  % Red=fail, Green=pass
    title(sprintf('Pass/Fail (Pass Rate: %.1f%%)', passRate));
    xlabel('X (mm)'); ylabel('Y (mm)');
    axis image;
    
    % Gamma histogram
    subplot(1, 3, 3);
    validGamma = gammaMap(gammaMap < 10 & ~isnan(gammaMap));
    histogram(validGamma, 50, 'FaceColor', [0.3 0.5 0.8]);
    xlabel('Gamma Index');
    ylabel('Voxels');
    title('Gamma Distribution');
    xline(1, 'r-', 'LineWidth', 2);
    grid on;
    
    % Add statistics text
    text(1.5, 0.8*max(histcounts(validGamma, 50)), ...
        {sprintf('Mean: %.2f', mean(validGamma)), ...
         sprintf('Median: %.2f', median(validGamma)), ...
         sprintf('Pass: %.1f%%', passRate)}, 'FontSize', 10);
    
    sgtitle(sprintf('Gamma Analysis: %d%%/%dmm (Cutoff: %d%%)', ...
        gammaDoseThreshold, gammaDTAThreshold, gammaDoseCutoff), 'FontSize', 14);
    
    if saveAllFigures
        saveFigure(fig6, fullfile(outputDir, '06_GammaAnalysis'), figureFormat, dpi);
    end
    
    fprintf('  Gamma pass rate (%.0f%%/%.0fmm): %.1f%%\n', ...
        gammaDoseThreshold, gammaDTAThreshold, passRate);
end

%% [8/9] MLC Aperture Visualization
fprintf('[8/9] Generating MLC aperture visualization...\n');

if ~isempty(segmentInfo)
    fig7 = figure('Name', 'MLC Apertures', 'Position', [50, 50, 1400, 800], 'Color', 'w');
    
    % Show MLC apertures for first beam's segments
    beam1 = segmentInfo.beams{1};
    numSegsToShow = min(6, length(beam1.segments));
    
    if beam1.numLeafPairs > 0 && ~isempty(beam1.leafBoundaries)
        leafBoundaries = beam1.leafBoundaries;
        
        for segIdx = 1:numSegsToShow
            subplot(2, 3, segIdx);
            
            seg = beam1.segments{segIdx};
            mlcPos = seg.mlcPositions;
            
            if ~isempty(mlcPos)
                numLeaves = beam1.numLeafPairs;
                leftBank = mlcPos(1:numLeaves);
                rightBank = mlcPos(numLeaves+1:end);
                
                % Draw leaf pairs
                hold on;
                for leafIdx = 1:numLeaves
                    yBottom = leafBoundaries(leafIdx);
                    yTop = leafBoundaries(leafIdx + 1);
                    
                    % Left leaf (closed region)
                    fill([-200, leftBank(leafIdx), leftBank(leafIdx), -200], ...
                         [yBottom, yBottom, yTop, yTop], ...
                         [0.5 0.5 0.5], 'EdgeColor', 'k');
                    
                    % Right leaf (closed region)
                    fill([rightBank(leafIdx), 200, 200, rightBank(leafIdx)], ...
                         [yBottom, yBottom, yTop, yTop], ...
                         [0.5 0.5 0.5], 'EdgeColor', 'k');
                    
                    % Open region
                    fill([leftBank(leafIdx), rightBank(leafIdx), rightBank(leafIdx), leftBank(leafIdx)], ...
                         [yBottom, yBottom, yTop, yTop], ...
                         [1 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                end
                
                % Draw jaws if available
                if ~isempty(seg.jawX)
                    xline(seg.jawX(1), 'r-', 'LineWidth', 2);
                    xline(seg.jawX(2), 'r-', 'LineWidth', 2);
                end
                if ~isempty(seg.jawY)
                    yline(seg.jawY(1), 'b-', 'LineWidth', 2);
                    yline(seg.jawY(2), 'b-', 'LineWidth', 2);
                end
                
                xlim([-150, 150]);
                ylim([min(leafBoundaries)-10, max(leafBoundaries)+10]);
                xlabel('X (mm)');
                ylabel('Y (mm)');
                title(sprintf('Segment %d (W=%.3f)', segIdx, seg.segmentWeight));
                axis equal;
                grid on;
            else
                text(0.5, 0.5, 'No MLC data', 'HorizontalAlignment', 'center');
                axis off;
            end
        end
        
        sgtitle(sprintf('MLC Apertures - Beam 1 (%s)', beam1.beamName), 'FontSize', 14);
    else
        text(0.5, 0.5, 'No MLC boundary data available', ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
        axis off;
    end
    
    if saveAllFigures
        saveFigure(fig7, fullfile(outputDir, '07_MLCApertures'), figureFormat, dpi);
    end
end

%% [9/9] Summary Statistics and Report
fprintf('[9/9] Generating summary report...\n');

fig8 = figure('Name', 'Summary Report', 'Position', [50, 50, 1000, 800], 'Color', 'w');
axis off;

% Build report text
reportLines = {
    '═══════════════════════════════════════════════════════════════';
    sprintf('  SEGMENT DOSE CALCULATION VERIFICATION REPORT');
    '═══════════════════════════════════════════════════════════════';
    '';
    sprintf('  Patient ID:        %s', patientID);
    sprintf('  Session:           %s', sessionName);
    sprintf('  Analysis Date:     %s', datestr(now));
    '';
    '───────────────────────────────────────────────────────────────';
    '  PLAN SUMMARY';
    '───────────────────────────────────────────────────────────────';
};

if ~isempty(segmentInfo)
    reportLines{end+1} = sprintf('  Number of Beams:       %d', segmentInfo.numBeams);
    reportLines{end+1} = sprintf('  Total Segments:        %d', segmentInfo.totalSegments);
    reportLines{end+1} = sprintf('  Segments Calculated:   %d', numSegments);
    
    totalMU = 0;
    for i = 1:length(segmentInfo.beams)
        totalMU = totalMU + segmentInfo.beams{i}.beamMeterset;
    end
    reportLines{end+1} = sprintf('  Total Meterset:        %.1f MU', totalMU);
end

reportLines{end+1} = '';
reportLines{end+1} = '───────────────────────────────────────────────────────────────';
reportLines{end+1} = '  DOSE STATISTICS';
reportLines{end+1} = '───────────────────────────────────────────────────────────────';
reportLines{end+1} = sprintf('  Calculated Max Dose:   %.4f Gy', max(totalDose(:)));
reportLines{end+1} = sprintf('  Calculated Mean Dose:  %.4f Gy', mean(totalDose(totalDose > 0)));

if ~isempty(referenceDose)
    reportLines{end+1} = sprintf('  Reference Max Dose:    %.4f Gy', max(referenceDose(:)));
    reportLines{end+1} = sprintf('  Reference Mean Dose:   %.4f Gy', mean(referenceDose(referenceDose > 0)));
    
    reportLines{end+1} = '';
    reportLines{end+1} = '───────────────────────────────────────────────────────────────';
    reportLines{end+1} = '  COMPARISON METRICS';
    reportLines{end+1} = '───────────────────────────────────────────────────────────────';
    
    if ~isempty(comparison)
        reportLines{end+1} = sprintf('  Mean Abs Difference:   %.4f Gy', comparison.metrics.meanAbsDiff);
        reportLines{end+1} = sprintf('  Max Difference:        %.4f Gy', comparison.metrics.maxDiff);
        reportLines{end+1} = sprintf('  RMS Difference:        %.4f Gy', comparison.metrics.rmsDiff);
        if isfield(comparison.metrics, 'meanRelativeDiff')
            reportLines{end+1} = sprintf('  Mean Rel Diff (>50%%):  %.2f%%', comparison.metrics.meanRelativeDiff);
        end
    end
    
    if exist('passRate', 'var')
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('  Gamma Pass Rate:       %.1f%% (%d%%/%dmm)', ...
            passRate, gammaDoseThreshold, gammaDTAThreshold);
    end
end

reportLines{end+1} = '';
reportLines{end+1} = '───────────────────────────────────────────────────────────────';
reportLines{end+1} = '  GRID INFORMATION';
reportLines{end+1} = '───────────────────────────────────────────────────────────────';
reportLines{end+1} = sprintf('  Grid Dimensions:       %d x %d x %d', gridSize(1), gridSize(2), gridSize(3));
reportLines{end+1} = sprintf('  Voxel Size:            %.2f x %.2f x %.2f mm', dx, dy, dz);

reportLines{end+1} = '';
reportLines{end+1} = '═══════════════════════════════════════════════════════════════';

% Display report
text(0.05, 0.98, reportLines, 'FontSize', 10, 'FontName', 'FixedWidth', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

if saveAllFigures
    saveFigure(fig8, fullfile(outputDir, '08_SummaryReport'), figureFormat, dpi);
end

% Also save report as text file
reportFile = fullfile(outputDir, 'VerificationReport.txt');
fid = fopen(reportFile, 'w');
for i = 1:length(reportLines)
    fprintf(fid, '%s\n', reportLines{i});
end
fclose(fid);
fprintf('  Report saved to: %s\n', reportFile);

%% Finish
fprintf('\n=======================================================\n');
fprintf('  Visualization Complete\n');
fprintf('=======================================================\n');
fprintf('  Output directory: %s\n', outputDir);
fprintf('  Figures saved: %d\n', 8);
fprintf('=======================================================\n\n');

%% Helper Functions

function [gammaMap, passRate] = calculateGamma2D(ref, calc, maxDose, doseCrit, dtaCrit, doseCutoff)
    % Simplified 2D gamma calculation
    % doseCrit: dose difference criterion as fraction (e.g., 0.03 for 3%)
    % dtaCrit: distance-to-agreement criterion in voxels
    % doseCutoff: minimum dose threshold as fraction of max
    
    [ny, nx] = size(ref);
    gammaMap = inf(ny, nx);
    
    % Create dose threshold mask
    doseMask = ref >= (doseCutoff * maxDose);
    
    % Search radius in voxels
    searchRadius = ceil(dtaCrit * 2);
    
    for i = 1:ny
        for j = 1:nx
            if ~doseMask(i, j)
                gammaMap(i, j) = NaN;
                continue;
            end
            
            refDose = ref(i, j);
            calcDose = calc(i, j);
            
            minGamma = inf;
            
            % Search neighboring points
            for di = -searchRadius:searchRadius
                for dj = -searchRadius:searchRadius
                    ni = i + di;
                    nj = j + dj;
                    
                    if ni < 1 || ni > ny || nj < 1 || nj > nx
                        continue;
                    end
                    
                    % Distance in voxels
                    dist = sqrt(di^2 + dj^2);
                    
                    % Dose difference
                    doseDiff = abs(calc(ni, nj) - refDose) / (doseCrit * maxDose);
                    
                    % Gamma
                    gamma = sqrt((dist/dtaCrit)^2 + doseDiff^2);
                    
                    if gamma < minGamma
                        minGamma = gamma;
                    end
                end
            end
            
            gammaMap(i, j) = minGamma;
        end
    end
    
    % Calculate pass rate
    validGamma = gammaMap(doseMask);
    passRate = 100 * sum(validGamma <= 1) / length(validGamma);
end

function saveFigure(fig, filename, format, dpi)
    % Save figure in specified format(s)
    switch format
        case 'png'
            print(fig, filename, '-dpng', sprintf('-r%d', dpi));
        case 'fig'
            savefig(fig, [filename, '.fig']);
        case 'both'
            print(fig, filename, '-dpng', sprintf('-r%d', dpi));
            savefig(fig, [filename, '.fig']);
    end
end
