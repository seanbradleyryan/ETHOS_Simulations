%% ETHOS Segment Dose Verification and Visualization
% Purpose: Create comprehensive visualizations of segment dose calculations
%          from calculate_field_dose.m output
%
% Prerequisites:
%   - calculate_field_dose.m must have been run successfully
%   - Output files must exist in the SegmentDoses directory
%
% Outputs:
%   - Multi-panel dose comparison figures
%   - Dose profiles (axial, coronal, sagittal)
%   - Beam-by-beam analysis
%   - Gamma analysis (if applicable)
%   - Summary statistics
%
% Author: Generated for ETHOS dose analysis
% Date: 2025

clear; clc; close all;

%% ==================== CONFIGURATION ====================
% Patient and session to analyze
patientID = '1194203';
sessionName = 'Session_1';

% Base directory (must match calculate_field_dose.m)
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';

% Output directory for figures
figureOutputDir = fullfile(wd, 'Figures', patientID, sessionName);

% Gamma analysis parameters
gammaDistCriteria = 3;    % mm (distance-to-agreement)
gammaDoseCriteria = 3;    % % (dose difference)
gammaThreshold = 10;      % % of max dose (exclude low dose regions)

% Plot settings
colormap_dose = 'jet';
figureVisible = 'on';       % 'on' or 'off' (off for batch processing)

%% ==================== INITIALIZATION ====================
fprintf('==========================================================\n');
fprintf('  ETHOS Segment Dose Verification\n');
fprintf('  Patient: %s, Session: %s\n', patientID, sessionName);
fprintf('==========================================================\n\n');

% Define paths
dataPath = fullfile(wd, 'SegmentDoses', patientID, sessionName);

% Create figure output directory
if ~exist(figureOutputDir, 'dir')
    mkdir(figureOutputDir);
    fprintf('Created figure output directory: %s\n', figureOutputDir);
end

%% ==================== LOAD DATA ====================
fprintf('[1/6] Loading data...\n');

% Check for required files
requiredFiles = {'segmentDoses.mat', 'segmentData.mat'};
for i = 1:length(requiredFiles)
    if ~exist(fullfile(dataPath, requiredFiles{i}), 'file')
        error('Required file not found: %s\nRun calculate_field_dose.m first.', requiredFiles{i});
    end
end

% Load main results
fprintf('  Loading segmentDoses.mat...\n');
load(fullfile(dataPath, 'segmentDoses.mat'));
fprintf('  Loading segmentData.mat...\n');
load(fullfile(dataPath, 'segmentData.mat'));

% Load comparison if available
if exist(fullfile(dataPath, 'doseComparison.mat'), 'file')
    fprintf('  Loading doseComparison.mat...\n');
    load(fullfile(dataPath, 'doseComparison.mat'));
    hasComparison = true;
else
    hasComparison = false;
    fprintf('  No doseComparison.mat found (reference dose comparison unavailable)\n');
end

% Load CT if available
if exist(fullfile(dataPath, 'ctResampled.mat'), 'file')
    fprintf('  Loading ctResampled.mat...\n');
    load(fullfile(dataPath, 'ctResampled.mat'));
    hasCT = true;
else
    hasCT = false;
    fprintf('  No ctResampled.mat found\n');
end

fprintf('  Data loaded successfully\n\n');

% Extract key variables
calcDose = totalDoseResampled;
if exist('referenceDose', 'var') && ~isempty(referenceDose)
    refDose = referenceDose;
    hasRefDose = true;
else
    hasRefDose = false;
end

doseSize = size(calcDose);
fprintf('  Dose grid size: [%d, %d, %d]\n', doseSize(1), doseSize(2), doseSize(3));
fprintf('  Dose grid resolution: [%.2f, %.2f, %.2f] mm\n', ...
    doseGrid.resolution(1), doseGrid.resolution(2), doseGrid.resolution(3));
fprintf('  Calculated max dose: %.4f Gy\n', max(calcDose(:)));
if hasRefDose
    fprintf('  Reference max dose: %.4f Gy\n', max(refDose(:)));
end
fprintf('  Number of beams: %d\n', segmentData.numBeams);
fprintf('  Total segments: %d\n', segmentData.totalSegments);

%% ==================== FIGURE 1: DOSE OVERVIEW ====================
fprintf('\n[2/6] Creating dose overview figure...\n');

% Find slice with maximum dose
[maxDose, maxIdx] = max(calcDose(:));
[maxI, maxJ, maxK] = ind2sub(doseSize, maxIdx);

fig1 = figure('Name', 'Dose Overview', 'Position', [50, 50, 1600, 900], ...
    'Visible', figureVisible, 'Color', 'w');

% Create coordinate vectors for plotting (in mm)
xVec = (0:doseSize(2)-1) * doseGrid.resolution(2);
yVec = (0:doseSize(1)-1) * doseGrid.resolution(1);
zVec = (0:doseSize(3)-1) * doseGrid.resolution(3);

if hasRefDose
    % 3x3 layout: Calc, Ref, Diff for each plane
    
    % Row 1: Axial slices (XY plane at max dose Z)
    subplot(3, 3, 1);
    imagesc(xVec, yVec, squeeze(calcDose(:, :, maxK)));
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Calculated - Axial (Z=%d)', maxK));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    subplot(3, 3, 2);
    imagesc(xVec, yVec, squeeze(refDose(:, :, maxK)));
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Reference - Axial (Z=%d)', maxK));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    subplot(3, 3, 3);
    diffSlice = squeeze(calcDose(:, :, maxK) - refDose(:, :, maxK));
    maxDiff = max(abs(diffSlice(:)));
    imagesc(xVec, yVec, diffSlice, [-maxDiff, maxDiff]);
    axis image; colorbar;
    title('Difference (Calc - Ref)');
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % Row 2: Coronal slices (XZ plane at max dose Y)
    subplot(3, 3, 4);
    imagesc(xVec, zVec, squeeze(calcDose(maxI, :, :))');
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Calculated - Coronal (Y=%d)', maxI));
    xlabel('X (mm)'); ylabel('Z (mm)');
    
    subplot(3, 3, 5);
    imagesc(xVec, zVec, squeeze(refDose(maxI, :, :))');
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Reference - Coronal (Y=%d)', maxI));
    xlabel('X (mm)'); ylabel('Z (mm)');
    
    subplot(3, 3, 6);
    diffSlice = squeeze(calcDose(maxI, :, :) - refDose(maxI, :, :))';
    maxDiff = max(abs(diffSlice(:)));
    imagesc(xVec, zVec, diffSlice, [-maxDiff, maxDiff]);
    axis image; colorbar;
    title('Difference (Calc - Ref)');
    xlabel('X (mm)'); ylabel('Z (mm)');
    
    % Row 3: Sagittal slices (YZ plane at max dose X)
    subplot(3, 3, 7);
    imagesc(yVec, zVec, squeeze(calcDose(:, maxJ, :))');
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Calculated - Sagittal (X=%d)', maxJ));
    xlabel('Y (mm)'); ylabel('Z (mm)');
    
    subplot(3, 3, 8);
    imagesc(yVec, zVec, squeeze(refDose(:, maxJ, :))');
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Reference - Sagittal (X=%d)', maxJ));
    xlabel('Y (mm)'); ylabel('Z (mm)');
    
    subplot(3, 3, 9);
    diffSlice = squeeze(calcDose(:, maxJ, :) - refDose(:, maxJ, :))';
    maxDiff = max(abs(diffSlice(:)));
    imagesc(yVec, zVec, diffSlice, [-maxDiff, maxDiff]);
    axis image; colorbar;
    title('Difference (Calc - Ref)');
    xlabel('Y (mm)'); ylabel('Z (mm)');
    
else
    % 1x3 layout: Just calculated dose
    subplot(1, 3, 1);
    imagesc(xVec, yVec, squeeze(calcDose(:, :, maxK)));
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Axial (Z=%d)', maxK));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    subplot(1, 3, 2);
    imagesc(xVec, zVec, squeeze(calcDose(maxI, :, :))');
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Coronal (Y=%d)', maxI));
    xlabel('X (mm)'); ylabel('Z (mm)');
    
    subplot(1, 3, 3);
    imagesc(yVec, zVec, squeeze(calcDose(:, maxJ, :))');
    axis image; colorbar; colormap(gca, colormap_dose);
    title(sprintf('Sagittal (X=%d)', maxJ));
    xlabel('Y (mm)'); ylabel('Z (mm)');
end

sgtitle(sprintf('Dose Distribution Overview - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(fig1, fullfile(figureOutputDir, 'dose_overview.png'));
saveas(fig1, fullfile(figureOutputDir, 'dose_overview.fig'));
fprintf('  Saved: dose_overview.png\n');

%% ==================== FIGURE 2: DOSE PROFILES ====================
fprintf('\n[3/6] Creating dose profile figure...\n');

fig2 = figure('Name', 'Dose Profiles', 'Position', [100, 100, 1400, 500], ...
    'Visible', figureVisible, 'Color', 'w');

% Profile through max dose point in each direction
% X profile (lateral)
subplot(1, 3, 1);
xProfile_calc = squeeze(calcDose(maxI, :, maxK));
plot(xVec, xProfile_calc, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if hasRefDose
    xProfile_ref = squeeze(refDose(maxI, :, maxK));
    plot(xVec, xProfile_ref, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('X Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('Lateral Profile (Y=%d, Z=%d)', maxI, maxK));
legend('Location', 'best');
grid on;

% Y profile (AP)
subplot(1, 3, 2);
yProfile_calc = squeeze(calcDose(:, maxJ, maxK));
plot(yVec, yProfile_calc, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if hasRefDose
    yProfile_ref = squeeze(refDose(:, maxJ, maxK));
    plot(yVec, yProfile_ref, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('Y Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('AP Profile (X=%d, Z=%d)', maxJ, maxK));
legend('Location', 'best');
grid on;

% Z profile (SI / depth)
subplot(1, 3, 3);
zProfile_calc = squeeze(calcDose(maxI, maxJ, :));
plot(zVec, zProfile_calc, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
hold on;
if hasRefDose
    zProfile_ref = squeeze(refDose(maxI, maxJ, :));
    plot(zVec, zProfile_ref, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
end
xlabel('Z Position (mm)');
ylabel('Dose (Gy)');
title(sprintf('SI Profile (X=%d, Y=%d)', maxJ, maxI));
legend('Location', 'best');
grid on;

sgtitle(sprintf('Dose Profiles Through Maximum - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig2, fullfile(figureOutputDir, 'dose_profiles.png'));
saveas(fig2, fullfile(figureOutputDir, 'dose_profiles.fig'));
fprintf('  Saved: dose_profiles.png\n');

%% ==================== FIGURE 3: BEAM CONTRIBUTIONS ====================
fprintf('\n[4/6] Creating beam contribution figure...\n');

numBeams = length(beamDosesResampled);
validBeams = find(~cellfun(@isempty, beamDosesResampled));
numValidBeams = length(validBeams);

% Determine subplot layout
nCols = min(4, numValidBeams);
nRows = ceil(numValidBeams / nCols);

fig3 = figure('Name', 'Beam Contributions', 'Position', [100, 100, 400*nCols, 350*nRows], ...
    'Visible', figureVisible, 'Color', 'w');

for idx = 1:numValidBeams
    beamIdx = validBeams(idx);
    beamData = beamDosesResampled{beamIdx};
    
    subplot(nRows, nCols, idx);
    
    % Show axial slice at max dose location
    beamSlice = squeeze(beamData.physicalDose(:, :, maxK));
    imagesc(xVec, yVec, beamSlice);
    axis image; colorbar; colormap(gca, colormap_dose);
    
    title(sprintf('Beam %d: %s\nGantry=%.0f°, Max=%.3f Gy', ...
        beamIdx, beamData.beamName, beamData.gantryAngle, beamData.maxDose));
    xlabel('X (mm)'); ylabel('Y (mm)');
end

sgtitle(sprintf('Individual Beam Contributions (Axial Z=%d) - Patient %s', maxK, patientID), ...
    'FontSize', 14, 'FontWeight', 'bold');

saveas(fig3, fullfile(figureOutputDir, 'beam_contributions.png'));
saveas(fig3, fullfile(figureOutputDir, 'beam_contributions.fig'));
fprintf('  Saved: beam_contributions.png\n');

%% ==================== FIGURE 4: BEAM STATISTICS ====================
fprintf('\n[5/6] Creating beam statistics figure...\n');

fig4 = figure('Name', 'Beam Statistics', 'Position', [100, 100, 1200, 800], ...
    'Visible', figureVisible, 'Color', 'w');

% Collect beam statistics
beamNames = cell(numValidBeams, 1);
beamMetersets = zeros(numValidBeams, 1);
beamMaxDoses = zeros(numValidBeams, 1);
beamGantryAngles = zeros(numValidBeams, 1);
beamNumSegments = zeros(numValidBeams, 1);

for idx = 1:numValidBeams
    beamIdx = validBeams(idx);
    beamData = beamDosesResampled{beamIdx};
    segData = segmentData.beams{beamIdx};
    
    beamNames{idx} = beamData.beamName;
    beamMetersets(idx) = beamData.beamMeterset;
    beamMaxDoses(idx) = beamData.maxDose;
    beamGantryAngles(idx) = beamData.gantryAngle;
    beamNumSegments(idx) = length(segData.segments);
end

% Subplot 1: Beam metersets (MU)
subplot(2, 2, 1);
bar(beamMetersets);
xlabel('Beam Index');
ylabel('Meterset (MU)');
title('Beam Metersets');
xticks(1:numValidBeams);
xticklabels(beamNames);
xtickangle(45);
grid on;

% Subplot 2: Max dose per beam
subplot(2, 2, 2);
bar(beamMaxDoses);
xlabel('Beam Index');
ylabel('Max Dose (Gy)');
title('Maximum Dose per Beam');
xticks(1:numValidBeams);
xticklabels(beamNames);
xtickangle(45);
grid on;

% Subplot 3: Gantry angles (polar plot)
subplot(2, 2, 3);
theta = deg2rad(beamGantryAngles);
rho = beamMetersets / max(beamMetersets);  % Normalize for visualization
polarplot(theta, rho, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
title('Beam Gantry Angles (radius = relative MU)');
thetalim([0 360]);

% Subplot 4: Segments per beam
subplot(2, 2, 4);
bar(beamNumSegments);
xlabel('Beam Index');
ylabel('Number of Segments');
title('Segments per Beam');
xticks(1:numValidBeams);
xticklabels(beamNames);
xtickangle(45);
grid on;

sgtitle(sprintf('Beam Statistics - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig4, fullfile(figureOutputDir, 'beam_statistics.png'));
saveas(fig4, fullfile(figureOutputDir, 'beam_statistics.fig'));
fprintf('  Saved: beam_statistics.png\n');

%% ==================== FIGURE 5: DOSE DIFFERENCE ANALYSIS ====================
if hasRefDose
    fprintf('\n[6/6] Creating dose difference analysis figure...\n');
    
    fig5 = figure('Name', 'Dose Difference Analysis', 'Position', [100, 100, 1400, 900], ...
        'Visible', figureVisible, 'Color', 'w');
    
    % Calculate difference
    doseDiff = calcDose - refDose;
    
    % Calculate relative difference (where reference > threshold)
    doseThreshold = 0.1 * max(refDose(:));  % 10% of max
    validMask = refDose > doseThreshold;
    relativeDiff = zeros(size(doseDiff));
    relativeDiff(validMask) = (doseDiff(validMask) ./ refDose(validMask)) * 100;
    
    % Subplot 1: Absolute difference histogram
    subplot(2, 3, 1);
    histogram(doseDiff(:), 100, 'FaceColor', 'b', 'EdgeColor', 'none');
    xlabel('Dose Difference (Gy)');
    ylabel('Voxel Count');
    title('Absolute Difference Distribution');
    xline(0, 'r--', 'LineWidth', 2);
    grid on;
    
    % Subplot 2: Relative difference histogram (high dose region only)
    subplot(2, 3, 2);
    histogram(relativeDiff(validMask), 100, 'FaceColor', 'g', 'EdgeColor', 'none');
    xlabel('Relative Difference (%)');
    ylabel('Voxel Count');
    title(sprintf('Relative Difference (>%.0f%% max dose)', 10));
    xline(0, 'r--', 'LineWidth', 2);
    grid on;
    
    % Subplot 3: Scatter plot - Calculated vs Reference
    subplot(2, 3, 3);
    % Subsample for performance
    subsampleFactor = max(1, round(numel(calcDose) / 50000));
    calcSub = calcDose(1:subsampleFactor:end);
    refSub = refDose(1:subsampleFactor:end);
    scatter(refSub(:), calcSub(:), 1, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    maxVal = max([max(calcSub(:)), max(refSub(:))]);
    plot([0, maxVal], [0, maxVal], 'r-', 'LineWidth', 2);
    xlabel('Reference Dose (Gy)');
    ylabel('Calculated Dose (Gy)');
    title('Calculated vs Reference');
    axis equal;
    xlim([0, maxVal]);
    ylim([0, maxVal]);
    grid on;
    
    % Calculate correlation
    validIdx = refSub > 0.01;
    if any(validIdx)
        R = corrcoef(refSub(validIdx), calcSub(validIdx));
        text(0.1*maxVal, 0.9*maxVal, sprintf('R = %.4f', R(1,2)), 'FontSize', 12);
    end
    
    % Subplot 4: Absolute difference map (axial)
    subplot(2, 3, 4);
    diffSlice = squeeze(doseDiff(:, :, maxK));
    maxAbsDiff = max(abs(diffSlice(:)));
    imagesc(xVec, yVec, diffSlice, [-maxAbsDiff, maxAbsDiff]);
    axis image; colorbar;
    title(sprintf('Absolute Diff - Axial (Z=%d)', maxK));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % Subplot 5: Relative difference map (axial)
    subplot(2, 3, 5);
    relDiffSlice = squeeze(relativeDiff(:, :, maxK));
    maxRelDiff = min(50, max(abs(relDiffSlice(:))));  % Cap at 50%
    imagesc(xVec, yVec, relDiffSlice, [-maxRelDiff, maxRelDiff]);
    axis image; colorbar;
    title(sprintf('Relative Diff (%%) - Axial (Z=%d)', maxK));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % Subplot 6: Statistics text box
    subplot(2, 3, 6);
    axis off;
    
    % Calculate statistics
    stats = struct();
    stats.meanAbsDiff = mean(abs(doseDiff(:)));
    stats.maxAbsDiff = max(abs(doseDiff(:)));
    stats.stdAbsDiff = std(doseDiff(:));
    stats.rmsDiff = sqrt(mean(doseDiff(:).^2));
    
    if any(validMask(:))
        stats.meanRelDiff = mean(abs(relativeDiff(validMask)));
        stats.maxRelDiff = max(abs(relativeDiff(validMask)));
    else
        stats.meanRelDiff = NaN;
        stats.maxRelDiff = NaN;
    end
    
    % Calculate percentage of voxels within tolerance
    tolerance3pct = sum(abs(relativeDiff(validMask)) <= 3) / sum(validMask(:)) * 100;
    tolerance5pct = sum(abs(relativeDiff(validMask)) <= 5) / sum(validMask(:)) * 100;
    
    statsText = {
        sprintf('DOSE DIFFERENCE STATISTICS')
        sprintf('================================')
        sprintf('')
        sprintf('Absolute Differences:')
        sprintf('  Mean: %.4f Gy', stats.meanAbsDiff)
        sprintf('  Max:  %.4f Gy', stats.maxAbsDiff)
        sprintf('  Std:  %.4f Gy', stats.stdAbsDiff)
        sprintf('  RMS:  %.4f Gy', stats.rmsDiff)
        sprintf('')
        sprintf('Relative Differences (>10%% max):')
        sprintf('  Mean: %.2f %%', stats.meanRelDiff)
        sprintf('  Max:  %.2f %%', stats.maxRelDiff)
        sprintf('')
        sprintf('Tolerance Analysis:')
        sprintf('  Within 3%%: %.1f %%', tolerance3pct)
        sprintf('  Within 5%%: %.1f %%', tolerance5pct)
    };
    
    text(0.1, 0.9, statsText, 'FontSize', 11, 'FontName', 'FixedWidth', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    sgtitle(sprintf('Dose Difference Analysis - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');
    
    saveas(fig5, fullfile(figureOutputDir, 'dose_difference_analysis.png'));
    saveas(fig5, fullfile(figureOutputDir, 'dose_difference_analysis.fig'));
    fprintf('  Saved: dose_difference_analysis.png\n');
    
else
    fprintf('\n[6/6] Skipping dose difference analysis (no reference dose)\n');
    stats = struct();
    tolerance3pct = NaN;
    tolerance5pct = NaN;
end

%% ==================== FIGURE 6: CT WITH DOSE OVERLAY ====================
if hasCT && isfield(ctResampled_struct, 'cubeHU') && ~isempty(ctResampled_struct.cubeHU{1})
    fprintf('\n[Bonus] Creating CT with dose overlay figure...\n');
    
    fig6 = figure('Name', 'CT with Dose Overlay', 'Position', [100, 100, 1200, 400], ...
        'Visible', figureVisible, 'Color', 'w');
    
    ctData = ctResampled_struct.cubeHU{1};
    
    % Normalize dose for overlay
    doseNorm = calcDose / max(calcDose(:));
    
    % Axial
    subplot(1, 3, 1);
    ctSlice = squeeze(ctData(:, :, maxK));
    imagesc(xVec, yVec, ctSlice, [-200, 200]);
    colormap(gca, 'gray');
    hold on;
    doseSlice = squeeze(doseNorm(:, :, maxK));
    h = imagesc(xVec, yVec, doseSlice);
    set(h, 'AlphaData', doseSlice * 0.6);
    axis image;
    title(sprintf('Axial (Z=%d)', maxK));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % Coronal
    subplot(1, 3, 2);
    ctSlice = squeeze(ctData(maxI, :, :))';
    imagesc(xVec, zVec, ctSlice, [-200, 200]);
    colormap(gca, 'gray');
    hold on;
    doseSlice = squeeze(doseNorm(maxI, :, :))';
    h = imagesc(xVec, zVec, doseSlice);
    set(h, 'AlphaData', doseSlice * 0.6);
    axis image;
    title(sprintf('Coronal (Y=%d)', maxI));
    xlabel('X (mm)'); ylabel('Z (mm)');
    
    % Sagittal
    subplot(1, 3, 3);
    ctSlice = squeeze(ctData(:, maxJ, :))';
    imagesc(yVec, zVec, ctSlice, [-200, 200]);
    colormap(gca, 'gray');
    hold on;
    doseSlice = squeeze(doseNorm(:, maxJ, :))';
    h = imagesc(yVec, zVec, doseSlice);
    set(h, 'AlphaData', doseSlice * 0.6);
    axis image;
    title(sprintf('Sagittal (X=%d)', maxJ));
    xlabel('Y (mm)'); ylabel('Z (mm)');
    
    sgtitle(sprintf('CT with Dose Overlay - Patient %s', patientID), 'FontSize', 14, 'FontWeight', 'bold');
    
    saveas(fig6, fullfile(figureOutputDir, 'ct_dose_overlay.png'));
    saveas(fig6, fullfile(figureOutputDir, 'ct_dose_overlay.fig'));
    fprintf('  Saved: ct_dose_overlay.png\n');
end

%% ==================== FIGURE 7: DVH COMPARISON ====================
if hasRefDose
    fprintf('\n[Bonus] Creating DVH comparison figure...\n');
    
    fig7 = figure('Name', 'DVH Comparison', 'Position', [100, 100, 800, 600], ...
        'Visible', figureVisible, 'Color', 'w');
    
    % Calculate cumulative DVH for whole volume
    maxDoseVal = max([max(calcDose(:)), max(refDose(:))]);
    doseBins = linspace(0, maxDoseVal, 100);
    
    dvh_calc = zeros(size(doseBins));
    dvh_ref = zeros(size(doseBins));
    totalVoxels = numel(calcDose);
    
    for i = 1:length(doseBins)
        dvh_calc(i) = sum(calcDose(:) >= doseBins(i)) / totalVoxels * 100;
        dvh_ref(i) = sum(refDose(:) >= doseBins(i)) / totalVoxels * 100;
    end
    
    plot(doseBins, dvh_calc, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
    hold on;
    plot(doseBins, dvh_ref, 'r--', 'LineWidth', 2, 'DisplayName', 'Reference');
    
    xlabel('Dose (Gy)');
    ylabel('Volume (%)');
    title(sprintf('Dose-Volume Histogram - Patient %s', patientID));
    legend('Location', 'best');
    grid on;
    xlim([0, maxDoseVal]);
    ylim([0, 100]);
    
    saveas(fig7, fullfile(figureOutputDir, 'dvh_comparison.png'));
    saveas(fig7, fullfile(figureOutputDir, 'dvh_comparison.fig'));
    fprintf('  Saved: dvh_comparison.png\n');
end

%% ==================== SUMMARY REPORT ====================
fprintf('\n==========================================================\n');
fprintf('  VERIFICATION COMPLETE\n');
fprintf('==========================================================\n\n');

fprintf('Patient: %s\n', patientID);
fprintf('Session: %s\n', sessionName);
fprintf('\n');

fprintf('DOSE SUMMARY:\n');
fprintf('  Calculated max dose: %.4f Gy\n', max(calcDose(:)));
if hasRefDose
    fprintf('  Reference max dose:  %.4f Gy\n', max(refDose(:)));
    fprintf('  Max absolute diff:   %.4f Gy\n', stats.maxAbsDiff);
    fprintf('  Mean absolute diff:  %.4f Gy\n', stats.meanAbsDiff);
    fprintf('  RMS difference:      %.4f Gy\n', stats.rmsDiff);
end
fprintf('\n');

fprintf('BEAM SUMMARY:\n');
fprintf('  Total beams: %d\n', numValidBeams);
fprintf('  Total segments: %d\n', sum(beamNumSegments));
fprintf('  Total MU: %.2f\n', sum(beamMetersets));
fprintf('\n');

fprintf('Figures saved to: %s\n', figureOutputDir);
fprintf('\n');

% Save summary to text file
summaryFile = fullfile(figureOutputDir, 'verification_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'ETHOS Segment Dose Verification Summary\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Patient: %s\n', patientID);
fprintf(fid, 'Session: %s\n', sessionName);
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, 'DOSE SUMMARY:\n');
fprintf(fid, '  Calculated max dose: %.4f Gy\n', max(calcDose(:)));
if hasRefDose
    fprintf(fid, '  Reference max dose:  %.4f Gy\n', max(refDose(:)));
    fprintf(fid, '  Max absolute diff:   %.4f Gy\n', stats.maxAbsDiff);
    fprintf(fid, '  Mean absolute diff:  %.4f Gy\n', stats.meanAbsDiff);
    fprintf(fid, '  RMS difference:      %.4f Gy\n', stats.rmsDiff);
    fprintf(fid, '  Within 3%%: %.1f %%\n', tolerance3pct);
    fprintf(fid, '  Within 5%%: %.1f %%\n', tolerance5pct);
end
fprintf(fid, '\nBEAM SUMMARY:\n');
for idx = 1:numValidBeams
    beamIdx = validBeams(idx);
    fprintf(fid, '  Beam %d: %s\n', beamIdx, beamNames{idx});
    fprintf(fid, '    Gantry: %.1f deg\n', beamGantryAngles(idx));
    fprintf(fid, '    MU: %.2f\n', beamMetersets(idx));
    fprintf(fid, '    Segments: %d\n', beamNumSegments(idx));
    fprintf(fid, '    Max dose: %.4f Gy\n', beamMaxDoses(idx));
end
fclose(fid);
fprintf('Summary saved to: %s\n', summaryFile);

fprintf('\n==========================================================\n');