%% RTDOSE File Inspector and Visualizer
% Purpose: Load and inspect RTDOSE files from ETHOS exports
% Displays detailed information and cross-sectional views
% Author: Generated for ETHOS dose analysis
% Date: 2025

clear; clc; close all;

%% Configuration
id = '1194203';
session = 'Session 1';

% Base directory
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
rawwd = fullfile(wd, 'EthosExports', id, 'Pancreas', session);
dicomPath = fullfile(rawwd, 'sct');

fprintf('========================================\n');
fprintf('RTDOSE File Inspector\n');
fprintf('========================================\n');
fprintf('Patient: %s, Session: %s\n', id, session);
fprintf('DICOM Path: %s\n\n', dicomPath);

%% Find and Load RTDOSE file
fprintf('[1/3] Searching for RTDOSE file...\n');

% Look for RTDOSE files
rtdoseFile = dir(fullfile(dicomPath, 'RD*.dcm'));
if isempty(rtdoseFile)
    rtdoseFile = dir(fullfile(dicomPath, '*RTDOSE*.dcm'));
end
if isempty(rtdoseFile)
    rtdoseFile = dir(fullfile(dicomPath, '*.dcm'));
    % Filter for RTDOSE by checking file headers
    for i = length(rtdoseFile):-1:1
        try
            info = dicominfo(fullfile(rtdoseFile(i).folder, rtdoseFile(i).name));
            if ~strcmp(info.Modality, 'RTDOSE')
                rtdoseFile(i) = [];
            end
        catch
            rtdoseFile(i) = [];
        end
    end
end

if isempty(rtdoseFile)
    error('No RTDOSE file found in directory: %s', dicomPath);
end

fprintf('  Found %d RTDOSE file(s):\n', length(rtdoseFile));
for i = 1:length(rtdoseFile)
    fprintf('    %d. %s\n', i, rtdoseFile(i).name);
end

% Use first RTDOSE file
rtdoseFilename = fullfile(rtdoseFile(1).folder, rtdoseFile(1).name);
fprintf('\n  Loading: %s\n', rtdoseFile(1).name);

%% Read RTDOSE Header Information
fprintf('\n[2/3] Reading RTDOSE header information...\n');

rtdoseInfo = dicominfo(rtdoseFilename);

fprintf('\n  === DICOM Header Information ===\n');
fprintf('  Patient ID: %s\n', rtdoseInfo.PatientID);
if isfield(rtdoseInfo, 'PatientName')
    fprintf('  Patient Name: %s\n', rtdoseInfo.PatientName.FamilyName);
end
fprintf('  Modality: %s\n', rtdoseInfo.Modality);
fprintf('  Manufacturer: %s\n', rtdoseInfo.Manufacturer);
if isfield(rtdoseInfo, 'ManufacturerModelName')
    fprintf('  Model: %s\n', rtdoseInfo.ManufacturerModelName);
end

fprintf('\n  === Dose Grid Information ===\n');
fprintf('  Grid Dimensions: %d x %d x %d\n', rtdoseInfo.Rows, rtdoseInfo.Columns, ...
    size(dicomread(rtdoseFilename), 3));
fprintf('  Pixel Spacing: [%.3f, %.3f] mm\n', rtdoseInfo.PixelSpacing(1), rtdoseInfo.PixelSpacing(2));

% Slice thickness calculation
if isfield(rtdoseInfo, 'GridFrameOffsetVector') && length(rtdoseInfo.GridFrameOffsetVector) > 1
    sliceSpacing = rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1);
    fprintf('  Slice Spacing: %.3f mm\n', sliceSpacing);
elseif isfield(rtdoseInfo, 'SliceThickness')
    fprintf('  Slice Thickness: %.3f mm\n', rtdoseInfo.SliceThickness);
else
    fprintf('  Slice Spacing: Not specified in header\n');
end

fprintf('  Dose Grid Scaling: %.6e\n', rtdoseInfo.DoseGridScaling);
fprintf('  Dose Units: %s\n', rtdoseInfo.DoseUnits);
fprintf('  Dose Type: %s\n', rtdoseInfo.DoseType);
fprintf('  Dose Summation Type: %s\n', rtdoseInfo.DoseSummationType);

if isfield(rtdoseInfo, 'ImagePositionPatient')
    fprintf('  Image Position: [%.2f, %.2f, %.2f] mm\n', ...
        rtdoseInfo.ImagePositionPatient(1), ...
        rtdoseInfo.ImagePositionPatient(2), ...
        rtdoseInfo.ImagePositionPatient(3));
end

%% Load and Scale Dose Data
fprintf('\n  === Loading Dose Data ===\n');

% Read raw dose data
doseRaw = dicomread(rtdoseFilename);
fprintf('  Raw data size: %d x %d x %d\n', size(doseRaw, 1), size(doseRaw, 2), size(doseRaw, 3));
fprintf('  Raw data class: %s\n', class(doseRaw));

% Apply scaling to get dose in Gy
dose = double(doseRaw) * rtdoseInfo.DoseGridScaling;

fprintf('\n  === Dose Statistics ===\n');
fprintf('  Minimum dose: %.4f Gy\n', min(dose(:)));
fprintf('  Maximum dose: %.4f Gy\n', max(dose(:)));
fprintf('  Mean dose: %.4f Gy\n', mean(dose(:)));
fprintf('  Median dose: %.4f Gy\n', median(dose(:)));
fprintf('  Std deviation: %.4f Gy\n', std(dose(:)));

% Find location of maximum dose
[maxDose, maxIdx] = max(dose(:));
[maxI, maxJ, maxK] = ind2sub(size(dose), maxIdx);
fprintf('  Max dose location: [%d, %d, %d]\n', maxI, maxJ, maxK);

%% Visualize Dose Distribution
fprintf('\n[3/3] Generating visualizations...\n');

% Create figure with multiple views
figure('Name', 'RTDOSE Cross-Sectional Views', 'Position', [100, 100, 1400, 500]);

% Find slice with maximum dose
centralSlice = maxK;

% Axial view (at max dose location)
subplot(1, 3, 1);
imagesc(dose(:, :, centralSlice)');
axis equal tight;
colorbar;
colormap(jet);
title(sprintf('Axial View (Slice %d/%d, Max Dose)', centralSlice, size(dose, 3)));
xlabel('i (rows)');
ylabel('j (columns)');
hold on;
plot(maxI, maxJ, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;

% Coronal view (at max dose location)
subplot(1, 3, 2);
coronalSlice = squeeze(dose(:, maxJ, :));
imagesc(coronalSlice');
axis equal tight;
colorbar;
colormap(jet);
title(sprintf('Coronal View (Column %d/%d)', maxJ, size(dose, 2)));
xlabel('i (rows)');
ylabel('k (slices)');
hold on;
plot(maxI, maxK, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;

% Sagittal view (at max dose location)
subplot(1, 3, 3);
sagittalSlice = squeeze(dose(maxI, :, :));
imagesc(sagittalSlice');
axis equal tight;
colorbar;
colormap(jet);
title(sprintf('Sagittal View (Row %d/%d)', maxI, size(dose, 1)));
xlabel('j (columns)');
ylabel('k (slices)');
hold on;
plot(maxJ, maxK, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
hold off;

% Add overall title
sgtitle(sprintf('RTDOSE: %s - Max Dose: %.2f Gy', rtdoseFile(1).name, maxDose), ...
    'Interpreter', 'none');

fprintf('  - Visualization complete\n');

%% Additional Analysis - Dose Profile Through Maximum
figure('Name', 'Dose Profiles Through Maximum', 'Position', [100, 650, 1400, 400]);

% Profile along i-axis (rows)
subplot(1, 3, 1);
profile_i = dose(:, maxJ, maxK);
plot(1:length(profile_i), profile_i, 'b-', 'LineWidth', 2);
hold on;
plot(maxI, maxDose, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Position along i (rows)');
ylabel('Dose (Gy)');
title('Dose Profile Along Rows (i-axis)');

% Profile along j-axis (columns)
subplot(1, 3, 2);
profile_j = dose(maxI, :, maxK);
plot(1:length(profile_j), profile_j, 'b-', 'LineWidth', 2);
hold on;
plot(maxJ, maxDose, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Position along j (columns)');
ylabel('Dose (Gy)');
title('Dose Profile Along Columns (j-axis)');

% Profile along k-axis (slices)
subplot(1, 3, 3);
profile_k = squeeze(dose(maxI, maxJ, :));
plot(1:length(profile_k), profile_k, 'b-', 'LineWidth', 2);
hold on;
plot(maxK, maxDose, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Position along k (slices)');
ylabel('Dose (Gy)');
title('Dose Profile Along Slices (k-axis)');

sgtitle('1D Dose Profiles Through Maximum Dose Point');

%% Summary
fprintf('\n========================================\n');
fprintf('RTDOSE Inspection Complete\n');
fprintf('========================================\n');
fprintf('Summary:\n');
fprintf('  - File: %s\n', rtdoseFile(1).name);
fprintf('  - Dimensions: %d x %d x %d\n', size(dose, 1), size(dose, 2), size(dose, 3));
fprintf('  - Resolution: [%.2f, %.2f, %.2f] mm\n', ...
    rtdoseInfo.PixelSpacing(1), rtdoseInfo.PixelSpacing(2), ...
    (isfield(rtdoseInfo, 'GridFrameOffsetVector') && length(rtdoseInfo.GridFrameOffsetVector) > 1) * ...
    abs(rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1)));
fprintf('  - Max dose: %.2f Gy at position [%d, %d, %d]\n', maxDose, maxI, maxJ, maxK);
fprintf('  - Dose Summation Type: %s\n', rtdoseInfo.DoseSummationType);
fprintf('========================================\n');