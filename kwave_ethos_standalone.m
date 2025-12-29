%% KWAVE_ETHOS_STANDALONE
% Self-contained k-Wave Simulation for ETHOS Field-by-Field Dose Reconstruction
%
% This script provides a complete workflow for:
% 1. Loading SCT and converting HU to material properties
% 2. Loading field-by-field dose grids
% 3. Converting dose to initial pressure
% 4. Running k-Wave forward simulation
% 5. Time reversal reconstruction
% 6. Converting back to dose and summing all fields
%
% Pressure-Dose Conversion: p0(r) = D(r) * Gamma(r) * rho(r)
%
% Author: Generated for ETHOS dose analysis
% Date: 2025

close all;
clear;
clc;

%% ========================================================================
%                           CONFIGURATION
% =========================================================================

fprintf('\n============================================================\n');
fprintf('    k-Wave ETHOS Dose Reconstruction - Standalone Version\n');
fprintf('============================================================\n\n');

% --- Path Configuration ---
% Adjust these paths to match your data location
dataPath = '/mnt/weka/home/80030361/ETHOS_Simulations';
patientID = 'PatientXX';  % Replace with your patient ID
sessionID = 'Session01';  % Replace with your session ID

% --- Dose per pulse ---
dosePerPulse_cGy = 0.16;  % cGy per pulse

% --- Simulation Parameters ---
scale = 1;                % Grid scaling (1 = full resolution, <1 for faster testing)
pmlSize = 10;             % Perfectly Matched Layer size (voxels)
useGPU = true;            % Use GPU if available
plotSimulation = false;   % Show simulation progress plots

% --- Time Reversal Parameters ---
numTRIterations = 10;     % Number of iterative TR refinements
positivityConstraint = true;  % Enforce positive pressure/dose

% --- Output Settings ---
saveResults = true;
outputPath = './kwave_ethos_results';

%% ========================================================================
%                    LOAD DATA (SCT AND DOSE GRIDS)
% =========================================================================
% 
% This section loads your data. Modify the loading code to match your
% actual data format and file locations.

fprintf('--- Loading Data ---\n');

% -------------------------------------------------------------------------
% OPTION 1: Load from DICOM files
% Uncomment and modify as needed
% -------------------------------------------------------------------------
% 
% % Load SCT (Simulation CT)
% sctPath = fullfile(dataPath, patientID, 'CT');
% [sctVolume, sctInfo] = loadDICOMVolume(sctPath);
% 
% % Load RTPLAN for MU/pulse information
% rtplanPath = fullfile(dataPath, patientID, 'RTPLAN');
% rtplanInfo = dicominfo(rtplanPath);
% 
% % Load field dose grids (from Step 1 output)
% dosePath = fullfile(dataPath, patientID, 'FieldDoses');
% doseFiles = dir(fullfile(dosePath, 'field_*.mat'));

% -------------------------------------------------------------------------
% OPTION 2: Load from pre-processed .mat files
% -------------------------------------------------------------------------

% For testing, you can create synthetic data or load your existing files
% Example: Load previously saved data

% Check if test data exists, otherwise create synthetic phantom
if exist('test_sct_data.mat', 'file')
    fprintf('  Loading existing test data...\n');
    load('test_sct_data.mat');
else
    fprintf('  Creating synthetic test phantom...\n');
    [sctVolume, fieldDoses, fieldNames, fieldMUs, sctInfo] = createTestPhantom();
end

%% ========================================================================
%                    CONVERT HU TO MATERIAL PROPERTIES
% =========================================================================

fprintf('\n--- Converting HU to Material Properties ---\n');

[density, soundSpeed, attenuation, gruneisen, materialMap] = ...
    hu2MaterialProps(sctVolume);

fprintf('  Density range: [%.1f, %.1f] kg/m³\n', min(density(:)), max(density(:)));
fprintf('  Sound speed range: [%.1f, %.1f] m/s\n', min(soundSpeed(:)), max(soundSpeed(:)));

%% ========================================================================
%                    CALCULATE NUMBER OF PULSES
% =========================================================================

fprintf('\n--- Calculating Pulses per Field ---\n');

numFields = length(fieldDoses);
dosePerPulse_Gy = dosePerPulse_cGy / 100;

% Calculate pulses based on MU weighting
totalMU = sum(fieldMUs);
fieldPulses = zeros(numFields, 1);

for i = 1:numFields
    % Estimate pulses based on maximum dose in field
    maxFieldDose = max(fieldDoses{i}(:));
    fieldPulses(i) = round(maxFieldDose / dosePerPulse_Gy);
    
    % Ensure at least 1 pulse
    fieldPulses(i) = max(fieldPulses(i), 1);
    
    fprintf('  Field %d (%s): %.1f MU -> %d pulses (max dose: %.4f Gy)\n', ...
        i, fieldNames{i}, fieldMUs(i), fieldPulses(i), maxFieldDose);
end

%% ========================================================================
%                    SETUP k-WAVE GRID
% =========================================================================

fprintf('\n--- Setting up k-Wave Grid ---\n');

[Nx, Ny, Nz] = size(sctVolume);

% Apply scaling if needed
if scale ~= 1
    Nx = round(Nx * scale);
    Ny = round(Ny * scale);
    Nz = round(Nz * scale);
    
    density = imresize3(density, [Nx, Ny, Nz], 'linear');
    soundSpeed = imresize3(soundSpeed, [Nx, Ny, Nz], 'linear');
    attenuation = imresize3(attenuation, [Nx, Ny, Nz], 'linear');
    gruneisen = imresize3(gruneisen, [Nx, Ny, Nz], 'linear');
    
    for i = 1:numFields
        fieldDoses{i} = imresize3(fieldDoses{i}, [Nx, Ny, Nz], 'linear');
    end
end

% Grid spacing (mm to m)
dx = sctInfo.PixelSpacing(1) / 1000 / scale;
dy = sctInfo.PixelSpacing(2) / 1000 / scale;
dz = sctInfo.SliceThickness / 1000 / scale;

% Create k-Wave grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% Time step (CFL-based)
cfl = 0.3;
c_max = max(soundSpeed(:));
dt_stability = cfl * min([dx, dy, dz]) / c_max;

% Simulation time (acoustic wave round-trip)
maxDist = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
c_min = min(soundSpeed(soundSpeed > 100));  % Exclude air
t_end = 2.5 * maxDist / c_min;

kgrid.dt = dt_stability;
kgrid.Nt = ceil(t_end / kgrid.dt);

fprintf('  Grid: [%d x %d x %d]\n', Nx, Ny, Nz);
fprintf('  Spacing: [%.3f, %.3f, %.3f] mm\n', dx*1000, dy*1000, dz*1000);
fprintf('  dt: %.3e s, Nt: %d\n', kgrid.dt, kgrid.Nt);
fprintf('  Simulation time: %.1f µs\n', t_end * 1e6);

%% ========================================================================
%                    SETUP MEDIUM
% =========================================================================

fprintf('\n--- Setting up Medium Properties ---\n');

medium.sound_speed = soundSpeed;
medium.density = density;
medium.alpha_coeff = attenuation;
medium.alpha_power = 1.5;

%% ========================================================================
%                    SETUP SENSOR (Planar, Inferior)
% =========================================================================

fprintf('\n--- Setting up Sensor ---\n');

% Sum all dose grids to find extent
allDose = zeros(Nx, Ny, Nz);
for i = 1:numFields
    allDose = allDose + fieldDoses{i};
end

% Find dose extent along Y (depth) axis
doseProfile_Y = squeeze(max(max(allDose, [], 1), [], 3));
doseThreshold = max(allDose(:)) * 0.01;

% Find inferior (largest Y) boundary of dose
inferiorIdx = find(doseProfile_Y > doseThreshold, 1, 'last');
if isempty(inferiorIdx)
    inferiorIdx = round(Ny * 0.9);
end

% Place sensor 5mm inferior to dose
sensorOffset = ceil(5 / (dy * 1000));
sensorY = min(inferiorIdx + sensorOffset, Ny - pmlSize - 1);

% Create planar sensor
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, sensorY, :) = 1;

fprintf('  Sensor Y position: %d (%.1f mm from origin)\n', sensorY, sensorY * dy * 1000);
fprintf('  Sensor elements: %d\n', sum(sensor.mask(:)));

%% ========================================================================
%                    RUN SIMULATIONS
% =========================================================================

fprintf('\n============================================================\n');
fprintf('           Running Field-by-Field Simulations\n');
fprintf('============================================================\n');

% GPU/CPU selection
if useGPU && gpuDeviceCount > 0
    dataCast = 'gpuArray-single';
    fprintf('Using GPU acceleration\n');
else
    dataCast = 'single';
    fprintf('Using CPU\n');
end

input_args = {'Smooth', false, ...
              'PMLInside', false, ...
              'PMLSize', pmlSize, ...
              'DataCast', dataCast, ...
              'PlotSim', plotSimulation};

% Storage for results
reconstructedDoses = cell(numFields, 1);

for fieldIdx = 1:numFields
    fprintf('\n--- Field %d/%d: %s ---\n', fieldIdx, numFields, fieldNames{fieldIdx});
    
    % Get dose per pulse
    doseGrid = fieldDoses{fieldIdx};
    dosePerPulse = doseGrid / fieldPulses(fieldIdx);
    
    % Convert to initial pressure: p0 = D * Gamma * rho
    p0 = dosePerPulse .* gruneisen .* density;
    
    % Ensure positive
    p0(p0 < 0) = 0;
    
    fprintf('  Dose range: [%.6f, %.6f] Gy/pulse\n', ...
        min(dosePerPulse(:)), max(dosePerPulse(:)));
    fprintf('  p0 range: [%.2f, %.2f] Pa\n', min(p0(:)), max(p0(:)));
    
    % Skip if no significant pressure
    if max(p0(:)) < 1e-6
        fprintf('  Skipping - no significant initial pressure\n');
        reconstructedDoses{fieldIdx} = zeros(Nx, Ny, Nz);
        continue;
    end
    
    % === Forward Simulation ===
    fprintf('  Forward simulation...\n');
    source.p0 = p0;
    
    tic;
    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    fprintf('    Completed in %.1f s\n', toc);
    
    % === Time Reversal ===
    fprintf('  Time reversal reconstruction...\n');
    
    source = rmfield(source, 'p0');
    source.p_mask = sensor.mask;
    source.p = fliplr(sensor_data);
    source.p_mode = 'dirichlet';
    sensor.record = {'p_final'};
    
    tic;
    p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    
    if positivityConstraint
        p0_recon.p_final = p0_recon.p_final .* (p0_recon.p_final > 0);
    end
    fprintf('    Initial TR completed in %.1f s\n', toc);
    
    % === Iterative Refinement ===
    if numTRIterations > 1
        fprintf('  Iterative refinement (%d iterations)...\n', numTRIterations - 1);
        
        for iter = 2:numTRIterations
            % Forward with current estimate
            src_iter.p0 = gather(p0_recon.p_final);
            sens_iter.mask = sensor.mask;
            
            data_iter = kspaceFirstOrder3D(kgrid, medium, src_iter, sens_iter, input_args{:});
            
            % Residual
            residual = sensor_data - data_iter;
            
            % TR of residual
            src_iter = rmfield(src_iter, 'p0');
            src_iter.p_mask = sensor.mask;
            src_iter.p = fliplr(residual);
            src_iter.p_mode = 'dirichlet';
            sens_iter.record = {'p_final'};
            
            p0_update = kspaceFirstOrder3D(kgrid, medium, src_iter, sens_iter, input_args{:});
            
            % Update
            p0_recon.p_final = p0_recon.p_final + p0_update.p_final;
            
            if positivityConstraint
                p0_recon.p_final = p0_recon.p_final .* (p0_recon.p_final > 0);
            end
            
            if mod(iter, 5) == 0
                fprintf('    Iteration %d/%d\n', iter, numTRIterations);
            end
        end
    end
    
    % === Convert Back to Dose ===
    p0_final = gather(p0_recon.p_final);
    
    % D = p0 / (Gamma * rho)
    denom = gruneisen .* density;
    denom(denom < 1e-6) = 1e-6;
    
    reconDose_perPulse = p0_final ./ denom;
    
    % Scale by pulses
    reconDose_total = reconDose_perPulse * fieldPulses(fieldIdx);
    
    reconstructedDoses{fieldIdx} = reconDose_total;
    
    fprintf('  Reconstructed dose: [%.6f, %.6f] Gy\n', ...
        min(reconDose_total(:)), max(reconDose_total(:)));
    
    % Clear GPU
    if useGPU && gpuDeviceCount > 0
        reset(gpuDevice);
    end
end

%% ========================================================================
%                    SUM ALL FIELDS
% =========================================================================

fprintf('\n============================================================\n');
fprintf('              Summing Field Doses\n');
fprintf('============================================================\n');

totalReconDose = zeros(Nx, Ny, Nz);
totalOrigDose = zeros(Nx, Ny, Nz);

for i = 1:numFields
    totalReconDose = totalReconDose + reconstructedDoses{i};
    totalOrigDose = totalOrigDose + fieldDoses{i};
end

fprintf('Original total dose: [%.6f, %.6f] Gy\n', ...
    min(totalOrigDose(:)), max(totalOrigDose(:)));
fprintf('Reconstructed total dose: [%.6f, %.6f] Gy\n', ...
    min(totalReconDose(:)), max(totalReconDose(:)));

%% ========================================================================
%                    VISUALIZATION
% =========================================================================

fprintf('\n--- Generating Visualizations ---\n');

% Find max dose location
[~, maxIdx] = max(totalReconDose(:));
[maxX, maxY, maxZ] = ind2sub(size(totalReconDose), maxIdx);

% Coordinates in cm
x_cm = (0:Nx-1) * dx * 100;
y_cm = (0:Ny-1) * dy * 100;
z_cm = (0:Nz-1) * dz * 100;

doseMax = max(max(totalOrigDose(:)), max(totalReconDose(:)));

% --- Figure 1: Side-by-side comparison ---
figure('Name', 'Dose Comparison', 'Position', [100, 100, 1400, 800]);

% Original
subplot(2,3,1);
imagesc(z_cm, x_cm, squeeze(totalOrigDose(:, maxY, :)));
axis image; colorbar; xlabel('Z (cm)'); ylabel('X (cm)');
title(sprintf('Original - Axial (Y=%.1f cm)', y_cm(maxY)));
caxis([0 doseMax]);

subplot(2,3,2);
imagesc(z_cm, y_cm, squeeze(totalOrigDose(maxX, :, :)));
axis image; colorbar; xlabel('Z (cm)'); ylabel('Y (cm)');
title('Original - Sagittal');
caxis([0 doseMax]);

subplot(2,3,3);
imagesc(y_cm, x_cm, squeeze(totalOrigDose(:, :, maxZ)));
axis image; colorbar; xlabel('Y (cm)'); ylabel('X (cm)');
title('Original - Coronal');
caxis([0 doseMax]);

% Reconstructed
subplot(2,3,4);
imagesc(z_cm, x_cm, squeeze(totalReconDose(:, maxY, :)));
axis image; colorbar; xlabel('Z (cm)'); ylabel('X (cm)');
title('Reconstructed - Axial');
caxis([0 doseMax]);

subplot(2,3,5);
imagesc(z_cm, y_cm, squeeze(totalReconDose(maxX, :, :)));
axis image; colorbar; xlabel('Z (cm)'); ylabel('Y (cm)');
title('Reconstructed - Sagittal');
caxis([0 doseMax]);

subplot(2,3,6);
imagesc(y_cm, x_cm, squeeze(totalReconDose(:, :, maxZ)));
axis image; colorbar; xlabel('Y (cm)'); ylabel('X (cm)');
title('Reconstructed - Coronal');
caxis([0 doseMax]);

colormap('jet');
sgtitle('Original vs Reconstructed Dose', 'FontSize', 14, 'FontWeight', 'bold');

% --- Figure 2: Difference ---
figure('Name', 'Dose Difference', 'Position', [150, 150, 1200, 400]);

doseDiff = totalReconDose - totalOrigDose;
diffMax = max(abs(doseDiff(:)));

subplot(1,3,1);
imagesc(z_cm, x_cm, squeeze(doseDiff(:, maxY, :)));
axis image; colorbar; xlabel('Z (cm)'); ylabel('X (cm)');
title('Difference - Axial');
caxis([-diffMax diffMax]);

subplot(1,3,2);
imagesc(z_cm, y_cm, squeeze(doseDiff(maxX, :, :)));
axis image; colorbar; xlabel('Z (cm)'); ylabel('Y (cm)');
title('Difference - Sagittal');
caxis([-diffMax diffMax]);

subplot(1,3,3);
imagesc(y_cm, x_cm, squeeze(doseDiff(:, :, maxZ)));
axis image; colorbar; xlabel('Y (cm)'); ylabel('X (cm)');
title('Difference - Coronal');
caxis([-diffMax diffMax]);

colormap('RdBu');
sgtitle('Dose Difference (Recon - Orig)', 'FontSize', 14);

% --- Figure 3: Profiles ---
figure('Name', 'Dose Profiles', 'Position', [200, 200, 1200, 400]);

% PDD
subplot(1,3,1);
pdd_orig = squeeze(totalOrigDose(maxX, :, maxZ));
pdd_recon = squeeze(totalReconDose(maxX, :, maxZ));
plot(y_cm, pdd_orig / max(pdd_orig) * 100, 'b-', 'LineWidth', 2); hold on;
plot(y_cm, pdd_recon / max(pdd_recon) * 100, 'r--', 'LineWidth', 2);
xlabel('Depth (cm)'); ylabel('Relative Dose (%)');
title('PDD'); legend('Original', 'Reconstructed'); grid on;

% Lateral X
subplot(1,3,2);
lat_x_orig = squeeze(totalOrigDose(:, maxY, maxZ));
lat_x_recon = squeeze(totalReconDose(:, maxY, maxZ));
plot(x_cm, lat_x_orig / max(lat_x_orig) * 100, 'b-', 'LineWidth', 2); hold on;
plot(x_cm, lat_x_recon / max(lat_x_recon) * 100, 'r--', 'LineWidth', 2);
xlabel('X (cm)'); ylabel('Relative Dose (%)');
title('Lateral Profile (X)'); legend('Original', 'Reconstructed'); grid on;

% Lateral Z
subplot(1,3,3);
lat_z_orig = squeeze(totalOrigDose(maxX, maxY, :));
lat_z_recon = squeeze(totalReconDose(maxX, maxY, :));
plot(z_cm, lat_z_orig / max(lat_z_orig) * 100, 'b-', 'LineWidth', 2); hold on;
plot(z_cm, lat_z_recon / max(lat_z_recon) * 100, 'r--', 'LineWidth', 2);
xlabel('Z (cm)'); ylabel('Relative Dose (%)');
title('Lateral Profile (Z)'); legend('Original', 'Reconstructed'); grid on;

sgtitle('Dose Profile Comparison', 'FontSize', 14);

%% ========================================================================
%                    SAVE RESULTS
% =========================================================================

if saveResults
    fprintf('\n--- Saving Results ---\n');
    
    if ~exist(outputPath, 'dir')
        mkdir(outputPath);
    end
    
    save(fullfile(outputPath, 'kwave_results.mat'), ...
        'totalOrigDose', 'totalReconDose', 'reconstructedDoses', ...
        'fieldNames', 'fieldPulses', 'kgrid', 'dx', 'dy', 'dz');
    
    saveas(figure(1), fullfile(outputPath, 'dose_comparison.png'));
    saveas(figure(2), fullfile(outputPath, 'dose_difference.png'));
    saveas(figure(3), fullfile(outputPath, 'dose_profiles.png'));
    
    fprintf('  Saved to: %s\n', outputPath);
end

fprintf('\n============================================================\n');
fprintf('                 Simulation Complete!\n');
fprintf('============================================================\n\n');

%% ========================================================================
%                    HELPER FUNCTIONS
% =========================================================================

function [density, soundSpeed, attenuation, gruneisen, materialMap] = ...
    hu2MaterialProps(huVolume)
%HU2MATERIALPROPS Convert HU to acoustic properties
%
%   HU Thresholds:
%       Air:         ≤ -950
%       Lung:        -950 to -500
%       Fat:         -500 to -50
%       Water:       -50 to 10
%       Blood:       10 to 55
%       Muscle:      55 to 100
%       Soft Tissue: 100 to 300
%       Bone:        300 to 3000
%       Metal:       > 3000

    [Nx, Ny, Nz] = size(huVolume);
    
    density = zeros(Nx, Ny, Nz);
    soundSpeed = zeros(Nx, Ny, Nz);
    attenuation = zeros(Nx, Ny, Nz);
    gruneisen = ones(Nx, Ny, Nz);  % Placeholder = 1
    materialMap = zeros(Nx, Ny, Nz, 'uint8');
    
    % Material properties: [density, sound_speed, attenuation, gruneisen]
    props.air = [1.2, 343, 0.0, 1.0];
    props.lung = [300, 600, 0.5, 1.0];
    props.fat = [920, 1450, 0.48, 1.0];
    props.water = [1000, 1500, 0.002, 1.0];
    props.blood = [1060, 1575, 0.14, 1.0];
    props.muscle = [1050, 1580, 0.57, 1.0];
    props.soft = [1040, 1540, 0.5, 1.0];
    props.bone = [1900, 3500, 4.0, 1.0];
    props.metal = [4500, 5000, 0.0, 1.0];
    
    % Segment and assign
    masks = struct();
    masks.air = huVolume <= -950;
    masks.lung = huVolume > -950 & huVolume <= -500;
    masks.fat = huVolume > -500 & huVolume <= -50;
    masks.water = huVolume > -50 & huVolume <= 10;
    masks.blood = huVolume > 10 & huVolume <= 55;
    masks.muscle = huVolume > 55 & huVolume <= 100;
    masks.soft = huVolume > 100 & huVolume <= 300;
    masks.bone = huVolume > 300 & huVolume <= 3000;
    masks.metal = huVolume > 3000;
    
    materials = fieldnames(masks);
    for i = 1:length(materials)
        mat = materials{i};
        m = masks.(mat);
        p = props.(mat);
        
        density(m) = p(1);
        soundSpeed(m) = p(2);
        attenuation(m) = p(3);
        gruneisen(m) = p(4);
        materialMap(m) = i;
    end
    
    % Handle bone with interpolation
    bone_mask = masks.bone;
    if any(bone_mask(:))
        bone_hu = huVolume(bone_mask);
        frac = (bone_hu - 300) / (3000 - 300);
        frac = max(0, min(1, frac));
        
        density(bone_mask) = 1100 + frac .* (1900 - 1100);
        soundSpeed(bone_mask) = 2500 + frac .* (4000 - 2500);
        attenuation(bone_mask) = 2.0 + frac .* (6.0 - 2.0);
    end
    
    % Ensure minimum values
    density = max(density, 100);
    soundSpeed = max(soundSpeed, 300);
    
    % Print summary
    total = numel(huVolume);
    fprintf('    Material Distribution:\n');
    for i = 1:length(materials)
        mat = materials{i};
        pct = 100 * sum(masks.(mat)(:)) / total;
        if pct > 0.01
            fprintf('      %-12s: %6.2f%%\n', mat, pct);
        end
    end
end

function [sctVolume, fieldDoses, fieldNames, fieldMUs, sctInfo] = createTestPhantom()
%CREATETESTPHANTOM Create synthetic phantom for testing

    % Grid size
    Nx = 64; Ny = 64; Nz = 64;
    
    % Create water phantom with embedded structures
    sctVolume = zeros(Nx, Ny, Nz);  % Water = 0 HU
    
    % Add tissue variations
    [X, Y, Z] = ndgrid(1:Nx, 1:Ny, 1:Nz);
    center = [Nx/2, Ny/2, Nz/2];
    
    % Outer body (muscle)
    bodyMask = sqrt((X-center(1)).^2 + (Y-center(2)).^2 + (Z-center(3)).^2) < 25;
    sctVolume(bodyMask) = 50;  % Muscle HU
    
    % Fat layer
    fatMask = sqrt((X-center(1)).^2 + (Y-center(2)).^2 + (Z-center(3)).^2) < 20;
    sctVolume(fatMask & ~bodyMask) = -100;
    
    % Target (soft tissue)
    targetMask = sqrt((X-center(1)).^2 + (Y-center(2)).^2 + (Z-center(3)).^2) < 8;
    sctVolume(targetMask) = 200;
    
    % Bone structure
    boneMask = (abs(X - center(1)) < 3) & (Y > center(2) + 5) & (Y < center(2) + 15);
    sctVolume(boneMask) = 1000;
    
    % Air outside
    sctVolume(~bodyMask) = -1000;
    
    % Info structure
    sctInfo.PixelSpacing = [2.0, 2.0];  % mm
    sctInfo.SliceThickness = 2.5;       % mm
    
    % Create test dose grids (3 fields, simple geometry)
    numFields = 3;
    fieldDoses = cell(numFields, 1);
    fieldNames = {'Field_AP', 'Field_LAO', 'Field_RAO'};
    fieldMUs = [100, 80, 80];
    
    for i = 1:numFields
        dose = zeros(Nx, Ny, Nz);
        
        % Simple dose distribution centered on target
        switch i
            case 1  % AP
                for y = 1:Ny
                    depth = abs(y - center(2));
                    pdd = exp(-0.05 * depth) * (1 - exp(-0.3 * depth));
                    dose(:, y, :) = targetMask(:, y, :) .* pdd;
                end
            case 2  % LAO
                for x = 1:Nx
                    depth = abs(x - center(1) - 10);
                    pdd = exp(-0.05 * depth) * (1 - exp(-0.3 * depth));
                    dose(x, :, :) = targetMask(x, :, :) .* pdd;
                end
            case 3  % RAO
                for x = 1:Nx
                    depth = abs(x - center(1) + 10);
                    pdd = exp(-0.05 * depth) * (1 - exp(-0.3 * depth));
                    dose(x, :, :) = targetMask(x, :, :) .* pdd;
                end
        end
        
        % Normalize to reasonable dose
        dose = dose / max(dose(:)) * 0.5 * (fieldMUs(i) / 100);  % Scale by MU
        fieldDoses{i} = dose;
    end
    
    fprintf('  Created %dx%dx%d test phantom with %d fields\n', Nx, Ny, Nz, numFields);
end
