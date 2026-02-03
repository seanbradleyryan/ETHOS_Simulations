%% ETHOS Field-by-Field Photoacoustic Dose Reconstruction
% Purpose: Simulate photoacoustic wave propagation using kWave for each
%          radiation field and reconstruct dose via time reversal
% 
% Workflow:
%   1. Load field doses from ethos_field_dosev3.m output (RTDOSE grid)
%   2. Load CT resampled to RTDOSE grid dimensions
%   3. Process CT to extract density and assign acoustic properties
%   4. Convert dose to initial pressure: p0 = D * Gamma * rho
%   5. Run kWave forward simulation
%   6. Perform time reversal reconstruction
%   7. Convert pressure back to dose
%   8. Sum all field contributions
%   9. Visualize results
%
% MODIFIED: Now loads CT resampled to RTDOSE grid from ctResampled.mat
%           instead of using original high-resolution CT
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
fieldDoseDir = fullfile(baseDir, 'RayStationFiles', patientID, sessionName);
dicomDir = fullfile(baseDir, 'EthosExports', patientID, 'Pancreas', sessionName, 'sct');
outputDir = fullfile(baseDir, 'kWaveReconstruction', patientID, sessionName);

% Dose per pulse (in cGy)
dosePerPulse_cGy = 0.16;
dosePerPulse_Gy = dosePerPulse_cGy / 100;  % Convert to Gy

% kWave simulation parameters
pmlSize = 10;                   % Perfectly matched layer size
useGPU = true;                  % Use GPU acceleration if available
numTimeReversalIterations = 1;  % Number of TR iterations (1 = single pass)

% Gruneisen coefficient (placeholder - to be filled with tissue-specific values later)
gruneisen_default = 1.0;

%% ======================== MATERIAL PROPERTIES ========================
% Define acoustic properties for each tissue type based on HU thresholds
% Format: [HU_min, HU_max, density (kg/m), speed_of_sound (m/s), alpha_coeff, alpha_power]

materialTable = struct();

% % Air: HU ~ -1000
% materialTable.air.HU_range = [-1024, -900];
% materialTable.air.density = 1000;
% materialTable.air.soundSpeed = 1480;
% materialTable.air.alphaCoeff = 0;
% materialTable.air.alphaPower = 1;
% materialTable.air.gruneisen = 0;  % No photoacoustic signal in air

% % Lung: HU -700 to -600
% materialTable.lung.HU_range = [-900, -500];
% materialTable.lung.density = 400;
% materialTable.lung.soundSpeed = 600;
% materialTable.lung.alphaCoeff = 0.5;
% materialTable.lung.alphaPower = 1.5;
% materialTable.lung.gruneisen = 0.5;

% Fat: HU -120 to -90
materialTable.fat.HU_range = [-199, -50];
materialTable.fat.density = 920;
materialTable.fat.soundSpeed = 1450;
materialTable.fat.alphaCoeff = 0.48;
materialTable.fat.alphaPower = 1.5;
materialTable.fat.gruneisen = 0.7;

% Water: HU ~ 0 override air
materialTable.water.HU_range = [-1000, -200];
materialTable.water.density = 1000;
materialTable.water.soundSpeed = 1480;
materialTable.water.alphaCoeff = 0.0022;
materialTable.water.alphaPower = 2;
materialTable.water.gruneisen = 0.11;

% Blood: HU +13 to +50
% materialTable.blood.HU_range = [13, 50];
% materialTable.blood.density = 1060;
% materialTable.blood.soundSpeed = 1575;
% materialTable.blood.alphaCoeff = 0.2;
% materialTable.blood.alphaPower = 1.3;
% materialTable.blood.gruneisen = 0.15;

% % Muscle: HU +35 to +55 (overlaps with blood, use higher priority)
% materialTable.muscle.HU_range = [51, 80];
% materialTable.muscle.density = 1050;
% materialTable.muscle.soundSpeed = 1580;
% materialTable.muscle.alphaCoeff = 0.5;
% materialTable.muscle.alphaPower = 1;
% materialTable.muscle.gruneisen = 0.2;

% Soft Tissue: HU +100 to +300
materialTable.softTissue.HU_range = [-50, 100];
materialTable.softTissue.density = 1040;
materialTable.softTissue.soundSpeed = 1540;
materialTable.softTissue.alphaCoeff = 0.5;
materialTable.softTissue.alphaPower = 1.1;
materialTable.softTissue.gruneisen = 1;

% Bone: HU +700 to +3000
materialTable.bone.HU_range = [100, 3000];
materialTable.bone.density = 1900;
materialTable.bone.soundSpeed = 3200;
materialTable.bone.alphaCoeff = 10;
materialTable.bone.alphaPower = 1;
materialTable.bone.gruneisen = 10;  % No signal propagation through bone

% Metal: HU > 3000
materialTable.metal.HU_range = [3001, 10000];
materialTable.metal.density = 7800;
materialTable.metal.soundSpeed = 5900;
materialTable.metal.alphaCoeff = 0;
materialTable.metal.alphaPower = 1;
materialTable.metal.gruneisen = 0;

%% ======================== VERIFY PATHS ========================
fprintf('=======================================================\n');
fprintf('  ETHOS Field-by-Field Photoacoustic Reconstruction\n');
fprintf('=======================================================\n\n');

fprintf('[1/9] Verifying paths and dependencies...\n');

% Check if field dose directory exists
if ~exist(fieldDoseDir, 'dir')
    error('Field dose directory does not exist: %s\nRun ethos_field_dosev3.m first.', fieldDoseDir);
end
fprintf('  - Field dose directory: OK\n');

% Check if DICOM directory exists
if ~exist(dicomDir, 'dir')
    error('DICOM directory does not exist: %s', dicomDir);
end
fprintf('  - DICOM directory: OK\n');

% Create output directory
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('  - Created output directory: %s\n', outputDir);
else
    fprintf('  - Output directory: OK\n');
end

% Check for kWave
if ~exist('kWaveGrid', 'file')
    error('kWave toolbox not found. Please add kWave to MATLAB path.');
end
fprintf('  - kWave toolbox: OK\n');

%% ======================== LOAD RESAMPLED CT AND FIELD DOSES ========================
fprintf('\n[2/9] Loading resampled CT and field doses from Step 1...\n');

% Check for resampled CT file (preferred)
ctResampledFile = fullfile(fieldDoseDir, 'ctResampled.mat');
fieldDoseFile = fullfile(fieldDoseDir, 'fieldDoses.mat');

if exist(ctResampledFile, 'file')
    fprintf('  Loading resampled CT from: ctResampled.mat\n');
    ctData = load(ctResampledFile);
    
    % Extract resampled CT structure
    if isfield(ctData, 'ctResampled_struct')
        ct = ctData.ctResampled_struct;
        fprintf('   Loaded resampled CT structure\n');
    else
        error('ctResampled.mat does not contain ctResampled_struct');
    end
    
    % Extract dose grid info
    if isfield(ctData, 'doseGrid')
        doseGrid = ctData.doseGrid;
    end
    if isfield(ctData, 'doseSpatial')
        doseSpatial = ctData.doseSpatial;
    end
    if isfield(ctData, 'ctSpatial')
        ctSpatial = ctData.ctSpatial;
    end
    
    usingResampledCT = true;
else
    fprintf('   Resampled CT file not found: %s\n', ctResampledFile);
    fprintf('  Falling back to original CT from fieldDoses.mat\n');
    usingResampledCT = false;
end

% Load field doses
if ~exist(fieldDoseFile, 'file')
    error('Field doses file not found: %s\nRun ethos_field_dosev3.m first.', fieldDoseFile);
end

loadedData = load(fieldDoseFile);

% Check which field dose format is available (resampled vs original)
if isfield(loadedData, 'fieldDosesResampled')
    fieldDoses = loadedData.fieldDosesResampled;
    fprintf('   Loaded resampled field doses (RTDOSE grid)\n');
else
    fieldDoses = loadedData.fieldDoses;
    fprintf('   Loaded original field doses (CT grid)\n');
end

% Load CT if not already loaded from ctResampled.mat
if ~usingResampledCT
    if isfield(loadedData, 'ctResampled_struct')
        ct = loadedData.ctResampled_struct;
        usingResampledCT = true;
        fprintf('   Found resampled CT in fieldDoses.mat\n');
    elseif isfield(loadedData, 'ct')
        ct = loadedData.ct;
        fprintf('   Using original CT (not resampled)\n');
    else
        error('No CT data found in loaded files');
    end
end

% Load additional data
if isfield(loadedData, 'pln')
    pln = loadedData.pln;
end
if isfield(loadedData, 'stf')
    stf = loadedData.stf;
end
if isfield(loadedData, 'referenceDose')
    referenceDose = loadedData.referenceDose;
end
if isfield(loadedData, 'doseGrid') && ~exist('doseGrid', 'var')
    doseGrid = loadedData.doseGrid;
end

% Count valid fields
numFields = length(fieldDoses);
validFields = find(~cellfun(@isempty, fieldDoses));
numValidFields = length(validFields);

fprintf('  - Loaded %d field doses (%d valid)\n', numFields, numValidFields);

% Display CT information
if isfield(ct, 'cubeDim')
    fprintf('  - CT dimensions: %d x %d x %d\n', ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));
elseif isfield(ct, 'cubeHU')
    ctSize = size(ct.cubeHU{1});
    fprintf('  - CT dimensions: %d x %d x %d\n', ctSize(1), ctSize(2), ctSize(3));
end

if isfield(ct, 'resolution')
    fprintf('  - CT resolution: [%.2f, %.2f, %.2f] mm\n', ...
        ct.resolution.x, ct.resolution.y, ct.resolution.z);
end

if usingResampledCT
    fprintf('   Using CT resampled to RTDOSE grid\n');
else
    fprintf('   Using original CT grid (may not match field dose dimensions)\n');
end

%% ======================== LOAD RTPLAN FOR PULSE CALCULATION ========================
fprintf('\n[3/9] Extracting pulse count from RTPLAN...\n');

% Find RTPLAN file
rtplanFiles = dir(fullfile(dicomDir, 'RP*.dcm'));
if isempty(rtplanFiles)
    rtplanFiles = dir(fullfile(dicomDir, '*RTPLAN*.dcm'));
end

if isempty(rtplanFiles)
    error('No RTPLAN file found in %s', dicomDir);
end

rtplanInfo = dicominfo(fullfile(rtplanFiles(1).folder, rtplanFiles(1).name));
fprintf('  - RTPLAN file: %s\n', rtplanFiles(1).name);

% Extract beam metersets and calculate pulses per beam
beamPulses = zeros(numFields, 1);
totalPlanDose_Gy = 0;

if isfield(rtplanInfo, 'FractionGroupSequence')
    fg = rtplanInfo.FractionGroupSequence.Item_1;
    
    if isfield(fg, 'ReferencedBeamSequence')
        numRefBeams = length(fieldnames(fg.ReferencedBeamSequence));
        
        for i = 1:numRefBeams
            refBeamField = sprintf('Item_%d', i);
            refBeam = fg.ReferencedBeamSequence.(refBeamField);
            
            if isfield(refBeam, 'BeamMeterset')
                % BeamMeterset is typically in MU (Monitor Units)
                % For most linacs, 1 MU  1 cGy at dmax for standard conditions
                % We use this as an approximation for the beam dose
                beamDose_cGy = refBeam.BeamMeterset;  % Approximate dose in cGy
                beamDose_Gy = beamDose_cGy / 100;
                
                % Calculate number of pulses for this beam
                beamPulses(i) = ceil(beamDose_cGy / dosePerPulse_cGy);
                totalPlanDose_Gy = totalPlanDose_Gy + beamDose_Gy;
                
                fprintf('  - Beam %d: Meterset=%.2f MU, Pulses=%d\n', ...
                    i, refBeam.BeamMeterset, beamPulses(i));
            end
        end
    end
end

totalPulses = sum(beamPulses);
fprintf('  - Total plan dose: %.2f Gy\n', totalPlanDose_Gy);
fprintf('  - Total pulses: %d\n', totalPulses);
fprintf('  - Dose per pulse: %.4f Gy\n', dosePerPulse_Gy);

%% ======================== PROCESS CT FOR ACOUSTIC PROPERTIES ========================
fprintf('\n[4/9] Processing CT for acoustic medium properties...\n');

% Get CT cube (Hounsfield Units)
ctCube = ct.cubeHU{1};
ctSize = size(ctCube);

fprintf('  - CT grid size: %d x %d x %d\n', ctSize(1), ctSize(2), ctSize(3));
fprintf('  - CT HU range: [%.0f, %.0f]\n', min(ctCube(:)), max(ctCube(:)));

% Initialize acoustic property arrays
density = zeros(ctSize);
soundSpeed = zeros(ctSize);
alphaCoeff = zeros(ctSize);
alphaPower = ones(ctSize);  % Default power law exponent
gruneisen = zeros(ctSize);

% Assign properties based on HU thresholds
% Process in order from lowest to highest HU (later assignments override earlier)
tissueTypes = {'air', 'lung', 'fat', 'water', 'blood', 'muscle', 'softTissue', 'bone', 'metal'};

for t = 1:length(tissueTypes)
    tissue = tissueTypes{t};
    mat = materialTable.(tissue);
    
    mask = (ctCube >= mat.HU_range(1)) & (ctCube <= mat.HU_range(2));
    numVoxels = sum(mask(:));
    
    if numVoxels > 0
        density(mask) = mat.density;
        soundSpeed(mask) = mat.soundSpeed;
        alphaCoeff(mask) = mat.alphaCoeff;
        alphaPower(mask) = mat.alphaPower;
        gruneisen(mask) = mat.gruneisen;
        
        fprintf('  - %s: %d voxels (%.1f%%)\n', tissue, numVoxels, 100*numVoxels/numel(ctCube));
    end
end

% Handle any unassigned voxels (default to water-like properties)
unassigned = (density == 0);
if any(unassigned(:))
    fprintf('  - Unassigned voxels: %d (setting to water)\n', sum(unassigned(:)));
    density(unassigned) = 1000;
    soundSpeed(unassigned) = 1480;
    alphaCoeff(unassigned) = 0.002;
    alphaPower(unassigned) = 2;
    gruneisen(unassigned) = 0.12;
end

% Ensure minimum values for numerical stability
density = max(density, 1);
soundSpeed = max(soundSpeed, 100);

fprintf('  - Density range: [%.0f, %.0f] kg/m\n', min(density(:)), max(density(:)));
fprintf('  - Sound speed range: [%.0f, %.0f] m/s\n', min(soundSpeed(:)), max(soundSpeed(:)));

%% ======================== SETUP KWAVE GRID ========================
fprintf('\n[5/9] Setting up kWave simulation grid...\n');

% Grid dimensions and spacing (convert from mm to m)
if isfield(ct, 'resolution')
    dx = ct.resolution.x / 1000;  % m
    dy = ct.resolution.y / 1000;  % m
    dz = ct.resolution.z / 1000;  % m
elseif exist('doseGrid', 'var')
    % Use dose grid resolution
    dx = doseGrid.resolution(2) / 1000;  % Column direction (X)
    dy = doseGrid.resolution(1) / 1000;  % Row direction (Y)
    dz = doseGrid.resolution(3) / 1000;  % Slice direction (Z)
else
    error('Cannot determine grid resolution');
end

Nx = ctSize(1);
Ny = ctSize(2);
Nz = ctSize(3);

fprintf('  - Grid size: %d x %d x %d\n', Nx, Ny, Nz);
fprintf('  - Grid spacing: [%.4f, %.4f, %.4f] mm\n', dx*1000, dy*1000, dz*1000);

% Create kWave grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% Calculate time step based on CFL condition
maxSoundSpeed = max(soundSpeed(:));
cfl = 0.3;  % CFL number for stability
dt = cfl * min([dx, dy, dz]) / maxSoundSpeed;

% Calculate simulation time (allow wave to traverse grid twice)
gridDiagonal = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
simTime = 2.5 * gridDiagonal / min(soundSpeed(:));
Nt = ceil(simTime / dt);

kgrid.dt = dt;
kgrid.Nt = Nt;

fprintf('  - Time step: %.2e s\n', dt);
fprintf('  - Number of time steps: %d\n', Nt);
fprintf('  - Simulation time: %.2e s\n', simTime);

%% ======================== SETUP MEDIUM PROPERTIES ========================
fprintf('\n[6/9] Configuring acoustic medium...\n');

medium = struct();
medium.density = density;
medium.sound_speed = soundSpeed;
medium.alpha_coeff = alphaCoeff;
medium.alpha_power = mean(alphaPower(:));  % kWave uses scalar alpha_power

fprintf('  - Medium properties assigned\n');
fprintf('  - Alpha power: %.2f\n', medium.alpha_power);

%% ======================== PROCESS EACH FIELD ========================
fprintf('\n[7/9] Processing radiation fields with kWave...\n');

% Initialize storage for reconstructed doses
reconstructedFieldDoses = cell(numFields, 1);
totalReconstructedDose = zeros(ctSize);

% kWave input arguments
if useGPU
    try
        gpuDevice;
        dataCast = 'gpuArray-single';
        fprintf('  - Using GPU acceleration\n');
    catch
        dataCast = 'single';
        fprintf('  - GPU not available, using CPU\n');
    end
else
    dataCast = 'single';
    fprintf('  - Using CPU\n');
end

inputArgs = {'Smooth', false, 'PMLInside', false, 'PMLSize', pmlSize, ...
             'DataCast', dataCast, 'PlotSim', true};

% Process each valid field
for fieldIdx = validFields'
    fprintf('\n  Processing Field %d/%d (Gantry: %.1f)...\n', ...
        fieldIdx, numFields, fieldDoses{fieldIdx}.gantryAngle);
    
    % Get field dose
    fieldDose = fieldDoses{fieldIdx}.physicalDose;
    
    % Check size compatibility with CT grid
    if ~isequal(size(fieldDose), ctSize)
        fprintf('     Field dose size mismatch with CT grid\n');
        fprintf('      Field dose: %d x %d x %d\n', size(fieldDose, 1), size(fieldDose, 2), size(fieldDose, 3));
        fprintf('      CT grid: %d x %d x %d\n', ctSize(1), ctSize(2), ctSize(3));
        
        % Attempt resampling if sizes don't match
        fprintf('    Attempting to resample field dose to CT grid...\n');
        [nRows, nCols, nSlices] = size(fieldDose);
        [iDose, jDose, kDose] = ndgrid(1:nRows, 1:nCols, 1:nSlices);
        
        scale_i = nRows / Nx;
        scale_j = nCols / Ny;
        scale_k = nSlices / Nz;
        
        iQuery = 1 + (0:Nx-1) * scale_i;
        jQuery = 1 + (0:Ny-1) * scale_j;
        kQuery = 1 + (0:Nz-1) * scale_k;
        
        [iQ, jQ, kQ] = ndgrid(iQuery, jQuery, kQuery);
        
        try
            fieldDose = interpn(iDose, jDose, kDose, fieldDose, iQ, jQ, kQ, 'linear', 0);
            fprintf('     Resampled to: %d x %d x %d\n', size(fieldDose));
        catch ME
            fprintf('     Resampling failed: %s\n', ME.message);
            continue;
        end
    end
    
    % Get number of pulses for this field
    numPulsesField = beamPulses(fieldIdx);
    if numPulsesField == 0
        numPulsesField = 1;  % Fallback
    end
    
    % Convert dose to dose per pulse
    dosePerPulseGrid = fieldDose / numPulsesField;
    
    fprintf('    Max dose: %.4f Gy, Per-pulse max: %.6f Gy\n', ...
        max(fieldDose(:)), max(dosePerPulseGrid(:)));
    
    % Convert dose per pulse to initial pressure
    % p0 = D * Gamma * rho
    % Using gruneisen = 1 for now (placeholder)
    p0 = dosePerPulseGrid .* gruneisen_default .* density;
    
    fprintf('    Initial pressure range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));
    
    % Find dose extent to place sensor
    doseThreshold = 0.01 * max(fieldDose(:));  % 1% of max dose
    doseMask = fieldDose > doseThreshold;
    
    if ~any(doseMask(:))
        fprintf('    WARNING: No significant dose found, skipping field\n');
        continue;
    end
    
    % Find inferior extent of dose (lowest Y index with dose)
    [~, yIndices, ~] = ind2sub(size(doseMask), find(doseMask));
    minDoseY = min(yIndices);
    maxDoseY = max(yIndices);
    doseCenter = round((minDoseY + maxDoseY) / 2);
    
    % Find X-Z extent of dose for sensor placement
    [xIndices, ~, zIndices] = ind2sub(size(doseMask), find(doseMask));
    minDoseX = max(1, min(xIndices) - 5);
    maxDoseX = min(Nx, max(xIndices) + 5);
    minDoseZ = max(1, min(zIndices) - 5);
    maxDoseZ = min(Nz, max(zIndices) + 5);
    
    % Place planar sensor inferior to dose (transverse plane)
    sensorY = max(1, minDoseY - 10);  % 10 voxels inferior to dose
    
    fprintf('    Dose Y extent: [%d, %d], Sensor at Y=%d\n', minDoseY, maxDoseY, sensorY);
    
    % Create sensor mask (planar sensor in transverse plane)
    sensor = struct();
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(minDoseX:maxDoseX, sensorY, minDoseZ:maxDoseZ) = 1;
    
    numSensorPoints = sum(sensor.mask(:));
    fprintf('    Sensor points: %d\n', numSensorPoints);
    
    % Setup source (initial pressure distribution)
    source = struct();
    source.p0 = p0;
    
    % Run forward simulation
    fprintf('    Running forward simulation...\n');
    try
        sensorData = kspaceFirstOrder3D(kgrid, medium, source, sensor, inputArgs{:});
        fprintf('    Forward simulation complete\n');
    catch ME
        fprintf('    ERROR in forward simulation: %s\n', ME.message);
        continue;
    end
    
    % Time reversal reconstruction
    fprintf('    Running time reversal reconstruction...\n');
    
    try
        % Clear source and setup for time reversal
        source = rmfield(source, 'p0');
        source.p_mask = sensor.mask;
        source.p = fliplr(sensorData);  % Time-reversed sensor data
        source.p_mode = 'dirichlet';
        
        % Record final pressure distribution
        sensor.record = {'p_final'};
        
        % Run time reversal
        p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, inputArgs{:});
        
        % Extract reconstructed pressure
        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = p0_recon.p_final;
        else
            reconPressure = p0_recon;
        end
        
        % Apply positivity constraint
        reconPressure = max(reconPressure, 0);
        
        fprintf('    Time reversal complete\n');
        fprintf('    Reconstructed pressure range: [%.2e, %.2e] Pa\n', ...
            min(reconPressure(:)), max(reconPressure(:)));
        
    catch ME
        fprintf('    ERROR in time reversal: %s\n', ME.message);
        continue;
    end
    
    % Convert pressure back to dose per pulse
    % D = p0 / (Gamma * rho)
    % Avoid division by zero
    conversionFactor = gruneisen_default .* density;
    conversionFactor(conversionFactor == 0) = 1;  % Prevent division by zero
    
    reconDosePerPulse = reconPressure ./ conversionFactor;
    
    % Multiply by number of pulses to get total field dose
    reconFieldDose = reconDosePerPulse * numPulsesField;
    
    fprintf('    Reconstructed dose range: [%.4f, %.4f] Gy\n', ...
        min(reconFieldDose(:)), max(reconFieldDose(:)));
    
    % Store reconstructed field dose
    reconstructedFieldDoses{fieldIdx} = struct();
    reconstructedFieldDoses{fieldIdx}.dose = reconFieldDose;
    reconstructedFieldDoses{fieldIdx}.gantryAngle = fieldDoses{fieldIdx}.gantryAngle;
    reconstructedFieldDoses{fieldIdx}.numPulses = numPulsesField;
    reconstructedFieldDoses{fieldIdx}.originalMaxDose = max(fieldDose(:));
    reconstructedFieldDoses{fieldIdx}.reconMaxDose = max(reconFieldDose(:));
    
    % Accumulate total dose
    totalReconstructedDose = totalReconstructedDose + reconFieldDose;
    
    fprintf('    Field %d complete\n', fieldIdx);
end

% Calculate total original dose for comparison
totalOriginalDose = zeros(ctSize);
for fieldIdx = validFields'
    if ~isempty(fieldDoses{fieldIdx})
        fieldDose = fieldDoses{fieldIdx}.physicalDose;
        
        % Resample if necessary
        if ~isequal(size(fieldDose), ctSize)
            [nRows, nCols, nSlices] = size(fieldDose);
            [iDose, jDose, kDose] = ndgrid(1:nRows, 1:nCols, 1:nSlices);
            
            scale_i = nRows / Nx;
            scale_j = nCols / Ny;
            scale_k = nSlices / Nz;
            
            iQuery = 1 + (0:Nx-1) * scale_i;
            jQuery = 1 + (0:Ny-1) * scale_j;
            kQuery = 1 + (0:Nz-1) * scale_k;
            
            [iQ, jQ, kQ] = ndgrid(iQuery, jQuery, kQuery);
            fieldDose = interpn(iDose, jDose, kDose, fieldDose, iQ, jQ, kQ, 'linear', 0);
        end
        
        totalOriginalDose = totalOriginalDose + fieldDose;
    end
end

fprintf('\n  Total dose reconstruction complete\n');
fprintf('  - Original total max dose: %.4f Gy\n', max(totalOriginalDose(:)));
fprintf('  - Reconstructed total max dose: %.4f Gy\n', max(totalReconstructedDose(:)));

%% ======================== SAVE RESULTS ========================
fprintf('\n[8/9] Saving results and generating visualization...\n');

% Save all results
results = struct();
results.reconstructedFieldDoses = reconstructedFieldDoses;
results.totalReconstructedDose = totalReconstructedDose;
results.totalOriginalDose = totalOriginalDose;
results.density = density;
results.soundSpeed = soundSpeed;
results.gruneisen = gruneisen;
results.beamPulses = beamPulses;
results.dosePerPulse_Gy = dosePerPulse_Gy;
results.gridSize = ctSize;
results.patientID = patientID;
results.sessionName = sessionName;
results.usingResampledCT = usingResampledCT;

% Store resolution info
if isfield(ct, 'resolution')
    results.ctResolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
elseif exist('doseGrid', 'var')
    results.ctResolution = doseGrid.resolution;
end

% Store reference dose if available
if exist('referenceDose', 'var') && ~isempty(referenceDose)
    results.referenceDose = referenceDose;
end

save(fullfile(outputDir, 'kwave_reconstruction_results.mat'), 'results', '-v7.3');
fprintf('  - Results saved to: kwave_reconstruction_results.mat\n');

%% ======================== VISUALIZATION ========================
fprintf('\n[9/9] Generating visualization...\n');

% Find slice with maximum dose for display
[maxDose, maxIdx] = max(totalOriginalDose(:));
[maxX, maxY, maxZ] = ind2sub(ctSize, maxIdx);

% Create figure with comparison plots
fig = figure('Position', [100, 100, 1400, 900], 'Color', 'w');
sgtitle(sprintf('ETHOS Photoacoustic Dose Reconstruction: Patient %s', patientID), 'FontSize', 14);

% Subplot 1: Original dose - Axial (XY) slice
subplot(2, 3, 1);
imagesc(squeeze(totalOriginalDose(:, :, maxZ)));
colorbar;
title(sprintf('Original Dose - Axial (Z=%d)', maxZ));
xlabel('Y'); ylabel('X');
axis image;

% Subplot 2: Reconstructed dose - Axial slice
subplot(2, 3, 2);
imagesc(squeeze(totalReconstructedDose(:, :, maxZ)));
colorbar;
title(sprintf('Reconstructed Dose - Axial (Z=%d)', maxZ));
xlabel('Y'); ylabel('X');
axis image;

% Subplot 3: Difference - Axial slice
subplot(2, 3, 3);
doseDiff = totalReconstructedDose - totalOriginalDose;
imagesc(squeeze(doseDiff(:, :, maxZ)));
colorbar;
title('Difference (Recon - Orig)');
xlabel('Y'); ylabel('X');
axis image;

% Subplot 4: Original dose - Coronal (XZ) slice
subplot(2, 3, 4);
imagesc(squeeze(totalOriginalDose(:, maxY, :)));
colorbar;
title(sprintf('Original Dose - Coronal (Y=%d)', maxY));
xlabel('Z'); ylabel('X');
axis image;

% Subplot 5: Reconstructed dose - Coronal slice
subplot(2, 3, 5);
imagesc(squeeze(totalReconstructedDose(:, maxY, :)));
colorbar;
title(sprintf('Reconstructed Dose - Coronal (Y=%d)', maxY));
xlabel('Z'); ylabel('X');
axis image;

% Subplot 6: CT with dose overlay
subplot(2, 3, 6);
ctSlice = squeeze(ctCube(:, maxY, :));
doseSlice = squeeze(totalOriginalDose(:, maxY, :));
imagesc(ctSlice);
colormap(gca, 'gray');
hold on;
doseOverlay = doseSlice / max(doseSlice(:));
h = imagesc(doseOverlay);
set(h, 'AlphaData', doseOverlay * 0.5);
colormap(gca, 'hot');
title('CT with Dose Overlay');
xlabel('Z'); ylabel('X');
axis image;

% Save figure
saveas(fig, fullfile(outputDir, 'dose_reconstruction_comparison.png'));
saveas(fig, fullfile(outputDir, 'dose_reconstruction_comparison.fig'));
fprintf('  - Visualization saved\n');

% Create dose profile comparison
fig2 = figure('Position', [100, 100, 1000, 400], 'Color', 'w');

% Central axis profile (Y direction through max dose)
subplot(1, 2, 1);
origProfile = squeeze(totalOriginalDose(maxX, :, maxZ));
reconProfile = squeeze(totalReconstructedDose(maxX, :, maxZ));

if isfield(ct, 'resolution')
    yAxis = (1:Ny) * ct.resolution.y;
elseif exist('doseGrid', 'var')
    yAxis = (1:Ny) * doseGrid.resolution(1);
else
    yAxis = 1:Ny;
end

plot(yAxis, origProfile, 'b-', 'LineWidth', 2, 'DisplayName', 'Original');
hold on;
plot(yAxis, reconProfile, 'r--', 'LineWidth', 2, 'DisplayName', 'Reconstructed');
xlabel('Depth (mm)');
ylabel('Dose (Gy)');
title('Central Axis Dose Profile');
legend('Location', 'best');
grid on;

% Lateral profile (X direction through max dose)
subplot(1, 2, 2);
origProfileX = squeeze(totalOriginalDose(:, maxY, maxZ));
reconProfileX = squeeze(totalReconstructedDose(:, maxY, maxZ));

if isfield(ct, 'resolution')
    xAxis = (1:Nx) * ct.resolution.x;
elseif exist('doseGrid', 'var')
    xAxis = (1:Nx) * doseGrid.resolution(2);
else
    xAxis = 1:Nx;
end

plot(xAxis, origProfileX, 'b-', 'LineWidth', 2, 'DisplayName', 'Original');
hold on;
plot(xAxis, reconProfileX, 'r--', 'LineWidth', 2, 'DisplayName', 'Reconstructed');
xlabel('Position (mm)');
ylabel('Dose (Gy)');
title('Lateral Dose Profile');
legend('Location', 'best');
grid on;

sgtitle('Dose Profile Comparison', 'FontSize', 12);

saveas(fig2, fullfile(outputDir, 'dose_profiles.png'));
saveas(fig2, fullfile(outputDir, 'dose_profiles.fig'));

% Calculate reconstruction metrics
fprintf('\n  Reconstruction Metrics:\n');

% Only consider voxels with significant dose
doseThreshold = 0.01 * maxDose;
validMask = totalOriginalDose > doseThreshold;

if any(validMask(:))
    origValid = totalOriginalDose(validMask);
    reconValid = totalReconstructedDose(validMask);
    diffValid = doseDiff(validMask);
    
    % Mean absolute error
    mae = mean(abs(diffValid));
    fprintf('    Mean Absolute Error: %.4f Gy\n', mae);
    
    % Root mean square error
    rmse = sqrt(mean(diffValid.^2));
    fprintf('    Root Mean Square Error: %.4f Gy\n', rmse);
    
    % Mean relative error (where original > threshold)
    mre = mean(abs(diffValid) ./ origValid) * 100;
    fprintf('    Mean Relative Error: %.2f%%\n', mre);
    
    % Correlation coefficient
    corrCoef = corrcoef(origValid, reconValid);
    fprintf('    Correlation Coefficient: %.4f\n', corrCoef(1, 2));
    
    % Save metrics
    results.metrics.mae = mae;
    results.metrics.rmse = rmse;
    results.metrics.mre = mre;
    results.metrics.corrCoef = corrCoef(1, 2);
    
    save(fullfile(outputDir, 'kwave_reconstruction_results.mat'), 'results', '-v7.3');
end

%% ======================== SUMMARY ========================
fprintf('\n=======================================================\n');
fprintf('  Reconstruction Complete\n');
fprintf('=======================================================\n');
fprintf('  Patient: %s\n', patientID);
fprintf('  Session: %s\n', sessionName);
fprintf('  Fields processed: %d/%d\n', sum(~cellfun(@isempty, reconstructedFieldDoses)), numFields);
fprintf('  Total pulses: %d\n', totalPulses);
fprintf('  Original max dose: %.4f Gy\n', max(totalOriginalDose(:)));
fprintf('  Reconstructed max dose: %.4f Gy\n', max(totalReconstructedDose(:)));
if usingResampledCT
    fprintf('   Used CT resampled to RTDOSE grid\n');
else
    fprintf('   Used original CT grid\n');
end
fprintf('  Output directory: %s\n', outputDir);
fprintf('=======================================================\n\n');