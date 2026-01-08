%% ETHOS IMRT Segment-by-Segment Dose Calculator (FIXED VERSION)
% Purpose: Calculate individual segment doses from ETHOS exported IMRT data
% Uses MATRAD for dose calculation
%
% FIXES IN THIS VERSION:
%   - Proper MLC position inheritance across control points
%   - Dual-layer MLC handling (combine MLCX1 and MLCX2)
%   - Memory management: clear segment doses after accumulating to beam dose
%   - Better debugging output to verify MLC aperture application
%
% KEY FEATURES:
%   - Calculates dose for each control point/segment within each beam
%   - IMRT plans have multiple segments per beam with different MLC configurations
%   - Each segment contributes a weighted portion of the total beam dose
%
% Author: Generated for ETHOS dose analysis
% Date: 2025

clear; clc; close all;

%% ==================== CONFIGURATION ====================
% Patient and session arrays (expandable for batch processing)
ids = {'1194203'};
sessions = {'Session_1'};

% Base directory
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
matradPath = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

%% ==================== OPTIMIZATION PARAMETERS ====================
% CT Downsampling Factor
% 1 = original resolution (slowest, most accurate)
% 2 = half resolution in each dimension (8x faster)
ctDownsampleFactor = 2;  % ADJUST THIS TO SPEED UP CALCULATIONS

% Caching Options
enableCaching = true;
cacheDir = fullfile(wd, 'Cache');

% Verbose output
verboseOutput = true;

% Debug MLC (print details about MLC aperture checking)
debugMLC = true;  % Set to true to see MLC position details

%% ==================== INITIALIZATION ====================
fprintf('==========================================================\n');
fprintf('  ETHOS IMRT Segment-by-Segment Dose Calculator\n');
fprintf('  FIXED VERSION - Proper MLC Handling\n');
fprintf('==========================================================\n\n');

fprintf('Optimization Settings:\n');
fprintf('  - CT Downsample Factor: %d (%.0fx speedup estimate)\n', ...
    ctDownsampleFactor, ctDownsampleFactor^3);
fprintf('  - Caching Enabled: %s\n', string(enableCaching));
fprintf('  - Cache Directory: %s\n', cacheDir);
fprintf('\n');

% Verify paths exist
fprintf('Verifying paths...\n');
if ~exist(wd, 'dir')
    error('Working directory does not exist: %s', wd);
else
    fprintf('  - Working directory OK: %s\n', wd);
end

if ~exist(matradPath, 'dir')
    error('matRad path does not exist: %s', matradPath);
else
    fprintf('  - matRad path OK: %s\n', matradPath);
end

% Create cache directory
if enableCaching && ~exist(cacheDir, 'dir')
    mkdir(cacheDir);
    fprintf('  - Created cache directory: %s\n', cacheDir);
end

% Add MATRAD to path
addpath(genpath(matradPath));
fprintf('  - matRad added to MATLAB path\n');

% Initialize MATRAD
fprintf('Initializing MATRAD...\n');
try
    matRad_rc;
    fprintf('  - MATRAD initialized successfully\n');
catch
    try
        matRad_cfg = MatRad_Config.instance();
        fprintf('  - MATRAD configured successfully\n');
    catch
        fprintf('  - Direct initialization not available, continuing...\n');
    end
end

%% ==================== MAIN PROCESSING LOOP ====================
for idxID = 1:length(ids)
    for idxSession = 1:length(sessions)
        
        currentID = ids{idxID};
        currentSession = sessions{idxSession};
        
        fprintf('\n========================================\n');
        fprintf('Processing Patient: %s, Session: %s\n', currentID, currentSession);
        fprintf('========================================\n');
        
        % Define paths
        rawwd = fullfile(wd, 'EthosExports', currentID, 'Pancreas', currentSession);
        dicomPath = fullfile(rawwd, 'sct');
        outputPath = fullfile(wd, 'SegmentDoses', currentID, currentSession);
        
        % Patient-specific cache directory (includes downsample factor)
        patientCacheDir = fullfile(cacheDir, currentID, currentSession, ...
            sprintf('ds%d', ctDownsampleFactor));
        if enableCaching && ~exist(patientCacheDir, 'dir')
            mkdir(patientCacheDir);
        end
        
        dicomCacheDir = fullfile(cacheDir, currentID, currentSession);
        if ~exist(dicomCacheDir, 'dir')
            mkdir(dicomCacheDir);
        end
        
        % Create output directory
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        
        %% Step 1: Import DICOM data (with caching)
        fprintf('\n[1/8] Importing DICOM data...\n');
        
        dicomCacheFile = fullfile(dicomCacheDir, 'dicom_import.mat');
        
        if enableCaching && exist(dicomCacheFile, 'file')
            fprintf('  Loading cached DICOM import...\n');
            load(dicomCacheFile, 'ct', 'cst', 'pln');
            fprintf('  - Loaded from cache\n');
        else
            try
                fprintf('  Creating matRad_DicomImporter object...\n');
                importer = matRad_DicomImporter(dicomPath);
                fprintf('  Importing DICOM files...\n');
                importer.matRad_importDicom();
                fprintf('  - Data imported successfully\n');
                
                % Cache the import (at original resolution)
                if enableCaching
                    save(dicomCacheFile, 'ct', 'cst', 'pln', '-v7.3');
                    fprintf('  - Cached DICOM import\n');
                end
                
            catch ME
                fprintf('Error importing DICOM data: %s\n', ME.message);
                continue;
            end
        end
        
        % Store original CT for reference
        ct_original = ct;
        cst_original = cst;
        originalCTSize = ct.cubeDim;
        
        fprintf('  - Original CT dimensions: %d x %d x %d\n', ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));
        fprintf('  - Original CT resolution: %.2f x %.2f x %.2f mm\n', ...
            ct.resolution.x, ct.resolution.y, ct.resolution.z);
        
        %% Step 1b: Downsample CT if requested
        if ctDownsampleFactor > 1
            fprintf('\n  Applying CT downsampling (factor=%d)...\n', ctDownsampleFactor);
            
            dsCacheFile = fullfile(patientCacheDir, 'ct_downsampled.mat');
            
            if enableCaching && exist(dsCacheFile, 'file')
                fprintf('    Loading cached downsampled CT...\n');
                load(dsCacheFile, 'ct', 'cst');
                fprintf('    Loaded from cache: [%d,%d,%d]\n', ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));
            else
                [ct, cst] = downsampleCTComplete(ct_original, cst_original, ctDownsampleFactor);
                
                if enableCaching
                    save(dsCacheFile, 'ct', 'cst', '-v7.3');
                    fprintf('    Cached downsampled CT\n');
                end
            end
        end
        
        %% Step 2: Load reference RTDOSE and spatial coordinates
        fprintf('\n[2/8] Loading reference RTDOSE and extracting spatial coordinates...\n');
        
        ctSpatial = struct();
        doseSpatial = struct();
        
        try
            % Extract CT spatial coordinates from ORIGINAL CT files
            ctFiles = dir(fullfile(dicomPath, 'CT*.dcm'));
            if isempty(ctFiles)
                ctFiles = dir(fullfile(dicomPath, '*CT*.dcm'));
            end
            
            if ~isempty(ctFiles)
                ctInfo = dicominfo(fullfile(ctFiles(1).folder, ctFiles(1).name));
                
                ctSpatial.ImagePositionPatient = ctInfo.ImagePositionPatient;
                ctSpatial.PixelSpacing = ctInfo.PixelSpacing;
                ctSpatial.Rows = ctInfo.Rows;
                ctSpatial.Columns = ctInfo.Columns;
                
                slicePositions = zeros(length(ctFiles), 1);
                for i = 1:length(ctFiles)
                    tempInfo = dicominfo(fullfile(ctFiles(i).folder, ctFiles(i).name));
                    slicePositions(i) = tempInfo.ImagePositionPatient(3);
                end
                slicePositions = sort(slicePositions);
                
                ctSpatial.SlicePositions = slicePositions;
                ctSpatial.SliceThickness = abs(slicePositions(2) - slicePositions(1));
                ctSpatial.NumSlices = length(slicePositions);
                
                fprintf('    CT spatial extent: X[%.1f,%.1f], Y[%.1f,%.1f], Z[%.1f,%.1f] mm\n', ...
                    ctSpatial.ImagePositionPatient(1), ...
                    ctSpatial.ImagePositionPatient(1) + (ctSpatial.Columns-1)*ctSpatial.PixelSpacing(2), ...
                    ctSpatial.ImagePositionPatient(2), ...
                    ctSpatial.ImagePositionPatient(2) + (ctSpatial.Rows-1)*ctSpatial.PixelSpacing(1), ...
                    min(slicePositions), max(slicePositions));
            end
            
            % Load RTDOSE
            rtdoseFile = dir(fullfile(dicomPath, 'RD*.dcm'));
            if isempty(rtdoseFile)
                rtdoseFile = dir(fullfile(dicomPath, '*RTDOSE*.dcm'));
            end
            
            if ~isempty(rtdoseFile)
                rtdoseInfo = dicominfo(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name));
                referenceDose = squeeze(double(dicomread(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name))));
                referenceDose = referenceDose * rtdoseInfo.DoseGridScaling;
                
                fprintf('  - Reference dose: %d x %d x %d, Max: %.2f Gy\n', ...
                    size(referenceDose,1), size(referenceDose,2), size(referenceDose,3), max(referenceDose(:)));
                
                if isfield(rtdoseInfo, 'DoseSummationType') && strcmp(rtdoseInfo.DoseSummationType, 'PLAN')
                    fprintf('  WARNING: PLAN dose detected, dividing by 10 fractions\n');
                    referenceDose = referenceDose / 10;
                end
                
                doseSpatial.ImagePositionPatient = rtdoseInfo.ImagePositionPatient;
                doseSpatial.PixelSpacing = rtdoseInfo.PixelSpacing;
                doseSpatial.Rows = rtdoseInfo.Rows;
                doseSpatial.Columns = rtdoseInfo.Columns;
                
                if isfield(rtdoseInfo, 'GridFrameOffsetVector') && length(rtdoseInfo.GridFrameOffsetVector) > 1
                    doseSpatial.GridFrameOffsetVector = rtdoseInfo.GridFrameOffsetVector;
                    doseSpatial.SliceThickness = abs(rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1));
                    doseSpatial.NumSlices = length(rtdoseInfo.GridFrameOffsetVector);
                else
                    doseSpatial.SliceThickness = 2.5;
                    doseSpatial.NumSlices = size(referenceDose, 3);
                    doseSpatial.GridFrameOffsetVector = (0:doseSpatial.NumSlices-1) * doseSpatial.SliceThickness;
                end
                
                doseGrid.resolution = [doseSpatial.PixelSpacing(1), doseSpatial.PixelSpacing(2), doseSpatial.SliceThickness];
                doseGrid.dimensions = [doseSpatial.Rows, doseSpatial.Columns, doseSpatial.NumSlices];
                doseGrid.ImagePositionPatient = doseSpatial.ImagePositionPatient;
                
            else
                fprintf('  WARNING: No RTDOSE file found\n');
                referenceDose = [];
                doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
                doseGrid.dimensions = ct.cubeDim;
            end
            
        catch ME
            fprintf('WARNING: Could not load RTDOSE: %s\n', ME.message);
            referenceDose = [];
            doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
            doseGrid.dimensions = ct.cubeDim;
        end
        
        %% Step 3: Extract segment information from RTPLAN with PROPER MLC HANDLING
        fprintf('\n[3/8] Extracting segment information from RTPLAN (with MLC inheritance)...\n');
        
        % Use a separate cache file for the fixed segment data
        segmentCacheFile = fullfile(dicomCacheDir, 'segmentData_fixed.mat');
        
        if enableCaching && exist(segmentCacheFile, 'file')
            fprintf('  Loading cached segment data...\n');
            load(segmentCacheFile, 'segmentData');
            fprintf('  - Loaded %d beams, %d total segments from cache\n', ...
                segmentData.numBeams, segmentData.totalSegments);
        else
            segmentData = struct();
            segmentData.beams = {};
            
            try
                rtplanFile = dir(fullfile(dicomPath, 'RP*.dcm'));
                if isempty(rtplanFile)
                    rtplanFile = dir(fullfile(dicomPath, '*RTPLAN*.dcm'));
                end
                
                if isempty(rtplanFile)
                    error('No RTPLAN file found');
                end
                
                rtplanInfo = dicominfo(fullfile(rtplanFile(1).folder, rtplanFile(1).name));
                
                % Get beam metersets
                beamMetersets = [];
                numFractions = 1;
                if isfield(rtplanInfo, 'FractionGroupSequence')
                    fg = rtplanInfo.FractionGroupSequence.Item_1;
                    if isfield(fg, 'NumberOfFractionsPlanned')
                        numFractions = fg.NumberOfFractionsPlanned;
                    end
                    if isfield(fg, 'ReferencedBeamSequence')
                        numRefBeams = length(fieldnames(fg.ReferencedBeamSequence));
                        beamMetersets = zeros(numRefBeams, 1);
                        for i = 1:numRefBeams
                            refBeam = fg.ReferencedBeamSequence.(sprintf('Item_%d', i));
                            if isfield(refBeam, 'BeamMeterset')
                                beamMetersets(i) = refBeam.BeamMeterset;
                            end
                        end
                    end
                end
                
                numBeams = length(fieldnames(rtplanInfo.BeamSequence));
                totalSegments = 0;
                
                for beamIdx = 1:numBeams
                    beamField = sprintf('Item_%d', beamIdx);
                    beam = rtplanInfo.BeamSequence.(beamField);
                    
                    beamData = struct();
                    beamData.beamIdx = beamIdx;
                    beamData.beamName = '';
                    beamData.gantryAngle = 0;
                    beamData.couchAngle = 0;
                    beamData.collimatorAngle = 0;
                    beamData.beamMeterset = 0;
                    beamData.segments = {};
                    beamData.leafBoundaries = [];
                    beamData.numLeafPairs = 0;
                    beamData.isDualLayerMLC = false;
                    
                    if isfield(beam, 'BeamName')
                        beamData.beamName = beam.BeamName;
                    else
                        beamData.beamName = sprintf('Beam_%d', beamIdx);
                    end
                    
                    if beamIdx <= length(beamMetersets)
                        beamData.beamMeterset = beamMetersets(beamIdx);
                    end
                    
                    % Get leaf boundaries from beam limiting device sequence
                    % Also detect dual-layer MLC
                    mlcDeviceNames = {};
                    if isfield(beam, 'BeamLimitingDeviceSequence')
                        numDevices = length(fieldnames(beam.BeamLimitingDeviceSequence));
                        for devIdx = 1:numDevices
                            device = beam.BeamLimitingDeviceSequence.(sprintf('Item_%d', devIdx));
                            if isfield(device, 'RTBeamLimitingDeviceType')
                                devType = device.RTBeamLimitingDeviceType;
                                if contains(devType, 'MLC', 'IgnoreCase', true)
                                    mlcDeviceNames{end+1} = devType;
                                    if isfield(device, 'LeafPositionBoundaries')
                                        beamData.leafBoundaries = device.LeafPositionBoundaries;
                                        beamData.numLeafPairs = length(beamData.leafBoundaries) - 1;
                                    end
                                end
                            end
                        end
                    end
                    
                    beamData.isDualLayerMLC = length(mlcDeviceNames) >= 2;
                    beamData.mlcDeviceNames = mlcDeviceNames;
                    
                    if beamData.isDualLayerMLC && debugMLC
                        fprintf('    Beam %d: Dual-layer MLC detected: %s\n', beamIdx, strjoin(mlcDeviceNames, ', '));
                    end
                    
                    % Process control points with PROPER INHERITANCE
                    if isfield(beam, 'ControlPointSequence')
                        numCP = length(fieldnames(beam.ControlPointSequence));
                        controlPoints = cell(numCP, 1);
                        
                        % Variables to track inherited values
                        inheritedMlcLayer1 = [];
                        inheritedMlcLayer2 = [];
                        inheritedJawX = [];
                        inheritedJawY = [];
                        inheritedGantryAngle = 0;
                        
                        for cpIdx = 1:numCP
                            cp = beam.ControlPointSequence.(sprintf('Item_%d', cpIdx));
                            
                            cpData = struct();
                            cpData.cpIdx = cpIdx;
                            cpData.cumulativeMetersetWeight = 0;
                            cpData.gantryAngle = [];
                            cpData.mlcLayer1 = [];  % First MLC layer
                            cpData.mlcLayer2 = [];  % Second MLC layer (if dual-layer)
                            cpData.mlcCombined = []; % Combined single-layer equivalent
                            cpData.jawX = [];
                            cpData.jawY = [];
                            
                            if isfield(cp, 'CumulativeMetersetWeight')
                                cpData.cumulativeMetersetWeight = cp.CumulativeMetersetWeight;
                            end
                            
                            % Gantry angle - inherit if not specified
                            if isfield(cp, 'GantryAngle')
                                inheritedGantryAngle = cp.GantryAngle;
                                beamData.gantryAngle = cp.GantryAngle;
                            end
                            cpData.gantryAngle = inheritedGantryAngle;
                            
                            if isfield(cp, 'PatientSupportAngle')
                                beamData.couchAngle = cp.PatientSupportAngle;
                            end
                            
                            % Extract MLC and jaw positions
                            if isfield(cp, 'BeamLimitingDevicePositionSequence')
                                numDevices = length(fieldnames(cp.BeamLimitingDevicePositionSequence));
                                
                                mlcLayerIdx = 0;  % Counter for MLC layers in this control point
                                
                                for devIdx = 1:numDevices
                                    device = cp.BeamLimitingDevicePositionSequence.(sprintf('Item_%d', devIdx));
                                    if ~isfield(device, 'RTBeamLimitingDeviceType') || ~isfield(device, 'LeafJawPositions')
                                        continue;
                                    end
                                    
                                    deviceType = device.RTBeamLimitingDeviceType;
                                    positions = device.LeafJawPositions;
                                    
                                    if contains(deviceType, 'ASYMX', 'IgnoreCase', true) || strcmp(deviceType, 'X')
                                        inheritedJawX = positions;
                                    elseif contains(deviceType, 'ASYMY', 'IgnoreCase', true) || strcmp(deviceType, 'Y')
                                        inheritedJawY = positions;
                                    elseif contains(deviceType, 'MLC', 'IgnoreCase', true)
                                        % Handle dual-layer MLC
                                        mlcLayerIdx = mlcLayerIdx + 1;
                                        if mlcLayerIdx == 1
                                            inheritedMlcLayer1 = positions;
                                        else
                                            inheritedMlcLayer2 = positions;
                                        end
                                    end
                                end
                            end
                            
                            % Apply inherited values
                            cpData.jawX = inheritedJawX;
                            cpData.jawY = inheritedJawY;
                            cpData.mlcLayer1 = inheritedMlcLayer1;
                            cpData.mlcLayer2 = inheritedMlcLayer2;
                            
                            % Combine dual-layer MLC into single-layer equivalent
                            if ~isempty(inheritedMlcLayer1)
                                if ~isempty(inheritedMlcLayer2) && length(inheritedMlcLayer1) == length(inheritedMlcLayer2)
                                    % Dual-layer: combine using max(left)/min(right) rule
                                    numLeaves = length(inheritedMlcLayer1) / 2;
                                    combined = zeros(size(inheritedMlcLayer1));
                                    
                                    % Left bank (A-bank): Take MAX (most restrictive for negative values)
                                    for i = 1:numLeaves
                                        combined(i) = max(inheritedMlcLayer1(i), inheritedMlcLayer2(i));
                                    end
                                    
                                    % Right bank (B-bank): Take MIN (most restrictive for positive values)
                                    for i = 1:numLeaves
                                        idx = numLeaves + i;
                                        combined(idx) = min(inheritedMlcLayer1(idx), inheritedMlcLayer2(idx));
                                    end
                                    
                                    cpData.mlcCombined = combined;
                                else
                                    % Single layer or mismatched: use layer 1
                                    cpData.mlcCombined = inheritedMlcLayer1;
                                end
                            end
                            
                            controlPoints{cpIdx} = cpData;
                        end
                        
                        % Create segments using combined MLC positions
                        numSegments = numCP - 1;
                        validSegmentCount = 0;
                        
                        for segIdx = 1:numSegments
                            cpStart = controlPoints{segIdx};
                            cpEnd = controlPoints{segIdx + 1};
                            
                            segmentWeight = cpEnd.cumulativeMetersetWeight - cpStart.cumulativeMetersetWeight;
                            
                            % Skip segments with zero weight
                            if segmentWeight <= 0
                                continue;
                            end
                            
                            validSegmentCount = validSegmentCount + 1;
                            
                            segment = struct();
                            segment.segmentIdx = segIdx;
                            segment.beamIdx = beamIdx;
                            segment.segmentWeight = segmentWeight;
                            segment.segmentMeterset = segmentWeight * beamData.beamMeterset;
                            segment.gantryAngle = cpStart.gantryAngle;
                            segment.mlcPositions = cpStart.mlcCombined;  % Use combined MLC
                            segment.jawX = cpStart.jawX;
                            segment.jawY = cpStart.jawY;
                            
                            beamData.segments{end+1} = segment;
                            totalSegments = totalSegments + 1;
                        end
                        
                        % Debug output for first beam
                        if beamIdx == 1 && debugMLC
                            fprintf('    Beam %d first segment MLC check:\n', beamIdx);
                            if ~isempty(beamData.segments) && ~isempty(beamData.segments{1}.mlcPositions)
                                mlc = beamData.segments{1}.mlcPositions;
                                numLeaves = length(mlc) / 2;
                                fprintf('      MLC positions available: YES (%d leaf pairs)\n', numLeaves);
                                fprintf('      Left bank range: [%.1f, %.1f] mm\n', min(mlc(1:numLeaves)), max(mlc(1:numLeaves)));
                                fprintf('      Right bank range: [%.1f, %.1f] mm\n', min(mlc(numLeaves+1:end)), max(mlc(numLeaves+1:end)));
                                
                                % Show aperture width for middle leaves
                                midLeaf = round(numLeaves/2);
                                aperture = mlc(numLeaves + midLeaf) - mlc(midLeaf);
                                fprintf('      Mid-leaf aperture: %.1f mm\n', aperture);
                            else
                                fprintf('      MLC positions available: NO\n');
                            end
                            if ~isempty(beamData.segments{1}.jawX)
                                fprintf('      Jaw X: [%.1f, %.1f] mm\n', beamData.segments{1}.jawX(1), beamData.segments{1}.jawX(2));
                            end
                            if ~isempty(beamData.segments{1}.jawY)
                                fprintf('      Jaw Y: [%.1f, %.1f] mm\n', beamData.segments{1}.jawY(1), beamData.segments{1}.jawY(2));
                            end
                        end
                    end
                    
                    segmentData.beams{beamIdx} = beamData;
                    fprintf('    Beam %d: %s, %.1f MU, %d segments\n', ...
                        beamIdx, beamData.beamName, beamData.beamMeterset, length(beamData.segments));
                end
                
                segmentData.totalSegments = totalSegments;
                segmentData.numBeams = numBeams;
                segmentData.numFractions = numFractions;
                
                % Cache segment data
                if enableCaching
                    save(segmentCacheFile, 'segmentData');
                    fprintf('  - Cached segment data\n');
                end
                
            catch ME
                fprintf('ERROR extracting segment data: %s\n', ME.message);
                for k = 1:length(ME.stack)
                    fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                end
                continue;
            end
        end
        
        %% Step 4: Configure dose calculation
        fprintf('\n[4/8] Configuring dose calculation...\n');
        
        machineDir = fullfile(matradPath, 'basedata');
        if exist(machineDir, 'dir')
            machineFiles = dir(fullfile(machineDir, 'photons*.mat'));
            if ~isempty(machineFiles)
                genericMachines = machineFiles(contains({machineFiles.name}, 'Generic', 'IgnoreCase', true));
                if ~isempty(genericMachines)
                    [~, selectedMachine, ~] = fileparts(genericMachines(1).name);
                else
                    [~, selectedMachine, ~] = fileparts(machineFiles(1).name);
                end
                pln.machine = selectedMachine;
            else
                pln.machine = 'Generic';
            end
        else
            pln.machine = 'Generic';
        end
        
        pln.propStf.bixelWidth = 5;
        pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x;
        pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y;
        pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z;
        
        if strcmp(pln.radiationMode, 'photons')
            pln.propDoseCalc.engine = 'pencilBeam';
        else
            pln.propDoseCalc.engine = 'matRad_pencilBeam';
        end
        
        fprintf('  - Machine: %s\n', pln.machine);
        fprintf('  - Dose engine: %s\n', pln.propDoseCalc.engine);
        fprintf('  - Dose grid resolution: [%.2f, %.2f, %.2f] mm\n', ...
            ct.resolution.x, ct.resolution.y, ct.resolution.z);
        fprintf('  - CT dimensions: [%d, %d, %d]\n', ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));
        
        %% Step 5: Generate steering information
        fprintf('\n[5/8] Generating steering information...\n');
        
        stfCacheFile = fullfile(patientCacheDir, 'stf.mat');
        
        if enableCaching && exist(stfCacheFile, 'file')
            fprintf('  Loading cached steering file...\n');
            load(stfCacheFile, 'stf');
            fprintf('  - Loaded %d beams from cache\n', length(stf));
        else
            try
                fprintf('  Generating stf (this may take a moment)...\n');
                stf = matRad_generateStf(ct, cst, pln);
                fprintf('  - Generated steering file for %d beams\n', length(stf));
                
                if enableCaching
                    save(stfCacheFile, 'stf');
                    fprintf('  - Cached steering file\n');
                end
            catch ME
                fprintf('Error generating steering file: %s\n', ME.message);
                for k = 1:length(ME.stack)
                    fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                end
                continue;
            end
        end
        
        for i = 1:length(stf)
            fprintf('    Beam %d: Gantry=%.1f, Rays=%d\n', i, stf(i).gantryAngle, stf(i).numOfRays);
        end
        
        %% Step 6: Calculate segment-by-segment doses with PROPER MLC HANDLING
        fprintf('\n[6/8] Calculating segment-by-segment doses...\n');
        
        calculatedGridSize = ct.cubeDim;
        totalDose = zeros(calculatedGridSize);
        beamDoses = cell(length(stf), 1);
        
        totalCalcTime = 0;
        cachedBeams = 0;
        calculatedBeams = 0;
        
        % Statistics for MLC aperture usage
        totalRaysChecked = 0;
        totalRaysInAperture = 0;
        
        for beamIdx = 1:length(stf)
            beamData = segmentData.beams{beamIdx};
            numSegments = length(beamData.segments);
            
            fprintf('\n  Beam %d/%d (%s): %d segments, %.1f MU\n', ...
                beamIdx, length(stf), beamData.beamName, numSegments, beamData.beamMeterset);
            
            % Check for cached beam dose
            beamDoseCacheFile = fullfile(patientCacheDir, sprintf('beam_%02d_dose_fixed.mat', beamIdx));
            
            if enableCaching && exist(beamDoseCacheFile, 'file')
                fprintf('    Loading cached beam dose...\n');
                load(beamDoseCacheFile, 'beamDoseResult');
                
                beamDoses{beamIdx} = beamDoseResult;
                totalDose = totalDose + beamDoseResult.physicalDose;
                cachedBeams = cachedBeams + 1;
                
                fprintf('    Loaded from cache (Max: %.4f Gy)\n', beamDoseResult.maxDose);
                continue;
            end
            
            % Calculate dij for this beam
            dijCacheFile = fullfile(patientCacheDir, sprintf('beam_%02d_dij.mat', beamIdx));
            
            plnSingle = pln;
            plnSingle.propStf.numOfBeams = 1;
            plnSingle.propStf.isoCenter = stf(beamIdx).isoCenter;
            stfSingle = stf(beamIdx);
            
            if enableCaching && exist(dijCacheFile, 'file')
                fprintf('    Loading cached dij matrix...\n');
                load(dijCacheFile, 'dij');
                fprintf('    Loaded dij (%d bixels) from cache\n', dij.totalNumOfBixels);
            else
                try
                    fprintf('    Calculating dij matrix...\n');
                    tic;
                    dij = matRad_calcDoseInfluence(ct, cst, stfSingle, plnSingle);
                    dijTime = toc;
                    totalCalcTime = totalCalcTime + dijTime;
                    fprintf('    Calculated dij (%d bixels) in %.1f sec\n', dij.totalNumOfBixels, dijTime);
                    
                    if enableCaching
                        fprintf('    Caching dij matrix...\n');
                        save(dijCacheFile, 'dij', '-v7.3');
                        fprintf('    Cached dij matrix\n');
                    end
                catch ME
                    fprintf('    ERROR calculating dij: %s\n', ME.message);
                    for k = 1:length(ME.stack)
                        fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                    end
                    continue;
                end
            end
            
            % Initialize beam dose accumulator
            beamDoseAccum = zeros(calculatedGridSize);
            beamRaysInAperture = 0;
            beamRaysTotal = 0;
            
            % Process each segment
            for segIdx = 1:numSegments
                segment = beamData.segments{segIdx};
                
                if segment.segmentWeight <= 0
                    continue;
                end
                
                % Create weight vector based on MLC aperture
                w = zeros(dij.totalNumOfBixels, 1);
                
                mlcPos = segment.mlcPositions;
                jawX = segment.jawX;
                jawY = segment.jawY;
                
                segRaysInAperture = 0;
                
                if ~isempty(mlcPos) && ~isempty(beamData.leafBoundaries)
                    numLeafPairs = beamData.numLeafPairs;
                    leafBoundaries = beamData.leafBoundaries;
                    
                    % Handle potential size mismatch
                    if length(mlcPos) < 2 * numLeafPairs
                        numLeafPairs = floor(length(mlcPos) / 2);
                    end
                    
                    leftBank = mlcPos(1:numLeafPairs);
                    rightBank = mlcPos(numLeafPairs+1:2*numLeafPairs);
                    
                    numRays = stfSingle.numOfRays;
                    bixelIdx = 0;
                    
                    for rayIdx = 1:numRays
                        ray = stfSingle.ray(rayIdx);
                        rayPosX = ray.rayPos_bev(1);
                        rayPosY = ray.rayPos_bev(2);
                        
                        % Check jaw constraints
                        inJawX = isempty(jawX) || (length(jawX) >= 2 && rayPosX >= jawX(1) && rayPosX <= jawX(2));
                        inJawY = isempty(jawY) || (length(jawY) >= 2 && rayPosY >= jawY(1) && rayPosY <= jawY(2));
                        
                        % Check MLC aperture
                        inMLC = false;
                        if inJawX && inJawY
                            % Find which leaf pair this ray falls into
                            for leafIdx = 1:numLeafPairs
                                if leafIdx <= length(leafBoundaries)-1
                                    if rayPosY >= leafBoundaries(leafIdx) && rayPosY < leafBoundaries(leafIdx + 1)
                                        % Check if ray is within the aperture for this leaf pair
                                        if rayPosX >= leftBank(leafIdx) && rayPosX <= rightBank(leafIdx)
                                            inMLC = true;
                                        end
                                        break;  % Found the correct leaf pair, no need to continue
                                    end
                                end
                            end
                        end
                        
                        numEnergies = length(ray.energy);
                        for energyIdx = 1:numEnergies
                            bixelIdx = bixelIdx + 1;
                            if inMLC
                                w(bixelIdx) = segment.segmentWeight * beamData.beamMeterset;
                                segRaysInAperture = segRaysInAperture + 1;
                            end
                        end
                    end
                    
                    beamRaysTotal = beamRaysTotal + dij.totalNumOfBixels;
                    beamRaysInAperture = beamRaysInAperture + segRaysInAperture;
                else
                    % WARNING: No MLC data - this should NOT happen with fixed code
                    fprintf('      WARNING: Seg %d has no MLC data!\n', segment.segmentIdx);
                    % DO NOT use fallback - this would cause overdose
                    % Instead, skip this segment
                    continue;
                end
                
                % Calculate dose using direct matrix multiplication
                try
                    if any(w > 0)
                        tic;
                        segmentDose = reshape(full(dij.physicalDose{1} * w), calculatedGridSize);
                        calcTime = toc;
                        totalCalcTime = totalCalcTime + calcTime;
                        
                        maxDose = max(segmentDose(:));
                        beamDoseAccum = beamDoseAccum + segmentDose;
                        
                        if verboseOutput
                            openFrac = segRaysInAperture / dij.totalNumOfBixels * 100;
                            fprintf('      Seg %d: Max=%.4f Gy, Open=%.1f%% (%.2fs)\n', ...
                                segment.segmentIdx, maxDose, openFrac, calcTime);
                        end
                    else
                        if verboseOutput
                            fprintf('      Seg %d: All rays blocked by MLC\n', segment.segmentIdx);
                        end
                    end
                catch ME
                    fprintf('      Seg %d ERROR: %s\n', segment.segmentIdx, ME.message);
                end
                
                % MEMORY MANAGEMENT: Clear segment-level data
                clear segmentDose w;
            end
            
            totalRaysChecked = totalRaysChecked + beamRaysTotal;
            totalRaysInAperture = totalRaysInAperture + beamRaysInAperture;
            
            % Store beam dose
            beamDoseResult = struct();
            beamDoseResult.beamIdx = beamIdx;
            beamDoseResult.beamName = beamData.beamName;
            beamDoseResult.gantryAngle = beamData.gantryAngle;
            beamDoseResult.couchAngle = beamData.couchAngle;
            beamDoseResult.beamMeterset = beamData.beamMeterset;
            beamDoseResult.numSegments = numSegments;
            beamDoseResult.physicalDose = beamDoseAccum;
            beamDoseResult.maxDose = max(beamDoseAccum(:));
            
            if beamRaysTotal > 0
                beamDoseResult.apertureOpenFraction = beamRaysInAperture / beamRaysTotal;
            else
                beamDoseResult.apertureOpenFraction = 0;
            end
            
            beamDoses{beamIdx} = beamDoseResult;
            totalDose = totalDose + beamDoseAccum;
            calculatedBeams = calculatedBeams + 1;
            
            % Cache complete beam result
            if enableCaching
                save(beamDoseCacheFile, 'beamDoseResult', '-v7.3');
            end
            
            fprintf('    Beam %d complete: Max=%.4f Gy, Avg open: %.1f%%\n', ...
                beamIdx, beamDoseResult.maxDose, beamDoseResult.apertureOpenFraction * 100);
            
            % MEMORY MANAGEMENT: Clear dij and beam accumulator
            clear dij beamDoseAccum;
        end
        
        fprintf('\n  Summary:\n');
        fprintf('    Beams processed: %d (cached: %d, calculated: %d)\n', ...
            length(stf), cachedBeams, calculatedBeams);
        fprintf('    Total calculation time: %.1f sec\n', totalCalcTime);
        fprintf('    Total dose max: %.4f Gy\n', max(totalDose(:)));
        if totalRaysChecked > 0
            fprintf('    Overall aperture open fraction: %.1f%%\n', totalRaysInAperture / totalRaysChecked * 100);
        end
        
        %% Step 7: Upsample doses back to original resolution and resample to RTDOSE grid
        fprintf('\n[7/8] Resampling doses...\n');
        
        % First upsample back to original CT resolution if we downsampled
        if ctDownsampleFactor > 1
            fprintf('  Upsampling doses to original CT resolution...\n');
            totalDose = imresize3(totalDose, originalCTSize, 'linear');
            
            for beamIdx = 1:length(beamDoses)
                if ~isempty(beamDoses{beamIdx})
                    beamDoses{beamIdx}.physicalDose = imresize3(beamDoses{beamIdx}.physicalDose, originalCTSize, 'linear');
                    beamDoses{beamIdx}.maxDose = max(beamDoses{beamIdx}.physicalDose(:));
                end
            end
            
            fprintf('  - Upsampled total and beam doses to [%d,%d,%d]\n', originalCTSize(1), originalCTSize(2), originalCTSize(3));
        end
        
        % Now resample to RTDOSE grid
        if ~isempty(referenceDose) && ~isempty(ctSpatial.ImagePositionPatient) && ~isempty(doseSpatial.ImagePositionPatient)
            fprintf('  Resampling to RTDOSE grid...\n');
            
            ct_x = double(ctSpatial.ImagePositionPatient(1)) + double((0:ctSpatial.Columns-1) * ctSpatial.PixelSpacing(2));
            ct_y = double(ctSpatial.ImagePositionPatient(2)) + double((0:ctSpatial.Rows-1) * ctSpatial.PixelSpacing(1));
            ct_z = double(ctSpatial.SlicePositions);
            
            dose_x = double(doseSpatial.ImagePositionPatient(1)) + double((0:doseSpatial.Columns-1) * doseSpatial.PixelSpacing(2));
            dose_y = double(doseSpatial.ImagePositionPatient(2)) + double((0:doseSpatial.Rows-1) * doseSpatial.PixelSpacing(1));
            dose_z = double(doseSpatial.ImagePositionPatient(3)) + double(doseSpatial.GridFrameOffsetVector);
            
            [CT_Y, CT_X, CT_Z] = ndgrid(ct_y, ct_x, ct_z);
            [DOSE_Y, DOSE_X, DOSE_Z] = ndgrid(dose_y, dose_x, dose_z);
            
            try
                totalDoseResampled = interpn(double(CT_Y), double(CT_X), double(CT_Z), totalDose, ...
                    DOSE_Y, DOSE_X, DOSE_Z, 'linear', 0);
                fprintf('  - Total dose resampled: Max=%.4f Gy\n', max(totalDoseResampled(:)));
            catch ME
                fprintf('  - ERROR resampling total dose: %s\n', ME.message);
                totalDoseResampled = totalDose;
            end
            
            beamDosesResampled = cell(size(beamDoses));
            for beamIdx = 1:length(beamDoses)
                if ~isempty(beamDoses{beamIdx})
                    try
                        doseResampled = interpn(double(CT_Y), double(CT_X), double(CT_Z), ...
                            beamDoses{beamIdx}.physicalDose, DOSE_Y, DOSE_X, DOSE_Z, 'linear', 0);
                        beamDosesResampled{beamIdx} = beamDoses{beamIdx};
                        beamDosesResampled{beamIdx}.physicalDose = doseResampled;
                        beamDosesResampled{beamIdx}.maxDose = max(doseResampled(:));
                    catch
                        beamDosesResampled{beamIdx} = beamDoses{beamIdx};
                    end
                end
            end
            
            ctCube = ct_original.cubeHU{1};
            try
                ctResampled = interpn(double(CT_Y), double(CT_X), double(CT_Z), ctCube, ...
                    DOSE_Y, DOSE_X, DOSE_Z, 'linear', -1000);
            catch
                ctResampled = [];
            end
            
            ctResampled_struct = struct();
            if ~isempty(ctResampled)
                ctResampled_struct.cubeHU = {ctResampled};
                ctResampled_struct.cubeDim = size(ctResampled);
                ctResampled_struct.resolution.x = doseGrid.resolution(2);
                ctResampled_struct.resolution.y = doseGrid.resolution(1);
                ctResampled_struct.resolution.z = doseGrid.resolution(3);
            end
            
        else
            totalDoseResampled = totalDose;
            beamDosesResampled = beamDoses;
            ctResampled_struct = [];
        end
        
        %% Step 8: Save results and compare
        fprintf('\n[8/8] Saving results...\n');
        
        save(fullfile(outputPath, 'segmentData.mat'), 'segmentData');
        fprintf('  - segmentData.mat\n');
        
        save(fullfile(outputPath, 'segmentDoses.mat'), 'beamDosesResampled', ...
             'totalDoseResampled', 'referenceDose', 'doseGrid', 'ctSpatial', 'doseSpatial', ...
             'stf', 'pln', 'ctDownsampleFactor', '-v7.3');
        fprintf('  - segmentDoses.mat\n');
        
        if ~isempty(ctResampled_struct) && isfield(ctResampled_struct, 'cubeHU')
            save(fullfile(outputPath, 'ctResampled.mat'), 'ctResampled_struct', 'doseGrid');
            fprintf('  - ctResampled.mat\n');
        end
        
        for beamIdx = 1:length(beamDosesResampled)
            if ~isempty(beamDosesResampled{beamIdx})
                beamFilename = sprintf('Beam_%02d.mat', beamIdx);
                beamDataSave = beamDosesResampled{beamIdx};
                save(fullfile(outputPath, beamFilename), 'beamDataSave');
            end
        end
        fprintf('  - Beam_XX.mat files\n');
        
        % Compare with reference
        if ~isempty(referenceDose) && isequal(size(totalDoseResampled), size(referenceDose))
            fprintf('\n  Comparison with reference RTDOSE:\n');
            doseDiff = totalDoseResampled - referenceDose;
            fprintf('    Calculated max: %.4f Gy\n', max(totalDoseResampled(:)));
            fprintf('    Reference max:  %.4f Gy\n', max(referenceDose(:)));
            fprintf('    Mean abs diff:  %.4f Gy\n', mean(abs(doseDiff(:))));
            fprintf('    Max diff:       %.4f Gy\n', max(abs(doseDiff(:))));
            fprintf('    RMS diff:       %.4f Gy\n', sqrt(mean(doseDiff(:).^2)));
            
            comparison = struct();
            comparison.calculated = totalDoseResampled;
            comparison.reference = referenceDose;
            comparison.difference = doseDiff;
            comparison.metrics.meanAbsDiff = mean(abs(doseDiff(:)));
            comparison.metrics.maxDiff = max(abs(doseDiff(:)));
            comparison.metrics.rmsDiff = sqrt(mean(doseDiff(:).^2));
            
            save(fullfile(outputPath, 'doseComparison.mat'), 'comparison');
            fprintf('  - doseComparison.mat\n');
        end
        
        fprintf('\n========================================\n');
        fprintf('Processing complete: %s - %s\n', currentID, currentSession);
        fprintf('Output: %s\n', outputPath);
        fprintf('Cache: %s\n', patientCacheDir);
        fprintf('========================================\n\n');
        
    end
end

fprintf('All processing complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function [ct_ds, cst_ds] = downsampleCTComplete(ct, cst, factor)
    if factor == 1
        ct_ds = ct;
        cst_ds = cst;
        return;
    end
    
    ct_ds = ct;
    originalSize = ct.cubeDim;
    newSize = max(1, round(originalSize / factor));
    
    fprintf('    Downsampling CT: [%d,%d,%d] -> [%d,%d,%d]\n', ...
        originalSize(1), originalSize(2), originalSize(3), ...
        newSize(1), newSize(2), newSize(3));
    
    origRes = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
    
    ct_ds.resolution.x = ct.resolution.x * (originalSize(2) / newSize(2));
    ct_ds.resolution.y = ct.resolution.y * (originalSize(1) / newSize(1));
    ct_ds.resolution.z = ct.resolution.z * (originalSize(3) / newSize(3));
    
    newRes = [ct_ds.resolution.x, ct_ds.resolution.y, ct_ds.resolution.z];
    fprintf('    Resolution: [%.2f,%.2f,%.2f] -> [%.2f,%.2f,%.2f] mm\n', ...
        origRes(1), origRes(2), origRes(3), newRes(1), newRes(2), newRes(3));
    
    ct_ds.cubeDim = newSize;
    
    if isfield(ct, 'cubeHU') && ~isempty(ct.cubeHU)
        for i = 1:length(ct.cubeHU)
            if ~isempty(ct.cubeHU{i})
                ct_ds.cubeHU{i} = imresize3(ct.cubeHU{i}, newSize, 'linear');
            end
        end
    end
    
    if isfield(ct, 'cube') && ~isempty(ct.cube)
        for i = 1:length(ct.cube)
            if ~isempty(ct.cube{i})
                ct_ds.cube{i} = imresize3(ct.cube{i}, newSize, 'linear');
            end
        end
    end
    
    if isfield(ct, 'x') && ~isempty(ct.x)
        xMin = min(ct.x);
        xMax = max(ct.x);
        ct_ds.x = linspace(xMin, xMax, newSize(2));
    end
    
    if isfield(ct, 'y') && ~isempty(ct.y)
        yMin = min(ct.y);
        yMax = max(ct.y);
        ct_ds.y = linspace(yMin, yMax, newSize(1));
    end
    
    if isfield(ct, 'z') && ~isempty(ct.z)
        zMin = min(ct.z);
        zMax = max(ct.z);
        ct_ds.z = linspace(zMin, zMax, newSize(3));
    end
    
    if isfield(ct, 'numOfCtScen')
        ct_ds.numOfCtScen = ct.numOfCtScen;
    end
    
    cst_ds = cst;
    scaleFactor = newSize ./ originalSize;
    
    for i = 1:size(cst, 1)
        if size(cst, 2) >= 4 && ~isempty(cst{i, 4})
            for scen = 1:length(cst{i, 4})
                if ~isempty(cst{i, 4}{scen})
                    originalIndices = cst{i, 4}{scen};
                    [r, c, s] = ind2sub(originalSize, originalIndices);
                    r_new = max(1, min(newSize(1), round(r * scaleFactor(1))));
                    c_new = max(1, min(newSize(2), round(c * scaleFactor(2))));
                    s_new = max(1, min(newSize(3), round(s * scaleFactor(3))));
                    newIndices = sub2ind(newSize, r_new, c_new, s_new);
                    cst_ds{i, 4}{scen} = unique(newIndices);
                end
            end
        end
    end
    
    fprintf('    CT downsampling complete\n');
end