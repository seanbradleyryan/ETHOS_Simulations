%% ETHOS IMRT Segment-by-Segment Dose Calculator (FIXED VERSION)
% Purpose: Calculate individual segment doses from ETHOS exported IMRT data
% Uses MATRAD for dose calculation
%
% FIXES APPLIED:
%   1. Memory management - Don't store individual segment doses, only beam totals
%   2. MLC position inheritance between control points
%   3. Dual-layer MLC handling for ETHOS/Halcyon
%   4. Proper weight normalization verification
%
% KEY FEATURES:
%   - Calculates dose for each control point/segment within each beam
%   - IMRT plans have multiple segments per beam with different MLC configurations
%   - Each segment contributes a weighted portion of the total beam dose
%
% Author: Generated for ETHOS dose analysis
% Date: 2025 (Fixed version)

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
% 1 = original resolution (requires more memory but most accurate)
% 2 = half resolution in each dimension (8x faster, moderate accuracy)
ctDownsampleFactor = 1;  % Try running at full resolution with fixed memory management

% Caching Options
enableCaching = true;           % Set to false to force recalculation
cacheDir = fullfile(wd, 'Cache');  % Directory to store cached results

% Verbose output (set to false to reduce console spam)
verboseOutput = true;

% Memory management - clear large arrays after use
aggressiveMemoryCleanup = true;

%% ==================== INITIALIZATION ====================
fprintf('==========================================================\n');
fprintf('  ETHOS IMRT Segment-by-Segment Dose Calculator\n');
fprintf('  FIXED VERSION - Memory + MLC Inheritance + Dual-Layer\n');
fprintf('==========================================================\n\n');

fprintf('Optimization Settings:\n');
fprintf('  - CT Downsample Factor: %d\n', ctDownsampleFactor);
fprintf('  - Caching Enabled: %s\n', string(enableCaching));
fprintf('  - Aggressive Memory Cleanup: %s\n', string(aggressiveMemoryCleanup));
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

% Verify MATRAD functions
essentialFunctions = {'matRad_calcDoseInfluence', 'matRad_generateStf', 'matRad_calcDoseForward'};
for i = 1:length(essentialFunctions)
    if exist(essentialFunctions{i}, 'file')
        fprintf('  - Found: %s\n', essentialFunctions{i});
    else
        fprintf('  - WARNING: %s not found\n', essentialFunctions{i});
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
        
        % Create output directory
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        
        %% Step 1: Import DICOM data (with caching)
        fprintf('\n[1/8] Importing DICOM data...\n');
        
        dicomCacheDir = fullfile(cacheDir, currentID, currentSession);
        if ~exist(dicomCacheDir, 'dir')
            mkdir(dicomCacheDir);
        end
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
            
            % Check for cached downsampled CT
            dsCacheFile = fullfile(patientCacheDir, 'ct_downsampled.mat');
            
            if enableCaching && exist(dsCacheFile, 'file')
                fprintf('    Loading cached downsampled CT...\n');
                load(dsCacheFile, 'ct', 'cst');
                fprintf('    Loaded from cache: [%d,%d,%d]\n', ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));
            else
                % Perform downsampling using helper function (defined at end of script)
                [ct, cst] = downsampleCTComplete(ct_original, cst_original, ctDownsampleFactor);
                
                % Cache downsampled CT
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
                
                % Check if this is a PLAN dose (sum of all fractions)
                if isfield(rtdoseInfo, 'DoseSummationType') && strcmp(rtdoseInfo.DoseSummationType, 'PLAN')
                    % Get number of fractions from RTPLAN
                    numFractionsForDose = 10;  % Default
                    rtplanFile = dir(fullfile(dicomPath, 'RP*.dcm'));
                    if ~isempty(rtplanFile)
                        rtplanInfoTemp = dicominfo(fullfile(rtplanFile(1).folder, rtplanFile(1).name));
                        if isfield(rtplanInfoTemp, 'FractionGroupSequence')
                            fg = rtplanInfoTemp.FractionGroupSequence.Item_1;
                            if isfield(fg, 'NumberOfFractionsPlanned')
                                numFractionsForDose = fg.NumberOfFractionsPlanned;
                            end
                        end
                    end
                    fprintf('  WARNING: PLAN dose detected, dividing by %d fractions\n', numFractionsForDose);
                    referenceDose = referenceDose / numFractionsForDose;
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
        
        %% Step 3: Extract segment information from RTPLAN (with FIXED MLC inheritance)
        fprintf('\n[3/8] Extracting segment information from RTPLAN (with MLC inheritance)...\n');
        
        % Segment data is independent of CT resolution, so cache at patient level
        segmentCacheFile = fullfile(dicomCacheDir, 'segmentData_v2_multilayer.mat');
        
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
                    beamData.hasDualLayerMLC = false;
                    
                    if isfield(beam, 'BeamName')
                        beamData.beamName = beam.BeamName;
                    else
                        beamData.beamName = sprintf('Beam_%d', beamIdx);
                    end
                    
                    if beamIdx <= length(beamMetersets)
                        beamData.beamMeterset = beamMetersets(beamIdx);
                    end
                    
                    % Get leaf boundaries and detect dual-layer MLC
                    % Store ALL MLC device info for proper dual-layer handling
                    beamData.mlcDevices = {};  % Cell array of MLC device info
                    
                    if isfield(beam, 'BeamLimitingDeviceSequence')
                        numDevices = length(fieldnames(beam.BeamLimitingDeviceSequence));
                        for devIdx = 1:numDevices
                            device = beam.BeamLimitingDeviceSequence.(sprintf('Item_%d', devIdx));
                            if isfield(device, 'RTBeamLimitingDeviceType')
                                devType = device.RTBeamLimitingDeviceType;
                                if contains(devType, 'MLC', 'IgnoreCase', true)
                                    mlcDeviceInfo = struct();
                                    mlcDeviceInfo.type = devType;
                                    mlcDeviceInfo.leafBoundaries = [];
                                    mlcDeviceInfo.numLeafPairs = 0;
                                    
                                    if isfield(device, 'LeafPositionBoundaries')
                                        mlcDeviceInfo.leafBoundaries = device.LeafPositionBoundaries;
                                        mlcDeviceInfo.numLeafPairs = length(device.LeafPositionBoundaries) - 1;
                                    elseif isfield(device, 'NumberOfLeafJawPairs')
                                        mlcDeviceInfo.numLeafPairs = device.NumberOfLeafJawPairs;
                                    end
                                    
                                    beamData.mlcDevices{end+1} = mlcDeviceInfo;
                                    
                                    % Keep first device as default for backward compatibility
                                    if isempty(beamData.leafBoundaries) && ~isempty(mlcDeviceInfo.leafBoundaries)
                                        beamData.leafBoundaries = mlcDeviceInfo.leafBoundaries;
                                        beamData.numLeafPairs = mlcDeviceInfo.numLeafPairs;
                                    end
                                end
                            end
                        end
                    end
                    
                    % Check for dual-layer MLC (ETHOS/Halcyon)
                    if length(beamData.mlcDevices) > 1
                        beamData.hasDualLayerMLC = true;
                        mlcTypes = cellfun(@(x) x.type, beamData.mlcDevices, 'UniformOutput', false);
                        fprintf('    Beam %d: Dual-layer MLC detected (%s)\n', beamIdx, strjoin(mlcTypes, ', '));
                        for mlcIdx = 1:length(beamData.mlcDevices)
                            fprintf('      Layer %d (%s): %d leaf pairs\n', mlcIdx, ...
                                beamData.mlcDevices{mlcIdx}.type, beamData.mlcDevices{mlcIdx}.numLeafPairs);
                        end
                    end
                    
                    % Process control points WITH INHERITANCE
                    if isfield(beam, 'ControlPointSequence')
                        numCP = length(fieldnames(beam.ControlPointSequence));
                        controlPoints = cell(numCP, 1);
                        
                        % Initialize inherited values
                        inheritedMlcPositions = {};  % Cell array for dual-layer
                        inheritedJawX = [];
                        inheritedJawY = [];
                        inheritedGantryAngle = 0;
                        
                        for cpIdx = 1:numCP
                            cp = beam.ControlPointSequence.(sprintf('Item_%d', cpIdx));
                            
                            cpData = struct();
                            cpData.cpIdx = cpIdx;
                            cpData.cumulativeMetersetWeight = 0;
                            cpData.gantryAngle = inheritedGantryAngle;
                            cpData.mlcPositions = inheritedMlcPositions;  % INHERIT from previous CP
                            cpData.jawX = inheritedJawX;  % INHERIT from previous CP
                            cpData.jawY = inheritedJawY;  % INHERIT from previous CP
                            
                            if isfield(cp, 'CumulativeMetersetWeight')
                                cpData.cumulativeMetersetWeight = cp.CumulativeMetersetWeight;
                            end
                            
                            if isfield(cp, 'GantryAngle')
                                cpData.gantryAngle = cp.GantryAngle;
                                inheritedGantryAngle = cp.GantryAngle;
                                beamData.gantryAngle = cp.GantryAngle;
                            end
                            
                            if isfield(cp, 'PatientSupportAngle')
                                beamData.couchAngle = cp.PatientSupportAngle;
                            end
                            
                            if isfield(cp, 'BeamLimitingDevicePositionSequence')
                                numDevices = length(fieldnames(cp.BeamLimitingDevicePositionSequence));
                                
                                % Collect ALL MLC layers for this CP
                                cpMlcLayers = {};
                                
                                for devIdx = 1:numDevices
                                    device = cp.BeamLimitingDevicePositionSequence.(sprintf('Item_%d', devIdx));
                                    if isfield(device, 'RTBeamLimitingDeviceType') && isfield(device, 'LeafJawPositions')
                                        deviceType = device.RTBeamLimitingDeviceType;
                                        positions = device.LeafJawPositions;
                                        
                                        if contains(deviceType, 'ASYMX', 'IgnoreCase', true) || strcmp(deviceType, 'X')
                                            cpData.jawX = positions;
                                            inheritedJawX = positions;
                                        elseif contains(deviceType, 'ASYMY', 'IgnoreCase', true) || strcmp(deviceType, 'Y')
                                            cpData.jawY = positions;
                                            inheritedJawY = positions;
                                        elseif contains(deviceType, 'MLC', 'IgnoreCase', true)
                                            % Store each MLC layer separately
                                            mlcLayer = struct();
                                            mlcLayer.type = deviceType;
                                            mlcLayer.positions = positions;
                                            cpMlcLayers{end+1} = mlcLayer;
                                        end
                                    end
                                end
                                
                                % If we found MLC data, update inherited values
                                if ~isempty(cpMlcLayers)
                                    cpData.mlcPositions = cpMlcLayers;
                                    inheritedMlcPositions = cpMlcLayers;
                                end
                            end
                            
                            controlPoints{cpIdx} = cpData;
                        end
                        
                        % Create segments
                        numSegments = numCP - 1;
                        for segIdx = 1:numSegments
                            cpStart = controlPoints{segIdx};
                            cpEnd = controlPoints{segIdx + 1};
                            
                            segment = struct();
                            segment.segmentIdx = segIdx;
                            segment.beamIdx = beamIdx;
                            segment.segmentWeight = cpEnd.cumulativeMetersetWeight - cpStart.cumulativeMetersetWeight;
                            segment.segmentMeterset = segment.segmentWeight * beamData.beamMeterset;
                            segment.gantryAngle = cpStart.gantryAngle;
                            segment.jawX = cpStart.jawX;
                            segment.jawY = cpStart.jawY;
                            
                            % Store ALL MLC layers - do NOT reduce here
                            % Reduction will happen during weight calculation where we check all layers
                            segment.mlcLayers = cpStart.mlcPositions;  % Cell array of layer structs
                            
                            beamData.segments{segIdx} = segment;
                            totalSegments = totalSegments + 1;
                        end
                    end
                    
                    segmentData.beams{beamIdx} = beamData;
                    fprintf('    Beam %d: %s, Gantry=%.1f°, %.1f MU, %d segments\n', ...
                        beamIdx, beamData.beamName, beamData.gantryAngle, beamData.beamMeterset, length(beamData.segments));
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
                fprintf('  Stack trace:\n');
                for k = 1:length(ME.stack)
                    fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                end
                continue;
            end
        end
        
        % Verify segment weights sum correctly
        fprintf('\n  Verifying segment weights...\n');
        for beamIdx = 1:segmentData.numBeams
            beamData = segmentData.beams{beamIdx};
            totalWeight = 0;
            for segIdx = 1:length(beamData.segments)
                totalWeight = totalWeight + beamData.segments{segIdx}.segmentWeight;
            end
            fprintf('    Beam %d: Total weight = %.6f (should be ~1.0)\n', beamIdx, totalWeight);
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
        
        % Set bixel width
        pln.propStf.bixelWidth = 5;
        
        % IMPORTANT: Set dose grid resolution to match the (possibly downsampled) CT
        pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x;
        pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y;
        pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z;
        
        % Set algorithm
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
        
        %% Step 5: Generate steering information (with caching)
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
                fprintf('  Error details: %s\n', ME.message);
                for k = 1:length(ME.stack)
                    fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                end
                continue;
            end
        end
        
        for i = 1:length(stf)
            fprintf('    Beam %d: Gantry=%.1f, Rays=%d\n', i, stf(i).gantryAngle, stf(i).numOfRays);
        end
        
        %% Step 6: Calculate segment-by-segment doses (FIXED MEMORY MANAGEMENT)
        fprintf('\n[6/8] Calculating segment-by-segment doses...\n');
        
        % Initialize storage - ONLY store beam totals, not individual segments
        calculatedGridSize = ct.cubeDim;
        totalDose = zeros(calculatedGridSize);
        beamDoses = cell(length(stf), 1);
        
        % Statistics tracking
        totalCalcTime = 0;
        cachedBeams = 0;
        calculatedBeams = 0;
        totalSegmentsProcessed = 0;
        
        for beamIdx = 1:length(stf)
            beamData = segmentData.beams{beamIdx};
            numSegments = length(beamData.segments);
            
            fprintf('\n  Beam %d/%d (%s): %d segments, %.1f MU\n', ...
                beamIdx, length(stf), beamData.beamName, numSegments, beamData.beamMeterset);
            
            % Check for cached beam dose
            beamDoseCacheFile = fullfile(patientCacheDir, sprintf('beam_%02d_dose_v2.mat', beamIdx));
            
            if enableCaching && exist(beamDoseCacheFile, 'file')
                fprintf('    Loading cached beam dose...\n');
                load(beamDoseCacheFile, 'beamDoseResult');
                
                beamDoses{beamIdx} = beamDoseResult;
                totalDose = totalDose + beamDoseResult.physicalDose;
                cachedBeams = cachedBeams + 1;
                totalSegmentsProcessed = totalSegmentsProcessed + numSegments;
                
                fprintf('    Loaded from cache (Max: %.4f Gy)\n', beamDoseResult.maxDose);
                continue;
            end
            
            % Calculate dij for this beam (with caching)
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
            segmentsWithDose = 0;
            segmentsSkipped = 0;
            
            % Track aperture statistics for debugging
            totalOpenBixels = 0;
            
            % Process each segment - ACCUMULATE ONLY, don't store individual segment doses
            tic;
            firstSegmentDiagnostics = true;  % Show diagnostics for first valid segment
            
            for segIdx = 1:numSegments
                segment = beamData.segments{segIdx};
                
                if segment.segmentWeight <= 0
                    segmentsSkipped = segmentsSkipped + 1;
                    continue;
                end
                
                % Create weight vector based on MLC aperture
                w = zeros(dij.totalNumOfBixels, 1);
                
                mlcLayers = segment.mlcLayers;  % Cell array of MLC layer data
                jawX = segment.jawX;
                jawY = segment.jawY;
                
                openBixelsThisSegment = 0;
                
                % Diagnostic output for first segment
                if firstSegmentDiagnostics && verboseOutput
                    fprintf('    First segment MLC diagnostics:\n');
                    if iscell(mlcLayers) && ~isempty(mlcLayers)
                        fprintf('      MLC layers: %d\n', length(mlcLayers));
                        for lIdx = 1:length(mlcLayers)
                            if ~isempty(mlcLayers{lIdx}) && isfield(mlcLayers{lIdx}, 'positions')
                                fprintf('        Layer %d: %d positions\n', lIdx, length(mlcLayers{lIdx}.positions));
                            end
                        end
                    else
                        fprintf('      WARNING: No MLC layer data!\n');
                    end
                    if isfield(beamData, 'mlcDevices')
                        fprintf('      MLC devices with boundaries: %d\n', length(beamData.mlcDevices));
                        for dIdx = 1:length(beamData.mlcDevices)
                            dev = beamData.mlcDevices{dIdx};
                            fprintf('        Device %d (%s): %d leaf pairs\n', dIdx, dev.type, dev.numLeafPairs);
                        end
                    end
                    firstSegmentDiagnostics = false;
                end
                
                % Check if we have MLC data and device info
                hasMlcData = ~isempty(mlcLayers) && iscell(mlcLayers) && ~isempty(mlcLayers);
                hasMlcDevices = isfield(beamData, 'mlcDevices') && ~isempty(beamData.mlcDevices);
                
                if hasMlcData && hasMlcDevices
                    numRays = stfSingle.numOfRays;
                    bixelIdx = 0;
                    
                    for rayIdx = 1:numRays
                        ray = stfSingle.ray(rayIdx);
                        rayPosX = ray.rayPos_bev(1);
                        rayPosY = ray.rayPos_bev(2);
                        
                        % Check jaw limits first
                        inJawX = isempty(jawX) || (length(jawX) >= 2 && rayPosX >= jawX(1) && rayPosX <= jawX(2));
                        inJawY = isempty(jawY) || (length(jawY) >= 2 && rayPosY >= jawY(1) && rayPosY <= jawY(2));
                        
                        % Check ALL MLC layers - ray must pass through ALL
                        inAllMlcLayers = inJawX && inJawY;  % Start assuming open if within jaws
                        
                        if inAllMlcLayers
                            % Check each MLC layer
                            for layerIdx = 1:length(mlcLayers)
                                mlcLayer = mlcLayers{layerIdx};
                                
                                if isempty(mlcLayer) || ~isfield(mlcLayer, 'positions')
                                    continue;  % Skip if no data for this layer
                                end
                                
                                mlcPos = mlcLayer.positions;
                                
                                % Get the corresponding device info for leaf boundaries
                                if layerIdx <= length(beamData.mlcDevices)
                                    deviceInfo = beamData.mlcDevices{layerIdx};
                                    leafBoundaries = deviceInfo.leafBoundaries;
                                    numLeafPairs = deviceInfo.numLeafPairs;
                                else
                                    % Fallback: try to infer from positions
                                    numLeafPairs = length(mlcPos) / 2;
                                    leafBoundaries = [];
                                end
                                
                                if isempty(leafBoundaries)
                                    % Can't check this layer without boundaries
                                    continue;
                                end
                                
                                % Extract left and right banks for this layer
                                if length(mlcPos) >= 2 * numLeafPairs
                                    leftBank = mlcPos(1:numLeafPairs);
                                    rightBank = mlcPos(numLeafPairs+1:2*numLeafPairs);
                                else
                                    % Incomplete MLC data for this layer
                                    continue;
                                end
                                
                                % Check if ray is within aperture for THIS layer
                                inThisLayer = false;
                                for leafIdx = 1:numLeafPairs
                                    if leafIdx < length(leafBoundaries)
                                        leafYMin = leafBoundaries(leafIdx);
                                        leafYMax = leafBoundaries(leafIdx + 1);
                                        
                                        if rayPosY >= leafYMin && rayPosY < leafYMax
                                            % This leaf covers the ray's Y position
                                            if rayPosX >= leftBank(leafIdx) && rayPosX <= rightBank(leafIdx)
                                                inThisLayer = true;
                                            end
                                            break;  % Found the covering leaf, no need to check others
                                        end
                                    end
                                end
                                
                                % If blocked by ANY layer, ray is blocked
                                if ~inThisLayer
                                    inAllMlcLayers = false;
                                    break;
                                end
                            end
                        end
                        
                        numEnergies = length(ray.energy);
                        for energyIdx = 1:numEnergies
                            bixelIdx = bixelIdx + 1;
                            if inAllMlcLayers
                                w(bixelIdx) = segment.segmentMeterset;  % Use segment MU directly
                                openBixelsThisSegment = openBixelsThisSegment + 1;
                            end
                        end
                    end
                elseif hasMlcData && ~isempty(beamData.leafBoundaries)
                    % Fallback: single layer mode using default boundaries
                    % This handles cases where mlcDevices wasn't populated
                    numLeafPairs = beamData.numLeafPairs;
                    leafBoundaries = beamData.leafBoundaries;
                    
                    % Get positions from first layer
                    mlcPos = [];
                    if iscell(mlcLayers) && ~isempty(mlcLayers{1})
                        if isfield(mlcLayers{1}, 'positions')
                            mlcPos = mlcLayers{1}.positions;
                        end
                    end
                    
                    if ~isempty(mlcPos) && length(mlcPos) >= 2 * numLeafPairs
                        leftBank = mlcPos(1:numLeafPairs);
                        rightBank = mlcPos(numLeafPairs+1:2*numLeafPairs);
                        
                        numRays = stfSingle.numOfRays;
                        bixelIdx = 0;
                        
                        for rayIdx = 1:numRays
                            ray = stfSingle.ray(rayIdx);
                            rayPosX = ray.rayPos_bev(1);
                            rayPosY = ray.rayPos_bev(2);
                            
                            inJawX = isempty(jawX) || (length(jawX) >= 2 && rayPosX >= jawX(1) && rayPosX <= jawX(2));
                            inJawY = isempty(jawY) || (length(jawY) >= 2 && rayPosY >= jawY(1) && rayPosY <= jawY(2));
                            
                            inMLC = false;
                            if inJawX && inJawY
                                for leafIdx = 1:numLeafPairs
                                    if rayPosY >= leafBoundaries(leafIdx) && rayPosY < leafBoundaries(leafIdx + 1)
                                        if rayPosX >= leftBank(leafIdx) && rayPosX <= rightBank(leafIdx)
                                            inMLC = true;
                                        end
                                        break;
                                    end
                                end
                            end
                            
                            numEnergies = length(ray.energy);
                            for energyIdx = 1:numEnergies
                                bixelIdx = bixelIdx + 1;
                                if inMLC
                                    w(bixelIdx) = segment.segmentMeterset;
                                    openBixelsThisSegment = openBixelsThisSegment + 1;
                                end
                            end
                        end
                    else
                        fprintf('      Seg %d: WARNING - Incomplete MLC data in fallback mode\n', segIdx);
                    end
                else
                    % No MLC data - this shouldn't happen with inheritance fix
                    fprintf('      Seg %d: WARNING - No MLC data, using uniform weights\n', segIdx);
                    w(:) = segment.segmentMeterset / dij.totalNumOfBixels;
                    openBixelsThisSegment = dij.totalNumOfBixels;
                end
                
                totalOpenBixels = totalOpenBixels + openBixelsThisSegment;
                
                % Calculate dose using direct matrix multiplication
                if any(w > 0)
                    segmentDose = reshape(full(dij.physicalDose{1} * w), calculatedGridSize);
                    beamDoseAccum = beamDoseAccum + segmentDose;
                    segmentsWithDose = segmentsWithDose + 1;
                    
                    if verboseOutput && mod(segIdx, 20) == 0
                        fprintf('      Seg %d: Open=%d bixels, Max=%.4f Gy\n', ...
                            segIdx, openBixelsThisSegment, max(segmentDose(:)));
                    end
                    
                    % Clear segment dose immediately to save memory
                    clear segmentDose;
                end
            end
            segmentCalcTime = toc;
            
            % Calculate average open bixels per segment
            avgOpenBixels = totalOpenBixels / max(1, segmentsWithDose);
            fprintf('    Processed %d segments (skipped %d) in %.1f sec\n', ...
                segmentsWithDose, segmentsSkipped, segmentCalcTime);
            fprintf('    Average open bixels per segment: %.1f / %d (%.1f%%)\n', ...
                avgOpenBixels, dij.totalNumOfBixels, 100*avgOpenBixels/dij.totalNumOfBixels);
            
            totalCalcTime = totalCalcTime + segmentCalcTime;
            totalSegmentsProcessed = totalSegmentsProcessed + segmentsWithDose;
            
            % Store beam dose (only the accumulated total)
            beamDoseResult = struct();
            beamDoseResult.beamIdx = beamIdx;
            beamDoseResult.beamName = beamData.beamName;
            beamDoseResult.gantryAngle = beamData.gantryAngle;
            beamDoseResult.couchAngle = beamData.couchAngle;
            beamDoseResult.beamMeterset = beamData.beamMeterset;
            beamDoseResult.numSegments = numSegments;
            beamDoseResult.segmentsWithDose = segmentsWithDose;
            beamDoseResult.avgOpenBixels = avgOpenBixels;
            beamDoseResult.physicalDose = beamDoseAccum;
            beamDoseResult.maxDose = max(beamDoseAccum(:));
            
            beamDoses{beamIdx} = beamDoseResult;
            totalDose = totalDose + beamDoseAccum;
            calculatedBeams = calculatedBeams + 1;
            
            % Cache beam result
            if enableCaching
                save(beamDoseCacheFile, 'beamDoseResult', '-v7.3');
            end
            
            fprintf('    Beam %d complete: Max=%.4f Gy\n', beamIdx, beamDoseResult.maxDose);
            
            % Clear dij to save memory before next beam
            if aggressiveMemoryCleanup
                clear dij beamDoseAccum;
            end
        end
        
        fprintf('\n  Summary:\n');
        fprintf('    Beams calculated: %d, from cache: %d\n', calculatedBeams, cachedBeams);
        fprintf('    Total segments processed: %d\n', totalSegmentsProcessed);
        fprintf('    Total calculation time: %.1f sec\n', totalCalcTime);
        fprintf('    Total dose max: %.4f Gy\n', max(totalDose(:)));
        
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
            
            % CT coordinate vectors (original resolution)
            ct_x = double(ctSpatial.ImagePositionPatient(1)) + double((0:ctSpatial.Columns-1) * ctSpatial.PixelSpacing(2));
            ct_y = double(ctSpatial.ImagePositionPatient(2)) + double((0:ctSpatial.Rows-1) * ctSpatial.PixelSpacing(1));
            ct_z = double(ctSpatial.SlicePositions);
            
            % RTDOSE coordinate vectors
            dose_x = double(doseSpatial.ImagePositionPatient(1)) + double((0:doseSpatial.Columns-1) * doseSpatial.PixelSpacing(2));
            dose_y = double(doseSpatial.ImagePositionPatient(2)) + double((0:doseSpatial.Rows-1) * doseSpatial.PixelSpacing(1));
            dose_z = double(doseSpatial.ImagePositionPatient(3)) + double(doseSpatial.GridFrameOffsetVector);
            
            [CT_Y, CT_X, CT_Z] = ndgrid(ct_y, ct_x, ct_z);
            [DOSE_Y, DOSE_X, DOSE_Z] = ndgrid(dose_y, dose_x, dose_z);
            
            % Resample total dose
            try
                totalDoseResampled = interpn(double(CT_Y), double(CT_X), double(CT_Z), totalDose, ...
                    DOSE_Y, DOSE_X, DOSE_Z, 'linear', 0);
                fprintf('  - Total dose resampled: Max=%.4f Gy\n', max(totalDoseResampled(:)));
            catch ME
                fprintf('  - ERROR resampling total dose: %s\n', ME.message);
                totalDoseResampled = totalDose;
            end
            
            % Resample beam doses
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
            
            % Resample CT
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
        
        % Save segment data
        save(fullfile(outputPath, 'segmentData.mat'), 'segmentData');
        fprintf('  - segmentData.mat\n');
        
        % Save main results (WITHOUT individual segment doses to save space)
        save(fullfile(outputPath, 'segmentDoses.mat'), 'beamDosesResampled', ...
             'totalDoseResampled', 'referenceDose', 'doseGrid', 'ctSpatial', 'doseSpatial', ...
             'stf', 'pln', 'ctDownsampleFactor', '-v7.3');
        fprintf('  - segmentDoses.mat\n');
        
        % Save CT resampled structure
        if ~isempty(ctResampled_struct) && isfield(ctResampled_struct, 'cubeHU')
            save(fullfile(outputPath, 'ctResampled.mat'), 'ctResampled_struct', 'doseGrid');
            fprintf('  - ctResampled.mat\n');
        end
        
        % Save individual beam doses
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
            
            % Calculate ratio for debugging
            calcMax = max(totalDoseResampled(:));
            refMax = max(referenceDose(:));
            fprintf('    Ratio (calc/ref): %.2f\n', calcMax/refMax);
            
            comparison = struct();
            comparison.calculated = totalDoseResampled;
            comparison.reference = referenceDose;
            comparison.difference = doseDiff;
            comparison.metrics.meanAbsDiff = mean(abs(doseDiff(:)));
            comparison.metrics.maxDiff = max(abs(doseDiff(:)));
            comparison.metrics.rmsDiff = sqrt(mean(doseDiff(:).^2));
            comparison.metrics.ratio = calcMax/refMax;
            
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
    % DOWNSAMPLECTCOMPLETE - Comprehensive CT downsampling for MATRAD compatibility
    %
    % This function downsamples a CT structure and adjusts all necessary fields
    % that MATRAD uses for dose calculation, including coordinate vectors.
    
    if factor == 1
        ct_ds = ct;
        cst_ds = cst;
        return;
    end
    
    ct_ds = ct;
    originalSize = ct.cubeDim;
    
    % Calculate new size (ensure at least 1 in each dimension)
    newSize = max(1, round(originalSize / factor));
    
    fprintf('    Downsampling CT: [%d,%d,%d] -> [%d,%d,%d]\n', ...
        originalSize(1), originalSize(2), originalSize(3), ...
        newSize(1), newSize(2), newSize(3));
    
    % Store original resolution for reference
    origRes = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
    
    % Update resolution based on actual size change
    ct_ds.resolution.x = ct.resolution.x * (originalSize(2) / newSize(2));
    ct_ds.resolution.y = ct.resolution.y * (originalSize(1) / newSize(1));
    ct_ds.resolution.z = ct.resolution.z * (originalSize(3) / newSize(3));
    
    newRes = [ct_ds.resolution.x, ct_ds.resolution.y, ct_ds.resolution.z];
    fprintf('    Resolution: [%.2f,%.2f,%.2f] -> [%.2f,%.2f,%.2f] mm\n', ...
        origRes(1), origRes(2), origRes(3), newRes(1), newRes(2), newRes(3));
    
    % Update dimensions
    ct_ds.cubeDim = newSize;
    
    % Downsample the HU cube(s)
    if isfield(ct, 'cubeHU') && ~isempty(ct.cubeHU)
        for i = 1:length(ct.cubeHU)
            if ~isempty(ct.cubeHU{i})
                ct_ds.cubeHU{i} = imresize3(ct.cubeHU{i}, newSize, 'linear');
            end
        end
    end
    
    % Downsample density cube if it exists
    if isfield(ct, 'cube') && ~isempty(ct.cube)
        for i = 1:length(ct.cube)
            if ~isempty(ct.cube{i})
                ct_ds.cube{i} = imresize3(ct.cube{i}, newSize, 'linear');
            end
        end
    end
    
    % Update coordinate vectors
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
    
    % Update numOfCtScen if it exists
    if isfield(ct, 'numOfCtScen')
        ct_ds.numOfCtScen = ct.numOfCtScen;
    end
    
    % Adjust CST structure indices for the new grid size
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