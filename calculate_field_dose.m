%% ETHOS IMRT Segment-by-Segment Dose Calculator
% Purpose: Calculate individual segment doses from ETHOS exported IMRT data
% Uses MATRAD for dose calculation
%
% KEY CHANGE FROM PREVIOUS VERSION:
%   - Calculates dose for each control point/segment within each beam
%   - IMRT plans have multiple segments per beam with different MLC configurations
%   - Each segment contributes a weighted portion of the total beam dose
%
% Workflow:
%   1. Import DICOM data using MATRAD
%   2. Load reference RTDOSE and extract spatial coordinates
%   3. Extract segment information from RTPLAN (MLC positions, weights)
%   4. Configure dose calculation
%   5. Generate steering information
%   6. For each beam and each segment:
%      - Calculate dose contribution for that segment
%      - Weight by segment's meterset fraction
%   7. Resample to RTDOSE grid
%   8. Save results and compare with reference
%
% Author: Generated for ETHOS dose analysis
% Date: 2025

clear; clc; close all;

%% Configuration
% Patient and session arrays (expandable for batch processing)
ids = {'1194203'};
sessions = {'Session_1'};

% Base directory
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
matradPath = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

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

% Add MATRAD to path
addpath(genpath(matradPath));
fprintf('  - matRad added to MATLAB path\n');

% Initialize MATRAD (required for some versions)
fprintf('Initializing MATRAD...\n');
try
    matRad_rc; % matRad run configuration/initialization
    fprintf('  - MATRAD initialized successfully\n');
catch
    fprintf('  - matRad_rc not found, trying matRad_cfg...\n');
    try
        matRad_cfg = MatRad_Config.instance();
        fprintf('  - MATRAD configured successfully\n');
    catch
        fprintf('  - Direct initialization not available, continuing...\n');
    end
end

% Verify MATRAD functions are available
fprintf('Verifying MATRAD installation...\n');
essentialFunctions = {'matRad_calcDoseInfluence', 'matRad_generateStf', 'matRad_calcDoseForward'};
for i = 1:length(essentialFunctions)
    if exist(essentialFunctions{i}, 'file')
        fprintf('  - Found: %s\n', essentialFunctions{i});
    else
        fprintf('  - WARNING: %s not found\n', essentialFunctions{i});
    end
end

%% Main Processing Loop
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
        
        % Create output directory
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        
        %% Step 1: Import DICOM data using MATRAD
        fprintf('\n[1/8] Importing DICOM data...\n');
        
        try
            % Use matRad_DicomImporter class (object-oriented approach)
            fprintf('  Creating matRad_DicomImporter object...\n');
            
            % Create importer instance
            importer = matRad_DicomImporter(dicomPath);
            
            % Import all DICOM data (CT, structures, plan)
            fprintf('  Importing DICOM files...\n');
            importer.matRad_importDicom();
            
            fprintf('  - Data imported successfully\n');
            
            % Display CT information
            fprintf('  - CT dimensions: %d x %d x %d\n', ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));
            fprintf('  - CT resolution: %.2f x %.2f x %.2f mm\n', ...
                ct.resolution.x, ct.resolution.y, ct.resolution.z);
            fprintf('  - Number of structures: %d\n', size(cst, 1));
            if isfield(pln, 'propStf')
                fprintf('  - Number of beams: %d\n', pln.propStf.numOfBeams);
            elseif isfield(pln, 'numOfBeams')
                fprintf('  - Number of beams: %d\n', pln.numOfBeams);
            end
            
        catch ME
            fprintf('Error importing DICOM data: %s\n', ME.message);
            fprintf('Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            continue;
        end
        
        %% Step 2: Load reference RTDOSE and extract spatial coordinates
        fprintf('\n[2/8] Loading reference RTDOSE and extracting spatial coordinates...\n');
        
        % Initialize spatial coordinate structures
        ctSpatial = struct();
        doseSpatial = struct();
        
        try
            % Extract CT spatial coordinates from DICOM files
            fprintf('  Extracting CT spatial coordinates...\n');
            ctFiles = dir(fullfile(dicomPath, 'CT*.dcm'));
            if isempty(ctFiles)
                ctFiles = dir(fullfile(dicomPath, '*CT*.dcm'));
            end
            
            if ~isempty(ctFiles)
                % Read first CT slice for ImagePositionPatient and PixelSpacing
                ctInfo = dicominfo(fullfile(ctFiles(1).folder, ctFiles(1).name));
                
                ctSpatial.ImagePositionPatient = ctInfo.ImagePositionPatient;
                ctSpatial.PixelSpacing = ctInfo.PixelSpacing;
                ctSpatial.Rows = ctInfo.Rows;
                ctSpatial.Columns = ctInfo.Columns;
                
                % Get all slice positions to determine Z coordinates
                slicePositions = zeros(length(ctFiles), 1);
                for i = 1:length(ctFiles)
                    tempInfo = dicominfo(fullfile(ctFiles(i).folder, ctFiles(i).name));
                    slicePositions(i) = tempInfo.ImagePositionPatient(3);
                end
                slicePositions = sort(slicePositions);
                
                ctSpatial.SlicePositions = slicePositions;
                ctSpatial.SliceThickness = abs(slicePositions(2) - slicePositions(1));
                ctSpatial.NumSlices = length(slicePositions);
                
                % Calculate CT bounding box in patient coordinates
                ctSpatial.XMin = ctSpatial.ImagePositionPatient(1);
                ctSpatial.XMax = ctSpatial.XMin + (ctSpatial.Columns - 1) * ctSpatial.PixelSpacing(2);
                ctSpatial.YMin = ctSpatial.ImagePositionPatient(2);
                ctSpatial.YMax = ctSpatial.YMin + (ctSpatial.Rows - 1) * ctSpatial.PixelSpacing(1);
                ctSpatial.ZMin = min(slicePositions);
                ctSpatial.ZMax = max(slicePositions);
                
                fprintf('    CT spatial extent:\n');
                fprintf('      X: [%.2f, %.2f] mm\n', ctSpatial.XMin, ctSpatial.XMax);
                fprintf('      Y: [%.2f, %.2f] mm\n', ctSpatial.YMin, ctSpatial.YMax);
                fprintf('      Z: [%.2f, %.2f] mm\n', ctSpatial.ZMin, ctSpatial.ZMax);
                fprintf('      Resolution: [%.2f, %.2f, %.2f] mm\n', ...
                    ctSpatial.PixelSpacing(1), ctSpatial.PixelSpacing(2), ctSpatial.SliceThickness);
            else
                fprintf('  WARNING: No CT DICOM files found for spatial extraction\n');
            end
            
            % Load RTDOSE and extract spatial coordinates
            rtdoseFile = dir(fullfile(dicomPath, 'RD*.dcm'));
            if isempty(rtdoseFile)
                rtdoseFile = dir(fullfile(dicomPath, '*RTDOSE*.dcm'));
            end
            
            if ~isempty(rtdoseFile)
                rtdoseInfo = dicominfo(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name));
                referenceDose = dicomread(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name));
                
                % RTDOSE is often 4D with singleton dimension - squeeze it
                referenceDose = squeeze(referenceDose);
                
                % Apply scaling
                referenceDose = double(referenceDose) * rtdoseInfo.DoseGridScaling;
                
                fprintf('  - Reference dose grid: %d x %d x %d\n', ...
                    size(referenceDose, 1), size(referenceDose, 2), size(referenceDose, 3));
                fprintf('  - Max dose: %.2f Gy\n', max(referenceDose(:)));
                
                % Check if this is a PLAN dose (cumulative) vs FRACTION dose
                if isfield(rtdoseInfo, 'DoseSummationType')
                    fprintf('  - DoseSummationType: %s\n', rtdoseInfo.DoseSummationType);
                    if strcmp(rtdoseInfo.DoseSummationType, 'PLAN')
                        fprintf('  WARNING: Entire planned dose detected. Dividing by 10 fractions\n');
                        referenceDose = referenceDose ./ 10;
                        fprintf('  - Adjusted max dose: %.2f Gy\n', max(referenceDose(:)));
                    end
                end
                
                % Extract RTDOSE spatial coordinates
                doseSpatial.ImagePositionPatient = rtdoseInfo.ImagePositionPatient;
                doseSpatial.PixelSpacing = rtdoseInfo.PixelSpacing;
                doseSpatial.Rows = rtdoseInfo.Rows;
                doseSpatial.Columns = rtdoseInfo.Columns;
                
                % Z resolution from GridFrameOffsetVector (multiframe DICOM)
                if isfield(rtdoseInfo, 'GridFrameOffsetVector') && length(rtdoseInfo.GridFrameOffsetVector) > 1
                    doseSpatial.GridFrameOffsetVector = rtdoseInfo.GridFrameOffsetVector;
                    doseSpatial.SliceThickness = abs(rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1));
                    doseSpatial.NumSlices = length(rtdoseInfo.GridFrameOffsetVector);
                    fprintf('  - Z resolution from GridFrameOffsetVector: %.3f mm\n', doseSpatial.SliceThickness);
                elseif isfield(rtdoseInfo, 'SliceThickness')
                    doseSpatial.SliceThickness = rtdoseInfo.SliceThickness;
                    doseSpatial.NumSlices = size(referenceDose, 3);
                    doseSpatial.GridFrameOffsetVector = (0:doseSpatial.NumSlices-1) * doseSpatial.SliceThickness;
                    fprintf('  - Z resolution from SliceThickness: %.3f mm\n', doseSpatial.SliceThickness);
                else
                    % Fallback
                    doseSpatial.SliceThickness = 2.5;  % Common default
                    doseSpatial.NumSlices = size(referenceDose, 3);
                    doseSpatial.GridFrameOffsetVector = (0:doseSpatial.NumSlices-1) * doseSpatial.SliceThickness;
                    fprintf('  - WARNING: Z resolution not found, using default: %.3f mm\n', doseSpatial.SliceThickness);
                end
                
                % Calculate RTDOSE bounding box in patient coordinates
                doseSpatial.XMin = doseSpatial.ImagePositionPatient(1);
                doseSpatial.XMax = doseSpatial.XMin + (doseSpatial.Columns - 1) * doseSpatial.PixelSpacing(2);
                doseSpatial.YMin = doseSpatial.ImagePositionPatient(2);
                doseSpatial.YMax = doseSpatial.YMin + (doseSpatial.Rows - 1) * doseSpatial.PixelSpacing(1);
                doseSpatial.ZMin = doseSpatial.ImagePositionPatient(3);
                doseSpatial.ZMax = doseSpatial.ZMin + doseSpatial.GridFrameOffsetVector(end);
                
                % Store dose grid parameters for later use
                doseGrid.resolution = [doseSpatial.PixelSpacing(1), doseSpatial.PixelSpacing(2), doseSpatial.SliceThickness];
                doseGrid.dimensions = [doseSpatial.Rows, doseSpatial.Columns, doseSpatial.NumSlices];
                doseGrid.ImagePositionPatient = doseSpatial.ImagePositionPatient;
                
                fprintf('    RTDOSE spatial extent:\n');
                fprintf('      X: [%.2f, %.2f] mm\n', doseSpatial.XMin, doseSpatial.XMax);
                fprintf('      Y: [%.2f, %.2f] mm\n', doseSpatial.YMin, doseSpatial.YMax);
                fprintf('      Z: [%.2f, %.2f] mm\n', doseSpatial.ZMin, doseSpatial.ZMax);
                fprintf('      Resolution: [%.3f, %.3f, %.3f] mm\n', ...
                    doseGrid.resolution(1), doseGrid.resolution(2), doseGrid.resolution(3));
                fprintf('      Dimensions: [%d, %d, %d]\n', ...
                    doseGrid.dimensions(1), doseGrid.dimensions(2), doseGrid.dimensions(3));
                
            else
                fprintf('  WARNING: No RTDOSE file found\n');
                doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
                doseGrid.dimensions = ct.cubeDim;
                referenceDose = [];
            end
            
        catch ME
            fprintf('WARNING: Could not load RTDOSE: %s\n', ME.message);
            doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
            doseGrid.dimensions = ct.cubeDim;
            referenceDose = [];
        end
        
        %% Step 3: Extract segment information from RTPLAN
        fprintf('\n[3/8] Extracting segment information from RTPLAN...\n');
        
        % Initialize segment data structure
        segmentData = struct();
        segmentData.beams = {};
        
        try
            % Find RTPLAN file
            rtplanFile = dir(fullfile(dicomPath, 'RP*.dcm'));
            if isempty(rtplanFile)
                rtplanFile = dir(fullfile(dicomPath, '*RTPLAN*.dcm'));
            end
            
            if isempty(rtplanFile)
                error('No RTPLAN file found in %s', dicomPath);
            end
            
            rtplanInfo = dicominfo(fullfile(rtplanFile(1).folder, rtplanFile(1).name));
            fprintf('  - RTPLAN file: %s\n', rtplanFile(1).name);
            
            % Get beam metersets from FractionGroupSequence
            beamMetersets = [];
            numFractions = 1;
            if isfield(rtplanInfo, 'FractionGroupSequence')
                fg = rtplanInfo.FractionGroupSequence.Item_1;
                
                if isfield(fg, 'NumberOfFractionsPlanned')
                    numFractions = fg.NumberOfFractionsPlanned;
                    fprintf('  - Number of fractions planned: %d\n', numFractions);
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
                    fprintf('  - Extracted beam metersets for %d beams\n', numRefBeams);
                end
            end
            
            % Process BeamSequence to extract control points
            if ~isfield(rtplanInfo, 'BeamSequence')
                error('No BeamSequence found in RTPLAN');
            end
            
            numBeams = length(fieldnames(rtplanInfo.BeamSequence));
            fprintf('  - Processing %d beams...\n\n', numBeams);
            
            totalSegments = 0;
            
            for beamIdx = 1:numBeams
                beamField = sprintf('Item_%d', beamIdx);
                beam = rtplanInfo.BeamSequence.(beamField);
                
                % Initialize beam data structure
                beamData = struct();
                beamData.beamIdx = beamIdx;
                beamData.beamName = '';
                beamData.gantryAngle = 0;
                beamData.couchAngle = 0;
                beamData.collimatorAngle = 0;
                beamData.beamMeterset = 0;
                beamData.segments = {};
                
                % Get beam name
                if isfield(beam, 'BeamName')
                    beamData.beamName = beam.BeamName;
                elseif isfield(beam, 'BeamDescription')
                    beamData.beamName = beam.BeamDescription;
                else
                    beamData.beamName = sprintf('Beam_%d', beamIdx);
                end
                
                % Get beam meterset
                if beamIdx <= length(beamMetersets)
                    beamData.beamMeterset = beamMetersets(beamIdx);
                end
                
                fprintf('  Beam %d: %s (Meterset: %.2f MU)\n', beamIdx, beamData.beamName, beamData.beamMeterset);
                
                % Get beam limiting device sequence for leaf boundary definitions
                leafBoundaries = [];
                numLeafPairs = 0;
                if isfield(beam, 'BeamLimitingDeviceSequence')
                    numDevices = length(fieldnames(beam.BeamLimitingDeviceSequence));
                    for devIdx = 1:numDevices
                        devField = sprintf('Item_%d', devIdx);
                        device = beam.BeamLimitingDeviceSequence.(devField);
                        if isfield(device, 'RTBeamLimitingDeviceType')
                            deviceType = device.RTBeamLimitingDeviceType;
                            if contains(deviceType, 'MLC', 'IgnoreCase', true)
                                if isfield(device, 'LeafPositionBoundaries')
                                    leafBoundaries = device.LeafPositionBoundaries;
                                    numLeafPairs = length(leafBoundaries) - 1;
                                    fprintf('    - MLC: %d leaf pairs\n', numLeafPairs);
                                elseif isfield(device, 'NumberOfLeafJawPairs')
                                    numLeafPairs = device.NumberOfLeafJawPairs;
                                    fprintf('    - MLC: %d leaf pairs (no boundaries)\n', numLeafPairs);
                                end
                            end
                        end
                    end
                end
                beamData.leafBoundaries = leafBoundaries;
                beamData.numLeafPairs = numLeafPairs;
                
                % Process control points to extract segments
                if ~isfield(beam, 'ControlPointSequence')
                    fprintf('    WARNING: No ControlPointSequence found\n');
                    segmentData.beams{beamIdx} = beamData;
                    continue;
                end
                
                numCP = length(fieldnames(beam.ControlPointSequence));
                fprintf('    - Control points: %d\n', numCP);
                
                % Extract cumulative meterset weights and MLC positions for all control points
                controlPoints = cell(numCP, 1);
                
                for cpIdx = 1:numCP
                    cpField = sprintf('Item_%d', cpIdx);
                    cp = beam.ControlPointSequence.(cpField);
                    
                    cpData = struct();
                    cpData.cpIdx = cpIdx;
                    cpData.cumulativeMetersetWeight = 0;
                    cpData.gantryAngle = [];
                    cpData.collimatorAngle = [];
                    cpData.mlcPositions = [];
                    cpData.jawX = [];
                    cpData.jawY = [];
                    
                    % Get cumulative meterset weight
                    if isfield(cp, 'CumulativeMetersetWeight')
                        cpData.cumulativeMetersetWeight = cp.CumulativeMetersetWeight;
                    end
                    
                    % Get gantry angle (first CP has it, subsequent may inherit)
                    if isfield(cp, 'GantryAngle')
                        cpData.gantryAngle = cp.GantryAngle;
                        beamData.gantryAngle = cp.GantryAngle;  % Store first gantry angle
                    end
                    
                    % Get collimator angle
                    if isfield(cp, 'BeamLimitingDeviceAngle')
                        cpData.collimatorAngle = cp.BeamLimitingDeviceAngle;
                        beamData.collimatorAngle = cp.BeamLimitingDeviceAngle;
                    end
                    
                    % Get couch angle
                    if isfield(cp, 'PatientSupportAngle')
                        beamData.couchAngle = cp.PatientSupportAngle;
                    end
                    
                    % Get beam limiting device positions (MLC and jaws)
                    if isfield(cp, 'BeamLimitingDevicePositionSequence')
                        numDevices = length(fieldnames(cp.BeamLimitingDevicePositionSequence));
                        
                        for devIdx = 1:numDevices
                            devField = sprintf('Item_%d', devIdx);
                            device = cp.BeamLimitingDevicePositionSequence.(devField);
                            
                            if ~isfield(device, 'RTBeamLimitingDeviceType') || ...
                               ~isfield(device, 'LeafJawPositions')
                                continue;
                            end
                            
                            deviceType = device.RTBeamLimitingDeviceType;
                            positions = device.LeafJawPositions;
                            
                            % Classify device type
                            if contains(deviceType, 'ASYMX', 'IgnoreCase', true) || strcmp(deviceType, 'X')
                                cpData.jawX = positions;
                            elseif contains(deviceType, 'ASYMY', 'IgnoreCase', true) || strcmp(deviceType, 'Y')
                                cpData.jawY = positions;
                            elseif contains(deviceType, 'MLC', 'IgnoreCase', true)
                                cpData.mlcPositions = positions;
                            end
                        end
                    end
                    
                    controlPoints{cpIdx} = cpData;
                end
                
                % Calculate segment weights from cumulative meterset weights
                % Segment i is between control point i and i+1
                % Weight for segment i = CumulativeMetersetWeight[i+1] - CumulativeMetersetWeight[i]
                
                numSegments = numCP - 1;  % Number of segments = number of control points - 1
                
                if numSegments < 1
                    fprintf('    WARNING: Not enough control points to form segments\n');
                    segmentData.beams{beamIdx} = beamData;
                    continue;
                end
                
                fprintf('    - Segments: %d\n', numSegments);
                
                % Create segments
                for segIdx = 1:numSegments
                    segment = struct();
                    segment.segmentIdx = segIdx;
                    segment.beamIdx = beamIdx;
                    
                    % Starting control point
                    cpStart = controlPoints{segIdx};
                    % Ending control point
                    cpEnd = controlPoints{segIdx + 1};
                    
                    % Calculate segment weight (fraction of beam dose)
                    segment.segmentWeight = cpEnd.cumulativeMetersetWeight - cpStart.cumulativeMetersetWeight;
                    
                    % Calculate segment meterset (MU)
                    segment.segmentMeterset = segment.segmentWeight * beamData.beamMeterset;
                    
                    % Get gantry angle for this segment
                    % For step-and-shoot IMRT, gantry is typically fixed
                    % For VMAT, we'd use the average or starting angle
                    if ~isempty(cpStart.gantryAngle)
                        segment.gantryAngle = cpStart.gantryAngle;
                    else
                        segment.gantryAngle = beamData.gantryAngle;
                    end
                    
                    % Get collimator angle
                    if ~isempty(cpStart.collimatorAngle)
                        segment.collimatorAngle = cpStart.collimatorAngle;
                    else
                        segment.collimatorAngle = beamData.collimatorAngle;
                    end
                    
                    % For step-and-shoot: use MLC positions from starting control point
                    % For sliding window/VMAT: you'd want to interpolate or use average
                    segment.mlcPositions = cpStart.mlcPositions;
                    segment.jawX = cpStart.jawX;
                    segment.jawY = cpStart.jawY;
                    
                    % Also store ending positions for reference (useful for sliding window)
                    segment.mlcPositionsEnd = cpEnd.mlcPositions;
                    
                    % Store control point indices
                    segment.cpStart = segIdx;
                    segment.cpEnd = segIdx + 1;
                    
                    beamData.segments{segIdx} = segment;
                    totalSegments = totalSegments + 1;
                end
                
                % Print segment summary
                fprintf('    Segment weights: ');
                weights = cellfun(@(s) s.segmentWeight, beamData.segments);
                fprintf('[%.4f', weights(1));
                for i = 2:min(5, length(weights))
                    fprintf(', %.4f', weights(i));
                end
                if length(weights) > 5
                    fprintf(', ... (%.4f total)', sum(weights));
                else
                    fprintf('] (%.4f total)', sum(weights));
                end
                fprintf('\n');
                
                segmentData.beams{beamIdx} = beamData;
            end
            
            fprintf('\n  Total segments across all beams: %d\n', totalSegments);
            segmentData.totalSegments = totalSegments;
            segmentData.numBeams = numBeams;
            segmentData.numFractions = numFractions;
            
        catch ME
            fprintf('ERROR extracting segment data: %s\n', ME.message);
            fprintf('Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            continue;
        end
        
        %% Step 4: Configure dose calculation
        fprintf('\n[4/8] Configuring dose calculation...\n');
        
        % Check and override machine specification
        if isfield(pln, 'machine')
            fprintf('  - Original machine: %s\n', pln.machine);
        end
        
        % Check available machine files
        machineDir = fullfile(matradPath, 'basedata');
        if exist(machineDir, 'dir')
            machineFiles = dir(fullfile(machineDir, 'photons*.mat'));
            if ~isempty(machineFiles)
                fprintf('  - Available machine files:\n');
                for i = 1:min(5, length(machineFiles))
                    [~, machineName, ~] = fileparts(machineFiles(i).name);
                    fprintf('    %d. %s\n', i, machineName);
                end
                
                % Use first generic machine or specific one
                genericMachines = machineFiles(contains({machineFiles.name}, 'Generic', 'IgnoreCase', true));
                if ~isempty(genericMachines)
                    [~, selectedMachine, ~] = fileparts(genericMachines(1).name);
                    fprintf('  - Using generic machine: %s\n', selectedMachine);
                    pln.machine = selectedMachine;
                else
                    % Use first available photon machine
                    [~, selectedMachine, ~] = fileparts(machineFiles(1).name);
                    fprintf('  - Using machine: %s\n', selectedMachine);
                    pln.machine = selectedMachine;
                end
            else
                fprintf('  - Warning: No machine files found, using default\n');
                pln.machine = 'Generic';
            end
        else
            fprintf('  - Warning: Machine directory not found, using default\n');
            pln.machine = 'Generic';
        end
        
        % Set up dose calculation parameters
        pln.propStf.bixelWidth = 5; % mm
        pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x;
        pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y;
        pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z;
        
        % Set algorithm
        if strcmp(pln.radiationMode, 'photons')
            pln.propDoseCalc.engine = 'pencilBeam';
        else
            pln.propDoseCalc.engine = 'matRad_pencilBeam';
        end
        
        fprintf('  - Radiation mode: %s\n', pln.radiationMode);
        fprintf('  - Machine: %s\n', pln.machine);
        fprintf('  - Dose engine: %s\n', pln.propDoseCalc.engine);
        fprintf('  - Calculation grid resolution: [%.2f, %.2f, %.2f] mm\n', ...
            ct.resolution.x, ct.resolution.y, ct.resolution.z);
        
        %% Step 5: Generate steering information
        fprintf('\n[5/8] Generating steering information...\n');
        
        try
            stf = matRad_generateStf(ct, cst, pln);
            fprintf('  - Steering file generated for %d beams\n', length(stf));
            
            % Display beam information
            for i = 1:length(stf)
                fprintf('    Beam %d: Gantry=%.1f, Couch=%.1f, Rays=%d, Bixels=%d\n', ...
                    i, stf(i).gantryAngle, stf(i).couchAngle, stf(i).numOfRays, ...
                    stf(i).numOfRays * length(stf(i).ray(1).energy));
            end
            
        catch ME
            fprintf('Error generating steering file: %s\n', ME.message);
            continue;
        end
        
        %% Step 6: Calculate segment-by-segment doses
        fprintf('\n[6/8] Calculating segment-by-segment doses...\n');
        
        % Initialize storage
        segmentDoses = {};  % Cell array to store all segment doses
        beamDoses = cell(length(stf), 1);  % Accumulated dose per beam
        totalDose = [];  % Total dose across all segments
        calculatedGridSize = [];
        
        segmentCounter = 0;
        
        for beamIdx = 1:length(stf)
            beamData = segmentData.beams{beamIdx};
            numSegments = length(beamData.segments);
            
            fprintf('\n  Processing Beam %d/%d (%s): %d segments\n', ...
                beamIdx, length(stf), beamData.beamName, numSegments);
            fprintf('    Gantry: %.1f°, Couch: %.1f°, Meterset: %.2f MU\n', ...
                beamData.gantryAngle, beamData.couchAngle, beamData.beamMeterset);
            
            % Create single-beam plan and steering file
            plnSingle = pln;
            plnSingle.propStf.numOfBeams = 1;
            plnSingle.propStf.isoCenter = stf(beamIdx).isoCenter;
            stfSingle = stf(beamIdx);
            
            % Calculate dose influence matrix for this beam (done once per beam)
            try
                fprintf('    Calculating dose influence matrix...\n');
                dij = matRad_calcDoseInfluence(ct, cst, stfSingle, plnSingle);
                fprintf('    - Total bixels in dij: %d\n', dij.totalNumOfBixels);
                
                % Store dij dimensions for reference
                if isempty(calculatedGridSize)
                    calculatedGridSize = ct.cubeDim;
                    totalDose = zeros(calculatedGridSize);
                end
                
            catch ME
                fprintf('    ERROR in dose influence calculation: %s\n', ME.message);
                fprintf('      %s\n', ME.message);
                continue;
            end
            
            % Initialize beam dose accumulator
            beamDoseAccum = zeros(calculatedGridSize);
            
            % Process each segment
            for segIdx = 1:numSegments
                segmentCounter = segmentCounter + 1;
                segment = beamData.segments{segIdx};
                
                fprintf('    Segment %d/%d (Weight: %.4f, MU: %.2f)...\n', ...
                    segIdx, numSegments, segment.segmentWeight, segment.segmentMeterset);
                
                % Skip segments with zero weight
                if segment.segmentWeight <= 0
                    fprintf('      Skipping (zero weight)\n');
                    continue;
                end
                
                % Create weight vector for this segment
                % For IMRT, we need to determine which bixels are "open" based on MLC
                % and weight them accordingly
                
                w = zeros(dij.totalNumOfBixels, 1);
                
                % Get MLC positions for this segment
                mlcPos = segment.mlcPositions;
                jawX = segment.jawX;
                jawY = segment.jawY;
                
                if ~isempty(mlcPos) && ~isempty(beamData.leafBoundaries)
                    % Determine which bixels are open based on MLC aperture
                    numLeafPairs = beamData.numLeafPairs;
                    leafBoundaries = beamData.leafBoundaries;
                    
                    % MLC positions: first half are left bank (A), second half are right bank (B)
                    leftBank = mlcPos(1:numLeafPairs);
                    rightBank = mlcPos(numLeafPairs+1:end);
                    
                    % For each ray in the steering file, check if it's within the aperture
                    numRays = stfSingle.numOfRays;
                    bixelIdx = 0;
                    
                    for rayIdx = 1:numRays
                        ray = stfSingle.ray(rayIdx);
                        
                        % Ray position at isocenter (in IEC beam coordinates)
                        rayPosX = ray.rayPos_bev(1);  % X position (left-right)
                        rayPosY = ray.rayPos_bev(2);  % Y position (sup-inf)
                        
                        % Check if ray is within jaw opening
                        inJawX = true;
                        inJawY = true;
                        
                        if ~isempty(jawX) && length(jawX) >= 2
                            inJawX = (rayPosX >= jawX(1)) && (rayPosX <= jawX(2));
                        end
                        
                        if ~isempty(jawY) && length(jawY) >= 2
                            inJawY = (rayPosY >= jawY(1)) && (rayPosY <= jawY(2));
                        end
                        
                        % Check if ray is within MLC opening
                        inMLC = false;
                        
                        if inJawX && inJawY
                            % Find which leaf pair this ray corresponds to
                            for leafIdx = 1:numLeafPairs
                                % Leaf boundaries define Y extent of each leaf pair
                                leafYMin = leafBoundaries(leafIdx);
                                leafYMax = leafBoundaries(leafIdx + 1);
                                
                                if rayPosY >= leafYMin && rayPosY < leafYMax
                                    % Ray is in this leaf pair's Y range
                                    % Check X position against MLC opening
                                    leftEdge = leftBank(leafIdx);
                                    rightEdge = rightBank(leafIdx);
                                    
                                    if rayPosX >= leftEdge && rayPosX <= rightEdge
                                        inMLC = true;
                                    end
                                    break;
                                end
                            end
                        end
                        
                        % Set bixel weights
                        numEnergies = length(ray.energy);
                        for energyIdx = 1:numEnergies
                            bixelIdx = bixelIdx + 1;
                            if inMLC
                                % Weight by segment weight (fraction of beam dose)
                                w(bixelIdx) = segment.segmentWeight * beamData.beamMeterset;
                            end
                        end
                    end
                    
                    numOpenBixels = nnz(w);
                    fprintf('      Open bixels: %d/%d (%.1f%%)\n', ...
                        numOpenBixels, dij.totalNumOfBixels, 100*numOpenBixels/dij.totalNumOfBixels);
                    
                else
                    % No MLC data - use uniform weights
                    fprintf('      WARNING: No MLC aperture data, using uniform weights\n');
                    w(:) = segment.segmentWeight * beamData.beamMeterset / dij.totalNumOfBixels;
                end
                
                % Calculate dose for this segment
                try
                    if any(w > 0)
                        resultGUI = matRad_calcDoseForward(ct, cst, stfSingle, plnSingle, w);
                        
                        if isfield(resultGUI, 'physicalDose')
                            segmentDose = resultGUI.physicalDose;
                            maxDose = max(segmentDose(:));
                            
                            fprintf('      Max dose: %.4f Gy\n', maxDose);
                            
                            % Store segment dose
                            segmentResult = struct();
                            segmentResult.beamIdx = beamIdx;
                            segmentResult.segmentIdx = segIdx;
                            segmentResult.segmentWeight = segment.segmentWeight;
                            segmentResult.segmentMeterset = segment.segmentMeterset;
                            segmentResult.gantryAngle = segment.gantryAngle;
                            segmentResult.physicalDose = segmentDose;
                            segmentResult.maxDose = maxDose;
                            segmentDoses{segmentCounter} = segmentResult;
                            
                            % Accumulate beam dose
                            beamDoseAccum = beamDoseAccum + segmentDose;
                            
                            % Accumulate total dose
                            totalDose = totalDose + segmentDose;
                        else
                            fprintf('      ERROR: No physicalDose in result\n');
                        end
                    else
                        fprintf('      Skipping (no open bixels)\n');
                    end
                    
                catch ME
                    fprintf('      ERROR in dose calculation: %s\n', ME.message);
                end
            end
            
            % Store accumulated beam dose
            beamDoses{beamIdx} = struct();
            beamDoses{beamIdx}.beamIdx = beamIdx;
            beamDoses{beamIdx}.beamName = beamData.beamName;
            beamDoses{beamIdx}.gantryAngle = beamData.gantryAngle;
            beamDoses{beamIdx}.couchAngle = beamData.couchAngle;
            beamDoses{beamIdx}.beamMeterset = beamData.beamMeterset;
            beamDoses{beamIdx}.numSegments = numSegments;
            beamDoses{beamIdx}.physicalDose = beamDoseAccum;
            beamDoses{beamIdx}.maxDose = max(beamDoseAccum(:));
            
            fprintf('    Beam %d total max dose: %.4f Gy\n', beamIdx, max(beamDoseAccum(:)));
        end
        
        numSuccessfulSegments = sum(~cellfun(@isempty, segmentDoses));
        fprintf('\n  Summary: %d/%d segments calculated successfully\n', ...
            numSuccessfulSegments, segmentData.totalSegments);
        fprintf('  Total dose max: %.4f Gy\n', max(totalDose(:)));
        
        %% Step 7: Resample to RTDOSE grid
        fprintf('\n[7/8] Resampling doses to RTDOSE grid...\n');
        
        if ~isempty(referenceDose) && ~isempty(ctSpatial.ImagePositionPatient) && ~isempty(doseSpatial.ImagePositionPatient)
            
            fprintf('  Computing spatial transformation from CT to RTDOSE grid...\n');
            
            % CT coordinate vectors (patient coordinates)
            ct_x = double(ctSpatial.ImagePositionPatient(1)) + double((0:ctSpatial.Columns-1) * ctSpatial.PixelSpacing(2));
            ct_y = double(ctSpatial.ImagePositionPatient(2)) + double((0:ctSpatial.Rows-1) * ctSpatial.PixelSpacing(1));
            ct_z = double(ctSpatial.SlicePositions);
            
            % RTDOSE coordinate vectors (patient coordinates)
            dose_x = double(doseSpatial.ImagePositionPatient(1)) + double((0:doseSpatial.Columns-1) * doseSpatial.PixelSpacing(2));
            dose_y = double(doseSpatial.ImagePositionPatient(2)) + double((0:doseSpatial.Rows-1) * doseSpatial.PixelSpacing(1));
            dose_z = double(doseSpatial.ImagePositionPatient(3)) + double(doseSpatial.GridFrameOffsetVector);
            
            % Create meshgrids
            [CT_Y, CT_X, CT_Z] = ndgrid(ct_y, ct_x, ct_z);
            CT_Y = double(CT_Y);
            CT_X = double(CT_X);
            CT_Z = double(CT_Z);
            
            [DOSE_Y, DOSE_X, DOSE_Z] = ndgrid(dose_y, dose_x, dose_z);
            
            % Resample CT
            fprintf('  Resampling CT to RTDOSE grid...\n');
            ctCube = ct.cubeHU{1};
            try
                ctResampled = interpn(CT_Y, CT_X, CT_Z, ctCube, DOSE_Y, DOSE_X, DOSE_Z, 'linear', -1000);
                fprintf('    CT resampled: %d x %d x %d\n', size(ctResampled));
            catch ME
                fprintf('    ERROR resampling CT: %s\n', ME.message);
                ctResampled = [];
            end
            
            % Resample total dose
            fprintf('  Resampling total dose to RTDOSE grid...\n');
            try
                totalDoseResampled = interpn(CT_Y, CT_X, CT_Z, totalDose, DOSE_Y, DOSE_X, DOSE_Z, 'linear', 0);
                fprintf('    Total dose resampled: %d x %d x %d, Max=%.4f Gy\n', ...
                    size(totalDoseResampled), max(totalDoseResampled(:)));
            catch ME
                fprintf('    ERROR resampling total dose: %s\n', ME.message);
                totalDoseResampled = totalDose;
            end
            
            % Resample beam doses
            fprintf('  Resampling beam doses...\n');
            beamDosesResampled = cell(length(beamDoses), 1);
            for beamIdx = 1:length(beamDoses)
                if ~isempty(beamDoses{beamIdx})
                    try
                        beamDoseResampled = interpn(CT_Y, CT_X, CT_Z, beamDoses{beamIdx}.physicalDose, ...
                            DOSE_Y, DOSE_X, DOSE_Z, 'linear', 0);
                        
                        beamDosesResampled{beamIdx} = beamDoses{beamIdx};
                        beamDosesResampled{beamIdx}.physicalDose = beamDoseResampled;
                        beamDosesResampled{beamIdx}.maxDose = max(beamDoseResampled(:));
                        
                        fprintf('    Beam %d: Max=%.4f Gy\n', beamIdx, max(beamDoseResampled(:)));
                    catch ME
                        fprintf('    ERROR resampling beam %d: %s\n', beamIdx, ME.message);
                        beamDosesResampled{beamIdx} = beamDoses{beamIdx};
                    end
                end
            end
            
            % Resample segment doses (optional - can be large)
            fprintf('  Resampling segment doses...\n');
            segmentDosesResampled = cell(size(segmentDoses));
            for i = 1:length(segmentDoses)
                if ~isempty(segmentDoses{i})
                    try
                        segDoseResampled = interpn(CT_Y, CT_X, CT_Z, segmentDoses{i}.physicalDose, ...
                            DOSE_Y, DOSE_X, DOSE_Z, 'linear', 0);
                        
                        segmentDosesResampled{i} = segmentDoses{i};
                        segmentDosesResampled{i}.physicalDose = segDoseResampled;
                        segmentDosesResampled{i}.maxDose = max(segDoseResampled(:));
                    catch
                        segmentDosesResampled{i} = segmentDoses{i};
                    end
                end
            end
            fprintf('    Segment doses resampled: %d segments\n', sum(~cellfun(@isempty, segmentDosesResampled)));
            
            % Create resampled CT structure
            ctResampled_struct = struct();
            ctResampled_struct.cubeHU = {ctResampled};
            ctResampled_struct.cubeDim = [size(ctResampled, 1), size(ctResampled, 2), size(ctResampled, 3)];
            ctResampled_struct.resolution.x = doseGrid.resolution(2);
            ctResampled_struct.resolution.y = doseGrid.resolution(1);
            ctResampled_struct.resolution.z = doseGrid.resolution(3);
            ctResampled_struct.x = dose_x;
            ctResampled_struct.y = dose_y;
            ctResampled_struct.z = dose_z;
            ctResampled_struct.ImagePositionPatient = doseSpatial.ImagePositionPatient;
            
        else
            fprintf('  Cannot resample - missing RTDOSE or spatial coordinates\n');
            totalDoseResampled = totalDose;
            beamDosesResampled = beamDoses;
            segmentDosesResampled = segmentDoses;
            ctResampled_struct = [];
        end
        
        %% Step 8: Save results and compare
        fprintf('\n[8/8] Saving results...\n');
        
        % Save segment data (extracted from RTPLAN)
        save(fullfile(outputPath, 'segmentData.mat'), 'segmentData');
        fprintf('  - Segment data saved: segmentData.mat\n');
        
        % Save segment doses (CT grid)
        save(fullfile(outputPath, 'segmentDoses_CTgrid.mat'), 'segmentDoses', 'beamDoses', ...
             'totalDose', 'stf', 'pln', 'ct', 'cst', '-v7.3');
        fprintf('  - Segment doses (CT grid) saved: segmentDoses_CTgrid.mat\n');
        
        % Save resampled doses (RTDOSE grid)
        if ~isempty(ctResampled_struct)
            save(fullfile(outputPath, 'segmentDoses.mat'), 'segmentDosesResampled', 'beamDosesResampled', ...
                 'totalDoseResampled', 'referenceDose', 'ctResampled_struct', 'doseGrid', ...
                 'ctSpatial', 'doseSpatial', 'stf', 'pln', '-v7.3');
            fprintf('  - Segment doses (RTDOSE grid) saved: segmentDoses.mat\n');
            
            save(fullfile(outputPath, 'ctResampled.mat'), 'ctResampled_struct', 'doseGrid', ...
                 'ctSpatial', 'doseSpatial');
            fprintf('  - Resampled CT saved: ctResampled.mat\n');
        end
        
        % Compare with reference
        if ~isempty(referenceDose) && ~isempty(totalDoseResampled)
            fprintf('\n  Comparing with reference RTDOSE...\n');
            
            if isequal(size(totalDoseResampled), size(referenceDose))
                doseDiff = totalDoseResampled - referenceDose;
                
                fprintf('    Grids match! Direct comparison:\n');
                fprintf('      Calculated max dose: %.4f Gy\n', max(totalDoseResampled(:)));
                fprintf('      Reference max dose:  %.4f Gy\n', max(referenceDose(:)));
                fprintf('      Mean absolute difference: %.4f Gy\n', mean(abs(doseDiff(:))));
                fprintf('      Max difference: %.4f Gy\n', max(abs(doseDiff(:))));
                fprintf('      RMS difference: %.4f Gy\n', sqrt(mean(doseDiff(:).^2)));
                
                % Calculate relative differences in high dose region
                highDoseThreshold = 0.5 * max(referenceDose(:));
                highDoseMask = referenceDose > highDoseThreshold;
                if any(highDoseMask(:))
                    relativeDiff = abs(doseDiff(highDoseMask)) ./ referenceDose(highDoseMask) * 100;
                    fprintf('      Mean relative diff (>50%% max): %.2f%%\n', mean(relativeDiff));
                end
                
                % Save comparison
                comparison = struct();
                comparison.calculated = totalDoseResampled;
                comparison.reference = referenceDose;
                comparison.difference = doseDiff;
                comparison.gridResolution = doseGrid.resolution;
                comparison.gridDimensions = doseGrid.dimensions;
                comparison.metrics.meanAbsDiff = mean(abs(doseDiff(:)));
                comparison.metrics.maxDiff = max(abs(doseDiff(:)));
                comparison.metrics.rmsDiff = sqrt(mean(doseDiff(:).^2));
                if any(highDoseMask(:))
                    comparison.metrics.meanRelativeDiff = mean(relativeDiff);
                end
                
                save(fullfile(outputPath, 'doseComparison.mat'), 'comparison');
                fprintf('    Comparison saved: doseComparison.mat\n');
            else
                fprintf('    WARNING: Grid size mismatch\n');
                fprintf('      Calculated: %d x %d x %d\n', size(totalDoseResampled));
                fprintf('      Reference:  %d x %d x %d\n', size(referenceDose));
            end
        end
        
        % Export individual beam doses to .mat files
        fprintf('\n  Exporting individual beam doses...\n');
        for beamIdx = 1:length(beamDosesResampled)
            if ~isempty(beamDosesResampled{beamIdx})
                beamFilename = sprintf('Beam_%02d.mat', beamIdx);
                beamData = beamDosesResampled{beamIdx};
                save(fullfile(outputPath, beamFilename), 'beamData');
            end
        end
        
        fprintf('\n========================================\n');
        fprintf('Processing complete for %s - %s\n', currentID, currentSession);
        fprintf('Results saved to: %s\n', outputPath);
        fprintf('========================================\n');
        
        fprintf('\n** SUMMARY **\n');
        fprintf('  Patient: %s\n', currentID);
        fprintf('  Session: %s\n', currentSession);
        fprintf('  Beams: %d\n', segmentData.numBeams);
        fprintf('  Total segments: %d\n', segmentData.totalSegments);
        fprintf('  Segments calculated: %d\n', numSuccessfulSegments);
        fprintf('  Total dose max: %.4f Gy\n', max(totalDoseResampled(:)));
        if ~isempty(referenceDose)
            fprintf('  Reference dose max: %.4f Gy\n', max(referenceDose(:)));
        end
        fprintf('\n');
        
        fprintf('** OUTPUT FILES **\n');
        fprintf('  - segmentData.mat: Segment information extracted from RTPLAN\n');
        fprintf('  - segmentDoses.mat: All segment and beam doses (RTDOSE grid)\n');
        fprintf('  - segmentDoses_CTgrid.mat: Doses on original CT grid\n');
        fprintf('  - ctResampled.mat: CT resampled to RTDOSE grid\n');
        fprintf('  - doseComparison.mat: Comparison with reference dose\n');
        fprintf('  - Beam_XX.mat: Individual beam doses\n\n');
        
    end
end

fprintf('\n\nAll processing complete!\n');
