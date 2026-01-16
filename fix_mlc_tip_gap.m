%% fix_mlc_tip_gap.m
% Script to fix minimum dynamic tip gap issues in RTPLAN DICOM files
% 
% This script:
% 1. Reads an RTPLAN DICOM file
% 2. Loops through all beams and control points (sorted by ControlPointIndex)
% 3. Identifies MLC leaf pairs that are DYNAMIC (move during beam delivery)
% 4. For dynamic leaves with gaps less than the minimum threshold, expands tips
% 5. Applies bounds checking to ensure positions stay within [-140, 140] mm
% 6. Verifies the corrections
% 7. Exports the corrected plan as a new DICOM file
%
% IMPORTANT NOTES:
% - All RTPLAN positions are in millimeters (mm)
% - Control points are sorted by ControlPointIndex to ensure correct processing
% - Only dynamic leaves are corrected; static leaves are left unchanged
% - MLC positions are clamped to physical limits if expansion would exceed bounds
%
% Author: Claude
% Date: January 2026

clear; clc;

%% ===================== USER PARAMETERS =====================
% Working directory
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';

% Patient IDs and sessions to process
patientIDs = {'1194203'};  % Add more IDs as needed
sessions = {'Session_1'};  % Add more sessions as needed

% Minimum tip gap threshold (mm) - from commissioning
% Note: RTPLAN uses millimeters for all positions
MIN_TIP_GAP = 0.6;  % mm (0.06 cm)

% Expansion amount per side when gap is too small (mm)
EXPANSION_PER_SIDE = 0.4;  % mm (0.04 cm - total expansion will be 0.8 mm)

% MLC position bounds (mm) - Halcyon physical limits
MLC_MIN_POSITION = -140;  % mm
MLC_MAX_POSITION = 140;   % mm

% Output filename suffix
OUTPUT_SUFFIX = '_adjusted_mlc';

%% ===================== LOOP OVER PATIENTS AND SESSIONS =====================
for p = 1:length(patientIDs)
    patientID = patientIDs{p};
    
    for s = 1:length(sessions)
        sessionName = sessions{s};
        
        fprintf('=============================================================\n');
        fprintf('MLC Tip Gap Correction Script\n');
        fprintf('Patient ID: %s, Session: %s\n', patientID, sessionName);
        fprintf('=============================================================\n\n');
        
        % Construct path to sct directory
        rawwd = fullfile(wd, 'EthosExports', patientID, 'Pancreas', sessionName);
        sctDir = fullfile(rawwd, 'sct');
        
        % Check if directory exists
        if ~exist(sctDir, 'dir')
            fprintf('WARNING: Directory does not exist: %s\n', sctDir);
            fprintf('Skipping patient %s, session %s\n\n', patientID, sessionName);
            continue;
        end
        
        % Find RTPLAN file in sct directory
        fprintf('Scanning directory: %s\n', sctDir);
        dcmFiles = dir(fullfile(sctDir, '*.dcm'));
        
        inputFile = '';
        for i = 1:length(dcmFiles)
            testFile = fullfile(sctDir, dcmFiles(i).name);
            try
                info = dicominfo(testFile);
                % Check if this is an RTPLAN file (Modality = RTPLAN)
                if isfield(info, 'Modality') && strcmp(info.Modality, 'RTPLAN')
                    inputFile = testFile;
                    fprintf('Found RTPLAN file: %s\n', dcmFiles(i).name);
                    break;
                end
            catch
                % Not a valid DICOM or can't read - skip
                continue;
            end
        end
        
        if isempty(inputFile)
            fprintf('WARNING: No RTPLAN file found in %s\n', sctDir);
            fprintf('Skipping patient %s, session %s\n\n', patientID, sessionName);
            continue;
        end

%% ===================== READ DICOM RTPLAN =====================
fprintf('=============================================================\n');
fprintf('MLC Tip Gap Correction Script\n');
fprintf('=============================================================\n\n');

fprintf('Reading RTPLAN file: %s\n', inputFile);

% Read the DICOM file
try
    rtplan = dicominfo(inputFile);
    fprintf('Successfully loaded RTPLAN: %s\n', rtplan.RTPlanLabel);
catch ME
    fprintf('ERROR: Failed to read DICOM file: %s\n', ME.message);
    fprintf('Skipping patient %s, session %s\n\n', patientID, sessionName);
    continue;
end

%% ===================== IDENTIFY MLC DEVICE TYPE =====================
% MLC device types in DICOM
MLC_TYPES = {'MLCX', 'MLCY', 'MLCX1', 'MLCX2', 'MLCY1', 'MLCY2'};

%% ===================== PROCESS EACH BEAM =====================
% Get number of beams
beamSeqFields = fieldnames(rtplan.BeamSequence);
numBeams = length(beamSeqFields);

fprintf('\nProcessing %d beams...\n', numBeams);
fprintf('Minimum tip gap threshold: %.2f mm\n', MIN_TIP_GAP);
fprintf('Expansion per side: %.2f mm\n', EXPANSION_PER_SIDE);
fprintf('MLC position bounds: [%.1f, %.1f] mm\n\n', MLC_MIN_POSITION, MLC_MAX_POSITION);

% Initialize counters for verification
totalControlPoints = 0;
totalLeafPairsChecked = 0;
totalCorrections = 0;
correctionLog = {};

% Loop through each beam
for beamIdx = 1:numBeams
    beamField = beamSeqFields{beamIdx};
    beam = rtplan.BeamSequence.(beamField);
    
    beamName = beam.BeamName;
    beamNumber = beam.BeamNumber;
    
    fprintf('Processing Beam %d: %s (Beam Number: %d)\n', beamIdx, beamName, beamNumber);
    
    % Get control point sequence
    if ~isfield(beam, 'ControlPointSequence')
        fprintf('  No ControlPointSequence found - skipping beam\n');
        continue;
    end
    
    cpSeqFields = fieldnames(beam.ControlPointSequence);
    numControlPoints = length(cpSeqFields);
    
    fprintf('  DEBUG: Beam %d has %d control points before sorting\n', beamNumber, numControlPoints);
    
    % Sort control point fields by their ControlPointIndex to ensure proper order
    cpIndices = zeros(numControlPoints, 1);
    for i = 1:numControlPoints
        cp = beam.ControlPointSequence.(cpSeqFields{i});
        if isfield(cp, 'ControlPointIndex')
            cpIndices(i) = cp.ControlPointIndex;
        else
            cpIndices(i) = i - 1;  % Default if not present
        end
    end
    
    fprintf('  DEBUG: Control point indices before sorting: ');
    fprintf('%d ', cpIndices);
    fprintf('\n');
    
    [sortedIndices, sortOrder] = sort(cpIndices);
    cpSeqFields = cpSeqFields(sortOrder);
    
    fprintf('  DEBUG: Control point indices after sorting: ');
    fprintf('%d ', sortedIndices);
    fprintf('\n');
    fprintf('  DEBUG: Sort order applied: ');
    fprintf('%d ', sortOrder);
    fprintf('\n');
    
    %% ===== STEP 1: IDENTIFY DYNAMIC VS STATIC LEAVES =====
    fprintf('  Step 1: Identifying dynamic and static leaves...\n');
    fprintf('  DEBUG: Starting dynamic leaf detection for beam %d\n', beamNumber);
    
    % Store ALL leaf positions across ALL control points for each MLC device
    % to determine which leaves are truly dynamic
    allPositionsMap = containers.Map();
    
    % First pass: Collect all positions from all control points
    fprintf('  DEBUG: First pass - collecting positions from %d control points\n', numControlPoints);
    for cpIdx = 1:numControlPoints
        cpField = cpSeqFields{cpIdx};
        controlPoint = beam.ControlPointSequence.(cpField);
        cpIndex = controlPoint.ControlPointIndex;
        
        if ~isfield(controlPoint, 'BeamLimitingDevicePositionSequence')
            fprintf('  DEBUG: CP %d - No BeamLimitingDevicePositionSequence, skipping\n', cpIndex);
            continue;
        end
        
        bldSeqFields = fieldnames(controlPoint.BeamLimitingDevicePositionSequence);
        fprintf('  DEBUG: CP %d - Found %d beam limiting devices\n', cpIndex, length(bldSeqFields));
        
        for devIdx = 1:length(bldSeqFields)
            devField = bldSeqFields{devIdx};
            device = controlPoint.BeamLimitingDevicePositionSequence.(devField);
            deviceType = device.RTBeamLimitingDeviceType;
            
            if ~any(strcmpi(deviceType, MLC_TYPES))
                fprintf('  DEBUG: CP %d, Device %s - Type %s is not MLC, skipping\n', cpIndex, devField, deviceType);
                continue;
            end
            
            leafPositions = device.LeafJawPositions;
            fprintf('  DEBUG: CP %d, Device %s (Type: %s) - Found %d leaf positions\n', ...
                cpIndex, devField, deviceType, length(leafPositions));
            
            % Initialize storage for this device if first encounter
            if ~allPositionsMap.isKey(devField)
                numLeaves = length(leafPositions);
                % Store positions from each control point as columns
                allPositionsMap(devField) = struct(...
                    'positions', leafPositions(:), ...  % Start with first CP as column
                    'numLeafPairs', numLeaves / 2);
                fprintf('  DEBUG: CP %d, Device %s - Initialized with %d leaves (%d pairs)\n', ...
                    cpIndex, devField, numLeaves, numLeaves/2);
            else
                % Append this control point's positions as a new column
                temp = allPositionsMap(devField);
                temp.positions = [temp.positions, leafPositions(:)];
                allPositionsMap(devField) = temp;
                fprintf('  DEBUG: CP %d, Device %s - Added positions (now have %d CPs)\n', ...
                    cpIndex, devField, size(temp.positions, 2));
            end
        end
    end
    
    % Second pass: Determine which leaves are dynamic (vary across control points)
    fprintf('  DEBUG: Second pass - determining dynamic leaves\n');
    dynamicLeavesMap = containers.Map();
    
    if ~isempty(allPositionsMap.keys)
        for devField = allPositionsMap.keys
            devInfo = allPositionsMap(devField{1});
            allPos = devInfo.positions;  % Each row is a leaf, each column is a CP
            numLeafPairs = devInfo.numLeafPairs;
            numLeaves = size(allPos, 1);
            numCPs = size(allPos, 2);
            
            fprintf('  DEBUG: Device %s - Analyzing %d leaves across %d control points\n', ...
                devField{1}, numLeaves, numCPs);
            
            % Show position data for first few leaves to help debug
            fprintf('  DEBUG: Sample positions for first 3 leaves across all CPs:\n');
            for leafIdx = 1:min(3, numLeaves)
                fprintf('  DEBUG:   Leaf %d positions: [', leafIdx);
                fprintf('%.4f ', allPos(leafIdx, :));
                fprintf(']\n');
            end
            
            % Check if each leaf varies across control points
            isDynamic = false(numLeaves, 1);
            numDynamic = 0;
            
            for leafIdx = 1:numLeaves
                % Get all positions for this leaf across all CPs
                leafPosAcrossCPs = allPos(leafIdx, :);
                % Check if there's any variation (max - min > tolerance)
                minPos = min(leafPosAcrossCPs);
                maxPos = max(leafPosAcrossCPs);
                variation = maxPos - minPos;
                
                if variation > 1e-4  % 0.0001 mm tolerance
                    isDynamic(leafIdx) = true;
                    numDynamic = numDynamic + 1;
                    if leafIdx <= 5 || leafIdx > numLeaves - 5  % Show first and last few
                        fprintf('  DEBUG:   Leaf %d - DYNAMIC (min=%.4f, max=%.4f, var=%.4f mm)\n', ...
                            leafIdx, minPos, maxPos, variation);
                    end
                else
                    if leafIdx <= 5 || leafIdx > numLeaves - 5  % Show first and last few
                        fprintf('  DEBUG:   Leaf %d - STATIC (min=%.4f, max=%.4f, var=%.4f mm)\n', ...
                            leafIdx, minPos, maxPos, variation);
                    end
                end
            end
            
            fprintf('  DEBUG: Device %s - Total dynamic leaves: %d out of %d\n', ...
                devField{1}, numDynamic, numLeaves);
            
            % Store dynamic information
            dynamicLeavesMap(devField{1}) = struct(...
                'isDynamic', isDynamic, ...
                'numLeafPairs', numLeafPairs);
            
            % Report statistics
            bankA_dynamic = isDynamic(1:numLeafPairs);
            bankB_dynamic = isDynamic(numLeafPairs+1:end);
            leafPair_dynamic = bankA_dynamic | bankB_dynamic;
            
            numDynamicPairs = sum(leafPair_dynamic);
            numStaticPairs = numLeafPairs - numDynamicPairs;
            
            fprintf('    Device %s: %d dynamic pairs, %d static pairs (total %d pairs)\n', ...
                devField{1}, numDynamicPairs, numStaticPairs, numLeafPairs);
        end
    end
    
    %% ===== STEP 2: CORRECT DYNAMIC LEAVES WITH SMALL GAPS =====
    fprintf('  Step 2: Correcting dynamic leaves with small gaps...\n');
    
    % Debug: Show what's in the dynamic leaves map before we start
    fprintf('  DEBUG: Dynamic leaves map contains %d devices\n', dynamicLeavesMap.Count);
    if dynamicLeavesMap.Count > 0
        fprintf('  DEBUG: Dynamic map keys: ');
        for key = dynamicLeavesMap.keys
            fprintf('%s ', key{1});
        end
        fprintf('\n');
    end
    
    beamCorrections = 0;
    
    % Loop through each control point for corrections
    for cpIdx = 1:numControlPoints
        cpField = cpSeqFields{cpIdx};
        controlPoint = beam.ControlPointSequence.(cpField);
        cpIndex = controlPoint.ControlPointIndex;
        
        fprintf('  DEBUG: Processing CP %d (loop index %d)\n', cpIndex, cpIdx);
        
        totalControlPoints = totalControlPoints + 1;
        
        if ~isfield(controlPoint, 'BeamLimitingDevicePositionSequence')
            fprintf('  DEBUG: CP %d - No BeamLimitingDevicePositionSequence, skipping\n', cpIndex);
            continue;
        end
        
        bldSeqFields = fieldnames(controlPoint.BeamLimitingDevicePositionSequence);
        fprintf('  DEBUG: CP %d - Found %d beam limiting devices\n', cpIndex, length(bldSeqFields));
        
        % Loop through each beam limiting device
        for devIdx = 1:length(bldSeqFields)
            devField = bldSeqFields{devIdx};
            device = controlPoint.BeamLimitingDevicePositionSequence.(devField);
            deviceType = device.RTBeamLimitingDeviceType;
            
            if ~any(strcmpi(deviceType, MLC_TYPES))
                fprintf('  DEBUG: CP %d, Device %s - Not MLC (type: %s), skipping\n', cpIndex, devField, deviceType);
                continue;
            end
            
            fprintf('  DEBUG: CP %d, Device %s - Processing MLC device (type: %s)\n', cpIndex, devField, deviceType);
            
            % Get leaf positions
            leafPositions = device.LeafJawPositions;
            numLeaves = length(leafPositions);
            numLeafPairs = numLeaves / 2;
            
            if mod(numLeaves, 2) ~= 0
                warning('Odd number of leaf positions at Beam %d, CP %d - skipping', beamNumber, cpIndex);
                continue;
            end
            
            % Get dynamic leaf information
            if ~dynamicLeavesMap.isKey(devField)
                fprintf('  DEBUG: CP %d, Device %s - WARNING: No dynamic leaf info found!\n', cpIndex, devField);
                continue;
            end
            devInfo = dynamicLeavesMap(devField);
            
            % Split into two banks
            bankA = leafPositions(1:numLeafPairs);
            bankB = leafPositions(numLeafPairs+1:end);
            
            % Get dynamic flags for each bank
            bankA_dynamic = devInfo.isDynamic(1:numLeafPairs);
            bankB_dynamic = devInfo.isDynamic(numLeafPairs+1:end);
            leafPair_dynamic = bankA_dynamic | bankB_dynamic;
            
            numDynamicInCP = sum(leafPair_dynamic);
            fprintf('  DEBUG: CP %d, Device %s - %d dynamic leaf pairs (out of %d total)\n', ...
                cpIndex, devField, numDynamicInCP, numLeafPairs);
            
            % Calculate gaps for each leaf pair
            gaps = bankB - bankA;
            
            totalLeafPairsChecked = totalLeafPairsChecked + numLeafPairs;
            
            % Find DYNAMIC leaf pairs with gap less than minimum
            smallGapIdx = find(gaps < MIN_TIP_GAP & leafPair_dynamic);
            
            if ~isempty(smallGapIdx)
                fprintf('  DEBUG: CP %d, Device %s - Found %d dynamic leaves with small gaps\n', ...
                    cpIndex, devField, length(smallGapIdx));
                
                for leafIdx = smallGapIdx'
                    originalGap = gaps(leafIdx);
                    originalA = bankA(leafIdx);
                    originalB = bankB(leafIdx);
                    
                    fprintf('  DEBUG:   CP %d, Leaf %d - Gap %.4f mm < threshold %.4f mm (A=%.4f, B=%.4f)\n', ...
                        cpIndex, leafIdx, originalGap, MIN_TIP_GAP, originalA, originalB);
                    
                    % Expand tips: move A left (more negative) and B right (more positive)
                    newA = originalA - EXPANSION_PER_SIDE;
                    newB = originalB + EXPANSION_PER_SIDE;
                    
                    % Apply bounds checking - clamp to valid MLC range
                    if newA < MLC_MIN_POSITION
                        newA = MLC_MIN_POSITION;
                        fprintf('    WARNING: Beam %d, CP %d, Leaf %d: Bank A clamped to %.1f mm\n', ...
                            beamNumber, cpIndex, leafIdx, MLC_MIN_POSITION);
                    end
                    if newB > MLC_MAX_POSITION
                        newB = MLC_MAX_POSITION;
                        fprintf('    WARNING: Beam %d, CP %d, Leaf %d: Bank B clamped to %.1f mm\n', ...
                            beamNumber, cpIndex, leafIdx, MLC_MAX_POSITION);
                    end
                    
                    newGap = newB - newA;
                    
                    fprintf('  DEBUG:   CP %d, Leaf %d - Corrected: new gap %.4f mm (A: %.4f->%.4f, B: %.4f->%.4f)\n', ...
                        cpIndex, leafIdx, newGap, originalA, newA, originalB, newB);
                    
                    % Update the positions
                    bankA(leafIdx) = newA;
                    bankB(leafIdx) = newB;
                    
                    % Log the correction
                    correctionLog{end+1} = sprintf(...
                        'Beam %d (%s), CP %d, Leaf %d: Gap %.2f -> %.2f mm (A: %.2f->%.2f, B: %.2f->%.2f)', ...
                        beamNumber, beamName, cpIndex, leafIdx, ...
                        originalGap, newGap, originalA, newA, originalB, newB);
                    
                    beamCorrections = beamCorrections + 1;
                    totalCorrections = totalCorrections + 1;
                end
                
                % Reconstruct leaf positions array and update in structure
                newLeafPositions = [bankA; bankB];
                rtplan.BeamSequence.(beamField).ControlPointSequence.(cpField)...
                    .BeamLimitingDevicePositionSequence.(devField).LeafJawPositions = newLeafPositions;
                
                fprintf('  DEBUG: CP %d, Device %s - Updated leaf positions in RTPLAN structure\n', cpIndex, devField);
            else
                if numDynamicInCP > 0
                    fprintf('  DEBUG: CP %d, Device %s - No dynamic leaves with small gaps (all gaps OK)\n', ...
                        cpIndex, devField);
                else
                    fprintf('  DEBUG: CP %d, Device %s - No dynamic leaves to check\n', cpIndex, devField);
                end
            end
        end
    end
    
    fprintf('  Processed %d control points, made %d corrections\n', ...
        numControlPoints, beamCorrections);
end

%% ===================== VERIFICATION =====================
fprintf('\n=============================================================\n');
fprintf('VERIFICATION\n');
fprintf('=============================================================\n\n');

fprintf('Total control points processed: %d\n', totalControlPoints);
fprintf('Total leaf pairs checked: %d\n', totalLeafPairsChecked);
fprintf('Total corrections made: %d\n\n', totalCorrections);

if totalCorrections > 0
    fprintf('Correction details:\n');
    fprintf('-------------------\n');
    for i = 1:length(correctionLog)
        fprintf('%s\n', correctionLog{i});
    end
    fprintf('\n');
end

% Verify all gaps are now >= minimum
fprintf('Verifying all DYNAMIC leaf gaps are now >= %.2f mm...\n', MIN_TIP_GAP);

verificationPassed = true;
remainingSmallGaps = 0;

for beamIdx = 1:numBeams
    beamField = beamSeqFields{beamIdx};
    beam = rtplan.BeamSequence.(beamField);
    
    if ~isfield(beam, 'ControlPointSequence')
        continue;
    end
    
    cpSeqFields = fieldnames(beam.ControlPointSequence);
    
    % Sort control points by ControlPointIndex
    cpIndices = zeros(length(cpSeqFields), 1);
    for i = 1:length(cpSeqFields)
        cp = beam.ControlPointSequence.(cpSeqFields{i});
        if isfield(cp, 'ControlPointIndex')
            cpIndices(i) = cp.ControlPointIndex;
        else
            cpIndices(i) = i - 1;
        end
    end
    [~, sortOrder] = sort(cpIndices);
    cpSeqFields = cpSeqFields(sortOrder);
    
    % Rebuild dynamic leaves map using same improved logic as main processing
    allPositionsMap = containers.Map();
    
    % First pass: Collect all positions from all control points
    for cpIdx = 1:length(cpSeqFields)
        cpField = cpSeqFields{cpIdx};
        controlPoint = beam.ControlPointSequence.(cpField);
        
        if ~isfield(controlPoint, 'BeamLimitingDevicePositionSequence')
            continue;
        end
        
        bldSeqFields = fieldnames(controlPoint.BeamLimitingDevicePositionSequence);
        
        for devIdx = 1:length(bldSeqFields)
            devField = bldSeqFields{devIdx};
            device = controlPoint.BeamLimitingDevicePositionSequence.(devField);
            deviceType = device.RTBeamLimitingDeviceType;
            
            if ~any(strcmpi(deviceType, MLC_TYPES))
                continue;
            end
            
            leafPositions = device.LeafJawPositions;
            
            if ~allPositionsMap.isKey(devField)
                numLeaves = length(leafPositions);
                allPositionsMap(devField) = struct(...
                    'positions', leafPositions(:), ...
                    'numLeafPairs', numLeaves / 2);
            else
                temp = allPositionsMap(devField);
                temp.positions = [temp.positions, leafPositions(:)];
                allPositionsMap(devField) = temp;
            end
        end
    end
    
    % Second pass: Determine which leaves are dynamic
    dynamicLeavesMap = containers.Map();
    
    if ~isempty(allPositionsMap.keys)
        for devField = allPositionsMap.keys
            devInfo = allPositionsMap(devField{1});
            allPos = devInfo.positions;
            numLeafPairs = devInfo.numLeafPairs;
            numLeaves = size(allPos, 1);
            
            isDynamic = false(numLeaves, 1);
            for leafIdx = 1:numLeaves
                leafPosAcrossCPs = allPos(leafIdx, :);
                variation = max(leafPosAcrossCPs) - min(leafPosAcrossCPs);
                if variation > 1e-4
                    isDynamic(leafIdx) = true;
                end
            end
            
            dynamicLeavesMap(devField{1}) = struct(...
                'isDynamic', isDynamic, ...
                'numLeafPairs', numLeafPairs);
        end
    end
    
    % Third pass: verify gaps only for dynamic leaves
    for cpIdx = 1:length(cpSeqFields)
        cpField = cpSeqFields{cpIdx};
        controlPoint = beam.ControlPointSequence.(cpField);
        
        if ~isfield(controlPoint, 'BeamLimitingDevicePositionSequence')
            continue;
        end
        
        bldSeqFields = fieldnames(controlPoint.BeamLimitingDevicePositionSequence);
        
        for devIdx = 1:length(bldSeqFields)
            devField = bldSeqFields{devIdx};
            device = controlPoint.BeamLimitingDevicePositionSequence.(devField);
            deviceType = device.RTBeamLimitingDeviceType;
            
            if ~any(strcmpi(deviceType, MLC_TYPES))
                continue;
            end
            
            leafPositions = device.LeafJawPositions;
            numLeafPairs = length(leafPositions) / 2;
            
            bankA = leafPositions(1:numLeafPairs);
            bankB = leafPositions(numLeafPairs+1:end);
            gaps = bankB - bankA;
            
            % Get dynamic leaf information
            if dynamicLeavesMap.isKey(devField)
                devInfo = dynamicLeavesMap(devField);
                bankA_dynamic = devInfo.isDynamic(1:numLeafPairs);
                bankB_dynamic = devInfo.isDynamic(numLeafPairs+1:end);
                leafPair_dynamic = bankA_dynamic | bankB_dynamic;
                
                % Only check dynamic leaves
                smallGaps = (gaps < MIN_TIP_GAP) & leafPair_dynamic;
            else
                % If no dynamic info, check all (shouldn't happen)
                smallGaps = gaps < MIN_TIP_GAP;
            end
            
            if any(smallGaps)
                verificationPassed = false;
                remainingSmallGaps = remainingSmallGaps + sum(smallGaps);
                cpIndex = controlPoint.ControlPointIndex;
                fprintf('WARNING: Small gap still found at Beam %d, CP %d\n', ...
                    beam.BeamNumber, cpIndex);
            end
        end
    end
end

if verificationPassed
    fprintf('VERIFICATION PASSED: All dynamic MLC gaps are now >= %.2f mm\n\n', MIN_TIP_GAP);
else
    fprintf('VERIFICATION FAILED: %d small gaps remaining\n\n', remainingSmallGaps);
end

%% ===================== EXPORT CORRECTED PLAN =====================
fprintf('=============================================================\n');
fprintf('EXPORTING CORRECTED RTPLAN\n');
fprintf('=============================================================\n\n');

% Generate output filename
[~, inputName, inputExt] = fileparts(inputFile);
outputFile = fullfile(sctDir, [inputName, OUTPUT_SUFFIX, inputExt]);

fprintf('Output file: %s\n', outputFile);

% Generate new SOP Instance UID for the modified plan
% (Important: modified DICOM files should have new UIDs)
newSOPInstanceUID = dicomuid;
rtplan.SOPInstanceUID = newSOPInstanceUID;
rtplan.MediaStorageSOPInstanceUID = newSOPInstanceUID;

% Update plan description to indicate modification
if isfield(rtplan, 'RTPlanDescription')
    rtplan.RTPlanDescription = [rtplan.RTPlanDescription, ' - MLC Gap Adjusted'];
else
    rtplan.RTPlanDescription = 'MLC Gap Adjusted';
end

% Update instance creation date/time
rtplan.InstanceCreationDate = datestr(now, 'yyyymmdd');
rtplan.InstanceCreationTime = datestr(now, 'HHMMSS');

% Write the DICOM file
try
    dicomwrite([], outputFile, rtplan, 'CreateMode', 'copy');
    fprintf('Successfully exported corrected RTPLAN!\n');
catch ME
    % Alternative method using dicomwrite with different syntax
    try
        % Some MATLAB versions require different approach
        status = dicomwrite([], outputFile, rtplan);
        fprintf('Successfully exported corrected RTPLAN!\n');
    catch ME2
        error('Failed to write DICOM file: %s\nAlternative method also failed: %s', ...
            ME.message, ME2.message);
    end
end

%% ===================== SUMMARY =====================
fprintf('\n=============================================================\n');
fprintf('SUMMARY\n');
fprintf('=============================================================\n');
fprintf('Patient ID:            %s\n', patientID);
fprintf('Session:               %s\n', sessionName);
fprintf('Input file:            %s\n', inputFile);
fprintf('Output file:           %s\n', outputFile);
fprintf('Minimum gap threshold: %.2f mm\n', MIN_TIP_GAP);
fprintf('Expansion per side:    %.2f mm\n', EXPANSION_PER_SIDE);
fprintf('MLC position bounds:   [%.1f, %.1f] mm\n', MLC_MIN_POSITION, MLC_MAX_POSITION);
fprintf('Total corrections:     %d\n', totalCorrections);
fprintf('Verification:          %s\n', conditional(verificationPassed, 'PASSED', 'FAILED'));
fprintf('=============================================================\n\n');

    end  % End session loop
end  % End patient loop

fprintf('=============================================================\n');
fprintf('ALL PATIENTS AND SESSIONS PROCESSED\n');
fprintf('=============================================================\n');

%% ===================== HELPER FUNCTIONS =====================
function result = conditional(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end