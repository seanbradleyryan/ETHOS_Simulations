%% fix_mlc_tip_gap.m
% Script to fix minimum dynamic tip gap issues in RTPLAN DICOM files
% 
% This script:
% 1. Reads an RTPLAN DICOM file
% 2. Loops through all beams and control points
% 3. Identifies MLC leaf pairs with gaps less than the minimum threshold
% 4. Expands the tips by a specified amount on each side
% 5. Verifies the corrections
% 6. Exports the corrected plan as a new DICOM file
%
% Author: Claude
% Date: January 2026

clear; clc;

%% ===================== USER PARAMETERS =====================
% Base data directory - modify this to your actual data path
BASE_DATA_DIR = 'path/to/data';  % <-- MODIFY THIS

% Patient IDs to process
PATIENT_IDS = {'1194203'};  % Add more patient IDs as needed

% Sessions to process
SESSIONS = {'Session_1'};  % Add more sessions as needed

% Minimum tip gap threshold (cm) - from commissioning
MIN_TIP_GAP = 0.06;  % cm

% Expansion amount per side when gap is too small (cm)
EXPANSION_PER_SIDE = 0.04;  % cm (total expansion will be 0.08 cm)

% MLC coordinate bounds (cm) - Varian standard
MLC_MIN_COORD = -140;  % cm
MLC_MAX_COORD = 140;   % cm

% Output filename suffix
OUTPUT_SUFFIX = '_adjusted_mlc';


%% ===================== MAIN PROCESSING LOOP =====================
fprintf('=============================================================\n');
fprintf('MLC Tip Gap Correction Script\n');
fprintf('=============================================================\n\n');
fprintf('Base data directory: %s\n', BASE_DATA_DIR);
fprintf('Patient IDs: %s\n', strjoin(PATIENT_IDS, ', '));
fprintf('Sessions: %s\n', strjoin(SESSIONS, ', '));
fprintf('Minimum tip gap threshold: %.4f cm\n', MIN_TIP_GAP);
fprintf('Expansion per side: %.4f cm\n\n', EXPANSION_PER_SIDE);

% Loop over patient IDs
for patIdx = 1:length(PATIENT_IDS)
    patientID = PATIENT_IDS{patIdx};
    
    % Loop over sessions
    for sessIdx = 1:length(SESSIONS)
        sessionName = SESSIONS{sessIdx};
        
        fprintf('\n#############################################################\n');
        fprintf('Processing Patient: %s, Session: %s\n', patientID, sessionName);
        fprintf('#############################################################\n\n');
        
        % Build path to sct directory
        sctDir = fullfile(BASE_DATA_DIR, patientID, sessionName, 'sct');
        
        % Check if directory exists
        if ~exist(sctDir, 'dir')
            fprintf('WARNING: Directory does not exist: %s\n', sctDir);
            fprintf('Skipping this patient/session combination.\n');
            continue;
        end
        
        % Search for RTPLAN files
        fprintf('Scanning directory: %s\n', sctDir);
        rtplanFiles = dir(fullfile(sctDir, 'RP*.dcm'));
        
        % Also try lowercase pattern
        if isempty(rtplanFiles)
            rtplanFiles = dir(fullfile(sctDir, 'rp*.dcm'));
        end
        
        % Also try general DICOM search and filter by modality
        if isempty(rtplanFiles)
            allDcmFiles = dir(fullfile(sctDir, '*.dcm'));
            rtplanFiles = [];
            for i = 1:length(allDcmFiles)
                try
                    tempInfo = dicominfo(fullfile(sctDir, allDcmFiles(i).name));
                    if isfield(tempInfo, 'Modality') && strcmpi(tempInfo.Modality, 'RTPLAN')
                        rtplanFiles = [rtplanFiles; allDcmFiles(i)];
                    end
                catch
                    % Skip files that can't be read
                end
            end
        end
        
        if isempty(rtplanFiles)
            fprintf('WARNING: No RTPLAN files found in %s\n', sctDir);
            fprintf('Skipping this patient/session combination.\n');
            continue;
        end
        
        fprintf('Found %d RTPLAN file(s)\n\n', length(rtplanFiles));
        
        % Process each RTPLAN file
        for fileIdx = 1:length(rtplanFiles)
            inputFile = fullfile(sctDir, rtplanFiles(fileIdx).name);
            
            fprintf('-------------------------------------------------------------\n');
            fprintf('Processing RTPLAN file %d of %d: %s\n', fileIdx, length(rtplanFiles), rtplanFiles(fileIdx).name);
            fprintf('-------------------------------------------------------------\n\n');
            
            % Call the processing function
            try
                processRTPlan(inputFile, MIN_TIP_GAP, EXPANSION_PER_SIDE, ...
                    MLC_MIN_COORD, MLC_MAX_COORD, OUTPUT_SUFFIX);
            catch ME
                fprintf('ERROR processing file %s: %s\n', inputFile, ME.message);
                fprintf('Continuing to next file...\n\n');
            end
        end
    end
end

fprintf('\n=============================================================\n');
fprintf('ALL PROCESSING COMPLETE\n');
fprintf('=============================================================\n');

%% ===================== PROCESSING FUNCTION =====================
function processRTPlan(inputFile, MIN_TIP_GAP, EXPANSION_PER_SIDE, ...
    MLC_MIN_COORD, MLC_MAX_COORD, OUTPUT_SUFFIX)
% Process a single RTPLAN file for MLC tip gap corrections

fprintf('Reading RTPLAN file: %s\n', inputFile);

% Check if file exists
if ~exist(inputFile, 'file')
    error('Input file does not exist: %s', inputFile);
end

% Read the DICOM file
try
    rtplan = dicominfo(inputFile);
    fprintf('Successfully loaded RTPLAN: %s\n', rtplan.RTPlanLabel);
catch ME
    error('Failed to read DICOM file: %s', ME.message);
end

%% ===================== IDENTIFY MLC DEVICE TYPE =====================
% MLC device types in DICOM
MLC_TYPES = {'MLCX', 'MLCY', 'MLCX1', 'MLCX2', 'MLCY1', 'MLCY2'};

%% ===================== PROCESS EACH BEAM =====================
% Get number of beams
beamSeqFields = fieldnames(rtplan.BeamSequence);
numBeams = length(beamSeqFields);

fprintf('\nProcessing %d beams...\n', numBeams);
fprintf('Minimum tip gap threshold: %.4f cm\n', MIN_TIP_GAP);
fprintf('Expansion per side: %.4f cm\n\n', EXPANSION_PER_SIDE);

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
    
    beamCorrections = 0;
    
    %% Identify dynamic vs static leaves for this beam
    % First, collect all MLC positions across all control points
    mlcDeviceData = struct();  % Will store data for each MLC device found
    
    fprintf('  Identifying dynamic vs static leaves...\n');
    
    for cpIdx = 1:numControlPoints
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
            
            % Only process MLC devices
            if ~any(strcmpi(deviceType, MLC_TYPES))
                continue;
            end
            
            % Create a unique key for this device type
            devKey = sprintf('dev_%s', deviceType);
            
            % Initialize storage for this device if first time seeing it
            if ~isfield(mlcDeviceData, devKey)
                mlcDeviceData.(devKey).deviceType = deviceType;
                mlcDeviceData.(devKey).positions = [];  % Will store all positions
                mlcDeviceData.(devKey).numLeafPairs = [];
            end
            
            % Store this control point's positions
            leafPositions = device.LeafJawPositions;
            mlcDeviceData.(devKey).positions = [mlcDeviceData.(devKey).positions; leafPositions'];
            mlcDeviceData.(devKey).numLeafPairs = length(leafPositions) / 2;
        end
    end
    
    % Determine which leaves are dynamic (change position) vs static
    devKeys = fieldnames(mlcDeviceData);
    for devKeyIdx = 1:length(devKeys)
        devKey = devKeys{devKeyIdx};
        positions = mlcDeviceData.(devKey).positions;  % Each row is a control point
        numLeafPairs = mlcDeviceData.(devKey).numLeafPairs;
        
        % Check if each leaf position changes across control points
        isDynamic = false(numLeafPairs * 2, 1);  % For both banks combined
        
        for leafIdx = 1:(numLeafPairs * 2)
            % Check if this leaf position varies across control points
            leafPositionsOverTime = positions(:, leafIdx);
            if max(leafPositionsOverTime) - min(leafPositionsOverTime) > 1e-6  % Tolerance for floating point
                isDynamic(leafIdx) = true;
            end
        end
        
        % Store the dynamic status for each leaf pair (both banks)
        mlcDeviceData.(devKey).isDynamic = isDynamic;
        
        % Count dynamic leaf pairs
        numDynamicPairs = sum(isDynamic(1:numLeafPairs) | isDynamic(numLeafPairs+1:end));
        fprintf('    Device %s: %d/%d leaf pairs are dynamic\n', ...
            mlcDeviceData.(devKey).deviceType, numDynamicPairs, numLeafPairs);
    end
    
    % Loop through each control point to make corrections
    for cpIdx = 1:numControlPoints
        cpField = cpSeqFields{cpIdx};
        controlPoint = beam.ControlPointSequence.(cpField);
        
        totalControlPoints = totalControlPoints + 1;
        
        % Check if BeamLimitingDevicePositionSequence exists
        if ~isfield(controlPoint, 'BeamLimitingDevicePositionSequence')
            % Some control points may not have this field (inherit from previous)
            continue;
        end
        
        bldSeqFields = fieldnames(controlPoint.BeamLimitingDevicePositionSequence);
        numDevices = length(bldSeqFields);
        
        % Loop through each beam limiting device
        for devIdx = 1:numDevices
            devField = bldSeqFields{devIdx};
            device = controlPoint.BeamLimitingDevicePositionSequence.(devField);
            
            % Check if this is an MLC device
            deviceType = device.RTBeamLimitingDeviceType;
            
            if ~any(strcmpi(deviceType, MLC_TYPES))
                % Not an MLC device (probably jaw) - skip
                continue;
            end
            
            % Get leaf positions
            leafPositions = device.LeafJawPositions;
            
            % LeafJawPositions contains positions for all leaves
            % First half: Bank A (or left bank)
            % Second half: Bank B (or right bank)
            numLeaves = length(leafPositions);
            numLeafPairs = numLeaves / 2;
            
            if mod(numLeaves, 2) ~= 0
                warning('Odd number of leaf positions at Beam %d, CP %d - skipping', ...
                    beamNumber, cpIdx);
                continue;
            end
            
            % Get dynamic status for this device type
            devKey = sprintf('dev_%s', deviceType);
            if ~isfield(mlcDeviceData, devKey)
                warning('No dynamic data found for device %s - skipping', deviceType);
                continue;
            end
            isDynamic = mlcDeviceData.(devKey).isDynamic;
            
            % Split into two banks
            bankA = leafPositions(1:numLeafPairs);       % Left bank (usually negative)
            bankB = leafPositions(numLeafPairs+1:end);   % Right bank (usually positive)
            
            % Calculate gaps for each leaf pair
            % Gap = bankB - bankA (right position - left position)
            gaps = bankB - bankA;
            
            totalLeafPairsChecked = totalLeafPairsChecked + numLeafPairs;
            
            % Find leaf pairs with gap less than minimum AND are dynamic
            smallGapIdx = find(gaps < MIN_TIP_GAP);
            
            if ~isempty(smallGapIdx)
                for leafIdx = smallGapIdx'
                    % Check if EITHER bank A or bank B for this leaf pair is dynamic
                    % If both are static, skip this leaf pair
                    if ~isDynamic(leafIdx) && ~isDynamic(leafIdx + numLeafPairs)
                        continue;  % Both banks static - skip
                    end
                    
                    originalGap = gaps(leafIdx);
                    originalA = bankA(leafIdx);
                    originalB = bankB(leafIdx);
                    
                    % Determine expansion strategy based on which banks are dynamic
                    bankAIsDynamic = isDynamic(leafIdx);
                    bankBIsDynamic = isDynamic(leafIdx + numLeafPairs);
                    
                    % Start with symmetric expansion
                    newA = originalA - EXPANSION_PER_SIDE;
                    newB = originalB + EXPANSION_PER_SIDE;
                    
                    % Apply boundary constraints
                    newA = max(newA, MLC_MIN_COORD);
                    newB = min(newB, MLC_MAX_COORD);
                    
                    % Check if after clamping we still have minimum gap
                    newGap = newB - newA;
                    if newGap < MIN_TIP_GAP
                        % Try asymmetric expansion to maintain minimum gap
                        if bankAIsDynamic && bankBIsDynamic
                            % Both dynamic - try to expand more on the side with room
                            totalExpansionNeeded = MIN_TIP_GAP - originalGap;
                            
                            % Calculate available room on each side
                            roomA = originalA - MLC_MIN_COORD;
                            roomB = MLC_MAX_COORD - originalB;
                            
                            if roomA + roomB >= totalExpansionNeeded
                                % Distribute expansion based on available room
                                if roomA >= totalExpansionNeeded / 2 && roomB >= totalExpansionNeeded / 2
                                    % Enough room on both sides for symmetric
                                    newA = originalA - totalExpansionNeeded / 2;
                                    newB = originalB + totalExpansionNeeded / 2;
                                elseif roomA >= roomB
                                    % More room on A side
                                    newB = min(originalB + roomB, MLC_MAX_COORD);
                                    newA = newB - MIN_TIP_GAP;
                                else
                                    % More room on B side
                                    newA = max(originalA - roomA, MLC_MIN_COORD);
                                    newB = newA + MIN_TIP_GAP;
                                end
                            else
                                % Not enough room even with asymmetric expansion
                                warning('Cannot achieve MIN_TIP_GAP for Beam %d, CP %d, Leaf %d due to bounds', ...
                                    beamNumber, cpIdx, leafIdx);
                                continue;
                            end
                        elseif bankAIsDynamic
                            % Only A is dynamic - expand only A
                            newA = originalB - MIN_TIP_GAP;
                            newB = originalB;  % Keep B static
                            if newA < MLC_MIN_COORD
                                warning('Cannot achieve MIN_TIP_GAP for Beam %d, CP %d, Leaf %d (A at boundary)', ...
                                    beamNumber, cpIdx, leafIdx);
                                continue;
                            end
                        elseif bankBIsDynamic
                            % Only B is dynamic - expand only B
                            newA = originalA;  % Keep A static
                            newB = originalA + MIN_TIP_GAP;
                            if newB > MLC_MAX_COORD
                                warning('Cannot achieve MIN_TIP_GAP for Beam %d, CP %d, Leaf %d (B at boundary)', ...
                                    beamNumber, cpIdx, leafIdx);
                                continue;
                            end
                        end
                        
                        newGap = newB - newA;
                    end
                    
                    % Only update dynamic banks
                    if bankAIsDynamic
                        bankA(leafIdx) = newA;
                    end
                    if bankBIsDynamic
                        bankB(leafIdx) = newB;
                    end
                    
                    % Log the correction
                    correctionLog{end+1} = sprintf(...
                        'Beam %d (%s), CP %d, Leaf %d: Gap %.4f -> %.4f cm (A: %.4f->%.4f%s, B: %.4f->%.4f%s)', ...
                        beamNumber, beamName, cpIdx, leafIdx, ...
                        originalGap, newGap, originalA, bankA(leafIdx), ...
                        conditional(bankAIsDynamic, '', '[static]'), ...
                        originalB, bankB(leafIdx), ...
                        conditional(bankBIsDynamic, '', '[static]'));
                    
                    beamCorrections = beamCorrections + 1;
                    totalCorrections = totalCorrections + 1;
                end
                
                % Reconstruct leaf positions array and update in structure
                newLeafPositions = [bankA; bankB];
                rtplan.BeamSequence.(beamField).ControlPointSequence.(cpField)...
                    .BeamLimitingDevicePositionSequence.(devField).LeafJawPositions = newLeafPositions;
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
fprintf('Verifying all gaps are now >= %.4f cm...\n', MIN_TIP_GAP);

verificationPassed = true;
remainingSmallGaps = 0;

for beamIdx = 1:numBeams
    beamField = beamSeqFields{beamIdx};
    beam = rtplan.BeamSequence.(beamField);
    
    if ~isfield(beam, 'ControlPointSequence')
        continue;
    end
    
    cpSeqFields = fieldnames(beam.ControlPointSequence);
    
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
            
            smallGaps = gaps < MIN_TIP_GAP;
            if any(smallGaps)
                verificationPassed = false;
                remainingSmallGaps = remainingSmallGaps + sum(smallGaps);
                fprintf('WARNING: Small gap still found at Beam %d, CP %d\n', ...
                    beam.BeamNumber, cpIdx);
            end
        end
    end
end

if verificationPassed
    fprintf('VERIFICATION PASSED: All MLC gaps are now >= %.4f cm\n\n', MIN_TIP_GAP);
else
    fprintf('VERIFICATION FAILED: %d small gaps remaining\n\n', remainingSmallGaps);
end

%% ===================== EXPORT CORRECTED PLAN =====================
fprintf('=============================================================\n');
fprintf('EXPORTING CORRECTED RTPLAN\n');
fprintf('=============================================================\n\n');

% Generate output filename
[inputPath, inputName, inputExt] = fileparts(inputFile);
outputFile = fullfile(inputPath, [inputName, OUTPUT_SUFFIX, inputExt]);

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
fprintf('Input file:  %s\n', inputFile);
fprintf('Output file: %s\n', outputFile);
fprintf('Minimum gap threshold: %.4f cm\n', MIN_TIP_GAP);
fprintf('Expansion per side:    %.4f cm\n', EXPANSION_PER_SIDE);
fprintf('Total corrections:     %d\n', totalCorrections);
fprintf('Verification:          %s\n', conditional(verificationPassed, 'PASSED', 'FAILED'));
fprintf('=============================================================\n\n');

end  % End of processRTPlan function

%% ===================== HELPER FUNCTIONS =====================
function result = conditional(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end