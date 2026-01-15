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
% Input RTPLAN file path - modify this to your actual file path
inputFile = 'path/to/your/rtplan.dcm';  % <-- MODIFY THIS

% Minimum tip gap threshold (cm) - from commissioning
MIN_TIP_GAP = 0.06;  % cm

% Expansion amount per side when gap is too small (cm)
EXPANSION_PER_SIDE = 0.04;  % cm (total expansion will be 0.08 cm)

% Output filename suffix
OUTPUT_SUFFIX = '_adjusted_mlc';

%% ===================== READ DICOM RTPLAN =====================
fprintf('=============================================================\n');
fprintf('MLC Tip Gap Correction Script\n');
fprintf('=============================================================\n\n');

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
    
    % Loop through each control point
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
            
            % Split into two banks
            bankA = leafPositions(1:numLeafPairs);       % Left bank (usually negative)
            bankB = leafPositions(numLeafPairs+1:end);   % Right bank (usually positive)
            
            % Calculate gaps for each leaf pair
            % Gap = bankB - bankA (right position - left position)
            gaps = bankB - bankA;
            
            totalLeafPairsChecked = totalLeafPairsChecked + numLeafPairs;
            
            % Find leaf pairs with gap less than minimum
            smallGapIdx = find(gaps < MIN_TIP_GAP);
            
            if ~isempty(smallGapIdx)
                for leafIdx = smallGapIdx'
                    originalGap = gaps(leafIdx);
                    originalA = bankA(leafIdx);
                    originalB = bankB(leafIdx);
                    
                    % Expand tips: move A left (more negative) and B right (more positive)
                    newA = originalA - EXPANSION_PER_SIDE;
                    newB = originalB + EXPANSION_PER_SIDE;
                    newGap = newB - newA;
                    
                    % Update the positions
                    bankA(leafIdx) = newA;
                    bankB(leafIdx) = newB;
                    
                    % Log the correction
                    correctionLog{end+1} = sprintf(...
                        'Beam %d (%s), CP %d, Leaf %d: Gap %.4f -> %.4f cm (A: %.4f->%.4f, B: %.4f->%.4f)', ...
                        beamNumber, beamName, cpIdx, leafIdx, ...
                        originalGap, newGap, originalA, newA, originalB, newB);
                    
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
fprintf('=============================================================\n');

%% ===================== HELPER FUNCTIONS =====================
function result = conditional(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end
