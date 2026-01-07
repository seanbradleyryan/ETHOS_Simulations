%% ETHOS Segment Dose Diagnostic Script
% Purpose: Verify that MLC data is being correctly extracted and applied
%          This helps debug dose calculation discrepancies
%
% Author: Generated for ETHOS dose analysis debugging
% Date: 2025

clear; clc; close all;

%% Configuration
patientID = '1194203';
sessionName = 'Session_1';
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
dicomPath = fullfile(wd, 'EthosExports', patientID, 'Pancreas', sessionName, 'sct');

fprintf('==========================================================\n');
fprintf('  ETHOS Segment Dose Diagnostics\n');
fprintf('  Patient: %s, Session: %s\n', patientID, sessionName);
fprintf('==========================================================\n\n');

%% Load RTPLAN and analyze structure
fprintf('[1] Loading RTPLAN...\n');

rtplanFile = dir(fullfile(dicomPath, 'RP*.dcm'));
if isempty(rtplanFile)
    rtplanFile = dir(fullfile(dicomPath, '*RTPLAN*.dcm'));
end

if isempty(rtplanFile)
    error('No RTPLAN file found');
end

rtplanInfo = dicominfo(fullfile(rtplanFile(1).folder, rtplanFile(1).name));
fprintf('  RTPLAN loaded: %s\n', rtplanFile(1).name);

%% Analyze beam limiting device sequence
fprintf('\n[2] Analyzing Beam Limiting Devices...\n');

numBeams = length(fieldnames(rtplanInfo.BeamSequence));
fprintf('  Found %d beams\n', numBeams);

for beamIdx = 1:min(3, numBeams)  % Just analyze first 3 beams
    beamField = sprintf('Item_%d', beamIdx);
    beam = rtplanInfo.BeamSequence.(beamField);
    
    beamName = '';
    if isfield(beam, 'BeamName')
        beamName = beam.BeamName;
    end
    
    fprintf('\n  Beam %d (%s):\n', beamIdx, beamName);
    
    % Analyze beam limiting device sequence
    if isfield(beam, 'BeamLimitingDeviceSequence')
        numDevices = length(fieldnames(beam.BeamLimitingDeviceSequence));
        fprintf('    Beam Limiting Devices: %d\n', numDevices);
        
        for devIdx = 1:numDevices
            device = beam.BeamLimitingDeviceSequence.(sprintf('Item_%d', devIdx));
            if isfield(device, 'RTBeamLimitingDeviceType')
                devType = device.RTBeamLimitingDeviceType;
                fprintf('      Device %d: %s', devIdx, devType);
                
                if isfield(device, 'NumberOfLeafJawPairs')
                    fprintf(', Pairs=%d', device.NumberOfLeafJawPairs);
                end
                if isfield(device, 'LeafPositionBoundaries')
                    fprintf(', Boundaries=%d', length(device.LeafPositionBoundaries));
                end
                fprintf('\n');
            end
        end
    end
    
    % Analyze control point sequence
    if isfield(beam, 'ControlPointSequence')
        numCP = length(fieldnames(beam.ControlPointSequence));
        fprintf('    Control Points: %d\n', numCP);
        
        % Check which CPs have device position data
        cpsWithDevicePos = 0;
        mlcDevicesPerCP = [];
        
        for cpIdx = 1:numCP
            cp = beam.ControlPointSequence.(sprintf('Item_%d', cpIdx));
            if isfield(cp, 'BeamLimitingDevicePositionSequence')
                cpsWithDevicePos = cpsWithDevicePos + 1;
                numDevs = length(fieldnames(cp.BeamLimitingDevicePositionSequence));
                mlcDevicesPerCP(end+1) = numDevs;
                
                % Detailed analysis for first few CPs
                if cpIdx <= 3
                    fprintf('      CP %d: Has %d device positions\n', cpIdx, numDevs);
                    for devIdx = 1:numDevs
                        device = cp.BeamLimitingDevicePositionSequence.(sprintf('Item_%d', devIdx));
                        if isfield(device, 'RTBeamLimitingDeviceType')
                            devType = device.RTBeamLimitingDeviceType;
                            if isfield(device, 'LeafJawPositions')
                                nPos = length(device.LeafJawPositions);
                                positions = device.LeafJawPositions;
                                fprintf('        %s: %d positions [%.1f to %.1f]\n', ...
                                    devType, nPos, min(positions), max(positions));
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('    CPs with explicit device positions: %d / %d (%.1f%%)\n', ...
            cpsWithDevicePos, numCP, 100*cpsWithDevicePos/numCP);
        
        if ~isempty(mlcDevicesPerCP)
            fprintf('    Devices per CP: min=%d, max=%d, mean=%.1f\n', ...
                min(mlcDevicesPerCP), max(mlcDevicesPerCP), mean(mlcDevicesPerCP));
        end
    end
end

%% Analyze cumulative meterset weights
fprintf('\n[3] Analyzing Cumulative Meterset Weights...\n');

beam = rtplanInfo.BeamSequence.Item_1;
numCP = length(fieldnames(beam.ControlPointSequence));

weights = zeros(numCP, 1);
for cpIdx = 1:numCP
    cp = beam.ControlPointSequence.(sprintf('Item_%d', cpIdx));
    if isfield(cp, 'CumulativeMetersetWeight')
        weights(cpIdx) = cp.CumulativeMetersetWeight;
    end
end

segmentWeights = diff(weights);
fprintf('  Beam 1 Cumulative Meterset Weights:\n');
fprintf('    Range: [%.6f, %.6f]\n', min(weights), max(weights));
fprintf('    Final weight: %.6f (should be 1.0)\n', weights(end));
fprintf('  Segment weights:\n');
fprintf('    Range: [%.6f, %.6f]\n', min(segmentWeights), max(segmentWeights));
fprintf('    Sum: %.6f (should be ~1.0)\n', sum(segmentWeights));
fprintf('    Mean: %.6f\n', mean(segmentWeights));
fprintf('    Segments with weight > 0: %d / %d\n', sum(segmentWeights > 0), length(segmentWeights));

%% Analyze MLC position patterns
fprintf('\n[4] Analyzing MLC Position Patterns...\n');

% Collect MLC positions for first beam
beam = rtplanInfo.BeamSequence.Item_1;
numCP = length(fieldnames(beam.ControlPointSequence));

mlcDataPresent = false(numCP, 1);
mlcPositionsFirstLayer = cell(numCP, 1);
mlcPositionsSecondLayer = cell(numCP, 1);

for cpIdx = 1:numCP
    cp = beam.ControlPointSequence.(sprintf('Item_%d', cpIdx));
    if isfield(cp, 'BeamLimitingDevicePositionSequence')
        numDevs = length(fieldnames(cp.BeamLimitingDevicePositionSequence));
        mlcCount = 0;
        
        for devIdx = 1:numDevs
            device = cp.BeamLimitingDevicePositionSequence.(sprintf('Item_%d', devIdx));
            if isfield(device, 'RTBeamLimitingDeviceType') && ...
               contains(device.RTBeamLimitingDeviceType, 'MLC', 'IgnoreCase', true)
                mlcCount = mlcCount + 1;
                if mlcCount == 1
                    mlcPositionsFirstLayer{cpIdx} = device.LeafJawPositions;
                elseif mlcCount == 2
                    mlcPositionsSecondLayer{cpIdx} = device.LeafJawPositions;
                end
            end
        end
        
        if mlcCount > 0
            mlcDataPresent(cpIdx) = true;
        end
    end
end

fprintf('  Control points with explicit MLC data: %d / %d\n', sum(mlcDataPresent), numCP);

% Find patterns of consecutive CPs without MLC data
consecutiveWithout = 0;
maxConsecutive = 0;
for cpIdx = 1:numCP
    if ~mlcDataPresent(cpIdx)
        consecutiveWithout = consecutiveWithout + 1;
        maxConsecutive = max(maxConsecutive, consecutiveWithout);
    else
        consecutiveWithout = 0;
    end
end
fprintf('  Max consecutive CPs without MLC data: %d\n', maxConsecutive);

% Check if MLC positions change between CPs that have data
cpsWithData = find(mlcDataPresent);
if length(cpsWithData) >= 2
    changeCount = 0;
    for i = 2:length(cpsWithData)
        cp1 = cpsWithData(i-1);
        cp2 = cpsWithData(i);
        
        if ~isempty(mlcPositionsFirstLayer{cp1}) && ~isempty(mlcPositionsFirstLayer{cp2})
            if ~isequal(mlcPositionsFirstLayer{cp1}, mlcPositionsFirstLayer{cp2})
                changeCount = changeCount + 1;
            end
        end
    end
    fprintf('  MLC position changes between data-present CPs: %d\n', changeCount);
end

%% Analyze aperture sizes
fprintf('\n[5] Analyzing MLC Aperture Sizes...\n');

% Get leaf boundaries
leafBoundaries = [];
if isfield(beam, 'BeamLimitingDeviceSequence')
    numDevices = length(fieldnames(beam.BeamLimitingDeviceSequence));
    for devIdx = 1:numDevices
        device = beam.BeamLimitingDeviceSequence.(sprintf('Item_%d', devIdx));
        if isfield(device, 'RTBeamLimitingDeviceType') && ...
           contains(device.RTBeamLimitingDeviceType, 'MLC', 'IgnoreCase', true)
            if isfield(device, 'LeafPositionBoundaries')
                leafBoundaries = device.LeafPositionBoundaries;
                break;
            end
        end
    end
end

if ~isempty(leafBoundaries) && ~isempty(cpsWithData)
    fprintf('  Leaf boundaries: %d values (for %d leaf pairs)\n', ...
        length(leafBoundaries), length(leafBoundaries)-1);
    
    % Calculate aperture opening for CPs with data
    apertureAreas = zeros(length(cpsWithData), 1);
    
    for i = 1:length(cpsWithData)
        cpIdx = cpsWithData(i);
        mlc1 = mlcPositionsFirstLayer{cpIdx};
        mlc2 = mlcPositionsSecondLayer{cpIdx};
        
        if ~isempty(mlc1)
            % Use first layer if only one present
            mlcPos = mlc1;
            numLeafPairs = length(mlcPos) / 2;
            leftBank = mlcPos(1:numLeafPairs);
            rightBank = mlcPos(numLeafPairs+1:end);
            
            % Reduce dual-layer if present
            if ~isempty(mlc2)
                leftBank2 = mlc2(1:numLeafPairs);
                rightBank2 = mlc2(numLeafPairs+1:end);
                leftBank = max(leftBank, leftBank2);  % More closed
                rightBank = min(rightBank, rightBank2);  % More closed
            end
            
            % Calculate aperture area
            totalArea = 0;
            for leafIdx = 1:numLeafPairs
                if leafIdx <= length(leafBoundaries) - 1
                    leafWidth = leafBoundaries(leafIdx+1) - leafBoundaries(leafIdx);
                    opening = max(0, rightBank(leafIdx) - leftBank(leafIdx));
                    totalArea = totalArea + opening * leafWidth;
                end
            end
            apertureAreas(i) = totalArea;
        end
    end
    
    fprintf('  Aperture areas (mm²) for CPs with data:\n');
    fprintf('    Min: %.1f, Max: %.1f, Mean: %.1f\n', ...
        min(apertureAreas), max(apertureAreas), mean(apertureAreas));
    
    % Compare to full field area
    if length(leafBoundaries) > 1
        fullFieldY = leafBoundaries(end) - leafBoundaries(1);
        mlcPos = mlcPositionsFirstLayer{cpsWithData(1)};
        if ~isempty(mlcPos)
            numLeafPairs = length(mlcPos) / 2;
            leftBank = mlcPos(1:numLeafPairs);
            rightBank = mlcPos(numLeafPairs+1:end);
            fullFieldX = max(rightBank) - min(leftBank);
            fullFieldArea = fullFieldX * fullFieldY;
            fprintf('    Full field area (approx): %.1f mm²\n', fullFieldArea);
            fprintf('    Average aperture as %% of full field: %.1f%%\n', ...
                100 * mean(apertureAreas) / fullFieldArea);
        end
    end
end

%% Summary and recommendations
fprintf('\n==========================================================\n');
fprintf('  DIAGNOSTIC SUMMARY\n');
fprintf('==========================================================\n');

% Check for dual-layer MLC
hasDualLayer = false;
beam = rtplanInfo.BeamSequence.Item_1;
if isfield(beam, 'BeamLimitingDeviceSequence')
    mlcCount = 0;
    numDevices = length(fieldnames(beam.BeamLimitingDeviceSequence));
    for devIdx = 1:numDevices
        device = beam.BeamLimitingDeviceSequence.(sprintf('Item_%d', devIdx));
        if isfield(device, 'RTBeamLimitingDeviceType') && ...
           contains(device.RTBeamLimitingDeviceType, 'MLC', 'IgnoreCase', true)
            mlcCount = mlcCount + 1;
        end
    end
    if mlcCount > 1
        hasDualLayer = true;
    end
end

fprintf('\n');
if hasDualLayer
    fprintf('  [!] DUAL-LAYER MLC DETECTED\n');
    fprintf('      The fixed script includes dual-layer reduction.\n');
end

needsInheritance = sum(~mlcDataPresent) > 0;
if needsInheritance
    fprintf('  [!] MLC INHERITANCE REQUIRED\n');
    fprintf('      %d CPs have no explicit MLC data and inherit from previous CP.\n', ...
        sum(~mlcDataPresent));
    fprintf('      The fixed script implements MLC position inheritance.\n');
end

fprintf('\n  EXPECTED BEHAVIOR AFTER FIX:\n');
fprintf('    - Each segment should have unique MLC positions (inherited or explicit)\n');
fprintf('    - Aperture area should vary significantly between segments\n');
fprintf('    - Calculated dose should be ~7 Gy max (matching reference)\n');
fprintf('\n');

% Plot aperture area distribution
if exist('apertureAreas', 'var') && ~isempty(apertureAreas)
    figure('Position', [100, 100, 800, 400]);
    subplot(1, 2, 1);
    histogram(apertureAreas, 30);
    xlabel('Aperture Area (mm²)');
    ylabel('Count');
    title('Distribution of MLC Aperture Areas');
    grid on;
    
    subplot(1, 2, 2);
    plot(1:length(apertureAreas), apertureAreas, '-o');
    xlabel('Control Point (with data)');
    ylabel('Aperture Area (mm²)');
    title('Aperture Area vs Control Point');
    grid on;
    
    saveas(gcf, fullfile(wd, 'Figures', 'mlc_diagnostics.png'));
    fprintf('  Saved aperture analysis figure to Figures/mlc_diagnostics.png\n');
end

fprintf('\n==========================================================\n');
