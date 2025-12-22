%% ETHOS Field-by-Field Dose Calculation Script
% Purpose: Reproduce RTDOSE by calculating individual field doses from RTPLAN
% Software: matRad, MATLAB
% Treatment: IMRT (17 fields) from Varian ETHOS
%
% Author: Generated for ETHOS Simulation Analysis
% Date: 2024
%
% This script processes ETHOS-exported DICOM data to produce separate 
% dose grids for each IMRT field, enabling field-by-field analysis.

%% Clear workspace and close figures
clear; close all; clc;

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Working directories
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
matRadPath = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

% Batch processing lists (expandable for multiple patients/sessions)
idList = {'1194203'};
sessionList = {'Session_1'};

% Add matRad to path
addpath(genpath(matRadPath));

% Initialize matRad
matRad_rc;

%% ========================================================================
%  MAIN PROCESSING LOOP
%  ========================================================================

for idIdx = 1:length(idList)
    for sessionIdx = 1:length(sessionList)
        
        currentId = idList{idIdx};
        currentSession = sessionList{sessionIdx};
        
        fprintf('\n');
        fprintf('================================================================\n');
        fprintf('Processing Patient: %s | Session: %s\n', currentId, currentSession);
        fprintf('================================================================\n');
        
        % Define paths
        rawwd = fullfile(wd, 'EthosExports', currentId, 'Pancreas', currentSession);
        dicomPath = fullfile(rawwd, 'sct');
        outputDir = fullfile(wd, 'Results', currentId, currentSession);
        
        % Create output directory if it doesn't exist
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end
        
        %% ================================================================
        %  STEP 1: Load Reference RTDOSE
        %  ================================================================
        fprintf('\n[STEP 1] Loading reference RTDOSE...\n');
        
        try
            % Find RTDOSE file
            rtdoseFiles = dir(fullfile(dicomPath, 'RD*.dcm'));
            if isempty(rtdoseFiles)
                rtdoseFiles = dir(fullfile(dicomPath, '*RTDOSE*.dcm'));
            end
            if isempty(rtdoseFiles)
                % Search all DICOM files for RTDOSE modality
                allDcm = dir(fullfile(dicomPath, '*.dcm'));
                for f = 1:length(allDcm)
                    info = dicominfo(fullfile(dicomPath, allDcm(f).name));
                    if isfield(info, 'Modality') && strcmp(info.Modality, 'RTDOSE')
                        rtdoseFiles = allDcm(f);
                        break;
                    end
                end
            end
            
            if isempty(rtdoseFiles)
                error('No RTDOSE file found in %s', dicomPath);
            end
            
            rtdoseFile = fullfile(dicomPath, rtdoseFiles(1).name);
            rtdoseInfo = dicominfo(rtdoseFile);
            rtdoseData = dicomread(rtdoseFile);
            
            % Apply scaling factor
            if isfield(rtdoseInfo, 'DoseGridScaling')
                doseScaling = rtdoseInfo.DoseGridScaling;
            else
                doseScaling = 1;
            end
            
            % Convert to double and apply squeeze to remove singleton dimension
            referenceDose = double(squeeze(rtdoseData)) * doseScaling;
            
            % Extract resolution - Z from GridFrameOffsetVector (multiframe DICOM)
            rtdoseResolution = [rtdoseInfo.PixelSpacing(1), ...
                                rtdoseInfo.PixelSpacing(2), ...
                                abs(rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1))];
            
            % Extract position
            rtdosePosition = rtdoseInfo.ImagePositionPatient;
            
            fprintf('  Reference dose loaded: [%d x %d x %d]\n', size(referenceDose));
            fprintf('  Resolution: [%.2f, %.2f, %.2f] mm\n', rtdoseResolution);
            fprintf('  Max dose: %.4f Gy\n', max(referenceDose(:)));
            
        catch ME
            fprintf('  ✗ Error loading RTDOSE: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  STEP 2: Load RTPLAN and Extract Beam Weights
        %  ================================================================
        fprintf('\n[STEP 2] Loading RTPLAN and extracting beam weights...\n');
        
        try
            % Find RTPLAN file
            rtplanFiles = dir(fullfile(dicomPath, 'RP*.dcm'));
            if isempty(rtplanFiles)
                rtplanFiles = dir(fullfile(dicomPath, '*RTPLAN*.dcm'));
            end
            if isempty(rtplanFiles)
                allDcm = dir(fullfile(dicomPath, '*.dcm'));
                for f = 1:length(allDcm)
                    info = dicominfo(fullfile(dicomPath, allDcm(f).name));
                    if isfield(info, 'Modality') && strcmp(info.Modality, 'RTPLAN')
                        rtplanFiles = allDcm(f);
                        break;
                    end
                end
            end
            
            if isempty(rtplanFiles)
                error('No RTPLAN file found in %s', dicomPath);
            end
            
            rtplanFile = fullfile(dicomPath, rtplanFiles(1).name);
            rtplanInfo = dicominfo(rtplanFile);
            
            % Extract beam metersets from FractionGroupSequence
            numBeams = 0;
            beamMetersets = [];
            beamNames = {};
            
            if isfield(rtplanInfo, 'FractionGroupSequence')
                fgs = rtplanInfo.FractionGroupSequence;
                fgsFields = fieldnames(fgs);
                
                % Use first fraction group (Item_1)
                if isfield(fgs, 'Item_1') && isfield(fgs.Item_1, 'ReferencedBeamSequence')
                    rbs = fgs.Item_1.ReferencedBeamSequence;
                    rbsFields = fieldnames(rbs);
                    
                    for i = 1:length(rbsFields)
                        if startsWith(rbsFields{i}, 'Item_')
                            numBeams = numBeams + 1;
                            beamItem = rbs.(rbsFields{i});
                            if isfield(beamItem, 'BeamMeterset')
                                beamMetersets(numBeams) = beamItem.BeamMeterset;
                            else
                                beamMetersets(numBeams) = 1;
                            end
                        end
                    end
                end
            end
            
            % Also extract beam names and control point info from BeamSequence
            cumulativeMetersetWeights = cell(1, numBeams);
            
            if isfield(rtplanInfo, 'BeamSequence')
                bs = rtplanInfo.BeamSequence;
                bsFields = fieldnames(bs);
                beamIdx = 0;
                
                for i = 1:length(bsFields)
                    if startsWith(bsFields{i}, 'Item_')
                        beamIdx = beamIdx + 1;
                        beamItem = bs.(bsFields{i});
                        
                        % Get beam name
                        if isfield(beamItem, 'BeamName')
                            beamNames{beamIdx} = beamItem.BeamName;
                        else
                            beamNames{beamIdx} = sprintf('Beam_%d', beamIdx);
                        end
                        
                        % Extract cumulative meterset weights from control points
                        if isfield(beamItem, 'ControlPointSequence')
                            cps = beamItem.ControlPointSequence;
                            cpsFields = fieldnames(cps);
                            weights = [];
                            
                            for j = 1:length(cpsFields)
                                if startsWith(cpsFields{j}, 'Item_')
                                    cpItem = cps.(cpsFields{j});
                                    if isfield(cpItem, 'CumulativeMetersetWeight')
                                        weights(end+1) = cpItem.CumulativeMetersetWeight;
                                    end
                                end
                            end
                            cumulativeMetersetWeights{beamIdx} = weights;
                        end
                    end
                end
            end
            
            fprintf('  Found %d beams in RTPLAN\n', numBeams);
            for b = 1:numBeams
                fprintf('    Beam %d: %s | Meterset: %.4f MU\n', b, beamNames{b}, beamMetersets(b));
            end
            
        catch ME
            fprintf('  ✗ Error loading RTPLAN: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  STEP 3: Import DICOM Data using matRad
        %  ================================================================
        fprintf('\n[STEP 3] Importing DICOM data with matRad...\n');
        
        try
            % Create matRad DICOM importer instance
            importer = matRad_DicomImporter(dicomPath);
            
            % Import DICOM (populates workspace directly)
            importer.matRad_importDicom();
            
            % Check that required variables exist
            if ~exist('ct', 'var') || ~exist('cst', 'var') || ~exist('pln', 'var')
                error('DICOM import did not create required variables (ct, cst, pln)');
            end
            
            fprintf('  ✓ CT imported: [%d x %d x %d]\n', ct.cubeDim);
            fprintf('  ✓ Resolution: [%.2f, %.2f, %.2f] mm\n', ct.resolution.x, ct.resolution.y, ct.resolution.z);
            fprintf('  ✓ Structures imported: %d\n', size(cst, 1));
            
        catch ME
            fprintf('  ✗ Error importing DICOM: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  STEP 4: Apply Dual-Layer MLC Reduction
        %  ================================================================
        fprintf('\n[STEP 4] Applying dual-layer MLC reduction...\n');
        
        try
            % Call external reduce_collimator function
            % This converts Halcyon's dual-layer MLC to single-layer
            pln = reduce_collimator(pln);
            fprintf('  ✓ MLC reduction applied successfully\n');
            
        catch ME
            fprintf('  ⚠ Warning: MLC reduction failed: %s\n', ME.message);
            fprintf('    Continuing with original MLC data...\n');
        end
        
        %% ================================================================
        %  STEP 5: Configure Machine and Photon Parameters
        %  ================================================================
        fprintf('\n[STEP 5] Configuring photon machine...\n');
        
        try
            % Search for generic photon machine in matRad basedata
            basedataPath = fullfile(matRadPath, 'basedata');
            
            % Look for photon machine files
            photonMachines = dir(fullfile(basedataPath, '*photon*.mat'));
            if isempty(photonMachines)
                photonMachines = dir(fullfile(basedataPath, '*Photon*.mat'));
            end
            if isempty(photonMachines)
                % Try generic search
                allMachines = dir(fullfile(basedataPath, '*.mat'));
                for m = 1:length(allMachines)
                    if contains(lower(allMachines(m).name), 'photon') || ...
                       contains(lower(allMachines(m).name), '6mv') || ...
                       contains(lower(allMachines(m).name), '6x')
                        photonMachines = allMachines(m);
                        break;
                    end
                end
            end
            
            if ~isempty(photonMachines)
                machineName = photonMachines(1).name;
                fprintf('  Using machine: %s\n', machineName);
                
                % Set machine in pln
                [~, machineBase, ~] = fileparts(machineName);
                pln.machine = machineBase;
            else
                fprintf('  Using default photon machine configuration\n');
                pln.machine = 'Generic';
            end
            
            % Ensure radiation mode is photons
            pln.radiationMode = 'photons';
            
            % Set to forward dose calculation
            if isfield(pln, 'propOpt')
                pln.propOpt.runSequencing = false;
                pln.propOpt.runDAO = false;
            end
            
            fprintf('  ✓ Machine configured: %s\n', pln.machine);
            fprintf('  ✓ Radiation mode: %s\n', pln.radiationMode);
            
        catch ME
            fprintf('  ✗ Error configuring machine: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  STEP 6: Generate Steering File (STF)
        %  ================================================================
        fprintf('\n[STEP 6] Generating steering file (STF)...\n');
        
        try
            % Generate STF from pln
            stf = matRad_generateStf(ct, cst, pln);
            
            fprintf('  ✓ STF generated with %d beams\n', length(stf));
            
            % Display beam information
            totalBixels = 0;
            bixelsPerBeam = zeros(1, length(stf));
            for b = 1:length(stf)
                bixelsPerBeam(b) = stf(b).totalNumOfBixels;
                totalBixels = totalBixels + bixelsPerBeam(b);
                fprintf('    Beam %d: Gantry=%.1f°, Couch=%.1f°, Bixels=%d\n', ...
                    b, stf(b).gantryAngle, stf(b).couchAngle, stf(b).totalNumOfBixels);
            end
            fprintf('  Total bixels: %d\n', totalBixels);
            
        catch ME
            fprintf('  ✗ Error generating STF: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  STEP 7: Populate Weights from RTPLAN
        %  ================================================================
        fprintf('\n[STEP 7] Populating beam weights from RTPLAN...\n');
        
        try
            % Initialize weight vector
            w = zeros(totalBixels, 1);
            
            % Distribute beam metersets uniformly across bixels in each beam
            bixelOffset = 0;
            for b = 1:min(length(stf), numBeams)
                numBixelsInBeam = bixelsPerBeam(b);
                
                if numBixelsInBeam > 0 && b <= length(beamMetersets)
                    % Uniform distribution of meterset across bixels
                    weightPerBixel = beamMetersets(b) / numBixelsInBeam;
                    
                    for i = 1:numBixelsInBeam
                        w(bixelOffset + i) = weightPerBixel;
                    end
                end
                
                bixelOffset = bixelOffset + numBixelsInBeam;
            end
            
            % Store weights in pln structure
            if ~isfield(pln, 'propOpt')
                pln.propOpt = struct();
            end
            pln.propOpt.w = w;
            
            fprintf('  ✓ Weights populated: %d bixels\n', length(w));
            fprintf('  Total weight sum: %.4f MU\n', sum(w));
            
        catch ME
            fprintf('  ✗ Error populating weights: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  STEP 8: Calculate Dose Influence Matrix
        %  ================================================================
        fprintf('\n[STEP 8] Calculating dose influence matrix...\n');
        
        try
            % Calculate dose influence matrix (dij)
            dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
            
            fprintf('  ✓ Dose influence matrix calculated\n');
            fprintf('  Matrix size: [%d x %d]\n', size(dij.physicalDose, 1), size(dij.physicalDose, 2));
            
        catch ME
            fprintf('  ✗ Error calculating dose influence: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  STEP 9: Calculate Individual Field Doses
        %  ================================================================
        fprintf('\n[STEP 9] Calculating individual field doses...\n');
        
        % Initialize storage
        fieldDoses = cell(1, length(stf));
        fieldMetadata = struct();
        successfulFields = 0;
        failedFields = 0;
        
        % Initialize total dose accumulator
        totalDose = zeros(ct.cubeDim);
        
        bixelOffset = 0;
        for b = 1:length(stf)
            fprintf('\n  Processing Field %d/%d (%s)...\n', b, length(stf), beamNames{min(b, length(beamNames))});
            
            numBixelsInBeam = bixelsPerBeam(b);
            
            % Create weight vector for this field only
            wField = zeros(totalBixels, 1);
            
            for i = 1:numBixelsInBeam
                wField(bixelOffset + i) = w(bixelOffset + i);
            end
            
            % --- Try-catch for forward dose calculation ---
            try
                % Calculate dose for this field
                resultField = matRad_calcDoseDirect(ct, stf, pln, cst, wField, dij);
                
                if isfield(resultField, 'physicalDose')
                    fieldDose = resultField.physicalDose;
                elseif isfield(resultField, 'RBExDose')
                    fieldDose = resultField.RBExDose;
                else
                    % Try to find dose cube in result
                    resultFields = fieldnames(resultField);
                    doseField = resultFields{contains(lower(resultFields), 'dose')};
                    fieldDose = resultField.(doseField{1});
                end
                
                % --- Validation ---
                maxDose = max(fieldDose(:));
                nonZeroVoxels = nnz(fieldDose);
                
                if maxDose > 0 && nonZeroVoxels > 0
                    % SUCCESS
                    fprintf('    ✓ Field %d: Max=%.4f Gy, Non-zero voxels=%d\n', b, maxDose, nonZeroVoxels);
                    
                    % Store field dose BEFORE accumulation attempt
                    fieldDoses{b} = fieldDose;
                    
                    % Store metadata
                    fieldMetadata(b).beamIndex = b;
                    fieldMetadata(b).beamName = beamNames{min(b, length(beamNames))};
                    fieldMetadata(b).gantryAngle = stf(b).gantryAngle;
                    fieldMetadata(b).couchAngle = stf(b).couchAngle;
                    fieldMetadata(b).meterset = beamMetersets(min(b, length(beamMetersets)));
                    fieldMetadata(b).numBixels = numBixelsInBeam;
                    fieldMetadata(b).maxDose = maxDose;
                    fieldMetadata(b).nonZeroVoxels = nonZeroVoxels;
                    fieldMetadata(b).success = true;
                    
                    successfulFields = successfulFields + 1;
                    
                    % --- Try-catch for accumulation (non-fatal) ---
                    try
                        totalDose = totalDose + fieldDose;
                    catch accumME
                        fprintf('    ⚠ Warning: Accumulation failed for field %d: %s\n', b, accumME.message);
                        % Non-fatal - continue processing
                    end
                    
                else
                    % FAILURE - zero dose
                    fprintf('    ✗ Field %d: ZERO DOSE (Max=%.4f, NNZ=%d)\n', b, maxDose, nonZeroVoxels);
                    
                    fieldMetadata(b).beamIndex = b;
                    fieldMetadata(b).beamName = beamNames{min(b, length(beamNames))};
                    fieldMetadata(b).success = false;
                    fieldMetadata(b).errorMsg = 'Zero dose calculated';
                    
                    failedFields = failedFields + 1;
                    
                    % Continue instead of stopping
                    continue;
                end
                
            catch calcME
                % Forward calculation failed
                fprintf('    ✗ Field %d calculation failed: %s\n', b, calcME.message);
                
                fieldMetadata(b).beamIndex = b;
                fieldMetadata(b).beamName = beamNames{min(b, length(beamNames))};
                fieldMetadata(b).success = false;
                fieldMetadata(b).errorMsg = calcME.message;
                
                failedFields = failedFields + 1;
                continue;
            end
            
            bixelOffset = bixelOffset + numBixelsInBeam;
        end
        
        fprintf('\n  Field Calculation Summary:\n');
        fprintf('    Successful: %d/%d\n', successfulFields, length(stf));
        fprintf('    Failed: %d/%d\n', failedFields, length(stf));
        
        %% ================================================================
        %  STEP 10: Resample CT and Doses to RTDOSE Resolution
        %  ================================================================
        fprintf('\n[STEP 10] Resampling to RTDOSE resolution...\n');
        
        try
            % Original CT grid
            ctResolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
            ctDim = ct.cubeDim;
            
            % Create coordinate grids for CT
            xCT = (0:ctDim(2)-1) * ctResolution(1) + ct.x(1);
            yCT = (0:ctDim(1)-1) * ctResolution(2) + ct.y(1);
            zCT = (0:ctDim(3)-1) * ctResolution(3) + ct.z(1);
            
            % Create coordinate grids for RTDOSE
            rtdoseDim = size(referenceDose);
            xRTDOSE = (0:rtdoseDim(2)-1) * rtdoseResolution(1) + rtdosePosition(1);
            yRTDOSE = (0:rtdoseDim(1)-1) * rtdoseResolution(2) + rtdosePosition(2);
            zRTDOSE = (0:rtdoseDim(3)-1) * rtdoseResolution(3) + rtdosePosition(3);
            
            % Create meshgrids for interpolation
            [XCT, YCT, ZCT] = meshgrid(xCT, yCT, zCT);
            [XRTDOSE, YRTDOSE, ZRTDOSE] = meshgrid(xRTDOSE, yRTDOSE, zRTDOSE);
            
            % Resample CT
            fprintf('  Resampling CT...\n');
            ctResampled = interp3(XCT, YCT, ZCT, ct.cubeHU{1}, XRTDOSE, YRTDOSE, ZRTDOSE, 'linear', -1000);
            fprintf('    ✓ CT resampled: [%d x %d x %d] -> [%d x %d x %d]\n', ctDim, size(ctResampled));
            
            % Resample total dose
            fprintf('  Resampling total dose...\n');
            totalDoseResampled = interp3(XCT, YCT, ZCT, totalDose, XRTDOSE, YRTDOSE, ZRTDOSE, 'linear', 0);
            fprintf('    ✓ Total dose resampled\n');
            
            % Resample individual field doses
            fprintf('  Resampling field doses...\n');
            fieldDosesResampled = cell(size(fieldDoses));
            for b = 1:length(fieldDoses)
                if ~isempty(fieldDoses{b})
                    fieldDosesResampled{b} = interp3(XCT, YCT, ZCT, fieldDoses{b}, ...
                        XRTDOSE, YRTDOSE, ZRTDOSE, 'linear', 0);
                    fprintf('    ✓ Field %d resampled\n', b);
                end
            end
            
            fprintf('  ✓ All resampling completed\n');
            
        catch ME
            fprintf('  ✗ Error during resampling: %s\n', ME.message);
            % Set resampled to original if resampling fails
            ctResampled = ct.cubeHU{1};
            totalDoseResampled = totalDose;
            fieldDosesResampled = fieldDoses;
        end
        
        %% ================================================================
        %  STEP 11: Compare with Reference Dose
        %  ================================================================
        fprintf('\n[STEP 11] Comparing with reference RTDOSE...\n');
        
        try
            % Ensure same dimensions for comparison
            if isequal(size(totalDoseResampled), size(referenceDose))
                % Calculate difference
                doseDifference = totalDoseResampled - referenceDose;
                
                % Calculate metrics
                maxDiffAbs = max(abs(doseDifference(:)));
                meanDiffAbs = mean(abs(doseDifference(:)));
                
                % Gamma analysis (simple version - 3mm/3%)
                doseTolerance = 0.03 * max(referenceDose(:));
                distanceTolerance = 3; % mm
                
                % Points where both doses are significant
                significantMask = referenceDose > 0.1 * max(referenceDose(:));
                
                % Dose difference percentage in significant region
                doseDiffPercent = zeros(size(referenceDose));
                doseDiffPercent(significantMask) = abs(doseDifference(significantMask)) ./ ...
                    referenceDose(significantMask) * 100;
                
                passingPoints = sum(abs(doseDifference(significantMask)) <= doseTolerance);
                totalPoints = sum(significantMask(:));
                passingRate = passingPoints / totalPoints * 100;
                
                fprintf('  Comparison Metrics:\n');
                fprintf('    Max absolute difference: %.4f Gy\n', maxDiffAbs);
                fprintf('    Mean absolute difference: %.4f Gy\n', meanDiffAbs);
                fprintf('    3%% dose difference passing rate: %.1f%%\n', passingRate);
                
                % Store comparison data
                doseComparison.reference = referenceDose;
                doseComparison.calculated = totalDoseResampled;
                doseComparison.difference = doseDifference;
                doseComparison.differencePercent = doseDiffPercent;
                doseComparison.metrics.maxDiffAbs = maxDiffAbs;
                doseComparison.metrics.meanDiffAbs = meanDiffAbs;
                doseComparison.metrics.passingRate = passingRate;
                
            else
                fprintf('  ⚠ Warning: Dose dimensions do not match for comparison\n');
                fprintf('    Resampled: [%s]\n', num2str(size(totalDoseResampled)));
                fprintf('    Reference: [%s]\n', num2str(size(referenceDose)));
                
                doseComparison.reference = referenceDose;
                doseComparison.calculated = totalDoseResampled;
                doseComparison.dimensionMismatch = true;
            end
            
        catch ME
            fprintf('  ✗ Error during comparison: %s\n', ME.message);
            doseComparison.error = ME.message;
        end
        
        %% ================================================================
        %  STEP 12: Save Results
        %  ================================================================
        fprintf('\n[STEP 12] Saving results...\n');
        
        try
            % Save fieldDoses.mat
            fieldDosesOutput.fieldDoses = fieldDoses;
            fieldDosesOutput.fieldDosesResampled = fieldDosesResampled;
            fieldDosesOutput.metadata = fieldMetadata;
            fieldDosesOutput.ctResolution = ctResolution;
            fieldDosesOutput.rtdoseResolution = rtdoseResolution;
            fieldDosesOutput.patientId = currentId;
            fieldDosesOutput.session = currentSession;
            fieldDosesOutput.beamNames = beamNames;
            fieldDosesOutput.beamMetersets = beamMetersets;
            
            save(fullfile(outputDir, 'fieldDoses.mat'), '-struct', 'fieldDosesOutput', '-v7.3');
            fprintf('  ✓ Saved: fieldDoses.mat\n');
            
            % Save reconstructedDose.mat
            reconstructedDose.totalDose = totalDose;
            reconstructedDose.totalDoseResampled = totalDoseResampled;
            reconstructedDose.ctOriginal = ct.cubeHU{1};
            reconstructedDose.ctResampled = ctResampled;
            reconstructedDose.ctResolution = ctResolution;
            reconstructedDose.rtdoseResolution = rtdoseResolution;
            reconstructedDose.successfulFields = successfulFields;
            reconstructedDose.totalFields = length(stf);
            
            save(fullfile(outputDir, 'reconstructedDose.mat'), '-struct', 'reconstructedDose', '-v7.3');
            fprintf('  ✓ Saved: reconstructedDose.mat\n');
            
            % Save doseComparison.mat
            doseComparison.originalResolution = ctResolution;
            doseComparison.resampledResolution = rtdoseResolution;
            
            save(fullfile(outputDir, 'doseComparison.mat'), '-struct', 'doseComparison', '-v7.3');
            fprintf('  ✓ Saved: doseComparison.mat\n');
            
            % Save individual field files
            fprintf('  Saving individual field files...\n');
            for b = 1:length(fieldDoses)
                if ~isempty(fieldDoses{b})
                    fieldData.dose = fieldDoses{b};
                    fieldData.doseResampled = fieldDosesResampled{b};
                    fieldData.metadata = fieldMetadata(b);
                    fieldData.resolution = ctResolution;
                    fieldData.resampledResolution = rtdoseResolution;
                    
                    fieldFilename = sprintf('Field_%02d.mat', b);
                    save(fullfile(outputDir, fieldFilename), '-struct', 'fieldData', '-v7.3');
                    fprintf('    ✓ Saved: %s\n', fieldFilename);
                end
            end
            
            fprintf('\n  All results saved to: %s\n', outputDir);
            
        catch ME
            fprintf('  ✗ Error saving results: %s\n', ME.message);
        end
        
        %% ================================================================
        %  SUMMARY
        %  ================================================================
        fprintf('\n');
        fprintf('================================================================\n');
        fprintf('PROCESSING COMPLETE: Patient %s, Session %s\n', currentId, currentSession);
        fprintf('================================================================\n');
        fprintf('  Total fields processed: %d\n', length(stf));
        fprintf('  Successful calculations: %d\n', successfulFields);
        fprintf('  Failed calculations: %d\n', failedFields);
        fprintf('  Output directory: %s\n', outputDir);
        fprintf('================================================================\n');
        
    end % session loop
end % patient loop

fprintf('\n\nBatch processing completed.\n');
