%% ETHOS Field-by-Field Dose Reconstruction using MATRAD
% =========================================================================
% Purpose: Reconstruct individual field doses from ETHOS IMRT treatments
%          using MATRAD dose calculation engine.
%
% This script performs independent dose calculation for validation of
% ETHOS treatment plans by separating individual beam contributions.
%
% Author: Generated for ETHOS dose analysis
% Date: 2025
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Batch processing lists (expandable for multiple patients/sessions)
patientIDs = {'1194203'};
sessions = {'Session_1'};

% Directory configuration
workingDir = '/mnt/weka/home/80030361/ETHOS_Simulations';
matRadPath = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

% Output configuration
saveIndividualFields = true;
generateVisualizations = true;
verboseOutput = true;

% Add MATRAD to path
addpath(genpath(matRadPath));

fprintf('=============================================================\n');
fprintf('  ETHOS Field-by-Field Dose Reconstruction\n');
fprintf('  Using MATRAD Dose Calculation Engine\n');
fprintf('=============================================================\n');
fprintf('  Working Directory: %s\n', workingDir);
fprintf('  MATRAD Path: %s\n', matRadPath);
fprintf('  Patients to process: %d\n', length(patientIDs));
fprintf('=============================================================\n\n');

%% ========================================================================
%  MAIN PROCESSING LOOP
%  ========================================================================

for patIdx = 1:length(patientIDs)
    for sessIdx = 1:length(sessions)
        
        patientID = patientIDs{patIdx};
        sessionName = sessions{sessIdx};
        
        fprintf('\n#############################################################\n');
        fprintf('  Processing Patient: %s | Session: %s\n', patientID, sessionName);
        fprintf('#############################################################\n');
        
        % Define paths
        rawDataDir = fullfile(workingDir, 'EthosExports', patientID, 'Pancreas', sessionName);
        dicomDir = fullfile(rawDataDir, 'sct');
        outputDir = fullfile(workingDir, 'Results', patientID, sessionName);
        
        % Verify directories exist
        if ~isfolder(dicomDir)
            fprintf('  ✗ DICOM directory not found: %s\n', dicomDir);
            fprintf('  Skipping this patient/session.\n');
            continue;
        end
        
        % Create output directory
        if ~isfolder(outputDir)
            mkdir(outputDir);
            fprintf('  Created output directory: %s\n', outputDir);
        end
        
        %% ================================================================
        %  PHASE 1: LOCATE DICOM FILES
        %  ================================================================
        fprintf('\n--- Phase 1: Locating DICOM Files ---\n');
        
        try
            % Find RTPLAN file
            rtplanFiles = dir(fullfile(dicomDir, '*RTPLAN*'));
            if isempty(rtplanFiles)
                rtplanFiles = dir(fullfile(dicomDir, 'RP*.dcm'));
            end
            if isempty(rtplanFiles)
                % Search all DCM files for RTPLAN modality
                allDcmFiles = dir(fullfile(dicomDir, '*.dcm'));
                for i = 1:length(allDcmFiles)
                    info = dicominfo(fullfile(dicomDir, allDcmFiles(i).name));
                    if isfield(info, 'Modality') && strcmp(info.Modality, 'RTPLAN')
                        rtplanFiles = allDcmFiles(i);
                        break;
                    end
                end
            end
            
            if isempty(rtplanFiles)
                error('RTPLAN file not found in %s', dicomDir);
            end
            rtplanFile = fullfile(dicomDir, rtplanFiles(1).name);
            fprintf('  ✓ RTPLAN: %s\n', rtplanFiles(1).name);
            
            % Find RTDOSE file
            rtdoseFiles = dir(fullfile(dicomDir, '*RTDOSE*'));
            if isempty(rtdoseFiles)
                rtdoseFiles = dir(fullfile(dicomDir, 'RD*.dcm'));
            end
            if isempty(rtdoseFiles)
                allDcmFiles = dir(fullfile(dicomDir, '*.dcm'));
                for i = 1:length(allDcmFiles)
                    info = dicominfo(fullfile(dicomDir, allDcmFiles(i).name));
                    if isfield(info, 'Modality') && strcmp(info.Modality, 'RTDOSE')
                        rtdoseFiles = allDcmFiles(i);
                        break;
                    end
                end
            end
            
            if isempty(rtdoseFiles)
                error('RTDOSE file not found in %s', dicomDir);
            end
            rtdoseFile = fullfile(dicomDir, rtdoseFiles(1).name);
            fprintf('  ✓ RTDOSE: %s\n', rtdoseFiles(1).name);
            
            % Find RTSTRUCT file
            rtstructFiles = dir(fullfile(dicomDir, '*RTSTRUCT*'));
            if isempty(rtstructFiles)
                rtstructFiles = dir(fullfile(dicomDir, 'RS*.dcm'));
            end
            if isempty(rtstructFiles)
                allDcmFiles = dir(fullfile(dicomDir, '*.dcm'));
                for i = 1:length(allDcmFiles)
                    info = dicominfo(fullfile(dicomDir, allDcmFiles(i).name));
                    if isfield(info, 'Modality') && strcmp(info.Modality, 'RTSTRUCT')
                        rtstructFiles = allDcmFiles(i);
                        break;
                    end
                end
            end
            
            if ~isempty(rtstructFiles)
                rtstructFile = fullfile(dicomDir, rtstructFiles(1).name);
                fprintf('  ✓ RTSTRUCT: %s\n', rtstructFiles(1).name);
            else
                rtstructFile = '';
                fprintf('  ⚠ RTSTRUCT: Not found (will proceed without structures)\n');
            end
            
        catch ME
            fprintf('  ✗ Error locating DICOM files: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  PHASE 2: LOAD REFERENCE RTDOSE
        %  ================================================================
        fprintf('\n--- Phase 2: Loading Reference RTDOSE ---\n');
        
        try
            rtdoseInfo = dicominfo(rtdoseFile);
            rtdoseData = squeeze(dicomread(rtdoseFile));  % Remove singleton dimensions
            
            % Get dose scaling factor
            if isfield(rtdoseInfo, 'DoseGridScaling')
                doseScaling = rtdoseInfo.DoseGridScaling;
            else
                doseScaling = 1;
            end
            
            % Convert to physical dose (Gy)
            referenceDose = double(rtdoseData) * doseScaling;
            
            % Extract RTDOSE grid information
            rtdoseGrid = struct();
            rtdoseGrid.PixelSpacing = rtdoseInfo.PixelSpacing;  % [row, col] spacing in mm
            
            % Z resolution from GridFrameOffsetVector (multiframe DICOM)
            if isfield(rtdoseInfo, 'GridFrameOffsetVector') && length(rtdoseInfo.GridFrameOffsetVector) > 1
                rtdoseGrid.ZSpacing = abs(rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1));
            elseif isfield(rtdoseInfo, 'SliceThickness')
                rtdoseGrid.ZSpacing = rtdoseInfo.SliceThickness;
            else
                rtdoseGrid.ZSpacing = rtdoseGrid.PixelSpacing(1);  % Fallback to pixel spacing
                fprintf('  ⚠ Z-spacing not found, using pixel spacing as fallback\n');
            end
            
            rtdoseGrid.ImagePositionPatient = rtdoseInfo.ImagePositionPatient;
            rtdoseGrid.Rows = rtdoseInfo.Rows;
            rtdoseGrid.Columns = rtdoseInfo.Columns;
            rtdoseGrid.NumberOfFrames = size(referenceDose, 3);
            rtdoseGrid.Resolution = [rtdoseGrid.PixelSpacing(2), rtdoseGrid.PixelSpacing(1), rtdoseGrid.ZSpacing];
            
            fprintf('  ✓ Reference dose loaded: [%d x %d x %d]\n', size(referenceDose));
            fprintf('    Max dose: %.4f Gy\n', max(referenceDose(:)));
            fprintf('    Resolution: [%.2f, %.2f, %.2f] mm\n', rtdoseGrid.Resolution);
            fprintf('    Origin: [%.2f, %.2f, %.2f] mm\n', rtdoseGrid.ImagePositionPatient);
            
        catch ME
            fprintf('  ✗ Error loading RTDOSE: %s\n', ME.message);
            continue;
        end
        
        %% ================================================================
        %  PHASE 3: IMPORT DICOM DATA USING MATRAD
        %  ================================================================
        fprintf('\n--- Phase 3: MATRAD DICOM Import ---\n');
        
        try
            % Create DICOM importer instance (object-oriented approach)
            importer = matRad_DicomImporter(dicomDir);
            
            % Execute import (populates workspace directly)
            importer.matRad_importDicom();
            
            % Retrieve imported data from base workspace
            ct = evalin('base', 'ct');
            cst = evalin('base', 'cst');
            pln = evalin('base', 'pln');
            
            fprintf('  ✓ CT imported: [%d x %d x %d]\n', size(ct.cubeHU{1}));
            fprintf('    Resolution: [%.2f, %.2f, %.2f] mm\n', ct.resolution.x, ct.resolution.y, ct.resolution.z);
            fprintf('  ✓ Structures imported: %d\n', size(cst, 1));
            fprintf('  ✓ Plan imported with %d beams\n', pln.propStf.numOfBeams);
            
        catch ME
            fprintf('  ✗ Error during MATRAD import: %s\n', ME.message);
            fprintf('    %s\n', ME.getReport('basic'));
            continue;
        end
        
        %% ================================================================
        %  PHASE 4: EXTRACT BEAM WEIGHTS FROM RTPLAN
        %  ================================================================
        fprintf('\n--- Phase 4: Extracting Beam Weights from RTPLAN ---\n');
        
        try
            rtplanInfo = dicominfo(rtplanFile);
            
            numBeams = pln.propStf.numOfBeams;
            beamWeights = zeros(numBeams, 1);
            beamMetersets = zeros(numBeams, 1);
            beamNames = cell(numBeams, 1);
            
            % Extract from FractionGroupSequence.Item_1.ReferencedBeamSequence
            if isfield(rtplanInfo, 'FractionGroupSequence') && ...
               isfield(rtplanInfo.FractionGroupSequence, 'Item_1') && ...
               isfield(rtplanInfo.FractionGroupSequence.Item_1, 'ReferencedBeamSequence')
                
                refBeamSeq = rtplanInfo.FractionGroupSequence.Item_1.ReferencedBeamSequence;
                refBeamFields = fieldnames(refBeamSeq);
                
                for i = 1:length(refBeamFields)
                    beamField = refBeamFields{i};
                    if startsWith(beamField, 'Item_')
                        beamData = refBeamSeq.(beamField);
                        beamNum = str2double(extractAfter(beamField, 'Item_'));
                        
                        if beamNum <= numBeams && isfield(beamData, 'BeamMeterset')
                            beamMetersets(beamNum) = beamData.BeamMeterset;
                        end
                    end
                end
                fprintf('  ✓ BeamMetersets extracted from FractionGroupSequence\n');
            else
                fprintf('  ⚠ FractionGroupSequence not found\n');
            end
            
            % Also extract beam names and cumulative weights from BeamSequence
            if isfield(rtplanInfo, 'BeamSequence')
                beamSeqFields = fieldnames(rtplanInfo.BeamSequence);
                
                for i = 1:length(beamSeqFields)
                    beamField = beamSeqFields{i};
                    if startsWith(beamField, 'Item_')
                        beamData = rtplanInfo.BeamSequence.(beamField);
                        beamNum = str2double(extractAfter(beamField, 'Item_'));
                        
                        if beamNum <= numBeams
                            % Get beam name
                            if isfield(beamData, 'BeamName')
                                beamNames{beamNum} = beamData.BeamName;
                            elseif isfield(beamData, 'BeamDescription')
                                beamNames{beamNum} = beamData.BeamDescription;
                            else
                                beamNames{beamNum} = sprintf('Beam_%d', beamNum);
                            end
                            
                            % Check ControlPointSequence for CumulativeMetersetWeight
                            if isfield(beamData, 'ControlPointSequence')
                                cpFields = fieldnames(beamData.ControlPointSequence);
                                lastCP = beamData.ControlPointSequence.(cpFields{end});
                                if isfield(lastCP, 'CumulativeMetersetWeight')
                                    beamWeights(beamNum) = lastCP.CumulativeMetersetWeight;
                                end
                            end
                        end
                    end
                end
            end
            
            % Normalize weights if needed
            totalMeterset = sum(beamMetersets);
            if totalMeterset > 0
                normalizedWeights = beamMetersets / totalMeterset;
            else
                normalizedWeights = ones(numBeams, 1) / numBeams;
            end
            
            fprintf('  Beam weights extracted:\n');
            for i = 1:numBeams
                fprintf('    Beam %2d (%s): Meterset = %.4f, Weight = %.4f\n', ...
                    i, beamNames{i}, beamMetersets(i), normalizedWeights(i));
            end
            
        catch ME
            fprintf('  ✗ Error extracting weights: %s\n', ME.message);
            fprintf('    Using uniform weights as fallback\n');
            normalizedWeights = ones(numBeams, 1) / numBeams;
            beamMetersets = ones(numBeams, 1);
            for i = 1:numBeams
                beamNames{i} = sprintf('Beam_%d', i);
            end
        end
        
        %% ================================================================
        %  PHASE 5: APPLY DUAL-LAYER MLC REDUCTION
        %  ================================================================
        fprintf('\n--- Phase 5: Dual-Layer MLC Reduction ---\n');
        
        try
            % Call the reduce_collimator function
            pln = reduce_collimator(pln, rtplanInfo, verboseOutput);
            
        catch ME
            fprintf('  ✗ Error in MLC reduction: %s\n', ME.message);
            fprintf('    Proceeding without MLC reduction\n');
        end
        
        %% ================================================================
        %  PHASE 6: FIND AND LOAD GENERIC PHOTON MACHINE
        %  ================================================================
        fprintf('\n--- Phase 6: Loading Machine Data ---\n');
        
        try
            % Search for generic photon machine in MATRAD basedata
            basedataPath = fullfile(matRadPath, 'basedata');
            machineFiles = dir(fullfile(basedataPath, '*photons*.mat'));
            
            if isempty(machineFiles)
                machineFiles = dir(fullfile(basedataPath, '*Photon*.mat'));
            end
            if isempty(machineFiles)
                machineFiles = dir(fullfile(basedataPath, '*6MV*.mat'));
            end
            if isempty(machineFiles)
                machineFiles = dir(fullfile(basedataPath, '*.mat'));
                % Filter for photon-related files
                photonIdx = cellfun(@(x) contains(lower(x), 'photon') || ...
                    contains(lower(x), '6mv') || contains(lower(x), 'generic'), ...
                    {machineFiles.name});
                machineFiles = machineFiles(photonIdx);
            end
            
            if ~isempty(machineFiles)
                machineFile = fullfile(basedataPath, machineFiles(1).name);
                machine = load(machineFile);
                
                % Extract machine name
                machineFieldNames = fieldnames(machine);
                if ~isempty(machineFieldNames)
                    machineName = machineFieldNames{1};
                    pln.machine = machineName;
                    fprintf('  ✓ Machine loaded: %s\n', machineName);
                    fprintf('    File: %s\n', machineFiles(1).name);
                end
            else
                fprintf('  ⚠ No photon machine file found\n');
                fprintf('    Using default MATRAD configuration\n');
                pln.machine = 'Generic';
            end
            
        catch ME
            fprintf('  ✗ Error loading machine: %s\n', ME.message);
            pln.machine = 'Generic';
        end
        
        %% ================================================================
        %  PHASE 7: GENERATE STEERING FILE (STF)
        %  ================================================================
        fprintf('\n--- Phase 7: Generating Steering File ---\n');
        
        try
            % Ensure required plan properties
            if ~isfield(pln, 'radiationMode')
                pln.radiationMode = 'photons';
            end
            
            % Generate steering file
            stf = matRad_generateStf(ct, cst, pln);
            
            fprintf('  ✓ Steering file generated\n');
            fprintf('    Number of beams: %d\n', length(stf));
            
            % Display beam information
            for i = 1:length(stf)
                fprintf('    Beam %d: Gantry=%.1f°, Couch=%.1f°, Rays=%d\n', ...
                    i, stf(i).gantryAngle, stf(i).couchAngle, stf(i).numOfRays);
            end
            
        catch ME
            fprintf('  ✗ Error generating steering file: %s\n', ME.message);
            fprintf('    %s\n', ME.getReport('basic'));
            continue;
        end
        
        %% ================================================================
        %  PHASE 8: CALCULATE DOSE INFLUENCE MATRIX
        %  ================================================================
        fprintf('\n--- Phase 8: Calculating Dose Influence Matrix ---\n');
        
        try
            % Calculate dose influence matrix (dij)
            dij = matRad_calcPhotonDose(ct, stf, pln, cst);
            
            fprintf('  ✓ Dose influence matrix calculated\n');
            fprintf('    Total bixels: %d\n', dij.totalNumOfBixels);
            fprintf('    Dose grid size: [%d x %d x %d]\n', dij.doseGrid.dimensions);
            
        catch ME
            fprintf('  ✗ Error calculating dose influence: %s\n', ME.message);
            fprintf('    %s\n', ME.getReport('basic'));
            continue;
        end
        
        %% ================================================================
        %  PHASE 9: FIELD-BY-FIELD DOSE CALCULATION
        %  ================================================================
        fprintf('\n--- Phase 9: Field-by-Field Dose Calculation ---\n');
        
        % Initialize storage
        fieldDoses = cell(numBeams, 1);
        fieldMetadata = struct();
        fieldMetadata.beamNames = beamNames;
        fieldMetadata.beamMetersets = beamMetersets;
        fieldMetadata.normalizedWeights = normalizedWeights;
        fieldMetadata.calculationSuccess = false(numBeams, 1);
        fieldMetadata.maxDose = zeros(numBeams, 1);
        fieldMetadata.nonZeroVoxels = zeros(numBeams, 1);
        
        totalDose = zeros(dij.doseGrid.dimensions);
        successfulFields = 0;
        
        for beamIdx = 1:numBeams
            fprintf('\n  [Beam %d/%d] %s\n', beamIdx, numBeams, beamNames{beamIdx});
            
            % === Dose Influence Extraction ===
            try
                % Find bixels belonging to this beam
                beamBixelStart = 1;
                beamBixelEnd = 0;
                
                for b = 1:beamIdx
                    beamBixelStart = beamBixelEnd + 1;
                    beamBixelEnd = beamBixelEnd + stf(b).totalNumOfBixels;
                end
                numBixels = stf(beamIdx).totalNumOfBixels;
                
                fprintf('    Bixels: %d to %d (%d total)\n', beamBixelStart, beamBixelEnd, numBixels);
                
            catch ME
                fprintf('    ✗ Error extracting beam bixels: %s\n', ME.message);
                continue;
            end
            
            % === Forward Dose Calculation ===
            try
                % Create weight vector for this beam only
                w = zeros(dij.totalNumOfBixels, 1);
                
                % Apply beam meterset uniformly across bixels
                if beamMetersets(beamIdx) > 0
                    bixelWeight = beamMetersets(beamIdx) / numBixels;
                else
                    bixelWeight = 1 / numBixels;
                end
                
                w(beamBixelStart:beamBixelEnd) = bixelWeight;
                
                % Calculate dose for this beam
                resultGUI = matRad_calcDoseDirect(ct, stf, pln, cst, w);
                
                % Extract dose cube
                if isfield(resultGUI, 'physicalDose')
                    beamDose = resultGUI.physicalDose;
                elseif isfield(resultGUI, 'RBExDose')
                    beamDose = resultGUI.RBExDose;
                else
                    error('No dose field found in result');
                end
                
            catch ME
                fprintf('    ✗ Error in dose calculation: %s\n', ME.message);
                continue;
            end
            
            % === Validation ===
            maxDoseVal = max(beamDose(:));
            nonZeroCount = nnz(beamDose);
            
            fieldMetadata.maxDose(beamIdx) = maxDoseVal;
            fieldMetadata.nonZeroVoxels(beamIdx) = nonZeroCount;
            
            if maxDoseVal > 0 && nonZeroCount > 0
                fieldDoses{beamIdx} = beamDose;
                fieldMetadata.calculationSuccess(beamIdx) = true;
                successfulFields = successfulFields + 1;
                
                fprintf('    ✓ Success: Max=%.4f Gy, NonZero=%d voxels\n', maxDoseVal, nonZeroCount);
            else
                fprintf('    ✗ Validation FAILED: Max=%.4f, NonZero=%d\n', maxDoseVal, nonZeroCount);
                fprintf('      Beam produced zero dose - check MLC/collimation settings\n');
            end
            
            % === Accumulation (non-fatal) ===
            try
                if fieldMetadata.calculationSuccess(beamIdx)
                    totalDose = totalDose + beamDose;
                end
            catch ME
                fprintf('    ⚠ Accumulation warning: %s\n', ME.message);
            end
        end
        
        fprintf('\n  Field calculation summary:\n');
        fprintf('    Successful: %d/%d beams\n', successfulFields, numBeams);
        fprintf('    Total accumulated dose max: %.4f Gy\n', max(totalDose(:)));
        
        %% ================================================================
        %  PHASE 10: RESAMPLE TO RTDOSE RESOLUTION
        %  ================================================================
        fprintf('\n--- Phase 10: Resampling to RTDOSE Grid ---\n');
        
        try
            % Get MATRAD dose grid info
            matradGrid = struct();
            matradGrid.x = ct.resolution.x;
            matradGrid.y = ct.resolution.y;
            matradGrid.z = ct.resolution.z;
            matradGrid.dimensions = size(totalDose);
            
            % Create coordinate vectors for MATRAD grid
            [matradX, matradY, matradZ] = ndgrid(...
                (0:matradGrid.dimensions(1)-1) * matradGrid.x, ...
                (0:matradGrid.dimensions(2)-1) * matradGrid.y, ...
                (0:matradGrid.dimensions(3)-1) * matradGrid.z);
            
            % Create coordinate vectors for RTDOSE grid
            [rtdoseX, rtdoseY, rtdoseZ] = ndgrid(...
                (0:rtdoseGrid.Columns-1) * rtdoseGrid.Resolution(1), ...
                (0:rtdoseGrid.Rows-1) * rtdoseGrid.Resolution(2), ...
                (0:rtdoseGrid.NumberOfFrames-1) * rtdoseGrid.Resolution(3));
            
            % Resample total dose
            totalDoseResampled = interp3(matradY, matradX, matradZ, totalDose, ...
                rtdoseY, rtdoseX, rtdoseZ, 'linear', 0);
            
            fprintf('  ✓ Total dose resampled: [%d x %d x %d] → [%d x %d x %d]\n', ...
                size(totalDose), size(totalDoseResampled));
            
            % Resample individual field doses
            fieldDosesResampled = cell(numBeams, 1);
            for beamIdx = 1:numBeams
                if fieldMetadata.calculationSuccess(beamIdx)
                    fieldDosesResampled{beamIdx} = interp3(matradY, matradX, matradZ, ...
                        fieldDoses{beamIdx}, rtdoseY, rtdoseX, rtdoseZ, 'linear', 0);
                end
            end
            
            fprintf('  ✓ Individual field doses resampled\n');
            
            % Resample CT for comparison
            ctResampled = interp3(matradY, matradX, matradZ, ct.cubeHU{1}, ...
                rtdoseY, rtdoseX, rtdoseZ, 'linear', -1000);
            
            fprintf('  ✓ CT resampled for comparison\n');
            
        catch ME
            fprintf('  ✗ Error during resampling: %s\n', ME.message);
            fprintf('    %s\n', ME.getReport('basic'));
            
            % Fallback: use original resolution
            totalDoseResampled = totalDose;
            fieldDosesResampled = fieldDoses;
            ctResampled = ct.cubeHU{1};
        end
        
        %% ================================================================
        %  PHASE 11: SAVE RESULTS
        %  ================================================================
        fprintf('\n--- Phase 11: Saving Results ---\n');
        
        try
            % Save individual field doses
            fieldDosesOutput = struct();
            fieldDosesOutput.doses = fieldDoses;
            fieldDosesOutput.dosesResampled = fieldDosesResampled;
            fieldDosesOutput.metadata = fieldMetadata;
            fieldDosesOutput.grid = matradGrid;
            fieldDosesOutput.rtdoseGrid = rtdoseGrid;
            fieldDosesOutput.patientID = patientID;
            fieldDosesOutput.session = sessionName;
            fieldDosesOutput.timestamp = datetime('now');
            
            save(fullfile(outputDir, 'fieldDoses.mat'), '-struct', 'fieldDosesOutput', '-v7.3');
            fprintf('  ✓ Saved: fieldDoses.mat\n');
            
            % Save reconstructed dose
            reconstructedDoseOutput = struct();
            reconstructedDoseOutput.totalDose = totalDose;
            reconstructedDoseOutput.totalDoseResampled = totalDoseResampled;
            reconstructedDoseOutput.grid = matradGrid;
            reconstructedDoseOutput.rtdoseGrid = rtdoseGrid;
            reconstructedDoseOutput.patientID = patientID;
            reconstructedDoseOutput.session = sessionName;
            reconstructedDoseOutput.timestamp = datetime('now');
            
            save(fullfile(outputDir, 'reconstructedDose.mat'), '-struct', 'reconstructedDoseOutput', '-v7.3');
            fprintf('  ✓ Saved: reconstructedDose.mat\n');
            
            % Save dose comparison
            doseComparisonOutput = struct();
            doseComparisonOutput.referenceDose = referenceDose;
            doseComparisonOutput.reconstructedDose = totalDose;
            doseComparisonOutput.reconstructedDoseResampled = totalDoseResampled;
            doseComparisonOutput.ctResampled = ctResampled;
            doseComparisonOutput.rtdoseGrid = rtdoseGrid;
            doseComparisonOutput.matradGrid = matradGrid;
            
            % Calculate difference metrics
            if isequal(size(totalDoseResampled), size(referenceDose))
                doseDiff = totalDoseResampled - referenceDose;
                doseComparisonOutput.doseDifference = doseDiff;
                doseComparisonOutput.maxDifference = max(abs(doseDiff(:)));
                doseComparisonOutput.meanDifference = mean(doseDiff(:));
                
                % Gamma analysis (simple 3%/3mm)
                % Note: Full gamma requires dedicated function
                fprintf('  Dose difference statistics:\n');
                fprintf('    Max difference: %.4f Gy\n', doseComparisonOutput.maxDifference);
                fprintf('    Mean difference: %.4f Gy\n', doseComparisonOutput.meanDifference);
            else
                fprintf('  ⚠ Grid size mismatch - skipping difference calculation\n');
                fprintf('    Reconstructed: [%d x %d x %d]\n', size(totalDoseResampled));
                fprintf('    Reference: [%d x %d x %d]\n', size(referenceDose));
            end
            
            doseComparisonOutput.patientID = patientID;
            doseComparisonOutput.session = sessionName;
            doseComparisonOutput.timestamp = datetime('now');
            
            save(fullfile(outputDir, 'doseComparison.mat'), '-struct', 'doseComparisonOutput', '-v7.3');
            fprintf('  ✓ Saved: doseComparison.mat\n');
            
            % Save individual field files
            if saveIndividualFields
                fieldOutputDir = fullfile(outputDir, 'IndividualFields');
                if ~isfolder(fieldOutputDir)
                    mkdir(fieldOutputDir);
                end
                
                for beamIdx = 1:numBeams
                    if fieldMetadata.calculationSuccess(beamIdx)
                        fieldOutput = struct();
                        fieldOutput.dose = fieldDoses{beamIdx};
                        fieldOutput.doseResampled = fieldDosesResampled{beamIdx};
                        fieldOutput.beamName = beamNames{beamIdx};
                        fieldOutput.beamMeterset = beamMetersets(beamIdx);
                        fieldOutput.beamIndex = beamIdx;
                        fieldOutput.maxDose = fieldMetadata.maxDose(beamIdx);
                        
                        filename = sprintf('Field_%02d.mat', beamIdx);
                        save(fullfile(fieldOutputDir, filename), '-struct', 'fieldOutput', '-v7.3');
                    end
                end
                fprintf('  ✓ Saved: %d individual field files\n', successfulFields);
            end
            
        catch ME
            fprintf('  ✗ Error saving results: %s\n', ME.message);
        end
        
        %% ================================================================
        %  PHASE 12: GENERATE VISUALIZATIONS (Optional)
        %  ================================================================
        if generateVisualizations
            fprintf('\n--- Phase 12: Generating Visualizations ---\n');
            
            try
                % Create visualization figure
                fig = figure('Position', [100, 100, 1600, 1000], 'Visible', 'off');
                
                % Find central slices
                [~, ~, ~, centralSlice] = max(sum(sum(totalDose, 1), 2));
                if isempty(centralSlice) || centralSlice == 0
                    centralSlice = round(size(totalDose, 3) / 2);
                end
                
                % 1. Total reconstructed dose (axial)
                subplot(2, 3, 1);
                imagesc(totalDose(:, :, centralSlice));
                colorbar; colormap(gca, 'jet');
                title(sprintf('Reconstructed Dose (Slice %d)', centralSlice));
                xlabel('X'); ylabel('Y');
                axis image;
                
                % 2. Reference dose (axial) - if sizes match
                subplot(2, 3, 2);
                if size(referenceDose, 3) >= centralSlice
                    imagesc(referenceDose(:, :, centralSlice));
                    colorbar; colormap(gca, 'jet');
                    title(sprintf('Reference RTDOSE (Slice %d)', centralSlice));
                else
                    text(0.5, 0.5, 'Size mismatch', 'HorizontalAlignment', 'center');
                    title('Reference RTDOSE');
                end
                xlabel('X'); ylabel('Y');
                axis image;
                
                % 3. Dose difference
                subplot(2, 3, 3);
                if isequal(size(totalDoseResampled), size(referenceDose))
                    diffSlice = totalDoseResampled(:, :, centralSlice) - referenceDose(:, :, centralSlice);
                    imagesc(diffSlice);
                    colorbar; colormap(gca, 'coolwarm');
                    caxis([-max(abs(diffSlice(:))), max(abs(diffSlice(:)))]);
                    title('Dose Difference');
                else
                    text(0.5, 0.5, 'Size mismatch', 'HorizontalAlignment', 'center');
                    title('Dose Difference');
                end
                xlabel('X'); ylabel('Y');
                axis image;
                
                % 4. Individual field contributions (bar chart)
                subplot(2, 3, 4);
                bar(fieldMetadata.maxDose);
                xlabel('Beam Number');
                ylabel('Max Dose (Gy)');
                title('Max Dose per Field');
                xticks(1:numBeams);
                
                % 5. Beam weight distribution
                subplot(2, 3, 5);
                bar(normalizedWeights);
                xlabel('Beam Number');
                ylabel('Normalized Weight');
                title('Beam Weight Distribution');
                xticks(1:numBeams);
                
                % 6. Dose profile comparison
                subplot(2, 3, 6);
                centerRow = round(size(totalDose, 1) / 2);
                profile_recon = squeeze(totalDose(centerRow, :, centralSlice));
                plot(profile_recon, 'b-', 'LineWidth', 1.5);
                hold on;
                if isequal(size(totalDoseResampled), size(referenceDose))
                    centerRowRef = round(size(referenceDose, 1) / 2);
                    profile_ref = squeeze(referenceDose(centerRowRef, :, centralSlice));
                    plot(profile_ref, 'r--', 'LineWidth', 1.5);
                    legend('Reconstructed', 'Reference');
                end
                xlabel('Position');
                ylabel('Dose (Gy)');
                title('Central Dose Profile');
                grid on;
                
                % Add overall title
                sgtitle(sprintf('Patient: %s | Session: %s', patientID, sessionName));
                
                % Save figure
                saveas(fig, fullfile(outputDir, 'dose_visualization.png'));
                saveas(fig, fullfile(outputDir, 'dose_visualization.fig'));
                close(fig);
                
                fprintf('  ✓ Visualizations saved\n');
                
            catch ME
                fprintf('  ✗ Error generating visualizations: %s\n', ME.message);
            end
        end
        
        %% ================================================================
        %  SUMMARY
        %  ================================================================
        fprintf('\n=============================================================\n');
        fprintf('  PROCESSING COMPLETE: %s / %s\n', patientID, sessionName);
        fprintf('=============================================================\n');
        fprintf('  Fields calculated: %d/%d successful\n', successfulFields, numBeams);
        fprintf('  Max reconstructed dose: %.4f Gy\n', max(totalDose(:)));
        fprintf('  Max reference dose: %.4f Gy\n', max(referenceDose(:)));
        fprintf('  Output directory: %s\n', outputDir);
        fprintf('=============================================================\n\n');
        
    end  % End session loop
end  % End patient loop

fprintf('\n\n*** ALL PROCESSING COMPLETE ***\n');
fprintf('Processed %d patient(s) across %d session(s)\n', ...
    length(patientIDs), length(sessions));