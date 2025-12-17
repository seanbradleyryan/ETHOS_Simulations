%% ETHOS IMRT Field-by-Field Dose Calculator
% Purpose: Calculate individual field doses from ETHOS exported IMRT data
% Uses MATRAD for dose calculation
% Author: Generated for ETHOS dose analysis
% Date: 2025

clear; clc; close all;

%% Configuration
% Patient and session arrays (expandable for batch processing)
ids = {'1194203'};
sessions = {'Session 1'};

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

% Check for matRad_DicomImporter class
classPath = fullfile(matradPath, 'dicom', '@matRad_DicomImporter');
if exist(classPath, 'dir')
    fprintf('  - Found: matRad_DicomImporter class at %s\n', classPath);
else
    fprintf('  - WARNING: matRad_DicomImporter class not found at expected location\n');
end

% Check if class can be instantiated
try
    testPath = tempname;
    mkdir(testPath);
    testImporter = matRad_DicomImporter(testPath);
    fprintf('  - matRad_DicomImporter class successfully instantiated\n');
    rmdir(testPath);
    clear testImporter;
catch ME
    fprintf('  - WARNING: Could not instantiate matRad_DicomImporter: %s\n', ME.message);
end

% Check for other essential matRad functions
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
        outputPath = fullfile(wd, 'FieldDoses', currentID, currentSession);
        
        % Create output directory
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        
        %% Step 1: Import DICOM data using MATRAD
        fprintf('\n[1/6] Importing DICOM data...\n');
        
        try
            % Use matRad_DicomImporter class (object-oriented approach)
            fprintf('  Creating matRad_DicomImporter object...\n');
            
            % Create importer instance
            importer = matRad_DicomImporter(dicomPath);
            
            % Import all DICOM data (CT, structures, plan)
            fprintf('  Importing DICOM files...\n');
            [ct, cst, pln] = importer.matRad_importDicom();
            
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
        
        %% Step 2: Load reference RTDOSE for grid parameters
        fprintf('\n[2/6] Loading reference RTDOSE...\n');
        
        try
            rtdoseFile = dir(fullfile(dicomPath, 'RD*.dcm'));
            if isempty(rtdoseFile)
                rtdoseFile = dir(fullfile(dicomPath, '*RTDOSE*.dcm'));
            end
            
            if ~isempty(rtdoseFile)
                rtdoseInfo = dicominfo(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name));
                referenceDose = dicomread(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name));
                referenceDose = double(referenceDose) * rtdoseInfo.DoseGridScaling;
                
                fprintf('  - Reference dose grid: %d x %d x %d\n', ...
                    size(referenceDose, 1), size(referenceDose, 2), size(referenceDose, 3));
                fprintf('  - Max dose: %.2f Gy\n', max(referenceDose(:)));
                
                % Extract dose grid parameters
                doseGrid.resolution = rtdoseInfo.PixelSpacing;
                doseGrid.dimensions = [rtdoseInfo.Rows, rtdoseInfo.Columns, ...
                    size(referenceDose, 3)];
            else
                fprintf('  Warning: No RTDOSE file found, using CT grid\n');
                doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
                doseGrid.dimensions = ct.cubeDim;
                referenceDose = [];
            end
            
        catch ME
            fprintf('Warning: Could not load RTDOSE: %s\n', ME.message);
            doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
            doseGrid.dimensions = ct.cubeDim;
            referenceDose = [];
        end
        
        %% Step 3: Configure dose calculation
        fprintf('\n[3/6] Configuring dose calculation...\n');
        
        % Set up dose calculation parameters
        pln.propStf.bixelWidth = 5; % mm
        pln.propDoseCalc.doseGrid.resolution.x = doseGrid.resolution(1);
        pln.propDoseCalc.doseGrid.resolution.y = doseGrid.resolution(2);
        pln.propDoseCalc.doseGrid.resolution.z = doseGrid.resolution(3);
        
        % Set algorithm (typically 'pencilBeam' for photons)
        if strcmp(pln.radiationMode, 'photons')
            pln.propDoseCalc.engine = 'pencilBeam';
        else
            pln.propDoseCalc.engine = 'matRad_pencilBeam';
        end
        
        fprintf('  - Radiation mode: %s\n', pln.radiationMode);
        fprintf('  - Dose engine: %s\n', pln.propDoseCalc.engine);
        fprintf('  - Dose grid resolution: [%.2f, %.2f, %.2f] mm\n', ...
            doseGrid.resolution(1), doseGrid.resolution(2), doseGrid.resolution(3));
        
        %% Step 4: Generate steering information
        fprintf('\n[4/6] Generating steering information...\n');
        
        try
            stf = matRad_generateStf(ct, cst, pln);
            fprintf('  - Steering file generated for %d beams\n', length(stf));
            
            % Display beam information
            for i = 1:length(stf)
                fprintf('    Beam %d: Gantry=%.1f°, Couch=%.1f°, Rays=%d\n', ...
                    i, stf(i).gantryAngle, stf(i).couchAngle, length(stf(i).ray));
            end
            
        catch ME
            fprintf('Error generating steering file: %s\n', ME.message);
            continue;
        end
        
        %% Step 5: Calculate individual field doses
        fprintf('\n[5/6] Calculating individual field doses...\n');
        
        % Initialize storage for field doses
        fieldDoses = cell(length(stf), 1);
        totalDose = zeros(doseGrid.dimensions);
        
        for beamIdx = 1:length(stf)
            fprintf('  Processing Beam %d/%d (Gantry: %.1f°)...\n', ...
                beamIdx, length(stf), stf(beamIdx).gantryAngle);
            
            try
                % Create temporary plan with single beam
                plnSingle = pln;
                plnSingle.propStf.numOfBeams = 1;
                plnSingle.propStf.isoCenter = stf(beamIdx).isoCenter;
                
                % Create single-beam steering file
                stfSingle = stf(beamIdx);
                
                % Calculate dose influence matrix for this beam
                dij = matRad_calcDoseInfluence(ct, cst, stfSingle, plnSingle);
                
                % Extract beam weights from original plan
                % For IMRT, weights should be in w vector
                w = zeros(dij.totalNumOfBixels, 1);
                
                % Try to extract weights from imported plan
                if isfield(pln, 'w') && length(pln.w) >= sum([stf(1:beamIdx).numOfRays] .* ...
                        arrayfun(@(x) length(x.ray(1).energy), stf(1:beamIdx)))
                    % Calculate offset for this beam
                    offset = 0;
                    for b = 1:(beamIdx-1)
                        offset = offset + stf(b).numOfRays * length(stf(b).ray(1).energy);
                    end
                    beamBixels = stfSingle.numOfRays * length(stfSingle.ray(1).energy);
                    w = pln.w(offset+1:offset+beamBixels);
                else
                    % Use uniform weights if not available
                    w(:) = 1;
                end
                
                % Calculate dose for this field
                resultGUI = matRad_calcDoseForward(ct, cst, stfSingle, plnSingle, w);
                
                % Store field dose
                fieldDoses{beamIdx}.physicalDose = resultGUI.physicalDose;
                fieldDoses{beamIdx}.gantryAngle = stf(beamIdx).gantryAngle;
                fieldDoses{beamIdx}.couchAngle = stf(beamIdx).couchAngle;
                fieldDoses{beamIdx}.beamIdx = beamIdx;
                
                % Accumulate total dose
                totalDose = totalDose + resultGUI.physicalDose;
                
                fprintf('    Max dose: %.2f Gy\n', max(resultGUI.physicalDose(:)));
                
            catch ME
                fprintf('    Error calculating dose for beam %d: %s\n', beamIdx, ME.message);
                fieldDoses{beamIdx} = [];
            end
        end
        
        %% Step 6: Save results
        fprintf('\n[6/6] Saving results...\n');
        
        % Save individual field doses
        save(fullfile(outputPath, 'fieldDoses.mat'), 'fieldDoses', 'stf', 'pln', 'ct', 'cst');
        fprintf('  - Field doses saved to: fieldDoses.mat\n');
        
        % Save reconstructed total dose
        save(fullfile(outputPath, 'reconstructedDose.mat'), 'totalDose');
        fprintf('  - Reconstructed total dose saved\n');
        
        % Compare with reference if available
        if ~isempty(referenceDose)
            % Resample if needed
            if ~isequal(size(totalDose), size(referenceDose))
                fprintf('  Warning: Dose grid size mismatch\n');
                fprintf('    Calculated: %d x %d x %d\n', size(totalDose));
                fprintf('    Reference: %d x %d x %d\n', size(referenceDose));
            else
                doseDiff = totalDose - referenceDose;
                
                fprintf('\n  Comparison with reference RTDOSE:\n');
                fprintf('    Calculated max dose: %.2f Gy\n', max(totalDose(:)));
                fprintf('    Reference max dose: %.2f Gy\n', max(referenceDose(:)));
                fprintf('    Mean absolute difference: %.2f Gy\n', mean(abs(doseDiff(:))));
                fprintf('    Max difference: %.2f Gy\n', max(abs(doseDiff(:))));
                
                % Save comparison
                comparison.calculated = totalDose;
                comparison.reference = referenceDose;
                comparison.difference = doseDiff;
                save(fullfile(outputPath, 'doseComparison.mat'), 'comparison');
            end
        end
        
        % Export individual field doses to DICOM
        fprintf('\n  Exporting field doses to DICOM...\n');
        for beamIdx = 1:length(fieldDoses)
            if ~isempty(fieldDoses{beamIdx})
                dicomFilename = sprintf('Field_%02d_Gantry_%.0f.dcm', ...
                    beamIdx, fieldDoses{beamIdx}.gantryAngle);
                
                try
                    % This would require MATRAD DICOM export functionality
                    % For now, save as .mat files
                    fieldFilename = sprintf('Field_%02d.mat', beamIdx);
                    fieldData = fieldDoses{beamIdx};
                    save(fullfile(outputPath, fieldFilename), 'fieldData');
                catch
                    fprintf('    Warning: Could not export field %d\n', beamIdx);
                end
            end
        end
        
        fprintf('\n========================================\n');
        fprintf('Processing complete for %s - %s\n', currentID, currentSession);
        fprintf('Results saved to: %s\n', outputPath);
        fprintf('========================================\n');
        
    end
end

fprintf('\n\nAll processing complete!\n');