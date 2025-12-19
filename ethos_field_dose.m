%% ETHOS IMRT Field-by-Field Dose Calculator
% Purpose: Calculate individual field doses from ETHOS exported IMRT data
% Uses MATRAD for dose calculation
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
            % Note: This populates variables directly in workspace without output args
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
                
                % RTDOSE is often 4D with singleton dimension - squeeze it
                referenceDose = squeeze(referenceDose);
                
                % Apply scaling
                referenceDose = double(referenceDose) * rtdoseInfo.DoseGridScaling;
                
                fprintf('  - Reference dose grid: %d x %d x %d\n', ...
                    size(referenceDose, 1), size(referenceDose, 2), size(referenceDose, 3));
                fprintf('  - Max dose: %.2f Gy\n', max(referenceDose(:)));
                
                % Extract dose grid parameters
                % X and Y resolution from PixelSpacing
                doseGrid.resolution(1) = rtdoseInfo.PixelSpacing(1);
                doseGrid.resolution(2) = rtdoseInfo.PixelSpacing(2);
                
                % Z resolution from GridFrameOffsetVector (multiframe DICOM)
                if isfield(rtdoseInfo, 'GridFrameOffsetVector') && length(rtdoseInfo.GridFrameOffsetVector) > 1
                    % Calculate slice spacing from frame offset vector
                    doseGrid.resolution(3) = abs(rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1));
                    fprintf('  - Z resolution calculated from GridFrameOffsetVector: %.3f mm\n', doseGrid.resolution(3));
                elseif isfield(rtdoseInfo, 'SliceThickness')
                    doseGrid.resolution(3) = rtdoseInfo.SliceThickness;
                    fprintf('  - Z resolution from SliceThickness: %.3f mm\n', doseGrid.resolution(3));
                else
                    % Fallback to CT resolution
                    doseGrid.resolution(3) = ct.resolution.z;
                    fprintf('  - Warning: Z resolution not found in RTDOSE, using CT resolution: %.3f mm\n', doseGrid.resolution(3));
                end
                
                doseGrid.dimensions = [rtdoseInfo.Rows, rtdoseInfo.Columns, size(referenceDose, 3)];
                
                fprintf('  - Dose grid resolution: [%.3f, %.3f, %.3f] mm\n', ...
                    doseGrid.resolution(1), doseGrid.resolution(2), doseGrid.resolution(3));
                
                % Resample RTDOSE to CT grid if dimensions don't match
                if ~isequal(size(referenceDose), ct.cubeDim)
                    fprintf('\n  Resampling RTDOSE to CT grid...\n');
                    fprintf('    Original RTDOSE: %d x %d x %d at [%.2f, %.2f, %.2f] mm\n', ...
                        size(referenceDose, 1), size(referenceDose, 2), size(referenceDose, 3), ...
                        doseGrid.resolution(1), doseGrid.resolution(2), doseGrid.resolution(3));
                    fprintf('    Target CT grid: %d x %d x %d at [%.2f, %.2f, %.2f] mm\n', ...
                        ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3), ...
                        ct.resolution.x, ct.resolution.y, ct.resolution.z);
                    
                    % Use index-based interpolation (matRad works in voxel space, not DICOM coordinates)
                    [nRows, nCols, nSlices] = size(referenceDose);
                    
                    % Create index grids for original RTDOSE
                    [iDose, jDose, kDose] = ndgrid(1:nRows, 1:nCols, 1:nSlices);
                    
                    % Create query points for CT grid (scaled by resolution ratios)
                    % This maps CT voxel indices to RTDOSE voxel indices based on resolution
                    scale_i = ct.resolution.y / doseGrid.resolution(1);  % Y direction
                    scale_j = ct.resolution.x / doseGrid.resolution(2);  % X direction
                    scale_k = ct.resolution.z / doseGrid.resolution(3);  % Z direction
                    
                    % Generate query grid in RTDOSE index space
                    iQuery = 1 + (0:ct.cubeDim(1)-1) * scale_i;
                    jQuery = 1 + (0:ct.cubeDim(2)-1) * scale_j;
                    kQuery = 1 + (0:ct.cubeDim(3)-1) * scale_k;
                    
                    [iQ, jQ, kQ] = ndgrid(iQuery, jQuery, kQuery);
                    
                    % Interpolate dose to CT grid
                    fprintf('    Interpolating dose values (index-based)...\n');
                    fprintf('    Resolution scaling: [%.3f, %.3f, %.3f]\n', scale_i, scale_j, scale_k);
                    
                    try
                        referenceDoseResampled = interpn(iDose, jDose, kDose, referenceDose, ...
                            iQ, jQ, kQ, 'linear', 0);
                        
                        fprintf('    ✓ Resampling complete\n');
                        fprintf('      Resampled grid: %d x %d x %d\n', ...
                            size(referenceDoseResampled, 1), size(referenceDoseResampled, 2), size(referenceDoseResampled, 3));
                        fprintf('      Max dose after resampling: %.2f Gy\n', max(referenceDoseResampled(:)));
                        fprintf('      Non-zero voxels: %d\n', nnz(referenceDoseResampled));
                        
                        % Verify resampled dose is not empty
                        if max(referenceDoseResampled(:)) == 0 || nnz(referenceDoseResampled) == 0
                            fprintf('    ✗ ERROR: Resampled dose is empty or all zeros!\n');
                            fprintf('      Using original RTDOSE grid (comparison will note size mismatch)\n');
                            referenceDoseOriginal = referenceDose;
                        else
                            % Store both versions
                            referenceDoseOriginal = referenceDose;
                            referenceDose = referenceDoseResampled;
                            fprintf('    ✓ Resampled dose validated - contains non-zero values\n');
                        end
                        
                    catch ME
                        fprintf('    ✗ ERROR during resampling: %s\n', ME.message);
                        fprintf('    Continuing with original RTDOSE grid (comparison will note size mismatch)\n');
                        referenceDoseOriginal = referenceDose;
                    end
                else
                    fprintf('  - RTDOSE and CT grids already match, no resampling needed\n');
                    referenceDoseOriginal = referenceDose;
                end
                
            else
                fprintf('  Warning: No RTDOSE file found, using CT grid\n');
                doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
                doseGrid.dimensions = ct.cubeDim;
                referenceDose = [];
                referenceDoseOriginal = [];
            end
            
        catch ME
            fprintf('Warning: Could not load RTDOSE: %s\n', ME.message);
            doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
            doseGrid.dimensions = ct.cubeDim;
            referenceDose = [];
            referenceDoseOriginal = [];
        end
        
        %% Step 3: Configure dose calculation
        fprintf('\n[3/6] Configuring dose calculation...\n');
        
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
        fprintf('  - Machine: %s\n', pln.machine);
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
            
            % Diagnose plan structure for weight extraction
            fprintf('\n  Diagnosing plan structure for weights:\n');
            if isfield(pln, 'w')
                fprintf('    - pln.w exists: %d values, sum=%.2f\n', length(pln.w), sum(pln.w));
            else
                fprintf('    - pln.w does NOT exist\n');
            end
            
            if isfield(pln, 'propStf')
                fprintf('    - pln.propStf exists\n');
                if isfield(pln.propStf, 'beam')
                    fprintf('      - pln.propStf.beam exists with %d beams\n', length(pln.propStf.beam));
                    if isfield(pln.propStf.beam, 'weight')
                        fprintf('      - pln.propStf.beam(1).weight exists\n');
                    end
                end
            end
            
            % Check stf structure for any weight information
            if isfield(stf, 'ray')
                if isfield(stf(1).ray, 'weight')
                    fprintf('    - stf(1).ray(1).weight exists: %.4f\n', stf(1).ray(1).weight);
                end
                if isfield(stf(1).ray, 'energy')
                    fprintf('    - stf(1).ray has %d energy levels\n', length(stf(1).ray(1).energy));
                end
            end
            
        catch ME
            fprintf('Error generating steering file: %s\n', ME.message);
            continue;
        end
        
        %% Step 5: Calculate individual field doses
        fprintf('\n[5/6] Calculating individual field doses...\n');
        
        % Check if weights are available in the plan
        if isfield(pln, 'w')
            fprintf('  - Found weights in plan: %d values\n', length(pln.w));
        else
            fprintf('  - WARNING: No weights found in pln.w, will use uniform weights\n');
        end
        
        % Initialize storage for field doses
        fieldDoses = cell(length(stf), 1);
        totalDose = [];  % Will be initialized after first successful calculation
        calculatedGridSize = [];  % Track the calculation grid size
        
        for beamIdx = 1:length(stf)
            fprintf('  Processing Beam %d/%d (Gantry: %.1f°)...\n', ...
                beamIdx, length(stf), stf(beamIdx).gantryAngle);
            
            % Create temporary plan with single beam
            plnSingle = pln;
            plnSingle.propStf.numOfBeams = 1;
            plnSingle.propStf.isoCenter = stf(beamIdx).isoCenter;
            
            % Create single-beam steering file
            stfSingle = stf(beamIdx);
            
            % Calculate dose influence matrix for this beam
            try
                fprintf('    Calculating dose influence matrix...\n');
                dij = matRad_calcDoseInfluence(ct, cst, stfSingle, plnSingle);
                fprintf('    - Total bixels in dij: %d\n', dij.totalNumOfBixels);
            catch ME
                fprintf('    ✗ ERROR in dose influence calculation: %s\n', ME.message);
                fieldDoses{beamIdx} = [];
                continue;
            end
            
            % Initialize weights vector
            w = ones(dij.totalNumOfBixels, 1);
            
            % Try to extract weights from imported plan
            if isfield(pln, 'w') && ~isempty(pln.w)
                % Calculate offset for this beam in the global weight vector
                offset = 0;
                for b = 1:(beamIdx-1)
                    numRays = stf(b).numOfRays;
                    numEnergies = length(stf(b).ray(1).energy);
                    offset = offset + numRays * numEnergies;
                end
                
                numBixelsThisBeam = stfSingle.numOfRays * length(stfSingle.ray(1).energy);
                fprintf('    - Extracting weights: offset=%d, num_bixels=%d\n', offset, numBixelsThisBeam);
                
                if offset + numBixelsThisBeam <= length(pln.w)
                    w = pln.w(offset+1:offset+numBixelsThisBeam);
                    fprintf('    - Using imported weights (sum=%.2f)\n', sum(w));
                else
                    fprintf('    - WARNING: Not enough weights in pln.w, using uniform weights\n');
                    w(:) = 1.0;
                end
            else
                fprintf('    - Using uniform weights (sum=%.2f)\n', sum(w));
            end
            
            % Calculate dose for this field
            try
                fprintf('    Calculating forward dose...\n');
                resultGUI = matRad_calcDoseForward(ct, cst, stfSingle, plnSingle, w);
            catch ME
                fprintf('    ✗ ERROR in forward dose calculation: %s\n', ME.message);
                fieldDoses{beamIdx} = [];
                continue;
            end
            
            % Validate and store the result
            if ~isfield(resultGUI, 'physicalDose')
                fprintf('    ✗ ERROR: No physicalDose field in result!\n');
                fieldDoses{beamIdx} = [];
                continue;
            end
            
            % Check dose values
            maxDoseField = max(resultGUI.physicalDose(:));
            nonzeroDose = nnz(resultGUI.physicalDose);
            doseSize = size(resultGUI.physicalDose);
            
            fprintf('    - Calculated dose grid: %d x %d x %d\n', ...
                doseSize(1), doseSize(2), doseSize(3));
            
            if maxDoseField == 0 || nonzeroDose == 0
                fprintf('    ✗ ERROR: Calculated dose is all zeros!\n');
                fieldDoses{beamIdx} = [];
                continue;
            end
            
            fprintf('    ✓ Dose calculated: Max = %.2f Gy, Non-zero voxels = %d\n', ...
                maxDoseField, nonzeroDose);
            
            % Store field dose (do this FIRST before any accumulation)
            fieldDoses{beamIdx} = struct();
            fieldDoses{beamIdx}.physicalDose = resultGUI.physicalDose;
            fieldDoses{beamIdx}.gantryAngle = stf(beamIdx).gantryAngle;
            fieldDoses{beamIdx}.couchAngle = stf(beamIdx).couchAngle;
            fieldDoses{beamIdx}.beamIdx = beamIdx;
            fieldDoses{beamIdx}.weights = w;
            fieldDoses{beamIdx}.maxDose = maxDoseField;
            
            fprintf('    ✓ Field dose stored successfully\n');
            
            % Initialize or accumulate total dose
            if isempty(totalDose)
                totalDose = zeros(doseSize);
                calculatedGridSize = doseSize;
                fprintf('    - Initialized totalDose: %d x %d x %d\n', ...
                    doseSize(1), doseSize(2), doseSize(3));
            end
            
            % Verify size compatibility before accumulation
            if isequal(doseSize, calculatedGridSize)
                try
                    totalDose = totalDose + resultGUI.physicalDose;
                    fprintf('    ✓ Added to totalDose (running max: %.2f Gy)\n', max(totalDose(:)));
                catch ME
                    fprintf('    ⚠ WARNING: Could not add to totalDose: %s\n', ME.message);
                    fprintf('      Field dose is still saved individually\n');
                end
            else
                fprintf('    ⚠ WARNING: Size mismatch with totalDose\n');
                fprintf('      Expected: %d x %d x %d\n', calculatedGridSize);
                fprintf('      Got: %d x %d x %d\n', doseSize);
                fprintf('      Field dose is still saved individually\n');
            end
            
            fprintf('\n');
        end
        
        % Final summary
        numSuccessful = sum(~cellfun(@isempty, fieldDoses));
        fprintf('  Summary: %d/%d fields calculated successfully\n', numSuccessful, length(stf));
        
        %% Step 6: Save results
        fprintf('\n[6/6] Saving results...\n');
        
        % Count successful calculations
        numSuccessful = sum(~cellfun(@isempty, fieldDoses));
        fprintf('  - Successfully calculated: %d/%d fields\n', numSuccessful, length(stf));
        
        if numSuccessful == 0
            fprintf('  ✗ ERROR: No fields were successfully calculated!\n');
            fprintf('  - Check MATRAD configuration and beam parameters\n');
            continue;
        end
        
        % Save individual field doses
        save(fullfile(outputPath, 'fieldDoses.mat'), 'fieldDoses', 'stf', 'pln', 'ct', 'cst');
        fprintf('  - Field doses saved to: fieldDoses.mat\n');
        
        % Save reconstructed total dose if available
        if ~isempty(totalDose) && max(totalDose(:)) > 0
            save(fullfile(outputPath, 'reconstructedDose.mat'), 'totalDose', 'calculatedGridSize');
            fprintf('  - Reconstructed total dose saved\n');
            fprintf('    Total max dose: %.2f Gy\n', max(totalDose(:)));
        else
            fprintf('  - WARNING: Total dose is empty or zero, not saved\n');
            fprintf('  - Individual field doses are still available\n');
        end
        
        % Compare with reference if available
        if ~isempty(referenceDose) && ~isempty(totalDose)
            fprintf('\n  Comparing with reference RTDOSE...\n');
            
            % Check size compatibility
            if ~isequal(size(totalDose), size(referenceDose))
                fprintf('    Dose grid size mismatch:\n');
                fprintf('      Calculated: %d x %d x %d\n', size(totalDose));
                fprintf('      Reference:  %d x %d x %d\n', size(referenceDose));
                fprintf('    Note: Grids have different resolutions\n');
                fprintf('      Calculated grid: [%.2f, %.2f, %.2f] mm\n', ...
                    pln.propDoseCalc.doseGrid.resolution.x, ...
                    pln.propDoseCalc.doseGrid.resolution.y, ...
                    pln.propDoseCalc.doseGrid.resolution.z);
                fprintf('      Reference grid: [%.2f, %.2f, %.2f] mm\n', ...
                    doseGrid.resolution(1), doseGrid.resolution(2), doseGrid.resolution(3));
                
                % Save without direct comparison
                comparison.calculated = totalDose;
                comparison.reference = referenceDose;
                if exist('referenceDoseOriginal', 'var')
                    comparison.referenceOriginal = referenceDoseOriginal;
                end
                comparison.note = 'Grid size mismatch - direct comparison not possible';
                comparison.calculatedGrid = calculatedGridSize;
                comparison.referenceGrid = size(referenceDose);
                save(fullfile(outputPath, 'doseComparison.mat'), 'comparison');
            else
                % Sizes match - can compare directly
                doseDiff = totalDose - referenceDose;
                
                fprintf('    Comparison statistics:\n');
                fprintf('      Calculated max dose: %.2f Gy\n', max(totalDose(:)));
                fprintf('      Reference max dose:  %.2f Gy\n', max(referenceDose(:)));
                fprintf('      Mean absolute difference: %.2f Gy\n', mean(abs(doseDiff(:))));
                fprintf('      Max difference: %.2f Gy\n', max(abs(doseDiff(:))));
                fprintf('      RMS difference: %.2f Gy\n', sqrt(mean(doseDiff(:).^2)));
                
                % Calculate relative differences in high dose region
                highDoseThreshold = 0.5 * max(referenceDose(:));
                highDoseMask = referenceDose > highDoseThreshold;
                if any(highDoseMask(:))
                    relativeDiff = abs(doseDiff(highDoseMask)) ./ referenceDose(highDoseMask) * 100;
                    fprintf('      Mean relative diff (>50%% max): %.2f%%\n', mean(relativeDiff));
                end
                
                % Save comparison
                comparison.calculated = totalDose;
                comparison.reference = referenceDose;
                if exist('referenceDoseOriginal', 'var')
                    comparison.referenceOriginal = referenceDoseOriginal;
                    comparison.note = 'Reference dose was resampled to CT grid for comparison';
                end
                comparison.difference = doseDiff;
                comparison.metrics.meanAbsDiff = mean(abs(doseDiff(:)));
                comparison.metrics.maxDiff = max(abs(doseDiff(:)));
                comparison.metrics.rmsDiff = sqrt(mean(doseDiff(:).^2));
                if any(highDoseMask(:))
                    comparison.metrics.meanRelativeDiff = mean(relativeDiff);
                end
                save(fullfile(outputPath, 'doseComparison.mat'), 'comparison');
                
                fprintf('    ✓ Comparison saved to doseComparison.mat\n');
            end
        elseif ~isempty(referenceDose) && isempty(totalDose)
            fprintf('  - Cannot compare: totalDose was not calculated\n');
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