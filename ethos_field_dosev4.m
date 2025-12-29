%% ETHOS IMRT Field-by-Field Dose Calculator
% Purpose: Calculate individual field doses from ETHOS exported IMRT data
%          and reproduce RTDOSE by producing individual fields from RTPLAN
% Uses MATRAD for dose calculation
% Author: Generated for ETHOS dose analysis
% Date: 2025
%
% Requirements:
%   - MATRAD toolkit (https://github.com/e0404/matRad)
%   - MATLAB Image Processing Toolbox (for DICOM functions)
%   - reduce_collimator.m (companion script for dual-layer MLC reduction)
%
% Output Files:
%   - fieldDoses.mat: All individual field doses with metadata
%   - reconstructedDose.mat: Sum of all fields (total plan dose)
%   - doseComparison.mat: Comparison metrics with reference RTDOSE
%   - Field_XX.mat: Individual field dose files
%   - croppedResampledData.mat: CT and doses cropped/resampled to RTDOSE grid

clear; clc; close all;

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Patient and session arrays (expandable for batch processing)
ids = {'1194203'};
sessions = {'Session_1'};

% Base directory
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
matradPath = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

% Add path to reduce_collimator.m if not in MATLAB path
% (assumes it's in the same directory as this script or in wd)
scriptPath = fileparts(mfilename('fullpath'));
if ~isempty(scriptPath)
    addpath(scriptPath);
end
addpath(wd);

%% ========================================================================
%  VERIFY PATHS AND INITIALIZE MATRAD
%  ========================================================================

fprintf('============================================================\n');
fprintf('  ETHOS IMRT Field-by-Field Dose Calculator\n');
fprintf('============================================================\n\n');

fprintf('Verifying paths...\n');
if ~exist(wd, 'dir')
    error('Working directory does not exist: %s', wd);
else
    fprintf('  ✓ Working directory OK: %s\n', wd);
end

if ~exist(matradPath, 'dir')
    error('matRad path does not exist: %s', matradPath);
else
    fprintf('  ✓ matRad path OK: %s\n', matradPath);
end

% Add MATRAD to path
addpath(genpath(matradPath));
fprintf('  ✓ matRad added to MATLAB path\n');

% Initialize MATRAD (required for some versions)
fprintf('\nInitializing MATRAD...\n');
try
    matRad_rc; % matRad run configuration/initialization
    fprintf('  ✓ MATRAD initialized successfully (matRad_rc)\n');
catch
    fprintf('  - matRad_rc not found, trying matRad_cfg...\n');
    try
        matRad_cfg = MatRad_Config.instance();
        fprintf('  ✓ MATRAD configured successfully (MatRad_Config)\n');
    catch
        fprintf('  ⚠ Direct initialization not available, continuing...\n');
    end
end

% Verify MATRAD functions are available
fprintf('\nVerifying MATRAD installation...\n');

% Check for matRad_DicomImporter class
classPath = fullfile(matradPath, 'dicom', '@matRad_DicomImporter');
if exist(classPath, 'dir')
    fprintf('  ✓ Found: matRad_DicomImporter class\n');
else
    fprintf('  ⚠ WARNING: matRad_DicomImporter class not found at expected location\n');
end

% Check if class can be instantiated
try
    testPath = tempname;
    mkdir(testPath);
    testImporter = matRad_DicomImporter(testPath);
    fprintf('  ✓ matRad_DicomImporter class successfully instantiated\n');
    rmdir(testPath);
    clear testImporter;
catch ME
    fprintf('  ⚠ WARNING: Could not instantiate matRad_DicomImporter: %s\n', ME.message);
end

% Check for other essential matRad functions
essentialFunctions = {'matRad_calcDoseInfluence', 'matRad_generateStf', 'matRad_calcDoseForward'};
for i = 1:length(essentialFunctions)
    if exist(essentialFunctions{i}, 'file')
        fprintf('  ✓ Found: %s\n', essentialFunctions{i});
    else
        fprintf('  ✗ WARNING: %s not found\n', essentialFunctions{i});
    end
end

% Check for reduce_collimator function
if exist('reduce_collimator', 'file')
    fprintf('  ✓ Found: reduce_collimator.m\n');
else
    fprintf('  ⚠ WARNING: reduce_collimator.m not found - dual-layer MLC reduction will be skipped\n');
end

%% ========================================================================
%  MAIN PROCESSING LOOP
%  ========================================================================

for idxID = 1:length(ids)
    for idxSession = 1:length(sessions)
        
        currentID = ids{idxID};
        currentSession = sessions{idxSession};
        
        fprintf('\n\n============================================================\n');
        fprintf('  Processing Patient: %s, Session: %s\n', currentID, currentSession);
        fprintf('============================================================\n');
        
        % Define paths
        rawwd = fullfile(wd, 'EthosExports', currentID, 'Pancreas', currentSession);
        dicomPath = fullfile(rawwd, 'sct');
        outputPath = fullfile(wd, 'FieldDoses', currentID, currentSession);
        
        % Verify DICOM path exists
        if ~exist(dicomPath, 'dir')
            fprintf('  ✗ ERROR: DICOM path does not exist: %s\n', dicomPath);
            fprintf('  Skipping this patient/session...\n');
            continue;
        end
        
        % Create output directory
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
            fprintf('  ✓ Created output directory: %s\n', outputPath);
        end
        
        %% ================================================================
        %  STEP 1: IMPORT DICOM DATA USING MATRAD
        %  ================================================================
        fprintf('\n[1/7] Importing DICOM data...\n');
        
        try
            % Use matRad_DicomImporter class (object-oriented approach)
            fprintf('  Creating matRad_DicomImporter object...\n');
            
            % Create importer instance
            importer = matRad_DicomImporter(dicomPath);
            
            % Import all DICOM data (CT, structures, plan)
            % Note: This populates variables directly in workspace without output args
            fprintf('  Importing DICOM files...\n');
            importer.matRad_importDicom();
            
            fprintf('  ✓ Data imported successfully\n');
            
            % Display CT information
            fprintf('  - CT dimensions: %d x %d x %d\n', ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));
            fprintf('  - CT resolution: %.2f x %.2f x %.2f mm\n', ...
                ct.resolution.x, ct.resolution.y, ct.resolution.z);
            fprintf('  - Number of structures: %d\n', size(cst, 1));
            
            if isfield(pln, 'propStf') && isfield(pln.propStf, 'numOfBeams')
                fprintf('  - Number of beams: %d\n', pln.propStf.numOfBeams);
            elseif isfield(pln, 'numOfBeams')
                fprintf('  - Number of beams: %d\n', pln.numOfBeams);
            end
            
        catch ME
            fprintf('  ✗ ERROR importing DICOM data: %s\n', ME.message);
            fprintf('  Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            continue;
        end
        
        %% ================================================================
        %  STEP 2: LOAD REFERENCE RTDOSE FOR GRID PARAMETERS
        %  ================================================================
        fprintf('\n[2/7] Loading reference RTDOSE...\n');
        
        referenceDose = [];
        referenceDoseOriginal = [];
        doseGrid = struct();
        rtdoseInfo = [];
        
        try
            rtdoseFile = dir(fullfile(dicomPath, 'RD*.dcm'));
            if isempty(rtdoseFile)
                rtdoseFile = dir(fullfile(dicomPath, '*RTDOSE*.dcm'));
            end
            
            if ~isempty(rtdoseFile)
                rtdoseInfo = dicominfo(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name));
                referenceDose = dicomread(fullfile(rtdoseFile(1).folder, rtdoseFile(1).name));
                
                fprintf('  - RTDOSE file: %s\n', rtdoseFile(1).name);
                
                % RTDOSE is often 4D with singleton dimension - squeeze it
                originalSize = size(referenceDose);
                referenceDose = squeeze(referenceDose);
                if ~isequal(originalSize, size(referenceDose))
                    fprintf('  - Squeezed RTDOSE from %s to %s\n', ...
                        mat2str(originalSize), mat2str(size(referenceDose)));
                end
                
                % Apply scaling
                referenceDose = double(referenceDose) * rtdoseInfo.DoseGridScaling;
                
                fprintf('  - Reference dose grid: %d x %d x %d\n', ...
                    size(referenceDose, 1), size(referenceDose, 2), size(referenceDose, 3));
                fprintf('  - Max dose: %.2f Gy\n', max(referenceDose(:)));
                fprintf('  - Non-zero voxels: %d\n', nnz(referenceDose));
                
                % Extract dose grid parameters
                % X and Y resolution from PixelSpacing
                doseGrid.resolution = zeros(1, 3);
                doseGrid.resolution(1) = rtdoseInfo.PixelSpacing(1);
                doseGrid.resolution(2) = rtdoseInfo.PixelSpacing(2);
                
                % Z resolution from GridFrameOffsetVector (multiframe DICOM)
                if isfield(rtdoseInfo, 'GridFrameOffsetVector') && length(rtdoseInfo.GridFrameOffsetVector) > 1
                    % Calculate slice spacing from frame offset vector
                    doseGrid.resolution(3) = abs(rtdoseInfo.GridFrameOffsetVector(2) - rtdoseInfo.GridFrameOffsetVector(1));
                    fprintf('  - Z resolution from GridFrameOffsetVector: %.3f mm\n', doseGrid.resolution(3));
                elseif isfield(rtdoseInfo, 'SliceThickness')
                    doseGrid.resolution(3) = rtdoseInfo.SliceThickness;
                    fprintf('  - Z resolution from SliceThickness: %.3f mm\n', doseGrid.resolution(3));
                else
                    % Fallback to CT resolution
                    doseGrid.resolution(3) = ct.resolution.z;
                    fprintf('  ⚠ Z resolution not found in RTDOSE, using CT resolution: %.3f mm\n', doseGrid.resolution(3));
                end
                
                doseGrid.dimensions = [rtdoseInfo.Rows, rtdoseInfo.Columns, size(referenceDose, 3)];
                
                % Store ImagePositionPatient for coordinate mapping
                if isfield(rtdoseInfo, 'ImagePositionPatient')
                    doseGrid.origin = rtdoseInfo.ImagePositionPatient;
                    fprintf('  - Dose grid origin: [%.2f, %.2f, %.2f] mm\n', ...
                        doseGrid.origin(1), doseGrid.origin(2), doseGrid.origin(3));
                end
                
                % Store GridFrameOffsetVector for Z coordinates
                if isfield(rtdoseInfo, 'GridFrameOffsetVector')
                    doseGrid.frameOffsets = rtdoseInfo.GridFrameOffsetVector;
                end
                
                fprintf('  ✓ Dose grid resolution: [%.3f, %.3f, %.3f] mm\n', ...
                    doseGrid.resolution(1), doseGrid.resolution(2), doseGrid.resolution(3));
                
                % Store original for later comparison
                referenceDoseOriginal = referenceDose;
                
            else
                fprintf('  ⚠ Warning: No RTDOSE file found, using CT grid\n');
                doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
                doseGrid.dimensions = ct.cubeDim;
            end
            
        catch ME
            fprintf('  ⚠ Warning: Could not load RTDOSE: %s\n', ME.message);
            doseGrid.resolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
            doseGrid.dimensions = ct.cubeDim;
        end
        
        %% ================================================================
        %  STEP 3: EXTRACT MLC DATA AND APPLY DUAL-LAYER REDUCTION
        %  ================================================================
        fprintf('\n[3/7] Extracting MLC data and applying dual-layer reduction...\n');
        
        rtplanInfo = [];
        
        try
            % Find RTPLAN file
            rtplanFile = dir(fullfile(dicomPath, 'RP*.dcm'));
            if isempty(rtplanFile)
                rtplanFile = dir(fullfile(dicomPath, '*RTPLAN*.dcm'));
            end
            
            if ~isempty(rtplanFile)
                rtplanInfo = dicominfo(fullfile(rtplanFile(1).folder, rtplanFile(1).name));
                fprintf('  - RTPLAN file: %s\n', rtplanFile(1).name);
                
                % Apply dual-layer MLC reduction using companion function
                if exist('reduce_collimator', 'file')
                    fprintf('  Applying dual-layer MLC reduction...\n');
                    pln = reduce_collimator(pln, rtplanInfo, true);
                else
                    fprintf('  ⚠ reduce_collimator.m not found, extracting MLC without reduction...\n');
                    
                    % Fallback: Basic MLC extraction without dual-layer reduction
                    pln = extract_mlc_basic(pln, rtplanInfo);
                end
            else
                fprintf('  ⚠ No RTPLAN file found\n');
            end
            
        catch ME
            fprintf('  ⚠ Error extracting MLC data: %s\n', ME.message);
            fprintf('    Will proceed with open field calculation\n');
        end
        
        %% ================================================================
        %  STEP 4: EXTRACT WEIGHTS FROM RTPLAN
        %  ================================================================
        fprintf('\n[4/7] Extracting weights from RTPLAN...\n');
        
        beamMetersets = [];
        
        try
            if ~isempty(rtplanInfo)
                % Method 1: Extract from FractionGroupSequence.ReferencedBeamSequence
                if isfield(rtplanInfo, 'FractionGroupSequence')
                    fg = rtplanInfo.FractionGroupSequence.Item_1;
                    if isfield(fg, 'ReferencedBeamSequence')
                        numRefBeams = length(fieldnames(fg.ReferencedBeamSequence));
                        beamMetersets = zeros(numRefBeams, 1);
                        
                        fprintf('  - Extracting metersets from FractionGroupSequence (%d beams)\n', numRefBeams);
                        
                        for i = 1:numRefBeams
                            refBeamField = sprintf('Item_%d', i);
                            refBeam = fg.ReferencedBeamSequence.(refBeamField);
                            if isfield(refBeam, 'BeamMeterset')
                                beamMetersets(i) = refBeam.BeamMeterset;
                                fprintf('    Beam %2d: Meterset = %.4f MU\n', i, beamMetersets(i));
                            end
                        end
                        
                        fprintf('  ✓ Total MU: %.2f\n', sum(beamMetersets));
                    end
                end
                
                % Method 2: Also extract CumulativeMetersetWeight from ControlPointSequence
                if isfield(rtplanInfo, 'BeamSequence')
                    numBeams = length(fieldnames(rtplanInfo.BeamSequence));
                    fprintf('\n  - Extracting control point weights from BeamSequence (%d beams)\n', numBeams);
                    
                    for beamIdx = 1:numBeams
                        beamField = sprintf('Item_%d', beamIdx);
                        beam = rtplanInfo.BeamSequence.(beamField);
                        
                        if isfield(beam, 'ControlPointSequence')
                            numCP = length(fieldnames(beam.ControlPointSequence));
                            
                            % Get final cumulative meterset weight
                            lastCPField = sprintf('Item_%d', numCP);
                            lastCP = beam.ControlPointSequence.(lastCPField);
                            
                            if isfield(lastCP, 'CumulativeMetersetWeight')
                                finalWeight = lastCP.CumulativeMetersetWeight;
                                fprintf('    Beam %2d: FinalCumulativeMetersetWeight = %.4f (%d control points)\n', ...
                                    beamIdx, finalWeight, numCP);
                            end
                        end
                    end
                end
            else
                fprintf('  ⚠ No RTPLAN info available for weight extraction\n');
            end
            
        catch ME
            fprintf('  ⚠ Error during weight extraction: %s\n', ME.message);
        end
        
        %% ================================================================
        %  STEP 5: CONFIGURE DOSE CALCULATION
        %  ================================================================
        fprintf('\n[5/7] Configuring dose calculation...\n');
        
        % Check and override machine specification
        if isfield(pln, 'machine')
            fprintf('  - Original machine: %s\n', pln.machine);
        end
        
        % Check available machine files in MATRAD basedata
        machineDir = fullfile(matradPath, 'basedata');
        if exist(machineDir, 'dir')
            machineFiles = dir(fullfile(machineDir, 'photons*.mat'));
            if ~isempty(machineFiles)
                fprintf('  - Available machine files:\n');
                for i = 1:min(5, length(machineFiles))
                    [~, machineName, ~] = fileparts(machineFiles(i).name);
                    fprintf('    %d. %s\n', i, machineName);
                end
                if length(machineFiles) > 5
                    fprintf('    ... and %d more\n', length(machineFiles) - 5);
                end
                
                % Use first generic machine or specific one
                genericMachines = machineFiles(contains({machineFiles.name}, 'Generic', 'IgnoreCase', true));
                if ~isempty(genericMachines)
                    [~, selectedMachine, ~] = fileparts(genericMachines(1).name);
                    fprintf('  ✓ Using generic machine: %s\n', selectedMachine);
                    pln.machine = selectedMachine;
                else
                    % Use first available photon machine
                    [~, selectedMachine, ~] = fileparts(machineFiles(1).name);
                    fprintf('  ✓ Using machine: %s\n', selectedMachine);
                    pln.machine = selectedMachine;
                end
            else
                fprintf('  ⚠ Warning: No machine files found, using default\n');
                pln.machine = 'Generic';
            end
        else
            fprintf('  ⚠ Warning: Machine directory not found, using default\n');
            pln.machine = 'Generic';
        end
        
        % Ensure propStf and propDoseCalc structures exist
        if ~isfield(pln, 'propStf')
            pln.propStf = struct();
        end
        if ~isfield(pln, 'propDoseCalc')
            pln.propDoseCalc = struct();
        end
        if ~isfield(pln.propDoseCalc, 'doseGrid')
            pln.propDoseCalc.doseGrid = struct();
        end
        if ~isfield(pln.propDoseCalc.doseGrid, 'resolution')
            pln.propDoseCalc.doseGrid.resolution = struct();
        end
        
        % Set up dose calculation parameters
        pln.propStf.bixelWidth = 5; % mm
        pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x;  % Use CT resolution for calculation
        pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y;
        pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z;
        
        % Set algorithm (typically 'pencilBeam' for photons)
        if isfield(pln, 'radiationMode') && strcmp(pln.radiationMode, 'photons')
            pln.propDoseCalc.engine = 'pencilBeam';
        else
            pln.propDoseCalc.engine = 'matRad_pencilBeam';
        end
        
        fprintf('  - Radiation mode: %s\n', pln.radiationMode);
        fprintf('  - Machine: %s\n', pln.machine);
        fprintf('  - Dose engine: %s\n', pln.propDoseCalc.engine);
        fprintf('  - Calculation grid resolution: [%.2f, %.2f, %.2f] mm\n', ...
            pln.propDoseCalc.doseGrid.resolution.x, ...
            pln.propDoseCalc.doseGrid.resolution.y, ...
            pln.propDoseCalc.doseGrid.resolution.z);
        
        %% ================================================================
        %  STEP 6: GENERATE STEERING FILE AND CALCULATE FIELD DOSES
        %  ================================================================
        fprintf('\n[6/7] Generating steering information and calculating field doses...\n');
        
        try
            % Generate steering file
            fprintf('  Generating steering file...\n');
            stf = matRad_generateStf(ct, cst, pln);
            fprintf('  ✓ Steering file generated for %d beams\n', length(stf));
            
            % Display beam information
            for i = 1:length(stf)
                fprintf('    Beam %2d: Gantry=%.1f°, Couch=%.1f°, Rays=%d\n', ...
                    i, stf(i).gantryAngle, stf(i).couchAngle, stf(i).numOfRays);
            end
            
            % Inject MLC shapes into stf structure if available
            fprintf('\n  Injecting MLC shapes into steering file...\n');
            stf = inject_mlc_to_stf(stf, pln);
            
        catch ME
            fprintf('  ✗ ERROR generating steering file: %s\n', ME.message);
            continue;
        end
        
        % Populate pln.w with weights
        fprintf('\n  Populating weight vector (pln.w)...\n');
        pln = populate_weights(pln, stf, beamMetersets);
        
        % Display weight diagnostics
        fprintf('\n  Weight diagnostics:\n');
        if isfield(pln, 'w') && ~isempty(pln.w)
            fprintf('    ✓ pln.w exists: %d values, sum=%.2f, mean=%.4f\n', ...
                length(pln.w), sum(pln.w), mean(pln.w));
        else
            fprintf('    ⚠ pln.w does NOT exist - will use uniform weights\n');
        end
        
        % Initialize storage for field doses
        fieldDoses = cell(length(stf), 1);
        totalDose = [];
        calculatedGridSize = [];
        
        fprintf('\n  Calculating individual field doses...\n');
        fprintf('  ----------------------------------------\n');
        
        for beamIdx = 1:length(stf)
            fprintf('\n  Processing Beam %d/%d (Gantry: %.1f°)...\n', ...
                beamIdx, length(stf), stf(beamIdx).gantryAngle);
            
            % Create temporary plan with single beam
            plnSingle = pln;
            plnSingle.propStf.numOfBeams = 1;
            plnSingle.propStf.isoCenter = stf(beamIdx).isoCenter;
            
            % Create single-beam steering file
            stfSingle = stf(beamIdx);
            
            % Calculate dose influence matrix for this beam
            dij = [];
            try
                fprintf('    Calculating dose influence matrix...\n');
                dij = matRad_calcDoseInfluence(ct, cst, stfSingle, plnSingle);
                fprintf('    ✓ Total bixels in dij: %d\n', dij.totalNumOfBixels);
            catch ME
                fprintf('    ✗ ERROR in dose influence calculation: %s\n', ME.message);
                fieldDoses{beamIdx} = [];
                continue;
            end
            
            % Initialize weights vector
            w = ones(dij.totalNumOfBixels, 1);
            
            % Extract weights for this beam from pln.w
            try
                if isfield(pln, 'w') && ~isempty(pln.w)
                    % Calculate offset for this beam in the global weight vector
                    offset = 0;
                    for b = 1:(beamIdx-1)
                        numRays = stf(b).numOfRays;
                        if isfield(stf(b).ray(1), 'energy')
                            numEnergies = length(stf(b).ray(1).energy);
                        else
                            numEnergies = 1;
                        end
                        offset = offset + numRays * numEnergies;
                    end
                    
                    if isfield(stfSingle.ray(1), 'energy')
                        numEnergies = length(stfSingle.ray(1).energy);
                    else
                        numEnergies = 1;
                    end
                    numBixelsThisBeam = stfSingle.numOfRays * numEnergies;
                    
                    fprintf('    - Extracting weights: offset=%d, num_bixels=%d\n', offset, numBixelsThisBeam);
                    
                    if offset + numBixelsThisBeam <= length(pln.w)
                        w = pln.w(offset+1:offset+numBixelsThisBeam);
                        fprintf('    ✓ Using extracted weights (sum=%.2f, mean=%.4f)\n', sum(w), mean(w));
                    else
                        fprintf('    ⚠ Not enough weights in pln.w, using uniform weights\n');
                        w(:) = 1.0;
                    end
                else
                    fprintf('    - Using uniform weights (sum=%.2f)\n', sum(w));
                end
            catch ME
                fprintf('    ⚠ Error extracting weights: %s\n', ME.message);
                fprintf('    - Using uniform weights\n');
                w = ones(dij.totalNumOfBixels, 1);
            end
            
            % Calculate dose for this field
            resultGUI = [];
            try
                fprintf('    Calculating forward dose...\n');
                resultGUI = matRad_calcDoseForward(ct, cst, stfSingle, plnSingle, w);
            catch ME
                fprintf('    ✗ ERROR in forward dose calculation: %s\n', ME.message);
                fieldDoses{beamIdx} = [];
                continue;
            end
            
            % Validate result
            if ~isfield(resultGUI, 'physicalDose')
                fprintf('    ✗ ERROR: No physicalDose field in result!\n');
                fieldDoses{beamIdx} = [];
                continue;
            end
            
            % Check dose values IMMEDIATELY
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
            
            fprintf('    ✓ Dose calculated: Max = %.4f Gy, Non-zero voxels = %d\n', ...
                maxDoseField, nonzeroDose);
            
            % STORE FIELD DOSE FIRST (before any accumulation)
            fieldDoses{beamIdx} = struct();
            fieldDoses{beamIdx}.physicalDose = resultGUI.physicalDose;
            fieldDoses{beamIdx}.gantryAngle = stf(beamIdx).gantryAngle;
            fieldDoses{beamIdx}.couchAngle = stf(beamIdx).couchAngle;
            fieldDoses{beamIdx}.beamIdx = beamIdx;
            fieldDoses{beamIdx}.weights = w;
            fieldDoses{beamIdx}.maxDose = maxDoseField;
            fieldDoses{beamIdx}.nonzeroVoxels = nonzeroDose;
            
            if ~isempty(beamMetersets) && beamIdx <= length(beamMetersets)
                fieldDoses{beamIdx}.meterset = beamMetersets(beamIdx);
            end
            
            fprintf('    ✓ Field dose stored successfully\n');
            
            % Initialize or accumulate total dose (in separate try-catch)
            try
                if isempty(totalDose)
                    totalDose = zeros(doseSize);
                    calculatedGridSize = doseSize;
                    fprintf('    - Initialized totalDose: %d x %d x %d\n', ...
                        doseSize(1), doseSize(2), doseSize(3));
                end
                
                % Verify size compatibility before accumulation
                if isequal(doseSize, calculatedGridSize)
                    totalDose = totalDose + resultGUI.physicalDose;
                    fprintf('    ✓ Added to totalDose (running max: %.4f Gy)\n', max(totalDose(:)));
                else
                    fprintf('    ⚠ WARNING: Size mismatch with totalDose\n');
                    fprintf('      Expected: %d x %d x %d\n', calculatedGridSize(1), calculatedGridSize(2), calculatedGridSize(3));
                    fprintf('      Got: %d x %d x %d\n', doseSize(1), doseSize(2), doseSize(3));
                    fprintf('      Field dose is saved individually (not added to total)\n');
                end
            catch ME
                fprintf('    ⚠ WARNING: Could not accumulate to totalDose: %s\n', ME.message);
                fprintf('      Field dose is still saved individually\n');
            end
        end
        
        % Field calculation summary
        numSuccessful = sum(~cellfun(@isempty, fieldDoses));
        fprintf('\n  ----------------------------------------\n');
        fprintf('  Field Calculation Summary: %d/%d fields calculated successfully\n', ...
            numSuccessful, length(stf));
        
        if numSuccessful == 0
            fprintf('  ✗ ERROR: No fields were successfully calculated!\n');
            fprintf('  - Check MATRAD configuration and beam parameters\n');
            continue;
        end
        
        %% ================================================================
        %  STEP 7: CROP CT TO RTDOSE EXTENT, RESAMPLE, AND SAVE RESULTS
        %  ================================================================
        fprintf('\n[7/7] Cropping to RTDOSE extent, resampling, and saving results...\n');
        
        % Save individual field doses (on CT grid)
        save(fullfile(outputPath, 'fieldDoses.mat'), 'fieldDoses', 'stf', 'pln', 'ct', 'cst', '-v7.3');
        fprintf('  ✓ Field doses saved to: fieldDoses.mat\n');
        
        % Save reconstructed total dose (on CT grid)
        if ~isempty(totalDose) && max(totalDose(:)) > 0
            save(fullfile(outputPath, 'reconstructedDose.mat'), 'totalDose', 'calculatedGridSize', '-v7.3');
            fprintf('  ✓ Reconstructed total dose saved: max = %.4f Gy\n', max(totalDose(:)));
        else
            fprintf('  ⚠ WARNING: Total dose is empty or zero, not saved\n');
        end
        
        % Save individual field files
        fprintf('\n  Saving individual field files...\n');
        for beamIdx = 1:length(fieldDoses)
            if ~isempty(fieldDoses{beamIdx})
                fieldFilename = sprintf('Field_%02d.mat', beamIdx);
                fieldData = fieldDoses{beamIdx};
                save(fullfile(outputPath, fieldFilename), 'fieldData', '-v7.3');
                fprintf('    ✓ Saved: %s (max=%.4f Gy)\n', fieldFilename, fieldData.maxDose);
            end
        end
        
        % Crop and resample all doses to match RTDOSE grid
        fprintf('\n  Cropping CT to RTDOSE extent and resampling...\n');
        
        croppedResampledData = struct();
        croppedResampledData.ctResolution = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
        croppedResampledData.doseResolution = doseGrid.resolution;
        croppedResampledData.doseGridDimensions = doseGrid.dimensions;
        
        if ~isempty(referenceDoseOriginal) && ~isempty(totalDose)
            try
                % Crop and resample totalDose to match RTDOSE grid
                fprintf('    Cropping and resampling totalDose to RTDOSE grid...\n');
                
                [totalDoseCropped, cropResampleInfo] = crop_and_resample_to_rtdose(...
                    totalDose, ct, doseGrid, rtdoseInfo);
                
                croppedResampledData.totalDose = totalDoseCropped;
                croppedResampledData.cropResampleInfo = cropResampleInfo;
                
                fprintf('    ✓ totalDose cropped/resampled: %s -> %s\n', ...
                    mat2str(size(totalDose)), mat2str(size(totalDoseCropped)));
                fprintf('      Max dose after crop/resample: %.4f Gy\n', max(totalDoseCropped(:)));
                fprintf('      Crop region (CT indices): rows [%d:%d], cols [%d:%d], slices [%d:%d]\n', ...
                    cropResampleInfo.cropIndices.rowStart, cropResampleInfo.cropIndices.rowEnd, ...
                    cropResampleInfo.cropIndices.colStart, cropResampleInfo.cropIndices.colEnd, ...
                    cropResampleInfo.cropIndices.sliceStart, cropResampleInfo.cropIndices.sliceEnd);
                
                % Crop and resample individual field doses
                fprintf('    Cropping and resampling individual field doses...\n');
                croppedResampledData.fieldDoses = cell(length(fieldDoses), 1);
                
                for beamIdx = 1:length(fieldDoses)
                    if ~isempty(fieldDoses{beamIdx})
                        [fieldCropped, ~] = crop_and_resample_to_rtdose(...
                            fieldDoses{beamIdx}.physicalDose, ct, doseGrid, rtdoseInfo);
                        
                        croppedResampledData.fieldDoses{beamIdx} = fieldDoses{beamIdx};
                        croppedResampledData.fieldDoses{beamIdx}.physicalDose = fieldCropped;
                        croppedResampledData.fieldDoses{beamIdx}.croppedResampledMaxDose = max(fieldCropped(:));
                        
                        fprintf('      Field %2d cropped/resampled: max=%.4f Gy\n', beamIdx, max(fieldCropped(:)));
                    end
                end
                
                % Crop and resample CT to RTDOSE grid
                fprintf('    Cropping and resampling CT to RTDOSE grid...\n');
                [ctCropped, ~] = crop_and_resample_to_rtdose(...
                    ct.cubeHU{1}, ct, doseGrid, rtdoseInfo);
                croppedResampledData.ct = ctCropped;
                fprintf('    ✓ CT cropped/resampled: %s -> %s\n', ...
                    mat2str(size(ct.cubeHU{1})), mat2str(size(ctCropped)));
                
                % Store reference dose for convenience
                croppedResampledData.referenceDose = referenceDoseOriginal;
                
                % Save cropped/resampled data
                save(fullfile(outputPath, 'croppedResampledData.mat'), 'croppedResampledData', '-v7.3');
                fprintf('  ✓ Cropped/resampled data saved to: croppedResampledData.mat\n');
                
            catch ME
                fprintf('    ⚠ Error during crop/resample: %s\n', ME.message);
                fprintf('    Stack trace:\n');
                for k = 1:length(ME.stack)
                    fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                end
                fprintf('    Continuing without cropped/resampled data...\n');
            end
        else
            fprintf('    ⚠ Cannot crop/resample: missing reference RTDOSE or totalDose\n');
        end
        
        % Compare with reference RTDOSE
        if ~isempty(referenceDoseOriginal)
            fprintf('\n  Comparing with reference RTDOSE...\n');
            
            comparison = struct();
            comparison.reference = referenceDoseOriginal;
            comparison.referenceGrid = doseGrid;
            
            if exist('totalDoseCropped', 'var') && ~isempty(totalDoseCropped)
                % Compare cropped/resampled calculated dose with reference
                if isequal(size(totalDoseCropped), size(referenceDoseOriginal))
                    doseDiff = totalDoseCropped - referenceDoseOriginal;
                    
                    fprintf('    Comparison statistics (on RTDOSE grid):\n');
                    fprintf('      Calculated max dose: %.4f Gy\n', max(totalDoseCropped(:)));
                    fprintf('      Reference max dose:  %.4f Gy\n', max(referenceDoseOriginal(:)));
                    fprintf('      Mean absolute difference: %.4f Gy\n', mean(abs(doseDiff(:))));
                    fprintf('      Max difference: %.4f Gy\n', max(abs(doseDiff(:))));
                    fprintf('      RMS difference: %.4f Gy\n', sqrt(mean(doseDiff(:).^2)));
                    
                    % Calculate relative differences in high dose region
                    highDoseThreshold = 0.5 * max(referenceDoseOriginal(:));
                    highDoseMask = referenceDoseOriginal > highDoseThreshold;
                    if any(highDoseMask(:))
                        relativeDiff = abs(doseDiff(highDoseMask)) ./ referenceDoseOriginal(highDoseMask) * 100;
                        fprintf('      Mean relative diff (>50%% max): %.2f%%\n', mean(relativeDiff));
                        comparison.metrics.meanRelativeDiff = mean(relativeDiff);
                    end
                    
                    comparison.calculated = totalDoseCropped;
                    comparison.difference = doseDiff;
                    comparison.metrics.meanAbsDiff = mean(abs(doseDiff(:)));
                    comparison.metrics.maxDiff = max(abs(doseDiff(:)));
                    comparison.metrics.rmsDiff = sqrt(mean(doseDiff(:).^2));
                    comparison.note = 'Comparison on RTDOSE grid (CT cropped to RTDOSE extent, then resampled)';
                    
                    fprintf('    ✓ Comparison complete\n');
                else
                    fprintf('    ⚠ Size mismatch after crop/resample:\n');
                    fprintf('      Calculated: %s\n', mat2str(size(totalDoseCropped)));
                    fprintf('      Reference: %s\n', mat2str(size(referenceDoseOriginal)));
                    comparison.calculated = totalDoseCropped;
                    comparison.note = 'Size mismatch - direct comparison not possible';
                end
            else
                fprintf('    ⚠ No cropped/resampled dose available for comparison\n');
                comparison.calculated = totalDose;
                comparison.calculatedGrid = calculatedGridSize;
                comparison.note = 'Different grids - data saved separately';
            end
            
            save(fullfile(outputPath, 'doseComparison.mat'), 'comparison', '-v7.3');
            fprintf('  ✓ Comparison saved to: doseComparison.mat\n');
        end
        
        %% ================================================================
        %  PROCESSING COMPLETE - SUMMARY
        %  ================================================================
        fprintf('\n============================================================\n');
        fprintf('  Processing complete for %s - %s\n', currentID, currentSession);
        fprintf('============================================================\n');
        fprintf('  Results saved to: %s\n\n', outputPath);
        
        fprintf('  ** IMPORTANT NOTES **\n\n');
        
        fprintf('  1. MLC Aperture Status:\n');
        if isfield(pln, 'propStf') && isfield(pln.propStf, 'collimation')
            if isfield(pln.propStf.collimation, 'dualLayerReduced') && pln.propStf.collimation.dualLayerReduced
                fprintf('     ✓ Dual-layer MLC successfully reduced to single-layer\n');
            else
                fprintf('     - Single-layer MLC used (no reduction needed)\n');
            end
        else
            fprintf('     ⚠ MLC data may not have been properly applied\n');
        end
        
        fprintf('\n  2. Weight Extraction:\n');
        if isfield(pln, 'w') && ~isempty(pln.w)
            fprintf('     ✓ Weights extracted from RTPLAN: %d values, sum=%.2f\n', ...
                length(pln.w), sum(pln.w));
        else
            fprintf('     ⚠ Uniform weights used\n');
        end
        
        fprintf('\n  3. Machine Commissioning:\n');
        fprintf('     - Using: %s\n', pln.machine);
        fprintf('     - Note: Generic machine data (not ETHOS/Halcyon-specific)\n');
        fprintf('     - Absolute dose values will differ from clinical plan\n');
        
        fprintf('\n  4. Output Files:\n');
        fprintf('     - fieldDoses.mat: All field doses + metadata (CT grid)\n');
        fprintf('     - reconstructedDose.mat: Sum of all fields (CT grid)\n');
        fprintf('     - doseComparison.mat: Comparison metrics\n');
        fprintf('     - croppedResampledData.mat: CT & doses cropped to RTDOSE extent\n');
        fprintf('     - Field_XX.mat: Individual field files\n');
        
        fprintf('\n  5. Crop/Resample Approach:\n');
        fprintf('     - CT cropped to RTDOSE spatial extent\n');
        fprintf('     - Cropped region resampled to match RTDOSE dimensions\n');
        fprintf('     - Enables direct voxel-by-voxel comparison\n');
        
        fprintf('\n  6. Validation Notes:\n');
        fprintf('     - %d/%d fields calculated successfully\n', numSuccessful, length(stf));
        if ~isempty(totalDose)
            fprintf('     - Total reconstructed max dose: %.4f Gy\n', max(totalDose(:)));
        end
        if ~isempty(referenceDoseOriginal)
            fprintf('     - Reference RTDOSE max: %.4f Gy\n', max(referenceDoseOriginal(:)));
        end
        fprintf('\n');
        
    end
end

fprintf('\n============================================================\n');
fprintf('  All processing complete!\n');
fprintf('============================================================\n');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function pln = extract_mlc_basic(pln, rtplanInfo)
%EXTRACT_MLC_BASIC Basic MLC extraction without dual-layer reduction
%   Fallback function when reduce_collimator.m is not available

    fprintf('  Performing basic MLC extraction...\n');
    
    if ~isfield(rtplanInfo, 'BeamSequence')
        return;
    end
    
    % Initialize structures
    if ~isfield(pln, 'propStf')
        pln.propStf = struct();
    end
    if ~isfield(pln.propStf, 'beam')
        pln.propStf.beam = struct();
    end
    
    numBeams = length(fieldnames(rtplanInfo.BeamSequence));
    
    for beamIdx = 1:numBeams
        beamField = sprintf('Item_%d', beamIdx);
        beam = rtplanInfo.BeamSequence.(beamField);
        
        if ~isfield(beam, 'ControlPointSequence')
            continue;
        end
        
        % Use first control point
        cp = beam.ControlPointSequence.Item_1;
        
        if ~isfield(cp, 'BeamLimitingDevicePositionSequence')
            continue;
        end
        
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
            
            if contains(deviceType, 'ASYMX', 'IgnoreCase', true) || strcmp(deviceType, 'X')
                pln.propStf.beam(beamIdx).jaw.x = positions;
            elseif contains(deviceType, 'ASYMY', 'IgnoreCase', true) || strcmp(deviceType, 'Y')
                pln.propStf.beam(beamIdx).jaw.y = positions;
            elseif contains(deviceType, 'MLC', 'IgnoreCase', true)
                pln.propStf.beam(beamIdx).shape.x = positions;
            end
        end
    end
    
    pln.propStf.collimation.type = 'mlc';
    pln.propStf.collimation.dualLayerReduced = false;
    
    fprintf('  ✓ Basic MLC extraction complete for %d beams\n', numBeams);
end

function stf = inject_mlc_to_stf(stf, pln)
%INJECT_MLC_TO_STF Inject MLC shapes from pln into stf structure

    if ~isfield(pln, 'propStf') || ~isfield(pln.propStf, 'beam')
        fprintf('    ⚠ No MLC data in pln.propStf.beam\n');
        return;
    end
    
    shapesInjected = 0;
    
    for beamIdx = 1:length(stf)
        if beamIdx > length(pln.propStf.beam)
            continue;
        end
        
        % Inject shape at beam level
        if isfield(pln.propStf.beam(beamIdx), 'shape')
            if isfield(pln.propStf.beam(beamIdx).shape, 'x')
                stf(beamIdx).shape.x = pln.propStf.beam(beamIdx).shape.x;
                
                % Also inject to each ray
                for rayIdx = 1:length(stf(beamIdx).ray)
                    stf(beamIdx).ray(rayIdx).shape.x = pln.propStf.beam(beamIdx).shape.x;
                end
                
                shapesInjected = shapesInjected + 1;
            end
        end
        
        % Inject jaw positions
        if isfield(pln.propStf.beam(beamIdx), 'jaw')
            stf(beamIdx).jaw = pln.propStf.beam(beamIdx).jaw;
        end
    end
    
    if shapesInjected > 0
        fprintf('    ✓ MLC shapes injected into %d beams\n', shapesInjected);
    else
        fprintf('    ⚠ No MLC shapes were injected\n');
    end
end

function pln = populate_weights(pln, stf, beamMetersets)
%POPULATE_WEIGHTS Create pln.w weight vector from beam metersets

    % Calculate total number of bixels
    totalBixels = 0;
    bixelsPerBeam = zeros(length(stf), 1);
    
    for beamIdx = 1:length(stf)
        numRays = stf(beamIdx).numOfRays;
        if isfield(stf(beamIdx).ray(1), 'energy')
            numEnergies = length(stf(beamIdx).ray(1).energy);
        else
            numEnergies = 1;
        end
        bixelsPerBeam(beamIdx) = numRays * numEnergies;
        totalBixels = totalBixels + bixelsPerBeam(beamIdx);
    end
    
    fprintf('    Total bixels across all beams: %d\n', totalBixels);
    
    % Initialize weight vector
    allWeights = zeros(totalBixels, 1);
    
    if ~isempty(beamMetersets) && length(beamMetersets) >= length(stf)
        % Distribute beam metersets uniformly across bixels
        offset = 0;
        for beamIdx = 1:length(stf)
            numBixels = bixelsPerBeam(beamIdx);
            
            % Distribute this beam's meterset uniformly
            weightPerBixel = beamMetersets(beamIdx) / numBixels;
            allWeights(offset+1:offset+numBixels) = weightPerBixel;
            
            fprintf('    Beam %2d: %d bixels, meterset=%.2f, weight/bixel=%.6f\n', ...
                beamIdx, numBixels, beamMetersets(beamIdx), weightPerBixel);
            
            offset = offset + numBixels;
        end
        
        pln.w = allWeights;
        fprintf('    ✓ Weight vector created: %d values, sum=%.2f\n', length(pln.w), sum(pln.w));
    else
        % Use uniform weights
        pln.w = ones(totalBixels, 1);
        fprintf('    ⚠ Using uniform weights: %d values\n', totalBixels);
    end
end

function [croppedResampled, info] = crop_and_resample_to_rtdose(data, ct, doseGrid, rtdoseInfo)
%CROP_AND_RESAMPLE_TO_RTDOSE Crop data to RTDOSE extent and resample to match dimensions
%
%   This function:
%   1. Computes the physical extent of the RTDOSE grid in patient coordinates
%   2. Finds the corresponding region in the CT coordinate system
%   3. Crops the input data (CT or dose on CT grid) to that region
%   4. Resamples the cropped region to match RTDOSE dimensions exactly
%
%   Inputs:
%       data        - 3D array on CT grid (dose or HU values)
%       ct          - matRad CT structure with resolution and coordinate info
%       doseGrid    - Structure with RTDOSE grid parameters
%       rtdoseInfo  - DICOM info structure from RTDOSE file
%
%   Outputs:
%       croppedResampled - 3D array matching RTDOSE dimensions
%       info             - Structure with crop/resample metadata

    info = struct();
    info.originalSize = size(data);
    info.method = 'crop_then_trilinear';
    
    % Get CT grid parameters
    ctRows = size(data, 1);
    ctCols = size(data, 2);
    ctSlices = size(data, 3);
    
    ctResX = ct.resolution.x;  % Column spacing (mm)
    ctResY = ct.resolution.y;  % Row spacing (mm)
    ctResZ = ct.resolution.z;  % Slice spacing (mm)
    
    % Get CT origin from dicomInfo if available, otherwise estimate from matRad structure
    if isfield(ct, 'dicomInfo') && isfield(ct.dicomInfo, 'ImagePositionPatient')
        ctOrigin = ct.dicomInfo.ImagePositionPatient;
    elseif isfield(ct, 'x') && isfield(ct, 'y') && isfield(ct, 'z')
        % matRad stores coordinate vectors
        ctOrigin = [ct.x(1), ct.y(1), ct.z(1)];
    else
        % Fallback: assume origin at (0,0,0) - this may need adjustment
        ctOrigin = [0, 0, 0];
        warning('CT origin not found, using [0,0,0]. Results may be misaligned.');
    end
    
    % Get RTDOSE grid parameters
    doseRows = rtdoseInfo.Rows;
    doseCols = rtdoseInfo.Columns;
    if isfield(rtdoseInfo, 'GridFrameOffsetVector')
        doseSlices = length(rtdoseInfo.GridFrameOffsetVector);
    else
        doseSlices = doseGrid.dimensions(3);
    end
    
    doseResX = doseGrid.resolution(2);  % Column spacing
    doseResY = doseGrid.resolution(1);  % Row spacing
    doseResZ = doseGrid.resolution(3);  % Slice spacing
    
    doseOrigin = doseGrid.origin;  % [x, y, z] of first voxel center
    
    % Compute RTDOSE physical extent (min/max coordinates)
    % DICOM convention: ImagePositionPatient is center of first voxel
    % X increases with column index, Y increases with row index
    doseXmin = doseOrigin(1);
    doseXmax = doseOrigin(1) + (doseCols - 1) * doseResX;
    doseYmin = doseOrigin(2);
    doseYmax = doseOrigin(2) + (doseRows - 1) * doseResY;
    
    % Z coordinates from GridFrameOffsetVector
    if isfield(rtdoseInfo, 'GridFrameOffsetVector')
        doseZmin = doseOrigin(3) + min(rtdoseInfo.GridFrameOffsetVector);
        doseZmax = doseOrigin(3) + max(rtdoseInfo.GridFrameOffsetVector);
    else
        doseZmin = doseOrigin(3);
        doseZmax = doseOrigin(3) + (doseSlices - 1) * doseResZ;
    end
    
    info.doseExtent.xRange = [doseXmin, doseXmax];
    info.doseExtent.yRange = [doseYmin, doseYmax];
    info.doseExtent.zRange = [doseZmin, doseZmax];
    
    % Compute CT physical coordinates
    % Note: In matRad, the CT cube is typically stored as (rows, cols, slices)
    % where rows correspond to Y, cols to X, slices to Z
    ctXmin = ctOrigin(1);
    ctXmax = ctOrigin(1) + (ctCols - 1) * ctResX;
    ctYmin = ctOrigin(2);
    ctYmax = ctOrigin(2) + (ctRows - 1) * ctResY;
    ctZmin = ctOrigin(3);
    ctZmax = ctOrigin(3) + (ctSlices - 1) * ctResZ;
    
    info.ctExtent.xRange = [ctXmin, ctXmax];
    info.ctExtent.yRange = [ctYmin, ctYmax];
    info.ctExtent.zRange = [ctZmin, ctZmax];
    
    % Find CT indices corresponding to RTDOSE extent
    % Add small margin to ensure we capture the full RTDOSE region
    margin = 0.5;  % Half voxel margin
    
    % Column indices (X direction)
    colStart = max(1, floor((doseXmin - ctXmin) / ctResX + 1 - margin));
    colEnd = min(ctCols, ceil((doseXmax - ctXmin) / ctResX + 1 + margin));
    
    % Row indices (Y direction)
    rowStart = max(1, floor((doseYmin - ctYmin) / ctResY + 1 - margin));
    rowEnd = min(ctRows, ceil((doseYmax - ctYmin) / ctResY + 1 + margin));
    
    % Slice indices (Z direction)
    sliceStart = max(1, floor((doseZmin - ctZmin) / ctResZ + 1 - margin));
    sliceEnd = min(ctSlices, ceil((doseZmax - ctZmin) / ctResZ + 1 + margin));
    
    info.cropIndices.rowStart = rowStart;
    info.cropIndices.rowEnd = rowEnd;
    info.cropIndices.colStart = colStart;
    info.cropIndices.colEnd = colEnd;
    info.cropIndices.sliceStart = sliceStart;
    info.cropIndices.sliceEnd = sliceEnd;
    
    % Crop the data
    croppedData = data(rowStart:rowEnd, colStart:colEnd, sliceStart:sliceEnd);
    info.croppedSize = size(croppedData);
    
    % Compute the physical coordinates of the cropped region
    croppedXmin = ctXmin + (colStart - 1) * ctResX;
    croppedYmin = ctYmin + (rowStart - 1) * ctResY;
    croppedZmin = ctZmin + (sliceStart - 1) * ctResZ;
    
    info.croppedOrigin = [croppedXmin, croppedYmin, croppedZmin];
    
    % Now resample from cropped CT grid to RTDOSE grid
    % Create coordinate vectors for cropped data
    croppedRows = size(croppedData, 1);
    croppedCols = size(croppedData, 2);
    croppedSlices = size(croppedData, 3);
    
    % Physical coordinates of cropped voxel centers
    croppedX = croppedXmin + (0:croppedCols-1) * ctResX;
    croppedY = croppedYmin + (0:croppedRows-1) * ctResY;
    croppedZ = croppedZmin + (0:croppedSlices-1) * ctResZ;
    
    % Physical coordinates of RTDOSE voxel centers
    doseX = doseOrigin(1) + (0:doseCols-1) * doseResX;
    doseY = doseOrigin(2) + (0:doseRows-1) * doseResY;
    
    if isfield(rtdoseInfo, 'GridFrameOffsetVector')
        doseZ = doseOrigin(3) + rtdoseInfo.GridFrameOffsetVector(:)';
    else
        doseZ = doseOrigin(3) + (0:doseSlices-1) * doseResZ;
    end
    
    % Create meshgrids for interpolation
    % Note: interpn expects dimensions in order (Y, X, Z) -> (row, col, slice)
    [croppedYgrid, croppedXgrid, croppedZgrid] = ndgrid(croppedY, croppedX, croppedZ);
    [doseYgrid, doseXgrid, doseZgrid] = ndgrid(doseY, doseX, doseZ);
    
    % Perform trilinear interpolation
    croppedResampled = interpn(croppedYgrid, croppedXgrid, croppedZgrid, ...
                               double(croppedData), ...
                               doseYgrid, doseXgrid, doseZgrid, ...
                               'linear', 0);  % 0 for extrapolation
    
    info.targetSize = [doseRows, doseCols, doseSlices];
    info.finalSize = size(croppedResampled);
    info.maxOriginal = max(data(:));
    info.maxCropped = max(croppedData(:));
    info.maxFinal = max(croppedResampled(:));
    
    % Verify output dimensions match RTDOSE
    if ~isequal(size(croppedResampled), [doseRows, doseCols, doseSlices])
        warning('Output size %s does not match expected RTDOSE size %s', ...
            mat2str(size(croppedResampled)), mat2str([doseRows, doseCols, doseSlices]));
    end
end