%% ETHOS RTPLAN Weight Diagnostics
% Purpose: Diagnose why beam weights are not being imported into matrad pln structure
% This script performs extensive investigation of RTPLAN DICOM structure
% Author: Diagnostic tool for weight import issues
% Date: 2025

clear; clc; close all;

%% Configuration
% Use same paths as main script
id = '1194203';
session = 'Session_1';

wd = '/mnt/weka/home/80030361/ETHOS_Simulations';
matradPath = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

rawwd = fullfile(wd, 'EthosExports', id, 'Pancreas', session);
dicomPath = fullfile(rawwd, 'sct');

% Add matrad to path
addpath(genpath(matradPath));

fprintf('========================================\n');
fprintf('ETHOS RTPLAN Weight Diagnostics\n');
fprintf('========================================\n\n');

%% Step 1: Locate all DICOM files
fprintf('[1/5] Locating DICOM files...\n');

% Find all DICOM files
allDicoms = dir(fullfile(dicomPath, '*.dcm'));
fprintf('  Found %d DICOM files total\n', length(allDicoms));

% Categorize by modality
rtplanFiles = {};
rtdoseFiles = {};
rtstructFiles = {};
ctFiles = {};

for i = 1:length(allDicoms)
    try
        info = dicominfo(fullfile(allDicoms(i).folder, allDicoms(i).name));
        modality = info.Modality;
        
        switch modality
            case 'RTPLAN'
                rtplanFiles{end+1} = fullfile(allDicoms(i).folder, allDicoms(i).name);
            case 'RTDOSE'
                rtdoseFiles{end+1} = fullfile(allDicoms(i).folder, allDicoms(i).name);
            case 'RTSTRUCT'
                rtstructFiles{end+1} = fullfile(allDicoms(i).folder, allDicoms(i).name);
            case 'CT'
                ctFiles{end+1} = fullfile(allDicoms(i).folder, allDicoms(i).name);
        end
    catch
        % Skip non-standard DICOM files
    end
end

fprintf('  - RTPLAN files: %d\n', length(rtplanFiles));
fprintf('  - RTDOSE files: %d\n', length(rtdoseFiles));
fprintf('  - RTSTRUCT files: %d\n', length(rtstructFiles));
fprintf('  - CT files: %d\n', length(ctFiles));

if isempty(rtplanFiles)
    error('No RTPLAN files found! Cannot proceed with weight diagnostics.');
end

%% Step 2: Deep dive into RTPLAN structure
fprintf('\n[2/5] Analyzing RTPLAN DICOM structure...\n\n');

for planIdx = 1:length(rtplanFiles)
    fprintf('================== RTPLAN #%d ==================\n', planIdx);
    rtplanFile = rtplanFiles{planIdx};
    [~, filename, ~] = fileparts(rtplanFile);
    fprintf('File: %s\n\n', filename);
    
    % Read RTPLAN metadata
    rtplan = dicominfo(rtplanFile);
    
    % Display high-level plan information
    fprintf('=== PLAN OVERVIEW ===\n');
    if isfield(rtplan, 'RTPlanLabel')
        fprintf('  Plan Label: %s\n', rtplan.RTPlanLabel);
    end
    if isfield(rtplan, 'RTPlanName')
        fprintf('  Plan Name: %s\n', rtplan.RTPlanName);
    end
    if isfield(rtplan, 'RTPlanDescription')
        fprintf('  Description: %s\n', rtplan.RTPlanDescription);
    end
    if isfield(rtplan, 'RTPlanGeometry')
        fprintf('  Geometry: %s\n', rtplan.RTPlanGeometry);
    end
    
    %% Check for Fraction Group Sequence (PRIMARY LOCATION FOR WEIGHTS)
    fprintf('\n=== FRACTION GROUP SEQUENCE ===\n');
    if isfield(rtplan, 'FractionGroupSequence')
        numFractionGroups = length(rtplan.FractionGroupSequence);
        fprintf('  Found %d Fraction Group(s)\n', numFractionGroups);
        
        for fgIdx = 1:numFractionGroups
            fg = rtplan.FractionGroupSequence(fgIdx);
            fprintf('\n  --- Fraction Group %d ---\n', fgIdx);
            
            if isfield(fg, 'FractionGroupNumber')
                fprintf('    Fraction Group Number: %d\n', fg.FractionGroupNumber);
            end
            if isfield(fg, 'NumberOfFractionsPlanned')
                fprintf('    Number of Fractions: %d\n', fg.NumberOfFractionsPlanned);
            end
            if isfield(fg, 'NumberOfBeams')
                fprintf('    Number of Beams: %d\n', fg.NumberOfBeams);
            end
            
            % Referenced Beam Sequence - THIS CONTAINS METERSET WEIGHTS
            if isfield(fg, 'ReferencedBeamSequence')
                numRefBeams = length(fg.ReferencedBeamSequence);
                fprintf('    Referenced Beam Sequence found: %d beams\n\n', numRefBeams);
                
                totalMeterset = 0;
                for rbIdx = 1:numRefBeams
                    rb = fg.ReferencedBeamSequence(rbIdx);
                    fprintf('    ** BEAM %d **\n', rbIdx);
                    
                    if isfield(rb, 'ReferencedBeamNumber')
                        fprintf('      Beam Number: %d\n', rb.ReferencedBeamNumber);
                    end
                    
                    % CRITICAL: BeamMeterset field
                    if isfield(rb, 'BeamMeterset')
                        meterset = rb.BeamMeterset;
                        totalMeterset = totalMeterset + meterset;
                        fprintf('      >>> BeamMeterset: %.6f MU <<<\n', meterset);
                    else
                        fprintf('      BeamMeterset: NOT FOUND\n');
                    end
                    
                    % Alternative: BeamDose
                    if isfield(rb, 'BeamDose')
                        fprintf('      BeamDose: %.6f Gy\n', rb.BeamDose);
                    end
                    
                    % Check for control point sequence in referenced beam
                    if isfield(rb, 'ReferencedControlPointSequence')
                        fprintf('      Referenced Control Points: %d\n', length(rb.ReferencedControlPointSequence));
                    end
                    
                    fprintf('\n');
                end
                
                fprintf('    Total Meterset across all beams: %.2f MU\n', totalMeterset);
            else
                fprintf('    WARNING: No ReferencedBeamSequence found!\n');
            end
        end
    else
        fprintf('  WARNING: No FractionGroupSequence found in RTPLAN!\n');
    end
    
    %% Check Beam Sequence (BEAM GEOMETRY AND PARAMETERS)
    fprintf('\n=== BEAM SEQUENCE ===\n');
    if isfield(rtplan, 'BeamSequence')
        numBeams = length(rtplan.BeamSequence);
        fprintf('  Found %d Beam(s) in BeamSequence\n', numBeams);
        
        for beamIdx = 1:numBeams
            beam = rtplan.BeamSequence(beamIdx);
            fprintf('\n  --- Beam %d ---\n', beamIdx);
            
            if isfield(beam, 'BeamNumber')
                fprintf('    Beam Number: %d\n', beam.BeamNumber);
            end
            if isfield(beam, 'BeamName')
                fprintf('    Beam Name: %s\n', beam.BeamName);
            end
            if isfield(beam, 'BeamType')
                fprintf('    Beam Type: %s\n', beam.BeamType);
            end
            if isfield(beam, 'RadiationType')
                fprintf('    Radiation Type: %s\n', beam.RadiationType);
            end
            if isfield(beam, 'TreatmentMachineName')
                fprintf('    Machine: %s\n', beam.TreatmentMachineName);
            end
            
            % Beam geometry
            if isfield(beam, 'ControlPointSequence')
                numCP = length(beam.ControlPointSequence);
                fprintf('    Control Points: %d\n', numCP);
                
                % Check first and last control points for meterset
                if numCP > 0
                    cp1 = beam.ControlPointSequence(1);
                    fprintf('    \n    Control Point 1:\n');
                    if isfield(cp1, 'ControlPointIndex')
                        fprintf('      Index: %d\n', cp1.ControlPointIndex);
                    end
                    if isfield(cp1, 'CumulativeMetersetWeight')
                        fprintf('      Cumulative Meterset Weight: %.6f\n', cp1.CumulativeMetersetWeight);
                    end
                    if isfield(cp1, 'GantryAngle')
                        fprintf('      Gantry Angle: %.2f°\n', cp1.GantryAngle);
                    end
                    if isfield(cp1, 'BeamLimitingDeviceAngle')
                        fprintf('      Collimator Angle: %.2f°\n', cp1.BeamLimitingDeviceAngle);
                    end
                    if isfield(cp1, 'PatientSupportAngle')
                        fprintf('      Couch Angle: %.2f°\n', cp1.PatientSupportAngle);
                    end
                    
                    % Check last control point
                    if numCP > 1
                        cpLast = beam.ControlPointSequence(numCP);
                        fprintf('    \n    Control Point %d (Last):\n', numCP);
                        if isfield(cpLast, 'CumulativeMetersetWeight')
                            fprintf('      Cumulative Meterset Weight: %.6f\n', cpLast.CumulativeMetersetWeight);
                            fprintf('      >>> This should equal 1.0 for normalized plan <<<\n');
                        end
                    end
                end
            else
                fprintf('    WARNING: No ControlPointSequence found!\n');
            end
            
            % Check for final cumulative meterset weight
            if isfield(beam, 'FinalCumulativeMetersetWeight')
                fprintf('    Final Cumulative Meterset Weight: %.6f\n', beam.FinalCumulativeMetersetWeight);
            end
            
            % Check for number of control points
            if isfield(beam, 'NumberOfControlPoints')
                fprintf('    Number of Control Points: %d\n', beam.NumberOfControlPoints);
            end
        end
    else
        fprintf('  WARNING: No BeamSequence found in RTPLAN!\n');
    end
    
    fprintf('\n');
end

%% Step 3: Import with matrad and check what gets populated
fprintf('\n[3/5] Importing with matRad_DicomImporter...\n');

try
    importer = matRad_DicomImporter(dicomPath);
    importer.matRad_importDicom();
    fprintf('  ✓ Import successful\n');
catch ME
    fprintf('  ✗ Import failed: %s\n', ME.message);
    return;
end

%% Step 4: Examine matrad pln structure for weights
fprintf('\n[4/5] Examining matRad plan structure...\n\n');

fprintf('=== PLN STRUCTURE ANALYSIS ===\n');

% Check top-level fields
fprintf('Top-level pln fields:\n');
plnFields = fieldnames(pln);
for i = 1:length(plnFields)
    fprintf('  - pln.%s\n', plnFields{i});
end

% Check for weight vector
fprintf('\n--- Weight Vector ---\n');
if isfield(pln, 'w')
    fprintf('  ✓ pln.w EXISTS\n');
    fprintf('    Size: %d elements\n', length(pln.w));
    fprintf('    Sum: %.6f\n', sum(pln.w));
    fprintf('    Min: %.6f\n', min(pln.w));
    fprintf('    Max: %.6f\n', max(pln.w));
    fprintf('    Non-zero elements: %d\n', nnz(pln.w));
    
    % Show first 20 values
    fprintf('    First 20 values:\n');
    for i = 1:min(20, length(pln.w))
        fprintf('      w(%d) = %.6f\n', i, pln.w(i));
    end
else
    fprintf('  ✗ pln.w DOES NOT EXIST\n');
end

% Check propStf structure
fprintf('\n--- propStf Structure ---\n');
if isfield(pln, 'propStf')
    fprintf('  ✓ pln.propStf EXISTS\n');
    propStfFields = fieldnames(pln.propStf);
    fprintf('    Fields:\n');
    for i = 1:length(propStfFields)
        fprintf('      - pln.propStf.%s\n', propStfFields{i});
    end
    
    % Check for beam array
    if isfield(pln.propStf, 'beam')
        fprintf('\n    ✓ pln.propStf.beam EXISTS\n');
        fprintf('      Number of beams: %d\n', length(pln.propStf.beam));
        
        % Examine first beam
        beamFields = fieldnames(pln.propStf.beam);
        fprintf('      Fields in beam structure:\n');
        for i = 1:length(beamFields)
            fprintf('        - pln.propStf.beam.%s\n', beamFields{i});
        end
        
        % Check for weight field
        if isfield(pln.propStf.beam, 'weight')
            fprintf('\n      ✓ pln.propStf.beam.weight EXISTS\n');
            for i = 1:length(pln.propStf.beam)
                fprintf('        Beam %d weight: %.6f\n', i, pln.propStf.beam(i).weight);
            end
        else
            fprintf('\n      ✗ pln.propStf.beam.weight DOES NOT EXIST\n');
        end
        
        % Check for meterset field
        if isfield(pln.propStf.beam, 'meterset')
            fprintf('\n      ✓ pln.propStf.beam.meterset EXISTS\n');
            for i = 1:length(pln.propStf.beam)
                fprintf('        Beam %d meterset: %.6f\n', i, pln.propStf.beam(i).meterset);
            end
        else
            fprintf('      ✗ pln.propStf.beam.meterset DOES NOT EXIST\n');
        end
    else
        fprintf('    ✗ pln.propStf.beam DOES NOT EXIST\n');
    end
else
    fprintf('  ✗ pln.propStf DOES NOT EXIST\n');
end

% Check for any meterset-related fields
fprintf('\n--- Searching for meterset/weight keywords in pln ---\n');
searchKeywords = {'weight', 'meterset', 'MU', 'dose'};
for kw = 1:length(searchKeywords)
    fprintf('  Searching for "%s"...\n', searchKeywords{kw});
    found = false;
    for i = 1:length(plnFields)
        if contains(lower(plnFields{i}), lower(searchKeywords{kw}))
            fprintf('    ✓ Found: pln.%s\n', plnFields{i});
            found = true;
        end
    end
    if ~found
        fprintf('    ✗ No matches found\n');
    end
end

%% Step 5: Generate steering file and compare
fprintf('\n[5/5] Generating steering file and checking for weights...\n');

try
    stf = matRad_generateStf(ct, cst, pln);
    fprintf('  ✓ Steering file generated\n');
    fprintf('    Number of beams: %d\n', length(stf));
    
    % Check stf structure for weights
    fprintf('\n  --- STF Structure Analysis ---\n');
    stfFields = fieldnames(stf);
    fprintf('  Top-level stf fields:\n');
    for i = 1:length(stfFields)
        fprintf('    - stf.%s\n', stfFields{i});
    end
    
    % Check first beam
    if isfield(stf, 'ray')
        fprintf('\n  Examining stf(1).ray structure:\n');
        rayFields = fieldnames(stf(1).ray);
        for i = 1:length(rayFields)
            fprintf('    - stf(1).ray.%s\n', rayFields{i});
        end
        
        % Check for weights in ray structure
        if isfield(stf(1).ray, 'weight')
            fprintf('\n  ✓ stf.ray.weight EXISTS\n');
            fprintf('    First 5 rays:\n');
            for i = 1:min(5, length(stf(1).ray))
                fprintf('      Ray %d weight: %.6f\n', i, stf(1).ray(i).weight);
            end
        else
            fprintf('\n  ✗ stf.ray.weight DOES NOT EXIST\n');
        end
        
        % Check energy structure
        if isfield(stf(1).ray, 'energy')
            fprintf('\n  Energy levels per ray: %d\n', length(stf(1).ray(1).energy));
            if length(stf(1).ray(1).energy) > 0
                energyFields = fieldnames(stf(1).ray(1).energy);
                fprintf('  Energy structure fields:\n');
                for i = 1:length(energyFields)
                    fprintf('    - stf(1).ray(1).energy.%s\n', energyFields{i});
                end
                
                if isfield(stf(1).ray(1).energy, 'weight')
                    fprintf('\n  ✓ stf.ray.energy.weight EXISTS\n');
                    fprintf('    Example: Ray 1, Energy 1 weight: %.6f\n', stf(1).ray(1).energy(1).weight);
                else
                    fprintf('\n  ✗ stf.ray.energy.weight DOES NOT EXIST\n');
                end
            end
        end
    end
    
catch ME
    fprintf('  ✗ Error generating steering file: %s\n', ME.message);
end

%% Summary and Recommendations
fprintf('\n========================================\n');
fprintf('DIAGNOSTIC SUMMARY\n');
fprintf('========================================\n\n');

fprintf('KEY FINDINGS:\n\n');

% Check if metersets were found in DICOM
metersetFound = false;
if exist('totalMeterset', 'var') && totalMeterset > 0
    metersetFound = true;
end

if metersetFound
    fprintf('✓ Meterset values found in RTPLAN DICOM files\n');
    fprintf('  Total meterset: %.2f MU\n', totalMeterset);
else
    fprintf('✗ No meterset values found in RTPLAN DICOM files\n');
end

% Check if weights were imported to pln
if isfield(pln, 'w')
    fprintf('✓ Weight vector exists in pln.w\n');
    fprintf('  Total weight: %.6f\n', sum(pln.w));
else
    fprintf('✗ No weight vector in pln.w\n');
end

fprintf('\nPOSSIBLE ISSUES:\n\n');

if ~metersetFound
    fprintf('1. RTPLAN files may not contain meterset information\n');
    fprintf('   - Check if this is a planning file vs. delivery record\n');
    fprintf('   - ETHOS may export plans without meterset weights\n\n');
end

if metersetFound && ~isfield(pln, 'w')
    fprintf('2. matRad importer may not be reading meterset values\n');
    fprintf('   - Check matRad_DicomImporter version/configuration\n');
    fprintf('   - May need to manually extract from DICOM and populate pln.w\n\n');
end

fprintf('RECOMMENDED ACTIONS:\n\n');
fprintf('1. Review RTPLAN file structure above for meterset locations\n');
fprintf('2. If metersets exist in DICOM but not in pln:\n');
fprintf('   - Manually extract from ReferencedBeamSequence\n');
fprintf('   - Populate pln.w before generating stf\n');
fprintf('3. If no metersets in DICOM:\n');
fprintf('   - Use uniform weights for initial calculation\n');
fprintf('   - Consider optimizing weights using matRad optimization\n');
fprintf('   - Check if separate delivery file exists with actual MU values\n');

fprintf('\n========================================\n');
fprintf('Diagnostics Complete\n');
fprintf('========================================\n');