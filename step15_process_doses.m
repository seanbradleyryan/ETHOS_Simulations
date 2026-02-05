function [field_doses, sct_resampled, total_rs_dose, metadata] = step15_process_doses(patient_id, session, config)
%% STEP15_PROCESS_DOSES - Process field doses and resample CT
%
%   [field_doses, sct_resampled, total_rs_dose, metadata] = step15_process_doses(patient_id, session, config)
%
%   PURPOSE:
%   Load all Raystation field dose DICOM files, extract dose grids with geometry
%   metadata, resample SCT to match dose grid, and save processed data. Each
%   field dose is saved as a separate file due to memory constraints.
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir    - Base directory path
%           .treatment_site - Subfolder name (default: 'Pancreas')
%
%   OUTPUTS:
%       field_doses     - Cell array of field dose structures (loaded from files)
%       sct_resampled   - Struct with CT resampled to dose grid:
%           .cubeHU       - 3D HU array
%           .cubeDensity  - 3D density array (kg/m³)
%           .origin       - [x, y, z] in mm
%           .spacing      - [dx, dy, dz] in mm
%           .dimensions   - [nx, ny, nz]
%       total_rs_dose   - Sum of all field doses (3D array in Gy)
%       metadata        - Struct with combined geometry info
%
%   FILES CREATED (in processed/ directory):
%       - field_dose_001.mat, field_dose_002.mat, ... (one per field)
%       - sct_resampled.mat
%       - total_rs_dose.mat
%       - metadata.mat
%
%   ALGORITHM:
%   1. Create processed/ subdirectory if not exists
%   2. Find all RD.*.dcm files in Raystation directory
%   3. Load RTPLAN to extract beam metadata (gantry angles, metersets)
%   4. Load first dose file to establish reference grid geometry
%   5. Process each dose file, save individually
%   6. Load SCT images, sort by z-position
%   7. Resample SCT to dose grid via 3D interpolation
%   8. Convert HU to density
%   9. Save all outputs and return
%
%   KEY TECHNICAL NOTES:
%   - Z-resolution MUST come from GridFrameOffsetVector, NOT PixelSpacing
%   - Use squeeze() to remove singleton dimensions in dose arrays
%   - Standard HU to density: ρ = 1000 + HU (approximate)
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       [field_doses, sct, total_dose, meta] = step15_process_doses('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (dicominfo, dicomread)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: step0_sort_dicom, step2_kwave_simulation

%% ======================== INPUT VALIDATION ========================

% Validate patient_id
if ~ischar(patient_id) && ~isstring(patient_id)
    error('step15_process_doses:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

% Validate session
if ~ischar(session) && ~isstring(session)
    error('step15_process_doses:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

% Validate config struct
if ~isstruct(config)
    error('step15_process_doses:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

% Validate required config fields
if ~isfield(config, 'working_dir')
    error('step15_process_doses:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

% Set default treatment_site if not provided
if ~isfield(config, 'treatment_site') || isempty(config.treatment_site)
    config.treatment_site = 'Pancreas';
    fprintf('  [INFO] Using default treatment_site: %s\n', config.treatment_site);
end

%% ======================== CONSTRUCT PATHS ========================

fprintf('\n========================================\n');
fprintf('  Step 1.5: Process Field Doses and Resample CT\n');
fprintf('  Patient: %s, Session: %s\n', patient_id, session);
fprintf('========================================\n');

% Raystation directory (contains RD.*.dcm field dose files)
rs_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session);

% SCT directory (contains CT images and RTPLAN)
sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session, 'sct');

% Processed output directory
processed_dir = fullfile(rs_dir, 'processed');

fprintf('  Raystation directory: %s\n', rs_dir);
fprintf('  SCT directory: %s\n', sct_dir);

%% ======================== VERIFY DIRECTORIES ========================

if ~isfolder(rs_dir)
    error('step15_process_doses:DirectoryNotFound', ...
        'Raystation directory does not exist: %s', rs_dir);
end

if ~isfolder(sct_dir)
    error('step15_process_doses:DirectoryNotFound', ...
        'SCT directory does not exist: %s', sct_dir);
end

% Create processed directory
if ~isfolder(processed_dir)
    mkdir(processed_dir);
    fprintf('  Created processed directory: %s\n', processed_dir);
else
    fprintf('  Processed directory exists: %s\n', processed_dir);
end

%% ======================== FIND FIELD DOSE FILES ========================

fprintf('\n[1/6] Finding field dose files...\n');

% Find all RD.*.dcm files (Raystation naming convention)
rd_files = dir(fullfile(rs_dir, 'RD.*.dcm'));

% Also check for RD*.dcm pattern as fallback
if isempty(rd_files)
    rd_files = dir(fullfile(rs_dir, 'RD*.dcm'));
end

if isempty(rd_files)
    error('step15_process_doses:NoFieldDoses', ...
        'No RD.*.dcm or RD*.dcm files found in: %s', rs_dir);
end

num_files = length(rd_files);
fprintf('  Found %d field dose file(s)\n', num_files);

for i = 1:num_files
    fprintf('    [%d] %s\n', i, rd_files(i).name);
end

%% ======================== LOAD RTPLAN FOR BEAM METADATA ========================

fprintf('\n[2/6] Loading RTPLAN for beam metadata...\n');

beam_metadata = loadRtplanMetadata(sct_dir);

if ~isempty(beam_metadata)
    fprintf('  Loaded metadata for %d beams from RTPLAN\n', length(beam_metadata));
else
    fprintf('  [WARNING] No RTPLAN metadata available, using defaults\n');
end

%% ======================== ESTABLISH REFERENCE GRID ========================

fprintf('\n[3/6] Establishing reference dose grid geometry...\n');

% Load first dose file to get reference geometry
ref_file = fullfile(rs_dir, rd_files(1).name);
ref_info = dicominfo(ref_file);
ref_dose = double(squeeze(dicomread(ref_file)));

% Apply dose grid scaling
if isfield(ref_info, 'DoseGridScaling')
    ref_dose = ref_dose * ref_info.DoseGridScaling;
end

% Extract reference geometry
% CRITICAL: Z-resolution from GridFrameOffsetVector, NOT PixelSpacing
ref_origin = ref_info.ImagePositionPatient(:);  % [x, y, z] in mm
ref_spacing = extractDoseSpacing(ref_info);      % [dx, dy, dz] in mm
ref_dims = size(ref_dose);                       % [rows, cols, slices]

fprintf('  Reference dose grid:\n');
fprintf('    Dimensions: [%d, %d, %d]\n', ref_dims(1), ref_dims(2), ref_dims(3));
fprintf('    Spacing (mm): [%.3f, %.3f, %.3f]\n', ref_spacing(1), ref_spacing(2), ref_spacing(3));
fprintf('    Origin (mm): [%.3f, %.3f, %.3f]\n', ref_origin(1), ref_origin(2), ref_origin(3));

% Initialize total dose accumulator
total_rs_dose = zeros(ref_dims);

% Initialize metadata structure
metadata = struct();
metadata.origin = ref_origin;
metadata.spacing = ref_spacing;
metadata.dimensions = ref_dims;
metadata.patient_id = patient_id;
metadata.session = session;
metadata.num_fields = num_files;
metadata.timestamp = datetime('now');
metadata.reference_file = rd_files(1).name;

%% ======================== PROCESS EACH FIELD DOSE ========================

fprintf('\n[4/6] Processing field doses (saving individually)...\n');

% Track which files were processed successfully
field_doses = cell(num_files, 1);
processed_count = 0;

for i = 1:num_files
    fprintf('  Processing field %d/%d: %s\n', i, num_files, rd_files(i).name);
    
    try
        % Load DICOM dose file
        dose_file = fullfile(rs_dir, rd_files(i).name);
        dose_info = dicominfo(dose_file);
        dose_data = double(squeeze(dicomread(dose_file)));
        
        % Apply dose grid scaling
        if isfield(dose_info, 'DoseGridScaling')
            dose_scaling = dose_info.DoseGridScaling;
            dose_data = dose_data * dose_scaling;
            fprintf('    Applied scaling: %e\n', dose_scaling);
        end
        
        % Extract geometry and verify it matches reference
        dose_origin = dose_info.ImagePositionPatient(:);
        dose_spacing = extractDoseSpacing(dose_info);
        dose_dims = size(dose_data);
        
        % Validate geometry matches reference
        [geom_match, geom_msg] = validateGeometry(dose_origin, dose_spacing, dose_dims, ...
            ref_origin, ref_spacing, ref_dims);
        
        if ~geom_match
            warning('step15_process_doses:GeometryMismatch', ...
                'Field %d geometry mismatch: %s', i, geom_msg);
            % Attempt to resample if dimensions don't match
            if ~isequal(dose_dims, ref_dims)
                fprintf('    Resampling to reference grid...\n');
                dose_data = resampleDoseToGrid(dose_data, dose_origin, dose_spacing, ...
                    ref_origin, ref_spacing, ref_dims);
                dose_dims = ref_dims;
            end
        end
        
        % Extract beam index from filename (e.g., RD.1.dcm -> beam 1)
        beam_index = extractBeamIndex(rd_files(i).name, i);
        
        % Get beam metadata if available
        [gantry_angle, meterset] = getBeamMetadata(beam_metadata, beam_index);
        
        % Create field dose structure
        field_dose = struct();
        field_dose.dose_Gy = dose_data;
        field_dose.origin = dose_origin;
        field_dose.spacing = dose_spacing;
        field_dose.dimensions = dose_dims;
        field_dose.beam_index = beam_index;
        field_dose.gantry_angle = gantry_angle;
        field_dose.meterset = meterset;
        field_dose.source_file = rd_files(i).name;
        field_dose.max_dose_Gy = max(dose_data(:));
        field_dose.mean_dose_Gy = mean(dose_data(dose_data > 0));
        
        % Add to total dose
        total_rs_dose = total_rs_dose + dose_data;
        
        % Save individual field dose file (MEMORY CONSTRAINT)
        field_filename = sprintf('field_dose_%03d.mat', i);
        field_filepath = fullfile(processed_dir, field_filename);
        save(field_filepath, 'field_dose', '-v7.3');
        fprintf('    Saved: %s (max: %.4f Gy, gantry: %.1f°)\n', ...
            field_filename, field_dose.max_dose_Gy, gantry_angle);
        
        % Store reference in output cell array (without full dose data for memory)
        field_doses{i} = struct();
        field_doses{i}.filepath = field_filepath;
        field_doses{i}.beam_index = beam_index;
        field_doses{i}.gantry_angle = gantry_angle;
        field_doses{i}.meterset = meterset;
        field_doses{i}.max_dose_Gy = field_dose.max_dose_Gy;
        field_doses{i}.source_file = rd_files(i).name;
        
        processed_count = processed_count + 1;
        
        % Clear field_dose to free memory
        clear field_dose dose_data;
        
    catch ME
        warning('step15_process_doses:FieldProcessingError', ...
            'Failed to process field %d (%s): %s', i, rd_files(i).name, ME.message);
    end
end

fprintf('  Successfully processed %d/%d field doses\n', processed_count, num_files);
fprintf('  Total dose max: %.4f Gy\n', max(total_rs_dose(:)));

% Update metadata
metadata.processed_count = processed_count;
metadata.total_dose_max_Gy = max(total_rs_dose(:));

%% ======================== SAVE TOTAL DOSE ========================

fprintf('\n  Saving total_rs_dose.mat...\n');
total_dose_file = fullfile(processed_dir, 'total_rs_dose.mat');
save(total_dose_file, 'total_rs_dose', '-v7.3');
fprintf('  Saved: total_rs_dose.mat\n');

%% ======================== LOAD AND RESAMPLE SCT ========================

fprintf('\n[5/6] Loading and resampling SCT to dose grid...\n');

% Load SCT images
[sct_hu, sct_origin, sct_spacing, sct_dims] = loadSctImages(sct_dir);

if isempty(sct_hu)
    error('step15_process_doses:NoSCT', ...
        'Failed to load SCT images from: %s', sct_dir);
end

fprintf('  Original SCT:\n');
fprintf('    Dimensions: [%d, %d, %d]\n', sct_dims(1), sct_dims(2), sct_dims(3));
fprintf('    Spacing (mm): [%.3f, %.3f, %.3f]\n', sct_spacing(1), sct_spacing(2), sct_spacing(3));
fprintf('    Origin (mm): [%.3f, %.3f, %.3f]\n', sct_origin(1), sct_origin(2), sct_origin(3));
fprintf('    HU range: [%.0f, %.0f]\n', min(sct_hu(:)), max(sct_hu(:)));

% Resample SCT to dose grid
fprintf('  Performing 3D interpolation (this may take a moment)...\n');

sct_hu_resampled = resampleSctToDoseGrid(sct_hu, sct_origin, sct_spacing, sct_dims, ...
    ref_origin, ref_spacing, ref_dims);

fprintf('  Resampled SCT dimensions: [%d, %d, %d]\n', ...
    size(sct_hu_resampled, 1), size(sct_hu_resampled, 2), size(sct_hu_resampled, 3));

% Convert HU to density
fprintf('  Converting HU to density...\n');
sct_density = huToDensity(sct_hu_resampled);

fprintf('  Density range: [%.0f, %.0f] kg/m³\n', min(sct_density(:)), max(sct_density(:)));

%% ======================== CREATE SCT RESAMPLED STRUCTURE ========================

sct_resampled = struct();
sct_resampled.cubeHU = sct_hu_resampled;
sct_resampled.cubeDensity = sct_density;
sct_resampled.origin = ref_origin;
sct_resampled.spacing = ref_spacing;
sct_resampled.dimensions = ref_dims;
sct_resampled.patient_id = patient_id;
sct_resampled.session = session;
sct_resampled.original_sct_dims = sct_dims;
sct_resampled.original_sct_spacing = sct_spacing;
sct_resampled.timestamp = datetime('now');

%% ======================== SAVE SCT RESAMPLED ========================

fprintf('\n[6/6] Saving processed data...\n');

sct_resampled_file = fullfile(processed_dir, 'sct_resampled.mat');
save(sct_resampled_file, 'sct_resampled', '-v7.3');
fprintf('  Saved: sct_resampled.mat\n');

% Save metadata
metadata_file = fullfile(processed_dir, 'metadata.mat');
save(metadata_file, 'metadata', '-v7.3');
fprintf('  Saved: metadata.mat\n');

%% ======================== RELOAD FIELD DOSES FOR OUTPUT ========================

% Optionally reload full field doses into cell array for return
% (Only if memory permits - otherwise caller should load from files)
fprintf('\n  Reloading field doses for output...\n');

try
    for i = 1:num_files
        if ~isempty(field_doses{i}) && isfield(field_doses{i}, 'filepath')
            loaded = load(field_doses{i}.filepath);
            field_doses{i} = loaded.field_dose;
        end
    end
    fprintf('  Field doses loaded into memory\n');
catch ME
    warning('step15_process_doses:MemoryWarning', ...
        'Could not reload all field doses into memory: %s\nAccess them from individual files.', ...
        ME.message);
end

%% ======================== SUMMARY ========================

fprintf('\n========================================\n');
fprintf('  Step 1.5 Complete\n');
fprintf('========================================\n');
fprintf('  Processed %d field doses\n', processed_count);
fprintf('  Dose grid: [%d x %d x %d]\n', ref_dims(1), ref_dims(2), ref_dims(3));
fprintf('  Spacing: [%.3f, %.3f, %.3f] mm\n', ref_spacing(1), ref_spacing(2), ref_spacing(3));
fprintf('  Total dose max: %.4f Gy\n', max(total_rs_dose(:)));
fprintf('  Output directory: %s\n', processed_dir);
fprintf('========================================\n\n');

end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function spacing = extractDoseSpacing(dose_info)
%EXTRACTDOSESPACING Extract dose grid spacing from DICOM info
%
%   CRITICAL: Z-resolution MUST come from GridFrameOffsetVector, NOT PixelSpacing
%   This is a common source of bugs in dose processing.

    % X and Y spacing from PixelSpacing
    dx = dose_info.PixelSpacing(1);  % Column spacing
    dy = dose_info.PixelSpacing(2);  % Row spacing
    
    % Z spacing from GridFrameOffsetVector
    if isfield(dose_info, 'GridFrameOffsetVector') && length(dose_info.GridFrameOffsetVector) >= 2
        dz = abs(dose_info.GridFrameOffsetVector(2) - dose_info.GridFrameOffsetVector(1));
    elseif isfield(dose_info, 'SliceThickness')
        dz = dose_info.SliceThickness;
        warning('extractDoseSpacing:NoGridFrameOffset', ...
            'GridFrameOffsetVector not available, using SliceThickness = %.3f mm', dz);
    else
        dz = dx;  % Fallback to isotropic assumption
        warning('extractDoseSpacing:NoZSpacing', ...
            'Cannot determine Z spacing, assuming isotropic: %.3f mm', dz);
    end
    
    spacing = [dx; dy; dz];
end


function beam_metadata = loadRtplanMetadata(sct_dir)
%LOADRTPLANMETADATA Load beam metadata from RTPLAN file
%
%   Extract gantry angles and metersets for each beam

    beam_metadata = [];
    
    % Find RTPLAN file
    rp_files = dir(fullfile(sct_dir, 'RP*.dcm'));
    
    if isempty(rp_files)
        return;
    end
    
    % Prefer adjusted MLC plan if available
    adjusted_idx = find(contains({rp_files.name}, 'adjusted_mlc'), 1);
    if ~isempty(adjusted_idx)
        rp_file = fullfile(sct_dir, rp_files(adjusted_idx).name);
    else
        rp_file = fullfile(sct_dir, rp_files(1).name);
    end
    
    try
        rtplan = dicominfo(rp_file);
        
        if ~isfield(rtplan, 'BeamSequence')
            return;
        end
        
        beam_fields = fieldnames(rtplan.BeamSequence);
        num_beams = length(beam_fields);
        
        beam_metadata = struct();
        
        for i = 1:num_beams
            beam = rtplan.BeamSequence.(beam_fields{i});
            
            % Beam number
            if isfield(beam, 'BeamNumber')
                beam_metadata(i).beam_number = beam.BeamNumber;
            else
                beam_metadata(i).beam_number = i;
            end
            
            % Beam name
            if isfield(beam, 'BeamName')
                beam_metadata(i).beam_name = beam.BeamName;
            else
                beam_metadata(i).beam_name = sprintf('Beam_%d', i);
            end
            
            % Gantry angle from first control point
            beam_metadata(i).gantry_angle = 0;
            if isfield(beam, 'ControlPointSequence')
                cp_fields = fieldnames(beam.ControlPointSequence);
                if ~isempty(cp_fields)
                    cp1 = beam.ControlPointSequence.(cp_fields{1});
                    if isfield(cp1, 'GantryAngle')
                        beam_metadata(i).gantry_angle = cp1.GantryAngle;
                    end
                end
            end
            
            % Meterset (from FractionGroupSequence)
            beam_metadata(i).meterset = 0;
        end
        
        % Extract metersets from FractionGroupSequence
        if isfield(rtplan, 'FractionGroupSequence')
            fg_fields = fieldnames(rtplan.FractionGroupSequence);
            fg = rtplan.FractionGroupSequence.(fg_fields{1});
            
            if isfield(fg, 'ReferencedBeamSequence')
                ref_beam_fields = fieldnames(fg.ReferencedBeamSequence);
                
                for i = 1:length(ref_beam_fields)
                    ref_beam = fg.ReferencedBeamSequence.(ref_beam_fields{i});
                    
                    if isfield(ref_beam, 'BeamMeterset') && isfield(ref_beam, 'ReferencedBeamNumber')
                        beam_num = ref_beam.ReferencedBeamNumber;
                        
                        % Find matching beam in metadata
                        for j = 1:length(beam_metadata)
                            if beam_metadata(j).beam_number == beam_num
                                beam_metadata(j).meterset = ref_beam.BeamMeterset;
                                break;
                            end
                        end
                    end
                end
            end
        end
        
    catch ME
        warning('loadRtplanMetadata:Error', ...
            'Failed to load RTPLAN metadata: %s', ME.message);
        beam_metadata = [];
    end
end


function [match, msg] = validateGeometry(dose_origin, dose_spacing, dose_dims, ...
    ref_origin, ref_spacing, ref_dims)
%VALIDATEGEOMETRY Validate that dose geometry matches reference grid

    match = true;
    msg = '';
    
    % Tolerance for floating point comparison
    tol = 0.01;  % mm
    
    % Check dimensions
    if ~isequal(dose_dims, ref_dims)
        match = false;
        msg = sprintf('Dimensions mismatch: [%d,%d,%d] vs [%d,%d,%d]', ...
            dose_dims(1), dose_dims(2), dose_dims(3), ...
            ref_dims(1), ref_dims(2), ref_dims(3));
        return;
    end
    
    % Check origin
    origin_diff = abs(dose_origin - ref_origin);
    if any(origin_diff > tol)
        match = false;
        msg = sprintf('Origin mismatch: max diff = %.3f mm', max(origin_diff));
        return;
    end
    
    % Check spacing
    spacing_diff = abs(dose_spacing - ref_spacing);
    if any(spacing_diff > tol)
        match = false;
        msg = sprintf('Spacing mismatch: max diff = %.3f mm', max(spacing_diff));
        return;
    end
end


function beam_index = extractBeamIndex(filename, default_index)
%EXTRACTBEAMINDEX Extract beam index from dose filename
%
%   Handles patterns like: RD.1.dcm, RD.2.dcm, RD.1.2.156.dcm, etc.

    beam_index = default_index;
    
    % Try to extract number after "RD."
    tokens = regexp(filename, 'RD\.(\d+)', 'tokens');
    
    if ~isempty(tokens) && ~isempty(tokens{1})
        beam_index = str2double(tokens{1}{1});
    end
end


function [gantry_angle, meterset] = getBeamMetadata(beam_metadata, beam_index)
%GETBEAMMETADATA Get gantry angle and meterset for specific beam

    gantry_angle = 0;
    meterset = 0;
    
    if isempty(beam_metadata)
        return;
    end
    
    % Try to match by beam number
    for i = 1:length(beam_metadata)
        if beam_metadata(i).beam_number == beam_index
            gantry_angle = beam_metadata(i).gantry_angle;
            meterset = beam_metadata(i).meterset;
            return;
        end
    end
    
    % Fallback: use index directly if within range
    if beam_index <= length(beam_metadata)
        gantry_angle = beam_metadata(beam_index).gantry_angle;
        meterset = beam_metadata(beam_index).meterset;
    end
end


function resampled_dose = resampleDoseToGrid(dose_data, dose_origin, dose_spacing, ...
    ref_origin, ref_spacing, ref_dims)
%RESAMPLEDOSETOGRID Resample dose array to reference grid

    dose_dims = size(dose_data);
    
    % Create coordinate grids for source dose
    [src_x, src_y, src_z] = ndgrid(...
        dose_origin(1) + (0:dose_dims(1)-1) * dose_spacing(1), ...
        dose_origin(2) + (0:dose_dims(2)-1) * dose_spacing(2), ...
        dose_origin(3) + (0:dose_dims(3)-1) * dose_spacing(3));
    
    % Create coordinate grids for target (reference)
    [tgt_x, tgt_y, tgt_z] = ndgrid(...
        ref_origin(1) + (0:ref_dims(1)-1) * ref_spacing(1), ...
        ref_origin(2) + (0:ref_dims(2)-1) * ref_spacing(2), ...
        ref_origin(3) + (0:ref_dims(3)-1) * ref_spacing(3));
    
    % Interpolate
    resampled_dose = interp3(src_y, src_x, src_z, dose_data, ...
        tgt_y, tgt_x, tgt_z, 'linear', 0);
end


function [sct_hu, origin, spacing, dims] = loadSctImages(sct_dir)
%LOADSCTIMAGES Load SCT DICOM images and extract geometry
%
%   Load all CT DICOM files, sort by z-position, convert to HU

    sct_hu = [];
    origin = [];
    spacing = [];
    dims = [];
    
    % Get all DICOM files that contain 'CT' in the name
    all_files = dir(fullfile(sct_dir, '*.dcm'));
    sct_files = [];
    
    for i = 1:length(all_files)
        if contains(all_files(i).name, 'CT', 'IgnoreCase', true) && ...
           ~contains(all_files(i).name, 'RTSTRUCT', 'IgnoreCase', true) && ...
           ~contains(all_files(i).name, 'RTPLAN', 'IgnoreCase', true) && ...
           ~contains(all_files(i).name, 'RTDOSE', 'IgnoreCase', true)
            sct_files = [sct_files; all_files(i)]; %#ok<AGROW>
        end
    end
    
    if isempty(sct_files)
        warning('loadSctImages:NoFiles', 'No CT DICOM files found in: %s', sct_dir);
        return;
    end
    
    fprintf('    Found %d CT image files\n', length(sct_files));
    
    % Get z-positions for sorting
    z_positions = zeros(length(sct_files), 1);
    
    for i = 1:length(sct_files)
        try
            info = dicominfo(fullfile(sct_dir, sct_files(i).name));
            z_positions(i) = info.ImagePositionPatient(3);
        catch
            z_positions(i) = i;  % Fallback
        end
    end
    
    % Sort by z-position
    [~, sort_idx] = sort(z_positions);
    
    % Load first slice for dimensions and metadata
    first_info = dicominfo(fullfile(sct_dir, sct_files(sort_idx(1)).name));
    first_img = dicomread(fullfile(sct_dir, sct_files(sort_idx(1)).name));
    
    % Initialize 3D array
    num_slices = length(sct_files);
    sct_data = zeros([size(first_img), num_slices], 'int16');
    
    % Load all slices in sorted order
    fprintf('    Loading %d SCT slices...\n', num_slices);
    
    for i = 1:num_slices
        idx = sort_idx(i);
        img = dicomread(fullfile(sct_dir, sct_files(idx).name));
        sct_data(:, :, i) = img;
    end
    
    % Extract geometry
    origin = first_info.ImagePositionPatient(:);
    
    % Get spacing
    dx = first_info.PixelSpacing(1);
    dy = first_info.PixelSpacing(2);
    
    % Z spacing from slice positions or SliceThickness
    if num_slices >= 2
        second_info = dicominfo(fullfile(sct_dir, sct_files(sort_idx(2)).name));
        dz = abs(second_info.ImagePositionPatient(3) - first_info.ImagePositionPatient(3));
    elseif isfield(first_info, 'SliceThickness')
        dz = first_info.SliceThickness;
    else
        dz = dx;  % Assume isotropic
    end
    
    spacing = [dx; dy; dz];
    dims = size(sct_data);
    
    % Convert to Hounsfield Units
    if isfield(first_info, 'RescaleSlope') && isfield(first_info, 'RescaleIntercept')
        sct_hu = double(sct_data) * first_info.RescaleSlope + first_info.RescaleIntercept;
    else
        sct_hu = double(sct_data);
    end
end


function sct_resampled = resampleSctToDoseGrid(sct_hu, sct_origin, sct_spacing, sct_dims, ...
    dose_origin, dose_spacing, dose_dims)
%RESAMPLESCTTODOSEGRID Resample SCT to dose grid using 3D interpolation
%
%   Uses interp3 with linear interpolation, -1000 HU for out-of-bounds

    % Create coordinate grids for SCT (source)
    % Note: meshgrid ordering is [Y, X, Z] but data is [rows, cols, slices]
    [sct_x, sct_y, sct_z] = meshgrid(...
        sct_origin(1) + (0:sct_dims(2)-1) * sct_spacing(1), ...   % X (columns)
        sct_origin(2) + (0:sct_dims(1)-1) * sct_spacing(2), ...   % Y (rows)
        sct_origin(3) + (0:sct_dims(3)-1) * sct_spacing(3));      % Z (slices)
    
    % Create coordinate grids for dose (target)
    [dose_x, dose_y, dose_z] = meshgrid(...
        dose_origin(1) + (0:dose_dims(2)-1) * dose_spacing(1), ...
        dose_origin(2) + (0:dose_dims(1)-1) * dose_spacing(2), ...
        dose_origin(3) + (0:dose_dims(3)-1) * dose_spacing(3));
    
    % Interpolate SCT to dose grid
    % Use -1000 HU (air) for extrapolated values
    sct_resampled = interp3(sct_x, sct_y, sct_z, sct_hu, ...
        dose_x, dose_y, dose_z, 'linear', -1000);
end


function density = huToDensity(hu)
%HUTODENSITY Convert Hounsfield Units to density (kg/m³)
%
%   Uses simplified linear conversion:
%   - Below -1000 HU (air): density = 1 kg/m³
%   - -1000 to 0 HU: linear interpolation from air (1) to water (1000)
%   - 0 to 1000 HU: linear from water (1000) to bone (~2000)
%   - Above 1000 HU: bone/metal region
%
%   Standard approximation: density = 1000 + HU (for soft tissue range)

    density = zeros(size(hu));
    
    % Air region (HU < -900)
    air_mask = hu < -900;
    density(air_mask) = 1.2;  % Air density
    
    % Lung region (-900 to -500 HU)
    lung_mask = (hu >= -900) & (hu < -500);
    density(lung_mask) = 400 + (hu(lung_mask) + 900) * (1000 - 400) / 400;
    
    % Soft tissue region (-500 to 100 HU)
    % Approximate: density ≈ 1000 + HU
    soft_mask = (hu >= -500) & (hu < 100);
    density(soft_mask) = 1000 + hu(soft_mask);
    
    % Bone region (HU >= 100)
    bone_mask = hu >= 100;
    % Linear interpolation from soft tissue to dense bone
    density(bone_mask) = 1100 + (hu(bone_mask) - 100) * (1900 - 1100) / 900;
    
    % Clamp to reasonable range
    density = max(density, 1);      % Minimum 1 kg/m³
    density = min(density, 7800);   % Maximum (metal)
end
