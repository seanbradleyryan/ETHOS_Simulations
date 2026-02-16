function [field_doses, sct_resampled, total_rs_dose, metadata] = step15_process_doses(patient_id, session, config)
%% STEP15_PROCESS_DOSES - Process field doses and resample CT
%
%   [field_doses, sct_resampled, total_rs_dose, metadata] = step15_process_doses(patient_id, session, config)
%
%   PURPOSE:
%   Load all Raystation field dose DICOM files, extract dose grids with geometry
%   metadata, resample SCT to match dose grid, extract tissue classifications
%   from RTSTRUCT, zero out couch regions, and save processed data. Each
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
%           .dose_Gy       - 3D dose array in Gy
%           .origin        - [x, y, z] in mm
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%           .beam_num      - Beam number from filename (n in Beam[n])
%           .seg_num       - Segment number from filename (m in Seg[m])
%           .field_num     - Field number from filename (o in Field[o])
%           .gantry_angle  - Gantry angle in degrees (from RTPLAN)
%           .meterset      - Monitor units (from RTPLAN, matched by field_num)
%       sct_resampled   - Struct with CT resampled to dose grid:
%           .cubeHU         - 3D HU array
%           .cubeDensity    - 3D density array (kg/m³)
%           .tissueMask     - 3D uint8 array with ROI labels (0 = unassigned)
%           .roiNames       - Cell array of ROI names (index matches label)
%           .bodyMask       - 3D logical array (true = inside body region)
%           .couchMask      - 3D logical array (true = couch region)
%           .origin         - [x, y, z] in mm
%           .spacing        - [dx, dy, dz] in mm
%           .dimensions     - [nx, ny, nz]
%       total_rs_dose   - Sum of all field doses (3D array in Gy), zeroed outside body and in couch
%       metadata        - Struct with combined geometry info
%
%   FILES CREATED (in processed/ directory):
%       - field_dose_001.mat, field_dose_002.mat, ... (one per field)
%       - sct_resampled.mat (includes tissueMask and couchMask)
%       - total_rs_dose.mat
%       - tissue_masks.mat (individual ROI masks)
%       - metadata.mat
%
%   ALGORITHM:
%   1. Create processed/ subdirectory if not exists
%   2. Find all Beam*_Seg*_Field*.dcm files in Raystation directory
%   3. Load RTPLAN to extract beam metadata (gantry angles, metersets)
%   4. Load first dose file to establish reference grid geometry
%   5. Process each dose file, match field_num to RTPLAN beam, save individually
%   6. Load SCT images, sort by z-position
%   7. Resample SCT to dose grid via 3D interpolation
%   8. Convert HU to density
%   9. Load RTSTRUCT and create tissue classification masks
%   10. Zero out dose in couch regions
%   11. Save all outputs and return
%
%   KEY TECHNICAL NOTES:
%   - Z-resolution MUST come from GridFrameOffsetVector, NOT PixelSpacing
%   - Use squeeze() to remove singleton dimensions in dose arrays
%   - Standard HU to density: ρ = 1000 + HU (approximate)
%   - Dose zeroed where: outside body OR inside couch
%   - Raystation files: Beam[n]_Seg[m]_Field [o].dcm pattern
%   - RTPLAN files: RTPLAN*.dcm pattern
%   - RTSTRUCT files: RTSTRUCT*.dcm pattern
%   - Meterset matching: field_num (o) from filename matches beam_number in RTPLAN
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       [field_doses, sct, total_dose, meta] = step15_process_doses('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (dicominfo, dicomread, poly2mask)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.1 (Added RTSTRUCT tissue classification and couch masking)
%
%   See also: load_processed_data, step0_sort_dicom, step2_kwave_simulation

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

% SCT directory (contains CT images, RTPLAN, and RTSTRUCT)
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

fprintf('\n[1/8] Finding field dose files...\n');

% Find all Beam*_Seg*_Field*.dcm files (Raystation naming convention)
rd_files = dir(fullfile(rs_dir, 'Beam*_Seg*_Field*.dcm'));

% Also check for alternative patterns as fallback
if isempty(rd_files)
    rd_files = dir(fullfile(rs_dir, 'Beam*.dcm'));
end

if isempty(rd_files)
    % Legacy pattern fallback
    rd_files = dir(fullfile(rs_dir, 'RD.*.dcm'));
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'RD*.dcm'));
    end
end

if isempty(rd_files)
    error('step15_process_doses:NoFieldDoses', ...
        'No field dose files found in: %s\nExpected pattern: Beam*_Seg*_Field*.dcm', rs_dir);
end

num_files = length(rd_files);
fprintf('  Found %d field dose file(s)\n', num_files);

for i = 1:num_files
    fprintf('    [%d] %s\n', i, rd_files(i).name);
end

%% ======================== LOAD RTPLAN FOR BEAM METADATA ========================

fprintf('\n[2/8] Loading RTPLAN for beam metadata...\n');

beam_metadata = loadRtplanMetadata(sct_dir);

if ~isempty(beam_metadata)
    fprintf('  Loaded metadata for %d beams from RTPLAN\n', length(beam_metadata));
else
    fprintf('  [WARNING] No RTPLAN metadata available, using defaults\n');
end

%% ======================== ESTABLISH REFERENCE GRID ========================

fprintf('\n[3/8] Establishing reference dose grid geometry...\n');

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

fprintf('\n[4/8] Processing field doses (saving individually)...\n');

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
        
        % Extract beam info from filename (e.g., Beam1_Seg0_Field 1.dcm)
        % field_num is used to match to RTPLAN beam metadata
        [beam_num, seg_num, field_num] = extractBeamInfo(rd_files(i).name, i);
        
        % Get beam metadata by matching field_num to beam_number in RTPLAN
        [gantry_angle, meterset] = getBeamMetadata(beam_metadata, field_num);
        
        % Create field dose structure
        field_dose = struct();
        field_dose.dose_Gy = dose_data;
        field_dose.origin = dose_origin;
        field_dose.spacing = dose_spacing;
        field_dose.dimensions = dose_dims;
        field_dose.beam_num = beam_num;         % Beam number from filename
        field_dose.seg_num = seg_num;           % Segment number from filename
        field_dose.field_num = field_num;       % Field number (matches RTPLAN beam)
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
        fprintf('    Saved: %s (field %d, max: %.4f Gy, gantry: %.1f°, MU: %.1f)\n', ...
            field_filename, field_num, field_dose.max_dose_Gy, gantry_angle, meterset);
        
        % Store reference in output cell array (without full dose data for memory)
        field_doses{i} = struct();
        field_doses{i}.filepath = field_filepath;
        field_doses{i}.beam_num = beam_num;
        field_doses{i}.seg_num = seg_num;
        field_doses{i}.field_num = field_num;
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
fprintf('  Total dose max (before couch masking): %.4f Gy\n', max(total_rs_dose(:)));

% Update metadata
metadata.processed_count = processed_count;
metadata.total_dose_max_Gy_before_masking = max(total_rs_dose(:));

%% ======================== LOAD AND RESAMPLE SCT ========================

fprintf('\n[5/8] Loading and resampling SCT to dose grid...\n');

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

%% ======================== LOAD RTSTRUCT AND CREATE TISSUE MASKS ========================

fprintf('\n[6/8] Loading RTSTRUCT and creating tissue classification masks...\n');

% Load RTSTRUCT and convert contours to masks on the dose grid
[tissue_mask, roi_names, roi_masks, body_mask, couch_mask] = loadRtstructAndCreateMasks(...
    sct_dir, ref_origin, ref_spacing, ref_dims);

if isempty(tissue_mask)
    warning('step15_process_doses:NoRTSTRUCT', ...
        'Could not create tissue masks from RTSTRUCT. Using empty masks.');
    tissue_mask = zeros(ref_dims, 'uint8');
    roi_names = {};
    roi_masks = struct();
    body_mask = false(ref_dims);
    couch_mask = false(ref_dims);
else
    fprintf('  Created masks for %d ROIs\n', length(roi_names));
    fprintf('  Body voxels identified: %d\n', sum(body_mask(:)));
    fprintf('  Couch voxels identified: %d\n', sum(couch_mask(:)));
    
    % List ROIs
    for i = 1:length(roi_names)
        if isfield(roi_masks, sprintf('ROI_%03d', i))
            mask_field = sprintf('ROI_%03d', i);
            num_voxels = sum(roi_masks.(mask_field)(:));
            fprintf('    [%d] %s: %d voxels\n', i, roi_names{i}, num_voxels);
        end
    end
end

% Save tissue masks separately (can be large)
fprintf('  Saving tissue_masks.mat...\n');
tissue_masks_file = fullfile(processed_dir, 'tissue_masks.mat');
save(tissue_masks_file, 'tissue_mask', 'roi_names', 'roi_masks', 'body_mask', 'couch_mask', '-v7.3');
fprintf('  Saved: tissue_masks.mat\n');

%% ======================== ZERO OUT DOSE OUTSIDE BODY AND IN COUCH ========================

fprintf('\n[7/8] Zeroing out dose outside body and in couch regions...\n');

% Create mask for valid dose region: inside body AND not in couch
valid_dose_mask = body_mask & ~couch_mask;
invalid_dose_mask = ~valid_dose_mask;

% Statistics before zeroing
dose_outside_body = sum(total_rs_dose(~body_mask));
dose_in_couch = sum(total_rs_dose(couch_mask));
dose_to_zero = sum(total_rs_dose(invalid_dose_mask));

num_voxels_outside_body = sum(~body_mask(:));
num_voxels_in_couch = sum(couch_mask(:));
num_voxels_zeroed = sum(invalid_dose_mask(:));

fprintf('  Voxels outside body: %d\n', num_voxels_outside_body);
fprintf('  Voxels in couch: %d\n', num_voxels_in_couch);
fprintf('  Total voxels to zero (outside body OR in couch): %d\n', num_voxels_zeroed);
fprintf('  Dose outside body before zeroing: %.4f Gy (sum)\n', dose_outside_body);
fprintf('  Dose in couch before zeroing: %.4f Gy (sum)\n', dose_in_couch);

% Zero out dose in invalid regions
total_rs_dose(invalid_dose_mask) = 0;

fprintf('  Total dose max (after masking): %.4f Gy\n', max(total_rs_dose(:)));

% Update metadata
metadata.total_dose_max_Gy = max(total_rs_dose(:));
metadata.body_voxels = sum(body_mask(:));
metadata.couch_voxels = sum(couch_mask(:));
metadata.voxels_zeroed = num_voxels_zeroed;
metadata.dose_outside_body_zeroed = dose_outside_body;
metadata.dose_in_couch_zeroed = dose_in_couch;

% Also update individual field dose files to zero invalid regions
fprintf('  Updating individual field doses to zero invalid regions...\n');
for i = 1:num_files
    if ~isempty(field_doses{i}) && isfield(field_doses{i}, 'filepath')
        try
            field_filepath = field_doses{i}.filepath;
            loaded = load(field_filepath);
            field_dose = loaded.field_dose;
            
            % Zero out dose outside body and in couch
            field_dose.dose_Gy(invalid_dose_mask) = 0;
            field_dose.max_dose_Gy = max(field_dose.dose_Gy(:));
            field_dose.body_masked = true;
            field_dose.couch_masked = true;
            
            % Re-save
            save(field_filepath, 'field_dose', '-v7.3');
            
            % Update reference
            field_doses{i}.max_dose_Gy = field_dose.max_dose_Gy;
            field_doses{i}.body_masked = true;
            field_doses{i}.couch_masked = true;
            
            clear field_dose;
        catch ME
            warning('step15_process_doses:MaskError', ...
                'Failed to update masks for field %d: %s', i, ME.message);
        end
    end
end
fprintf('  Updated %d field dose files\n', num_files);

%% ======================== SAVE TOTAL DOSE ========================

fprintf('\n  Saving total_rs_dose.mat...\n');
total_dose_file = fullfile(processed_dir, 'total_rs_dose.mat');
save(total_dose_file, 'total_rs_dose', '-v7.3');
fprintf('  Saved: total_rs_dose.mat\n');

%% ======================== CREATE SCT RESAMPLED STRUCTURE ========================

sct_resampled = struct();
sct_resampled.cubeHU = sct_hu_resampled;
sct_resampled.cubeDensity = sct_density;
sct_resampled.tissueMask = tissue_mask;
sct_resampled.roiNames = roi_names;
sct_resampled.bodyMask = body_mask;
sct_resampled.couchMask = couch_mask;
sct_resampled.origin = ref_origin;
sct_resampled.spacing = ref_spacing;
sct_resampled.dimensions = ref_dims;
sct_resampled.patient_id = patient_id;
sct_resampled.session = session;
sct_resampled.original_sct_dims = sct_dims;
sct_resampled.original_sct_spacing = sct_spacing;
sct_resampled.timestamp = datetime('now');

%% ======================== SAVE SCT RESAMPLED ========================

fprintf('\n[8/8] Saving processed data...\n');

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
fprintf('  Tissue ROIs: %d\n', length(roi_names));
fprintf('  Body voxels: %d\n', sum(body_mask(:)));
fprintf('  Couch voxels: %d\n', sum(couch_mask(:)));
fprintf('  Voxels zeroed (outside body OR in couch): %d\n', num_voxels_zeroed);
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
    
    % Find RTPLAN file (RTPLAN*.dcm naming convention)
    rp_files = dir(fullfile(sct_dir, 'RTPLAN*.dcm'));
    
    if isempty(rp_files)
        % Try alternative naming patterns
        rp_files = dir(fullfile(sct_dir, 'RP*.dcm'));
    end
    
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


function [beam_num, seg_num, field_num] = extractBeamInfo(filename, default_index)
%EXTRACTBEAMINFO Extract beam, segment, and field numbers from dose filename
%
%   Handles patterns like: Beam1_Seg0_Field 1.dcm, Beam2_Seg0_Field 2.dcm, etc.
%   Also handles legacy patterns: RD.1.dcm, RD.2.dcm
%
%   OUTPUTS:
%       beam_num  - Beam number from filename (n in Beam[n])
%       seg_num   - Segment number from filename (m in Seg[m])
%       field_num - Field number from filename (o in Field [o])
%                   This is used to match to RTPLAN beam metadata

    beam_num = default_index;
    seg_num = 0;
    field_num = default_index;
    
    % Try new pattern: Beam[n]_Seg[m]_Field [o].dcm
    % Note: Field number may have a space before the number
    tokens = regexp(filename, 'Beam(\d+)_Seg(\d+)_Field\s*(\d+)', 'tokens');
    
    if ~isempty(tokens) && ~isempty(tokens{1})
        beam_num = str2double(tokens{1}{1});
        seg_num = str2double(tokens{1}{2});
        field_num = str2double(tokens{1}{3});
        return;
    end
    
    % Try alternative pattern without space: Beam[n]_Seg[m]_Field[o].dcm
    tokens = regexp(filename, 'Beam(\d+)_Seg(\d+)_Field(\d+)', 'tokens');
    
    if ~isempty(tokens) && ~isempty(tokens{1})
        beam_num = str2double(tokens{1}{1});
        seg_num = str2double(tokens{1}{2});
        field_num = str2double(tokens{1}{3});
        return;
    end
    
    % Legacy pattern: RD.[n].dcm
    tokens = regexp(filename, 'RD\.(\d+)', 'tokens');
    
    if ~isempty(tokens) && ~isempty(tokens{1})
        beam_num = str2double(tokens{1}{1});
        field_num = beam_num;
    end
end


function [gantry_angle, meterset] = getBeamMetadata(beam_metadata, field_num)
%GETBEAMMETADATA Get gantry angle and meterset for specific field
%
%   Matches field_num (from filename) to beam_number in RTPLAN metadata

    gantry_angle = 0;
    meterset = 0;
    
    if isempty(beam_metadata)
        return;
    end
    
    % Match field_num to beam_number in RTPLAN
    for i = 1:length(beam_metadata)
        if beam_metadata(i).beam_number == field_num
            gantry_angle = beam_metadata(i).gantry_angle;
            meterset = beam_metadata(i).meterset;
            return;
        end
    end
    
    % Fallback: use field_num as index if within range
    if field_num <= length(beam_metadata)
        gantry_angle = beam_metadata(field_num).gantry_angle;
        meterset = beam_metadata(field_num).meterset;
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


function [tissue_mask, roi_names, roi_masks, body_mask, couch_mask] = loadRtstructAndCreateMasks(...
    sct_dir, dose_origin, dose_spacing, dose_dims)
%LOADRTSTRUCTANDCREATEMASKS Load RTSTRUCT and create tissue classification masks
%
%   [tissue_mask, roi_names, roi_masks, body_mask, couch_mask] = loadRtstructAndCreateMasks(...)
%
%   Loads the RTSTRUCT file, extracts all ROI contours, and converts them
%   to 3D binary masks on the dose grid. Identifies body and couch regions.
%
%   INPUTS:
%       sct_dir      - Path to directory containing RTSTRUCT
%       dose_origin  - [x, y, z] origin of dose grid (mm)
%       dose_spacing - [dx, dy, dz] spacing of dose grid (mm)
%       dose_dims    - [rows, cols, slices] dimensions of dose grid
%
%   OUTPUTS:
%       tissue_mask  - 3D uint8 array with ROI labels (0 = unassigned)
%       roi_names    - Cell array of ROI names (index matches label in tissue_mask)
%       roi_masks    - Struct with individual ROI binary masks (ROI_001, ROI_002, ...)
%       body_mask    - 3D logical array (true = inside body region)
%       couch_mask   - 3D logical array (true = couch region)

    tissue_mask = [];
    roi_names = {};
    roi_masks = struct();
    body_mask = false(dose_dims);
    couch_mask = false(dose_dims);
    
    % Find RTSTRUCT file (RTSTRUCT*.dcm naming convention)
    rs_files = dir(fullfile(sct_dir, 'RTSTRUCT*.dcm'));
    
    if isempty(rs_files)
        % Try alternative naming patterns
        rs_files = dir(fullfile(sct_dir, 'RS*.dcm'));
    end
    
    if isempty(rs_files)
        warning('loadRtstructAndCreateMasks:NoRTSTRUCT', ...
            'No RTSTRUCT file found in: %s', sct_dir);
        return;
    end
    
    % Use first RTSTRUCT file found
    rs_file = fullfile(sct_dir, rs_files(1).name);
    fprintf('    Loading RTSTRUCT: %s\n', rs_files(1).name);
    
    try
        rtstruct = dicominfo(rs_file);
    catch ME
        warning('loadRtstructAndCreateMasks:LoadError', ...
            'Failed to load RTSTRUCT: %s', ME.message);
        return;
    end
    
    % Check for required sequences
    if ~isfield(rtstruct, 'StructureSetROISequence') || ...
       ~isfield(rtstruct, 'ROIContourSequence')
        warning('loadRtstructAndCreateMasks:MissingSequence', ...
            'RTSTRUCT missing StructureSetROISequence or ROIContourSequence');
        return;
    end
    
    % Extract ROI names from StructureSetROISequence
    roi_seq_fields = fieldnames(rtstruct.StructureSetROISequence);
    num_rois = length(roi_seq_fields);
    
    roi_info = struct();
    for i = 1:num_rois
        roi = rtstruct.StructureSetROISequence.(roi_seq_fields{i});
        roi_info(i).number = roi.ROINumber;
        roi_info(i).name = roi.ROIName;
    end
    
    fprintf('    Found %d ROIs in StructureSetROISequence\n', num_rois);
    
    % Initialize tissue mask
    tissue_mask = zeros(dose_dims, 'uint8');
    roi_names = cell(num_rois, 1);
    
    % Body region names to identify (case-insensitive matching)
    body_patterns = {'body', 'external', 'patient', 'skin', 'outer contour', ...
                     'body contour', 'external contour'};
    
    % Couch region names to identify (case-insensitive matching)
    couch_patterns = {'couch exterior', 'couch interior', 'couchexterior', ...
                      'couchinterior', 'couch_exterior', 'couch_interior', ...
                      'couch', 'table'};
    
    % Calculate z-coordinates for each slice in dose grid
    dose_z_coords = dose_origin(3) + (0:dose_dims(3)-1) * dose_spacing(3);
    
    % Process each ROI in ROIContourSequence
    contour_seq_fields = fieldnames(rtstruct.ROIContourSequence);
    
    for c_idx = 1:length(contour_seq_fields)
        try
            contour_item = rtstruct.ROIContourSequence.(contour_seq_fields{c_idx});
            
            % Get referenced ROI number
            if ~isfield(contour_item, 'ReferencedROINumber')
                continue;
            end
            ref_roi_num = contour_item.ReferencedROINumber;
            
            % Find matching ROI info
            roi_idx = find([roi_info.number] == ref_roi_num, 1);
            if isempty(roi_idx)
                continue;
            end
            
            roi_name = roi_info(roi_idx).name;
            roi_names{roi_idx} = roi_name;
            
            fprintf('    Processing ROI %d: %s\n', roi_idx, roi_name);
            
            % Check if this is a body region (case-insensitive)
            roi_name_lower = lower(roi_name);
            is_body = any(strcmpi(roi_name_lower, body_patterns)) || ...
                      any(contains(roi_name_lower, body_patterns));
            
            % Check if this is a couch region (case-insensitive)
            is_couch = any(strcmpi(roi_name_lower, couch_patterns)) || ...
                       any(contains(roi_name_lower, couch_patterns));
            
            % Initialize mask for this ROI
            roi_mask = false(dose_dims);
            
            % Get contour sequence
            if ~isfield(contour_item, 'ContourSequence')
                continue;
            end
            
            contour_fields = fieldnames(contour_item.ContourSequence);
            num_contours = length(contour_fields);
            
            % Process each contour (typically one per slice)
            for j = 1:num_contours
                contour = contour_item.ContourSequence.(contour_fields{j});
                
                % Get contour data (x, y, z triplets)
                if ~isfield(contour, 'ContourData') || isempty(contour.ContourData)
                    continue;
                end
                
                contour_data = contour.ContourData;
                num_points = length(contour_data) / 3;
                
                if num_points < 3
                    continue;  % Need at least 3 points for a polygon
                end
                
                % Reshape to [N x 3]
                points = reshape(contour_data, 3, num_points)';
                
                % Get contour z-coordinate
                contour_z = points(1, 3);  % All points should have same z
                
                % Find matching slice in dose grid
                [min_diff, slice_idx] = min(abs(dose_z_coords - contour_z));
                
                % Check if within tolerance (half slice thickness)
                if min_diff > dose_spacing(3)
                    continue;  % Contour not on this grid
                end
                
                % Convert contour points to pixel coordinates
                % X -> column, Y -> row
                col_coords = (points(:, 1) - dose_origin(1)) / dose_spacing(1) + 1;
                row_coords = (points(:, 2) - dose_origin(2)) / dose_spacing(2) + 1;
                
                % Create polygon mask for this slice
                try
                    slice_mask = poly2mask(col_coords, row_coords, dose_dims(1), dose_dims(2));
                    roi_mask(:, :, slice_idx) = roi_mask(:, :, slice_idx) | slice_mask;
                catch
                    % poly2mask can fail for degenerate polygons
                    continue;
                end
            end
            
            % Store individual ROI mask
            mask_field = sprintf('ROI_%03d', roi_idx);
            roi_masks.(mask_field) = roi_mask;
            
            % Add to tissue mask (later ROIs overwrite earlier ones in overlap)
            tissue_mask(roi_mask) = uint8(roi_idx);
            
            % Add to body mask if applicable
            if is_body
                body_mask = body_mask | roi_mask;
                fprintf('      -> Identified as BODY region (%d voxels)\n', sum(roi_mask(:)));
            end
            
            % Add to couch mask if applicable
            if is_couch
                couch_mask = couch_mask | roi_mask;
                fprintf('      -> Identified as COUCH region (%d voxels)\n', sum(roi_mask(:)));
            end
            
        catch ME
            warning('loadRtstructAndCreateMasks:ContourError', ...
                'Error processing contour %d: %s', c_idx, ME.message);
            continue;
        end
    end
    
    % Fill in empty ROI names for consistency
    for i = 1:length(roi_names)
        if isempty(roi_names{i})
            roi_names{i} = sprintf('ROI_%d_NoContour', i);
        end
    end
    
    % Count valid ROIs
    valid_count = sum(~contains(roi_names, 'NoContour'));
    fprintf('    Tissue mask created with %d labeled ROIs\n', valid_count);
    fprintf('    Body mask: %d voxels\n', sum(body_mask(:)));
    fprintf('    Couch mask: %d voxels\n', sum(couch_mask(:)));
end