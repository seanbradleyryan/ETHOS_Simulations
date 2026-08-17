function [field_doses, cbct_resampled, total_rs_dose, metadata] = step15_process_doses(patient_id, session, config)
%% STEP15_PROCESS_DOSES - Process field doses and resample per-CBCT geometry
%
%   [field_doses, cbct_resampled, total_rs_dose, metadata] = step15_process_doses(patient_id, session, config)
%
%   PURPOSE:
%   Load all Raystation field dose files, extract dose grids with geometry
%   metadata, resample BOTH CBCTs (CT_1 and CT_3) found in the RayStation
%   directory to the dose grid, build per-CBCT tissue/body/couch masks from
%   each CBCT's own RTSTRUCT, zero out couch regions per field, and save
%   processed data. Each field dose is saved as a separate file due to
%   memory constraints.
%
%   The earlier-acquired CBCT (by SeriesDate+SeriesTime) becomes CT_1, the
%   later becomes CT_3. Field doses without a CT_1/CT_3 label in their
%   filename are rejected — per-CBCT masking has no sensible default.
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir        - Base directory path
%           .treatment_site     - Subfolder name (default: 'Pancreas')
%           .apply_dose_masking - Boolean, zero dose outside body/in couch
%                                 (default: true, set false for debugging)
%
%   OUTPUTS:
%       field_doses     - Cell array of field dose structures (loaded from files)
%           .dose_Gy       - 3D dose array in Gy
%           .origin        - [x, y, z] in mm
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%           .beam_num      - Beam number from filename (n in B[n])
%           .seg_num       - Segment number from filename
%           .field_num     - Field number (= beam_num, used to match RTPLAN)
%           .plan_type     - 'adapted' or 'reference' (from filename)
%           .gantry_angle  - Gantry angle in degrees (from RTPLAN)
%           .meterset      - Monitor units (from RTPLAN, matched by beam_num)
%       cbct_resampled  - Struct with two fields, CT_1 and CT_3, each a
%                         per-CBCT resampled struct with:
%           .cubeHU              - 3D HU array
%           .cubeDensity         - 3D density array (kg/mÂ³)
%           .tissueMask          - 3D uint8 array with ROI labels (0 = unassigned)
%           .roiNames            - Cell array of ROI names (index matches label)
%           .bodyMask            - 3D logical array (true = inside body region)
%           .couchMask           - 3D logical array (true = couch region)
%           .origin/.spacing     - dose-grid geometry (mm)
%           .dimensions          - [nx, ny, nz]
%           .series_uid          - source CBCT SeriesInstanceUID
%           .series_datetime     - numeric SeriesDate+SeriesTime
%           .ct_label            - 'CT_1' or 'CT_3'
%       total_rs_dose   - Sum of all field doses (3D array in Gy), zeroed via
%                         the union of CT_1/CT_3 body masks
%       metadata        - Struct with combined geometry info
%
%   FILES CREATED (in processed/ directory):
%       - dose_[id]_[session]_[adapted|reference]_{CT_n}_B[n]_[seg].mat (per field)
%       - CBCT1_resampled.mat  (variable: CBCT1_resampled)
%       - CBCT3_resampled.mat  (variable: CBCT3_resampled)
%       - total_rs_dose.mat
%       - total_dose_[CT_label].mat (per-CT-image total; e.g. total_dose_CT_1.mat, total_dose_CT_3.mat)
%       - tissue_masks.mat (per-CBCT masks: *_ct1 and *_ct3 variables)
%       - metadata.mat
%
%   ALGORITHM:
%   1. Create processed/ subdirectory if not exists
%   2. Find all dose_*.{mat,dcm} files in Raystation directory
%   3. Load RTPLAN to extract beam metadata (gantry angles, metersets)
%   4. Load first dose file to establish reference grid geometry
%   5. Discover the two CBCT series + their RTSTRUCTs in RayStationFiles,
%      pair by SeriesInstanceUID, sort by datetime -> CT_1 / CT_3
%   6. Resample each CBCT to the dose grid; build per-CBCT body/couch/tissue
%      masks from the matching RTSTRUCT; precompute per-CBCT invalid masks
%   7. Process each dose file: scale -> validate -> zero invalid regions
%      (per the field's ct_label) -> save
%   8. Zero combined total_rs_dose using union-body / intersection-couch;
%      zero each per-CT total with its own mask; save all outputs
%
%   KEY TECHNICAL NOTES:
%   - Z-resolution MUST come from GridFrameOffsetVector, NOT PixelSpacing
%   - Use squeeze() to remove singleton dimensions in dose arrays
%   - Standard HU to density: Ï = 1000 + HU (approximate)
%   - Dose zeroed where: outside body OR inside couch (unless disabled)
%   - Set config.apply_dose_masking = false to skip dose zeroing (debugging)
%   - RayStation files: dose_[patientid]_[session]_[adapted|reference]_B[n]_[seg].dcm
%       e.g. dose_1885729_Session_4_adapted_B6_103.dcm
%       beam_num = B[n], seg_num = [seg], field_num = beam_num (matches RTPLAN)
%       plan_type = 'adapted' or 'reference'
%   - RTPLAN files: RTPLAN*.dcm pattern
%   - RTSTRUCT files: RTSTRUCT*.dcm pattern
%   - Meterset matching: beam_num (n) from filename matches beam_number in RTPLAN
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       config.apply_dose_masking = true;  % Set false to disable dose zeroing (debugging)
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

% Set default for dose masking (for debugging, can disable zeroing outside body/couch)
if ~isfield(config, 'apply_dose_masking')
    config.apply_dose_masking = true;  % Default: apply masking
end
if ~config.apply_dose_masking
    fprintf('  [INFO] Dose masking DISABLED (debugging mode)\n');
end

% Set default batch size for field dose processing (memory management)
if ~isfield(config, 'batch_size') || isempty(config.batch_size)
    config.batch_size = 1000;
end

% Resume-friendly behavior: detect which outputs in processed/ already
% exist and skip the corresponding work. Set false for a forced full re-run.
if ~isfield(config, 'skip_completed')
    config.skip_completed = true;
end

% Store dose/mask arrays as sparse 2D matrices for compressibility.
% Dose grids are mostly zero after body/couch masking.
% Reconstruct on load: reshape(full(var_sp), var_dims)
if ~isfield(config, 'use_sparse_storage')
    config.use_sparse_storage = true;
end
if config.use_sparse_storage
    fprintf('  [INFO] Sparse storage ENABLED — dose/mask arrays saved as sparse 2D\n');
    fprintf('         To reconstruct: reshape(full(var_sp), var_dims)\n');
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

% --- Input format priority ---
%   1. dose_*.mat   (preferred, written by step14_npz_to_mat from RayStation NPZ)
%   2. dose_*.dcm   (legacy RayStation DICOM export)
%   3. older legacy patterns (Plan_Field*, RD.*, ...)
%
% input_format selects the per-field loader branch later in this function.
input_format       = '';
rd_files_mat       = dir(fullfile(rs_dir, 'dose_*.mat'));

if ~isempty(rd_files_mat)
    rd_files = rd_files_mat;
    input_format = 'mat';
    fprintf('  Using converted .mat input (dose_*.mat): %d file(s) found\n', numel(rd_files));
else
    rd_files_preferred = dir(fullfile(rs_dir, 'dose_*.dcm'));
end

if isempty(input_format) && ~isempty(rd_files_preferred)
    rd_files = rd_files_preferred;
    input_format = 'dicom';
    fprintf('  Using preferred format (dose_*.dcm): %d file(s) found\n', numel(rd_files));
elseif isempty(input_format)
    input_format = 'dicom';
    % --- Legacy fallback patterns ---
    rd_files = dir(fullfile(rs_dir, 'Plan_Field*_Beam*_B*_S*.dcm'));

    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'Beam*_Seg*_Field*.dcm'));
    end

    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'Beam*.dcm'));
    end

    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'RD.*.dcm'));
        if isempty(rd_files)
            rd_files = dir(fullfile(rs_dir, 'RD*.dcm'));
        end
    end

    if ~isempty(rd_files)
        fprintf('  Preferred format not found; using legacy pattern: %d file(s)\n', numel(rd_files));
    end
end

if isempty(rd_files)
    error('step15_process_doses:NoFieldDoses', ...
        ['No field dose files found in: %s\n' ...
         'Searched patterns (in priority order):\n' ...
         '  1. dose_*.mat  (preferred, from step14_npz_to_mat)\n' ...
         '  2. dose_*.dcm  (legacy RayStation DICOM export)\n' ...
         '  3. Plan_Field*_Beam*_B*_S*.dcm\n' ...
         '  4. Beam*_Seg*_Field*.dcm\n' ...
         '  5. Beam*.dcm\n' ...
         '  6. RD.*.dcm / RD*.dcm'], rs_dir);
end

num_files = length(rd_files);
fprintf('  Found %d field dose file(s)\n', num_files);

for i = 1:num_files
    fprintf('    [%d] %s\n', i, rd_files(i).name);
end

%% ======================== DETECT EXISTING PROCESSED OUTPUTS ========================

fprintf('\n[1.5/8] Detecting existing processed outputs (skip_completed=%d)...\n', ...
    config.skip_completed);

[skip_status, expected_field_outputs] = detectProcessedOutputs( ...
    processed_dir, rd_files, patient_id, session);

fprintf('  CBCT1_resampled.mat:   %s\n', presentStr(skip_status.cbct1_done));
fprintf('  CBCT3_resampled.mat:   %s\n', presentStr(skip_status.cbct3_done));
fprintf('  tissue_masks.mat:      %s\n', presentStr(skip_status.masks_done));
fprintf('  total_rs_dose.mat:     %s\n', presentStr(skip_status.total_rs_done));
fprintf('  total_dose_CT_1.mat:   %s\n', presentStr(skip_status.ct1_total_done));
fprintf('  total_dose_CT_3.mat:   %s\n', presentStr(skip_status.ct3_total_done));
fprintf('  metadata.mat:          %s\n', presentStr(skip_status.metadata_done));
fprintf('  Per-field outputs:     %d / %d present\n', ...
    skip_status.num_fields_done, num_files);

cbct_ready    = config.skip_completed && skip_status.cbct1_done && skip_status.cbct3_done;
masks_ready   = config.skip_completed && skip_status.masks_done;
fields_ready  = config.skip_completed && (skip_status.num_fields_done == num_files);
totals_ready  = config.skip_completed && skip_status.total_rs_done && ...
                skip_status.ct1_total_done && skip_status.ct3_total_done && ...
                skip_status.metadata_done;

if cbct_ready && masks_ready && fields_ready && totals_ready
    fprintf('\n  All Step 1.5 outputs already exist — loading and returning early.\n');
    [field_doses, cbct_resampled, total_rs_dose, metadata] = ...
        loadExistingStep15Outputs(processed_dir, expected_field_outputs);
    fprintf('  Loaded %d field dose entries, total_rs_dose max=%.4f Gy\n', ...
        numel(field_doses), max(total_rs_dose(:)));
    return;
end

% Whether the per-field loop needs to accumulate into running totals.
% If every field output is on disk AND every totals file is on disk, we
% can skip total accumulation entirely and load existing totals at the end.
need_total_accum = ~(fields_ready && totals_ready);

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

% Load first dose file to get reference geometry. Both .mat and DICOM
% inputs are pre-converted to the same (origin_mm, spacing_mm, dimensions)
% triple before the per-field loop runs.
ref_file = fullfile(rs_dir, rd_files(1).name);
switch input_format
    case 'mat'
        ref_loaded  = load(ref_file, 'raw_field_dose');
        ref_rf      = ref_loaded.raw_field_dose;
        ref_dose    = ref_rf.dose_Gy;
        ref_origin  = ref_rf.origin(:);                  % [x; y; z] in mm
        ref_spacing = ref_rf.spacing(:);                 % [dx; dy; dz] in mm
        ref_dims    = size(ref_dose);                    % [rows, cols, slices]
    otherwise  % 'dicom'
        ref_info = dicominfo(ref_file);
        ref_dose = double(squeeze(dicomread(ref_file)));
        % Apply dose grid scaling
        if isfield(ref_info, 'DoseGridScaling')
            ref_dose = ref_dose * ref_info.DoseGridScaling;
        end
        % CRITICAL: Z-resolution from GridFrameOffsetVector, NOT PixelSpacing
        ref_origin  = ref_info.ImagePositionPatient(:);  % [x, y, z] in mm
        ref_spacing = extractDoseSpacing(ref_info);      % [dx, dy, dz] in mm
        ref_dims    = size(ref_dose);                    % [rows, cols, slices]
end

fprintf('  Reference dose grid:\n');
fprintf('    Dimensions: [%d, %d, %d]\n', ref_dims(1), ref_dims(2), ref_dims(3));
fprintf('    Spacing (mm): [%.3f, %.3f, %.3f]\n', ref_spacing(1), ref_spacing(2), ref_spacing(3));
fprintf('    Origin (mm): [%.3f, %.3f, %.3f]\n', ref_origin(1), ref_origin(2), ref_origin(3));

% Initialize total dose accumulator
total_rs_dose = zeros(ref_dims);

% Per-CT-label accumulators; keys added dynamically from filenames (e.g. 'CT_1', 'CT_3')
ct_dose_accum = struct();

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
metadata.beam_metadata = beam_metadata;  % Includes isocenter + jaw data for sensor placement

%% ======================== DISCOVER, LOAD, AND RESAMPLE CBCT1 / CBCT3 ========================

if cbct_ready
    fprintf('\n[4/8] CBCT*_resampled.mat already exist — loading cached resampled CBCTs.\n');

    cbct1_cache = load(fullfile(processed_dir, 'CBCT1_resampled.mat'), 'CBCT1_resampled');
    cbct3_cache = load(fullfile(processed_dir, 'CBCT3_resampled.mat'), 'CBCT3_resampled');
    CBCT1_cached = cbct1_cache.CBCT1_resampled;
    CBCT3_cached = cbct3_cache.CBCT3_resampled;

    ct1_hu_resampled = CBCT1_cached.cubeHU;
    ct1_density      = CBCT1_cached.cubeDensity;
    ct1_dims         = CBCT1_cached.original_cbct_dims;
    ct1_spacing      = CBCT1_cached.original_cbct_spacing;
    cbct_ct1_meta    = struct('series_uid', CBCT1_cached.series_uid, ...
                              'datetime',   CBCT1_cached.series_datetime, ...
                              'label',      'CT_1');

    ct3_hu_resampled = CBCT3_cached.cubeHU;
    ct3_density      = CBCT3_cached.cubeDensity;
    ct3_dims         = CBCT3_cached.original_cbct_dims;
    ct3_spacing      = CBCT3_cached.original_cbct_spacing;
    cbct_ct3_meta    = struct('series_uid', CBCT3_cached.series_uid, ...
                              'datetime',   CBCT3_cached.series_datetime, ...
                              'label',      'CT_3');

    fprintf('    CT_1 cached: HU=[%.0f %.0f]  density=[%.0f %.0f] kg/m^3\n', ...
        min(ct1_hu_resampled(:)), max(ct1_hu_resampled(:)), ...
        min(ct1_density(:)), max(ct1_density(:)));
    fprintf('    CT_3 cached: HU=[%.0f %.0f]  density=[%.0f %.0f] kg/m^3\n', ...
        min(ct3_hu_resampled(:)), max(ct3_hu_resampled(:)), ...
        min(ct3_density(:)), max(ct3_density(:)));

    clear cbct1_cache cbct3_cache CBCT1_cached CBCT3_cached;
else
    fprintf('\n[4/8] Discovering CBCTs and RTSTRUCTs in RayStation directory...\n');

    [cbct_ct1_meta, cbct_ct3_meta] = discoverCbctSeries(rs_dir);

    fprintf('\n  Loading and resampling CBCT1 (CT_1, earlier) to dose grid...\n');
    [ct1_hu, ct1_origin, ct1_spacing, ct1_dims] = loadCbctImagesFromFiles(cbct_ct1_meta.files);
    if isempty(ct1_hu)
        error('step15_process_doses:NoCBCT', 'Failed to load CBCT1 slices.');
    end
    fprintf('    Original CT_1: dims=[%d %d %d]  spacing=[%.3f %.3f %.3f]  origin=[%.3f %.3f %.3f]  HU=[%.0f %.0f]\n', ...
        ct1_dims(1), ct1_dims(2), ct1_dims(3), ct1_spacing(1), ct1_spacing(2), ct1_spacing(3), ...
        ct1_origin(1), ct1_origin(2), ct1_origin(3), min(ct1_hu(:)), max(ct1_hu(:)));
    ct1_hu_resampled = resampleSctToDoseGrid(ct1_hu, ct1_origin, ct1_spacing, ct1_dims, ...
        ref_origin, ref_spacing, ref_dims);
    ct1_density = huToDensity(ct1_hu_resampled);
    fprintf('    CT_1 resampled HU range: [%.0f, %.0f]  density range: [%.0f, %.0f] kg/m^3\n', ...
        min(ct1_hu_resampled(:)), max(ct1_hu_resampled(:)), min(ct1_density(:)), max(ct1_density(:)));

    fprintf('\n  Loading and resampling CBCT3 (CT_3, later) to dose grid...\n');
    [ct3_hu, ct3_origin, ct3_spacing, ct3_dims] = loadCbctImagesFromFiles(cbct_ct3_meta.files);
    if isempty(ct3_hu)
        error('step15_process_doses:NoCBCT', 'Failed to load CBCT3 slices.');
    end
    fprintf('    Original CT_3: dims=[%d %d %d]  spacing=[%.3f %.3f %.3f]  origin=[%.3f %.3f %.3f]  HU=[%.0f %.0f]\n', ...
        ct3_dims(1), ct3_dims(2), ct3_dims(3), ct3_spacing(1), ct3_spacing(2), ct3_spacing(3), ...
        ct3_origin(1), ct3_origin(2), ct3_origin(3), min(ct3_hu(:)), max(ct3_hu(:)));
    ct3_hu_resampled = resampleSctToDoseGrid(ct3_hu, ct3_origin, ct3_spacing, ct3_dims, ...
        ref_origin, ref_spacing, ref_dims);
    ct3_density = huToDensity(ct3_hu_resampled);
    fprintf('    CT_3 resampled HU range: [%.0f, %.0f]  density range: [%.0f, %.0f] kg/m^3\n', ...
        min(ct3_hu_resampled(:)), max(ct3_hu_resampled(:)), min(ct3_density(:)), max(ct3_density(:)));

    % Free raw CBCT volumes; we only need the dose-grid-resampled cubes from here on.
    clear ct1_hu ct3_hu;
end

%% ======================== LOAD RTSTRUCTS AND CREATE PER-CBCT MASKS ========================

tissue_masks_file = fullfile(processed_dir, 'tissue_masks.mat');

% When we used the cached CBCT cubes we never called discoverCbctSeries and
% therefore don't have an .rtstruct path. Fill it in now if we are about to
% rebuild the masks.
if cbct_ready && ~masks_ready && (~isfield(cbct_ct1_meta, 'rtstruct') || isempty(cbct_ct1_meta.rtstruct))
    fprintf('  Cached CBCT cubes lack RTSTRUCT path — running discoverCbctSeries...\n');
    [tmp_meta1, tmp_meta3] = discoverCbctSeries(rs_dir);
    cbct_ct1_meta.rtstruct = tmp_meta1.rtstruct;
    cbct_ct3_meta.rtstruct = tmp_meta3.rtstruct;
    clear tmp_meta1 tmp_meta3;
end

if masks_ready
    fprintf('\n[5/8] tissue_masks.mat already exists — loading cached masks.\n');
    [tissue_mask_ct1, roi_names_ct1, roi_masks_ct1, body_mask_ct1, couch_mask_ct1, ...
     tissue_mask_ct3, roi_names_ct3, roi_masks_ct3, body_mask_ct3, couch_mask_ct3] = ...
        loadCachedTissueMasks(tissue_masks_file, ref_dims);
    fprintf('    CT_1 ROIs=%d  body voxels=%d  couch voxels=%d\n', ...
        length(roi_names_ct1), sum(body_mask_ct1(:)), sum(couch_mask_ct1(:)));
    fprintf('    CT_3 ROIs=%d  body voxels=%d  couch voxels=%d\n', ...
        length(roi_names_ct3), sum(body_mask_ct3(:)), sum(couch_mask_ct3(:)));
else
    fprintf('\n[5/8] Loading per-CBCT RTSTRUCTs and creating tissue classification masks...\n');

    fprintf('  CT_1 RTSTRUCT...\n');
    [tissue_mask_ct1, roi_names_ct1, roi_masks_ct1, body_mask_ct1, couch_mask_ct1] = ...
        loadRtstructAndCreateMasksFromFile(cbct_ct1_meta.rtstruct, ref_origin, ref_spacing, ref_dims);
    if isempty(tissue_mask_ct1)
        error('step15_process_doses:NoRTSTRUCT', ...
            'Could not create CT_1 masks from RTSTRUCT: %s', cbct_ct1_meta.rtstruct);
    end
    fprintf('    CT_1 ROIs=%d  body voxels=%d  couch voxels=%d\n', ...
        length(roi_names_ct1), sum(body_mask_ct1(:)), sum(couch_mask_ct1(:)));

    fprintf('  CT_3 RTSTRUCT...\n');
    [tissue_mask_ct3, roi_names_ct3, roi_masks_ct3, body_mask_ct3, couch_mask_ct3] = ...
        loadRtstructAndCreateMasksFromFile(cbct_ct3_meta.rtstruct, ref_origin, ref_spacing, ref_dims);
    if isempty(tissue_mask_ct3)
        error('step15_process_doses:NoRTSTRUCT', ...
            'Could not create CT_3 masks from RTSTRUCT: %s', cbct_ct3_meta.rtstruct);
    end
    fprintf('    CT_3 ROIs=%d  body voxels=%d  couch voxels=%d\n', ...
        length(roi_names_ct3), sum(body_mask_ct3(:)), sum(couch_mask_ct3(:)));

    % Save both mask sets in a single tissue_masks.mat
    fprintf('  Saving tissue_masks.mat (both CT_1 and CT_3)...\n');
if config.use_sparse_storage
    tissue_mask_ct1_dims = size(tissue_mask_ct1);
    body_mask_ct1_dims   = size(body_mask_ct1);
    couch_mask_ct1_dims  = size(couch_mask_ct1);
    tissue_mask_ct1_sp   = sparse(reshape(double(tissue_mask_ct1), [], tissue_mask_ct1_dims(end)));
    body_mask_ct1_sp     = sparse(reshape(double(body_mask_ct1),   [], body_mask_ct1_dims(end)));
    couch_mask_ct1_sp    = sparse(reshape(double(couch_mask_ct1),  [], couch_mask_ct1_dims(end)));

    tissue_mask_ct3_dims = size(tissue_mask_ct3);
    body_mask_ct3_dims   = size(body_mask_ct3);
    couch_mask_ct3_dims  = size(couch_mask_ct3);
    tissue_mask_ct3_sp   = sparse(reshape(double(tissue_mask_ct3), [], tissue_mask_ct3_dims(end)));
    body_mask_ct3_sp     = sparse(reshape(double(body_mask_ct3),   [], body_mask_ct3_dims(end)));
    couch_mask_ct3_sp    = sparse(reshape(double(couch_mask_ct3),  [], couch_mask_ct3_dims(end)));

    save(tissue_masks_file, ...
        'tissue_mask_ct1_sp', 'tissue_mask_ct1_dims', ...
        'body_mask_ct1_sp',   'body_mask_ct1_dims', ...
        'couch_mask_ct1_sp',  'couch_mask_ct1_dims', ...
        'roi_names_ct1', 'roi_masks_ct1', ...
        'tissue_mask_ct3_sp', 'tissue_mask_ct3_dims', ...
        'body_mask_ct3_sp',   'body_mask_ct3_dims', ...
        'couch_mask_ct3_sp',  'couch_mask_ct3_dims', ...
        'roi_names_ct3', 'roi_masks_ct3', ...
        '-v7.3');
else
    save(tissue_masks_file, ...
        'tissue_mask_ct1', 'roi_names_ct1', 'roi_masks_ct1', 'body_mask_ct1', 'couch_mask_ct1', ...
        'tissue_mask_ct3', 'roi_names_ct3', 'roi_masks_ct3', 'body_mask_ct3', 'couch_mask_ct3', ...
        '-v7.3');
    end
    fprintf('  Saved: tissue_masks.mat\n');
end  % if masks_ready / else

% Pre-compute invalid dose mask PER CBCT (selected per field via ct_label)
if config.apply_dose_masking
    invalid_dose_mask_ct1 = ~(body_mask_ct1 & ~couch_mask_ct1);
    invalid_dose_mask_ct3 = ~(body_mask_ct3 & ~couch_mask_ct3);
    fprintf('  Invalid-dose masks: CT_1 zeros %d voxels, CT_3 zeros %d voxels\n', ...
        sum(invalid_dose_mask_ct1(:)), sum(invalid_dose_mask_ct3(:)));
else
    invalid_dose_mask_ct1 = false(ref_dims);
    invalid_dose_mask_ct3 = false(ref_dims);
    fprintf('  Dose masking DISABLED (debugging mode)\n');
end

%% ======================== PROCESS EACH FIELD DOSE ========================

fprintf('\n[6/8] Processing field doses (masking before export, batch_size=%d)...\n', config.batch_size);

% Track which files were processed successfully
field_doses    = cell(num_files, 1);
processed_count = 0;
save_count      = 0;

num_batches = ceil(num_files / config.batch_size);
fprintf('  Total files: %d | Batches: %d\n', num_files, num_batches);

for batch_idx = 1:num_batches
    batch_start = (batch_idx - 1) * config.batch_size + 1;
    batch_end   = min(batch_idx * config.batch_size, num_files);

    fprintf('\n  --- Batch %d/%d (files %d-%d) ---\n', ...
        batch_idx, num_batches, batch_start, batch_end);

    % Accumulate dose contribution for this batch only, then fold into total
    batch_total_dose = zeros(ref_dims);

    for i = batch_start:batch_end
        fprintf('  Processing field %d/%d: %s\n', i, num_files, rd_files(i).name);

        try
            % ----- Skip if the processed output for this field already exists.
            % If totals also need to be rebuilt we still load the cached dose
            % and fold it into the running totals; otherwise we only populate
            % the lightweight field_doses{i} entry.
            cached_out = expected_field_outputs{i};
            if config.skip_completed && ~isempty(cached_out) && isfile(cached_out)
                cached = load(cached_out, 'field_dose');
                cfd    = cached.field_dose;

                if need_total_accum
                    if isfield(cfd, 'is_sparse') && cfd.is_sparse
                        dense_dose = reshape(full(cfd.dose_Gy), cfd.dose_dims);
                    else
                        dense_dose = cfd.dose_Gy;
                    end
                    batch_total_dose = batch_total_dose + dense_dose;
                    if isfield(cfd, 'ct_label') && ~isempty(cfd.ct_label)
                        ck = strrep(cfd.ct_label, '-', '_');
                        if ~isfield(ct_dose_accum, ck)
                            ct_dose_accum.(ck) = zeros(ref_dims);
                        end
                        ct_dose_accum.(ck) = ct_dose_accum.(ck) + dense_dose;
                    end
                    clear dense_dose;
                end

                field_doses{i} = struct( ...
                    'filepath',     cached_out, ...
                    'beam_num',     cfd.beam_num, ...
                    'seg_num',      cfd.seg_num, ...
                    'field_num',    cfd.field_num, ...
                    'plan_type',    cfd.plan_type, ...
                    'gantry_angle', cfd.gantry_angle, ...
                    'meterset',     cfd.meterset, ...
                    'max_dose_Gy',  cfd.max_dose_Gy, ...
                    'source_file',  cfd.source_file, ...
                    'isocenter',    cfd.isocenter, ...
                    'jaw_x',        cfd.jaw_x, ...
                    'jaw_y',        cfd.jaw_y, ...
                    'body_masked',  cfd.body_masked, ...
                    'couch_masked', cfd.couch_masked);

                processed_count = processed_count + 1;
                fprintf('    [Skip] %s already processed (max: %.4f Gy)\n', ...
                    rd_files(i).name, cfd.max_dose_Gy);
                clear cached cfd;
                continue;
            end

            dose_file = fullfile(rs_dir, rd_files(i).name);
            switch input_format
                case 'mat'
                    % Load pre-converted .mat (from step14_npz_to_mat). Geometry
                    % is already in mm and the dose array is already double in
                    % MATLAB (row=Y, col=X, slice=Z) order.
                    loaded       = load(dose_file, 'raw_field_dose');
                    rf           = loaded.raw_field_dose;
                    dose_data    = rf.dose_Gy;
                    dose_origin  = rf.origin(:);
                    dose_spacing = rf.spacing(:);
                    dose_dims    = size(dose_data);
                otherwise  % 'dicom'
                    % Load DICOM dose file
                    dose_info = dicominfo(dose_file);
                    dose_data = double(squeeze(dicomread(dose_file)));

                    % Apply dose grid scaling
                    if isfield(dose_info, 'DoseGridScaling')
                        dose_scaling = dose_info.DoseGridScaling;
                        dose_data = dose_data * dose_scaling;
                        fprintf('    Applied scaling: %e\n', dose_scaling);
                    end

                    % Extract geometry and verify it matches reference
                    dose_origin  = dose_info.ImagePositionPatient(:);
                    dose_spacing = extractDoseSpacing(dose_info);
                    dose_dims    = size(dose_data);
            end

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

            % Extract beam info from filename
            % e.g. dose_1885729_Session_4_adapted_B6_103.dcm
            %   -> beam_num=6, seg_num=103, field_num=6, plan_type='adapted'
            % field_num (= beam_num) is used to match to RTPLAN beam metadata
            [beam_num, seg_num, field_num, plan_type] = extractBeamInfo(rd_files(i).name, i);

            % Extract CT label (e.g. 'CT_1' / 'CT_3'). Required: every field
            % dose must carry a CT label so we can pick the matching CBCT
            % geometry + mask for simulation and masking.
            ct_tokens = regexp(rd_files(i).name, ...
                '_(?:adapted|reference)_(CT_\d+)_B\d+_\d+\.(?:dcm|mat)$', ...
                'tokens', 'once', 'ignorecase');
            if ~isempty(ct_tokens)
                ct_label = ct_tokens{1};
            else
                ct_label = '';
            end

            % Select per-CBCT invalid-dose mask. Reject anything that doesn't
            % map to a supported CBCT — there is no sensible default now that
            % masks are per-CT.
            switch ct_label
                case 'CT_1'
                    invalid_for_field = invalid_dose_mask_ct1;
                case 'CT_3'
                    invalid_for_field = invalid_dose_mask_ct3;
                otherwise
                    error('step15_process_doses:UnsupportedCtLabel', ...
                        'Field %s has ct_label="%s"; only CT_1 and CT_3 are supported.', ...
                        rd_files(i).name, ct_label);
            end

            % Get beam metadata by matching field_num to beam_number in RTPLAN
            [gantry_angle, meterset] = getBeamMetadata(beam_metadata, field_num);

            % Propagate isocenter and jaw data from beam_metadata
            [iso, jx, jy] = getBeamGeometry(beam_metadata, field_num);

            % Zero out invalid regions (outside body or in couch) BEFORE saving
            if config.apply_dose_masking
                dose_data(invalid_for_field) = 0;
            end

            % Create field dose structure with masking already applied
            field_dose = struct();
            field_dose.dose_Gy = dose_data;
            field_dose.origin = dose_origin;
            field_dose.spacing = dose_spacing;
            field_dose.dimensions = dose_dims;
            field_dose.beam_num = beam_num;         % Beam number from filename (B[n])
            field_dose.seg_num = seg_num;           % Segment number from filename
            field_dose.field_num = field_num;       % Field number (= beam_num, matches RTPLAN)
            field_dose.plan_type = plan_type;       % 'adapted' or 'reference'
            field_dose.ct_label  = ct_label;        % '' for legacy DICOM, 'CT_n' for NPZ-derived
            field_dose.gantry_angle = gantry_angle;
            field_dose.meterset = meterset;
            field_dose.source_file = rd_files(i).name;
            field_dose.max_dose_Gy = max(dose_data(:));
            field_dose.mean_dose_Gy = mean(dose_data(dose_data > 0));
            field_dose.isocenter = iso;
            field_dose.jaw_x = jx;
            field_dose.jaw_y = jy;
            field_dose.body_masked = config.apply_dose_masking;
            field_dose.couch_masked = config.apply_dose_masking;

            % Accumulate into batch subtotal (masking already applied)
            batch_total_dose = batch_total_dose + dose_data;

            % Accumulate per-CT-label total (NPZ-derived inputs only)
            if ~isempty(ct_label)
                ct_key = strrep(ct_label, '-', '_');  % ensure valid struct field name
                if ~isfield(ct_dose_accum, ct_key)
                    ct_dose_accum.(ct_key) = zeros(ref_dims);
                end
                ct_dose_accum.(ct_key) = ct_dose_accum.(ct_key) + dose_data;
            end

            % Save individual field dose file — name mirrors source.
            % Format (legacy DICOM):  dose_[id]_[session]_[plan_type]_B[beam]_[seg].mat
            % Format (NPZ-derived):   dose_[id]_[session]_[plan_type]_[ct_label]_B[beam]_[seg].mat
            if isempty(ct_label)
                field_filename = sprintf('dose_%s_%s_%s_B%d_%d.mat', ...
                    patient_id, session, plan_type, beam_num, seg_num);
            else
                field_filename = sprintf('dose_%s_%s_%s_%s_B%d_%d.mat', ...
                    patient_id, session, plan_type, ct_label, beam_num, seg_num);
            end
            field_filepath = fullfile(processed_dir, field_filename);

            % Convert 3D dose to sparse 2D [nRows*nCols, nSlices] before saving.
            % Field doses are mostly zero outside the treated volume after masking.
            % Reconstruct: reshape(full(field_dose.dose_Gy), field_dose.dose_dims)
            if config.use_sparse_storage
                field_dose.dose_dims = size(field_dose.dose_Gy);
                field_dose.dose_Gy   = sparse(reshape(field_dose.dose_Gy, [], field_dose.dose_dims(end)));
                field_dose.is_sparse = true;
            end
            save(field_filepath, 'field_dose', '-v7.3');
            save_count = save_count + 1;
            fprintf('    [Save %d] %s (B%d S%d [%s], max: %.4f Gy, gantry: %.1f deg, MU: %.1f)\n', ...
                save_count, field_filename, beam_num, seg_num, plan_type, ...
                field_dose.max_dose_Gy, gantry_angle, meterset);

            % Store reference in output cell array (without full dose data for memory)
            field_doses{i} = struct();
            field_doses{i}.filepath = field_filepath;
            field_doses{i}.beam_num = beam_num;
            field_doses{i}.seg_num = seg_num;
            field_doses{i}.field_num = field_num;
            field_doses{i}.plan_type = plan_type;
            field_doses{i}.gantry_angle = gantry_angle;
            field_doses{i}.meterset = meterset;
            field_doses{i}.max_dose_Gy = field_dose.max_dose_Gy;
            field_doses{i}.source_file = rd_files(i).name;
            field_doses{i}.isocenter = iso;
            field_doses{i}.jaw_x = jx;
            field_doses{i}.jaw_y = jy;
            field_doses{i}.body_masked = config.apply_dose_masking;
            field_doses{i}.couch_masked = config.apply_dose_masking;

            processed_count = processed_count + 1;

            % Clear per-file variables immediately to free memory
            clear field_dose dose_data dose_info dose_origin dose_spacing dose_dims;

        catch ME
            warning('step15_process_doses:FieldProcessingError', ...
                'Failed to process field %d (%s): %s', i, rd_files(i).name, ME.message);
        end
    end  % inner file loop

    % Fold batch subtotal into running total, then clear batch arrays
    total_rs_dose = total_rs_dose + batch_total_dose;
    clear batch_total_dose;

    fprintf('  [Batch %d/%d] Done. Running total max: %.4f Gy. Memory cleared.\n', ...
        batch_idx, num_batches, max(total_rs_dose(:)));
end  % batch loop

fprintf('  Successfully processed %d/%d field doses\n', processed_count, num_files);
fprintf('  Total dose max: %.4f Gy\n', max(total_rs_dose(:)));

% Update metadata
metadata.processed_count = processed_count;
metadata.total_dose_max_Gy_before_masking = max(total_rs_dose(:));

%% ======================== ZERO OUT DOSE OUTSIDE BODY AND IN COUCH ========================

if ~need_total_accum
    fprintf('\n[7/8] Total dose / metadata already on disk — loading for return.\n');

    td = load(fullfile(processed_dir, 'total_rs_dose.mat'));
    if isfield(td, 'total_rs_dose_sparse')
        total_rs_dose = reshape(full(td.total_rs_dose_sparse), td.total_rs_dose_dims);
    else
        total_rs_dose = td.total_rs_dose;
    end
    clear td;

    md = load(fullfile(processed_dir, 'metadata.mat'), 'metadata');
    metadata = md.metadata;
    clear md;

    % Skip the masking + save blocks entirely; jump to CBCT struct construction.
    num_voxels_zeroed = NaN;  %#ok<NASGU>  (kept for summary print compatibility)
else

% Combined-anatomy masks for the cross-CBCT total dose: keep a voxel that
% is body in EITHER CBCT and not couch in BOTH.
body_mask_union  = body_mask_ct1  | body_mask_ct3;
couch_mask_inter = couch_mask_ct1 & couch_mask_ct3;
invalid_dose_mask_union = ~(body_mask_union & ~couch_mask_inter);

if config.apply_dose_masking
    fprintf('\n[7/8] Applying final mask to total dose...\n');

    % Apply union mask to total_rs_dose (aggregates both CBCTs)
    num_voxels_outside_body = sum(~body_mask_union(:));
    num_voxels_in_couch     = sum(couch_mask_inter(:));
    num_voxels_zeroed       = sum(invalid_dose_mask_union(:));

    total_rs_dose(invalid_dose_mask_union) = 0;

    % Apply per-CT mask to each per-CT-label accumulator
    ct_keys = fieldnames(ct_dose_accum);
    for k = 1:numel(ct_keys)
        switch ct_keys{k}
            case 'CT_1'
                ct_dose_accum.(ct_keys{k})(invalid_dose_mask_ct1) = 0;
            case 'CT_3'
                ct_dose_accum.(ct_keys{k})(invalid_dose_mask_ct3) = 0;
            otherwise
                error('step15_process_doses:UnsupportedCtLabel', ...
                    'Accumulator key "%s" is not CT_1 or CT_3.', ct_keys{k});
        end
    end

    fprintf('  Voxels outside body (union): %d\n', num_voxels_outside_body);
    fprintf('  Voxels in couch (intersection): %d\n', num_voxels_in_couch);
    fprintf('  Total voxels zeroed: %d\n', num_voxels_zeroed);
    fprintf('  Total dose max (after masking): %.4f Gy\n', max(total_rs_dose(:)));

    metadata.total_dose_max_Gy = max(total_rs_dose(:));
    metadata.body_voxels_ct1   = sum(body_mask_ct1(:));
    metadata.body_voxels_ct3   = sum(body_mask_ct3(:));
    metadata.couch_voxels_ct1  = sum(couch_mask_ct1(:));
    metadata.couch_voxels_ct3  = sum(couch_mask_ct3(:));
    metadata.voxels_zeroed     = num_voxels_zeroed;
    metadata.dose_masking_applied = true;

else
    fprintf('\n[7/8] Dose masking SKIPPED (config.apply_dose_masking = false)\n');

    metadata.total_dose_max_Gy = max(total_rs_dose(:));
    metadata.body_voxels_ct1   = sum(body_mask_ct1(:));
    metadata.body_voxels_ct3   = sum(body_mask_ct3(:));
    metadata.couch_voxels_ct1  = sum(couch_mask_ct1(:));
    metadata.couch_voxels_ct3  = sum(couch_mask_ct3(:));
    metadata.voxels_zeroed     = 0;
    metadata.dose_masking_applied = false;

    fprintf('  Total dose max (unmasked): %.4f Gy\n', max(total_rs_dose(:)));
    fprintf('  Body voxels CT_1: %d   CT_3: %d\n', sum(body_mask_ct1(:)), sum(body_mask_ct3(:)));
    fprintf('  Couch voxels CT_1: %d   CT_3: %d\n', sum(couch_mask_ct1(:)), sum(couch_mask_ct3(:)));
end

%% ======================== SAVE TOTAL DOSE ========================

fprintf('\n  Saving total_rs_dose.mat...\n');
total_dose_file = fullfile(processed_dir, 'total_rs_dose.mat');
if config.use_sparse_storage
    total_rs_dose_dims   = size(total_rs_dose);
    total_rs_dose_sparse = sparse(reshape(total_rs_dose, [], total_rs_dose_dims(end)));
    save(total_dose_file, 'total_rs_dose_sparse', 'total_rs_dose_dims', '-v7.3');
else
    save(total_dose_file, 'total_rs_dose', '-v7.3');
end
fprintf('  Saved: total_rs_dose.mat\n');

% Save per-CT-label total doses
ct_keys = fieldnames(ct_dose_accum);
if ~isempty(ct_keys)
    fprintf('\n  Saving per-CT total doses (%d CT label(s) found)...\n', numel(ct_keys));
    for k = 1:numel(ct_keys)
        ct_key  = ct_keys{k};
        ct_total = ct_dose_accum.(ct_key);
        ct_dose_file = fullfile(processed_dir, sprintf('total_dose_%s.mat', ct_key));
        if config.use_sparse_storage
            ct_total_dims   = size(ct_total);
            ct_total_sparse = sparse(reshape(ct_total, [], ct_total_dims(end)));
            save(ct_dose_file, 'ct_total_sparse', 'ct_total_dims', '-v7.3');
        else
            save(ct_dose_file, 'ct_total', '-v7.3');
        end
        fprintf('    Saved: total_dose_%s.mat (max: %.4f Gy)\n', ct_key, max(ct_total(:)));
    end
end

end  % if need_total_accum

%% ======================== CREATE PER-CBCT RESAMPLED STRUCTS ========================

CBCT1_resampled = struct();
CBCT1_resampled.cubeHU                = ct1_hu_resampled;
CBCT1_resampled.cubeDensity           = ct1_density;
CBCT1_resampled.tissueMask            = tissue_mask_ct1;
CBCT1_resampled.roiNames              = roi_names_ct1;
CBCT1_resampled.bodyMask              = body_mask_ct1;
CBCT1_resampled.couchMask             = couch_mask_ct1;
CBCT1_resampled.origin                = ref_origin;
CBCT1_resampled.spacing               = ref_spacing;
CBCT1_resampled.dimensions            = ref_dims;
CBCT1_resampled.patient_id            = patient_id;
CBCT1_resampled.session               = session;
CBCT1_resampled.original_cbct_dims    = ct1_dims;
CBCT1_resampled.original_cbct_spacing = ct1_spacing;
CBCT1_resampled.series_uid            = cbct_ct1_meta.series_uid;
CBCT1_resampled.series_datetime       = cbct_ct1_meta.datetime;
CBCT1_resampled.ct_label              = 'CT_1';
CBCT1_resampled.timestamp             = datetime('now');

CBCT3_resampled = struct();
CBCT3_resampled.cubeHU                = ct3_hu_resampled;
CBCT3_resampled.cubeDensity           = ct3_density;
CBCT3_resampled.tissueMask            = tissue_mask_ct3;
CBCT3_resampled.roiNames              = roi_names_ct3;
CBCT3_resampled.bodyMask              = body_mask_ct3;
CBCT3_resampled.couchMask             = couch_mask_ct3;
CBCT3_resampled.origin                = ref_origin;
CBCT3_resampled.spacing               = ref_spacing;
CBCT3_resampled.dimensions            = ref_dims;
CBCT3_resampled.patient_id            = patient_id;
CBCT3_resampled.session               = session;
CBCT3_resampled.original_cbct_dims    = ct3_dims;
CBCT3_resampled.original_cbct_spacing = ct3_spacing;
CBCT3_resampled.series_uid            = cbct_ct3_meta.series_uid;
CBCT3_resampled.series_datetime       = cbct_ct3_meta.datetime;
CBCT3_resampled.ct_label              = 'CT_3';
CBCT3_resampled.timestamp             = datetime('now');

% Return bundle (kept in the 2nd output slot in place of the old sct_resampled)
cbct_resampled = struct('CT_1', CBCT1_resampled, 'CT_3', CBCT3_resampled);

%% ======================== SAVE PER-CBCT RESAMPLED FILES ========================

fprintf('\n[8/8] Saving processed data...\n');

if ~cbct_ready
    cbct1_file = fullfile(processed_dir, 'CBCT1_resampled.mat');
    save(cbct1_file, 'CBCT1_resampled', '-v7.3');
    fprintf('  Saved: CBCT1_resampled.mat\n');

    cbct3_file = fullfile(processed_dir, 'CBCT3_resampled.mat');
    save(cbct3_file, 'CBCT3_resampled', '-v7.3');
    fprintf('  Saved: CBCT3_resampled.mat\n');
else
    fprintf('  CBCT1/CBCT3 .mat already on disk — not re-saving.\n');
end

% Save metadata (only when totals were rebuilt; otherwise the on-disk
% metadata reflects the current totals and we leave it alone)
if need_total_accum
    metadata_file = fullfile(processed_dir, 'metadata.mat');
    save(metadata_file, 'metadata', '-v7.3');
    fprintf('  Saved: metadata.mat\n');
else
    fprintf('  metadata.mat already on disk — not re-saving.\n');
end

%% ======================== SUMMARY ========================

fprintf('\n========================================\n');
fprintf('  Step 1.5 Complete\n');
fprintf('========================================\n');
fprintf('  Processed %d field doses\n', processed_count);
fprintf('  .mat files saved: %d\n', save_count);
ct_keys = fieldnames(ct_dose_accum);
if ~isempty(ct_keys)
    fprintf('  Per-CT total doses saved: %s\n', strjoin(ct_keys, ', '));
end
fprintf('  Dose grid: [%d x %d x %d]\n', ref_dims(1), ref_dims(2), ref_dims(3));
fprintf('  Spacing: [%.3f, %.3f, %.3f] mm\n', ref_spacing(1), ref_spacing(2), ref_spacing(3));
fprintf('  Total dose max: %.4f Gy\n', max(total_rs_dose(:)));
fprintf('  CT_1 tissue ROIs: %d   body voxels: %d   couch voxels: %d\n', ...
    length(roi_names_ct1), sum(body_mask_ct1(:)), sum(couch_mask_ct1(:)));
fprintf('  CT_3 tissue ROIs: %d   body voxels: %d   couch voxels: %d\n', ...
    length(roi_names_ct3), sum(body_mask_ct3(:)), sum(couch_mask_ct3(:)));
if config.apply_dose_masking
    fprintf('  Dose masking: ENABLED (applied before export)\n');
    fprintf('    Total voxels zeroed: %d\n', num_voxels_zeroed);
else
    fprintf('  Dose masking: DISABLED (debugging mode)\n');
end
fprintf('  Output directory: %s\n', processed_dir);
fprintf('========================================\n\n');

end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function s = presentStr(tf)
%PRESENTSTR Human-readable 'present'/'MISSING' label for status logs.
    if tf, s = 'present'; else, s = 'MISSING'; end
end


function [status, expected_outputs] = detectProcessedOutputs(processed_dir, rd_files, patient_id, session)
%DETECTPROCESSEDOUTPUTS Inspect processed/ and report which artifacts exist.
%
%   Returns:
%     status            - struct of booleans + counts
%     expected_outputs  - cell array (1 per rd_files entry) with the expected
%                         per-field .mat path under processed/ (empty if the
%                         input filename can't be parsed)

    status = struct();
    status.cbct1_done      = false;
    status.cbct3_done      = false;
    status.masks_done      = false;
    status.total_rs_done   = false;
    status.ct1_total_done  = false;
    status.ct3_total_done  = false;
    status.metadata_done   = false;
    status.num_fields_done = 0;

    n = numel(rd_files);
    expected_outputs = cell(n, 1);

    if ~isfolder(processed_dir)
        return;
    end

    status.cbct1_done     = isfile(fullfile(processed_dir, 'CBCT1_resampled.mat'));
    status.cbct3_done     = isfile(fullfile(processed_dir, 'CBCT3_resampled.mat'));
    status.masks_done     = isfile(fullfile(processed_dir, 'tissue_masks.mat'));
    status.total_rs_done  = isfile(fullfile(processed_dir, 'total_rs_dose.mat'));
    status.ct1_total_done = isfile(fullfile(processed_dir, 'total_dose_CT_1.mat'));
    status.ct3_total_done = isfile(fullfile(processed_dir, 'total_dose_CT_3.mat'));
    status.metadata_done  = isfile(fullfile(processed_dir, 'metadata.mat'));

    for i = 1:n
        out_path = buildFieldOutputPath(processed_dir, rd_files(i).name, patient_id, session);
        expected_outputs{i} = out_path;
        if ~isempty(out_path) && isfile(out_path)
            status.num_fields_done = status.num_fields_done + 1;
        end
    end
end


function out_path = buildFieldOutputPath(processed_dir, rd_name, patient_id, session)
%BUILDFIELDOUTPUTPATH Reproduce the per-field output filename used by step15.
%
%   Mirrors the same parsing + sprintf format used in the main loop so we
%   can check "was this input already processed?" without re-running it.

    out_path = '';
    [beam_num, seg_num, ~, plan_type] = extractBeamInfo(rd_name, NaN);
    if isnan(beam_num) || isnan(seg_num) || strcmp(plan_type, 'unknown')
        return;
    end

    ct_tokens = regexp(rd_name, ...
        '_(?:adapted|reference)_(CT_\d+)_B\d+_\d+\.(?:dcm|mat)$', ...
        'tokens', 'once', 'ignorecase');
    if isempty(ct_tokens)
        fname = sprintf('dose_%s_%s_%s_B%d_%d.mat', ...
            patient_id, session, plan_type, beam_num, seg_num);
    else
        fname = sprintf('dose_%s_%s_%s_%s_B%d_%d.mat', ...
            patient_id, session, plan_type, ct_tokens{1}, beam_num, seg_num);
    end
    out_path = fullfile(processed_dir, fname);
end


function [field_doses, cbct_resampled, total_rs_dose, metadata] = ...
    loadExistingStep15Outputs(processed_dir, expected_outputs)
%LOADEXISTINGSTEP15OUTPUTS Reconstruct return values from cached processed/.
%
%   Used when every output in processed/ is already present and
%   skip_completed=true. The per-field cell array mirrors the lightweight
%   summary struct that the normal path produces (filepath + beam/jaw
%   metadata, no dose array).

    % Per-field summaries
    n = numel(expected_outputs);
    field_doses = cell(n, 1);
    for i = 1:n
        fpath = expected_outputs{i};
        if isempty(fpath) || ~isfile(fpath), continue; end
        s = load(fpath, 'field_dose');
        fd = s.field_dose;
        field_doses{i} = struct( ...
            'filepath',     fpath, ...
            'beam_num',     fd.beam_num, ...
            'seg_num',      fd.seg_num, ...
            'field_num',    fd.field_num, ...
            'plan_type',    fd.plan_type, ...
            'gantry_angle', fd.gantry_angle, ...
            'meterset',     fd.meterset, ...
            'max_dose_Gy',  fd.max_dose_Gy, ...
            'source_file',  fd.source_file, ...
            'isocenter',    fd.isocenter, ...
            'jaw_x',        fd.jaw_x, ...
            'jaw_y',        fd.jaw_y, ...
            'body_masked',  fd.body_masked, ...
            'couch_masked', fd.couch_masked);
    end

    % CBCT bundle
    c1 = load(fullfile(processed_dir, 'CBCT1_resampled.mat'), 'CBCT1_resampled');
    c3 = load(fullfile(processed_dir, 'CBCT3_resampled.mat'), 'CBCT3_resampled');
    cbct_resampled = struct('CT_1', c1.CBCT1_resampled, 'CT_3', c3.CBCT3_resampled);

    % Total dose (sparse 2D or dense 3D)
    td = load(fullfile(processed_dir, 'total_rs_dose.mat'));
    if isfield(td, 'total_rs_dose_sparse')
        total_rs_dose = reshape(full(td.total_rs_dose_sparse), td.total_rs_dose_dims);
    else
        total_rs_dose = td.total_rs_dose;
    end

    md = load(fullfile(processed_dir, 'metadata.mat'), 'metadata');
    metadata = md.metadata;
end


function [tm_ct1, rn_ct1, rm_ct1, bm_ct1, cm_ct1, ...
          tm_ct3, rn_ct3, rm_ct3, bm_ct3, cm_ct3] = ...
    loadCachedTissueMasks(tissue_masks_file, ref_dims)
%LOADCACHEDTISSUEMASKS Restore CT_1/CT_3 mask variables from tissue_masks.mat
%
%   Handles both storage layouts: sparse 2D ([nVoxPerSlice x nSlices], with
%   the *_dims companion variable) and dense 3D logical arrays.

    s = load(tissue_masks_file);

    % CT_1
    if isfield(s, 'tissue_mask_ct1_sp')
        tm_ct1 = uint8(reshape(full(s.tissue_mask_ct1_sp), s.tissue_mask_ct1_dims));
        bm_ct1 = logical(reshape(full(s.body_mask_ct1_sp),   s.body_mask_ct1_dims));
        cm_ct1 = logical(reshape(full(s.couch_mask_ct1_sp),  s.couch_mask_ct1_dims));
    else
        tm_ct1 = s.tissue_mask_ct1;
        bm_ct1 = s.body_mask_ct1;
        cm_ct1 = s.couch_mask_ct1;
    end
    rn_ct1 = s.roi_names_ct1;
    rm_ct1 = s.roi_masks_ct1;

    % CT_3
    if isfield(s, 'tissue_mask_ct3_sp')
        tm_ct3 = uint8(reshape(full(s.tissue_mask_ct3_sp), s.tissue_mask_ct3_dims));
        bm_ct3 = logical(reshape(full(s.body_mask_ct3_sp),   s.body_mask_ct3_dims));
        cm_ct3 = logical(reshape(full(s.couch_mask_ct3_sp),  s.couch_mask_ct3_dims));
    else
        tm_ct3 = s.tissue_mask_ct3;
        bm_ct3 = s.body_mask_ct3;
        cm_ct3 = s.couch_mask_ct3;
    end
    rn_ct3 = s.roi_names_ct3;
    rm_ct3 = s.roi_masks_ct3;

    % Defensive: verify shapes line up with the live reference grid
    if ~isequal(size(bm_ct1), ref_dims) || ~isequal(size(bm_ct3), ref_dims)
        error('step15_process_doses:CachedMasksShape', ...
            'Cached tissue_masks.mat shape does not match ref_dims [%d %d %d].', ...
            ref_dims(1), ref_dims(2), ref_dims(3));
    end
end


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
%   Extract gantry angles, metersets, isocenter positions, and jaw positions
%   for each beam. Isocenter and jaw data are required by determine_sensor_mask
%   for computing beam field exclusion zones on the patient surface.
%
%   FIELDS EXTRACTED:
%       .beam_number   - Beam number from BeamSequence
%       .beam_name     - Beam name string
%       .gantry_angle  - Gantry angle from first ControlPoint (degrees)
%       .meterset      - Monitor units from FractionGroupSequence
%       .isocenter     - [x, y, z] mm from ControlPointSequence.Item_1.IsocenterPosition
%       .jaw_x         - [x1, x2] mm at isocenter (ASYMX / X jaw positions)
%       .jaw_y         - [y1, y2] mm at isocenter (ASYMY / Y jaw positions)

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
            
            % Initialize defaults for all fields
            beam_metadata(i).gantry_angle = 0;
            beam_metadata(i).meterset = 0;
            beam_metadata(i).isocenter = [];
            beam_metadata(i).jaw_x = [];
            beam_metadata(i).jaw_y = [];
            
            % Extract from first ControlPoint
            if isfield(beam, 'ControlPointSequence')
                cp_fields = fieldnames(beam.ControlPointSequence);
                if ~isempty(cp_fields)
                    cp1 = beam.ControlPointSequence.(cp_fields{1});
                    
                    % Gantry angle
                    if isfield(cp1, 'GantryAngle')
                        beam_metadata(i).gantry_angle = cp1.GantryAngle;
                    end
                    
                    % Isocenter position [x, y, z] in mm (DICOM patient coords)
                    if isfield(cp1, 'IsocenterPosition')
                        beam_metadata(i).isocenter = cp1.IsocenterPosition(:)';
                    end
                    
                    % Jaw positions from BeamLimitingDevicePositionSequence
                    if isfield(cp1, 'BeamLimitingDevicePositionSequence')
                        bld_seq = cp1.BeamLimitingDevicePositionSequence;
                        bld_fields = fieldnames(bld_seq);
                        
                        for j = 1:length(bld_fields)
                            bld_item = bld_seq.(bld_fields{j});
                            
                            if ~isfield(bld_item, 'RTBeamLimitingDeviceType') || ...
                               ~isfield(bld_item, 'LeafJawPositions')
                                continue;
                            end
                            
                            dev_type = upper(bld_item.RTBeamLimitingDeviceType);
                            positions = bld_item.LeafJawPositions;
                            
                            % ASYMX or X jaws -> jaw_x
                            if strcmp(dev_type, 'ASYMX') || strcmp(dev_type, 'X')
                                if length(positions) >= 2
                                    beam_metadata(i).jaw_x = positions(1:2)';
                                end
                            end
                            
                            % ASYMY or Y jaws -> jaw_y
                            if strcmp(dev_type, 'ASYMY') || strcmp(dev_type, 'Y')
                                if length(positions) >= 2
                                    beam_metadata(i).jaw_y = positions(1:2)';
                                end
                            end
                        end
                    end
                end
            end
            
            % Log extracted geometry
            if ~isempty(beam_metadata(i).isocenter)
                fprintf('    Beam %d: gantry=%.1f deg, iso=[%.1f, %.1f, %.1f] mm', ...
                    beam_metadata(i).beam_number, beam_metadata(i).gantry_angle, ...
                    beam_metadata(i).isocenter(1), beam_metadata(i).isocenter(2), ...
                    beam_metadata(i).isocenter(3));
            else
                fprintf('    Beam %d: gantry=%.1f deg, iso=N/A', ...
                    beam_metadata(i).beam_number, beam_metadata(i).gantry_angle);
            end
            
            if ~isempty(beam_metadata(i).jaw_x)
                fprintf(', jaw_x=[%.1f, %.1f]', beam_metadata(i).jaw_x(1), beam_metadata(i).jaw_x(2));
            end
            if ~isempty(beam_metadata(i).jaw_y)
                fprintf(', jaw_y=[%.1f, %.1f]', beam_metadata(i).jaw_y(1), beam_metadata(i).jaw_y(2));
            end
            fprintf('\n');
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


function [beam_num, seg_num, field_num, plan_type] = extractBeamInfo(filename, default_index)
%EXTRACTBEAMINFO Extract beam, segment, field numbers and plan type from dose filename
%
%   Supported patterns (checked in priority order):
%     1. dose_[id]_[session]_(adapted|reference)[_CT_k]_B[n]_[seg].(dcm|mat)
%        e.g. dose_1885729_Session_4_adapted_B6_103.dcm        (legacy DICOM)
%             dose_1194203_Session_1_adapted_CT_1_B13_00.mat   (RayStation NPZ -> mat)
%        -> beam_num  = 6 or 13  (B[n])
%        -> seg_num   = 103 or 0
%        -> field_num = beam_num (matches RTPLAN beam_number)
%        -> plan_type = 'adapted'
%     2. Plan_Field [n]_Beam[m]_B[n]_S[m].dcm  (legacy step06 format)
%     3. Beam[n]_Seg[m]_Field [o].dcm           (legacy RayStation export)
%     4. Beam[n]_Seg[m]_Field[o].dcm            (legacy, no space)
%     5. RD.[n].dcm                              (legacy DICOM default)
%
%   OUTPUTS:
%       beam_num  - Beam number (B[n] from filename)
%       seg_num   - Segment number within that beam
%       field_num - Field number used to match RTPLAN beam metadata (= beam_num)
%       plan_type - 'adapted', 'reference', or 'unknown'

    beam_num  = default_index;
    seg_num   = 0;
    field_num = default_index;
    plan_type = 'unknown';

    % --- Pattern 1 (preferred): dose_*_(adapted|reference)[_CT_k]_B[n]_[seg].(dcm|mat) ---
    tokens = regexp(filename, ...
        '_(adapted|reference)(?:_CT_\d+)?_B(\d+)_(\d+)\.(?:dcm|mat)$', ...
        'tokens', 'ignorecase');
    if ~isempty(tokens) && ~isempty(tokens{1})
        plan_type = lower(tokens{1}{1});        % 'adapted' or 'reference'
        beam_num  = str2double(tokens{1}{2});   % B[n]
        seg_num   = str2double(tokens{1}{3});   % segment number
        field_num = beam_num;                   % field_num matches RTPLAN beam_number
        return;
    end

    % --- Pattern 2: Plan_Field [n]_Beam[m]_B[n]_S[m].dcm ---
    tokens = regexp(filename, 'Plan_Field\s*(\d+)_Beam(\d+)_B(\d+)_S(\d+)', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        field_num = str2double(tokens{1}{1});
        beam_num  = str2double(tokens{1}{3});
        seg_num   = str2double(tokens{1}{4});
        return;
    end

    % --- Pattern 3: Beam[n]_Seg[m]_Field [o].dcm ---
    tokens = regexp(filename, 'Beam(\d+)_Seg(\d+)_Field\s*(\d+)', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        beam_num  = str2double(tokens{1}{1});
        seg_num   = str2double(tokens{1}{2});
        field_num = str2double(tokens{1}{3});
        return;
    end

    % --- Pattern 4: Beam[n]_Seg[m]_Field[o].dcm (no space) ---
    tokens = regexp(filename, 'Beam(\d+)_Seg(\d+)_Field(\d+)', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        beam_num  = str2double(tokens{1}{1});
        seg_num   = str2double(tokens{1}{2});
        field_num = str2double(tokens{1}{3});
        return;
    end

    % --- Pattern 5 (legacy): RD.[n].dcm ---
    tokens = regexp(filename, 'RD\.(\d+)', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        beam_num  = str2double(tokens{1}{1});
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


function [isocenter, jaw_x, jaw_y] = getBeamGeometry(beam_metadata, field_num)
%GETBEAMGEOMETRY Get isocenter and jaw positions for a beam by field number
%
%   Matches field_num to beam_number in beam_metadata struct array.
%   Returns empty arrays if fields are not found.

    isocenter = [];
    jaw_x = [];
    jaw_y = [];
    
    if isempty(beam_metadata) || ~isstruct(beam_metadata)
        return;
    end
    
    % Match field_num to beam_number in RTPLAN
    for j = 1:length(beam_metadata)
        if beam_metadata(j).beam_number == field_num
            if isfield(beam_metadata(j), 'isocenter')
                isocenter = beam_metadata(j).isocenter;
            end
            if isfield(beam_metadata(j), 'jaw_x')
                jaw_x = beam_metadata(j).jaw_x;
            end
            if isfield(beam_metadata(j), 'jaw_y')
                jaw_y = beam_metadata(j).jaw_y;
            end
            return;
        end
    end
    
    % Fallback: use field_num as index if within range
    if field_num <= length(beam_metadata)
        if isfield(beam_metadata(field_num), 'isocenter')
            isocenter = beam_metadata(field_num).isocenter;
        end
        if isfield(beam_metadata(field_num), 'jaw_x')
            jaw_x = beam_metadata(field_num).jaw_x;
        end
        if isfield(beam_metadata(field_num), 'jaw_y')
            jaw_y = beam_metadata(field_num).jaw_y;
        end
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

    % Get all DICOM files whose SeriesDescription is 'sct' (case-insensitive)
    all_files = dir(fullfile(sct_dir, '*.dcm'));
    sct_files = [];
    z_positions = [];

    for i = 1:length(all_files)
        try
            info = dicominfo(fullfile(sct_dir, all_files(i).name));
            if isfield(info, 'SeriesDescription') && strcmpi(info.SeriesDescription, 'sct')
                sct_files = [sct_files; all_files(i)]; %#ok<AGROW>
                z_positions = [z_positions; info.ImagePositionPatient(3)]; %#ok<AGROW>
            end
        catch
            % Skip unreadable files
        end
    end

    if isempty(sct_files)
        warning('loadSctImages:NoFiles', 'No SCT DICOM files found in: %s', sct_dir);
        return;
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


function [cbct_hu, origin, spacing, dims] = loadCbctImagesFromFiles(file_list)
%LOADCBCTIMAGESFROMFILES Load DICOM CT slices from an explicit file list
%
%   Same return contract as loadSctImages but operates on an explicit cell
%   array of DICOM filepaths (e.g. dicomCollection().Filenames{i}) instead
%   of globbing a directory by SeriesDescription.

    cbct_hu = [];
    origin  = [];
    spacing = [];
    dims    = [];

    if iscell(file_list) && ~isempty(file_list) && iscell(file_list{1})
        % Some dicomCollection variants nest filenames in a single cell
        file_list = file_list{1};
    end

    n = numel(file_list);
    if n == 0
        warning('loadCbctImagesFromFiles:NoFiles', 'Empty CBCT file list.');
        return;
    end

    z_positions = nan(n, 1);
    for i = 1:n
        try
            info = dicominfo(file_list{i});
            z_positions(i) = info.ImagePositionPatient(3);
        catch
            % leave NaN; will sort to end
        end
    end

    [~, sort_idx] = sort(z_positions);

    first_info = dicominfo(file_list{sort_idx(1)});
    first_img  = dicomread(file_list{sort_idx(1)});

    cbct_data = zeros([size(first_img), n], 'int16');
    fprintf('    Loading %d CBCT slices...\n', n);
    for i = 1:n
        idx = sort_idx(i);
        cbct_data(:, :, i) = dicomread(file_list{idx});
    end

    origin = first_info.ImagePositionPatient(:);
    dx = first_info.PixelSpacing(1);
    dy = first_info.PixelSpacing(2);

    if n >= 2
        second_info = dicominfo(file_list{sort_idx(2)});
        dz = abs(second_info.ImagePositionPatient(3) - first_info.ImagePositionPatient(3));
    elseif isfield(first_info, 'SliceThickness')
        dz = first_info.SliceThickness;
    else
        dz = dx;
    end

    spacing = [dx; dy; dz];
    dims    = size(cbct_data);

    if isfield(first_info, 'RescaleSlope') && isfield(first_info, 'RescaleIntercept')
        cbct_hu = double(cbct_data) * first_info.RescaleSlope + first_info.RescaleIntercept;
    else
        cbct_hu = double(cbct_data);
    end
end


function [cbct_ct1, cbct_ct3] = discoverCbctSeries(rs_dir)
%DISCOVERCBCTSERIES Find both CBCT image series + RTSTRUCTs in a RayStation dir
%
%   Scans rs_dir with dicomCollection, picks out the two CT series (the
%   CBCTs), pairs each with its RTSTRUCT by SeriesInstanceUID, sorts the
%   pair by SeriesDate+SeriesTime, and returns:
%       cbct_ct1 - earlier-acquired CBCT (treated as CT_1)
%       cbct_ct3 - later-acquired CBCT  (treated as CT_3)
%   Each struct has fields:
%       .files      - cell array of DICOM file paths for the CT series
%       .rtstruct   - filepath to the paired RTSTRUCT
%       .series_uid - SeriesInstanceUID of the CT series
%       .datetime   - numeric SeriesDate+SeriesTime
%       .label      - 'CT_1' or 'CT_3'

    cbct_ct1 = [];
    cbct_ct3 = [];

    try
        info = dicomCollection(rs_dir);
    catch ME
        error('step15_process_doses:CBCTDiscoveryFailed', ...
            'dicomCollection failed on %s: %s', rs_dir, ME.message);
    end

    if isempty(info) || height(info) == 0
        error('step15_process_doses:CBCTDiscoveryFailed', ...
            'No DICOM series found in: %s', rs_dir);
    end

    if ~ismember('Modality', info.Properties.VariableNames)
        error('step15_process_doses:CBCTDiscoveryFailed', ...
            'dicomCollection missing Modality column for: %s', rs_dir);
    end

    ct_rows     = strcmpi(info.Modality, 'CT');
    rtstr_rows  = strcmpi(info.Modality, 'RTSTRUCT');

    ct_info    = info(ct_rows, :);
    rtstr_info = info(rtstr_rows, :);

    if height(ct_info) < 2
        error('step15_process_doses:CBCTDiscoveryFailed', ...
            'Expected at least 2 CT series in %s, found %d.', rs_dir, height(ct_info));
    end

    % Read SeriesInstanceUID + SeriesDate/Time per CT series
    nCt = height(ct_info);
    ct_series_uid = strings(nCt, 1);
    ct_datetime   = nan(nCt, 1);
    ct_files      = cell(nCt, 1);
    for i = 1:nCt
        files_i = ct_info.Filenames{i};
        if iscell(files_i) && ~isempty(files_i) && iscell(files_i{1})
            files_i = files_i{1};
        end
        ct_files{i} = files_i;
        if isempty(files_i)
            continue;
        end
        try
            meta = dicominfo(files_i{1});
            if isfield(meta, 'SeriesInstanceUID')
                ct_series_uid(i) = string(meta.SeriesInstanceUID);
            end
            dateStr = '';
            timeStr = '';
            if isfield(meta, 'SeriesDate'), dateStr = strtrim(meta.SeriesDate); end
            if isfield(meta, 'SeriesTime'), timeStr = strtrim(meta.SeriesTime); end
            if ~isempty(dateStr) || ~isempty(timeStr)
                ct_datetime(i) = str2double([dateStr, timeStr]);
            end
        catch
            % leave as default
        end
    end

    % If more than two CT series are present, keep the two earliest
    if nCt > 2
        [~, dt_sort] = sort(ct_datetime);
        keep_idx = dt_sort(1:2);
        ct_series_uid = ct_series_uid(keep_idx);
        ct_datetime   = ct_datetime(keep_idx);
        ct_files      = ct_files(keep_idx);
        nCt = 2;
        fprintf('    More than 2 CT series found; keeping the two earliest by datetime.\n');
    end

    % Read referenced-SeriesInstanceUID per RTSTRUCT
    nRs = height(rtstr_info);
    rs_ref_uid = strings(nRs, 1);
    rs_path    = strings(nRs, 1);
    for i = 1:nRs
        files_i = rtstr_info.Filenames{i};
        if iscell(files_i) && ~isempty(files_i) && iscell(files_i{1})
            files_i = files_i{1};
        end
        if isempty(files_i)
            continue;
        end
        rs_path(i) = string(files_i{1});
        try
            meta = dicominfo(files_i{1});
            ref = meta.ReferencedFrameOfReferenceSequence.Item_1 ...
                      .RTReferencedStudySequence.Item_1 ...
                      .RTReferencedSeriesSequence.Item_1 ...
                      .SeriesInstanceUID;
            rs_ref_uid(i) = string(ref);
        catch
            % leave empty; will fail the match check below
        end
    end

    % Pair each CT series with the matching RTSTRUCT
    matched_rtstruct = strings(nCt, 1);
    for i = 1:nCt
        match = find(rs_ref_uid == ct_series_uid(i), 1);
        if isempty(match)
            error('step15_process_doses:CBCTDiscoveryFailed', ...
                'No RTSTRUCT references CT SeriesInstanceUID %s in %s.', ...
                ct_series_uid(i), rs_dir);
        end
        matched_rtstruct(i) = rs_path(match);
    end

    % Sort the two CBCTs by datetime: earlier -> CT_1, later -> CT_3
    [~, sort_order] = sort(ct_datetime);
    earlier = sort_order(1);
    later   = sort_order(2);

    cbct_ct1 = struct( ...
        'files',      {ct_files{earlier}}, ...
        'rtstruct',   char(matched_rtstruct(earlier)), ...
        'series_uid', char(ct_series_uid(earlier)), ...
        'datetime',   ct_datetime(earlier), ...
        'label',      'CT_1');

    cbct_ct3 = struct( ...
        'files',      {ct_files{later}}, ...
        'rtstruct',   char(matched_rtstruct(later)), ...
        'series_uid', char(ct_series_uid(later)), ...
        'datetime',   ct_datetime(later), ...
        'label',      'CT_3');

    fprintf('    CT_1 (earlier): UID=%s  datetime=%s  files=%d  RTSTRUCT=%s\n', ...
        cbct_ct1.series_uid, num2str(cbct_ct1.datetime), ...
        numel(cbct_ct1.files), cbct_ct1.rtstruct);
    fprintf('    CT_3 (later):   UID=%s  datetime=%s  files=%d  RTSTRUCT=%s\n', ...
        cbct_ct3.series_uid, num2str(cbct_ct3.datetime), ...
        numel(cbct_ct3.files), cbct_ct3.rtstruct);
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
%HUTODENSITY Convert Hounsfield Units to density (kg/mÂ³)
%
%   Uses simplified linear conversion:
%   - Below -1000 HU (air): density = 1 kg/mÂ³
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
    % Approximate: density â‰ˆ 1000 + HU
    soft_mask = (hu >= -500) & (hu < 100);
    density(soft_mask) = 1000 + hu(soft_mask);
    
    % Bone region (HU >= 100)
    bone_mask = hu >= 100;
    % Linear interpolation from soft tissue to dense bone
    density(bone_mask) = 1100 + (hu(bone_mask) - 100) * (1900 - 1100) / 900;
    
    % Clamp to reasonable range
    density = max(density, 1);      % Minimum 1 kg/mÂ³
    density = min(density, 7800);   % Maximum (metal)
end


function [tissue_mask, roi_names, roi_masks, body_mask, couch_mask] = loadRtstructAndCreateMasks(...
    sct_dir, dose_origin, dose_spacing, dose_dims)
%LOADRTSTRUCTANDCREATEMASKS Find RTSTRUCT in a directory and create masks
%
%   Thin wrapper around loadRtstructAndCreateMasksFromFile: globs
%   RTSTRUCT*.dcm (then RS*.dcm) in sct_dir and uses the first match.

    tissue_mask = [];
    roi_names   = {};
    roi_masks   = struct();
    body_mask   = false(dose_dims);
    couch_mask  = false(dose_dims);

    rs_files = dir(fullfile(sct_dir, 'RTSTRUCT*.dcm'));
    if isempty(rs_files)
        rs_files = dir(fullfile(sct_dir, 'RS*.dcm'));
    end
    if isempty(rs_files)
        warning('loadRtstructAndCreateMasks:NoRTSTRUCT', ...
            'No RTSTRUCT file found in: %s', sct_dir);
        return;
    end

    rs_path = fullfile(sct_dir, rs_files(1).name);
    [tissue_mask, roi_names, roi_masks, body_mask, couch_mask] = ...
        loadRtstructAndCreateMasksFromFile(rs_path, dose_origin, dose_spacing, dose_dims);
end


function [tissue_mask, roi_names, roi_masks, body_mask, couch_mask] = loadRtstructAndCreateMasksFromFile(...
    rs_path, dose_origin, dose_spacing, dose_dims)
%LOADRTSTRUCTANDCREATEMASKSFROMFILE Load a specific RTSTRUCT file and create masks
%
%   Same outputs as loadRtstructAndCreateMasks but takes an explicit
%   RTSTRUCT filepath rather than a directory to glob.

    tissue_mask = [];
    roi_names = {};
    roi_masks = struct();
    body_mask = false(dose_dims);
    couch_mask = false(dose_dims);

    if isempty(rs_path) || ~isfile(rs_path)
        warning('loadRtstructAndCreateMasksFromFile:NoRTSTRUCT', ...
            'RTSTRUCT file not found: %s', rs_path);
        return;
    end

    [~, rs_name, rs_ext] = fileparts(rs_path);
    fprintf('    Loading RTSTRUCT: %s%s\n', rs_name, rs_ext);

    try
        rtstruct = dicominfo(rs_path);
    catch ME
        warning('loadRtstructAndCreateMasksFromFile:LoadError', ...
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
    fprintf('    Body mask (before gap fill): %d voxels\n', sum(body_mask(:)));
    fprintf('    Couch mask (before gap fill): %d voxels\n', sum(couch_mask(:)));
    
    % ===== FIX GAPS IN BODY AND COUCH MASKS =====
    % RTSTRUCT contours may not align with dose grid slices, causing gaps.
    % Fill gaps by interpolating along z-direction.
    fprintf('    Filling z-direction gaps in masks...\n');
    body_mask = fillMaskZGaps(body_mask);
    couch_mask = fillMaskZGaps(couch_mask);
    fprintf('    Body mask (after gap fill): %d voxels\n', sum(body_mask(:)));
    fprintf('    Couch mask (after gap fill): %d voxels\n', sum(couch_mask(:)));
end


function mask_filled = fillMaskZGaps(mask)
%FILLMASKZGAPS Fill gaps in a 3D binary mask along the z-direction
%
%   For each (x,y) column, finds gaps (zeros between ones) and fills them.
%   This fixes the common issue where RTSTRUCT contour z-coordinates don't
%   align with dose grid slice positions, causing alternating filled/empty slices.
%
%   Uses vectorized operations for efficiency on large volumes.
%
%   INPUT:
%       mask - 3D logical array
%
%   OUTPUT:
%       mask_filled - 3D logical array with z-gaps filled

    [nRows, nCols, nSlices] = size(mask);
    
    % Reshape to 2D: each column is a z-profile for one (row,col) position
    % Shape: [nSlices x (nRows*nCols)]
    mask_2d = reshape(mask, [], nSlices)';  % [nSlices x nPixels]
    
    % For each pixel column, find first and last true indices
    % We'll use cumsum tricks for vectorization
    
    % Method: For each column, create a mask that is true from first_true to last_true
    % 1. Forward cumsum: marks everything from first true onward
    % 2. Backward cumsum: marks everything up to last true
    % 3. AND them together
    
    % Forward pass: cumsum > 0 means we've seen at least one true
    forward_cumsum = cumsum(mask_2d, 1);
    from_first = forward_cumsum > 0;
    
    % Backward pass: flip, cumsum, flip back
    backward_cumsum = flipud(cumsum(flipud(mask_2d), 1));
    to_last = backward_cumsum > 0;
    
    % Fill: true where we're between first and last (inclusive)
    mask_2d_filled = from_first & to_last;
    
    % Reshape back to 3D
    mask_filled = reshape(mask_2d_filled', nRows, nCols, nSlices);
end
