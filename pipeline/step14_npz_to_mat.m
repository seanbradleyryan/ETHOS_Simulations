function n_converted = step14_npz_to_mat(patient_id, session, config)
%STEP14_NPZ_TO_MAT Convert RayStation NPZ field doses into per-field .mat files.
%
%   n_converted = step14_npz_to_mat(patient_id, session, config)
%
%   PURPOSE:
%   Bridges the RayStation NPZ export (calc_beam_plan_doses.py) and the
%   MATLAB pipeline. Each NPZ in
%       <config.working_dir>/RayStationFiles/<patient_id>/<session>/
%   is unpacked and converted into a per-field .mat file with the same
%   stem (extension swapped). The .mat files are the input format expected
%   by the .mat-branch in step15_process_doses.m.
%
%   INPUTS:
%       patient_id - char/string, patient identifier (e.g. '1194203')
%       session    - char/string, session name (e.g. 'Session_1')
%       config     - struct with at least:
%           .working_dir  - base pipeline directory (char)
%
%   OUTPUT:
%       n_converted - number of NPZ files successfully converted
%
%   NPZ SCHEMA (from calc_beam_plan_doses.py):
%       dose          : float32, shape (nz, ny, nx)  -- C-order
%       voxel_size_cm : float32, shape (3,)          -- [vx, vy, vz] in cm
%       corner_cm     : float32, shape (3,)          -- [x,  y,  z] in cm
%       nx, ny, nz    : int32 scalars
%
%   OUTPUT .MAT SCHEMA (variable name: raw_field_dose):
%       .dose_Gy    : double, shape [ny, nx, nz]  (MATLAB (row=Y, col=X, slice=Z))
%       .origin     : double 3x1 column vector, mm (DICOM patient coords)
%       .spacing    : double 3x1 column vector, mm [dx; dy; dz]
%       .dimensions : double row vector [ny, nx, nz]
%       .ct_label   : char, parsed from filename (e.g. 'CT_1') or ''
%       .source_npz : original NPZ filename
%
%   See also: step15_process_doses, pipeline_compress

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('step14_npz_to_mat:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

if ~ischar(session) && ~isstring(session)
    error('step14_npz_to_mat:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

if ~isstruct(config)
    error('step14_npz_to_mat:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end
if ~isfield(config, 'working_dir')
    error('step14_npz_to_mat:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

% Skip per-file when its target .mat already exists. Set false to force
% full re-conversion.
if ~isfield(config, 'skip_completed')
    config.skip_completed = true;
end

%% ======================== LOCATE NPZ FILES ========================

rs_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session);

fprintf('\n========================================\n');
fprintf('  Step 1.4: Convert NPZ Field Doses to .mat\n');
fprintf('  Patient: %s, Session: %s\n', patient_id, session);
fprintf('========================================\n');
fprintf('  Raystation directory: %s\n', rs_dir);

if ~isfolder(rs_dir)
    error('step14_npz_to_mat:DirectoryNotFound', ...
        'Raystation directory does not exist: %s', rs_dir);
end

npz_files = dir(fullfile(rs_dir, 'dose_*.npz'));
if isempty(npz_files)
    fprintf('  No dose_*.npz files found. Nothing to convert.\n');
    n_converted = 0;
    return;
end

fprintf('  Found %d NPZ file(s)\n', numel(npz_files));

%% ======================== CONVERT EACH NPZ ========================

n_converted = 0;
n_skipped   = 0;
for i = 1:numel(npz_files)
    npz_name = npz_files(i).name;
    npz_path = fullfile(rs_dir, npz_name);
    [~, stem, ~] = fileparts(npz_name);
    mat_name = [stem '.mat'];
    mat_path = fullfile(rs_dir, mat_name);

    if config.skip_completed && isfile(mat_path)
        fprintf('  [%d/%d] %s -> %s already exists, skipping.\n', ...
            i, numel(npz_files), npz_name, mat_name);
        n_skipped = n_skipped + 1;
        continue;
    end

    fprintf('  [%d/%d] %s\n', i, numel(npz_files), npz_name);

    tmp_dir = '';
    try
        tmp_dir = tempname;
        mkdir(tmp_dir);
        extracted = local_unpack_npz(npz_path, tmp_dir);

        dose_npy_path     = local_find_member(extracted, 'dose.npy');
        voxel_npy_path    = local_find_member(extracted, 'voxel_size_cm.npy');
        corner_npy_path   = local_find_member(extracted, 'corner_cm.npy');

        dose_npz       = local_read_npy(dose_npy_path);       % (nz, ny, nx)
        voxel_size_cm  = local_read_npy(voxel_npy_path);      % [vx; vy; vz]
        corner_cm      = local_read_npy(corner_npy_path);     % [x;  y;  z ]

        % NPZ stores dose in C-order (nz, ny, nx). MATLAB step15 expects
        % (rows=Y, cols=X, slices=Z) = (ny, nx, nz).
        dose_mat = double(permute(dose_npz, [2 3 1]));

        % Convert geometry: cm -> mm.
        spacing_mm = double(voxel_size_cm(:)) * 10;   % 3x1 column
        origin_mm  = double(corner_cm(:))     * 10;   % 3x1 column
        dims       = size(dose_mat);                  % [ny, nx, nz]

        % Best-effort extract of CT label from filename
        ct_label = local_extract_ct_label(npz_name);

        raw_field_dose = struct();
        raw_field_dose.dose_Gy    = dose_mat;
        raw_field_dose.origin     = origin_mm;
        raw_field_dose.spacing    = spacing_mm;
        raw_field_dose.dimensions = dims;
        raw_field_dose.ct_label   = ct_label;
        raw_field_dose.source_npz = npz_name;

        save(mat_path, 'raw_field_dose', '-v7.3');

        fprintf('    -> %s  (size=[%d %d %d], spacing=[%.3f %.3f %.3f] mm)\n', ...
            mat_name, dims(1), dims(2), dims(3), ...
            spacing_mm(1), spacing_mm(2), spacing_mm(3));

        n_converted = n_converted + 1;
    catch ME
        warning('step14_npz_to_mat:ConvertFailed', ...
            'Failed to convert %s: %s', npz_name, ME.message);
    end

    if ~isempty(tmp_dir) && isfolder(tmp_dir)
        rmdir(tmp_dir, 's');
    end
end

fprintf('  Converted %d/%d NPZ file(s). Skipped (already converted): %d\n', ...
    n_converted, numel(npz_files), n_skipped);

end


%% =========================================================================
%  PRIVATE HELPERS
%% =========================================================================

function extracted = local_unpack_npz(npz_path, dest_dir)
%LOCAL_UNPACK_NPZ Extract an NPZ file (which is a ZIP of .npy members).
    extracted = unzip(npz_path, dest_dir);
    if isempty(extracted)
        error('local_unpack_npz:Empty', 'No members extracted from %s', npz_path);
    end
end


function p = local_find_member(extracted, name)
%LOCAL_FIND_MEMBER Find an extracted member path by basename.
    for k = 1:numel(extracted)
        [~, base, ext] = fileparts(extracted{k});
        if strcmpi([base ext], name)
            p = extracted{k};
            return;
        end
    end
    error('local_find_member:NotFound', 'NPZ missing member: %s', name);
end


function ct_label = local_extract_ct_label(npz_name)
%LOCAL_EXTRACT_CT_LABEL Best-effort parse of CT_n label from NPZ filename.
%
%   Filenames follow: dose_{id}_{session}_{plan_type}_{ct_label}_{origbeam}_{seg:02d}.npz
%   Plan_type is adapted|reference; we capture the segment that follows it
%   before B<digits>.
    ct_label = '';
    tokens = regexp(npz_name, '_(?:adapted|reference)_(.+?)_B\d+_\d+\.npz$', ...
        'tokens', 'once', 'ignorecase');
    if ~isempty(tokens)
        ct_label = tokens{1};
    end
end


function arr = local_read_npy(npy_path)
%LOCAL_READ_NPY Minimal NPY 1.0/2.0 reader (sufficient for RayStation NPZ).
%
%   Supports dtypes: <f4, <f8, <i4, <i8, <u4, <u1, |i1, |u1
%   Returns a MATLAB array with the logical shape declared in the NPY header
%   (C-order arrays are permuted back to the original shape after reshape).

    fid = fopen(npy_path, 'rb');
    if fid < 0
        error('local_read_npy:OpenFail', 'Cannot open %s', npy_path);
    end
    cleaner = onCleanup(@() fclose(fid));

    magic = fread(fid, 6, '*uint8')';
    if ~isequal(magic, uint8([147 78 85 77 80 89]))   % \x93NUMPY
        error('local_read_npy:BadMagic', 'Not a valid NPY file: %s', npy_path);
    end

    ver_major = fread(fid, 1, 'uint8');
    fread(fid, 1, 'uint8');  % minor (ignored)

    if ver_major == 1
        header_len = double(fread(fid, 1, 'uint16'));
    elseif ver_major >= 2
        header_len = double(fread(fid, 1, 'uint32'));
    else
        error('local_read_npy:BadVersion', 'Unsupported NPY version %d', ver_major);
    end

    header = char(fread(fid, header_len, '*uint8')');

    descr         = local_parse_descr(header);
    fortran_order = local_parse_fortran(header);
    shape         = local_parse_shape(header);

    switch descr
        case '<f4',          precision = 'single';
        case '<f8',          precision = 'double';
        case '<i4',          precision = 'int32';
        case '<i8',          precision = 'int64';
        case '<u4',          precision = 'uint32';
        case {'<u1','|u1'},  precision = 'uint8';
        case '|i1',          precision = 'int8';
        otherwise
            error('local_read_npy:UnsupportedDtype', ...
                'Unsupported dtype "%s" in %s', descr, npy_path);
    end

    if isempty(shape)
        n_elements = 1;
    else
        n_elements = prod(shape);
    end

    data = fread(fid, n_elements, ['*' precision]);

    if isempty(shape)
        arr = data;
    elseif numel(shape) == 1
        arr = data;  % 1-D array; column vector is fine for our consumers
    elseif fortran_order
        arr = reshape(data, shape);
    else
        % C-order on disk: reshape with reversed dims, then permute back.
        arr = permute(reshape(data, flip(shape)), numel(shape):-1:1);
    end
end


function descr = local_parse_descr(header)
    t = regexp(header, '''descr''\s*:\s*''([^'']+)''', 'tokens', 'once');
    if isempty(t)
        error('local_parse_descr:Missing', 'descr field not found in NPY header');
    end
    descr = t{1};
end


function fo = local_parse_fortran(header)
    t = regexp(header, '''fortran_order''\s*:\s*(True|False)', 'tokens', 'once');
    if isempty(t)
        error('local_parse_fortran:Missing', 'fortran_order field not found in NPY header');
    end
    fo = strcmp(t{1}, 'True');
end


function sh = local_parse_shape(header)
    t = regexp(header, '''shape''\s*:\s*\(([^)]*)\)', 'tokens', 'once');
    if isempty(t)
        error('local_parse_shape:Missing', 'shape field not found in NPY header');
    end
    inner = strtrim(t{1});
    if isempty(inner)
        sh = [];
        return;
    end
    parts = strsplit(inner, ',');
    sh = [];
    for k = 1:numel(parts)
        p = strtrim(parts{k});
        if isempty(p)
            continue;
        end
        sh(end+1) = str2double(p); %#ok<AGROW>
    end
end
