function [field_index, metadata] = list_processed_field_doses(patient_id, session, config)
%LIST_PROCESSED_FIELD_DOSES Enumerate processed field-dose files (no heavy load).
%
%   [field_index, metadata] = list_processed_field_doses(patient_id, session, config)
%
%   PURPOSE:
%   Build a lightweight index of the processed dose_*.mat files for a
%   patient/session WITHOUT loading the large dose_Gy arrays. The simulation
%   loop uses this to load one dose per worker on demand (via
%   load_field_dose_file) instead of holding every dose in a cell array and
%   broadcasting it into every parfor worker. Only small per-field metadata
%   stays resident.
%
%   The returned struct array is sorted by beam then segment, matching the
%   ordering used by load_processed_data, so field indices remain stable.
%
%   INPUTS:
%       patient_id - Char/string, patient identifier (e.g., '1194203').
%       session    - Char/string, session name (e.g., 'Session_1').
%       config     - Struct with:
%           .working_dir - Base directory path.
%
%   OUTPUTS:
%       field_index - 1xN struct array (one per dose file), each with:
%           .file                - Full path to the dose_*.mat file
%           .source_mat_filename - Basename (with extension), for output naming
%           .beam_index          - Beam number parsed from _B(\d+)_
%           .segment             - Segment number parsed from _(\d+)\.mat
%           .gantry_angle        - From metadata.beam_metadata matched by
%                                  beam_number, else NaN (cosmetic; the
%                                  authoritative value comes from the loaded
%                                  dose inside the worker)
%       metadata    - Struct loaded from metadata.mat (dimensions, spacing,
%                     num_fields, beam_metadata, ...).
%
%   ALGORITHM:
%       1. Load metadata.mat (cheap).
%       2. dir() the dose_*.mat files; parse beam/segment from each name.
%       3. Sort by [beam, segment]; assemble the struct array in that order.
%       4. Attach gantry_angle from beam_metadata when matchable.
%
%   EXAMPLE:
%       [fidx, meta] = list_processed_field_doses('1194203', 'Session_1', config);
%       fd = load_field_dose_file(fidx(1).file);
%
%   DEPENDENCIES:
%       - Processed data from step15_process_doses (metadata.mat + dose_*.mat).
%
%   See also: load_processed_data, load_field_dose_file, pipeline_simulate

    %% ======================== INPUT VALIDATION ========================

    if ~ischar(patient_id) && ~isstring(patient_id)
        error('list_processed_field_doses:InvalidInput', ...
            'patient_id must be a string or character array. Received: %s', ...
            class(patient_id));
    end
    patient_id = char(patient_id);

    if ~ischar(session) && ~isstring(session)
        error('list_processed_field_doses:InvalidInput', ...
            'session must be a string or character array. Received: %s', ...
            class(session));
    end
    session = char(session);

    if ~isstruct(config) || ~isfield(config, 'working_dir')
        error('list_processed_field_doses:MissingConfig', ...
            'config must be a struct containing working_dir.');
    end

    %% ======================== CONSTRUCT PATHS ========================

    processed_dir = fullfile(config.working_dir, 'RayStationFiles', ...
        patient_id, session, 'processed');

    if ~isfolder(processed_dir)
        error('list_processed_field_doses:DirectoryNotFound', ...
            ['Processed directory does not exist: %s\n' ...
             'Run step15_process_doses first.'], processed_dir);
    end

    %% ======================== LOAD METADATA ========================

    metadata_file = fullfile(processed_dir, 'metadata.mat');
    if ~isfile(metadata_file)
        error('list_processed_field_doses:FileNotFound', ...
            'metadata.mat not found in: %s', processed_dir);
    end
    loaded   = load(metadata_file);
    metadata = loaded.metadata;

    %% ======================== ENUMERATE + SORT DOSE FILES ========================

    field_files = dir(fullfile(processed_dir, 'dose_*.mat'));
    if isempty(field_files)
        error('list_processed_field_doses:FileNotFound', ...
            'No dose_*.mat files found in: %s', processed_dir);
    end

    num_files = numel(field_files);
    beam_nums = zeros(num_files, 1);
    seg_nums  = zeros(num_files, 1);
    for i = 1:num_files
        tokens = regexp(field_files(i).name, '_B(\d+)_(\d+)\.mat$', 'tokens');
        if ~isempty(tokens) && ~isempty(tokens{1})
            beam_nums(i) = str2double(tokens{1}{1});
            seg_nums(i)  = str2double(tokens{1}{2});
        else
            beam_nums(i) = i;
            seg_nums(i)  = 0;
        end
    end

    [~, sort_order] = sortrows([beam_nums, seg_nums]);

    %% ======================== BUILD LIGHTWEIGHT INDEX ========================

    % Optional gantry lookup table from beam_metadata (matched by beam_number).
    have_beam_meta = isfield(metadata, 'beam_metadata') && ~isempty(metadata.beam_metadata);

    field_index = repmat(struct( ...
        'file', '', 'source_mat_filename', '', ...
        'beam_index', NaN, 'segment', NaN, 'gantry_angle', NaN), 1, num_files);

    for i = 1:num_files
        idx  = sort_order(i);
        name = field_files(idx).name;

        field_index(i).file                = fullfile(processed_dir, name);
        field_index(i).source_mat_filename = name;
        field_index(i).beam_index          = beam_nums(idx);
        field_index(i).segment             = seg_nums(idx);
        field_index(i).gantry_angle        = lookup_gantry( ...
            have_beam_meta, metadata, beam_nums(idx));
    end

    fprintf(['[list_processed_field_doses] Indexed %d field dose file(s) for ' ...
             '%s / %s (no dose arrays loaded).\n'], num_files, patient_id, session);
end


function gantry = lookup_gantry(have_beam_meta, metadata, beam_number)
%LOOKUP_GANTRY Best-effort gantry angle for a beam number from beam_metadata.
%  Returns NaN when no match is available; the real per-field value is read
%  from the loaded dose inside the worker.
    gantry = NaN;
    if ~have_beam_meta || isnan(beam_number)
        return;
    end
    bm = metadata.beam_metadata;
    for k = 1:numel(bm)
        if isfield(bm(k), 'beam_number') && ~isempty(bm(k).beam_number) && ...
                bm(k).beam_number == beam_number
            if isfield(bm(k), 'gantry_angle') && ~isempty(bm(k).gantry_angle)
                gantry = bm(k).gantry_angle;
            end
            return;
        end
    end
end
