function [field_doses, sct_resampled, total_rs_dose, metadata] = load_processed_data(patient_id, session, config)
%% LOAD_PROCESSED_DATA - Load previously processed field doses and SCT data
%
%   [field_doses, sct_resampled, total_rs_dose, metadata] = load_processed_data(patient_id, session, config)
%
%   PURPOSE:
%   Load all processed data files created by step15_process_doses from the
%   processed/ directory. Returns data in the same format as the processing
%   function for seamless pipeline integration.
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir    - Base directory path
%
%   OUTPUTS:
%       field_doses     - Cell array of field dose structures, each containing:
%           .dose_Gy       - 3D dose array in Gy
%           .origin        - [x, y, z] in mm
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%           .beam_index    - Beam number
%           .gantry_angle  - Gantry angle in degrees
%           .meterset      - Monitor units
%           .source_file   - Original DICOM filename
%       sct_resampled   - Struct with CT resampled to dose grid:
%           .cubeHU        - 3D HU array
%           .cubeDensity   - 3D density array (kg/m³)
%           .origin        - [x, y, z] in mm
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%       total_rs_dose   - Sum of all field doses (3D array in Gy)
%       metadata        - Struct with combined geometry info
%
%   EXPECTED FILES (in processed/ directory):
%       - field_dose_001.mat, field_dose_002.mat, ... (one per field)
%       - sct_resampled.mat
%       - total_rs_dose.mat
%       - metadata.mat
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       [field_doses, sct, total_dose, meta] = load_processed_data('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       - Processed data from step15_process_doses
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: step15_process_doses, step2_kwave_simulation

%% ======================== INPUT VALIDATION ========================

% Validate patient_id
if ~ischar(patient_id) && ~isstring(patient_id)
    error('load_processed_data:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

% Validate session
if ~ischar(session) && ~isstring(session)
    error('load_processed_data:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

% Validate config struct
if ~isstruct(config)
    error('load_processed_data:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

% Validate required config fields
if ~isfield(config, 'working_dir')
    error('load_processed_data:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

%% ======================== CONSTRUCT PATHS ========================

fprintf('\n[load_processed_data] Loading processed data for %s / %s\n', patient_id, session);

% Processed data directory
processed_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session, 'processed');

%% ======================== VERIFY DIRECTORY EXISTS ========================

if ~isfolder(processed_dir)
    error('load_processed_data:DirectoryNotFound', ...
        'Processed directory does not exist: %s\nRun step15_process_doses first.', processed_dir);
end

%% ======================== LOAD METADATA ========================

fprintf('  Loading metadata...\n');

metadata_file = fullfile(processed_dir, 'metadata.mat');

if ~isfile(metadata_file)
    error('load_processed_data:FileNotFound', ...
        'metadata.mat not found in: %s', processed_dir);
end

loaded = load(metadata_file);
metadata = loaded.metadata;

fprintf('    Grid dimensions: [%d, %d, %d]\n', ...
    metadata.dimensions(1), metadata.dimensions(2), metadata.dimensions(3));
fprintf('    Number of fields: %d\n', metadata.num_fields);

%% ======================== LOAD TOTAL DOSE ========================

fprintf('  Loading total_rs_dose...\n');

total_dose_file = fullfile(processed_dir, 'total_rs_dose.mat');

if ~isfile(total_dose_file)
    error('load_processed_data:FileNotFound', ...
        'total_rs_dose.mat not found in: %s', processed_dir);
end

loaded = load(total_dose_file);
total_rs_dose = loaded.total_rs_dose;

fprintf('    Total dose max: %.4f Gy\n', max(total_rs_dose(:)));

%% ======================== LOAD SCT RESAMPLED ========================

fprintf('  Loading sct_resampled...\n');

sct_file = fullfile(processed_dir, 'sct_resampled.mat');

if ~isfile(sct_file)
    error('load_processed_data:FileNotFound', ...
        'sct_resampled.mat not found in: %s', processed_dir);
end

loaded = load(sct_file);
sct_resampled = loaded.sct_resampled;

fprintf('    SCT HU range: [%.0f, %.0f]\n', ...
    min(sct_resampled.cubeHU(:)), max(sct_resampled.cubeHU(:)));
fprintf('    SCT density range: [%.0f, %.0f] kg/m³\n', ...
    min(sct_resampled.cubeDensity(:)), max(sct_resampled.cubeDensity(:)));

%% ======================== LOAD FIELD DOSES ========================

fprintf('  Loading field doses...\n');

% Find all field_dose_XXX.mat files
field_files = dir(fullfile(processed_dir, 'field_dose_*.mat'));

if isempty(field_files)
    error('load_processed_data:FileNotFound', ...
        'No field_dose_*.mat files found in: %s', processed_dir);
end

num_files = length(field_files);
fprintf('    Found %d field dose files\n', num_files);

% Sort files by index number to ensure correct order
file_indices = zeros(num_files, 1);
for i = 1:num_files
    tokens = regexp(field_files(i).name, 'field_dose_(\d+)\.mat', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        file_indices(i) = str2double(tokens{1}{1});
    else
        file_indices(i) = i;
    end
end

[~, sort_order] = sort(file_indices);

% Initialize cell array
field_doses = cell(num_files, 1);

% Load each field dose file
loaded_count = 0;
for i = 1:num_files
    idx = sort_order(i);
    field_filepath = fullfile(processed_dir, field_files(idx).name);
    
    try
        loaded = load(field_filepath);
        field_doses{i} = loaded.field_dose;
        loaded_count = loaded_count + 1;
        
        fprintf('    [%d] %s (gantry: %.1f°, max: %.4f Gy)\n', ...
            i, field_files(idx).name, ...
            field_doses{i}.gantry_angle, field_doses{i}.max_dose_Gy);
        
    catch ME
        warning('load_processed_data:LoadError', ...
            'Failed to load %s: %s', field_files(idx).name, ME.message);
        field_doses{i} = [];
    end
end

%% ======================== SUMMARY ========================

fprintf('\n[load_processed_data] Complete\n');
fprintf('  Loaded %d/%d field doses\n', loaded_count, num_files);
fprintf('  Grid: [%d x %d x %d], Spacing: [%.2f, %.2f, %.2f] mm\n', ...
    metadata.dimensions(1), metadata.dimensions(2), metadata.dimensions(3), ...
    metadata.spacing(1), metadata.spacing(2), metadata.spacing(3));

end
