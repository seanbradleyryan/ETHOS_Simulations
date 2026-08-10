function fd = load_field_dose_file(file_path)
%LOAD_FIELD_DOSE_FILE Load a single processed field-dose .mat from disk.
%
%   fd = load_field_dose_file(file_path)
%
%   PURPOSE:
%   Load exactly ONE processed dose file (dose_*.mat produced by
%   step15_process_doses) and return its field_dose struct with a dense
%   dose_Gy array and source_mat_filename populated. This lets a parfor worker
%   load just its own dose on demand instead of broadcasting the entire
%   field_doses cell array (which copies all doses into every worker and blows
%   up host RAM). Mirrors the per-file loading block in load_processed_data so
%   both paths share one implementation of sparse reconstruction + naming.
%
%   INPUTS:
%       file_path - Char/string, absolute path to a dose_*.mat file containing
%                   a variable named 'field_dose'. The field_dose struct may be
%                   stored in sparse 2D form (field_dose.is_sparse == true with
%                   field_dose.dose_dims giving the original [nx, ny, nz]).
%
%   OUTPUTS:
%       fd - field_dose struct with:
%           .dose_Gy             - 3D dose array (Gy), densified if needed
%           .source_mat_filename - basename of file_path (with extension), used
%                                  downstream to name derived outputs
%           (all other fields from the saved struct passed through unchanged)
%
%   ALGORITHM:
%       1. load() the .mat and pull out field_dose.
%       2. If saved sparse, reshape full(dose_Gy) back to dose_dims.
%       3. Set source_mat_filename from the file basename.
%
%   EXAMPLE:
%       fd = load_field_dose_file('/mnt/.../processed/dose_1194203_Session_1_adapted_B1_1.mat');
%
%   DEPENDENCIES:
%       - Processed dose files from step15_process_doses / pipeline_compress.
%
%   See also: load_processed_data, list_processed_field_doses,
%             run_single_field_simulation

    if ~ischar(file_path) && ~isstring(file_path)
        error('load_field_dose_file:InvalidInput', ...
            'file_path must be a string or character array. Received: %s', ...
            class(file_path));
    end
    file_path = char(file_path);

    if ~isfile(file_path)
        error('load_field_dose_file:FileNotFound', ...
            'Field dose file not found: %s', file_path);
    end

    loaded = load(file_path);
    if ~isfield(loaded, 'field_dose')
        error('load_field_dose_file:NoFieldDose', ...
            'File does not contain a ''field_dose'' variable: %s', file_path);
    end
    fd = loaded.field_dose;

    % Reconstruct dense 3D dose_Gy if it was saved in sparse 2D format.
    if isfield(fd, 'is_sparse') && fd.is_sparse
        fd.dose_Gy = reshape(full(fd.dose_Gy), fd.dose_dims);
    end

    % Expose the source .mat filename so downstream consumers can name derived
    % outputs after the input (e.g. <basename>_sim.mat).
    [~, nm, ext] = fileparts(file_path);
    fd.source_mat_filename = [nm, ext];
end
