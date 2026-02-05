function [adjusted_rtplan_path, num_corrections] = step05_fix_mlc_gaps(patient_id, session, config)
%% STEP05_FIX_MLC_GAPS - Correct small MLC gaps in RTPLAN for Halcyon dual-layer MLC
%
%   [adjusted_rtplan_path, num_corrections] = step05_fix_mlc_gaps(patient_id, session, config)
%
%   PURPOSE:
%   Correct small MLC leaf gaps (<0.5mm) in RTPLAN files from ETHOS exports
%   for Halcyon LINAC with dual-layer MLC. Small gaps can cause Raystation
%   import failures. This function expands gaps symmetrically to minimum
%   threshold while respecting leaf position boundaries.
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir          - Base directory path
%           .treatment_site       - Subfolder name (default: 'Pancreas')
%           .mlc_min_gap_mm       - Minimum allowed gap (default: 0.5)
%           .mlc_expansion_mm     - Expansion per side (default: 0.4)
%           .mlc_position_range   - Valid leaf position range (default: [-140, 140])
%           .mlc_position_tol     - Tolerance for dynamic leaf detection (default: 1e-6)
%
%   OUTPUTS:
%       adjusted_rtplan_path - String, full path to corrected RTPLAN file
%                              (*_adjusted_mlc.dcm) or empty if failed
%       num_corrections      - Integer, total number of leaf gap corrections made
%
%   ALGORITHM:
%   1. Find RTPLAN in sct directory
%   2. Read DICOM metadata
%   3. For each beam:
%      a. Identify MLC devices (MLCX, MLCX1, MLCX2)
%      b. Classify leaves as dynamic/static based on position changes
%      c. For each control point, for each dynamic leaf pair:
%         - Calculate gap = bank_b - bank_a
%         - If gap < threshold: expand symmetrically
%         - Handle boundary conditions (clamp to valid range)
%   4. Generate new SOPInstanceUID
%   5. Update RTPlanLabel with '_adj' suffix
%   6. Write modified DICOM
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       config.mlc_min_gap_mm = 0.5;
%       config.mlc_expansion_mm = 0.4;
%       [path, n] = step05_fix_mlc_gaps('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (dicominfo, dicomwrite)
%       - dicomuid() for generating new UIDs
%
%   NOTES:
%       - Halcyon has dual-layer MLC (MLCX1 and MLCX2)
%       - Bank A = negative side (leaves 1:N), Bank B = positive side (leaves N+1:2N)
%       - Gap = Bank B position - Bank A position (positive when open)
%       - Only dynamic leaves (those that move during delivery) are corrected
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 2.0 (Refactored for master pipeline integration)
%
%   See also: step0_sort_dicom, dicominfo, dicomwrite, dicomuid

%% ======================== INPUT VALIDATION ========================

% Initialize outputs
adjusted_rtplan_path = '';
num_corrections = 0;

% Validate patient_id
if ~ischar(patient_id) && ~isstring(patient_id)
    error('step05_fix_mlc_gaps:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

% Validate session
if ~ischar(session) && ~isstring(session)
    error('step05_fix_mlc_gaps:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

% Validate config struct
if ~isstruct(config)
    error('step05_fix_mlc_gaps:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

% Validate required config fields
if ~isfield(config, 'working_dir')
    error('step05_fix_mlc_gaps:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

% Validate working directory exists
if ~isfolder(config.working_dir)
    error('step05_fix_mlc_gaps:DirectoryNotFound', ...
        'Working directory does not exist: %s', config.working_dir);
end

%% ======================== SET DEFAULT CONFIG VALUES ========================

% Set default treatment_site
if ~isfield(config, 'treatment_site') || isempty(config.treatment_site)
    config.treatment_site = 'Pancreas';
    fprintf('  [INFO] Using default treatment_site: %s\n', config.treatment_site);
end

% Set default MLC parameters
if ~isfield(config, 'mlc_min_gap_mm') || isempty(config.mlc_min_gap_mm)
    config.mlc_min_gap_mm = 0.5;
    fprintf('  [INFO] Using default mlc_min_gap_mm: %.2f mm\n', config.mlc_min_gap_mm);
end

if ~isfield(config, 'mlc_expansion_mm') || isempty(config.mlc_expansion_mm)
    config.mlc_expansion_mm = 0.4;
    fprintf('  [INFO] Using default mlc_expansion_mm: %.2f mm\n', config.mlc_expansion_mm);
end

if ~isfield(config, 'mlc_position_range') || isempty(config.mlc_position_range)
    config.mlc_position_range = [-140, 140];
    fprintf('  [INFO] Using default mlc_position_range: [%.0f, %.0f] mm\n', ...
        config.mlc_position_range(1), config.mlc_position_range(2));
end

if ~isfield(config, 'mlc_position_tol') || isempty(config.mlc_position_tol)
    config.mlc_position_tol = 1e-6;
end

% Extract parameters for cleaner code
MIN_GAP_THRESHOLD = config.mlc_min_gap_mm;
EXPANSION_AMOUNT = config.mlc_expansion_mm;
LEAF_POS_MIN = config.mlc_position_range(1);
LEAF_POS_MAX = config.mlc_position_range(2);
POSITION_TOLERANCE = config.mlc_position_tol;

%% ======================== CONSTRUCT PATHS ========================

% Build path to sct directory containing sorted DICOM files
sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session, 'sct');

fprintf('  Processing: Patient %s, %s\n', patient_id, session);
fprintf('  SCT directory: %s\n', sct_dir);

%% ======================== VERIFY SCT DIRECTORY ========================

if ~isfolder(sct_dir)
    warning('step05_fix_mlc_gaps:DirectoryNotFound', ...
        'SCT directory not found: %s\nRun step0_sort_dicom first.', sct_dir);
    return;
end

%% ======================== FIND RTPLAN FILE ========================

fprintf('  Searching for RTPLAN file...\n');

rtplan_file = findRtplanFile(sct_dir);

if isempty(rtplan_file)
    warning('step05_fix_mlc_gaps:NoRTPLAN', ...
        'No RTPLAN file found in: %s', sct_dir);
    return;
end

rtplan_filepath = fullfile(sct_dir, rtplan_file);
fprintf('  Found RTPLAN: %s\n', rtplan_file);

%% ======================== CHECK IF ALREADY PROCESSED ========================

% Check if an adjusted version already exists
[~, rtplan_name, ~] = fileparts(rtplan_file);
adjusted_filename = [rtplan_name '_adjusted_mlc.dcm'];
adjusted_rtplan_path = fullfile(sct_dir, adjusted_filename);

if exist(adjusted_rtplan_path, 'file')
    fprintf('  [INFO] Adjusted RTPLAN already exists: %s\n', adjusted_filename);
    fprintf('  [INFO] To reprocess, delete existing file first.\n');
    
    % Count existing corrections from the file (optional - read and report)
    % For now, return with path but 0 corrections (already done)
    return;
end

%% ======================== READ RTPLAN DICOM ========================

fprintf('  Reading RTPLAN DICOM metadata...\n');

try
    rtplan = dicominfo(rtplan_filepath);
catch ME
    error('step05_fix_mlc_gaps:DicomReadError', ...
        'Failed to read RTPLAN DICOM: %s\nError: %s', rtplan_filepath, ME.message);
end

% Verify this is an RTPLAN
if ~isfield(rtplan, 'Modality') || ~strcmpi(rtplan.Modality, 'RTPLAN')
    error('step05_fix_mlc_gaps:InvalidModality', ...
        'File is not an RTPLAN (Modality: %s)', rtplan.Modality);
end

% Verify BeamSequence exists
if ~isfield(rtplan, 'BeamSequence')
    error('step05_fix_mlc_gaps:NoBeamSequence', ...
        'RTPLAN does not contain BeamSequence');
end

%% ======================== GET BEAM INFORMATION ========================

beam_fields = fieldnames(rtplan.BeamSequence);
num_beams = length(beam_fields);

fprintf('  Number of beams: %d\n', num_beams);
fprintf('  MLC correction parameters:\n');
fprintf('    - Minimum gap threshold: %.2f mm\n', MIN_GAP_THRESHOLD);
fprintf('    - Expansion per side: %.2f mm\n', EXPANSION_AMOUNT);
fprintf('    - Valid position range: [%.0f, %.0f] mm\n', LEAF_POS_MIN, LEAF_POS_MAX);

%% ======================== PROCESS EACH BEAM ========================

total_corrections = 0;
beam_corrections = zeros(num_beams, 1);

for b = 1:num_beams
    beam_field = beam_fields{b};
    beam = rtplan.BeamSequence.(beam_field);
    
    % Get beam name for logging
    if isfield(beam, 'BeamName')
        beam_name = beam.BeamName;
    else
        beam_name = sprintf('Beam_%d', b);
    end
    
    fprintf('\n  --- Processing %s (Beam %d/%d) ---\n', beam_name, b, num_beams);
    
    % Get control point field names
    if ~isfield(beam, 'ControlPointSequence')
        fprintf('    [WARN] No ControlPointSequence found. Skipping beam.\n');
        continue;
    end
    
    cp_fields = fieldnames(beam.ControlPointSequence);
    num_cp = length(cp_fields);
    fprintf('    Number of control points: %d\n', num_cp);
    
    %% Identify MLC devices in BeamLimitingDeviceSequence
    mlc_info = identifyMlcDevices(beam);
    
    if isempty(mlc_info)
        fprintf('    [WARN] No MLC devices found in beam. Skipping.\n');
        continue;
    end
    
    fprintf('    Found %d MLC device(s)\n', length(mlc_info));
    
    %% Process each MLC device (Halcyon has dual-layer: MLCX1 and MLCX2)
    for mlc_idx = 1:length(mlc_info)
        mlc_device = mlc_info(mlc_idx);
        fprintf('    Processing MLC: %s\n', mlc_device.type);
        
        %% Step 1: Identify dynamic vs static leaf pairs
        % Dynamic leaves change position across control points
        % We only correct gaps on dynamic leaves
        [dynamic_leaves, num_pairs] = identifyDynamicLeaves(...
            beam, cp_fields, mlc_device, POSITION_TOLERANCE);
        
        if num_pairs == 0
            fprintf('      [WARN] Could not determine leaf pair count. Skipping.\n');
            continue;
        end
        
        num_dynamic = sum(dynamic_leaves);
        fprintf('      Leaf pairs: %d, Dynamic: %d, Static: %d\n', ...
            num_pairs, num_dynamic, num_pairs - num_dynamic);
        
        %% Step 2: Correct small gaps for dynamic leaves
        mlc_corrections = 0;
        
        for cp_idx = 1:num_cp
            cp_field = cp_fields{cp_idx};
            
            % Get BeamLimitingDevicePositionSequence for this control point
            if ~isfield(beam.ControlPointSequence.(cp_field), 'BeamLimitingDevicePositionSequence')
                continue;
            end
            
            bld_pos_seq = beam.ControlPointSequence.(cp_field).BeamLimitingDevicePositionSequence;
            bld_pos_fields = fieldnames(bld_pos_seq);
            
            % Find the matching MLC device in this control point
            mlc_bld_idx = findMlcInControlPoint(bld_pos_seq, bld_pos_fields, mlc_device.type);
            
            if mlc_bld_idx == 0
                continue;
            end
            
            mlc_bld_field = bld_pos_fields{mlc_bld_idx};
            positions = bld_pos_seq.(mlc_bld_field).LeafJawPositions;
            
            % Ensure positions is a column vector
            positions = positions(:);
            
            % Verify we have correct number of positions
            if length(positions) ~= 2 * num_pairs
                fprintf('      [WARN] CP %d: Position count mismatch (expected %d, got %d)\n', ...
                    cp_idx, 2 * num_pairs, length(positions));
                continue;
            end
            
            % Split into Bank A (negative side) and Bank B (positive side)
            % Convention: first half is Bank A, second half is Bank B
            bank_a = positions(1:num_pairs);
            bank_b = positions(num_pairs+1:end);
            
            % Correct gaps where needed
            [bank_a_corrected, bank_b_corrected, cp_corrections] = correctMlcGaps(...
                bank_a, bank_b, dynamic_leaves, ...
                MIN_GAP_THRESHOLD, EXPANSION_AMOUNT, LEAF_POS_MIN, LEAF_POS_MAX);
            
            if cp_corrections > 0
                % Update positions in the RTPLAN structure
                corrected_positions = [bank_a_corrected; bank_b_corrected];
                rtplan.BeamSequence.(beam_field).ControlPointSequence.(cp_field)...
                    .BeamLimitingDevicePositionSequence.(mlc_bld_field).LeafJawPositions = corrected_positions;
                
                mlc_corrections = mlc_corrections + cp_corrections;
            end
        end
        
        fprintf('      Corrections made for %s: %d\n', mlc_device.type, mlc_corrections);
        beam_corrections(b) = beam_corrections(b) + mlc_corrections;
    end
    
    total_corrections = total_corrections + beam_corrections(b);
    fprintf('    Total corrections for %s: %d\n', beam_name, beam_corrections(b));
end

num_corrections = total_corrections;

%% ======================== UPDATE DICOM METADATA ========================

fprintf('\n  Updating DICOM metadata...\n');

% Generate new SOP Instance UID for modified plan
% This is required because we're creating a new/modified object
rtplan.SOPInstanceUID = dicomuid;
rtplan.MediaStorageSOPInstanceUID = rtplan.SOPInstanceUID;

% Update plan label to indicate modification
if isfield(rtplan, 'RTPlanLabel')
    original_label = rtplan.RTPlanLabel;
    new_label = [original_label '_adj'];
    % DICOM RTPlanLabel has max length of 16 characters
    if length(new_label) > 16
        new_label = new_label(1:16);
    end
    rtplan.RTPlanLabel = new_label;
    fprintf('    RTPlanLabel: %s -> %s\n', original_label, new_label);
else
    rtplan.RTPlanLabel = 'adj_plan';
end

% Update plan description if present
if isfield(rtplan, 'RTPlanDescription')
    rtplan.RTPlanDescription = [rtplan.RTPlanDescription ' - MLC gaps adjusted'];
else
    rtplan.RTPlanDescription = 'MLC gaps adjusted';
end

%% ======================== WRITE MODIFIED DICOM ========================

fprintf('  Writing adjusted RTPLAN...\n');
fprintf('    Output: %s\n', adjusted_filename);

try
    % dicomwrite requires empty first argument when writing from metadata struct
    dicomwrite([], adjusted_rtplan_path, rtplan, 'CreateMode', 'Copy');
    fprintf('    Write successful.\n');
catch ME
    error('step05_fix_mlc_gaps:DicomWriteError', ...
        'Failed to write adjusted RTPLAN: %s\nError: %s', adjusted_rtplan_path, ME.message);
end

%% ======================== SUMMARY ========================

fprintf('\n  ========================================\n');
fprintf('  Step 0.5 Complete\n');
fprintf('  ========================================\n');
fprintf('  Patient: %s\n', patient_id);
fprintf('  Session: %s\n', session);
fprintf('  Total beams processed: %d\n', num_beams);
fprintf('  Total MLC gap corrections: %d\n', num_corrections);
fprintf('  Output file: %s\n', adjusted_rtplan_path);

if num_corrections == 0
    fprintf('  [INFO] No corrections were needed.\n');
end

fprintf('  ========================================\n\n');

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function rtplan_file = findRtplanFile(directory)
%FINDRTPLANFILE Find RTPLAN DICOM file in directory
%
%   rtplan_file = findRtplanFile(directory)
%
%   Search for RTPLAN file by checking DICOM Modality field.
%   Returns filename (not full path) or empty string if not found.
%   Skips files that already have '_adjusted' in the name.

    rtplan_file = '';
    
    % Get all files in directory
    files = dir(fullfile(directory, '*'));
    
    for i = 1:length(files)
        % Skip directories
        if files(i).isdir
            continue;
        end
        
        % Skip files that are already adjusted versions
        if contains(files(i).name, '_adjusted', 'IgnoreCase', true)
            continue;
        end
        
        filepath = fullfile(directory, files(i).name);
        
        try
            % Try to read DICOM info
            info = dicominfo(filepath);
            
            % Check if it's an RTPLAN
            % Method 1: Check Modality field
            if isfield(info, 'Modality') && strcmpi(info.Modality, 'RTPLAN')
                rtplan_file = files(i).name;
                return;
            end
            
            % Method 2: Check SOPClassUID (RTPLAN = 1.2.840.10008.5.1.4.1.1.481.5)
            if isfield(info, 'SOPClassUID') && ...
               contains(info.SOPClassUID, '1.2.840.10008.5.1.4.1.1.481.5')
                rtplan_file = files(i).name;
                return;
            end
        catch
            % Not a valid DICOM file or unreadable, skip
            continue;
        end
    end
end


function mlc_info = identifyMlcDevices(beam)
%IDENTIFYMLCDEVICES Identify MLC devices in BeamLimitingDeviceSequence
%
%   mlc_info = identifyMlcDevices(beam)
%
%   Returns structure array with MLC device information:
%       .type           - Device type string (e.g., 'MLCX', 'MLCX1', 'MLCX2')
%       .bld_seq_idx    - Index in BeamLimitingDeviceSequence
%       .bld_seq_field  - Field name in BeamLimitingDeviceSequence
%       .num_pairs      - Number of leaf pairs
%       .boundaries     - Leaf position boundaries (if available)

    mlc_info = [];
    
    if ~isfield(beam, 'BeamLimitingDeviceSequence')
        return;
    end
    
    bld_seq = beam.BeamLimitingDeviceSequence;
    bld_fields = fieldnames(bld_seq);
    
    mlc_count = 0;
    
    for i = 1:length(bld_fields)
        device = bld_seq.(bld_fields{i});
        
        if ~isfield(device, 'RTBeamLimitingDeviceType')
            continue;
        end
        
        device_type = device.RTBeamLimitingDeviceType;
        
        % Check if this is an MLC device
        % Common types: MLCX, MLCY, MLCX1, MLCX2 (Halcyon dual-layer)
        if contains(device_type, 'MLC', 'IgnoreCase', true)
            mlc_count = mlc_count + 1;
            mlc_info(mlc_count).type = device_type; %#ok<AGROW>
            mlc_info(mlc_count).bld_seq_idx = i; %#ok<AGROW>
            mlc_info(mlc_count).bld_seq_field = bld_fields{i}; %#ok<AGROW>
            
            % Get number of leaf pairs
            if isfield(device, 'NumberOfLeafJawPairs')
                mlc_info(mlc_count).num_pairs = device.NumberOfLeafJawPairs; %#ok<AGROW>
            else
                mlc_info(mlc_count).num_pairs = 0; %#ok<AGROW> % Will be determined later
            end
            
            % Get leaf boundaries if available
            if isfield(device, 'LeafPositionBoundaries')
                mlc_info(mlc_count).boundaries = device.LeafPositionBoundaries; %#ok<AGROW>
            else
                mlc_info(mlc_count).boundaries = []; %#ok<AGROW>
            end
        end
    end
end


function mlc_bld_idx = findMlcInControlPoint(bld_pos_seq, bld_pos_fields, target_type)
%FINDMLCINCONTROLPOINT Find index of specific MLC type in control point
%
%   mlc_bld_idx = findMlcInControlPoint(bld_pos_seq, bld_pos_fields, target_type)
%
%   Returns the index (1-based) of the MLC device matching target_type
%   in the BeamLimitingDevicePositionSequence, or 0 if not found.

    mlc_bld_idx = 0;
    
    for i = 1:length(bld_pos_fields)
        device = bld_pos_seq.(bld_pos_fields{i});
        
        if isfield(device, 'RTBeamLimitingDeviceType')
            if strcmpi(device.RTBeamLimitingDeviceType, target_type)
                mlc_bld_idx = i;
                return;
            end
        end
    end
end


function [dynamic_leaves, num_pairs] = identifyDynamicLeaves(beam, cp_fields, mlc_device, tolerance)
%IDENTIFYDYNAMICLEAVES Identify which leaf pairs are dynamic vs static
%
%   [dynamic_leaves, num_pairs] = identifyDynamicLeaves(beam, cp_fields, mlc_device, tolerance)
%
%   A leaf pair is considered "dynamic" if either leaf (Bank A or Bank B)
%   changes position across control points by more than the tolerance.
%
%   Returns:
%       dynamic_leaves - Logical array (num_pairs x 1), true if dynamic
%       num_pairs      - Number of leaf pairs

    num_cp = length(cp_fields);
    num_pairs = mlc_device.num_pairs;
    
    % Collect all leaf positions across control points
    leaf_positions_all = cell(num_cp, 1);
    
    for cp_idx = 1:num_cp
        cp = beam.ControlPointSequence.(cp_fields{cp_idx});
        
        if ~isfield(cp, 'BeamLimitingDevicePositionSequence')
            continue;
        end
        
        bld_pos_seq = cp.BeamLimitingDevicePositionSequence;
        bld_pos_fields = fieldnames(bld_pos_seq);
        mlc_bld_idx = findMlcInControlPoint(bld_pos_seq, bld_pos_fields, mlc_device.type);
        
        if mlc_bld_idx > 0
            positions = bld_pos_seq.(bld_pos_fields{mlc_bld_idx}).LeafJawPositions;
            leaf_positions_all{cp_idx} = positions(:);
        end
    end
    
    % Remove empty entries
    valid_positions = leaf_positions_all(~cellfun(@isempty, leaf_positions_all));
    
    if isempty(valid_positions)
        dynamic_leaves = false(num_pairs, 1);
        num_pairs = 0;
        return;
    end
    
    % Determine num_pairs from data if not set
    first_positions = valid_positions{1};
    total_leaves = length(first_positions);
    actual_num_pairs = total_leaves / 2;
    
    if num_pairs == 0
        num_pairs = actual_num_pairs;
    end
    
    % Stack all positions into matrix (num_leaves x num_valid_cp)
    positions_matrix = zeros(total_leaves, length(valid_positions));
    for i = 1:length(valid_positions)
        positions_matrix(:, i) = valid_positions{i};
    end
    
    % Check each leaf position for changes across control points
    leaf_is_dynamic = false(total_leaves, 1);
    
    for leaf = 1:total_leaves
        leaf_pos = positions_matrix(leaf, :);
        pos_range = max(leaf_pos) - min(leaf_pos);
        
        if pos_range > tolerance
            leaf_is_dynamic(leaf) = true;
        end
    end
    
    % A leaf PAIR is dynamic if either Bank A or Bank B leaf is dynamic
    bank_a_dynamic = leaf_is_dynamic(1:num_pairs);
    bank_b_dynamic = leaf_is_dynamic(num_pairs+1:end);
    
    dynamic_leaves = bank_a_dynamic | bank_b_dynamic;
end


function [bank_a_out, bank_b_out, num_corrected] = correctMlcGaps(...
    bank_a, bank_b, dynamic_leaves, min_gap, expansion, pos_min, pos_max)
%CORRECTMLCGAPS Correct small MLC gaps for dynamic leaf pairs
%
%   [bank_a_out, bank_b_out, num_corrected] = correctMlcGaps(...)
%
%   For dynamic leaf pairs with gaps smaller than min_gap, expand the gap
%   symmetrically by moving Bank A more negative and Bank B more positive.
%
%   Inputs:
%       bank_a          - Positions of Bank A (negative/left side)
%       bank_b          - Positions of Bank B (positive/right side)
%       dynamic_leaves  - Logical array indicating which pairs are dynamic
%       min_gap         - Minimum gap threshold (gaps smaller than this are corrected)
%       expansion       - Amount to expand on each side
%       pos_min, pos_max - Allowed position range
%
%   Outputs:
%       bank_a_out, bank_b_out - Corrected positions
%       num_corrected          - Number of corrections made

    bank_a_out = bank_a;
    bank_b_out = bank_b;
    num_corrected = 0;
    
    num_pairs = length(bank_a);
    
    for i = 1:num_pairs
        % Skip static leaves - only correct dynamic ones
        if ~dynamic_leaves(i)
            continue;
        end
        
        % Calculate gap
        % Convention: Bank B (positive side) - Bank A (negative side)
        % Positive gap = open, negative gap = overlapping (error state)
        gap = bank_b(i) - bank_a(i);
        
        % Check if gap is below threshold and positive (closed or nearly closed)
        % We only fix small positive gaps, not negative (overlapping) ones
        if gap >= 0 && gap < min_gap
            % Need to expand the gap symmetrically
            % Try expanding both sides by 'expansion' amount
            new_bank_a = bank_a(i) - expansion;  % Move Bank A further negative
            new_bank_b = bank_b(i) + expansion;  % Move Bank B further positive
            
            % Check boundary conditions
            a_within_bounds = (new_bank_a >= pos_min);
            b_within_bounds = (new_bank_b <= pos_max);
            
            if a_within_bounds && b_within_bounds
                % Both sides can be expanded normally
                bank_a_out(i) = new_bank_a;
                bank_b_out(i) = new_bank_b;
                
            elseif ~a_within_bounds && b_within_bounds
                % Bank A hits boundary, compensate by expanding Bank B more
                bank_a_out(i) = pos_min;
                bank_b_out(i) = bank_b(i) + 2 * expansion;
                % Clamp Bank B if it exceeds max
                if bank_b_out(i) > pos_max
                    bank_b_out(i) = pos_max;
                end
                
            elseif a_within_bounds && ~b_within_bounds
                % Bank B hits boundary, compensate by expanding Bank A more
                bank_a_out(i) = bank_a(i) - 2 * expansion;
                bank_b_out(i) = pos_max;
                % Clamp Bank A if it exceeds min
                if bank_a_out(i) < pos_min
                    bank_a_out(i) = pos_min;
                end
                
            else
                % Both at boundaries - maximize gap as much as possible
                bank_a_out(i) = pos_min;
                bank_b_out(i) = pos_max;
            end
            
            num_corrected = num_corrected + 1;
        end
    end
end
