%% fix_mlc_gaps.m
% Script to correct small MLC gaps in RTPLAN files from Ethos exports
% For Halcyon LINAC with dual-layer MLC
%
% This script:
% 1. Identifies dynamic vs static MLC leaf pairs per beam
% 2. Expands gaps that are < 0.5mm for dynamic leaves by 0.4mm on each side
% 3. Handles boundary conditions (leaf positions must stay within [-140, 140])
% 4. Exports corrected plan with "adjusted_mlc" suffix
%
% Author: Generated for ETHOS Simulations project
% Date: January 2025

clear; clc;

%% Configuration
% Working directory base path
wd = '/mnt/weka/home/80030361/ETHOS_Simulations';

% Patient and session lists (expandable)
patient_ids = {'1194203'};
sessions = {'Session_1'};

% MLC gap correction parameters
MIN_GAP_THRESHOLD = 0.5;    % mm - gaps smaller than this will be corrected
EXPANSION_AMOUNT = 0.4;     % mm - amount to expand on each side
LEAF_POS_MIN = -140;        % mm - minimum allowed leaf position
LEAF_POS_MAX = 140;         % mm - maximum allowed leaf position
POSITION_TOLERANCE = 1e-6;  % mm - tolerance for determining if leaves moved

%% Main Processing Loop
for pid = 1:length(patient_ids)
    for sid = 1:length(sessions)
        
        patient_id = patient_ids{pid};
        session = sessions{sid};
        
        fprintf('\n========================================\n');
        fprintf('Processing Patient: %s, Session: %s\n', patient_id, session);
        fprintf('========================================\n');
        
        % Construct path to raw data
        rawwd = fullfile(wd, 'EthosExports', patient_id, 'Pancreas', session, 'sct');
        
        % Check if directory exists
        if ~isfolder(rawwd)
            warning('Directory not found: %s\nSkipping...', rawwd);
            continue;
        end
        
        %% Find RTPLAN file in directory
        rtplan_file = find_rtplan_file(rawwd);
        
        if isempty(rtplan_file)
            warning('No RTPLAN file found in: %s\nSkipping...', rawwd);
            continue;
        end
        
        fprintf('Found RTPLAN file: %s\n', rtplan_file);
        
        %% Read RTPLAN DICOM
        fprintf('Reading RTPLAN DICOM...\n');
        rtplan = dicominfo(fullfile(rawwd, rtplan_file));
        
        %% Get beam sequence field names
        beam_fields = fieldnames(rtplan.BeamSequence);
        num_beams = length(beam_fields);
        fprintf('Number of beams: %d\n', num_beams);
        
        % Initialize statistics
        total_corrections = 0;
        
        %% Process each beam
        for b = 1:num_beams
            beam_field = beam_fields{b};
            beam = rtplan.BeamSequence.(beam_field);
            
            fprintf('\n--- Beam %d: %s ---\n', b, beam.BeamName);
            
            % Get control point field names
            cp_fields = fieldnames(beam.ControlPointSequence);
            num_cp = length(cp_fields);
            fprintf('Number of control points: %d\n', num_cp);
            
            %% Identify MLC devices in BeamLimitingDeviceSequence
            mlc_info = identify_mlc_devices(beam);
            
            if isempty(mlc_info)
                fprintf('No MLC devices found in beam. Skipping...\n');
                continue;
            end
            
            %% Process each MLC device (Halcyon has dual-layer)
            for mlc_idx = 1:length(mlc_info)
                mlc_device = mlc_info(mlc_idx);
                fprintf('\nProcessing MLC: %s (Index: %d)\n', mlc_device.type, mlc_device.bld_seq_idx);
                
                %% Step 1: Identify dynamic vs static leaf pairs
                [dynamic_leaves, leaf_positions_all] = identify_dynamic_leaves(...
                    beam, cp_fields, mlc_device, POSITION_TOLERANCE);
                
                num_leaf_pairs = mlc_device.num_pairs;
                num_dynamic = sum(dynamic_leaves);
                fprintf('Leaf pairs: %d, Dynamic: %d, Static: %d\n', ...
                    num_leaf_pairs, num_dynamic, num_leaf_pairs - num_dynamic);
                
                %% Step 2: Correct small gaps for dynamic leaves
                beam_corrections = 0;
                
                for cp_idx = 1:num_cp
                    cp_field = cp_fields{cp_idx};
                    
                    % Get BeamLimitingDevicePositionSequence for this control point
                    if ~isfield(beam.ControlPointSequence.(cp_field), 'BeamLimitingDevicePositionSequence')
                        continue;
                    end
                    
                    bld_pos_seq = beam.ControlPointSequence.(cp_field).BeamLimitingDevicePositionSequence;
                    bld_pos_fields = fieldnames(bld_pos_seq);
                    
                    % Find the matching MLC device in this control point
                    mlc_bld_idx = find_mlc_in_control_point(bld_pos_seq, bld_pos_fields, mlc_device.type);
                    
                    if mlc_bld_idx == 0
                        continue;
                    end
                    
                    mlc_bld_field = bld_pos_fields{mlc_bld_idx};
                    positions = bld_pos_seq.(mlc_bld_field).LeafJawPositions;
                    
                    % Ensure positions is a column vector
                    positions = positions(:);
                    
                    % Split into Bank A (negative side) and Bank B (positive side)
                    % Convention: first half is Bank A, second half is Bank B
                    bank_a = positions(1:num_leaf_pairs);
                    bank_b = positions(num_leaf_pairs+1:end);
                    
                    % Calculate gaps and correct where needed
                    [bank_a_corrected, bank_b_corrected, num_corrected] = correct_mlc_gaps(...
                        bank_a, bank_b, dynamic_leaves, ...
                        MIN_GAP_THRESHOLD, EXPANSION_AMOUNT, LEAF_POS_MIN, LEAF_POS_MAX);
                    
                    if num_corrected > 0
                        % Update positions in the RTPLAN
                        corrected_positions = [bank_a_corrected; bank_b_corrected];
                        rtplan.BeamSequence.(beam_field).ControlPointSequence.(cp_field)...
                            .BeamLimitingDevicePositionSequence.(mlc_bld_field).LeafJawPositions = corrected_positions;
                        
                        beam_corrections = beam_corrections + num_corrected;
                    end
                end
                
                fprintf('Corrections made for %s: %d\n', mlc_device.type, beam_corrections);
                total_corrections = total_corrections + beam_corrections;
            end
        end
        
        fprintf('\n========================================\n');
        fprintf('Total corrections for patient %s: %d\n', patient_id, total_corrections);
        
        %% Generate new SOP Instance UID for modified plan
        rtplan.SOPInstanceUID = dicomuid;
        rtplan.MediaStorageSOPInstanceUID = rtplan.SOPInstanceUID;
        
        % Update plan label and description to indicate modification
        original_label = rtplan.RTPlanLabel;
        rtplan.RTPlanLabel = [original_label '_adj'];
        if length(rtplan.RTPlanLabel) > 16
            rtplan.RTPlanLabel = rtplan.RTPlanLabel(1:16);
        end
        
        if isfield(rtplan, 'RTPlanDescription')
            rtplan.RTPlanDescription = [rtplan.RTPlanDescription ' - MLC gaps adjusted'];
        end
        
        %% Export corrected RTPLAN
        [~, original_name, ~] = fileparts(rtplan_file);
        output_filename = [original_name '_adjusted_mlc.dcm'];
        output_path = fullfile(rawwd, output_filename);
        
        fprintf('\nExporting corrected plan to:\n%s\n', output_path);
        
        % Write DICOM file
        dicomwrite([], output_path, rtplan, 'CreateMode', 'Copy');
        
        fprintf('Export complete!\n');
    end
end

fprintf('\n========================================\n');
fprintf('All processing complete!\n');
fprintf('========================================\n');

%% ========================================================================
%  Helper Functions
%% ========================================================================

function rtplan_file = find_rtplan_file(directory)
    % Find RTPLAN DICOM file in directory by checking DICOM modality
    % Returns filename or empty string if not found
    
    rtplan_file = '';
    
    % Get all files in directory
    files = dir(fullfile(directory, '*'));
    
    for i = 1:length(files)
        if files(i).isdir
            continue;
        end
        
        filepath = fullfile(directory, files(i).name);
        
        try
            % Try to read DICOM info
            info = dicominfo(filepath);
            
            % Check if it's an RTPLAN (Modality = 'RTPLAN' or check SOPClassUID)
            % RTPLAN SOP Class UID: 1.2.840.10008.5.1.4.1.1.481.5
            if isfield(info, 'Modality') && strcmpi(info.Modality, 'RTPLAN')
                rtplan_file = files(i).name;
                return;
            end
            
            if isfield(info, 'SOPClassUID') && contains(info.SOPClassUID, '1.2.840.10008.5.1.4.1.1.481.5')
                rtplan_file = files(i).name;
                return;
            end
        catch
            % Not a valid DICOM file, skip
            continue;
        end
    end
end

function mlc_info = identify_mlc_devices(beam)
    % Identify MLC devices in the BeamLimitingDeviceSequence
    % Returns structure array with MLC device information
    
    mlc_info = [];
    
    if ~isfield(beam, 'BeamLimitingDeviceSequence')
        return;
    end
    
    bld_seq = beam.BeamLimitingDeviceSequence;
    bld_fields = fieldnames(bld_seq);
    
    mlc_count = 0;
    
    for i = 1:length(bld_fields)
        device = bld_seq.(bld_fields{i});
        device_type = device.RTBeamLimitingDeviceType;
        
        % Check if this is an MLC (MLCX, MLCY, MLCX1, MLCX2, etc.)
        if contains(device_type, 'MLC', 'IgnoreCase', true)
            mlc_count = mlc_count + 1;
            mlc_info(mlc_count).type = device_type;
            mlc_info(mlc_count).bld_seq_idx = i;
            mlc_info(mlc_count).bld_seq_field = bld_fields{i};
            
            % Get number of leaf pairs
            if isfield(device, 'NumberOfLeafJawPairs')
                mlc_info(mlc_count).num_pairs = device.NumberOfLeafJawPairs;
            else
                % Try to infer from leaf positions in first control point
                mlc_info(mlc_count).num_pairs = 0; % Will be set later
            end
            
            % Get leaf boundaries if available
            if isfield(device, 'LeafPositionBoundaries')
                mlc_info(mlc_count).boundaries = device.LeafPositionBoundaries;
            else
                mlc_info(mlc_count).boundaries = [];
            end
        end
    end
end

function mlc_bld_idx = find_mlc_in_control_point(bld_pos_seq, bld_pos_fields, target_type)
    % Find the index of specific MLC type in BeamLimitingDevicePositionSequence
    
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

function [dynamic_leaves, leaf_positions_all] = identify_dynamic_leaves(beam, cp_fields, mlc_device, tolerance)
    % Identify which leaf pairs are dynamic (change position) vs static
    %
    % Returns:
    %   dynamic_leaves: logical array (num_pairs x 1), true if leaf pair is dynamic
    %   leaf_positions_all: cell array of all leaf positions for each control point
    
    num_cp = length(cp_fields);
    num_pairs = mlc_device.num_pairs;
    
    % If num_pairs not set, determine from first control point
    if num_pairs == 0
        for cp_idx = 1:num_cp
            cp = beam.ControlPointSequence.(cp_fields{cp_idx});
            if isfield(cp, 'BeamLimitingDevicePositionSequence')
                bld_pos_seq = cp.BeamLimitingDevicePositionSequence;
                bld_pos_fields = fieldnames(bld_pos_seq);
                mlc_bld_idx = find_mlc_in_control_point(bld_pos_seq, bld_pos_fields, mlc_device.type);
                if mlc_bld_idx > 0
                    positions = bld_pos_seq.(bld_pos_fields{mlc_bld_idx}).LeafJawPositions;
                    num_pairs = length(positions) / 2;
                    break;
                end
            end
        end
    end
    
    % Initialize storage
    leaf_positions_all = cell(num_cp, 1);
    
    % Collect all leaf positions across control points
    for cp_idx = 1:num_cp
        cp = beam.ControlPointSequence.(cp_fields{cp_idx});
        
        if ~isfield(cp, 'BeamLimitingDevicePositionSequence')
            continue;
        end
        
        bld_pos_seq = cp.BeamLimitingDevicePositionSequence;
        bld_pos_fields = fieldnames(bld_pos_seq);
        mlc_bld_idx = find_mlc_in_control_point(bld_pos_seq, bld_pos_fields, mlc_device.type);
        
        if mlc_bld_idx > 0
            positions = bld_pos_seq.(bld_pos_fields{mlc_bld_idx}).LeafJawPositions;
            leaf_positions_all{cp_idx} = positions(:);
        end
    end
    
    % Remove empty entries
    valid_positions = leaf_positions_all(~cellfun(@isempty, leaf_positions_all));
    
    if isempty(valid_positions)
        dynamic_leaves = false(num_pairs, 1);
        return;
    end
    
    % Stack all positions into matrix (num_leaves x num_valid_cp)
    positions_matrix = cell2mat(cellfun(@(x) x', valid_positions, 'UniformOutput', false))';
    
    % Check each leaf position for changes across control points
    % Total leaves = 2 * num_pairs (Bank A + Bank B)
    total_leaves = size(positions_matrix, 1);
    actual_num_pairs = total_leaves / 2;
    
    % Update num_pairs if we determined it from data
    if num_pairs == 0
        num_pairs = actual_num_pairs;
    end
    
    % Initialize dynamic tracking for each leaf
    leaf_is_dynamic = false(total_leaves, 1);
    
    % Check each leaf for position changes
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

function [bank_a_out, bank_b_out, num_corrected] = correct_mlc_gaps(...
    bank_a, bank_b, dynamic_leaves, min_gap, expansion, pos_min, pos_max)
    % Correct small MLC gaps for dynamic leaf pairs
    %
    % Inputs:
    %   bank_a: positions of Bank A (negative/left side)
    %   bank_b: positions of Bank B (positive/right side)
    %   dynamic_leaves: logical array indicating which pairs are dynamic
    %   min_gap: minimum gap threshold (gaps smaller than this are corrected)
    %   expansion: amount to expand on each side
    %   pos_min, pos_max: allowed position range
    %
    % Outputs:
    %   bank_a_out, bank_b_out: corrected positions
    %   num_corrected: number of corrections made
    
    bank_a_out = bank_a;
    bank_b_out = bank_b;
    num_corrected = 0;
    
    num_pairs = length(bank_a);
    
    for i = 1:num_pairs
        % Skip static leaves
        if ~dynamic_leaves(i)
            continue;
        end
        
        % Calculate gap (Bank B - Bank A should be positive for open leaves)
        % Convention: Bank A is negative side (moves in -X), Bank B is positive side (moves in +X)
        gap = bank_b(i) - bank_a(i);
        
        % Check if gap is below threshold (and positive - closed or nearly closed)
        if gap >= 0 && gap < min_gap
            % Need to expand the gap
            % Try expanding both sides by 'expansion' amount
            new_bank_a = bank_a(i) - expansion;  % Move Bank A further negative
            new_bank_b = bank_b(i) + expansion;  % Move Bank B further positive
            
            % Check boundary conditions
            a_ok = (new_bank_a >= pos_min);
            b_ok = (new_bank_b <= pos_max);
            
            if a_ok && b_ok
                % Both sides can be expanded
                bank_a_out(i) = new_bank_a;
                bank_b_out(i) = new_bank_b;
            elseif ~a_ok && b_ok
                % Bank A hits boundary, expand Bank B by full amount
                bank_a_out(i) = pos_min;
                bank_b_out(i) = bank_b(i) + 2 * expansion;
                % Check if Bank B now exceeds max
                if bank_b_out(i) > pos_max
                    bank_b_out(i) = pos_max;
                end
            elseif a_ok && ~b_ok
                % Bank B hits boundary, expand Bank A by full amount
                bank_a_out(i) = bank_a(i) - 2 * expansion;
                bank_b_out(i) = pos_max;
                % Check if Bank A now exceeds min
                if bank_a_out(i) < pos_min
                    bank_a_out(i) = pos_min;
                end
            else
                % Both at boundaries - can't expand further
                bank_a_out(i) = pos_min;
                bank_b_out(i) = pos_max;
            end
            
            num_corrected = num_corrected + 1;
        end
    end
end