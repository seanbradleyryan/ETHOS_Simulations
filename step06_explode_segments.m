function output_paths = step06_explode_segments(patient_id, session, config)
%STEP06_EXPLODE_SEGMENTS  Per-beam plan explosion for both REFERENCE and ADAPTED plans.
%
%  Takes the MLC-gap-corrected RTPLAN files produced by step05_fix_mlc_gaps
%  (RP_reference_adjusted_mlc.dcm and RP_adapted_adjusted_mlc.dcm) and
%  produces ONE output RTPLAN file PER ORIGINAL BEAM for each plan type.
%  Within each output file every segment of that beam is promoted to its
%  own independent two-control-point beam, so:
%
%    original plan with 17 beams, 165 segments each
%      -> 17 output RTPLAN files
%         each file has 165 beams (one per segment), one FractionGroup
%
%  OUTPUT FILENAME FORMAT:
%    RTPLAN_{patient_id}_{session}_{reference|adapted}_B{orig_beam}.dcm
%    e.g.  RTPLAN_1194203_Session_1_reference_B10.dcm
%          RTPLAN_1194203_Session_1_adapted_B10.dcm
%
%  INTRA-DICOM NAMING:
%    RTPlanLabel (max 16 chars): {session} {plan_type} B{orig_beam_number}  (truncated if needed)
%    RTPlanName:                 {session} {plan_type} B{orig_beam_number}
%    RTPlanDescription:          {session} {plan_type} B{orig} ({beam_name}) exploded segments
%
%  Each output file contains:
%    - BeamSequence with N_seg beams (one per segment, BeamNumber 1..N_seg)
%    - FractionGroupSequence with N_seg referenced beams and per-segment MUs
%
%  RayStation imports these plans; it exports one dose file per beam, which
%  step15 loads via the preferred pattern:  Plan_Field*_Beam*_B*_S*.dcm
%
%  INPUTS:
%    patient_id  - char,   e.g. '1194203'
%    session     - char,   e.g. 'Session_1'
%    config      - struct with fields:
%                    .working_dir      (char)  root workspace directory
%                    .treatment_site   (char)  e.g. 'Pancreas'
%
%  OUTPUT:
%    output_paths - struct with fields:
%                    .reference  - cell array of char, full paths for reference plan files
%                    .adapted    - cell array of char, full paths for adapted plan files
%                    .all        - cell array of char, all written paths combined
%
%  STANDALONE USAGE (no arguments  uses built-in defaults):
%    step06_explode_segments();
%
%  PIPELINE USAGE (called by pipeline_setup.m):
%    output_paths = step06_explode_segments(patient_id, session, CONFIG);

% -----------------------------------------------------------------------
% Standalone defaults  (used only when called with no arguments)
% -----------------------------------------------------------------------
if nargin == 0
    patient_id = '1194203';
    session    = 'Session_1';
    config     = struct( ...
        'working_dir',    '/mnt/weka/home/80030361/ETHOS_Simulations', ...
        'treatment_site', 'Pancreas');
end

fprintf('=== Step 0.6: Segment Explosion  1 Plan per Beam (Reference + Adapted) ===\n');
fprintf('Patient: %s | Session: %s\n\n', patient_id, session);

% -----------------------------------------------------------------------
% Paths
% -----------------------------------------------------------------------
sct_dir    = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session, 'sct');

output_dir = fullfile(config.working_dir, 'Raystation_Input', patient_id, session);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('SCT dir:    %s\n', sct_dir);
fprintf('Output dir: %s\n\n', output_dir);

% -----------------------------------------------------------------------
% Copy supporting files from sct_dir and sim_ct_dir into output_dir
% -----------------------------------------------------------------------
fprintf('--- Copying supporting files to RayStation input directory ---\n');

% sCT image slices
copy_files_to_dir(sct_dir, 'CT*.dcm', output_dir, 'sCT images');

% Reference RTSTRUCT only
copy_files_to_dir(sct_dir, 'RS_reference.dcm', output_dir, 'RTSTRUCT (reference)');

% Image registration (REG) files
copy_files_to_dir(sct_dir, 'REG_*.dcm', output_dir, 'Image registrations');

% CBCT series (up to two, sorted by sort_CBCT)
copy_files_to_dir(sct_dir, 'CBCT*.dcm', output_dir, 'CBCT series');

fprintf('\n');

% -----------------------------------------------------------------------
% Plan types to process
% -----------------------------------------------------------------------
plan_types  = {'reference'};
input_files = { ...
    fullfile(sct_dir, 'RP_reference_adjusted_mlc.dcm')};

output_paths.reference = {};

% -----------------------------------------------------------------------
% Main loop  one pass per plan type
% -----------------------------------------------------------------------
for pt = 1:numel(plan_types)

    plan_type    = plan_types{pt};
    input_rtplan = input_files{pt};

    fprintf('--- Processing %s plan ---\n', upper(plan_type));
    fprintf('Input: %s\n', input_rtplan);

    if ~isfile(input_rtplan)
        warning('step06:fileNotFound', ...
            'Input file not found for %s plan, skipping:\n  %s', plan_type, input_rtplan);
        fprintf('[SKIP] %s plan  file missing.\n\n', upper(plan_type));
        continue;
    end

    % -------------------------------------------------------------------
    % Load and validate RTPLAN
    % -------------------------------------------------------------------
    fprintf('Loading RTPLAN...\n');
    rtplan = dicominfo(input_rtplan, 'UseDictionaryVR', true);

    if ~isfield(rtplan, 'Modality') || ~strcmp(rtplan.Modality, 'RTPLAN')
        error('step06:wrongModality', ...
            '[%s] File Modality is not RTPLAN: %s', plan_type, rtplan.Modality);
    end
    if ~isfield(rtplan, 'BeamSequence')
        error('step06:noBeams', '[%s] RTPLAN has no BeamSequence.', plan_type);
    end

    beam_fields_orig   = fieldnames(rtplan.BeamSequence);
    num_original_beams = numel(beam_fields_orig);
    fprintf('RTPlanLabel:         %s\n', rtplan.RTPlanLabel);
    fprintf('Original beam count: %d\n\n', num_original_beams);

    % -------------------------------------------------------------------
    % Build MU map: beam_number -> total_MU
    % -------------------------------------------------------------------
    mu_map = containers.Map('KeyType', 'int32', 'ValueType', 'double');

    if ~isfield(rtplan, 'FractionGroupSequence')
        error('step06:noFractionGroup', ...
            '[%s] RTPLAN has no FractionGroupSequence.', plan_type);
    end
    frac_group      = rtplan.FractionGroupSequence.Item_1;
    ref_beam_fields = fieldnames(frac_group.ReferencedBeamSequence);

    for k = 1:numel(ref_beam_fields)
        rb        = frac_group.ReferencedBeamSequence.(ref_beam_fields{k});
        bnum      = int32(rb.ReferencedBeamNumber);
        bmu       = double(rb.BeamMeterset);
        mu_map(bnum) = bmu;
    end
    fprintf('MU map loaded for %d beams.\n\n', mu_map.Count);

    % -------------------------------------------------------------------
    % Per-beam MU accumulators for validation
    % -------------------------------------------------------------------
    mu_sum_per_orig  = containers.Map('KeyType', 'int32', 'ValueType', 'double');
    mu_orig_per_beam = containers.Map('KeyType', 'int32', 'ValueType', 'double');

    plan_output_paths = {};
    any_warn_global   = false;

    % -------------------------------------------------------------------
    % Beam loop  each original beam produces ONE output RTPLAN file
    %             containing all its segments promoted to beams
    % -------------------------------------------------------------------
    for b = 1:num_original_beams

        orig_beam            = rtplan.BeamSequence.(beam_fields_orig{b});
        original_beam_number = int32(orig_beam.BeamNumber);
        original_beam_name   = orig_beam.BeamName;

        if ~isKey(mu_map, original_beam_number)
            error('step06:muNotFound', ...
                '[%s] No MU found in FractionGroup for beam %d (%s).', ...
                plan_type, original_beam_number, original_beam_name);
        end
        total_mu = mu_map(original_beam_number);
        mu_orig_per_beam(original_beam_number) = total_mu;
        mu_sum_per_orig(original_beam_number)  = 0;

        % Sort control point fields numerically
        cp_fields = fieldnames(orig_beam.ControlPointSequence);
        cp_nums   = cellfun(@(f) sscanf(f, 'Item_%d'), cp_fields);
        [~, sort_idx] = sort(cp_nums);
        cp_fields = cp_fields(sort_idx);
        N_cp  = numel(cp_fields);
        N_seg = N_cp - 1;

        if N_seg < 1
            warning('step06:noSegments', ...
                '[%s] Beam %d (%s) has only %d control points  skipping.', ...
                plan_type, original_beam_number, original_beam_name, N_cp);
            continue;
        end

        % First CP of original beam (geometry back-fill source)
        cp1_orig = orig_beam.ControlPointSequence.(cp_fields{1});

        if isfield(cp1_orig, 'GantryAngle')
            gantry_str = sprintf('%.1f deg', double(cp1_orig.GantryAngle));
        else
            gantry_str = 'N/A';
        end
        fprintf('  B%d (%s, gantry %s, %d segments, MU=%.4f)...\n', ...
            original_beam_number, original_beam_name, gantry_str, N_seg, total_mu);

        % -----------------------------------------------------------
        % Accumulators for this beam's output plan
        %   new_beam_seq:     BeamSequence of the output RTPLAN
        %   new_ref_beam_seq: FractionGroup.ReferencedBeamSequence
        %   local_beam_idx:   running beam counter within this plan (1..N_seg)
        % -----------------------------------------------------------
        new_beam_seq     = struct();
        new_ref_beam_seq = struct();
        local_beam_idx   = 0;

        % -----------------------------------------------------------
        % Segment loop  each segment becomes one beam in the output plan
        % -----------------------------------------------------------
        for s = 1:N_seg

            local_beam_idx = local_beam_idx + 1;

            % Deep-copy original beam
            new_beam                       = orig_beam;
            new_beam.BeamNumber            = local_beam_idx;   % sequential within this plan
            new_beam.BeamName              = sprintf('B%d_S%03d', original_beam_number, s - 1);
            new_beam.NumberOfControlPoints = 2;

            % Bounding control points
            cp_entry_orig = orig_beam.ControlPointSequence.(cp_fields{s});
            cp_exit_orig  = orig_beam.ControlPointSequence.(cp_fields{s + 1});

            % CumulativeMetersetWeight values
            if isfield(cp_entry_orig, 'CumulativeMetersetWeight')
                cmw_entry = double(cp_entry_orig.CumulativeMetersetWeight);
            else
                cmw_entry = (s - 1) / N_seg;
            end
            if isfield(cp_exit_orig, 'CumulativeMetersetWeight')
                cmw_exit = double(cp_exit_orig.CumulativeMetersetWeight);
            else
                cmw_exit = s / N_seg;
            end

            % Segment MU
            seg_mu = (cmw_exit - cmw_entry) * total_mu;
            mu_sum_per_orig(original_beam_number) = ...
                mu_sum_per_orig(original_beam_number) + seg_mu;

            % --- Entry CP (Item_1): all geometry fields present ---
            cp_entry_new = cp_entry_orig;
            cp_entry_new.CumulativeMetersetWeight = 0;
            cp_entry_new.ControlPointIndex        = 0;

            % Back-fill geometry fields absent from this intermediate CP
            cp1_varying_fields = { ...
                'BeamLimitingDevicePositionSequence', ...
                'CumulativeMetersetWeight', ...
                'ControlPointIndex'};

            cp1_all_fields = fieldnames(cp1_orig);
            for gf = 1:numel(cp1_all_fields)
                fld = cp1_all_fields{gf};
                if ~isfield(cp_entry_new, fld) && ~ismember(fld, cp1_varying_fields)
                    cp_entry_new.(fld) = cp1_orig.(fld);
                end
            end

            % BeamLimitingDevicePositionSequence for entry CP:
            %   Item_1: X jaw  (from original beam CP_1)
            %   Item_2: Y jaw  (from original beam CP_1)
            %   Item_3: MLC1   (from this segment's entry CP)
            %   Item_4: MLC2   (from this segment's entry CP)
            bldps_cp1   = cp1_orig.BeamLimitingDevicePositionSequence;
            bldps_entry = cp_entry_orig.BeamLimitingDevicePositionSequence;

            n_bldps_entry = numel(fieldnames(bldps_entry));
            if n_bldps_entry == 4
                bldps_entry_mlc1 = bldps_entry.Item_3;
                bldps_entry_mlc2 = bldps_entry.Item_4;
            else
                bldps_entry_mlc1 = bldps_entry.Item_1;
                bldps_entry_mlc2 = bldps_entry.Item_2;
            end

            bldps_new_entry        = struct();
            bldps_new_entry.Item_1 = bldps_cp1.Item_1;    % X jaw
            bldps_new_entry.Item_2 = bldps_cp1.Item_2;    % Y jaw
            bldps_new_entry.Item_3 = bldps_entry_mlc1;    % MLC1
            bldps_new_entry.Item_4 = bldps_entry_mlc2;    % MLC2
            cp_entry_new.BeamLimitingDevicePositionSequence = bldps_new_entry;

            % --- Exit CP (Item_2): MLC-only ---
            cp_exit_new = cp_exit_orig;
            cp_exit_new.CumulativeMetersetWeight = 1;
            cp_exit_new.ControlPointIndex        = 1;

            if ~isfield(cp_exit_new, 'BeamLimitingDevicePositionSequence')
                bldps_exit_fb        = struct();
                bldps_exit_fb.Item_1 = bldps_entry_mlc1;
                bldps_exit_fb.Item_2 = bldps_entry_mlc2;
                cp_exit_new.BeamLimitingDevicePositionSequence = bldps_exit_fb;
            end

            % 2-item ControlPointSequence
            new_cp_seq        = struct();
            new_cp_seq.Item_1 = cp_entry_new;
            new_cp_seq.Item_2 = cp_exit_new;
            new_beam.ControlPointSequence         = new_cp_seq;
            new_beam.FinalCumulativeMetersetWeight = 1;

            % -----------------------------------------------------------
            % Accumulate this segment-beam into the plan's BeamSequence
            % and FractionGroup.ReferencedBeamSequence
            % -----------------------------------------------------------
            new_beam_seq.(sprintf('Item_%d', local_beam_idx)) = new_beam;

            ref_entry                      = struct();
            ref_entry.ReferencedBeamNumber = local_beam_idx;
            ref_entry.BeamMeterset         = seg_mu;
            new_ref_beam_seq.(sprintf('Item_%d', local_beam_idx)) = ref_entry;

            % Progress every 50 new beams within this plan
            if mod(local_beam_idx, 50) == 0
                fprintf('    Created beam %d of ~%d (B%d S%d, MU=%.4f)\n', ...
                    local_beam_idx, N_seg, original_beam_number, s - 1, seg_mu);
            end

        end  % segment loop

        fprintf('    -> %d segment-beams assembled for B%d\n', local_beam_idx, original_beam_number);

        % ---------------------------------------------------------------
        % MU validation for this beam
        % ---------------------------------------------------------------
        orig_mu   = mu_orig_per_beam(original_beam_number);
        summed    = mu_sum_per_orig(original_beam_number);
        delta_abs = abs(summed - orig_mu);
        delta_pct = 0;
        if orig_mu > 0
            delta_pct = 100 * delta_abs / orig_mu;
        end
        if delta_pct > 0.5
            mu_status       = '[WARN >0.5%]';
            any_warn_global = true;
        else
            mu_status = '[OK]';
        end
        fprintf('    MU: orig=%.4f  sum=%.4f  delta=%.4f%%  %s\n', ...
            orig_mu, summed, delta_pct, mu_status);

        % ---------------------------------------------------------------
        % Assemble output RTPLAN for this beam
        %   - BeamSequence         : N_seg beams (one per segment)
        %   - FractionGroupSequence: N_seg referenced beams with per-segment MUs
        %   - NumberOfBeams        : N_seg
        % ---------------------------------------------------------------
        rtplan_out = rtplan;   % deep copy of original header

        rtplan_out.BeamSequence = new_beam_seq;
        rtplan_out.FractionGroupSequence.Item_1.ReferencedBeamSequence = new_ref_beam_seq;
        rtplan_out.FractionGroupSequence.Item_1.NumberOfBeams = local_beam_idx;

        % Unique SOP Instance UID per file
        new_uid                               = dicomuid();
        rtplan_out.SOPInstanceUID             = new_uid;
        rtplan_out.MediaStorageSOPInstanceUID = new_uid;

        % -----------------------------------------------------------
        % Naming  (plan describes the original beam and session)
        %
        %   File:              RTPLAN_{patient_id}_{session}_{plan_type}_B{orig}.dcm
        %   RTPlanLabel:       {session} {plan_type} B{orig}  (truncated to 16 chars, LO VR)
        %   RTPlanName:        {session} {plan_type} B{orig_beam_number}
        %   RTPlanDescription: {session} {plan_type} B{orig} ({beam_name}) exploded segments
        % -----------------------------------------------------------
        % RTPlanLabel is max 16 chars (DICOM LO VR); include session prefix truncated to fit
        label_full = sprintf('%s %s B%d', session, plan_type, original_beam_number);
        rtplan_out.RTPlanLabel       = label_full(1:min(16, end));
        if isfield(rtplan_out, 'RTPlanName')
            rtplan_out.RTPlanName    = sprintf('%s %s B%d', session, plan_type, original_beam_number);
        end
        rtplan_out.RTPlanDescription = sprintf('%s %s B%d (%s) exploded segments', ...
            session, plan_type, original_beam_number, original_beam_name);

        % -----------------------------------------------------------
        % Write output file   one file per original beam
        % -----------------------------------------------------------
        out_filename = sprintf('RTPLAN_%s_%s_%s_B%d.dcm', ...
            patient_id, session, plan_type, original_beam_number);
        out_filepath = fullfile(output_dir, out_filename);

        try
            dicomwrite([], out_filepath, rtplan_out, 'CreateMode', 'Copy');
        catch ME
            error('step06:writeFailed', ...
                'Failed to write output RTPLAN:\n  %s\nError: %s', out_filepath, ME.message);
        end

        plan_output_paths{end+1} = out_filepath; %#ok<AGROW>
        fprintf('    Written: %s\n\n', out_filename);

    end  % beam loop

    % -------------------------------------------------------------------
    % Per-plan-type summary
    % -------------------------------------------------------------------
    fprintf('%s plan: %d beam file(s) written.\n', upper(plan_type), numel(plan_output_paths));
    print_file_list(plan_output_paths);

    if any_warn_global
        warning('step06:muMismatch', ...
            '[%s] One or more beams had MU sum deviating >0.5%% from original.', plan_type);
    end

    output_paths.(plan_type) = plan_output_paths;
    fprintf('\n');

end  % plan type loop

% -----------------------------------------------------------------------
% Combine all paths for convenience
% -----------------------------------------------------------------------
output_paths.all = output_paths.reference;

total = numel(output_paths.all);
fprintf('=== Step 0.6 complete ===\n');
fprintf('Total beam files written: %d (reference plan only)\n', total);
fprintf('Output directory: %s\n\n', output_dir);

end  % function step06_explode_segments


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function copy_files_to_dir(src_dir, pattern, dst_dir, label)
%COPY_FILES_TO_DIR Copy files matching pattern from src_dir to dst_dir.
%   Skips files already present. Reports counts.
    files = dir(fullfile(src_dir, pattern));
    if isempty(files)
        warning('step06:noFiles', 'No %s files (%s) found in:\n  %s', ...
            label, pattern, src_dir);
        return;
    end
    n_copied  = 0;
    n_skipped = 0;
    for k = 1:numel(files)
        src = fullfile(src_dir, files(k).name);
        dst = fullfile(dst_dir, files(k).name);
        if isfile(dst)
            n_skipped = n_skipped + 1;
        else
            copyfile(src, dst);
            n_copied = n_copied + 1;
        end
    end
    fprintf('  %-20s %3d copied, %3d already present.\n', ...
        [label ':'], n_copied, n_skipped);
end


function print_file_list(paths)
%PRINT_FILE_LIST Print a summary list of written file names.
    n = numel(paths);
    if n == 0, return; end
    if n <= 20
        for k = 1:n
            [~, fn, ext] = fileparts(paths{k});
            fprintf('  %s%s\n', fn, ext);
        end
    else
        for k = 1:5
            [~, fn, ext] = fileparts(paths{k});
            fprintf('  %s%s\n', fn, ext);
        end
        fprintf('  ... (%d more files) ...\n', n - 10);
        for k = n-4 : n
            [~, fn, ext] = fileparts(paths{k});
            fprintf('  %s%s\n', fn, ext);
        end
    end
end