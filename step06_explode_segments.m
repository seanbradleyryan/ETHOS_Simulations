%% step06_explode_segments.m
% Per-Segment Beam Explosion Script
%
% Takes MLC-gap-corrected RTPLAN files (_adjusted_mlc.dcm) for BOTH the
% reference and adapted plans and produces one RTPLAN file PER SEGMENT,
% each containing a single independent two-control-point beam.
%
% Processes two input files from the sct directory:
%   RP_reference_adjusted_mlc.dcm   -> "reference" plan set
%   RP_adapted_adjusted_mlc.dcm     -> "adapted"   plan set
%
% Output filename format:
%   plan_{patient_id}_{session}_{plantype}_B{orig_beam}_S{seg}.dcm
%   e.g.  plan_1194203_Session_1_reference_B10_S3.dcm
%         plan_1194203_Session_1_adapted_B10_S3.dcm
%
% Intra-DICOM naming (to distinguish plans in RayStation):
%   RTPlanLabel (max 16 chars): {short}_B{orig}_S{seg}
%                                where short = 'ref' | 'adp'
%                                e.g.  ref_B10_S3    adp_B10_S3
%   RTPlanName:                 {session}_{plantype}_B{orig}_S{seg}
%   RTPlanDescription:          {plantype} B{orig} ({beam_name}) segment {seg}
%
% Each output file has:
%   - BeamSequence with one beam  (the single 2-CP exploded segment)
%   - FractionGroupSequence with one referenced beam and its segment MU
%
% RayStation then exports one dose file per plan, which step15 loads via
% the preferred pattern:  Plan_Field*_Beam*_B*_S*.dcm
%
% Usage: run as a standalone script (no arguments).

clc;
fprintf('=== Step 0.6: Segment Explosion (Reference + Adapted) ===\n');

%% -----------------------------------------------------------------------
% Hardcoded defaults
% -----------------------------------------------------------------------
patient_id     = '1194203';
session        = 'Session_1';
working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
treatment_site = 'Pancreas';

fprintf('Patient: %s | Session: %s\n\n', patient_id, session);

%% -----------------------------------------------------------------------
% Paths
% -----------------------------------------------------------------------
sct_dir    = fullfile(working_dir, 'EthosExports', patient_id, treatment_site, session, 'sct');
output_dir = fullfile(working_dir, 'Raystation_Input', patient_id, session);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('SCT dir:    %s\n', sct_dir);
fprintf('Output dir: %s\n\n', output_dir);

%% -----------------------------------------------------------------------
% Define the two plans to process
%   label       - used in filenames and RTPlanDescription
%   short_label - used in RTPlanLabel (<=16 char budget)
%   filename    - expected _adjusted_mlc.dcm file in sct_dir
% -----------------------------------------------------------------------
plan_specs(1).label       = 'reference';
plan_specs(1).short_label = 'ref';
plan_specs(1).filename    = 'RP_reference_adjusted_mlc.dcm';

plan_specs(2).label       = 'adapted';
plan_specs(2).short_label = 'adp';
plan_specs(2).filename    = 'RP_adapted_adjusted_mlc.dcm';

%% -----------------------------------------------------------------------
% Grand totals
% -----------------------------------------------------------------------
all_output_paths = {};
any_warn_global  = false;

%% -----------------------------------------------------------------------
% Outer loop — one iteration per plan type
% -----------------------------------------------------------------------
for p = 1:numel(plan_specs)

    plan_label       = plan_specs(p).label;
    plan_short       = plan_specs(p).short_label;
    input_filename   = plan_specs(p).filename;
    input_rtplan     = fullfile(sct_dir, input_filename);

    fprintf('------------------------------------------------------------\n');
    fprintf('Plan type : %s\n', upper(plan_label));
    fprintf('Input file: %s\n', input_filename);

    % --- Check existence ---
    if ~isfile(input_rtplan)
        warning('step06:fileNotFound', ...
            'Input file not found for %s plan — skipping:\n  %s', ...
            plan_label, input_rtplan);
        fprintf('[SKIPPED] %s plan not found.\n\n', plan_label);
        continue;
    end

    fprintf('Filename pattern: plan_%s_%s_%s_B<orig>_S<seg>.dcm\n\n', ...
        patient_id, session, plan_label);

    % --- Load DICOM ---
    fprintf('Loading RTPLAN...\n');
    rtplan = dicominfo(input_rtplan, 'UseDictionaryVR', true);

    % Validate
    if ~isfield(rtplan, 'Modality') || ~strcmp(rtplan.Modality, 'RTPLAN')
        warning('step06:wrongModality', ...
            'File Modality is not RTPLAN for %s plan: %s — skipping.', ...
            plan_label, rtplan.Modality);
        continue;
    end
    if ~isfield(rtplan, 'BeamSequence')
        warning('step06:noBeams', ...
            'RTPLAN has no BeamSequence for %s plan — skipping.', plan_label);
        continue;
    end

    beam_fields_orig   = fieldnames(rtplan.BeamSequence);
    num_original_beams = numel(beam_fields_orig);
    fprintf('RTPlanLabel:         %s\n', rtplan.RTPlanLabel);
    fprintf('Original beam count: %d\n\n', num_original_beams);

    % --- Extract MU map  (beam_number -> total_MU) ---
    if ~isfield(rtplan, 'FractionGroupSequence')
        warning('step06:noFractionGroup', ...
            'RTPLAN has no FractionGroupSequence for %s plan — skipping.', plan_label);
        continue;
    end

    mu_map = containers.Map('KeyType', 'int32', 'ValueType', 'double');
    frac_group      = rtplan.FractionGroupSequence.Item_1;
    ref_beam_fields = fieldnames(frac_group.ReferencedBeamSequence);

    for k = 1:numel(ref_beam_fields)
        rb   = frac_group.ReferencedBeamSequence.(ref_beam_fields{k});
        bnum = int32(rb.ReferencedBeamNumber);
        bmu  = double(rb.BeamMeterset);
        mu_map(bnum) = bmu;
    end
    fprintf('MU map loaded for %d beams.\n\n', mu_map.Count);

    % --- Per-plan output accumulators ---
    plan_output_paths   = {};
    mu_sum_per_orig     = containers.Map('KeyType', 'int32', 'ValueType', 'double');
    mu_orig_per_beam    = containers.Map('KeyType', 'int32', 'ValueType', 'double');

    % -----------------------------------------------------------------------
    % Beam loop — explode each beam into per-segment RTPLAN files
    % -----------------------------------------------------------------------
    for b = 1:num_original_beams

        orig_beam            = rtplan.BeamSequence.(beam_fields_orig{b});
        original_beam_number = int32(orig_beam.BeamNumber);
        original_beam_name   = orig_beam.BeamName;

        if ~isKey(mu_map, original_beam_number)
            error('step06:muNotFound', ...
                'No MU found in FractionGroup for %s plan beam %d (%s).', ...
                plan_label, original_beam_number, original_beam_name);
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
                '%s plan beam %d (%s) has only %d control points — skipping.', ...
                plan_label, original_beam_number, original_beam_name, N_cp);
            continue;
        end

        cp1_orig = orig_beam.ControlPointSequence.(cp_fields{1});

        if isfield(cp1_orig, 'GantryAngle')
            gantry_str = sprintf('%.1f deg', double(cp1_orig.GantryAngle));
        else
            gantry_str = 'N/A';
        end
        fprintf('Processing [%s] B%d (%s, gantry %s, %d segments, MU=%.4f)...\n', ...
            plan_label, original_beam_number, original_beam_name, ...
            gantry_str, N_seg, total_mu);

        % ---------------------------------------------------------------
        % Segment loop
        % ---------------------------------------------------------------
        for s = 1:N_seg

            % Deep-copy original beam
            new_beam            = orig_beam;
            new_beam.BeamNumber = 1;
            new_beam.BeamName   = sprintf('B%d_S%03d', original_beam_number, s - 1);
            new_beam.NumberOfControlPoints = 2;

            cp_entry_orig = orig_beam.ControlPointSequence.(cp_fields{s});
            cp_exit_orig  = orig_beam.ControlPointSequence.(cp_fields{s + 1});

            % CMW values
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

            % Build entry CP (Item_1): all geometry fields present
            cp_entry_new = cp_entry_orig;
            cp_entry_new.CumulativeMetersetWeight = 0;
            cp_entry_new.ControlPointIndex        = 0;

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

            % BeamLimitingDevicePositionSequence for entry CP
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
            bldps_new_entry.Item_1 = bldps_cp1.Item_1;   % X jaw
            bldps_new_entry.Item_2 = bldps_cp1.Item_2;   % Y jaw
            bldps_new_entry.Item_3 = bldps_entry_mlc1;   % MLC1
            bldps_new_entry.Item_4 = bldps_entry_mlc2;   % MLC2
            cp_entry_new.BeamLimitingDevicePositionSequence = bldps_new_entry;

            % Build exit CP (Item_2): MLC-only
            cp_exit_new = cp_exit_orig;
            cp_exit_new.CumulativeMetersetWeight = 1;
            cp_exit_new.ControlPointIndex        = 1;

            if ~isfield(cp_exit_new, 'BeamLimitingDevicePositionSequence')
                bldps_exit_fb        = struct();
                bldps_exit_fb.Item_1 = bldps_entry_mlc1;
                bldps_exit_fb.Item_2 = bldps_entry_mlc2;
                cp_exit_new.BeamLimitingDevicePositionSequence = bldps_exit_fb;
            end

            % Assemble 2-item ControlPointSequence
            new_cp_seq        = struct();
            new_cp_seq.Item_1 = cp_entry_new;
            new_cp_seq.Item_2 = cp_exit_new;
            new_beam.ControlPointSequence         = new_cp_seq;
            new_beam.FinalCumulativeMetersetWeight = 1;

            % ------------------------------------------------------------------
            % Assemble per-segment output plan
            % ------------------------------------------------------------------
            rtplan_out = rtplan;

            rtplan_out.BeamSequence = struct('Item_1', new_beam);

            ref_entry                      = struct();
            ref_entry.ReferencedBeamNumber = 1;
            ref_entry.BeamMeterset         = seg_mu;

            rtplan_out.FractionGroupSequence.Item_1.ReferencedBeamSequence = ...
                struct('Item_1', ref_entry);
            rtplan_out.FractionGroupSequence.Item_1.NumberOfBeams = 1;

            % Unique SOP Instance UID
            new_uid = dicomuid();
            rtplan_out.SOPInstanceUID             = new_uid;
            rtplan_out.MediaStorageSOPInstanceUID = new_uid;

            % ------------------------------------------------------------------
            % Naming
            %
            %   File:             plan_{pid}_{session}_{plantype}_B{orig}_S{seg}.dcm
            %   RTPlanLabel:      {short}_B{orig}_S{seg}         (max 16 chars)
            %                       e.g. ref_B10_S3  /  adp_B10_S3
            %   RTPlanName:       {session}_{plantype}_B{orig}_S{seg}
            %   RTPlanDescription:{plantype} B{orig} ({beam_name}) segment {seg}
            % ------------------------------------------------------------------
            seg_label_short_str = sprintf('%s_B%d_S%d', plan_short, ...
                original_beam_number, s - 1);
            seg_label_short_str = seg_label_short_str(1:min(16, end));

            seg_label_long_str  = sprintf('%s_%s_B%d_S%d', session, plan_label, ...
                original_beam_number, s - 1);

            rtplan_out.RTPlanLabel       = seg_label_short_str;
            if isfield(rtplan_out, 'RTPlanName')
                rtplan_out.RTPlanName    = seg_label_long_str;
            end
            rtplan_out.RTPlanDescription = sprintf('%s B%d (%s) segment %d', ...
                plan_label, original_beam_number, original_beam_name, s - 1);

            % ------------------------------------------------------------------
            % Write output file
            % ------------------------------------------------------------------
            out_filename = sprintf('plan_%s_%s_%s_B%d_S%d.dcm', ...
                patient_id, session, plan_label, original_beam_number, s - 1);
            out_filepath = fullfile(output_dir, out_filename);

            try
                dicomwrite([], out_filepath, rtplan_out, 'CreateMode', 'Copy');
            catch ME
                error('step06:writeFailed', ...
                    'Failed to write output RTPLAN:\n  %s\nError: %s', ...
                    out_filepath, ME.message);
            end

            plan_output_paths{end+1} = out_filepath; %#ok<AGROW>

            % Progress every 50 new files
            total_written = numel(all_output_paths) + numel(plan_output_paths);
            if mod(total_written, 50) == 0
                fprintf('  Written %d files so far (all plans combined)...\n', total_written);
            end

        end % segment loop

        fprintf('  -> %d segment file(s) written for [%s] B%d\n', ...
            N_seg, plan_label, original_beam_number);

        % MU validation for this beam
        orig_mu   = mu_orig_per_beam(original_beam_number);
        summed    = mu_sum_per_orig(original_beam_number);
        delta_abs = abs(summed - orig_mu);
        delta_pct = 0;
        if orig_mu > 0
            delta_pct = 100 * delta_abs / orig_mu;
        end
        if delta_pct > 0.5
            mu_status = '[WARN >0.5%]';
            any_warn_global = true;
        else
            mu_status = '[OK]';
        end
        fprintf('  MU: orig=%.4f  sum=%.4f  delta=%.4f%%  %s\n\n', ...
            orig_mu, summed, delta_pct, mu_status);

    end % beam loop

    fprintf('[%s] plan complete: %d segment files written.\n\n', ...
        upper(plan_label), numel(plan_output_paths));

    all_output_paths = [all_output_paths, plan_output_paths]; %#ok<AGROW>

end % plan type loop

%% -----------------------------------------------------------------------
% Final summary
% -----------------------------------------------------------------------
fprintf('============================================================\n');
fprintf('=== Step 0.6 complete ===\n');
fprintf('Total segment files written: %d\n', numel(all_output_paths));
fprintf('Output directory: %s\n\n', output_dir);

if numel(all_output_paths) <= 20
    for k = 1:numel(all_output_paths)
        [~, fn, ext] = fileparts(all_output_paths{k});
        fprintf('  %s%s\n', fn, ext);
    end
else
    for k = 1:5
        [~, fn, ext] = fileparts(all_output_paths{k});
        fprintf('  %s%s\n', fn, ext);
    end
    fprintf('  ... (%d files total) ...\n', numel(all_output_paths));
    for k = numel(all_output_paths)-4:numel(all_output_paths)
        [~, fn, ext] = fileparts(all_output_paths{k});
        fprintf('  %s%s\n', fn, ext);
    end
end

if any_warn_global
    warning('step06:muMismatch', ...
        'One or more beams had MU sum deviating >0.5%% from original. Check output carefully.');
end