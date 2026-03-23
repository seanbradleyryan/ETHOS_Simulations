%% step06_explode_segments.m
% Per-Segment Beam Explosion Script
%
% Takes a single-fraction RTPLAN that has already been MLC-gap-corrected
% (_adjusted_mlc.dcm) and produces one RTPLAN file PER SEGMENT, each
% containing a single independent two-control-point beam.
%
% Output filename format:  plan_{patient_id}_{session}_B{orig_beam}_S{seg}.dcm
%   e.g.  plan_1194203_Session_1_B10_S3.dcm
%
% Intra-DICOM naming (omits 'plan' and patient_id to reduce clutter):
%   RTPlanLabel (max 16 chars): {session}_B{orig}_S{seg}  (truncated if needed)
%   RTPlanName:                 {session}_B{orig}_S{seg}
%   RTPlanDescription:          B{orig} ({beam_name}) segment {seg}
%
% Each output file has:
%   - BeamSequence with one beam  (the single 2-CP exploded segment)
%   - FractionGroupSequence with one referenced beam and its segment MU
%
% RayStation then exports one dose file per plan, which step15 loads via
% the preferred pattern:  plan_{id}_{session}_B{orig}_S{seg}.dcm
%
% Usage: run as a standalone script (no arguments).

clc;
fprintf('=== Step 0.6: Segment Explosion ===\n');

%% -----------------------------------------------------------------------
% Hardcoded defaults
% -----------------------------------------------------------------------
patient_id     = '1194203';
session        = 'Session_1';
working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
treatment_site = 'Pancreas';

fprintf('Patient: %s | Session: %s\n', patient_id, session);

%% -----------------------------------------------------------------------
% Step 1 — Derive paths and load the adjusted RTPLAN
% -----------------------------------------------------------------------
sct_dir = fullfile(working_dir, 'EthosExports', patient_id, treatment_site, session, 'sct');

% Find the _adjusted_mlc.dcm file
listing = dir(fullfile(sct_dir, '*_adjusted_mlc.dcm'));
if isempty(listing)
    error('step06:fileNotFound', ...
        'No *_adjusted_mlc.dcm file found in:\n  %s', sct_dir);
end
if numel(listing) > 1
    warning('step06:multipleFiles', ...
        'Multiple *_adjusted_mlc.dcm files found; using first: %s', listing(1).name);
end

input_rtplan = fullfile(sct_dir, listing(1).name);
[~, base_name, ~] = fileparts(listing(1).name);

% Output directory for RayStation import
output_dir = fullfile(working_dir, 'Raystation_Input', patient_id, session);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Input:      %s\n', input_rtplan);
fprintf('Output dir: %s\n', output_dir);
fprintf('Filename pattern: plan_%s_%s_B<orig>_S<seg>.dcm\n\n', patient_id, session);

% Load DICOM
fprintf('Loading RTPLAN...\n');
rtplan = dicominfo(input_rtplan, 'UseDictionaryVR', true);

% Validate
if ~isfield(rtplan, 'Modality') || ~strcmp(rtplan.Modality, 'RTPLAN')
    error('step06:wrongModality', 'File Modality is not RTPLAN: %s', rtplan.Modality);
end
if ~isfield(rtplan, 'BeamSequence')
    error('step06:noBeams', 'RTPLAN has no BeamSequence.');
end

beam_fields_orig    = fieldnames(rtplan.BeamSequence);
num_original_beams  = numel(beam_fields_orig);
fprintf('RTPlanLabel:        %s\n', rtplan.RTPlanLabel);
fprintf('Original beam count: %d\n\n', num_original_beams);

%% -----------------------------------------------------------------------
% Step 2 — Extract original beam MU map  (beam_number -> total_MU)
% -----------------------------------------------------------------------
mu_map = containers.Map('KeyType', 'int32', 'ValueType', 'double');

if ~isfield(rtplan, 'FractionGroupSequence')
    error('step06:noFractionGroup', 'RTPLAN has no FractionGroupSequence.');
end
frac_group      = rtplan.FractionGroupSequence.Item_1;
ref_beam_fields = fieldnames(frac_group.ReferencedBeamSequence);

for k = 1:numel(ref_beam_fields)
    rb   = frac_group.ReferencedBeamSequence.(ref_beam_fields{k});
    bnum = int32(rb.ReferencedBeamNumber);
    bmu  = double(rb.BeamMeterset);
    mu_map(bnum) = bmu;
end

fprintf('MU map loaded for %d beams.\n\n', mu_map.Count);

%% -----------------------------------------------------------------------
% Step 3 — Explode each beam; one output RTPLAN file per segment
% -----------------------------------------------------------------------
output_paths     = {};   % collect for final summary
any_warn_global  = false;

% Per-original-beam MU validation accumulators
mu_sum_per_orig  = containers.Map('KeyType', 'int32', 'ValueType', 'double');
mu_orig_per_beam = containers.Map('KeyType', 'int32', 'ValueType', 'double');

for b = 1:num_original_beams

    orig_beam            = rtplan.BeamSequence.(beam_fields_orig{b});
    original_beam_number = int32(orig_beam.BeamNumber);
    original_beam_name   = orig_beam.BeamName;

    % Retrieve total MU for this beam
    if ~isKey(mu_map, original_beam_number)
        error('step06:muNotFound', ...
            'No MU found in FractionGroup for beam %d (%s).', ...
            original_beam_number, original_beam_name);
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
            'Beam %d (%s) has only %d control points — skipping.', ...
            original_beam_number, original_beam_name, N_cp);
        continue;
    end

    % First control point of original beam (geometry back-fill source)
    cp1_orig = orig_beam.ControlPointSequence.(cp_fields{1});

    % Gantry angle for progress display
    if isfield(cp1_orig, 'GantryAngle')
        gantry_str = sprintf('%.1f deg', double(cp1_orig.GantryAngle));
    else
        gantry_str = 'N/A';
    end
    fprintf('Processing beam B%d (%s, gantry %s, %d segments, total MU=%.4f)...\n', ...
        original_beam_number, original_beam_name, gantry_str, N_seg, total_mu);

    % ------------------------------------------------------------------
    % Loop over segments — each segment becomes its own output RTPLAN
    % ------------------------------------------------------------------
    for s = 1:N_seg

        % --- Deep-copy original beam ---
        new_beam            = orig_beam;
        new_beam.BeamNumber = 1;   % single beam in each plan -> always 1
        new_beam.BeamName   = sprintf('B%d_S%03d', original_beam_number, s - 1);
        new_beam.NumberOfControlPoints = 2;

        % --- Retrieve the two bounding control points ---
        cp_entry_orig = orig_beam.ControlPointSequence.(cp_fields{s});
        cp_exit_orig  = orig_beam.ControlPointSequence.(cp_fields{s + 1});

        % --- CMW values ---
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

        % --- Segment MU ---
        seg_mu = (cmw_exit - cmw_entry) * total_mu;
        mu_sum_per_orig(original_beam_number) = ...
            mu_sum_per_orig(original_beam_number) + seg_mu;

        % --- Build entry CP (Item_1): all geometry fields present ---
        cp_entry_new = cp_entry_orig;
        cp_entry_new.CumulativeMetersetWeight = 0;
        cp_entry_new.ControlPointIndex        = 0;

        % Back-fill any geometry fields absent from this intermediate CP
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

        % --- BeamLimitingDevicePositionSequence for entry CP ---
        % Item_1: X jaw  (from original beam CP_1)
        % Item_2: Y jaw  (from original beam CP_1)
        % Item_3: MLC1   (from this segment's actual entry CP)
        % Item_4: MLC2   (from this segment's actual entry CP)
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

        % --- Build exit CP (Item_2): MLC-only ---
        cp_exit_new = cp_exit_orig;
        cp_exit_new.CumulativeMetersetWeight = 1;
        cp_exit_new.ControlPointIndex        = 1;

        if ~isfield(cp_exit_new, 'BeamLimitingDevicePositionSequence')
            bldps_exit_fb        = struct();
            bldps_exit_fb.Item_1 = bldps_entry_mlc1;
            bldps_exit_fb.Item_2 = bldps_entry_mlc2;
            cp_exit_new.BeamLimitingDevicePositionSequence = bldps_exit_fb;
        end

        % --- Assign 2-item ControlPointSequence ---
        new_cp_seq        = struct();
        new_cp_seq.Item_1 = cp_entry_new;
        new_cp_seq.Item_2 = cp_exit_new;
        new_beam.ControlPointSequence           = new_cp_seq;
        new_beam.FinalCumulativeMetersetWeight   = 1;

        % ------------------------------------------------------------------
        % Assemble per-segment output plan
        % ------------------------------------------------------------------
        rtplan_out = rtplan;   % deep copy from original

        % BeamSequence: single beam
        rtplan_out.BeamSequence = struct('Item_1', new_beam);

        % FractionGroupSequence: single referenced beam
        ref_entry                     = struct();
        ref_entry.ReferencedBeamNumber = 1;
        ref_entry.BeamMeterset          = seg_mu;

        rtplan_out.FractionGroupSequence.Item_1.ReferencedBeamSequence = ...
            struct('Item_1', ref_entry);
        rtplan_out.FractionGroupSequence.Item_1.NumberOfBeams = 1;

        % --- Unique SOP Instance UID ---
        new_uid = dicomuid();
        rtplan_out.SOPInstanceUID             = new_uid;
        rtplan_out.MediaStorageSOPInstanceUID = new_uid;

        % ------------------------------------------------------------------
        % Naming
        %   File:             plan_{patient_id}_{session}_B{orig}_S{seg}.dcm
        %   RTPlanLabel:      {session}_B{orig}_S{seg}  (max 16 chars)
        %   RTPlanName:       {session}_B{orig}_S{seg}
        %   RTPlanDescription: B{orig} ({beam_name}) segment {seg}
        % ------------------------------------------------------------------
        seg_label_long  = sprintf('%s_B%d_S%d', session, original_beam_number, s - 1);
        seg_label_short = seg_label_long(1:min(16, end));  % LO max 16 chars

        rtplan_out.RTPlanLabel       = seg_label_short;
        if isfield(rtplan_out, 'RTPlanName')
            rtplan_out.RTPlanName    = seg_label_long;
        end
        rtplan_out.RTPlanDescription = sprintf('B%d (%s) segment %d', ...
            original_beam_number, original_beam_name, s - 1);

        % ------------------------------------------------------------------
        % Write output file
        % ------------------------------------------------------------------
        out_filename = sprintf('plan_%s_%s_B%d_S%d.dcm', ...
            patient_id, session, original_beam_number, s - 1);
        out_filepath = fullfile(output_dir, out_filename);

        try
            dicomwrite([], out_filepath, rtplan_out, 'CreateMode', 'Copy');
        catch ME
            error('step06:writeFailed', ...
                'Failed to write output RTPLAN:\n  %s\nError: %s', out_filepath, ME.message);
        end

        output_paths{end+1} = out_filepath; %#ok<AGROW>

        % Progress every 50 new files
        if mod(numel(output_paths), 50) == 0
            fprintf('  Written %d files so far...\n', numel(output_paths));
        end

    end % segment loop

    fprintf('  -> %d segment file(s) written for beam B%d\n', N_seg, original_beam_number);

    %% -------------------------------------------------------------------
    % MU validation for this original beam
    % -------------------------------------------------------------------
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

%% -----------------------------------------------------------------------
% Final summary
% -----------------------------------------------------------------------
fprintf('=== Step 0.6 complete ===\n');
fprintf('Total segment files written: %d\n', numel(output_paths));
fprintf('Output directory: %s\n\n', output_dir);

if numel(output_paths) <= 20
    % Print all paths for small plans
    for k = 1:numel(output_paths)
        [~, fn, ext] = fileparts(output_paths{k});
        fprintf('  %s%s\n', fn, ext);
    end
else
    % Print first and last few for large plans
    for k = 1:5
        [~, fn, ext] = fileparts(output_paths{k});
        fprintf('  %s%s\n', fn, ext);
    end
    fprintf('  ... (%d files) ...\n', numel(output_paths) - 10);
    for k = numel(output_paths)-4:numel(output_paths)
        [~, fn, ext] = fileparts(output_paths{k});
        fprintf('  %s%s\n', fn, ext);
    end
end

if any_warn_global
    warning('step06:muMismatch', ...
        'One or more beams had MU sum deviating >0.5%% from original. Check output carefully.');
end
