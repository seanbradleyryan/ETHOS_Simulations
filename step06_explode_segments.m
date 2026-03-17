%% step06_explode_segments.m
% Per-Segment Beam Explosion Script
%
% Takes a single-fraction RTPLAN that has already been MLC-gap-corrected
% (_adjusted_mlc.dcm) and produces a new RTPLAN in which every segment of
% every beam has been promoted into its own independent single-segment beam.
%
% All new beams are placed in one plan so RayStation can export one dose
% file per beam.
%
% Usage: run as a standalone script (no arguments).

clc;
fprintf('=== Step 0.6: Segment Explosion ===\n');

%% -----------------------------------------------------------------------
% Hardcoded defaults
% -----------------------------------------------------------------------
patient_id    = '1194203';
session       = 'Session_1';
working_dir   = '/mnt/weka/home/80030361/ETHOS_Simulations';
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

% Derive base name for per-beam output filenames
[~, base_name, ~] = fileparts(listing(1).name);

fprintf('Input:  %s\n', input_rtplan);
fprintf('Output: one file per beam -> %s_B<N>_exploded_segments.dcm\n\n', base_name);

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

beam_fields_orig = fieldnames(rtplan.BeamSequence);
num_original_beams = numel(beam_fields_orig);
fprintf('RTPlanLabel: %s\n', rtplan.RTPlanLabel);
fprintf('Original beam count: %d\n\n', num_original_beams);

%% -----------------------------------------------------------------------
% Step 2 — Extract original beam MU map  (beam_number -> total_MU)
% -----------------------------------------------------------------------
mu_map = containers.Map('KeyType', 'int32', 'ValueType', 'double');

if ~isfield(rtplan, 'FractionGroupSequence')
    error('step06:noFractionGroup', 'RTPLAN has no FractionGroupSequence.');
end
frac_group = rtplan.FractionGroupSequence.Item_1;
ref_beam_fields = fieldnames(frac_group.ReferencedBeamSequence);

for k = 1:numel(ref_beam_fields)
    rb = frac_group.ReferencedBeamSequence.(ref_beam_fields{k});
    bnum = int32(rb.ReferencedBeamNumber);
    bmu  = double(rb.BeamMeterset);
    mu_map(bnum) = bmu;
end

fprintf('MU map loaded for %d beams.\n\n', mu_map.Count);

%% -----------------------------------------------------------------------
% Step 3 — Explode each beam; one output RTPLAN file per original beam
% -----------------------------------------------------------------------
% For MU validation: accumulate per-original-beam sums
mu_sum_per_orig  = containers.Map('KeyType', 'int32', 'ValueType', 'double');
mu_orig_per_beam = containers.Map('KeyType', 'int32', 'ValueType', 'double');

output_paths = {};   % collect for final summary
any_warn_global = false;

for b = 1:num_original_beams

    % Per-beam accumulators (reset each iteration)
    new_beam_seq     = struct();
    new_ref_beam_seq = struct();
    local_beam_idx   = 0;
    orig_beam = rtplan.BeamSequence.(beam_fields_orig{b});

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

    % Control point fields (sorted so CP_1, CP_2 ... are in order)
    cp_fields = fieldnames(orig_beam.ControlPointSequence);
    cp_fields = sort(cp_fields);   % Item_1, Item_10, Item_100 ... need numeric sort
    % Numeric sort by the trailing integer
    cp_nums   = cellfun(@(f) sscanf(f, 'Item_%d'), cp_fields);
    [~, sort_idx] = sort(cp_nums);
    cp_fields = cp_fields(sort_idx);
    N_cp = numel(cp_fields);
    N_seg = N_cp - 1;

    if N_seg < 1
        warning('step06:noSegments', ...
            'Beam %d (%s) has only %d control points — skipping.', ...
            original_beam_number, original_beam_name, N_cp);
        continue;
    end

    % First control point of original beam (used for geometry back-fill)
    cp1_orig = orig_beam.ControlPointSequence.(cp_fields{1});

    % Gantry angle for progress display
    if isfield(cp1_orig, 'GantryAngle')
        gantry_str = sprintf('%.1f deg', double(cp1_orig.GantryAngle));
    else
        gantry_str = 'N/A';
    end
    fprintf('Processing beam B%d (%s, gantry %s, %d segments, total MU=%.4f)...\n', ...
        original_beam_number, original_beam_name, gantry_str, N_seg, total_mu);

    for s = 1:N_seg
        local_beam_idx = local_beam_idx + 1;

        % --- Deep-copy original beam (MATLAB structs are value types) ---
        new_beam = orig_beam;

        % Update identifying fields
        new_beam.BeamNumber = local_beam_idx;
        new_beam.BeamName   = sprintf('B%d_S%03d', original_beam_number, s - 1);
        new_beam.NumberOfControlPoints = 2;

        % --- Retrieve the two bounding control points ---
        cp_entry_orig = orig_beam.ControlPointSequence.(cp_fields{s});
        cp_exit_orig  = orig_beam.ControlPointSequence.(cp_fields{s + 1});

        % --- CMW values (guard for missing field) ---
        if isfield(cp_entry_orig, 'CumulativeMetersetWeight')
            cmw_entry = double(cp_entry_orig.CumulativeMetersetWeight);
        else
            % Reconstruct from linear spacing (fallback)
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

        % --- Build Item_1 (entry CP): ensure all geometry fields present ---
        cp_entry_new = cp_entry_orig;
        cp_entry_new.CumulativeMetersetWeight = 0;
        cp_entry_new.ControlPointIndex = 0;

        % Copy ALL fields from the original beam's first CP that are absent
        % from this entry CP, EXCEPT fields that legitimately vary between
        % control points and must come from the actual segment CP.
        cp1_varying_fields = { ...
            'BeamLimitingDevicePositionSequence', ...  % handled separately below
            'CumulativeMetersetWeight', ...            % already set to 0 above
            'ControlPointIndex'};                      % already set to 0 above

        cp1_all_fields = fieldnames(cp1_orig);
        for gf = 1:numel(cp1_all_fields)
            fld = cp1_all_fields{gf};
            if ~isfield(cp_entry_new, fld) && ~ismember(fld, cp1_varying_fields)
                cp_entry_new.(fld) = cp1_orig.(fld);
            end
        end

        % --- Build BeamLimitingDevicePositionSequence for entry CP ---
        % Entry CP (Item_1 of new beam) must have all four items:
        %   Item_1: X jaw  (from original beam CP_1 — jaws don't change between CPs)
        %   Item_2: Y jaw  (from original beam CP_1)
        %   Item_3: MLC1   (from this segment's actual entry CP)
        %   Item_4: MLC2   (from this segment's actual entry CP)
        %
        % The segment CPs only carry MLC positions (Item_1=MLC1, Item_2=MLC2),
        % so we splice the jaw items in from CP_1 of the original beam.
        bldps_cp1   = cp1_orig.BeamLimitingDevicePositionSequence;
        bldps_entry = cp_entry_orig.BeamLimitingDevicePositionSequence;

        % Determine which items in the entry CP's BLDPS are the MLCs.
        % CP_1 of the original beam has 4 items: jaw, jaw, MLC1, MLC2.
        % All other CPs have 2 items: MLC1, MLC2.
        n_bldps_entry = numel(fieldnames(bldps_entry));
        if n_bldps_entry == 4
            % Entry CP is the original beam's first CP — MLCs are Item_3 & Item_4
            bldps_entry_mlc1 = bldps_entry.Item_3;
            bldps_entry_mlc2 = bldps_entry.Item_4;
        else
            % Entry CP is an intermediate CP — MLCs are Item_1 & Item_2
            bldps_entry_mlc1 = bldps_entry.Item_1;
            bldps_entry_mlc2 = bldps_entry.Item_2;
        end

        bldps_new_entry = struct();
        bldps_new_entry.Item_1 = bldps_cp1.Item_1;  % X jaw (always from original CP_1)
        bldps_new_entry.Item_2 = bldps_cp1.Item_2;  % Y jaw (always from original CP_1)
        bldps_new_entry.Item_3 = bldps_entry_mlc1;  % MLC1 from this segment's entry
        bldps_new_entry.Item_4 = bldps_entry_mlc2;  % MLC2 from this segment's entry
        cp_entry_new.BeamLimitingDevicePositionSequence = bldps_new_entry;

        % --- Build Item_2 (exit CP) ---
        % Exit CP only needs MLC positions (Item_1=MLC1, Item_2=MLC2).
        % Jaws are not repeated on intermediate/exit CPs.
        cp_exit_new = cp_exit_orig;
        cp_exit_new.CumulativeMetersetWeight = 1;
        cp_exit_new.ControlPointIndex = 1;

        % Ensure MLC positions are present in exit CP (backfill from entry if missing)
        if ~isfield(cp_exit_new, 'BeamLimitingDevicePositionSequence')
            bldps_exit_fallback = struct();
            bldps_exit_fallback.Item_1 = bldps_entry_mlc1;
            bldps_exit_fallback.Item_2 = bldps_entry_mlc2;
            cp_exit_new.BeamLimitingDevicePositionSequence = bldps_exit_fallback;
        end

        % --- Assign new 2-item ControlPointSequence ---
        new_cp_seq = struct();
        new_cp_seq.Item_1 = cp_entry_new;
        new_cp_seq.Item_2 = cp_exit_new;
        new_beam.ControlPointSequence = new_cp_seq;

        % --- FinalCumulativeMetersetWeight ---
        % Should equal the exit CMW of the original beam at this level;
        % conventionally set to 1 for a standalone beam.
        new_beam.FinalCumulativeMetersetWeight = 1;

        % --- Store new beam ---
        new_beam_seq.(sprintf('Item_%d', local_beam_idx)) = new_beam;

        % --- Build FractionGroup reference entry ---
        ref_entry = struct();
        ref_entry.ReferencedBeamNumber = local_beam_idx;
        ref_entry.BeamMeterset         = seg_mu;
        new_ref_beam_seq.(sprintf('Item_%d', local_beam_idx)) = ref_entry;

        % Progress report every 50 new beams
        if mod(local_beam_idx, 50) == 0
            fprintf('  Created beam %d (orig B%d S%d, MU=%.4f)\n', ...
                local_beam_idx, original_beam_number, s - 1, seg_mu);
        end
    end % segment loop

    fprintf('  -> %d segments exploded for beam B%d\n', local_beam_idx, original_beam_number);

    %% -------------------------------------------------------------------
    % Step 4 — Assemble per-beam output plan
    % -------------------------------------------------------------------
    rtplan_out = rtplan;   % deep copy

    rtplan_out.BeamSequence = new_beam_seq;
    rtplan_out.FractionGroupSequence.Item_1.ReferencedBeamSequence = new_ref_beam_seq;
    rtplan_out.FractionGroupSequence.Item_1.NumberOfBeams = local_beam_idx;

    % New SOP Instance UID (unique per file)
    new_uid = dicomuid();
    rtplan_out.SOPInstanceUID             = new_uid;
    rtplan_out.MediaStorageSOPInstanceUID = new_uid;

    % Plan label: e.g. 'exp_B1' — max 16 chars
    rtplan_out.RTPlanLabel       = sprintf('exp_B%d', original_beam_number);
    rtplan_out.RTPlanDescription = sprintf('Exploded segments beam B%d', original_beam_number);

    %% -------------------------------------------------------------------
    % Step 5 — Validate MU total for this beam
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
    fprintf('  MU: orig=%.4f  sum=%.4f  delta=%.4f%%  %s\n', ...
        orig_mu, summed, delta_pct, mu_status);

    %% -------------------------------------------------------------------
    % Step 6 — Write per-beam output file
    % -------------------------------------------------------------------
    beam_output_path = fullfile(sct_dir, ...
        sprintf('%s_B%d_exploded_segments.dcm', base_name, original_beam_number));

    fprintf('  Writing: %s\n', beam_output_path);
    try
        dicomwrite([], beam_output_path, rtplan_out, 'CreateMode', 'Copy');
    catch ME
        error('step06:writeFailed', ...
            'Failed to write output RTPLAN:\n  %s\nError: %s', beam_output_path, ME.message);
    end
    output_paths{end+1} = beam_output_path; %#ok<AGROW>

end % beam loop

%% -----------------------------------------------------------------------
% Final summary
% -----------------------------------------------------------------------
fprintf('\n=== Step 0.6 complete ===\n');
fprintf('Output files written (%d plans):\n', numel(output_paths));
for k = 1:numel(output_paths)
    fprintf('  %s\n', output_paths{k});
end
if any_warn_global
    warning('step06:muMismatch', ...
        'One or more beams had MU sum deviating >0.5%% from original. Check output carefully.');
end
