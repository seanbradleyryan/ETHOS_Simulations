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

% Derive output path
[~, base_name, ~] = fileparts(listing(1).name);
output_rtplan = fullfile(sct_dir, [base_name '_exploded_segments.dcm']);

fprintf('Input:  %s\n', input_rtplan);
fprintf('Output: %s\n\n', output_rtplan);

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
% Step 3 — Explode each beam into N single-segment beams
% -----------------------------------------------------------------------
new_beam_seq     = struct();
new_ref_beam_seq = struct();
global_beam_idx  = 0;

% For MU validation: accumulate per-original-beam sums
mu_sum_per_orig = containers.Map('KeyType', 'int32', 'ValueType', 'double');
mu_orig_per_beam = containers.Map('KeyType', 'int32', 'ValueType', 'double');

for b = 1:num_original_beams
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
        global_beam_idx = global_beam_idx + 1;

        % --- Deep-copy original beam (MATLAB structs are value types) ---
        new_beam = orig_beam;

        % Update identifying fields
        new_beam.BeamNumber = global_beam_idx;
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
            'BeamLimitingDevicePositionSequence', ...  % MLC/jaw positions
            'CumulativeMetersetWeight', ...            % already set to 0 above
            'ControlPointIndex'};                      % already set to 0 above

        cp1_all_fields = fieldnames(cp1_orig);
        for gf = 1:numel(cp1_all_fields)
            fld = cp1_all_fields{gf};
            if ~isfield(cp_entry_new, fld) && ~ismember(fld, cp1_varying_fields)
                cp_entry_new.(fld) = cp1_orig.(fld);
            end
        end

        % --- Build Item_2 (exit CP) ---
        cp_exit_new = cp_exit_orig;
        cp_exit_new.CumulativeMetersetWeight = 1;
        cp_exit_new.ControlPointIndex = 1;

        % Ensure MLC positions are present in exit CP (backfill if missing)
        if ~isfield(cp_exit_new, 'BeamLimitingDevicePositionSequence') && ...
                isfield(cp_entry_new, 'BeamLimitingDevicePositionSequence')
            cp_exit_new.BeamLimitingDevicePositionSequence = ...
                cp_entry_new.BeamLimitingDevicePositionSequence;
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
        new_beam_seq.(sprintf('Item_%d', global_beam_idx)) = new_beam;

        % --- Build FractionGroup reference entry ---
        ref_entry = struct();
        ref_entry.ReferencedBeamNumber = global_beam_idx;
        ref_entry.BeamMeterset         = seg_mu;
        new_ref_beam_seq.(sprintf('Item_%d', global_beam_idx)) = ref_entry;

        % Progress report every 50 new beams
        if mod(global_beam_idx, 50) == 0
            fprintf('  Created beam %d (orig B%d S%d, MU=%.4f)\n', ...
                global_beam_idx, original_beam_number, s - 1, seg_mu);
        end
    end % segment loop
end % beam loop

fprintf('\nTotal new beams created: %d\n\n', global_beam_idx);

%% -----------------------------------------------------------------------
% Step 4 — Assemble the output plan
% -----------------------------------------------------------------------
rtplan_out = rtplan;   % deep copy (struct assignment in MATLAB)

rtplan_out.BeamSequence = new_beam_seq;
rtplan_out.FractionGroupSequence.Item_1.ReferencedBeamSequence = new_ref_beam_seq;
rtplan_out.FractionGroupSequence.Item_1.NumberOfBeams = global_beam_idx;

% New SOP Instance UID
new_uid = dicomuid();
rtplan_out.SOPInstanceUID           = new_uid;
rtplan_out.MediaStorageSOPInstanceUID = new_uid;

% Plan label / description (LO max 16 chars)
rtplan_out.RTPlanLabel       = 'exploded_segs';   % 13 chars — safe
rtplan_out.RTPlanDescription = 'Per-segment explosion for RS import';

%% -----------------------------------------------------------------------
% Step 5 — Validate MU totals
% -----------------------------------------------------------------------
fprintf('MU Validation:\n');
any_warn = false;
orig_beam_nums = keys(mu_orig_per_beam);
for k = 1:numel(orig_beam_nums)
    bnum     = orig_beam_nums{k};
    orig_mu  = mu_orig_per_beam(bnum);
    summed   = mu_sum_per_orig(bnum);
    delta_abs = abs(summed - orig_mu);
    if orig_mu > 0
        delta_pct = 100 * delta_abs / orig_mu;
    else
        delta_pct = 0;
    end
    if delta_pct > 0.5
        status = '[WARN >0.5%]';
        any_warn = true;
    else
        status = '[OK]';
    end
    fprintf('  B%d: orig=%.4f  sum=%.4f  delta=%.4f%%  %s\n', ...
        bnum, orig_mu, summed, delta_pct, status);
end
if any_warn
    warning('step06:muMismatch', ...
        'One or more beams have MU sum deviating >0.5%% from original. Check output carefully.');
end

%% -----------------------------------------------------------------------
% Step 6 — Write output
% -----------------------------------------------------------------------
fprintf('\nWriting output DICOM...\n');
try
    dicomwrite([], output_rtplan, rtplan_out, 'CreateMode', 'Copy');
catch ME
    error('step06:writeFailed', ...
        'Failed to write output RTPLAN:\n  %s\nError: %s', output_rtplan, ME.message);
end

fprintf('Write complete: %s\n', output_rtplan);
fprintf('\n=== Step 0.6 complete. Total new beams: %d ===\n', global_beam_idx);
