function swapped = step055_swap_imageset(patient_id, session, config)
%% STEP055_SWAP_IMAGESET - Re-associate MLC-corrected RTPLAN onto ICBCT and VCBCT
%
%   swapped = step055_swap_imageset(patient_id, session, config)
%
%   PURPOSE:
%   Takes the MLC-gap-corrected RTPLAN produced by step05_fix_mlc_gaps (planned
%   on the synthetic CT, SCT) and produces TWO new RTPLAN files re-associated
%   with the ICBCT and VCBCT image sets respectively.  Beam-spatial coordinates
%   (isocenters, surface-entry points, beam-dose-specification points) are
%   transformed via the spatial registration object (SRO) so that RayStation,
%   on import, treats each plan as a nominal-dose-capable plan rooted on the
%   target CBCT image set.  Step 0.6 then explodes both swapped plans and
%   skips the original SCT plan.
%
%   INPUTS:
%       patient_id  - char | string, patient identifier (e.g., '1194203')
%       session     - char | string, session name (e.g., 'Session_1')
%       config      - struct with configuration parameters:
%           .working_dir       - Base directory path
%           .treatment_site    - Subfolder name (default: 'Pancreas')
%           .swap_source_plan  - (optional) filename of source RTPLAN in
%                                sct_dir.  Default: 'RTPLAN_reference_adjusted_mlc.dcm'
%           .swap_rotation_tol_deg - (optional) Euler-angle tolerance for
%                                rotational-registration warning. Default: 1.0
%           .swap_coord_warn_mm    - (optional) magnitude threshold for the
%                                transformed-coordinate sanity warning.
%                                Default: 5000
%
%   OUTPUTS:
%       swapped - struct with fields:
%           .icbct      - char, full path to the ICBCT-associated RTPLAN
%                         (empty if generation failed)
%           .vcbct      - char, full path to the VCBCT-associated RTPLAN
%                         (empty if generation failed)
%           .all        - cell array of char, the plans step 0.6 should explode
%                         {icbct, vcbct} (only successful outputs included)
%           .sct_source - char, full path to the original SCT-rooted RTPLAN
%                         (the plan step 0.6 must IGNORE)
%
%   EXPECTED INPUT FILES (resolved from patient_id, session, config):
%       <working_dir>/EthosExports/<patient_id>/<treatment_site>/<session>/sct/
%           RTPLAN_reference_adjusted_mlc.dcm  - source plan (step 0.5 output)
%           CBCT1_*.dcm                        - ICBCT image slices
%           CBCT2_*.dcm                        - VCBCT image slices
%           RTSTRUCT_ICBCT.dcm                 - ICBCT structure set
%           RTSTRUCT_VCBCT.dcm                 - VCBCT structure set
%           REG_*.dcm                          - one or more SROs (step 0 copy)
%
%   OUTPUT FILES (written to step 0.6's input directory):
%       <working_dir>/Raystation_Input/<patient_id>/<session>/
%           RTPLAN_<patient>_<session>_reference_on_ICBCT.dcm
%           RTPLAN_<patient>_<session>_reference_on_VCBCT.dcm
%
%   ALGORITHM:
%   1. Resolve paths for source plan, target CBCT image sets, target RTSTRUCTs
%      and REG_*.dcm SRO files.
%   2. Read source plan and target frame-of-reference UIDs from the first
%      slice of each CBCT series.
%   3. For each REG file, parse RegistrationSequence items; build a list of
%      (FrameOfReferenceUID -> 4x4) entries.  Within each SRO one item has
%      the identity transform (fixed frame); the other item carries the
%      moving-frame transform.  Per the DICOM C.20.2 standard the matrix is
%      16 doubles, row-major; in MATLAB use reshape(v,4,4).' (default reshape
%      is column-major - the transpose is required).
%   4. VCBCT output: select the SRO whose MOVING FrameOfReferenceUID equals
%      the SCT FoR and whose FIXED FrameOfReferenceUID equals the VCBCT FoR.
%      Apply that 4x4 directly.
%   5. ICBCT output: compare the SCT and ICBCT FoR UIDs.
%        - If identical (ICBCT is typically the source for the synthetic CT,
%          so they often share a FoR), use the 4x4 identity (UID rewrite only).
%        - Otherwise derive T_sct_to_icbct = inv(T_icbct_to_vcbct) * T_sct_to_vcbct
%          using a second SRO; if a second SRO is not present, error out.
%   6. Apply the chosen 4x4 to every beam:
%        - ControlPointSequence.Item_N.IsocenterPosition
%        - ControlPointSequence.Item_N.SurfaceEntryPoint (if present)
%        - FractionGroupSequence.Item_N.ReferencedBeamSequence.Item_M
%              .BeamDoseSpecificationPoint (if present)
%      Recursively scan for any other *Position / *Point 3-vec under
%      BeamSequence and warn (do NOT silently transform).
%      Leave SourceToSurfaceDistance, LeafJawPositions, and gantry / collimator
%      / patient-support angles untouched.
%   7. Rewrite UIDs and references to root the plan on the target CT:
%        - Fresh via dicomuid: SOPInstanceUID, MediaStorageSOPInstanceUID,
%          SeriesInstanceUID
%        - From target CT first slice: StudyInstanceUID, FrameOfReferenceUID,
%          ReferencedFrameOfReferenceSequence.Item_1.FrameOfReferenceUID,
%          nested ...RTReferencedStudySequence.Item_1
%                       .RTReferencedSeriesSequence.Item_1.SeriesInstanceUID
%        - From target RTSTRUCT:
%          ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID
%        - RTPlanLabel: append '_on_ICBCT' / '_on_VCBCT' suffix (truncated
%          to DICOM 16-char limit)
%        - InstanceCreationDate / InstanceCreationTime: now
%        - PatientID / PatientName / StudyID are LEFT ALONE
%   8. Write with dicomwrite(s, path, 'CreateMode','Copy', 'WritePrivate', true).
%      CreateMode='Copy' is critical - the default ('Create') strips unknown
%      RTPLAN private fields and breaks RayStation import.
%
%   COORDINATES:
%   RTPLAN spatial fields are DICOM patient-coordinate (x,y,z) in mm.  No
%   MATLAB-array (row, col, slice) transposition applies here - we are not
%   reading a volume, only transforming sparse 3-vectors.  The CLAUDE.md
%   array-axis convention is irrelevant to this step.
%
%   LIMITATIONS:
%   - Rigid registration only.  Deformable-registration objects (DRO) are
%     not resampled - if a REG file carries a DeformableRegistrationSequence
%     with no MatrixRegistrationSequence it will be skipped with a warning.
%   - Gantry / collimator / patient-support angles are NOT rotated even when
%     the registration matrix carries a rotation.  For dose-on-CBCT this is
%     conventionally correct (the linac frame is unchanged); reviewer
%     confirms acceptance prior to clinical use.
%   - No dosimetric validation - RayStation re-computes dose on import.
%
%   EXAMPLE:
%       config.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       swapped = step055_swap_imageset('1194203', 'Session_1', config);
%       % swapped.icbct, swapped.vcbct -- ready for step06
%       % step06 should iterate over swapped.all and skip swapped.sct_source
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (dicominfo, dicomwrite)
%       - dicomuid() for generating new UIDs
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: May 2026
%   VERSION: 1.0
%
%   See also: step05_fix_mlc_gaps, step06_explode_segments, dicominfo,
%             dicomwrite, dicomuid

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('step055_swap_imageset:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

if ~ischar(session) && ~isstring(session)
    error('step055_swap_imageset:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

if ~isstruct(config)
    error('step055_swap_imageset:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

if ~isfield(config, 'working_dir')
    error('step055_swap_imageset:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

if ~isfolder(config.working_dir)
    error('step055_swap_imageset:DirectoryNotFound', ...
        'Working directory does not exist: %s', config.working_dir);
end

%% ======================== SET DEFAULT CONFIG VALUES ========================

if ~isfield(config, 'treatment_site') || isempty(config.treatment_site)
    config.treatment_site = 'Pancreas';
end

if ~isfield(config, 'swap_source_plan') || isempty(config.swap_source_plan)
    config.swap_source_plan = 'RTPLAN_reference_adjusted_mlc.dcm';
end

if ~isfield(config, 'swap_rotation_tol_deg') || isempty(config.swap_rotation_tol_deg)
    config.swap_rotation_tol_deg = 1.0;
end

if ~isfield(config, 'swap_coord_warn_mm') || isempty(config.swap_coord_warn_mm)
    config.swap_coord_warn_mm = 5000;
end

ROT_TOL_DEG    = config.swap_rotation_tol_deg;
COORD_WARN_MM  = config.swap_coord_warn_mm;

%% ======================== CONSTRUCT PATHS ========================

sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session, 'sct');

output_dir = fullfile(config.working_dir, 'Raystation_Input', patient_id, session);

fprintf('\n  [STEP 0.55] Image-set swap: SCT -> ICBCT + VCBCT\n');
fprintf('  [STEP 0.55] Patient: %s | Session: %s\n', patient_id, session);
fprintf('  [STEP 0.55] SCT dir   : %s\n', sct_dir);
fprintf('  [STEP 0.55] Output dir: %s\n', output_dir);

if ~isfolder(sct_dir)
    error('step055_swap_imageset:DirectoryNotFound', ...
        'SCT directory not found: %s\nRun step0_sort_dicom first.', sct_dir);
end

if ~isfolder(output_dir)
    mkdir(output_dir);
    fprintf('  [STEP 0.55] Created output dir.\n');
end

%% ======================== INITIALIZE OUTPUT ========================

source_plan_path = fullfile(sct_dir, config.swap_source_plan);

swapped = struct( ...
    'icbct',      '', ...
    'vcbct',      '', ...
    'all',        {{}}, ...
    'sct_source', source_plan_path);

%% ======================== LOCATE INPUTS ========================

if ~isfile(source_plan_path)
    error('step055_swap_imageset:MissingSourcePlan', ...
        'Source RTPLAN not found: %s\nRun step05_fix_mlc_gaps first.', source_plan_path);
end
fprintf('  [STEP 0.55] Source plan : %s\n', config.swap_source_plan);

% Target CT first-slice paths (used to read FoR / Study / Series UIDs)
icbct_ct_path = locateFirstSlice(sct_dir, 'CBCT1_*.dcm', 'ICBCT (CBCT1_)');
vcbct_ct_path = locateFirstSlice(sct_dir, 'CBCT2_*.dcm', 'VCBCT (CBCT2_)');

% Target RTSTRUCTs
icbct_rs_path = fullfile(sct_dir, 'RTSTRUCT_ICBCT.dcm');
vcbct_rs_path = fullfile(sct_dir, 'RTSTRUCT_VCBCT.dcm');
if ~isfile(icbct_rs_path)
    error('step055_swap_imageset:MissingRTSTRUCT', ...
        'ICBCT RTSTRUCT not found: %s', icbct_rs_path);
end
if ~isfile(vcbct_rs_path)
    error('step055_swap_imageset:MissingRTSTRUCT', ...
        'VCBCT RTSTRUCT not found: %s', vcbct_rs_path);
end

% Registration files
reg_listing = dir(fullfile(sct_dir, 'REG_*.dcm'));
if isempty(reg_listing)
    error('step055_swap_imageset:MissingREG', ...
        'No REG_*.dcm spatial registration files found in %s', sct_dir);
end
fprintf('  [STEP 0.55] Found %d REG_*.dcm file(s):\n', numel(reg_listing));
for i = 1:numel(reg_listing)
    fprintf('              - %s\n', reg_listing(i).name);
end

%% ======================== READ FoR UIDs ========================

src_plan_meta = dicominfo(source_plan_path);
sct_for_uid   = getFieldOrEmpty(src_plan_meta, 'FrameOfReferenceUID');
if isempty(sct_for_uid)
    error('step055_swap_imageset:NoFoR', ...
        'Source RTPLAN has no FrameOfReferenceUID: %s', source_plan_path);
end

icbct_ct_meta = dicominfo(icbct_ct_path);
vcbct_ct_meta = dicominfo(vcbct_ct_path);
icbct_for_uid = getFieldOrEmpty(icbct_ct_meta, 'FrameOfReferenceUID');
vcbct_for_uid = getFieldOrEmpty(vcbct_ct_meta, 'FrameOfReferenceUID');

if isempty(icbct_for_uid) || isempty(vcbct_for_uid)
    error('step055_swap_imageset:NoFoR', ...
        'Target CT slices are missing FrameOfReferenceUID.');
end

fprintf('  [STEP 0.55] SCT   FrameOfReferenceUID: %s\n', sct_for_uid);
fprintf('  [STEP 0.55] ICBCT FrameOfReferenceUID: %s\n', icbct_for_uid);
fprintf('  [STEP 0.55] VCBCT FrameOfReferenceUID: %s\n', vcbct_for_uid);

%% ======================== PARSE ALL REGISTRATIONS ========================

% Aggregate every (moving_for, fixed_for, 4x4) tuple found across all REG files.
reg_entries = struct('moving', {}, 'fixed', {}, 'matrix', {}, 'source_file', {});

for ri = 1:numel(reg_listing)
    reg_path = fullfile(sct_dir, reg_listing(ri).name);
    try
        new_entries = parseRegistration(reg_path);
    catch ME
        warning('step055_swap_imageset:RegParseError', ...
            'Failed to parse %s: %s', reg_listing(ri).name, ME.message);
        continue;
    end
    for ei = 1:numel(new_entries)
        new_entries(ei).source_file = reg_listing(ri).name; %#ok<AGROW>
        reg_entries(end+1) = new_entries(ei); %#ok<AGROW>
        fprintf('  [STEP 0.55] REG entry: moving=%s  fixed=%s  (from %s)\n', ...
            new_entries(ei).moving, new_entries(ei).fixed, reg_listing(ri).name);
    end
end

if isempty(reg_entries)
    error('step055_swap_imageset:NoUsableRegistration', ...
        'No usable (non-identity) MatrixRegistrationSequence entries found.');
end

%% ======================== BUILD T_SCT_TO_VCBCT ========================

T_sct_to_vcbct = findTransform(reg_entries, sct_for_uid, vcbct_for_uid);
if isempty(T_sct_to_vcbct)
    error('step055_swap_imageset:NoSCTtoVCBCT', ...
        ['No registration found mapping SCT FoR (%s) -> VCBCT FoR (%s). ', ...
         'Check that the correct SRO is present in %s.'], ...
        sct_for_uid, vcbct_for_uid, sct_dir);
end

%% ======================== BUILD T_SCT_TO_ICBCT ========================

if strcmp(sct_for_uid, icbct_for_uid)
    fprintf('  [STEP 0.55] SCT and ICBCT share FrameOfReferenceUID - using identity for ICBCT swap.\n');
    T_sct_to_icbct = eye(4);
else
    fprintf('  [STEP 0.55] SCT and ICBCT FoR differ - deriving T_sct_to_icbct from second SRO.\n');
    T_icbct_to_vcbct = findTransform(reg_entries, icbct_for_uid, vcbct_for_uid);
    if isempty(T_icbct_to_vcbct)
        error('step055_swap_imageset:MissingICBCTtoVCBCT', ...
            ['SCT and ICBCT have different FoRs and no ICBCT->VCBCT registration ', ...
             'is available; a second SRO is required to derive T_sct_to_icbct.']);
    end
    T_sct_to_icbct = T_icbct_to_vcbct \ T_sct_to_vcbct;   % inv(A) * B
end

%% ======================== GENERATE BOTH OUTPUTS ========================

icbct_out = fullfile(output_dir, sprintf('RTPLAN_%s_%s_reference_on_ICBCT.dcm', ...
    patient_id, session));
vcbct_out = fullfile(output_dir, sprintf('RTPLAN_%s_%s_reference_on_VCBCT.dcm', ...
    patient_id, session));

ok_icbct = writeSwappedPlan(src_plan_meta, T_sct_to_icbct, ...
    icbct_ct_meta, icbct_rs_path, 'ICBCT', icbct_out, ...
    ROT_TOL_DEG, COORD_WARN_MM);

ok_vcbct = writeSwappedPlan(src_plan_meta, T_sct_to_vcbct, ...
    vcbct_ct_meta, vcbct_rs_path, 'VCBCT', vcbct_out, ...
    ROT_TOL_DEG, COORD_WARN_MM);

%% ======================== POPULATE OUTPUT STRUCT ========================

if ok_icbct
    swapped.icbct = icbct_out;
    swapped.all{end+1} = icbct_out;
end
if ok_vcbct
    swapped.vcbct = vcbct_out;
    swapped.all{end+1} = vcbct_out;
end

%% ======================== PER-OUTPUT SUMMARY ========================

fprintf('\n  [STEP 0.55] ============= Summary =============\n');
fprintf('  [STEP 0.55] SCT source plan (step 0.6 must IGNORE):\n');
fprintf('              %s\n', swapped.sct_source);
fprintf('  [STEP 0.55] ICBCT output: %s\n', ternaryStr(ok_icbct, 'SUCCESS', 'FAILED'));
if ok_icbct
    fprintf('              %s\n', swapped.icbct);
end
fprintf('  [STEP 0.55] VCBCT output: %s\n', ternaryStr(ok_vcbct, 'SUCCESS', 'FAILED'));
if ok_vcbct
    fprintf('              %s\n', swapped.vcbct);
end
fprintf('  [STEP 0.55] ====================================\n\n');

if ~ok_icbct && ~ok_vcbct
    error('step055_swap_imageset:AllOutputsFailed', ...
        'Both ICBCT and VCBCT swapped plans failed to write; see warnings above.');
end

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function pth = locateFirstSlice(directory, glob_pattern, label)
%LOCATEFIRSTSLICE Return the first matching file (alphabetically) or error.

    listing = dir(fullfile(directory, glob_pattern));
    if isempty(listing)
        error('step055_swap_imageset:MissingTargetCT', ...
            '%s slices not found (pattern "%s") in %s', ...
            label, glob_pattern, directory);
    end
    names = sort({listing.name});
    pth   = fullfile(directory, names{1});
end


function v = getFieldOrEmpty(s, fname)
%GETFIELDOREMPTY Safe struct field accessor.

    if isstruct(s) && isfield(s, fname)
        v = s.(fname);
    else
        v = '';
    end
end


function entries = parseRegistration(reg_path)
%PARSEREGISTRATION Extract (moving_for, fixed_for, 4x4) tuples from an SRO.
%
%   A spatial registration object contains a RegistrationSequence with one
%   item per registered frame-of-reference.  Per the DICOM standard, one
%   item carries the identity matrix (the "fixed" frame); the other carries
%   the actual transform from that item's FrameOfReferenceUID (the moving
%   frame) into the fixed frame.  This function returns one entry per
%   non-identity item, paired with the fixed item's FrameOfReferenceUID.

    entries = struct('moving', {}, 'fixed', {}, 'matrix', {}, 'source_file', {});

    meta = dicominfo(reg_path);

    if ~isfield(meta, 'RegistrationSequence')
        warning('parseRegistration:NoRegSeq', ...
            'File has no RegistrationSequence: %s', reg_path);
        return;
    end

    reg_seq = meta.RegistrationSequence;
    items   = fieldnames(reg_seq);

    % Collect FoR UID and 4x4 for each item.
    n = numel(items);
    item_for = cell(n, 1);
    item_mat = cell(n, 1);
    item_is_identity = false(n, 1);

    for i = 1:n
        item = reg_seq.(items{i});
        if ~isfield(item, 'FrameOfReferenceUID')
            warning('parseRegistration:NoItemFoR', ...
                'RegistrationSequence.%s missing FrameOfReferenceUID', items{i});
            continue;
        end
        item_for{i} = item.FrameOfReferenceUID;

        T = extractMatrixFromItem(item);
        if isempty(T)
            warning('parseRegistration:NoMatrix', ...
                'RegistrationSequence.%s missing MatrixRegistrationSequence; skipping (deformable?)', ...
                items{i});
            continue;
        end
        item_mat{i} = T;
        item_is_identity(i) = isApproxIdentity(T);
    end

    % Identify the fixed FoR (the item with identity matrix).
    fixed_idx = find(item_is_identity, 1, 'first');
    if isempty(fixed_idx)
        warning('parseRegistration:NoFixedFrame', ...
            'No identity-matrix item found in %s; using item 1 as fixed.', reg_path);
        fixed_idx = 1;
    end
    fixed_for = item_for{fixed_idx};

    % Each non-identity item with a valid matrix is a (moving -> fixed) entry.
    for i = 1:n
        if i == fixed_idx, continue; end
        if isempty(item_mat{i}) || isempty(item_for{i}), continue; end

        T = item_mat{i};
        if abs(T(4,1)) > 1e-9 || abs(T(4,2)) > 1e-9 || ...
           abs(T(4,3)) > 1e-9 || abs(T(4,4) - 1) > 1e-9
            warning('parseRegistration:BadBottomRow', ...
                'Bottom row of 4x4 in %s item %s is not [0 0 0 1]; got [%g %g %g %g]', ...
                reg_path, items{i}, T(4,1), T(4,2), T(4,3), T(4,4));
        end

        new_entry = struct('moving', item_for{i}, 'fixed', fixed_for, ...
            'matrix', T, 'source_file', '');
        entries(end+1) = new_entry; %#ok<AGROW>
    end
end


function T = extractMatrixFromItem(item)
%EXTRACTMATRIXFROMITEM Pull the 4x4 FrameOfReferenceTransformationMatrix.
%
%   The 16 doubles in FrameOfReferenceTransformationMatrix are stored
%   row-major per DICOM.  MATLAB's reshape is column-major, so we reshape
%   into (4,4) and transpose.

    T = [];
    if ~isfield(item, 'MatrixRegistrationSequence'), return; end
    mrs = item.MatrixRegistrationSequence;
    if ~isfield(mrs, 'Item_1'), return; end
    mrs1 = mrs.Item_1;
    if ~isfield(mrs1, 'MatrixSequence'), return; end
    ms = mrs1.MatrixSequence;
    if ~isfield(ms, 'Item_1'), return; end
    ms1 = ms.Item_1;
    if ~isfield(ms1, 'FrameOfReferenceTransformationMatrix'), return; end

    v = ms1.FrameOfReferenceTransformationMatrix;
    v = double(v(:));
    if numel(v) ~= 16
        warning('extractMatrixFromItem:BadLength', ...
            'FrameOfReferenceTransformationMatrix has %d elements (expected 16)', numel(v));
        return;
    end
    T = reshape(v, 4, 4).';   % row-major -> 4x4
end


function tf = isApproxIdentity(T)
%ISAPPROXIDENTITY Test 4x4 matrix against identity within tight tolerance.

    tf = all(abs(T(:) - reshape(eye(4), [], 1)) < 1e-9);
end


function T = findTransform(reg_entries, moving_for, fixed_for)
%FINDTRANSFORM Return the 4x4 matching (moving_for -> fixed_for), or [].

    T = [];
    for i = 1:numel(reg_entries)
        if strcmp(reg_entries(i).moving, moving_for) && ...
           strcmp(reg_entries(i).fixed,  fixed_for)
            T = reg_entries(i).matrix;
            return;
        end
    end
end


function ok = writeSwappedPlan(src_meta, T, target_ct_meta, target_rs_path, ...
        target_label, out_path, rot_tol_deg, coord_warn_mm)
%WRITESWAPPEDPLAN Apply transform, rewrite references, and dicomwrite.

    ok = false;

    fprintf('\n  [STEP 0.55] --- Building %s-rooted plan ---\n', target_label);

    % Log the 4x4 (4 decimal places)
    fprintf('  [STEP 0.55] Applied 4x4 (row-major):\n');
    for r = 1:4
        fprintf('              [ %10.4f %10.4f %10.4f %10.4f ]\n', ...
            T(r,1), T(r,2), T(r,3), T(r,4));
    end

    % Rotational-deviation warning
    checkRotation(T, rot_tol_deg, target_label);

    src_for_uid    = getFieldOrEmpty(src_meta, 'FrameOfReferenceUID');
    target_for_uid = getFieldOrEmpty(target_ct_meta, 'FrameOfReferenceUID');
    fprintf('  [STEP 0.55] Source FoR: %s\n', src_for_uid);
    fprintf('  [STEP 0.55] Target FoR: %s\n', target_for_uid);

    % Transform all spatial coordinates
    [out_meta, n_xformed, first_beam_log] = transformAllCoords( ...
        src_meta, T, coord_warn_mm);

    if ~isempty(first_beam_log)
        fprintf('  [STEP 0.55] Beam 1 isocenter: [%.3f %.3f %.3f]  -> [%.3f %.3f %.3f]  delta=[%.3f %.3f %.3f]\n', ...
            first_beam_log.orig(1), first_beam_log.orig(2), first_beam_log.orig(3), ...
            first_beam_log.new(1),  first_beam_log.new(2),  first_beam_log.new(3), ...
            first_beam_log.delta(1), first_beam_log.delta(2), first_beam_log.delta(3));
    end
    fprintf('  [STEP 0.55] Total coordinates transformed: %d\n', n_xformed);

    % Rewrite UIDs / references to root on target CT + RTSTRUCT
    try
        out_meta = rewriteReferences(out_meta, target_ct_meta, target_rs_path, target_label);
    catch ME
        warning('step055_swap_imageset:RewriteFailed', ...
            '[%s] Reference rewrite failed: %s', target_label, ME.message);
        return;
    end

    % Write
    try
        dicomwrite([], out_path, out_meta, 'CreateMode', 'Copy', 'WritePrivate', true);
        fprintf('  [STEP 0.55] Wrote: %s\n', out_path);
        ok = true;
    catch ME
        warning('step055_swap_imageset:WriteFailed', ...
            '[%s] dicomwrite failed: %s', target_label, ME.message);
    end
end


function [meta, n_xformed, first_beam_log] = transformAllCoords(meta, T, coord_warn_mm)
%TRANSFORMALLCOORDS Apply T to known RTPLAN spatial fields.
%
%   Transforms:
%     - BeamSequence.Item_b.ControlPointSequence.Item_c.IsocenterPosition
%     - BeamSequence.Item_b.ControlPointSequence.Item_c.SurfaceEntryPoint
%     - FractionGroupSequence.Item_g.ReferencedBeamSequence.Item_r
%           .BeamDoseSpecificationPoint
%
%   Leaves untouched (distances / leaf-frame / angles):
%     SourceToSurfaceDistance, LeafJawPositions, GantryAngle,
%     BeamLimitingDeviceAngle, PatientSupportAngle.
%
%   Recursively warns on other *Position / *Point 3-vec fields found inside
%   ControlPointSequence that are NOT in the known-transform list.

    n_xformed       = 0;
    first_beam_log  = [];

    KNOWN_CP_FIELDS = {'IsocenterPosition', 'SurfaceEntryPoint'};
    SKIP_LIKE_POS   = {'TableTopEccentricRotationDirection'};  %#ok<NASGU> % placeholder

    if isfield(meta, 'BeamSequence')
        beam_fields = fieldnames(meta.BeamSequence);
        for bi = 1:numel(beam_fields)
            bf   = beam_fields{bi};
            beam = meta.BeamSequence.(bf);
            if ~isfield(beam, 'ControlPointSequence'), continue; end

            cp_fields = fieldnames(beam.ControlPointSequence);
            for ci = 1:numel(cp_fields)
                cpf = cp_fields{ci};
                cp  = beam.ControlPointSequence.(cpf);

                % Known transformable fields
                for ki = 1:numel(KNOWN_CP_FIELDS)
                    fname = KNOWN_CP_FIELDS{ki};
                    if isfield(cp, fname) && isnumeric(cp.(fname)) && numel(cp.(fname)) == 3
                        p_old = double(cp.(fname)(:)).';
                        p_new = applyAffine(T, p_old);
                        checkCoordMagnitude(p_new, coord_warn_mm, fname);
                        cp.(fname) = reshape(p_new, size(cp.(fname)));
                        n_xformed  = n_xformed + 1;

                        if bi == 1 && ci == 1 && strcmp(fname, 'IsocenterPosition')
                            first_beam_log = struct( ...
                                'orig',  p_old, ...
                                'new',   p_new, ...
                                'delta', p_new - p_old);
                        end
                    end
                end

                % Recursive scan for any other *Position/*Point 3-vec - warn only
                warnUnknownPositionFields(cp, KNOWN_CP_FIELDS, ...
                    sprintf('BeamSequence.%s.ControlPointSequence.%s', bf, cpf));

                beam.ControlPointSequence.(cpf) = cp;
            end

            meta.BeamSequence.(bf) = beam;
        end
    end

    if isfield(meta, 'FractionGroupSequence')
        fg_fields = fieldnames(meta.FractionGroupSequence);
        for fi = 1:numel(fg_fields)
            fgf = fg_fields{fi};
            fg  = meta.FractionGroupSequence.(fgf);
            if ~isfield(fg, 'ReferencedBeamSequence'), continue; end

            rb_fields = fieldnames(fg.ReferencedBeamSequence);
            for ri = 1:numel(rb_fields)
                rbf = rb_fields{ri};
                rb  = fg.ReferencedBeamSequence.(rbf);
                if isfield(rb, 'BeamDoseSpecificationPoint') && ...
                   isnumeric(rb.BeamDoseSpecificationPoint) && ...
                   numel(rb.BeamDoseSpecificationPoint) == 3
                    p_old = double(rb.BeamDoseSpecificationPoint(:)).';
                    p_new = applyAffine(T, p_old);
                    checkCoordMagnitude(p_new, coord_warn_mm, 'BeamDoseSpecificationPoint');
                    rb.BeamDoseSpecificationPoint = reshape(p_new, ...
                        size(rb.BeamDoseSpecificationPoint));
                    n_xformed = n_xformed + 1;
                    fg.ReferencedBeamSequence.(rbf) = rb;
                end
            end
            meta.FractionGroupSequence.(fgf) = fg;
        end
    end
end


function p_new = applyAffine(T, p_old)
%APPLYAFFINE Apply 4x4 to a 3-vec (row).

    h = [p_old(:); 1];
    h = T * h;
    p_new = h(1:3).';
end


function checkCoordMagnitude(p, threshold, label)
    if any(abs(p) > threshold)
        warning('step055_swap_imageset:CoordOutOfRange', ...
            'Transformed %s = [%g %g %g] exceeds %g mm (unit bug?)', ...
            label, p(1), p(2), p(3), threshold);
    end
end


function warnUnknownPositionFields(s, known_list, prefix)
%WARNUNKNOWNPOSITIONFIELDS Shallow scan for *Position / *Point fields
%   (numeric, length 3) not in known_list.

    if ~isstruct(s), return; end
    fns = fieldnames(s);
    for i = 1:numel(fns)
        fn = fns{i};
        if any(strcmp(fn, known_list)), continue; end
        if ~(endsWithCI(fn, 'Position') || endsWithCI(fn, 'Point')), continue; end
        v = s.(fn);
        if isnumeric(v) && numel(v) == 3
            warning('step055_swap_imageset:UnknownPositionField', ...
                'Found 3-vec field %s.%s not in known-transform list (NOT transformed)', ...
                prefix, fn);
        end
    end
end


function tf = endsWithCI(str, suffix)
    str    = char(str);
    suffix = char(suffix);
    if numel(suffix) > numel(str)
        tf = false; return;
    end
    tf = strcmpi(str(end-numel(suffix)+1:end), suffix);
end


function checkRotation(T, tol_deg, label)
%CHECKROTATION Warn if rotational part deviates from identity > tol on any axis.

    R = T(1:3, 1:3);
    % Extract Euler angles (ZYX intrinsic); guard against gimbal singularity.
    sy = sqrt(R(1,1)^2 + R(2,1)^2);
    if sy < 1e-6
        rx = atan2(-R(2,3), R(2,2));
        ry = atan2(-R(3,1), sy);
        rz = 0;
    else
        rx = atan2(R(3,2), R(3,3));
        ry = atan2(-R(3,1), sy);
        rz = atan2(R(2,1), R(1,1));
    end
    eul_deg = [rx, ry, rz] * 180 / pi;

    if any(abs(eul_deg) > tol_deg)
        warning('step055_swap_imageset:LargeRotation', ...
            ['[%s] Rotation part deviates from identity: Euler ZYX = [%.3f %.3f %.3f] deg ', ...
             '(threshold %.2f). Gantry/collimator angles are NOT rotated - reviewer must confirm.'], ...
            label, eul_deg(3), eul_deg(2), eul_deg(1), tol_deg);
    end
end


function meta = rewriteReferences(meta, target_ct_meta, target_rs_path, target_label)
%REWRITEREFERENCES Update UIDs and reference chain to root the plan on the
%   target CT + RTSTRUCT.  PatientID / PatientName / StudyID untouched.

    % --- Target CT-derived UIDs ---
    target_study_uid = getFieldOrEmpty(target_ct_meta, 'StudyInstanceUID');
    target_for_uid   = getFieldOrEmpty(target_ct_meta, 'FrameOfReferenceUID');
    target_series_uid = getFieldOrEmpty(target_ct_meta, 'SeriesInstanceUID');

    if isempty(target_study_uid) || isempty(target_for_uid) || isempty(target_series_uid)
        error('step055_swap_imageset:TargetCTMissingUIDs', ...
            'Target CT slice is missing StudyInstanceUID / FrameOfReferenceUID / SeriesInstanceUID');
    end

    % --- Target RTSTRUCT SOPInstanceUID ---
    target_rs_meta = dicominfo(target_rs_path);
    target_rs_sop  = getFieldOrEmpty(target_rs_meta, 'SOPInstanceUID');
    if isempty(target_rs_sop)
        error('step055_swap_imageset:TargetRSMissingSOP', ...
            'Target RTSTRUCT %s has no SOPInstanceUID', target_rs_path);
    end

    % --- Fresh UIDs for the new plan ---
    new_sop    = dicomuid;
    new_series = dicomuid;

    meta.SOPInstanceUID             = new_sop;
    meta.MediaStorageSOPInstanceUID = new_sop;
    meta.SeriesInstanceUID          = new_series;

    % --- Copy from target CT ---
    meta.StudyInstanceUID    = target_study_uid;
    meta.FrameOfReferenceUID = target_for_uid;

    % ReferencedFrameOfReferenceSequence.Item_1.FrameOfReferenceUID
    %   and nested RTReferencedStudySequence -> RTReferencedSeriesSequence.SeriesInstanceUID
    if isfield(meta, 'ReferencedFrameOfReferenceSequence') && ...
       isfield(meta.ReferencedFrameOfReferenceSequence, 'Item_1')
        item1 = meta.ReferencedFrameOfReferenceSequence.Item_1;
        item1.FrameOfReferenceUID = target_for_uid;

        if isfield(item1, 'RTReferencedStudySequence') && ...
           isfield(item1.RTReferencedStudySequence, 'Item_1')
            studyItem = item1.RTReferencedStudySequence.Item_1;
            if isfield(studyItem, 'RTReferencedSeriesSequence') && ...
               isfield(studyItem.RTReferencedSeriesSequence, 'Item_1')
                seriesItem = studyItem.RTReferencedSeriesSequence.Item_1;
                seriesItem.SeriesInstanceUID = target_series_uid;
                studyItem.RTReferencedSeriesSequence.Item_1 = seriesItem;
            else
                warning('step055_swap_imageset:NoRTRefSeriesSeq', ...
                    '[%s] Plan has no RTReferencedSeriesSequence.Item_1 to update', target_label);
            end
            item1.RTReferencedStudySequence.Item_1 = studyItem;
        else
            warning('step055_swap_imageset:NoRTRefStudySeq', ...
                '[%s] Plan has no RTReferencedStudySequence.Item_1 to update', target_label);
        end

        meta.ReferencedFrameOfReferenceSequence.Item_1 = item1;
    else
        warning('step055_swap_imageset:NoRefFoRSeq', ...
            '[%s] Plan has no ReferencedFrameOfReferenceSequence.Item_1', target_label);
    end

    % ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID
    if isfield(meta, 'ReferencedStructureSetSequence') && ...
       isfield(meta.ReferencedStructureSetSequence, 'Item_1')
        meta.ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID = target_rs_sop;
    else
        warning('step055_swap_imageset:NoRefSSSeq', ...
            '[%s] Plan has no ReferencedStructureSetSequence.Item_1', target_label);
    end

    % --- Label / description / timestamps ---
    label_suffix = sprintf('_on_%s', target_label);
    if isfield(meta, 'RTPlanLabel')
        new_label = [meta.RTPlanLabel, label_suffix];
    else
        new_label = ['plan', label_suffix];
    end
    if numel(new_label) > 16   % DICOM SH max 16 chars
        new_label = new_label(1:16);
    end
    meta.RTPlanLabel = new_label;

    if isfield(meta, 'RTPlanDescription')
        meta.RTPlanDescription = [meta.RTPlanDescription, ' - swapped onto ', target_label];
    else
        meta.RTPlanDescription = ['Swapped onto ', target_label];
    end

    now_dt = datetime('now');
    meta.InstanceCreationDate = datestr(now_dt, 'yyyymmdd');
    meta.InstanceCreationTime = datestr(now_dt, 'HHMMSS');
end


function s = ternaryStr(cond, t_str, f_str)
    if cond, s = t_str; else, s = f_str; end
end
