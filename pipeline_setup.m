%% =========================================================================
%  PIPELINE_SETUP.m
%  ETHOS Photoacoustic Pipeline  Pre-RayStation Steps
%  =========================================================================
%
%  PURPOSE:
%  Runs all steps that must be completed BEFORE exporting field doses from
%  RayStation.  After this script completes, the user must:
%    1. Import the per-beam exploded-segment RTPLAN files (from
%       Raystation_Input/<PatientID>/<Session>/) into RayStation
%    2. Recalculate dose for each exploded-segment beam
%    3. Export field doses as Plan_Field*_Beam*_B*_S*.dcm files into
%       RayStationFiles/<PatientID>/<Session>/
%    4. Run pipeline_simulate.m
%
%  STEPS EXECUTED:
%    Step 0    Sort DICOM files (SCT + matched RT files)
%    Step 0.5  Fix Halcyon dual-layer MLC minimum gaps in RTPLAN
%    Step 0.6  Explode each beam's segments into individual 2-CP beams
%               (one output RTPLAN file per original beam)
%
%  PREREQUISITES:
%    - MATLAB R2022a or later
%    - Image Processing Toolbox (dicominfo, dicomwrite, dicomuid)
%    - Raw ETHOS DICOM exports in EthosExports/<PatientID>/Pancreas/<Session>/
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: March 2026
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION =================================

% --- Patient and Session Selection ---
CONFIG.patients        = {'1194203'};
CONFIG.sessions        = {'Session_1'};
CONFIG.treatment_site  = 'Pancreas';

% --- Directory Paths ---
CONFIG.working_dir  = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.matrad_path  = '/mnt/weka/home/80030361/MATLAB/Addons/matRad';

% --- MLC Gap Correction Parameters (Step 0.5) ---
CONFIG.mlc_min_gap_mm       = 0.5;             % Minimum allowed leaf gap (mm)
CONFIG.mlc_expansion_mm     = 0.4;             % Symmetric expansion per side (mm)
CONFIG.mlc_position_range   = [-140, 140];     % Valid Halcyon leaf position range (mm)

% --- Pipeline Control Flags ---
CONFIG.run_step0   = true;   % Step 0 : Sort DICOM files
CONFIG.run_step05  = true;   % Step 0.5: Fix MLC gaps
CONFIG.run_step06  = true;   % Step 0.6: Explode beam segments

%% ========================= INITIALIZATION ================================

fprintf('=========================================================\n');
fprintf('  ETHOS Pipeline  Setup (Steps 0 / 0.5 / 0.6)\n');
fprintf('=========================================================\n');
fprintf('  Started: %s\n', datetime('now'));
fprintf('  Working directory: %s\n', CONFIG.working_dir);
fprintf('=========================================================\n\n');

% Add required paths
addpath(genpath(CONFIG.matrad_path));
addpath(genpath(fullfile(CONFIG.working_dir, 'PipelineScripts')));

% Initialize results structure
RESULTS = struct();
RESULTS.timestamp = datetime('now');
RESULTS.config    = CONFIG;
RESULTS.patients  = struct();

%% ========================= MAIN PROCESSING LOOP ==========================

for p_idx = 1:length(CONFIG.patients)
    patient_id = CONFIG.patients{p_idx};

    for s_idx = 1:length(CONFIG.sessions)
        session = CONFIG.sessions{s_idx};

        fprintf('\n=== Processing: Patient %s, %s ===\n', patient_id, session);

        result_key = make_result_key(patient_id, session);
        RESULTS.patients.(result_key) = init_patient_result(patient_id, session);

        try

            %% ============================================================
            %  STEP 0: Sort DICOM Files
            %% ============================================================
            if CONFIG.run_step0
                fprintf('\n[STEP 0] Sorting DICOM files...\n');

                sct_dir = step0_sort_dicom(patient_id, session, CONFIG);

                RESULTS.patients.(result_key).sct_dir = sct_dir;
                fprintf('[STEP 0] Complete. Output: %s\n', sct_dir);
            else
                sct_dir = get_sct_directory(patient_id, session, CONFIG);
                fprintf('[STEP 0] Skipped. Using existing: %s\n', sct_dir);
            end

            if isempty(sct_dir) || ~isfolder(sct_dir)
                error('SCT directory not found after Step 0: %s', sct_dir);
            end

            %% ============================================================
            %  STEP 0.5: Fix MLC Gaps in RTPLAN
            %% ============================================================
            if CONFIG.run_step05
                fprintf('\n[STEP 0.5] Fixing MLC gaps...\n');

                [adjusted_rtplan_path, num_corrections] = ...
                    step05_fix_mlc_gaps(patient_id, session, CONFIG);

                RESULTS.patients.(result_key).adjusted_rtplan = adjusted_rtplan_path;
                RESULTS.patients.(result_key).mlc_corrections = num_corrections;
                fprintf('[STEP 0.5] Complete. Corrections made: Hopefully a lot?');
            else
                fprintf('\n[STEP 0.5] Skipped.\n');
                % Verify the adjusted plan exists so Step 0.6 can proceed
                adjusted_rtplan_path = find_adjusted_rtplan(sct_dir);
                if isempty(adjusted_rtplan_path)
                    error('No *_adjusted_mlc.dcm found in %s  run Step 0.5 first.', sct_dir);
                end
                fprintf('[STEP 0.5] Using existing adjusted RTPLAN: %s\n', ...
                    adjusted_rtplan_path);
            end

            %% ============================================================
            %  STEP 0.6: Explode Beam Segments
            %% ============================================================
            if CONFIG.run_step06
                fprintf('\n[STEP 0.6] Exploding beam segments...\n');

                output_paths = step06_explode_segments(patient_id, session, CONFIG);

                RESULTS.patients.(result_key).exploded_plans = output_paths;
                RESULTS.patients.(result_key).num_exploded_plans = numel(output_paths);
                fprintf('[STEP 0.6] Complete. %d per-beam RTPLAN file(s) written.\n', ...
                    numel(output_paths));

                fprintf('\n--- RayStation Import Instructions ---\n');
                rs_input_dir = fullfile(CONFIG.working_dir, 'Raystation_Input', ...
                    patient_id, session);
                fprintf('  Import the following files into RayStation:\n');
                fprintf('    %s\n', rs_input_dir);
                fprintf('  Then:\n');
                fprintf('    1. Recalculate dose for each exploded-segment beam\n');
                fprintf('    2. Export field doses as Plan_Field*_Beam*_B*_S*.dcm into:\n');
                fprintf('       %s\n', fullfile(CONFIG.working_dir, ...
                    'RayStationFiles', patient_id, session));
                fprintf('    3. Run pipeline_simulate.m\n');
                fprintf('--------------------------------------\n');
            else
                fprintf('\n[STEP 0.6] Skipped.\n');
            end

            RESULTS.patients.(result_key).status = 'complete';
            fprintf('\n=== %s/%s: SETUP COMPLETE ===\n', patient_id, session);

        catch ME
            fprintf('\n[ERROR] %s\n', ME.message);
            for k = 1:length(ME.stack)
                fprintf('        In %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            RESULTS.patients.(result_key).status        = 'error';
            RESULTS.patients.(result_key).error.message = ME.message;
            RESULTS.patients.(result_key).error.stack   = ME.stack;
            fprintf('\n=== %s/%s: SETUP FAILED ===\n', patient_id, session);
            continue;
        end

    end  % session loop
end  % patient loop

%% ========================= FINALIZE ======================================

fprintf('\n=========================================================\n');
fprintf('  Setup Pipeline Complete\n');
fprintf('=========================================================\n');
generate_setup_summary(RESULTS);
fprintf('\n=========================================================\n\n');


%% =========================================================================
%  STEP 0.6 IMPLEMENTATION
%  (Calls dicominfo / dicomwrite directly  no external function file needed)
%% =========================================================================

function output_paths = step06_explode_segments(patient_id, session, config)
%STEP06_EXPLODE_SEGMENTS Explode each beam's segments into individual 2-CP beams.
%
%  Takes the *_adjusted_mlc.dcm RTPLAN and produces one output RTPLAN per
%  original beam.  Each new beam contains N_seg independent 2-control-point
%  beams (one per segment).  Output files are written to:
%    {working_dir}/Raystation_Input/{patient_id}/{session}/
%  named {patient_id}_{session}_B{N}.dcm.

    output_paths = {};

    sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
        config.treatment_site, session, 'sct');

    output_dir = fullfile(config.working_dir, 'Raystation_Input', patient_id, session);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % --- Find *_adjusted_mlc.dcm ---
    listing = dir(fullfile(sct_dir, '*_adjusted_mlc.dcm'));
    if isempty(listing)
        error('step06_explode_segments:fileNotFound', ...
            'No *_adjusted_mlc.dcm file found in:\n  %s', sct_dir);
    end
    if numel(listing) > 1
        warning('step06_explode_segments:multipleFiles', ...
            'Multiple *_adjusted_mlc.dcm files found; using first: %s', listing(1).name);
    end

    input_rtplan = fullfile(sct_dir, listing(1).name);
    fprintf('  Input:      %s\n', input_rtplan);
    fprintf('  Output dir: %s\n', output_dir);

    % --- Load RTPLAN ---
    fprintf('  Loading RTPLAN...\n');
    rtplan = dicominfo(input_rtplan, 'UseDictionaryVR', true);
    if ~isfield(rtplan, 'Modality') || ~strcmp(rtplan.Modality, 'RTPLAN')
        error('step06_explode_segments:wrongModality', ...
            'File Modality is not RTPLAN: %s', rtplan.Modality);
    end
    if ~isfield(rtplan, 'BeamSequence')
        error('step06_explode_segments:noBeams', 'RTPLAN has no BeamSequence.');
    end

    beam_fields_orig   = fieldnames(rtplan.BeamSequence);
    num_original_beams = numel(beam_fields_orig);
    fprintf('  RTPlanLabel: %s\n', rtplan.RTPlanLabel);
    fprintf('  Original beams: %d\n\n', num_original_beams);

    % --- Build MU map: beam_number -> total_MU ---
    mu_map = containers.Map('KeyType','int32','ValueType','double');
    if ~isfield(rtplan, 'FractionGroupSequence')
        error('step06_explode_segments:noFractionGroup', ...
            'RTPLAN has no FractionGroupSequence.');
    end
    frac_group      = rtplan.FractionGroupSequence.Item_1;
    ref_beam_fields = fieldnames(frac_group.ReferencedBeamSequence);
    for k = 1:numel(ref_beam_fields)
        rb   = frac_group.ReferencedBeamSequence.(ref_beam_fields{k});
        bnum = int32(rb.ReferencedBeamNumber);
        bmu  = double(rb.BeamMeterset);
        mu_map(bnum) = bmu;
    end
    fprintf('  MU map loaded for %d beams.\n\n', mu_map.Count);

    any_warn_global = false;

    for b = 1:num_original_beams

        orig_beam            = rtplan.BeamSequence.(beam_fields_orig{b});
        original_beam_number = int32(orig_beam.BeamNumber);
        original_beam_name   = orig_beam.BeamName;

        if ~isKey(mu_map, original_beam_number)
            error('step06_explode_segments:muNotFound', ...
                'No MU found for beam %d (%s).', original_beam_number, original_beam_name);
        end
        total_mu = mu_map(original_beam_number);

        % Sort control points by their numeric index
        cp_fields = fieldnames(orig_beam.ControlPointSequence);
        cp_nums   = cellfun(@(f) sscanf(f, 'Item_%d'), cp_fields);
        [~, si]   = sort(cp_nums);
        cp_fields = cp_fields(si);
        N_cp      = numel(cp_fields);
        N_seg     = N_cp - 1;

        if N_seg < 1
            warning('step06_explode_segments:noSegments', ...
                'Beam %d (%s) has only %d control points  skipping.', ...
                original_beam_number, original_beam_name, N_cp);
            continue;
        end

        cp1_orig = orig_beam.ControlPointSequence.(cp_fields{1});
        if isfield(cp1_orig, 'GantryAngle')
            gantry_str = sprintf('%.1f deg', double(cp1_orig.GantryAngle));
        else
            gantry_str = 'N/A';
        end
        fprintf('  Beam B%d (%s, gantry %s, %d segments, MU=%.4f)...\n', ...
            original_beam_number, original_beam_name, gantry_str, N_seg, total_mu);

        new_beam_seq     = struct();
        new_ref_beam_seq = struct();
        mu_sum           = 0;

        for s = 1:N_seg
            new_beam = orig_beam;
            new_beam.BeamNumber            = s;
            new_beam.BeamName              = sprintf('B%d_S%03d', original_beam_number, s-1);
            new_beam.NumberOfControlPoints = 2;

            cp_entry_orig = orig_beam.ControlPointSequence.(cp_fields{s});
            cp_exit_orig  = orig_beam.ControlPointSequence.(cp_fields{s+1});

            % CumulativeMetersetWeight values
            if isfield(cp_entry_orig, 'CumulativeMetersetWeight')
                cmw_entry = double(cp_entry_orig.CumulativeMetersetWeight);
            else
                cmw_entry = (s-1) / N_seg;
            end
            if isfield(cp_exit_orig, 'CumulativeMetersetWeight')
                cmw_exit = double(cp_exit_orig.CumulativeMetersetWeight);
            else
                cmw_exit = s / N_seg;
            end
            seg_mu = (cmw_exit - cmw_entry) * total_mu;
            mu_sum = mu_sum + seg_mu;

            % --- Build entry CP (Item_1) ---
            cp_entry_new = cp_entry_orig;
            cp_entry_new.CumulativeMetersetWeight = 0;
            cp_entry_new.ControlPointIndex        = 0;

            % Back-fill any geometry fields from original CP_1 that are absent
            cp1_skip = {'BeamLimitingDevicePositionSequence', ...
                        'CumulativeMetersetWeight', 'ControlPointIndex'};
            cp1_all = fieldnames(cp1_orig);
            for gf = 1:numel(cp1_all)
                fld = cp1_all{gf};
                if ~isfield(cp_entry_new, fld) && ~ismember(fld, cp1_skip)
                    cp_entry_new.(fld) = cp1_orig.(fld);
                end
            end

            % Build 4-item BLDPS for entry CP (jaw + jaw + MLC1 + MLC2)
            bldps_cp1   = cp1_orig.BeamLimitingDevicePositionSequence;
            bldps_entry = cp_entry_orig.BeamLimitingDevicePositionSequence;
            n_bldps     = numel(fieldnames(bldps_entry));
            if n_bldps == 4
                mlc1 = bldps_entry.Item_3;
                mlc2 = bldps_entry.Item_4;
            else
                mlc1 = bldps_entry.Item_1;
                mlc2 = bldps_entry.Item_2;
            end
            bldps_new_entry        = struct();
            bldps_new_entry.Item_1 = bldps_cp1.Item_1;   % X jaw
            bldps_new_entry.Item_2 = bldps_cp1.Item_2;   % Y jaw
            bldps_new_entry.Item_3 = mlc1;               % MLC1
            bldps_new_entry.Item_4 = mlc2;               % MLC2
            cp_entry_new.BeamLimitingDevicePositionSequence = bldps_new_entry;

            % --- Build exit CP (Item_2) ---
            cp_exit_new = cp_exit_orig;
            cp_exit_new.CumulativeMetersetWeight = 1;
            cp_exit_new.ControlPointIndex        = 1;
            if ~isfield(cp_exit_new, 'BeamLimitingDevicePositionSequence')
                bldps_exit_fallback        = struct();
                bldps_exit_fallback.Item_1 = mlc1;
                bldps_exit_fallback.Item_2 = mlc2;
                cp_exit_new.BeamLimitingDevicePositionSequence = bldps_exit_fallback;
            end

            new_beam.ControlPointSequence              = struct();
            new_beam.ControlPointSequence.Item_1       = cp_entry_new;
            new_beam.ControlPointSequence.Item_2       = cp_exit_new;
            new_beam.FinalCumulativeMetersetWeight      = 1;

            new_beam_seq.(sprintf('Item_%d', s)) = new_beam;

            ref_entry                         = struct();
            ref_entry.ReferencedBeamNumber    = s;
            ref_entry.BeamMeterset            = seg_mu;
            new_ref_beam_seq.(sprintf('Item_%d', s)) = ref_entry;
        end  % segment loop

        % MU validation
        delta_pct = 0;
        if total_mu > 0
            delta_pct = 100 * abs(mu_sum - total_mu) / total_mu;
        end
        if delta_pct > 0.5
            warning('step06_explode_segments:muMismatch', ...
                'Beam B%d MU sum deviation = %.4f%% (>0.5%%)', ...
                original_beam_number, delta_pct);
            any_warn_global = true;
            mu_status = '[WARN]';
        else
            mu_status = '[OK]';
        end
        fprintf('    MU: orig=%.4f  sum=%.4f  delta=%.4f%%  %s\n', ...
            total_mu, mu_sum, delta_pct, mu_status);

        % --- Assemble per-beam output plan ---
        rtplan_out = rtplan;
        rtplan_out.BeamSequence = new_beam_seq;
        rtplan_out.FractionGroupSequence.Item_1.ReferencedBeamSequence = new_ref_beam_seq;
        rtplan_out.FractionGroupSequence.Item_1.NumberOfBeams = N_seg;

        new_uid = dicomuid();
        rtplan_out.SOPInstanceUID             = new_uid;
        rtplan_out.MediaStorageSOPInstanceUID = new_uid;

        rtplan_out.RTPlanLabel       = original_beam_name(1:min(16,end));
        if isfield(rtplan_out, 'RTPlanName')
            rtplan_out.RTPlanName    = original_beam_name;
        end
        rtplan_out.RTPlanDescription = sprintf('Beam %s (%s %s) exploded segments', ...
            original_beam_name, patient_id, session);

        beam_output_path = fullfile(output_dir, ...
            sprintf('%s_%s_B%d.dcm', patient_id, session, original_beam_number));

        fprintf('    Writing: %s\n', beam_output_path);
        dicomwrite([], beam_output_path, rtplan_out, 'CreateMode', 'Copy');

        output_paths{end+1} = beam_output_path; %#ok<AGROW>

    end  % beam loop

    if any_warn_global
        warning('step06_explode_segments:globalMuWarn', ...
            'One or more beams had MU sum deviation >0.5%%  verify output.');
    end
    fprintf('  %d per-beam RTPLAN files written to:\n    %s\n', ...
        numel(output_paths), output_dir);
end


%% =========================================================================
%  HELPER FUNCTIONS
%% =========================================================================

function sct_dir = get_sct_directory(patient_id, session, config)
    sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
        config.treatment_site, session, 'sct');
end

function adjusted_path = find_adjusted_rtplan(sct_dir)
    listing = dir(fullfile(sct_dir, '*_adjusted_mlc.dcm'));
    if isempty(listing)
        adjusted_path = '';
    else
        adjusted_path = fullfile(sct_dir, listing(1).name);
    end
end

function key = make_result_key(patient_id, session)
    key = sprintf('P%s_%s', patient_id, strrep(session, '_', ''));
end

function result = init_patient_result(patient_id, session)
    result            = struct();
    result.patient_id = patient_id;
    result.session    = session;
    result.status     = 'started';
    result.start_time = datetime('now');
end

function generate_setup_summary(results)
    fprintf('\n--- Setup Summary ---\n');
    pf = fieldnames(results.patients);
    for i = 1:length(pf)
        p = results.patients.(pf{i});
        fprintf('\n%s / %s:\n', p.patient_id, p.session);
        fprintf('  Status: %s\n', p.status);
        if strcmp(p.status, 'complete')
            if isfield(p, 'sct_dir') && ~isempty(p.sct_dir)
                fprintf('  SCT dir: %s\n', p.sct_dir);
            end
            if isfield(p, 'mlc_corrections')
%                 fprintf('  MLC corrections: %d\n', p.mlc_corrections);
            end
            if isfield(p, 'num_exploded_plans')
                fprintf('  Exploded plans written: %d\n', p.num_exploded_plans);
            end
        elseif strcmp(p.status, 'error')
            fprintf('  Error: %s\n', p.error.message);
        end
    end
end
