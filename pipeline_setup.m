%% =========================================================================
%  PIPELINE_SETUP.m
%  ETHOS Pipeline  Pre-RayStation Steps
%  =========================================================================
%
%  PURPOSE:
%  Runs all steps that must be completed BEFORE exporting field doses from
%  RayStation. 
%
%  STEPS EXECUTED:
%    Step 0    Sort DICOM files (SCT + RT + REG + all CBCTs) and moves
%    them into Raystation_input (non-SCT RTSTRUCTs go to
%    Raystation_input/<pid>/<session>/extra_structs)
%    Step 0.5  Raystation is configured to reject plans without a minimum
%    gap between mlc teeth. This fixes Halcyon dual-layer MLC minimum gaps in RTPLAN
%    Step 0.6  Raystation does not have an object for the dose from a single segment which is needed for IRAI, but it does store the dose from a beam. 
% This explodes each beam's segments into individual 2-CP beams.
%               (one output RTPLAN file per original beam)
%
%  PREREQUISITES:
%    - Raw ETHOS DICOM exports in EthosExports/<PatientID>/Pancreas/<Session>/
%
%  PLATFORM: Windows machine (ideally F: drive for Moffitt users.) 
%
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION =================================

% --- Patient and Session Selection ---
CONFIG.patients        = {'1194203'};
CONFIG.sessions = {'Session_1','Session_2'}; 
%CONFIG.sessions        = {'Session_1','Session_2','Session_3','Session_4','Session_5','Session_6','Session_7','Session_8','Session_9','Session_10'};
CONFIG.treatment_site  = 'Pancreas';

% --- Directory Paths ---
% working_dir = the git repo root (found by utils/get_repo_root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));
CONFIG.working_dir  = get_repo_root();

% --- MLC Gap Correction Parameters (Step 0.5) ---
CONFIG.mlc_min_gap_mm       = 0.5*2;             % Minimum allowed leaf gap (mm)
CONFIG.mlc_expansion_mm     = 0.4*2;             % Symmetric expansion per side (mm)
CONFIG.mlc_position_range   = [-140, 140];     % Valid Halcyon leaf position range (mm)

% --- Pipeline Control Flags ---
CONFIG.run_step0   = true;   % Step 0 : Sort DICOM files
CONFIG.run_step05  = true;   % Step 0.5: Fix MLC gaps
CONFIG.run_step06  = true;   % Step 0.6: Explode beam segments
CONFIG.skip_completed = true; % Skip a patient/session if all setup outputs already exist

%% ========================= INITIALIZATION ================================

fprintf('=========================================================\n');
fprintf('  ETHOS Pipeline  Setup (Steps 0 / 0.5 / 0.6)\n');
fprintf('=========================================================\n');
fprintf('  Started: %s\n', datetime('now'));
fprintf('  Working directory: %s\n', CONFIG.working_dir);
fprintf('=========================================================\n\n');

% Add required paths
addpath(genpath(fullfile(CONFIG.working_dir, 'PipelineScripts')));
addpath(genpath(fullfile(CONFIG.working_dir, 'pipeline')));  % step* functions live here
addpath(genpath(fullfile(CONFIG.working_dir, 'utils')));

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

        if CONFIG.skip_completed && setup_outputs_exist(patient_id, session, CONFIG)
            RESULTS.patients.(result_key).status  = 'skipped';
            RESULTS.patients.(result_key).sct_dir = get_sct_directory(patient_id, session, CONFIG);
            fprintf('  All setup outputs already exist. Skipping.\n');
            continue;
        end

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
            %  STEP 0.6: Explode Beam Segments (Reference only)
            %% ============================================================
            if CONFIG.run_step06
                fprintf('\n[STEP 0.6] Exploding beam segments (reference only)...\n');

                exploded = step06_explode_segments(patient_id, session, CONFIG);

                RESULTS.patients.(result_key).exploded_plans          = exploded;
                RESULTS.patients.(result_key).num_exploded_reference  = numel(exploded.reference);
                RESULTS.patients.(result_key).num_exploded_total      = numel(exploded.all);

                fprintf('[STEP 0.6] Complete. %d segment RTPLAN file(s) written (reference plan).\n', ...
                    numel(exploded.all));

                fprintf('\n--- RayStation Import Instructions ---\n');
                rs_input_dir = fullfile(CONFIG.working_dir, 'Raystation_Input', ...
                    patient_id, session);
                fprintf('  Import REFERENCE plan RTPLAN files from:\n');
                fprintf('    %s\n', rs_input_dir);
                fprintf('  Files follow the convention:\n');
                fprintf('    RTPLAN_{patient}_{session}_reference_B<N>.dcm\n');
                fprintf('  Then:\n');
                fprintf('    1. Recalculate dose for each exploded-segment plan\n');
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
%  HELPER FUNCTIONS
%% =========================================================================

function sct_dir = get_sct_directory(patient_id, session, config)
    sct_dir = fullfile(config.working_dir, 'EthosExports', patient_id, ...
        config.treatment_site, session, 'sct');
end

function done = setup_outputs_exist(patient_id, session, config)
% True when Steps 0 / 0.5 / 0.6 outputs are all on disk:
%   Step 0  : sCT slices in sct/
%   Step 0.5: RTPLAN_reference_adjusted_mlc.dcm in sct/
%   Step 0.6: plan copy in RayStationFiles/, and one exploded RTPLAN per
%             beam of the adjusted plan in Raystation_Input/
    done = false;
    sct_dir      = get_sct_directory(patient_id, session, config);
    adjusted     = fullfile(sct_dir, 'RTPLAN_reference_adjusted_mlc.dcm');
    rs_input_dir = fullfile(config.working_dir, 'Raystation_Input', patient_id, session);
    rs_files_dir = fullfile(config.working_dir, 'RayStationFiles', patient_id, session);

    if isempty(dir(fullfile(sct_dir, 'CT*.dcm'))) || ~isfile(adjusted)
        return;
    end
    if ~isfile(fullfile(rs_files_dir, 'RTPLAN_reference_adjusted_mlc.dcm'))
        return;
    end

    % Compare exploded file count to the beam count of the adjusted plan
    exploded = dir(fullfile(rs_input_dir, ...
        sprintf('RTPLAN_%s_%s_reference_B*.dcm', patient_id, session)));
    rtplan = dicominfo(adjusted, 'UseDictionaryVR', true);
    numBeams = numel(fieldnames(rtplan.BeamSequence));
    done = numel(exploded) >= numBeams;
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
            if isfield(p, 'num_exploded_total')
                fprintf('  Exploded plans written: %d (reference plan)\n', ...
                    p.num_exploded_total);
            end
        elseif strcmp(p.status, 'error')
            fprintf('  Error: %s\n', p.error.message);
        end
    end
end