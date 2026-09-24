%% =========================================================================
%  PIPELINE_COMPRESS.m
%  ETHOS Photoacoustic Pipeline — Dose Processing (Steps 1.4 / 1.5)
%  =========================================================================
%
%  PURPOSE:
%  Runs on the Linux cluster after RayStation field dose NPZ export.
%  Converts NPZ field doses to .mat (Step 1.4) and processes them into
%  per-field / total / SCT-resampled outputs (Step 1.5). The resulting
%  .mat files in RayStationFiles/<id>/<session>/processed/ are the direct
%  input to pipeline_simulate.m — no upload-packaging step is performed.
%
%  Run this BEFORE pipeline_simulate.m.
%
%  STEPS EXECUTED:
%    Step 1.4 — Convert RayStation NPZ field doses to .mat
%    Step 1.5 — Process field doses and resample SCT to dose grid
%
%  PREREQUISITES:
%    - MATLAB R2022a or later
%    - Image Processing Toolbox
%    - pipeline_setup.m must have been run successfully
%    - RayStation field dose NPZ files placed under
%      <working_dir>/RayStationFiles/[PatientID]/[Session]/
%
%  PLATFORM: Linux cluster. working_dir defaults to the directory this
%            script lives in (fileparts(mfilename('fullpath'))).
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: April 2026
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION =================================

% --- Patient and Session Selection ---
CONFIG.patients        = {'1194203'};
CONFIG.sessions        = {'Session_1'};
CONFIG.treatment_site  = 'Pancreas';

% --- Directory Paths ---
% Default working_dir = the directory this script lives in. Works on any
% host (Linux cluster, Windows laptop) as long as the pipeline scripts and
% the RayStationFiles/EthosExports trees are checked out side-by-side.
CONFIG.working_dir  = fileparts(mfilename('fullpath'));

% --- Dose Masking (Step 1.5) ---
% Set false to disable zeroing outside body / in couch (for debugging only)
CONFIG.apply_dose_masking = true;

% --- Storage Format ---
% Dose arrays are mostly zero after body/couch masking;  2D format
% greatly reduces file size. Set false only for debugging/compatibility.
% Reconstruct 3D on load: reshape(full(dose_Gy), dose_dims)
CONFIG.use_sparse_storage  = true;

% --- Pipeline Control Flags ---
CONFIG.run_step14       = true;   % Step 1.4: Convert NPZ field doses to .mat
CONFIG.run_step15       = true;   % Step 1.5: Process doses and resample CT

% --- Resume / Skip-Completed Behavior ---
% When true, each step inspects its expected outputs under processed/ and
% skips any sub-task that is already complete (per-NPZ conversion in 1.4,
% CBCT resample / mask creation / per-field processing / total accumulation
% in 1.5). Set false to force a full re-run.
CONFIG.skip_completed   = true;

%% ========================= INITIALIZATION ================================

fprintf('=========================================================\n');
fprintf('  ETHOS Pipeline — Step 1.4 / 1.5 Dose Processing\n');
fprintf('=========================================================\n');
fprintf('  Started: %s\n', datetime('now'));
fprintf('  Working directory: %s\n', CONFIG.working_dir);
fprintf('=========================================================\n\n');

% Add pipeline scripts to path (run_*/pipeline_* live at the repo root next to
% this file; step* functions and helpers now live in pipeline/ and utils/).
addpath(CONFIG.working_dir);
addpath(genpath(fullfile(CONFIG.working_dir, 'pipeline')));  % step* functions live here
addpath(genpath(fullfile(CONFIG.working_dir, 'utils')));

% Initialize results structure
RESULTS           = struct();
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
            %  Check RayStation files are present before doing anything
            %% ============================================================
            rs_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', patient_id, session);
            [rs_files_exist, num_dose_files] = check_raystation_files(rs_dir);

            if ~rs_files_exist
                fprintf('[STEP 1] WAITING: No RayStation field dose files found.\n');
                fprintf('         Expected location: %s\n', rs_dir);
                fprintf('         Please complete the RayStation export, then re-run.\n');
                RESULTS.patients.(result_key).status = 'awaiting_raystation';
                continue;
            end
            fprintf('[STEP 1] Found %d dose file(s) in %s\n', num_dose_files, rs_dir);

            %% ============================================================
            %  Report which processed outputs are already on disk
            %% ============================================================
            processed_dir = fullfile(rs_dir, 'processed');
            report_processed_status(processed_dir, rs_dir, CONFIG.skip_completed);

            %% ============================================================
            %  STEP 1.4: Convert RayStation NPZ Field Doses to .mat
            %% ============================================================
            if CONFIG.run_step14
                fprintf('\n[STEP 1.4] Converting NPZ field doses to .mat...\n');
                n_converted = step14_npz_to_mat(patient_id, session, CONFIG);
                RESULTS.patients.(result_key).num_npz_converted = n_converted;
                fprintf('[STEP 1.4] Converted %d NPZ file(s).\n', n_converted);
            else
                fprintf('\n[STEP 1.4] Skipped (CONFIG.run_step14 = false).\n');
            end

            %% ============================================================
            %  STEP 1.5: Process Field Doses and Resample CT
            %% ============================================================
            if CONFIG.run_step15
                fprintf('\n[STEP 1.5] Processing field doses and resampling CT...\n');

                [field_doses, sct_resampled, total_rs_dose, dose_metadata] = ...
                    step15_process_doses(patient_id, session, CONFIG);

                num_valid = sum(~cellfun(@isempty, field_doses));
                RESULTS.patients.(result_key).num_fields        = num_valid;
                RESULTS.patients.(result_key).dose_grid_size    = dose_metadata.dimensions;
                RESULTS.patients.(result_key).total_rs_dose_max = max(total_rs_dose(:));

                fprintf('[STEP 1.5] Complete. %d fields processed.\n', num_valid);
                fprintf('           Dose grid: [%d x %d x %d]\n', dose_metadata.dimensions);
                fprintf('           Max RS dose: %.4f Gy\n', max(total_rs_dose(:)));
            else
                fprintf('\n[STEP 1.5] Skipped (CONFIG.run_step15 = false).\n');
            end

            RESULTS.patients.(result_key).status = 'complete';
            fprintf('\n=== %s/%s: PROCESSING COMPLETE ===\n', patient_id, session);

        catch ME
            fprintf('\n[ERROR] %s\n', ME.message);
            for k = 1:length(ME.stack)
                fprintf('        In %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            RESULTS.patients.(result_key).status        = 'error';
            RESULTS.patients.(result_key).error.message = ME.message;
            RESULTS.patients.(result_key).error.stack   = ME.stack;
            fprintf('\n=== %s/%s: PROCESSING FAILED ===\n', patient_id, session);
            continue;
        end

    end  % session loop
end  % patient loop

%% ========================= FINALIZE ======================================

fprintf('\n=========================================================\n');
fprintf('  Compress Pipeline Complete\n');
fprintf('=========================================================\n\n');

generate_compress_summary(RESULTS);


%% =========================================================================
%  HELPER FUNCTIONS
%% =========================================================================

function [exists, num_files] = check_raystation_files(rs_dir)
    exists    = false;
    num_files = 0;
    if ~exist(rs_dir, 'dir'), return; end
    % Preferred input format on Linux: NPZ from calc_beam_plan_doses.py
    rd_files = dir(fullfile(rs_dir, 'dose_*.npz'));
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'dose_*.dcm'));
    end
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'Plan_Field*_Beam*_B*_S*.dcm'));
    end
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'RD.*.dcm'));
    end
    if isempty(rd_files)
        rd_files = dir(fullfile(rs_dir, 'RD*.dcm'));
    end
    num_files = numel(rd_files);
    exists    = (num_files > 0);
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

function report_processed_status(processed_dir, rs_dir, skip_completed)
%REPORT_PROCESSED_STATUS Print which Step 1.4/1.5 outputs are already on disk.
%
%   Purely informational — actual skipping is performed inside the daughter
%   functions (step14_npz_to_mat / step15_process_doses), which honor
%   config.skip_completed. This summary lets the user see what will be
%   reused before processing begins.

    fprintf('\n[STATUS] Processed output check (skip_completed=%d)\n', skip_completed);

    if ~isfolder(processed_dir)
        fprintf('         No processed/ directory yet — full run expected.\n');
        return;
    end

    cbct1     = isfile(fullfile(processed_dir, 'CBCT1_resampled.mat'));
    cbct3     = isfile(fullfile(processed_dir, 'CBCT3_resampled.mat'));
    masks     = isfile(fullfile(processed_dir, 'tissue_masks.mat'));
    total_rs  = isfile(fullfile(processed_dir, 'total_rs_dose.mat'));
    meta      = isfile(fullfile(processed_dir, 'metadata.mat'));
    ct1_total = isfile(fullfile(processed_dir, 'total_dose_CT_1.mat'));
    ct3_total = isfile(fullfile(processed_dir, 'total_dose_CT_3.mat'));

    field_outs = dir(fullfile(processed_dir, 'dose_*.mat'));
    n_field_out = numel(field_outs);

    % NPZ vs converted .mat counts in rs_dir (Step 1.4 progress)
    npz_files  = dir(fullfile(rs_dir, 'dose_*.npz'));
    mat_inputs = dir(fullfile(rs_dir, 'dose_*.mat'));

    fprintf('         Step 1.4 (NPZ -> .mat):  %d NPZ / %d .mat in %s\n', ...
        numel(npz_files), numel(mat_inputs), rs_dir);
    fprintf('         CBCT1_resampled.mat:     %s\n', present(cbct1));
    fprintf('         CBCT3_resampled.mat:     %s\n', present(cbct3));
    fprintf('         tissue_masks.mat:        %s\n', present(masks));
    fprintf('         total_rs_dose.mat:       %s\n', present(total_rs));
    fprintf('         total_dose_CT_1.mat:     %s\n', present(ct1_total));
    fprintf('         total_dose_CT_3.mat:     %s\n', present(ct3_total));
    fprintf('         metadata.mat:            %s\n', present(meta));
    fprintf('         Per-field dose .mat:     %d in processed/\n', n_field_out);
end

function s = present(tf)
    if tf, s = 'present'; else, s = 'MISSING'; end
end

function generate_compress_summary(results)
    fprintf('\n--- Compress Summary ---\n');
    pf = fieldnames(results.patients);
    for i = 1:numel(pf)
        p = results.patients.(pf{i});
        fprintf('\n%s / %s:\n', p.patient_id, p.session);
        fprintf('  Status: %s\n', p.status);
        if strcmp(p.status, 'complete')
            if isfield(p, 'num_fields')
                fprintf('  Fields processed: %d\n', p.num_fields);
            end
            if isfield(p, 'total_rs_dose_max')
                fprintf('  Max RS dose: %.4f Gy\n', p.total_rs_dose_max);
            end
        elseif strcmp(p.status, 'error')
            fprintf('  Error: %s\n', p.error.message);
        end
    end
end
