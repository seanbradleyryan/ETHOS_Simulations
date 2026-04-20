%% =========================================================================
%  PIPELINE_COMPRESS.m
%  ETHOS Photoacoustic Pipeline — Dose Processing and Upload Preparation
%  =========================================================================
%
%  PURPOSE:
%  Runs on the Windows work laptop after RayStation field dose export.
%  Processes raw field dose DICOMs (Step 1.5) and prepares the processed
%  .mat files for upload to the Linux cluster (prepare_uploads).
%
%  Run this BEFORE pipeline_simulate.m on the cluster.
%
%  STEPS EXECUTED:
%    Step 1.5 — Process field doses and resample SCT to dose grid
%    Upload   — prepare_uploads: package processed .mat files for transfer
%
%  PREREQUISITES:
%    - MATLAB R2022a or later
%    - Image Processing Toolbox
%    - pipeline_setup.m must have been run successfully on the cluster
%    - RayStation field dose DICOMs exported to
%      F:\RayStationFiles\[PatientID]\[Session]\
%
%  PLATFORM: Windows work laptop only.
%            pipeline_setup.m and pipeline_simulate.m run on the Linux cluster.
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

% --- Directory Paths (Windows work laptop) ---
CONFIG.working_dir  = 'C:/Users/80030361/Documents/ETHOS_Simulations';

% --- Dose Masking (Step 1.5) ---
% Set false to disable zeroing outside body / in couch (for debugging only)
CONFIG.apply_dose_masking = true;

% --- Pipeline Control Flags ---
CONFIG.run_step15       = true;   % Step 1.5: Process doses and resample CT
CONFIG.run_prepare      = true;   % Upload  : Prepare .mat files for cluster

%% ========================= INITIALIZATION ================================

fprintf('=========================================================\n');
fprintf('  ETHOS Pipeline — Compress & Prepare (Steps 1.5 / Upload)\n');
fprintf('=========================================================\n');
fprintf('  Started: %s\n', datetime('now'));
fprintf('  Working directory: %s\n', CONFIG.working_dir);
fprintf('=========================================================\n\n');

% Add pipeline scripts to path
addpath(genpath(fullfile(CONFIG.working_dir, 'PipelineScripts')));

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

            %% ============================================================
            %  UPLOAD PREPARATION: Package processed files for cluster
            %% ============================================================
            if CONFIG.run_prepare
                fprintf('\n[UPLOAD] Preparing files for cluster transfer...\n');

                prepare_uploads(patient_id, session, CONFIG);

                fprintf('[UPLOAD] Complete.\n');
            else
                fprintf('\n[UPLOAD] Skipped (CONFIG.run_prepare = false).\n');
            end

            RESULTS.patients.(result_key).status = 'complete';
            fprintf('\n=== %s/%s: COMPRESS COMPLETE ===\n', patient_id, session);

        catch ME
            fprintf('\n[ERROR] %s\n', ME.message);
            for k = 1:length(ME.stack)
                fprintf('        In %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            RESULTS.patients.(result_key).status        = 'error';
            RESULTS.patients.(result_key).error.message = ME.message;
            RESULTS.patients.(result_key).error.stack   = ME.stack;
            fprintf('\n=== %s/%s: COMPRESS FAILED ===\n', patient_id, session);
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
    rd_files = dir(fullfile(rs_dir, 'dose_*.dcm'));
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
