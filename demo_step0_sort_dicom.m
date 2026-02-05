%% DEMO_STEP0_SORT_DICOM.m
%  Demonstration script showing how to use step0_sort_dicom
%  
%  This script can be run standalone or the function can be called
%  from the master pipeline.
%
%  Author: ETHOS Pipeline Team
%  Date: February 2026

clear; clc;

fprintf('=========================================================\n');
fprintf('  Step 0: Sort DICOM Files - Demo Script\n');
fprintf('=========================================================\n\n');

%% ======================== CONFIGURATION ========================
% Configuration struct matching master pipeline format

config = struct();

% Directory paths
config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
config.treatment_site = 'Pancreas';

% Patient and session lists (for batch processing)
config.patients = {'1194203'};
config.sessions = {'Session_1', 'Session_2'};

% Pipeline control flags (for master pipeline integration)
config.run_step0 = true;

%% ======================== STANDALONE USAGE ========================
% Example: Process a single patient/session

fprintf('[Standalone Usage Example]\n');
fprintf('----------------------------\n');

patient_id = '1194203';
session = 'Session_1';

fprintf('Processing single patient/session: %s / %s\n\n', patient_id, session);

try
    sct_dir = step0_sort_dicom(patient_id, session, config);
    
    if ~isempty(sct_dir)
        fprintf('\nSUCCESS: Files sorted to:\n  %s\n', sct_dir);
        
        % List files in output directory
        fprintf('\nFiles in sct directory:\n');
        files = dir(fullfile(sct_dir, '*.dcm'));
        for i = 1:length(files)
            fprintf('  %s\n', files(i).name);
        end
    else
        fprintf('\nWARNING: No sct_dir returned (check warnings above)\n');
    end
    
catch ME
    fprintf('\nERROR: %s\n', ME.message);
    fprintf('  In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end

%% ======================== BATCH PROCESSING ========================
% Example: Process multiple patients/sessions (as in master pipeline)

fprintf('\n\n[Batch Processing Example]\n');
fprintf('----------------------------\n');

% Results storage
results = struct();
results.timestamp = datetime('now');
results.patients = struct();

for p_idx = 1:length(config.patients)
    patient_id = config.patients{p_idx};
    
    for s_idx = 1:length(config.sessions)
        session = config.sessions{s_idx};
        
        fprintf('\n=== Patient: %s, Session: %s ===\n', patient_id, session);
        
        % Create result key
        result_key = sprintf('P%s_%s', patient_id, strrep(session, '_', ''));
        results.patients.(result_key) = struct();
        results.patients.(result_key).patient_id = patient_id;
        results.patients.(result_key).session = session;
        
        if config.run_step0
            try
                sct_dir = step0_sort_dicom(patient_id, session, config);
                
                results.patients.(result_key).sct_dir = sct_dir;
                results.patients.(result_key).status = 'success';
                
                if ~isempty(sct_dir)
                    % Count files
                    dcm_files = dir(fullfile(sct_dir, '*.dcm'));
                    results.patients.(result_key).num_files = length(dcm_files);
                end
                
            catch ME
                results.patients.(result_key).status = 'error';
                results.patients.(result_key).error = ME.message;
                fprintf('  ERROR: %s\n', ME.message);
            end
        else
            fprintf('  Step 0 skipped (config.run_step0 = false)\n');
            results.patients.(result_key).status = 'skipped';
        end
    end
end

%% ======================== SUMMARY ========================

fprintf('\n=========================================================\n');
fprintf('  Processing Summary\n');
fprintf('=========================================================\n');

patient_fields = fieldnames(results.patients);
for i = 1:length(patient_fields)
    p = results.patients.(patient_fields{i});
    fprintf('\n%s / %s:\n', p.patient_id, p.session);
    fprintf('  Status: %s\n', p.status);
    
    if strcmp(p.status, 'success') && isfield(p, 'sct_dir') && ~isempty(p.sct_dir)
        fprintf('  Output: %s\n', p.sct_dir);
        if isfield(p, 'num_files')
            fprintf('  Files: %d\n', p.num_files);
        end
    elseif strcmp(p.status, 'error')
        fprintf('  Error: %s\n', p.error);
    end
end

fprintf('\n=========================================================\n');
fprintf('  Demo Complete\n');
fprintf('=========================================================\n');
