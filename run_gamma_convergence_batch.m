%% =========================================================================
%  RUN_GAMMA_CONVERGENCE_BATCH.m
%  [+ethos] NOT migrated: already a thin batch that orchestrates
%  study_gamma_index_convergence (itself not migrated -- needs per-iteration
%  recon). See CLAUDE-SIMULATION_CONTEXT.md migration table.
%  Batch driver: run study_gamma_index_convergence over a RANDOM selection of
%  beam/segment field-dose files.
%
%  For each selected CT_1 field dose (with an existing CT_3 counterpart) this
%  calls study_gamma_index_convergence, which produces, per field:
%     - Study A / Study B truth-vs-recon dose + gamma-index panels (iters 1 & 10)
%     - gamma pass-rate vs TR-iteration convergence figures (A, B, overlaid)
%  Figures are saved to BATCH.figure_dir and closed (headless-friendly); a
%  per-field RESULTS .mat is written, plus a batch summary at the end.
%
%  NOTE: HIPAA / remote-execution - WRITTEN here, RUN on the remote device.
%  Do not execute locally.
%% =========================================================================

clear; clc; close all;
batch_timer = tic;

%% ========================= BATCH CONFIGURATION ==========================

BATCH.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
BATCH.patient_id  = '1194203';
BATCH.session     = 'Session_1';

BATCH.num_files   = 10;      % how many beam/segment combinations to sample
BATCH.rng_seed    = 42;      % fixed seed => reproducible random selection ([] = shuffle)

% Only CT_1 reference field doses are candidates (each must have a CT_3 twin).
% Adjust the glob if you want to include 'adapted' plans, other CTs, etc.
BATCH.dose_glob   = 'dose_*_CT_1_B*_*.mat';

% Overrides forwarded to study_gamma_index_convergence for every field. Anything
% left unset falls back to that function's default_config().
BATCH.config = struct();
BATCH.config.working_dir  = BATCH.working_dir;
BATCH.config.patient_id   = BATCH.patient_id;
BATCH.config.session      = BATCH.session;
BATCH.config.plot_results = true;
BATCH.config.save_figures = true;   % save each figure as PNG
BATCH.config.close_figures = true;  % close after saving (no display pile-up)
BATCH.config.save_results  = true;
% Common output location for figures + per-field RESULTS .mat files.
BATCH.config.figure_dir = fullfile(BATCH.working_dir, 'AnalysisResults', ...
    BATCH.patient_id, BATCH.session, 'gamma_convergence');

%% ========================= DISCOVER CANDIDATE FILES =====================

processed_dir = fullfile(BATCH.working_dir, 'RayStationFiles', ...
    BATCH.patient_id, BATCH.session, 'processed');

if ~isfolder(processed_dir)
    error('run_gamma_convergence_batch:NoProcessedDir', ...
        'Processed dir not found: %s', processed_dir);
end

listing = dir(fullfile(processed_dir, BATCH.dose_glob));
if isempty(listing)
    error('run_gamma_convergence_batch:NoDoseFiles', ...
        'No dose files matching "%s" in %s', BATCH.dose_glob, processed_dir);
end

% Keep only CT_1 files whose derived CT_3 counterpart also exists on disk (the
% study needs both to run).
keep = false(numel(listing), 1);
for i = 1:numel(listing)
    ct1_name = listing(i).name;
    ct3_name = strrep(ct1_name, '_CT_1_', '_CT_3_');
    if ~strcmp(ct3_name, ct1_name) && isfile(fullfile(processed_dir, ct3_name))
        keep(i) = true;
    end
end
listing = listing(keep);

nAvail = numel(listing);
if nAvail == 0
    error('run_gamma_convergence_batch:NoPairs', ...
        'No CT_1 dose files with a matching CT_3 counterpart in %s', processed_dir);
end

%% ========================= RANDOM SELECTION =============================

if ~isempty(BATCH.rng_seed)
    rng(BATCH.rng_seed);
else
    rng('shuffle');
end

nPick   = min(BATCH.num_files, nAvail);
sel_idx = randperm(nAvail, nPick);
selected = listing(sel_idx);

fprintf('=========================================================\n');
fprintf('  Gamma Convergence Batch\n');
fprintf('=========================================================\n');
fprintf('  Patient/session:   %s / %s\n', BATCH.patient_id, BATCH.session);
fprintf('  Candidate pairs:   %d (CT_1 with CT_3 twin)\n', nAvail);
fprintf('  Randomly selected: %d  (seed = %s)\n', nPick, mat2str(BATCH.rng_seed));
fprintf('  Output dir:        %s\n', BATCH.config.figure_dir);
fprintf('---------------------------------------------------------\n');
for i = 1:nPick
    fprintf('    %2d.  %s\n', i, selected(i).name);
end
fprintf('=========================================================\n\n');

%% ========================= RUN OVER SELECTION ===========================

summary = struct('name', {}, 'ok', {}, 'pass_A_final', {}, 'pass_B_final', {}, ...
    'pass_N_final', {}, 'runtime_s', {}, 'error', {});

for i = 1:nPick
    name      = selected(i).name;
    dose_path = fullfile(selected(i).folder, name);

    fprintf('\n########## [%d/%d] %s ##########\n', i, nPick, name);
    field_tic = tic;

    summary(i).name = name;
    try
        RESULTS = study_gamma_index_convergence(dose_path, BATCH.config);

        summary(i).ok           = true;
        summary(i).pass_A_final = RESULTS.pass_rate_A(end);
        summary(i).pass_B_final = RESULTS.pass_rate_B(end);
        if isfield(RESULTS, 'pass_rate_N') && ~isempty(RESULTS.pass_rate_N)
            summary(i).pass_N_final = RESULTS.pass_rate_N(end);   % noise floor
        else
            summary(i).pass_N_final = NaN;
        end
        summary(i).error        = '';
    catch ME
        summary(i).ok           = false;
        summary(i).pass_A_final = NaN;
        summary(i).pass_B_final = NaN;
        summary(i).pass_N_final = NaN;
        summary(i).error        = ME.message;
        fprintf('  [ERROR] %s failed: %s\n', name, ME.message);
    end
    summary(i).runtime_s = toc(field_tic);
end

%% ========================= BATCH SUMMARY ================================

batch_runtime = toc(batch_timer);

fprintf('\n==================== BATCH SUMMARY ====================\n');
fprintf('  %-42s  %-6s  %-10s  %-10s  %-12s  %-8s\n', ...
    'dose file', 'ok', 'A final %', 'B final %', 'Noise flr %', 'time (s)');
for i = 1:numel(summary)
    okstr = 'yes';
    if ~summary(i).ok, okstr = 'NO'; end
    fprintf('  %-42s  %-6s  %-10.2f  %-10.2f  %-12.2f  %-8.1f\n', ...
        summary(i).name, okstr, summary(i).pass_A_final, ...
        summary(i).pass_B_final, summary(i).pass_N_final, summary(i).runtime_s);
end
fprintf('------------------------------------------------------\n');
fprintf('  Total: %d files, %d ok, %.1f s (%.2f min).\n', ...
    numel(summary), sum([summary.ok]), batch_runtime, batch_runtime / 60);
fprintf('======================================================\n');

% Persist the batch manifest + summary alongside the figures.
if ~isfolder(BATCH.config.figure_dir)
    mkdir(BATCH.config.figure_dir);
end
summary_file = fullfile(BATCH.config.figure_dir, 'gamma_convergence_batch_summary.mat');
save(summary_file, 'summary', 'BATCH', 'selected', 'batch_runtime', '-v7.3');
fprintf('\nBatch summary saved to: %s\n', summary_file);
