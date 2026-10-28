%% =========================================================================
%  MIGRATE_RESULTS_LAYOUT.m
%  One-off: move results written before the storage overhaul to the new layout
%  =========================================================================
%
%  PURPOSE:
%  The pipeline now keeps simulation data beside the processed doses and names
%  analysis folders after the config instead of the hash. This script moves
%  results that were written in the old layout:
%    1. SimulationResults/<id>/<session>/<method>/simulation_debug/<plot>_<hash>.png
%         -> AnalysisResults/<id>/<session>/<label>/verify_pipeline_simulate/
%    2. SimulationResults/<id>/<session>/<method>/*  (recons, totals, masks,
%       manifests, config_registry.json, simulation_log.txt, field_status/,
%       Step 2.5 / 2.05 files)
%         -> RayStationFiles/<id>/<session>/processed/simulation_results/<method>/
%    3. AnalysisResults/<id>/<session>/noise_ensemble_<id>_<session>_<hash>.mat
%         -> the simulation_results/<method>/ folder of that hash
%    4. AnalysisResults/<id>/<session>/<method>/<hash>/  (Step 3 outputs)
%         -> AnalysisResults/<id>/<session>/<label>/step3_analysis/
%    5. AnalysisResults/<id>/<session>/compress_verification/
%         -> AnalysisResults/<id>/<session>/verify_pipeline_compress_output/
%  <label> = config_label of the hash's config in config_registry.json
%  ('canonical' or the fields that differ from get_default_config).
%
%  Left in place: AnalysisResults/pipeline_results_*.mat and timing_report_*.txt
%  (one file covers several sessions), and hashes missing from every registry.
%  Nothing is overwritten: an item whose destination exists is skipped and listed.
%  Emptied old folders are removed.
%
%  INPUTS:  CONFIG.dry_run (true = only print what would move). Run once with
%           true, read the output, then set false and run again.
%  OUTPUTS: Files moved; console report.
%
%  DEPENDENCIES: utils/ (get_repo_root, get_simulation_dir, config_label,
%                compute_sim_config_hash, get_default_config)
%  =========================================================================

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'utils')));

%% ========================= CONFIGURATION =================================
CONFIG.dry_run     = true;              % true = only print; set false to move
CONFIG.working_dir = get_repo_root();

%% ========================= MIGRATION =====================================
workingDir   = CONFIG.working_dir;
oldSimRoot   = fullfile(workingDir, 'SimulationResults');
analysisRoot = fullfile(workingDir, 'AnalysisResults');
if CONFIG.dry_run
    fprintf('DRY RUN: nothing is moved.\n');
end

% Every patient/session that has an old SimulationResults or AnalysisResults folder
sessions = cell(0, 2);   % rows of {patient_id, session}
for root = {oldSimRoot, analysisRoot}
    patients = list_subfolders(root{1});
    for p = 1:numel(patients)
        sessionNames = list_subfolders(fullfile(root{1}, patients{p}));
        for s = 1:numel(sessionNames)
            sessions(end + 1, :) = {patients{p}, sessionNames{s}}; %#ok<SAGROW>
        end
    end
end
[~, keep] = unique(strcat(sessions(:, 1), '/', sessions(:, 2)));
sessions  = sessions(sort(keep), :);

for k = 1:size(sessions, 1)
    pid     = sessions{k, 1};
    session = sessions{k, 2};
    fprintf('\n=== %s / %s ===\n', pid, session);
    oldSessionDir      = fullfile(oldSimRoot, pid, session);
    analysisSessionDir = fullfile(analysisRoot, pid, session);

    % Hash -> gruneisen method and label, from the old and new registries
    regFiles = [dir(fullfile(oldSessionDir, '*', 'config_registry.json')); ...
                dir(fullfile(workingDir, 'RayStationFiles', pid, session, 'processed', ...
                    'simulation_results', '*', 'config_registry.json'))];
    hashes = {};
    methods = {};
    labels = {};
    for r = 1:numel(regFiles)
        [~, method] = fileparts(regFiles(r).folder);
        registry = jsondecode(fileread(fullfile(regFiles(r).folder, regFiles(r).name)));
        keys = fieldnames(registry);
        for i = 1:numel(keys)
            hashes{end + 1}  = keys{i}(3:end);   %#ok<SAGROW> key is 'h_<hash>'
            methods{end + 1} = method;           %#ok<SAGROW>
            labels{end + 1}  = config_label(registry.(keys{i}).config); %#ok<SAGROW>
            fprintf('  hash %s (%s) -> %s\n', hashes{end}, method, labels{end});
        end
    end

    % 1 + 2. Old simulation folders: debug plots to AnalysisResults, the rest to processed/
    oldMethods = list_subfolders(oldSessionDir);
    for m = 1:numel(oldMethods)
        oldMethodDir = fullfile(oldSessionDir, oldMethods{m});
        debugDir     = fullfile(oldMethodDir, 'simulation_debug');
        pngs = dir(fullfile(debugDir, '*.png'));
        for i = 1:numel(pngs)
            tok = regexp(pngs(i).name, '_([0-9a-f]{8})\.png$', 'tokens', 'once');
            label = '';
            if ~isempty(tok)
                label = lookup_label(tok{1}, hashes, labels);
            end
            if isempty(label)
                fprintf('  SKIP (hash not in a registry): %s\n', fullfile(debugDir, pngs(i).name));
                continue;
            end
            move_item(fullfile(debugDir, pngs(i).name), fullfile(analysisSessionDir, label, ...
                'verify_pipeline_simulate', pngs(i).name), CONFIG.dry_run);
        end
        remove_if_empty(debugDir, CONFIG.dry_run);

        newSimDir = get_simulation_dir(workingDir, pid, session, oldMethods{m});
        items = dir(oldMethodDir);
        items = items(~ismember({items.name}, {'.', '..', 'simulation_debug'}));
        fprintf('  %d item(s): %s\n    -> %s\n', numel(items), oldMethodDir, newSimDir);
        for i = 1:numel(items)
            src = fullfile(oldMethodDir, items(i).name);
            dst = fullfile(newSimDir, items(i).name);
            if strcmp(items(i).name, 'config_registry.json') && isfile(dst)
                merge_registries(src, dst, CONFIG.dry_run);   % pipeline already ran in the new layout
            else
                move_item(src, dst, CONFIG.dry_run, true);
            end
        end
        remove_if_empty(oldMethodDir, CONFIG.dry_run);
    end
    remove_if_empty(oldSessionDir, CONFIG.dry_run);
    remove_if_empty(fileparts(oldSessionDir), CONFIG.dry_run);

    % 3. Noise-only ensemble caches -> beside the recons of their hash
    caches = dir(fullfile(analysisSessionDir, 'noise_ensemble_*.mat'));
    for i = 1:numel(caches)
        tok = regexp(caches(i).name, '_([0-9a-f]{8})\.mat$', 'tokens', 'once');
        idx = [];
        if ~isempty(tok)
            idx = find(strcmp(hashes, tok{1}), 1);
        end
        if isempty(idx)
            fprintf('  SKIP (hash not in a registry): %s\n', caches(i).name);
            continue;
        end
        move_item(fullfile(analysisSessionDir, caches(i).name), fullfile(get_simulation_dir( ...
            workingDir, pid, session, methods{idx}), caches(i).name), CONFIG.dry_run);
    end

    % 4. Step 3: <method>/<hash>/ -> <label>/step3_analysis/ (only registry methods,
    %    so study folders such as gamma_sensitivity/<hash>/ are not touched)
    uniqueMethods = unique(methods);
    for m = 1:numel(uniqueMethods)
        methodDir = fullfile(analysisSessionDir, uniqueMethods{m});
        hashDirs  = list_subfolders(methodDir);
        for i = 1:numel(hashDirs)
            label = lookup_label(hashDirs{i}, hashes, labels);
            if isempty(label)
                fprintf('  SKIP (not a known hash): %s\n', fullfile(methodDir, hashDirs{i}));
                continue;
            end
            move_item(fullfile(methodDir, hashDirs{i}), ...
                fullfile(analysisSessionDir, label, 'step3_analysis'), CONFIG.dry_run);
        end
        remove_if_empty(methodDir, CONFIG.dry_run);
    end

    % 5. Compress verification folder rename
    move_item(fullfile(analysisSessionDir, 'compress_verification'), ...
        fullfile(analysisSessionDir, 'verify_pipeline_compress_output'), CONFIG.dry_run);
end
remove_if_empty(oldSimRoot, CONFIG.dry_run);
fprintf('\nDone.\n');


%% =========================================================================
%  HELPERS
%% =========================================================================

function names = list_subfolders(folder)
    % Names of the subfolders of folder ({} when it does not exist).
    names = {};
    if ~isfolder(folder), return; end
    d = dir(folder);
    names = {d([d.isdir] & ~ismember({d.name}, {'.', '..'})).name};
end

function label = lookup_label(hash8, hashes, labels)
    % Label of a hash, '' when no registry lists it.
    idx   = find(strcmp(hashes, hash8), 1);
    label = '';
    if ~isempty(idx)
        label = labels{idx};
    end
end

function move_item(src, dst, dryRun, quiet)
    % Move a file or folder; never overwrite. quiet = only print skips.
    if nargin < 4, quiet = false; end
    if ~exist(src, 'file'), return; end
    if exist(dst, 'file')
        fprintf('  SKIP (destination exists): %s\n', dst);
        return;
    end
    if ~quiet
        fprintf('  %s\n    -> %s\n', src, dst);
    end
    if dryRun, return; end
    if ~isfolder(fileparts(dst)), mkdir(fileparts(dst)); end
    movefile(src, dst);
end

function merge_registries(src, dst, dryRun)
    % Add the old registry's hashes to the new one (an entry already in the new
    % registry is kept), then delete the old file.
    oldReg = jsondecode(fileread(src));
    newReg = jsondecode(fileread(dst));
    keys = fieldnames(oldReg);
    for i = 1:numel(keys)
        if ~isfield(newReg, keys{i})
            newReg.(keys{i}) = oldReg.(keys{i});
        end
    end
    fprintf('  merge %s\n    -> %s\n', src, dst);
    if dryRun, return; end
    fid = fopen(dst, 'w');
    fprintf(fid, '%s', jsonencode(newReg, 'PrettyPrint', true));
    fclose(fid);
    delete(src);
end

function remove_if_empty(folder, dryRun)
    % rmdir without 's' only removes an empty folder; the output argument keeps
    % it from erroring when the folder still has files.
    if dryRun || ~isfolder(folder), return; end
    [~] = rmdir(folder);
end
