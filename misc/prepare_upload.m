%% =========================================================================
%  PREPARE_UPLOAD.m
%  ETHOS Photoacoustic Pipeline — Pre-Upload Data Preparation
%  =========================================================================
%
%  PURPOSE:
%  Dramatically reduces upload size by replacing the full RayStation DICOM
%  export with only the processed .mat files needed by pipeline_simulate.m.
%
%  WHAT THIS DOES:
%    1. Runs step15_process_doses if processed/ outputs don't already exist
%       (or if CONFIG.force_rerun = true)
%    2. Collects all .mat files from the processed/ directory
%    3. Optionally also includes RTPLAN DICOMs (small, sometimes useful)
%    4. Creates a single compressed archive ready to upload
%    5. Reports size reduction and estimated upload time savings
%
%  WHAT YOU UPLOAD vs. WHAT YOU USED TO UPLOAD:
%    Before: Full RayStation export (CT slices + RTDOSE fields + RTSTRUCT)
%            Typically 1–5 GB
%    After:  Processed .mat archive only
%            Typically 50–300 MB
%
%  RECIPIENT USAGE:
%    1. Extract archive into RayStationFiles/<patient_id>/<session>/processed/
%    2. In pipeline_simulate.m, set CONFIG.run_step15 = false
%    3. Run pipeline_simulate.m normally — it loads from processed/ directly
%
%  PREREQUISITES:
%    - step15_process_doses.m on MATLAB path
%    - load_processed_data.m on MATLAB path
%    - RayStation DICOM field dose files present in RayStationFiles/ dir
%      (only needed if step15 has not been run yet)
%    - Image Processing Toolbox (for step15 DICOM reading)
%    - tar must be available on system PATH (Linux/Mac); on Windows uses zip
%
%  AUTHOR: ETHOS Pipeline
%  DATE: April 2026
%  =========================================================================

clear; clc;

% Script lives one level below the repo root; add utils/ and pipeline/ from there.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'utils')));
addpath(genpath(fullfile(repoRoot, 'pipeline')));

%% ========================= CONFIGURATION =================================

% --- Patient / Session ---
CONFIG.patient_id      = '1194203';
CONFIG.session         = 'Session_1';
CONFIG.treatment_site  = 'Pancreas';

% --- Paths ---
CONFIG.working_dir     = '/mnt/weka/home/80030361/ETHOS_Simulations';

% --- Step 15 options (passed through if step15 needs to run) ---
CONFIG.apply_dose_masking = true;

% --- Upload prep options ---

% If processed/ .mat files already exist, skip re-running step15.
% Set to true to force step15 to re-run even if outputs exist.
CONFIG.force_rerun = false;

% Include the RTPLAN DICOM files in the archive (small, ~100-500 KB each).
% Useful if the recipient might need to re-run step15 with different settings.
% RTSTRUCT and CT DICOMs are NOT included — they are large and not needed
% once step15 has run.
CONFIG.include_rtplan_dicoms = true;

% Archive format: 'tar.gz' (recommended, ~30% smaller) or 'zip' (Windows-safe)
% 'tar.gz' requires tar on system PATH. Falls back to 'zip' automatically.
CONFIG.archive_format = 'tar.gz';

% Upload speed estimate in Mbps (used only for time estimate printout).
% Adjust to your actual connection speed.
CONFIG.upload_mbps = 5;

%% ========================= DERIVED PATHS =================================

rs_dir        = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
                    CONFIG.patient_id, CONFIG.session);
processed_dir = fullfile(rs_dir, 'processed');

% Archive will be placed in working_dir for easy access
timestamp_str  = datestr(now, 'yyyymmdd_HHMM');
archive_name   = sprintf('upload_%s_%s_%s', ...
                    CONFIG.patient_id, CONFIG.session, timestamp_str);
archive_dir    = CONFIG.working_dir;

fprintf('==========================================================\n');
fprintf('  PREPARE UPLOAD\n');
fprintf('  Patient:  %s\n', CONFIG.patient_id);
fprintf('  Session:  %s\n', CONFIG.session);
fprintf('  RS dir:   %s\n', rs_dir);
fprintf('  Output:   %s\n', archive_dir);
fprintf('==========================================================\n\n');

%% ========================= VERIFY RS DIRECTORY ===========================

if ~isfolder(rs_dir)
    error('prepare_upload:DirectoryNotFound', ...
        'RayStation directory not found:\n  %s\n', rs_dir);
end

%% ========================= RUN STEP15 IF NEEDED ==========================

% Check whether processed outputs already exist
processed_mat_files = dir(fullfile(processed_dir, '*.mat'));
step15_already_done = ~isempty(processed_mat_files) && ...
    isfile(fullfile(processed_dir, 'metadata.mat')) && ...
    isfile(fullfile(processed_dir, 'sct_resampled.mat')) && ...
    isfile(fullfile(processed_dir, 'total_rs_dose.mat'));

if step15_already_done && ~CONFIG.force_rerun
    fprintf('[Step 1.5] Processed outputs already exist (%d .mat files found).\n', ...
        numel(processed_mat_files));
    fprintf('           Skipping step15. Set CONFIG.force_rerun = true to re-run.\n\n');
else
    if CONFIG.force_rerun && step15_already_done
        fprintf('[Step 1.5] Force re-run requested — overwriting existing outputs.\n\n');
    else
        fprintf('[Step 1.5] No processed outputs found. Running step15_process_doses...\n\n');
    end

    step15_config = struct( ...
        'working_dir',        CONFIG.working_dir, ...
        'treatment_site',     CONFIG.treatment_site, ...
        'apply_dose_masking', CONFIG.apply_dose_masking);

    [~, ~, ~, ~] = step15_process_doses( ...
        CONFIG.patient_id, CONFIG.session, step15_config);

    fprintf('\n[Step 1.5] Complete.\n\n');
end

%% ========================= COLLECT FILES FOR ARCHIVE =====================

fprintf('[Collect] Gathering files to include in archive...\n');

files_to_archive = {};

% --- All .mat files from processed/ ---
mat_files = dir(fullfile(processed_dir, '*.mat'));
if isempty(mat_files)
    error('prepare_upload:NoMatFiles', ...
        'No .mat files found in processed/ directory:\n  %s\n', processed_dir);
end

for i = 1:numel(mat_files)
    files_to_archive{end+1} = fullfile(processed_dir, mat_files(i).name); %#ok<SAGROW>
end

fprintf('  .mat files:  %d file(s) from processed/\n', numel(mat_files));

% --- RTPLAN DICOMs (optional, small) ---
rtplan_count = 0;
if CONFIG.include_rtplan_dicoms
    % Look in both the RS dir and the sct dir
    rtplan_search_dirs = { ...
        rs_dir, ...
        fullfile(CONFIG.working_dir, 'EthosExports', CONFIG.patient_id, ...
            CONFIG.treatment_site, CONFIG.session, 'sct') ...
    };

    for d = 1:numel(rtplan_search_dirs)
        if ~isfolder(rtplan_search_dirs{d}), continue; end
        rp_files = dir(fullfile(rtplan_search_dirs{d}, 'RP*.dcm'));
        rp_files = [rp_files; dir(fullfile(rtplan_search_dirs{d}, 'RTPLAN*.dcm'))]; %#ok<AGROW>
        for i = 1:numel(rp_files)
            fp = fullfile(rp_files(i).folder, rp_files(i).name);
            if ~ismember(fp, files_to_archive)
                files_to_archive{end+1} = fp; %#ok<SAGROW>
                rtplan_count = rtplan_count + 1;
            end
        end
    end
    fprintf('  RTPLAN DCMs: %d file(s)\n', rtplan_count);
end

fprintf('  Total files: %d\n\n', numel(files_to_archive));

%% ========================= MEASURE SIZES =================================

fprintf('[Size]   Calculating file sizes...\n');

% Size of files being archived
archive_bytes = 0;
for i = 1:numel(files_to_archive)
    info = dir(files_to_archive{i});
    if ~isempty(info)
        archive_bytes = archive_bytes + info(1).bytes;
    end
end

% Size of original RayStation directory (full export)
rs_total_bytes = dirSize(rs_dir);
% Size of SCT directory
sct_dir_path = fullfile(CONFIG.working_dir, 'EthosExports', CONFIG.patient_id, ...
    CONFIG.treatment_site, CONFIG.session, 'sct');
sct_total_bytes = 0;
if isfolder(sct_dir_path)
    sct_total_bytes = dirSize(sct_dir_path);
end
original_bytes = rs_total_bytes + sct_total_bytes;

fprintf('  Original export size:   %s  (RayStation dir + SCT dir)\n', ...
    formatBytes(original_bytes));
fprintf('  Files to archive:       %s  (processed .mat outputs)\n', ...
    formatBytes(archive_bytes));
if original_bytes > 0
    fprintf('  Raw reduction:          %.1f%%\n\n', ...
        (1 - archive_bytes / original_bytes) * 100);
end

%% ========================= CREATE ARCHIVE ================================

fprintf('[Archive] Creating compressed archive...\n');

% Determine format — fall back to zip if tar not available on Windows
use_tar = false;
if strcmpi(CONFIG.archive_format, 'tar.gz')
    [tar_status, ~] = system('tar --version');
    if tar_status == 0
        use_tar = true;
    else
        fprintf('  [WARN] tar not found on PATH — falling back to zip format.\n');
    end
end

if use_tar
    archive_path = fullfile(archive_dir, [archive_name '.tar.gz']);

    % Build file list string for tar (quote each path)
    file_list_str = '';
    for i = 1:numel(files_to_archive)
        file_list_str = [file_list_str ' "' files_to_archive{i} '"']; %#ok<AGROW>
    end

    tar_cmd = sprintf('tar -czf "%s"%s', archive_path, file_list_str);
    fprintf('  Running: tar -czf %s [%d files]\n', [archive_name '.tar.gz'], ...
        numel(files_to_archive));

    [status, output] = system(tar_cmd);
    if status ~= 0
        error('prepare_upload:TarFailed', ...
            'tar command failed:\n%s\n', output);
    end

else
    archive_path = fullfile(archive_dir, [archive_name '.zip']);
    fprintf('  Running: zip -> %s\n', [archive_name '.zip']);
    zip(archive_path, files_to_archive);
end

%% ========================= FINAL REPORT ==================================

archive_info = dir(archive_path);
compressed_bytes = archive_info(1).bytes;
compression_ratio = archive_bytes / max(compressed_bytes, 1);

fprintf('\n==========================================================\n');
fprintf('  UPLOAD PACKAGE READY\n');
fprintf('==========================================================\n');
fprintf('  Archive:            %s\n', archive_path);
fprintf('  Uncompressed size:  %s  (processed .mat files)\n', formatBytes(archive_bytes));
fprintf('  Compressed size:    %s  (%.1fx compression)\n', ...
    formatBytes(compressed_bytes), compression_ratio);
if original_bytes > 0
    fprintf('  vs. full export:    %s  ->  %s  (%.1f%% reduction)\n', ...
        formatBytes(original_bytes), formatBytes(compressed_bytes), ...
        (1 - compressed_bytes / original_bytes) * 100);
end
fprintf('\n');

% Upload time estimate
mbytes = compressed_bytes / 1e6;
upload_sec  = (mbytes * 8) / CONFIG.upload_mbps;
upload_min  = upload_sec / 60;
orig_upload_sec  = (original_bytes / 1e6 * 8) / CONFIG.upload_mbps;
orig_upload_min  = orig_upload_sec / 60;

fprintf('  Estimated upload @ %g Mbps\n', CONFIG.upload_mbps);
if orig_upload_min < 60
    fprintf('    Before: %.0f min\n', orig_upload_min);
else
    fprintf('    Before: %.1f hr\n', orig_upload_min / 60);
end
if upload_min < 60
    fprintf('    After:  %.0f min\n', upload_min);
else
    fprintf('    After:  %.1f hr\n', upload_min / 60);
end

fprintf('\n--- Contents breakdown ---\n');
fprintf('  %-40s  %10s\n', 'File', 'Size');
fprintf('  %s\n', repmat('-', 1, 55));
for i = 1:numel(files_to_archive)
    finfo = dir(files_to_archive{i});
    if isempty(finfo), continue; end
    [~, fname, fext] = fileparts(files_to_archive{i});
    fprintf('  %-40s  %10s\n', [fname fext], formatBytes(finfo(1).bytes));
end

fprintf('\n--- Recipient instructions ---\n');
fprintf('  1. Extract archive to:\n');
fprintf('       RayStationFiles/%s/%s/\n', CONFIG.patient_id, CONFIG.session);
fprintf('     The processed/ folder must land at:\n');
fprintf('       RayStationFiles/%s/%s/processed/*.mat\n', ...
    CONFIG.patient_id, CONFIG.session);
fprintf('  2. In pipeline_simulate.m, set:\n');
fprintf('       CONFIG.run_step15 = false;\n');
fprintf('  3. Run pipeline_simulate.m normally.\n');
fprintf('==========================================================\n');

%% ========================= LOCAL HELPER FUNCTIONS ========================

function total = dirSize(d)
%DIRSIZE  Recursively sum byte sizes of all files under directory d.
    total = 0;
    if ~isfolder(d), return; end
    listing = dir(fullfile(d, '**', '*'));
    for k = 1:numel(listing)
        if ~listing(k).isdir
            total = total + listing(k).bytes;
        end
    end
end

function s = formatBytes(b)
%FORMATBYTES  Return a human-readable string for a byte count.
    if b >= 1e9
        s = sprintf('%.2f GB', b / 1e9);
    elseif b >= 1e6
        s = sprintf('%.1f MB', b / 1e6);
    elseif b >= 1e3
        s = sprintf('%.1f KB', b / 1e3);
    else
        s = sprintf('%d B', round(b));
    end
end
