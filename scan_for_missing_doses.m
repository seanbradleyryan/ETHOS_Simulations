%% =========================================================================
%  SCAN_FOR_MISSING_DOSES.m
%  ETHOS Photoacoustic Pipeline — RayStation Dose Completeness Check
%  =========================================================================
%
%  PURPOSE:
%  Scan the RayStationFiles directory for a given patient/session and report
%  any missing segments from each beamset. Segments are expected to run
%  continuously from CONFIG.seg_range(1) to CONFIG.seg_range(2) inclusive.
%
%  INPUTS (via CONFIG struct below):
%    patient_id  - Patient identifier string  (e.g. '1194203')
%    session     - Session string             (e.g. 'Session_1')
%    working_dir - Repository root            (default: directory of this file)
%    seg_range   - [min max] segment range    (default: [0 164])
%
%  OUTPUT:
%    Console report listing each beamset, missing segments, and a summary.
%
%  EXAMPLE:
%    % Edit the CONFIG block below, then run:
%    scan_for_missing_doses
%
%  DEPENDENCIES: none
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: May 2026
%  =========================================================================

clear; clc;

%% ========================= CONFIGURATION =================================

CONFIG.patient_id  = '1194203';
CONFIG.session     = 'Session_1';
CONFIG.working_dir = fileparts(mfilename('fullpath'));
CONFIG.seg_range   = [0, 164];

%% ========================= SETUP =========================================

patient_id = char(CONFIG.patient_id);
session    = char(CONFIG.session);
seg_min    = CONFIG.seg_range(1);
seg_max    = CONFIG.seg_range(2);
n_expected = seg_max - seg_min + 1;

dose_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', patient_id, session);

fprintf('=========================================================\n');
fprintf('  SCAN_FOR_MISSING_DOSES\n');
fprintf('  Patient: %s  |  Session: %s\n', patient_id, session);
fprintf('=========================================================\n');
fprintf('  Directory: %s\n\n', dose_dir);

if ~isfolder(dose_dir)
    error('scan_for_missing_doses:DirNotFound', ...
        'Directory not found:\n  %s', dose_dir);
end

%% ========================= SCAN FILES ====================================

files = dir(fullfile(dose_dir, 'dose_*.dcm'));

if isempty(files)
    fprintf('  No dose_*.dcm files found.\n');
    fprintf('=========================================================\n');
    return
end

% Parse filenames: dose_{id}_{session}_{plan_type}_{beam}_{seg}.dcm
% Capture groups: (plan_type) (B\d+) (seg)
pattern = '_([^_]+)_(B\d+)_(\d+)\.dcm$';

beamset_keys = {};   % unique 'plan_type_beam' strings
beamset_segs = {};   % cell of segment number vectors per beamset

n_parsed = 0;
n_skipped = 0;

for i = 1:numel(files)
    tokens = regexp(files(i).name, pattern, 'tokens');
    if isempty(tokens)
        n_skipped = n_skipped + 1;
        continue
    end
    plan_type = tokens{1}{1};
    beam      = tokens{1}{2};
    seg_num   = str2double(tokens{1}{3});

    key = sprintf('%s_%s', plan_type, beam);
    idx = find(strcmp(beamset_keys, key), 1);
    if isempty(idx)
        beamset_keys{end+1} = key; %#ok<AGROW>
        beamset_segs{end+1} = seg_num; %#ok<AGROW>
    else
        beamset_segs{idx}(end+1) = seg_num;
    end
    n_parsed = n_parsed + 1;
end

fprintf('  Found %d dose files across %d beamset(s).', n_parsed, numel(beamset_keys));
if n_skipped > 0
    fprintf('  (%d file(s) skipped — name did not match pattern)', n_skipped);
end
fprintf('\n\n');

%% ========================= REPORT ========================================

n_incomplete = 0;
expected_segs = seg_min:seg_max;

% Sort beamsets alphabetically for consistent output
[beamset_keys, sort_idx] = sort(beamset_keys);
beamset_segs = beamset_segs(sort_idx);

for i = 1:numel(beamset_keys)
    key      = beamset_keys{i};
    present  = unique(beamset_segs{i});
    missing  = setdiff(expected_segs, present);

    fprintf('  Beamset: %-20s  [%d/%d segments present]\n', ...
        key, numel(present), n_expected);

    if isempty(missing)
        fprintf('    OK — all %d segments present (%d–%d)\n\n', ...
            n_expected, seg_min, seg_max);
    else
        n_incomplete = n_incomplete + 1;
        % Print missing segments in compact ranges
        fprintf('    MISSING %d segment(s): %s\n\n', ...
            numel(missing), format_seg_list(missing));
    end
end

fprintf('=========================================================\n');
if n_incomplete == 0
    fprintf('  Summary: all %d beamset(s) complete.\n', numel(beamset_keys));
else
    fprintf('  Summary: %d/%d beamset(s) incomplete.\n', ...
        n_incomplete, numel(beamset_keys));
end
fprintf('=========================================================\n');


%% =========================================================================
%  HELPER
%% =========================================================================

function str = format_seg_list(segs)
% Compact representation: consecutive runs shown as a-b, singles as a.
    segs = sort(segs(:)');
    parts = {};
    i = 1;
    while i <= numel(segs)
        j = i;
        while j < numel(segs) && segs(j+1) == segs(j) + 1
            j = j + 1;
        end
        if j > i
            parts{end+1} = sprintf('%d-%d', segs(i), segs(j)); %#ok<AGROW>
        else
            parts{end+1} = sprintf('%d', segs(i)); %#ok<AGROW>
        end
        i = j + 1;
    end
    str = strjoin(parts, ', ');
end
