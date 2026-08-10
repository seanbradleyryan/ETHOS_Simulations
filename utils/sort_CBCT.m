function cbct_paths = sort_CBCT(patient_id, session, config)
%% SORT_CBCT - Sort CBCT DICOM series into the sct directory
%
%   cbct_paths = sort_CBCT(patient_id, session, config)
%
%   PURPOSE:
%   Identify all CT series in the raw ETHOS DICOM export whose
%   SeriesDescription contains 'CBCT' and copy their files into the
%   patient/session sct directory (alongside the SCT, RT, and image
%   registration files sorted by step0_sort_dicom).  When more than two
%   CBCT series are present, the two earliest by SeriesDate/SeriesTime
%   are selected.
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir    - Base directory path
%           .treatment_site - Subfolder name (default: 'Pancreas')
%
%   OUTPUTS:
%       cbct_paths  - Cell array of destination file paths copied into the
%                     sct directory.  Empty if no CBCT series were found.
%
%   ALGORITHM:
%   1. Scan EthosExports/<patient>/<treatment_site>/<session>/ via
%      dicomCollection.
%   2. Select rows where Modality is 'CT' and SeriesDescription contains
%      'CBCT' (case-insensitive).
%   3. If more than 2 CBCT series are found, sort by
%      SeriesDate+SeriesTime and keep the two earliest.
%   4. Copy each selected series' files into the sct directory using a
%      'CBCT<n>_' prefix on the original filename to avoid collisions
%      with the SCT 'CT*.dcm' filenames.
%
%   FILE NAMING CONVENTION:
%       CBCT1_<original_name>.dcm
%       CBCT2_<original_name>.dcm
%
%   EXAMPLE:
%       config.working_dir    = 'C:/Users/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       cbct_paths = sort_CBCT('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (dicomCollection, dicominfo)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: April 2026
%
%   See also: step0_sort_dicom, step06_explode_segments, dicomCollection

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('sort_CBCT:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

if ~ischar(session) && ~isstring(session)
    error('sort_CBCT:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

if ~isstruct(config)
    error('sort_CBCT:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

if ~isfield(config, 'working_dir')
    error('sort_CBCT:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

if ~isfield(config, 'treatment_site') || isempty(config.treatment_site)
    config.treatment_site = 'Pancreas';
end

%% ======================== CONSTRUCT PATHS ========================

rawwd      = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session);
sct_dir    = fullfile(rawwd, 'sct');
cbct_paths = {};

fprintf('  Sorting CBCT series for: Patient %s, %s\n', patient_id, session);

if ~isfolder(rawwd)
    warning('sort_CBCT:DirectoryNotFound', ...
        'Raw directory not found for patient %s, %s: %s', ...
        patient_id, session, rawwd);
    return;
end

if ~isfolder(sct_dir)
    mkdir(sct_dir);
    fprintf('    Created sct directory: %s\n', sct_dir);
end

%% ======================== SCAN DICOM COLLECTION ========================

try
    ctInfo = dicomCollection(rawwd);
catch ME
    error('sort_CBCT:DicomScanFailed', ...
        'Failed to scan DICOM directory: %s\nError: %s', rawwd, ME.message);
end

if isempty(ctInfo) || height(ctInfo) == 0
    warning('sort_CBCT:EmptyCollection', 'No DICOM files found in: %s', rawwd);
    return;
end

if ~ismember('Modality', ctInfo.Properties.VariableNames) || ...
   ~ismember('SeriesDescription', ctInfo.Properties.VariableNames)
    warning('sort_CBCT:MissingColumns', ...
        'Modality or SeriesDescription column missing from DICOM collection');
    return;
end

%% ======================== SELECT CBCT SERIES ========================

isCT     = strcmpi(ctInfo.Modality, 'CT');
descUpper = upper(string(ctInfo.SeriesDescription));
isCBCT   = contains(descUpper, 'CBCT');
matches  = isCT & isCBCT;

if ~any(matches)
    fprintf('    No CBCT series found for %s/%s.\n', patient_id, session);
    return;
end

cbctInfo = ctInfo(matches, :);
nFound   = height(cbctInfo);
fprintf('    Found %d CBCT series\n', nFound);

%% ======================== ENFORCE MAX 2 (BY DATETIME) ========================

if nFound > 2
    fprintf('    More than 2 CBCT series found; selecting the two earliest by datetime.\n');

    dtValues = NaN(nFound, 1);
    for i = 1:nFound
        fileCell = cbctInfo.Filenames{i};
        if isempty(fileCell) || isempty(fileCell{1}), continue; end
        try
            meta    = dicominfo(fileCell{1});
            dateStr = '';
            timeStr = '';
            if isfield(meta, 'SeriesDate'), dateStr = strtrim(meta.SeriesDate); end
            if isfield(meta, 'SeriesTime'), timeStr = strtrim(meta.SeriesTime); end
            if ~isempty(dateStr) || ~isempty(timeStr)
                dt = str2double([dateStr, timeStr]);
                if ~isnan(dt) && dt > 0
                    dtValues(i) = dt;
                end
            end
        catch
            % leave as NaN
        end
    end

    [~, sortIdx] = sort(dtValues);   % NaNs sort to the end
    keepIdx      = sortIdx(1:2);
    fprintf('    Keeping CBCT series indices [%d, %d] with datetimes [%s, %s]\n', ...
        keepIdx(1), keepIdx(2), num2str(dtValues(keepIdx(1))), num2str(dtValues(keepIdx(2))));
    cbctInfo = cbctInfo(keepIdx, :);
end

%% ======================== COPY FILES TO SCT DIR ========================

totalCopied = 0;
for i = 1:height(cbctInfo)
    fileCell = cbctInfo.Filenames{i};
    if isempty(fileCell), continue; end

    desc = '';
    seriesDate = '';
    seriesTime = '';
    try
        if ~isempty(fileCell{1})
            meta = dicominfo(fileCell{1});
            if isfield(meta, 'SeriesDescription')
                desc = strtrim(meta.SeriesDescription);
            end
            if isfield(meta, 'SeriesDate'), seriesDate = strtrim(meta.SeriesDate); end
            if isfield(meta, 'SeriesTime'), seriesTime = strtrim(meta.SeriesTime); end
        end
    catch
    end
    fprintf('    CBCT series %d: "%s"  Date/Time: %s / %s  (%d files)\n', ...
        i, desc, seriesDate, seriesTime, numel(fileCell));

    prefix   = sprintf('CBCT%d_', i);
    nCopied  = 0;
    nSkipped = 0;

    for k = 1:numel(fileCell)
        srcFile = fileCell{k};
        if isempty(srcFile) || ~isfile(srcFile), continue; end

        [~, name, ext] = fileparts(srcFile);
        destFile = fullfile(sct_dir, [prefix, name, ext]);

        if isfile(destFile)
            nSkipped = nSkipped + 1;
            cbct_paths{end+1} = destFile; %#ok<AGROW>
            continue;
        end

        try
            copyfile(srcFile, destFile);
            cbct_paths{end+1} = destFile; %#ok<AGROW>
            nCopied = nCopied + 1;
        catch ME
            warning('sort_CBCT:CopyError', ...
                'Failed to copy %s: %s', name, ME.message);
        end
    end
    fprintf('      Copied %d, skipped %d (already present)\n', nCopied, nSkipped);
    totalCopied = totalCopied + nCopied;
end

fprintf('  CBCT sorting complete. Total CBCT files copied: %d\n', totalCopied);

end
