function [sct_dir, sim_ct_dir] = step0_sort_dicom(patient_id, session, config)
%% STEP0_SORT_DICOM - Sort DICOM files for ETHOS pipeline
%
%   sct_dir = step0_sort_dicom(patient_id, session, config)
%
%   PURPOSE:
%   Organize raw ETHOS DICOM export by identifying SCT (synthetic CT) series
%   and matching RT files (RTSTRUCT, RTPLAN, RTDOSE) based on DICOM reference
%   chains. Copy matched files to a clean 'sct' subdirectory for downstream
%   processing.
%
%   INPUTS:
%       patient_id  - String, patient identifier (e.g., '1194203')
%       session     - String, session name (e.g., 'Session_1')
%       config      - Struct with configuration parameters:
%           .working_dir    - Base directory path
%           .treatment_site - Subfolder name (default: 'Pancreas')
%
%   OUTPUTS:
%       sct_dir     - String, path to directory containing sorted files:
%                     - CT*.dcm files from SCT series
%                     - Matched RS*.dcm (RTSTRUCT)
%                     - Matched RP*.dcm (RTPLAN)
%                     - Matched RD*.dcm (RTDOSE)
%       sim_ct_dir  - String, path to directory containing the most recent
%                     simulation CT series (original planning CT, not SCT).
%                     Searches for CT series whose SeriesDescription is NOT
%                     'sct', preferring descriptions that contain 'sim'.
%                     Returns '' if no simulation CT is found.
%
%   ALGORITHM:
%   1. Create 'sct' subdirectory if not exists
%   2. Use dicomCollection() to scan directory
%   3. Find SCT series (SeriesDescription = 'sct')
%   4. Extract SCT SeriesInstanceUID
%   5. Find RTSTRUCT with ReferencedFrameOfReferenceSequence matching SCT
%   6. Find RTPLAN with ReferencedStructureSetSequence matching RTSTRUCT
%   7. Find RTDOSE with ReferencedRTPlanSequence matching RTPLAN
%   8. Copy matched files to sct directory
%   9. Return sct_dir path
%
%   EXAMPLE:
%       config.working_dir = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       config.treatment_site = 'Pancreas';
%       sct_dir = step0_sort_dicom('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (dicomCollection, dicominfo)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 2.0 (Refactored for master pipeline integration)
%
%   See also: dicomCollection, dicominfo, step05_fix_mlc_gaps

%% ======================== INPUT VALIDATION ========================

% Validate patient_id
if ~ischar(patient_id) && ~isstring(patient_id)
    error('step0_sort_dicom:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);  % Ensure char array for consistency

% Validate session
if ~ischar(session) && ~isstring(session)
    error('step0_sort_dicom:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

% Validate config struct
if ~isstruct(config)
    error('step0_sort_dicom:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

% Validate required config fields
if ~isfield(config, 'working_dir')
    error('step0_sort_dicom:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

% Set default treatment_site if not provided
if ~isfield(config, 'treatment_site') || isempty(config.treatment_site)
    config.treatment_site = 'Pancreas';
    fprintf('  [INFO] Using default treatment_site: %s\n', config.treatment_site);
end

% Validate working directory exists
if ~isfolder(config.working_dir)
    error('step0_sort_dicom:DirectoryNotFound', ...
        'Working directory does not exist: %s', config.working_dir);
end

%% ======================== CONSTRUCT PATHS ========================

% Build path to raw DICOM directory
rawwd = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session);

% Define output directories
sct_dir    = fullfile(rawwd, 'sct');
sim_ct_dir = '';   % populated later if a simulation CT is found

fprintf('  Processing: Patient %s, %s\n', patient_id, session);
fprintf('  Raw directory: %s\n', rawwd);

%% ======================== VERIFY RAW DIRECTORY ========================

if ~isfolder(rawwd)
    warning('step0_sort_dicom:DirectoryNotFound', ...
        'Raw directory not found for patient %s, %s: %s', ...
        patient_id, session, rawwd);
    sct_dir    = '';
    sim_ct_dir = '';
    return;
end

%% ======================== SCAN DICOM COLLECTION ========================

fprintf('  Scanning DICOM collection...\n');

try
    ctInfo = dicomCollection(rawwd);
catch ME
    error('step0_sort_dicom:DicomScanFailed', ...
        'Failed to scan DICOM directory: %s\nError: %s', rawwd, ME.message);
end

if isempty(ctInfo) || height(ctInfo) == 0
    warning('step0_sort_dicom:EmptyCollection', ...
        'No DICOM files found in: %s', rawwd);
    sct_dir    = '';
    sim_ct_dir = '';
    return;
end

fprintf('  Found %d DICOM series\n', height(ctInfo));

%% ======================== PRINT DIAGNOSTIC INFO ========================

printCollectionInfo(ctInfo);

%% ======================== CREATE SCT DIRECTORY ========================

if ~isfolder(sct_dir)
    mkdir(sct_dir);
    fprintf('  Created sct directory: %s\n', sct_dir);
else
    fprintf('  sct directory exists: %s\n', sct_dir);
end

%% ======================== SORT SCT FILES ========================

fprintf('  Sorting SCT files...\n');
sctSeriesUID = sortSctFiles(ctInfo, rawwd, sct_dir);

if isempty(sctSeriesUID)
    warning('step0_sort_dicom:NoSCT', ...
        'No SCT series found for patient %s, %s', patient_id, session);
end

%% ======================== SORT RT FILES ========================

fprintf('  Sorting RT files (RTSTRUCT, RTPLAN, RTDOSE)...\n');
sortRTFiles(ctInfo, rawwd, sct_dir, sctSeriesUID);

%% ======================== SORT SIMULATION CT FILES ========================

fprintf('  Searching for simulation CT series...\n');
sim_ct_dir = fullfile(rawwd, 'sim_ct');
sim_ct_dir = sortSimCtFiles(ctInfo, rawwd, sim_ct_dir, sctSeriesUID);

%% ======================== VERIFY OUTPUT ========================

% Count files in sct directory
sctFiles = dir(fullfile(sct_dir, '*.dcm'));
fprintf('  Sorting complete. %d files in sct directory.\n', length(sctFiles));

% Verify critical files exist
hasRTSTRUCT = ~isempty(dir(fullfile(sct_dir, 'RS*.dcm')));
hasRTPLAN   = ~isempty(dir(fullfile(sct_dir, 'RP*.dcm')));
hasRTDOSE   = ~isempty(dir(fullfile(sct_dir, 'RD*.dcm')));
hasCT       = ~isempty(dir(fullfile(sct_dir, 'CT*.dcm'))) || ...
              ~isempty(dir(fullfile(sct_dir, '*CT*.dcm')));

if ~hasRTSTRUCT
    warning('step0_sort_dicom:MissingFile', 'No RTSTRUCT file found in sct directory');
end
if ~hasRTPLAN
    warning('step0_sort_dicom:MissingFile', 'No RTPLAN file found in sct directory');
end
if ~hasRTDOSE
    warning('step0_sort_dicom:MissingFile', 'No RTDOSE file found in sct directory');
end
if ~hasCT
    warning('step0_sort_dicom:MissingFile', 'No CT files found in sct directory');
end

% Report sim_ct result
if ~isempty(sim_ct_dir)
    simCtFiles = dir(fullfile(sim_ct_dir, '*.dcm'));
    fprintf('  Simulation CT: %d files in %s\n', length(simCtFiles), sim_ct_dir);
else
    warning('step0_sort_dicom:NoSimCT', ...
        'No simulation CT series found for patient %s, %s', patient_id, session);
end

fprintf('  Step 0 complete for %s/%s\n', patient_id, session);

end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function sctSeriesUID = sortSctFiles(ctInfo, sourceDir, destDir)
%SORTSCTFILES Sort SCT DICOM files and extract SeriesInstanceUID
%
%   sctSeriesUID = sortSctFiles(ctInfo, sourceDir, destDir)
%
%   Find SCT series by SeriesDescription, move files to destination,
%   and return the SeriesInstanceUID for RT file matching.

    sctSeriesUID = '';
    
    % Find rows with SeriesDescription = 'sct'
    if ~ismember('SeriesDescription', ctInfo.Properties.VariableNames)
        warning('sortSctFiles:NoSeriesDescription', ...
            'SeriesDescription column not found in DICOM collection');
        return;
    end
    
    rowIndex = strcmp(ctInfo.SeriesDescription, 'sct');
    sctInfo = ctInfo(rowIndex, :);
    
    if height(sctInfo) == 0
        % Try case-insensitive search
        rowIndex = strcmpi(ctInfo.SeriesDescription, 'sct');
        sctInfo = ctInfo(rowIndex, :);
    end
    
    if height(sctInfo) == 0
        warning('sortSctFiles:NoSCT', 'No SCT series found in collection');
        return;
    end
    
    fprintf('    Found %d SCT series\n', height(sctInfo));
    
    % Get files from first (or only) SCT series
    sctFiles = sctInfo.Filenames{1};
    
    if isempty(sctFiles)
        warning('sortSctFiles:EmptySeries', 'SCT series has no files');
        return;
    end
    
    % Extract SeriesInstanceUID from first file
    firstFile = sctFiles{1};
    try
        sctMetadata = dicominfo(firstFile);
        sctSeriesUID = sctMetadata.SeriesInstanceUID;
        fprintf('    SCT SeriesInstanceUID: %s\n', sctSeriesUID);
        fprintf('    SCT Series Date/Time: %s / %s\n', ...
            sctMetadata.SeriesDate, sctMetadata.SeriesTime);
    catch ME
        warning('sortSctFiles:MetadataError', ...
            'Failed to read SCT metadata: %s', ME.message);
        return;
    end
    
    % Move SCT files to destination
    numMoved = 0;
    numSkipped = 0;
    
    for k = 1:length(sctFiles)
        srcFile = sctFiles{k};
        [~, name, ext] = fileparts(srcFile);
        destFile = fullfile(destDir, strcat(name, ext));
        
        if exist(srcFile, 'file')
            if exist(destFile, 'file')
                numSkipped = numSkipped + 1;
            else
                try
                    movefile(srcFile, destFile);
                    numMoved = numMoved + 1;
                catch ME
                    warning('sortSctFiles:MoveError', ...
                        'Failed to move file %s: %s', name, ME.message);
                end
            end
        end
    end
    
    fprintf('    SCT files: %d moved, %d already existed\n', numMoved, numSkipped);
end


function sim_ct_dir = sortSimCtFiles(ctInfo, sourceDir, sim_ct_dir, sctSeriesUID)
%SORTSIMCTFILES Find and move the most recent simulation CT to sim_ct folder
%
%   sim_ct_dir = sortSimCtFiles(ctInfo, sourceDir, sim_ct_dir, sctSeriesUID)
%
%   Searches for a CT series whose SeriesDescription is NOT 'sct'.
%   Priority:
%     1. CT series with SeriesDescription containing 'sim' (case-insensitive)
%     2. Most recent remaining CT series
%   Moves all files from the selected series to sim_ct_dir and returns
%   that path.  Returns '' if no eligible series is found.

    sim_ct_dir = '';

    if ~ismember('Modality', ctInfo.Properties.VariableNames) || ...
       ~ismember('SeriesDescription', ctInfo.Properties.VariableNames)
        warning('sortSimCtFiles:MissingColumns', ...
            'Modality or SeriesDescription column not found in DICOM collection');
        return;
    end

    % Find all CT rows that are NOT the SCT series
    isCT    = strcmpi(ctInfo.Modality, 'CT');
    isSct   = strcmpi(ctInfo.SeriesDescription, 'sct');
    eligible = isCT & ~isSct;

    if ~any(eligible)
        fprintf('    No non-SCT CT series found.\n');
        return;
    end

    eligibleInfo = ctInfo(eligible, :);
    fprintf('    Found %d non-SCT CT series\n', height(eligibleInfo));

    % Prefer series whose description contains 'sim'
    hasSim = contains(lower(eligibleInfo.SeriesDescription), 'sim');
    if any(hasSim)
        candidateInfo = eligibleInfo(hasSim, :);
        fprintf('    Found %d series with ''sim'' in description\n', height(candidateInfo));
    else
        candidateInfo = eligibleInfo;
    end

    % Selection priority among candidates:
    %   1. Series with an empty/missing SeriesDate or SeriesTime — ETHOS
    %      sometimes exports the simulation CT without a date stamp.
    %   2. Oldest series by SeriesDate/SeriesTime (simulation CT is acquired
    %      before the adaptive fractions, so it has the earliest timestamp).
    %
    % Build a datetime value for every candidate series.
    % NaN is used as a sentinel for missing/empty datetime (highest priority).
    numCandidates = height(candidateInfo);
    dtValues      = zeros(numCandidates, 1);   % 0 = placeholder

    for i = 1:numCandidates
        fileCell = candidateInfo.Filenames{i};
        if isempty(fileCell) || isempty(fileCell{1})
            dtValues(i) = NaN;   % treat unreadable file same as missing date
            continue;
        end
        try
            meta     = dicominfo(fileCell{1});
            dateStr  = '';
            timeStr  = '';
            if isfield(meta, 'SeriesDate'), dateStr = strtrim(meta.SeriesDate); end
            if isfield(meta, 'SeriesTime'), timeStr = strtrim(meta.SeriesTime); end

            if isempty(dateStr) && isempty(timeStr)
                % Empty datetime — highest priority candidate
                dtValues(i) = NaN;
            else
                dt = str2double([dateStr, timeStr]);
                if isnan(dt) || dt == 0
                    dtValues(i) = NaN;   % treat unparseable as missing
                else
                    dtValues(i) = dt;
                end
            end
        catch
            dtValues(i) = NaN;   % unreadable metadata → treat as missing date
        end
    end

    % Choose: first prefer any NaN (empty datetime) entry, then the oldest
    nanIdx = find(isnan(dtValues));
    if ~isempty(nanIdx)
        % Multiple NaN entries — just take the first
        bestIdx = nanIdx(1);
        fprintf('    Selecting simulation CT with empty/missing datetime (index %d)\n', bestIdx);
    else
        % All have valid timestamps — pick the oldest (smallest value)
        [~, bestIdx] = min(dtValues);
        fprintf('    Selecting oldest simulation CT by SeriesDate/Time\n');
    end

    selectedSeries = candidateInfo(bestIdx, :);
    simFiles       = selectedSeries.Filenames{1};

    if isempty(simFiles)
        fprintf('    Selected simulation CT series has no files.\n');
        return;
    end

    % Log selected series info
    try
        meta = dicominfo(simFiles{1});
        fprintf('    Simulation CT: "%s"  Date: %s  Files: %d\n', ...
            meta.SeriesDescription, meta.SeriesDate, length(simFiles));
    catch
        fprintf('    Simulation CT: %d files (metadata unavailable)\n', length(simFiles));
    end

    % Create destination directory (rawwd/sim_ct, passed in as sim_ct_dir arg)
    destDir = fullfile(sourceDir, 'sim_ct');

    if ~isfolder(destDir)
        mkdir(destDir);
        fprintf('    Created sim_ct directory: %s\n', destDir);
    else
        fprintf('    sim_ct directory exists: %s\n', destDir);
    end

    % Move files
    numMoved   = 0;
    numSkipped = 0;

    for k = 1:length(simFiles)
        srcFile = simFiles{k};
        [~, name, ext] = fileparts(srcFile);
        destFile = fullfile(destDir, [name, ext]);

        if ~exist(srcFile, 'file')
            continue;
        end
        if exist(destFile, 'file')
            numSkipped = numSkipped + 1;
        else
            try
                movefile(srcFile, destFile);
                numMoved = numMoved + 1;
            catch ME
                warning('sortSimCtFiles:MoveError', ...
                    'Failed to move %s: %s', name, ME.message);
            end
        end
    end

    fprintf('    Simulation CT files: %d moved, %d already existed\n', ...
        numMoved, numSkipped);

    sim_ct_dir = destDir;
end


function sortRTFiles(ctInfo, sourceDir, destDir, sctSeriesUID)
%SORTRTFILES Sort RT DICOM files based on reference chain
%
%   sortRTFiles(ctInfo, sourceDir, destDir, sctSeriesUID)
%
%   Find and copy RTSTRUCT, RTPLAN, and RTDOSE files that form a
%   complete reference chain matching the SCT series.

    modalityList = {'RTSTRUCT', 'RTPLAN', 'RTDOSE'};
    
    % Build structure of RT file information with metadata
    rtFileInfo = struct();
    
    for modIdx = 1:length(modalityList)
        modality = modalityList{modIdx};
        
        % Find rows matching this modality
        if ismember('Modality', ctInfo.Properties.VariableNames)
            rowIndices = strcmp(ctInfo.Modality, modality);
        else
            continue;
        end
        
        if ~any(rowIndices)
            fprintf('    No %s files found\n', modality);
            continue;
        end
        
        selectedRows = ctInfo(rowIndices, :);
        rtFileInfo.(modality) = {};
        
        % Extract metadata for each file
        for fileIdx = 1:height(selectedRows)
            fileCell = selectedRows.Filenames{fileIdx};
            
            if ~isempty(fileCell) && ~isempty(fileCell{1})
                filePath = fileCell{1};
                
                try
                    metadata = dicominfo(filePath);
                    
                    fileInfo = struct();
                    fileInfo.FilePath = filePath;
                    fileInfo.SeriesInstanceUID = metadata.SeriesInstanceUID;
                    fileInfo.SeriesDate = metadata.SeriesDate;
                    fileInfo.SeriesTime = metadata.SeriesTime;
                    
                    % Extract SeriesDescription if available
                    if isfield(metadata, 'SeriesDescription')
                        fileInfo.SeriesDescription = metadata.SeriesDescription;
                    else
                        fileInfo.SeriesDescription = 'N/A';
                    end
                    
                    % Extract reference UIDs based on modality
                    fileInfo = extractReferenceUIDs(fileInfo, metadata, modality);
                    
                    rtFileInfo.(modality){end+1} = fileInfo;
                    
                catch ME
                    warning('sortRTFiles:MetadataError', ...
                        'Failed to read %s metadata from %s: %s', ...
                        modality, filePath, ME.message);
                end
            end
        end
        
        fprintf('    Found %d %s file(s)\n', length(rtFileInfo.(modality)), modality);
    end
    
    % Find best matching set of RT files
    [selectedRTStruct, selectedRTPlan, selectedRTDose] = ...
        findMatchingRTSet(rtFileInfo, sctSeriesUID);
    
    % Copy selected files to destination
    if ~isempty(selectedRTStruct)
        copyRTFile(selectedRTStruct, destDir, 'RTSTRUCT');
    end
    
    if ~isempty(selectedRTPlan)
        copyRTFile(selectedRTPlan, destDir, 'RTPLAN');
    end
    
    if ~isempty(selectedRTDose)
        copyRTFile(selectedRTDose, destDir, 'RTDOSE');
    end
end


function fileInfo = extractReferenceUIDs(fileInfo, metadata, modality)
%EXTRACTREFERENCEUIDS Extract reference UIDs from DICOM metadata
%
%   Extract the appropriate reference UIDs based on modality:
%   - RTSTRUCT references CT Series
%   - RTPLAN references RTSTRUCT
%   - RTDOSE references RTPLAN

    switch modality
        case 'RTSTRUCT'
            % RTSTRUCT references CT series via ReferencedFrameOfReferenceSequence
            if isfield(metadata, 'ReferencedFrameOfReferenceSequence')
                try
                    refSeq = metadata.ReferencedFrameOfReferenceSequence;
                    if isfield(refSeq, 'Item_1') && ...
                       isfield(refSeq.Item_1, 'RTReferencedStudySequence')
                        studySeq = refSeq.Item_1.RTReferencedStudySequence;
                        if isfield(studySeq, 'Item_1') && ...
                           isfield(studySeq.Item_1, 'RTReferencedSeriesSequence')
                            seriesSeq = studySeq.Item_1.RTReferencedSeriesSequence;
                            if isfield(seriesSeq, 'Item_1') && ...
                               isfield(seriesSeq.Item_1, 'SeriesInstanceUID')
                                fileInfo.ReferencedSeriesUID = seriesSeq.Item_1.SeriesInstanceUID;
                            end
                        end
                    end
                catch
                    % Reference extraction failed, continue without it
                end
            end
            
        case 'RTPLAN'
            % RTPLAN references RTSTRUCT via ReferencedStructureSetSequence
            if isfield(metadata, 'ReferencedStructureSetSequence')
                try
                    refSeq = metadata.ReferencedStructureSetSequence;
                    if isfield(refSeq, 'Item_1') && ...
                       isfield(refSeq.Item_1, 'ReferencedSOPInstanceUID')
                        fileInfo.ReferencedSOPInstanceUID = refSeq.Item_1.ReferencedSOPInstanceUID;
                    end
                catch
                    % Reference extraction failed
                end
            end
            
        case 'RTDOSE'
            % RTDOSE references RTPLAN via ReferencedRTPlanSequence
            if isfield(metadata, 'ReferencedRTPlanSequence')
                try
                    refSeq = metadata.ReferencedRTPlanSequence;
                    if isfield(refSeq, 'Item_1') && ...
                       isfield(refSeq.Item_1, 'ReferencedSOPInstanceUID')
                        fileInfo.ReferencedSOPInstanceUID = refSeq.Item_1.ReferencedSOPInstanceUID;
                    end
                catch
                    % Reference extraction failed
                end
            end
    end
end


function [selectedStruct, selectedPlan, selectedDose] = findMatchingRTSet(rtFileInfo, sctSeriesUID)
%FINDMATCHINGRTSET Find a matching set of RTSTRUCT, RTPLAN, RTDOSE
%
%   Priority order:
%   1. RTSTRUCT that references the SCT SeriesInstanceUID
%   2. Most recent complete set
%   3. Most recent individual files

    selectedStruct = [];
    selectedPlan = [];
    selectedDose = [];
    
    % Get lists (may be empty)
    rtStructs = {};
    rtPlans = {};
    rtDoses = {};
    
    if isfield(rtFileInfo, 'RTSTRUCT')
        rtStructs = rtFileInfo.RTSTRUCT;
    end
    if isfield(rtFileInfo, 'RTPLAN')
        rtPlans = rtFileInfo.RTPLAN;
    end
    if isfield(rtFileInfo, 'RTDOSE')
        rtDoses = rtFileInfo.RTDOSE;
    end
    
    % Strategy 1: Find RTSTRUCT that references SCT
    if ~isempty(sctSeriesUID) && ~isempty(rtStructs)
        for i = 1:length(rtStructs)
            if isfield(rtStructs{i}, 'ReferencedSeriesUID') && ...
               strcmp(rtStructs{i}.ReferencedSeriesUID, sctSeriesUID)
                selectedStruct = rtStructs{i};
                fprintf('    Found RTSTRUCT referencing SCT\n');
                break;
            end
        end
    end
    
    % Fallback: Most recent RTSTRUCT
    if isempty(selectedStruct) && ~isempty(rtStructs)
        selectedStruct = getMostRecentFile(rtStructs);
        fprintf('    Selected most recent RTSTRUCT\n');
    end
    
    % Find RTPLAN that references selected RTSTRUCT
    if ~isempty(selectedStruct) && ~isempty(rtPlans)
        try
            structSOPInstanceUID = dicominfo(selectedStruct.FilePath).SOPInstanceUID;
            for i = 1:length(rtPlans)
                if isfield(rtPlans{i}, 'ReferencedSOPInstanceUID') && ...
                   strcmp(rtPlans{i}.ReferencedSOPInstanceUID, structSOPInstanceUID)
                    selectedPlan = rtPlans{i};
                    fprintf('    Found RTPLAN referencing selected RTSTRUCT\n');
                    break;
                end
            end
        catch
            % Could not read RTSTRUCT SOPInstanceUID
        end
    end
    
    % Fallback: Most recent RTPLAN
    if isempty(selectedPlan) && ~isempty(rtPlans)
        selectedPlan = getMostRecentFile(rtPlans);
        fprintf('    Selected most recent RTPLAN\n');
    end
    
    % Find RTDOSE that references selected RTPLAN
    if ~isempty(selectedPlan) && ~isempty(rtDoses)
        try
            planSOPInstanceUID = dicominfo(selectedPlan.FilePath).SOPInstanceUID;
            for i = 1:length(rtDoses)
                if isfield(rtDoses{i}, 'ReferencedSOPInstanceUID') && ...
                   strcmp(rtDoses{i}.ReferencedSOPInstanceUID, planSOPInstanceUID)
                    selectedDose = rtDoses{i};
                    fprintf('    Found RTDOSE referencing selected RTPLAN\n');
                    break;
                end
            end
        catch
            % Could not read RTPLAN SOPInstanceUID
        end
    end
    
    % Fallback: Most recent RTDOSE
    if isempty(selectedDose) && ~isempty(rtDoses)
        selectedDose = getMostRecentFile(rtDoses);
        fprintf('    Selected most recent RTDOSE\n');
    end
end


function mostRecent = getMostRecentFile(fileList)
%GETMOSTRECENTFILE Get most recent file based on SeriesDate/Time
%
%   mostRecent = getMostRecentFile(fileList)

    mostRecent = [];
    
    if isempty(fileList)
        return;
    end
    
    mostRecent = fileList{1};
    mostRecentDateTime = parseDateTime(mostRecent.SeriesDate, mostRecent.SeriesTime);
    
    for i = 2:length(fileList)
        currentDateTime = parseDateTime(fileList{i}.SeriesDate, fileList{i}.SeriesTime);
        if currentDateTime > mostRecentDateTime
            mostRecent = fileList{i};
            mostRecentDateTime = currentDateTime;
        end
    end
end


function dt = parseDateTime(dateStr, timeStr)
%PARSEDATETIME Parse DICOM date/time strings to numeric value
%
%   dt = parseDateTime(dateStr, timeStr)

    try
        dt = str2double([dateStr, timeStr]);
    catch
        dt = 0;
    end
    
    if isnan(dt)
        dt = 0;
    end
end


function copyRTFile(fileInfo, destDir, modality)
%COPYRTFILE Copy RT file to destination directory
%
%   copyRTFile(fileInfo, destDir, modality)

    [~, name, ext] = fileparts(fileInfo.FilePath);
    destFile = fullfile(destDir, strcat(name, ext));
    
    if exist(destFile, 'file')
        fprintf('    %s already exists in destination (skipping)\n', modality);
        return;
    end
    
    try
        copyfile(fileInfo.FilePath, destFile);
        fprintf('    Copied %s (Date: %s)\n', modality, fileInfo.SeriesDate);
    catch ME
        warning('copyRTFile:CopyError', ...
            'Failed to copy %s: %s', modality, ME.message);
    end
end


function printCollectionInfo(ctInfo)
%PRINTCOLLECTIONINFO Print diagnostic information about DICOM collection
%
%   printCollectionInfo(ctInfo)

    fprintf('\n  --- DICOM Collection Summary ---\n');
    fprintf('  Total series: %d\n', height(ctInfo));
    
    % Group by modality
    if ismember('Modality', ctInfo.Properties.VariableNames)
        modalities = unique(ctInfo.Modality);
        for i = 1:length(modalities)
            modality = modalities{i};
            count = sum(strcmp(ctInfo.Modality, modality));
            fprintf('    %s: %d series\n', modality, count);
        end
    end
    
    % Print SCT info
    fprintf('\n  --- SCT Series ---\n');
    if ismember('SeriesDescription', ctInfo.Properties.VariableNames)
        sctRows = strcmp(ctInfo.SeriesDescription, 'sct');
        if any(sctRows)
            sctInfo = ctInfo(sctRows, :);
            fileCell = sctInfo.Filenames{1};
            if ~isempty(fileCell) && ~isempty(fileCell{1})
                try
                    metadata = dicominfo(fileCell{1});
                    fprintf('    Series: %s\n', metadata.SeriesDescription);
                    fprintf('    Date/Time: %s / %s\n', metadata.SeriesDate, metadata.SeriesTime);
                    fprintf('    SeriesInstanceUID: %s\n', metadata.SeriesInstanceUID);
                    fprintf('    Number of images: %d\n', length(fileCell));
                catch
                    fprintf('    (metadata unavailable)\n');
                end
            end
        else
            fprintf('    No SCT series found\n');
        end
    end
    
    % Print RT file summary
    fprintf('\n  --- RT File Summary ---\n');
    rtModalities = {'RTSTRUCT', 'RTPLAN', 'RTDOSE'};
    
    for i = 1:length(rtModalities)
        modality = rtModalities{i};
        if ismember('Modality', ctInfo.Properties.VariableNames)
            rowIndices = strcmp(ctInfo.Modality, modality);
            if any(rowIndices)
                fprintf('    %s: %d found\n', modality, sum(rowIndices));
            end
        end
    end
    
    fprintf('  ------------------------------\n\n');
end
