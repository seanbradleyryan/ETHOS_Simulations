function sct_dir = step0_sort_dicom(patient_id, session, config)
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

% Define output sct directory
sct_dir = fullfile(rawwd, 'sct');

fprintf('  Processing: Patient %s, %s\n', patient_id, session);
fprintf('  Raw directory: %s\n', rawwd);

%% ======================== VERIFY RAW DIRECTORY ========================

if ~isfolder(rawwd)
    warning('step0_sort_dicom:DirectoryNotFound', ...
        'Raw directory not found for patient %s, %s: %s', ...
        patient_id, session, rawwd);
    sct_dir = '';
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
    sct_dir = '';
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

%% ======================== VERIFY OUTPUT ========================

% Count files in sct directory
sctFiles = dir(fullfile(sct_dir, '*.dcm'));
fprintf('  Sorting complete. %d files in sct directory.\n', length(sctFiles));

% Verify critical files exist
hasRTSTRUCT = ~isempty(dir(fullfile(sct_dir, 'RS*.dcm')));
hasRTPLAN = ~isempty(dir(fullfile(sct_dir, 'RP*.dcm')));
hasRTDOSE = ~isempty(dir(fullfile(sct_dir, 'RD*.dcm')));
hasCT = ~isempty(dir(fullfile(sct_dir, 'CT*.dcm'))) || ...
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
