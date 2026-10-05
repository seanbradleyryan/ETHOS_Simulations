function sct_dir = step0_sort_dicom(patient_id, session, config)
%% STEP0_SORT_DICOM - Sort DICOM files for ETHOS pipeline
%
%   sct_dir = step0_sort_dicom(patient_id, session, config)
%
%   PURPOSE:
%   Organize raw ETHOS DICOM export by identifying SCT (synthetic CT) series
%   and classifying RT plan files as REFERENCE or ADAPTED based on the
%   RTPlanRelationship field inside ReferenceRTPlanSequence.Item_1.
%   Each plan type's associated RTSTRUCT and RTDOSE are traced through the
%   standard DICOM reference chain and copied with standardized filenames.
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
%                     - RTSTRUCT_SCT.dcm   (plan structure set on the SCT;
%                                           _2 if the plans use two)
%                     - RTSTRUCT_extra_<n>.dcm (every other RTSTRUCT,
%                                           e.g. CBCT structure sets)
%                     - CBCT<n>_*.dcm      (every CBCT series, via sort_CBCT)
%                     - RTPLAN_reference.dcm  (reference RTPLAN)
%                     - RTPLAN_adapted.dcm    (adapted RTPLAN)
%                     - RTDOSE_reference.dcm  (RTDOSE for reference plan)
%                     - RTDOSE_adapted.dcm    (RTDOSE for adapted plan)
%
%   ALGORITHM:
%   1. Sort SCT series files (three-tier priority)
%   2. Scan all RTPLAN files; for each, read ReferenceRTPlanSequence.Item_1.RTPlanRelationship
%      - Contains 'REFERENCE'   reference plan
%      - Contains 'ADAPTED'     adapted plan
%   3. For each classified plan:
%      a. Trace ReferencedStructureSetSequence  find matching RTSTRUCT by SOPInstanceUID
%      b. Keep the RTSTRUCT as RTSTRUCT_SCT.dcm only if it references the SCT
%      c. Trace RTDOSE ReferencedRTPlanSequence  find RTDOSE referencing this plan
%   4. Copy all files with standardized names into sct_dir
%   5. Copy every CBCT series (sort_CBCT) and every non-SCT RTSTRUCT
%      (RTSTRUCT_extra_<n>.dcm) into sct_dir
%
%   FILE NAMING CONVENTION:
%       RTPLAN_reference.dcm   RTPLAN_adapted.dcm
%       RTSTRUCT_SCT.dcm       RTSTRUCT_extra_<n>.dcm
%       RTDOSE_reference.dcm   RTDOSE_adapted.dcm
%
%   EXAMPLE:
%       config.working_dir = get_repo_root();
%       config.treatment_site = 'Pancreas';
%       sct_dir = step0_sort_dicom('1194203', 'Session_1', config);
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (dicomCollection, dicominfo)
%
%   DATE: April 2026
%   VERSION: 3.0 (REFERENCE/ADAPTED plan classification)

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('step0_sort_dicom:InvalidInput', ...
        'patient_id must be a string or character array. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

if ~ischar(session) && ~isstring(session)
    error('step0_sort_dicom:InvalidInput', ...
        'session must be a string or character array. Received: %s', class(session));
end
session = char(session);

if ~isstruct(config)
    error('step0_sort_dicom:InvalidInput', ...
        'config must be a struct. Received: %s', class(config));
end

if ~isfield(config, 'working_dir')
    error('step0_sort_dicom:MissingConfig', ...
        'config.working_dir is required but not provided.');
end

if ~isfield(config, 'treatment_site') || isempty(config.treatment_site)
    config.treatment_site = 'Pancreas';
    fprintf('  [INFO] Using default treatment_site: %s\n', config.treatment_site);
end

if ~isfolder(config.working_dir)
    error('step0_sort_dicom:DirectoryNotFound', ...
        'Working directory does not exist: %s', config.working_dir);
end

%% ======================== CONSTRUCT PATHS ========================

rawwd = fullfile(config.working_dir, 'EthosExports', patient_id, ...
    config.treatment_site, session);

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

fprintf('  Classifying RTPLAN files (REFERENCE / ADAPTED)...\n');
sortRTFiles(ctInfo, rawwd, sct_dir, sctSeriesUID);

%% ======================== SORT REG (IMAGE REGISTRATION) FILES ========================

fprintf('  Sorting image registration (REG) files...\n');
sortREGFiles(ctInfo, rawwd, sct_dir);

%% ======================== SORT CBCT FILES ========================

fprintf('  Sorting CBCT series...\n');
try
    sort_CBCT(patient_id, session, config);
catch ME
    warning('step0_sort_dicom:CBCTFailed', ...
        'sort_CBCT failed: %s', ME.message);
end

%% ======================== SORT EXTRA (NON-SCT) RTSTRUCT FILES ========================

fprintf('  Sorting extra (non-SCT) RTSTRUCT files...\n');
sortExtraStructFiles(ctInfo, sct_dir, sctSeriesUID);

%% ======================== PATCH PATIENT NAMES ========================

%fprintf('  Patching PatientName to append "research" in sorted files...\n');
%appendResearchToPatientName(sct_dir);

%% ======================== VERIFY OUTPUT ========================

sctFiles = dir(fullfile(sct_dir, '*.dcm'));
fprintf('  Sorting complete. %d files in sct directory.\n', length(sctFiles));

hasCT              = ~isempty(dir(fullfile(sct_dir, 'CT*.dcm')));
nRTSTRUCT          = numel(dir(fullfile(sct_dir, 'RTSTRUCT_SCT*.dcm')));
hasRTSTRUCT        = nRTSTRUCT > 0;
nExtraStructs      = numel(dir(fullfile(sct_dir, 'RTSTRUCT_extra_*.dcm')));
hasRTPLANref       = isfile(fullfile(sct_dir, 'RTPLAN_reference.dcm'));
hasRTPLANadp       = isfile(fullfile(sct_dir, 'RTPLAN_adapted.dcm'));
hasRTDOSEref       = isfile(fullfile(sct_dir, 'RTDOSE_reference.dcm'));
hasRTDOSEadp       = isfile(fullfile(sct_dir, 'RTDOSE_adapted.dcm'));
hasREG             = ~isempty(dir(fullfile(sct_dir, 'REG_*.dcm')));
nCBCTfiles         = length(dir(fullfile(sct_dir, 'CBCT*.dcm')));
hasCBCT            = nCBCTfiles > 0;

if ~hasCT
    warning('step0_sort_dicom:MissingFile', 'No CT files found in sct directory');
end
if ~hasRTPLANref
    warning('step0_sort_dicom:MissingFile', 'No RTPLAN_reference.dcm found in sct directory');
end
if ~hasRTPLANadp
    fprintf('  [INFO] RTPLAN_adapted.dcm not found (not required for RayStation import).\n');
end
if ~hasRTSTRUCT
    warning('step0_sort_dicom:MissingFile', 'No RTSTRUCT_SCT*.dcm files found in sct directory');
end
if ~hasRTDOSEref
    warning('step0_sort_dicom:MissingFile', 'No RTDOSE_reference.dcm found in sct directory');
end
if ~hasRTDOSEadp
    fprintf('  [INFO] RTDOSE_adapted.dcm not found (not required for RayStation import).\n');
end
fprintf('\n  --- Sorted file summary ---\n');
fprintf('    CT slices            : %s\n', tf2str(hasCT));
fprintf('    RTPLAN_reference     : %s\n', tf2str(hasRTPLANref));
fprintf('    RTSTRUCT_SCT* files  : %s (%d found)\n', tf2str(hasRTSTRUCT), nRTSTRUCT);
fprintf('    RTSTRUCT_extra_*     : %d found\n', nExtraStructs);
fprintf('    RTDOSE_reference     : %s\n', tf2str(hasRTDOSEref));
fprintf('    RTPLAN_adapted(opt.) : %s\n', tf2str(hasRTPLANadp));
fprintf('    RTDOSE_adapted(opt.) : %s\n', tf2str(hasRTDOSEadp));
fprintf('    REG (image reg.)     : %s\n', tf2str(hasREG));
fprintf('    CBCT files           : %s (%d files)\n', tf2str(hasCBCT), nCBCTfiles);
fprintf('  ---------------------------\n');

fprintf('  Step 0 complete for %s/%s\n', patient_id, session);

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function s = tf2str(val)
    if val, s = 'FOUND'; else, s = 'MISSING'; end
end


function sctSeriesUID = sortSctFiles(ctInfo, sourceDir, destDir) %#ok<INUSL>
%SORTSCTFILES Sort SCT DICOM files and extract SeriesInstanceUID
%
%   sctSeriesUID = sortSctFiles(ctInfo, sourceDir, destDir)
%
%   Three-tier priority:
%     1. SeriesDescription exactly 'sct' (case-insensitive)
%     2. CT series with empty SeriesDate/SeriesTime
%     3. Oldest CT timestamp

    sctSeriesUID = '';

    if ~ismember('SeriesDescription', ctInfo.Properties.VariableNames)
        warning('sortSctFiles:NoSeriesDescription', ...
            'SeriesDescription column not found in DICOM collection');
        return;
    end

    rowIndex = strcmpi(ctInfo.SeriesDescription, 'sct');
    sctInfo  = ctInfo(rowIndex, :);

    if height(sctInfo) == 0
        warning('sortSctFiles:NoSCT', 'No SCT series found in collection');
        return;
    end

    fprintf('    Found %d SCT series\n', height(sctInfo));

    sctFiles = sctInfo.Filenames{1};

    if isempty(sctFiles)
        warning('sortSctFiles:EmptySeries', 'SCT series has no files');
        return;
    end

    firstFile = sctFiles{1};
    try
        sctMetadata  = dicominfo(firstFile);
        sctSeriesUID = sctMetadata.SeriesInstanceUID;
        fprintf('    SCT SeriesInstanceUID: %s\n', sctSeriesUID);
        fprintf('    SCT Series Date/Time : %s / %s\n', ...
            sctMetadata.SeriesDate, sctMetadata.SeriesTime);
    catch ME
        warning('sortSctFiles:MetadataError', ...
            'Failed to read SCT metadata: %s', ME.message);
        return;
    end

    numMoved   = 0;
    numSkipped = 0;

    for k = 1:length(sctFiles)
        srcFile = sctFiles{k};
        [~, name, ext] = fileparts(srcFile);
        destFile = fullfile(destDir, [name, ext]);

        if exist(srcFile, 'file')
            if exist(destFile, 'file')
                numSkipped = numSkipped + 1;
            else
                try
                    copyfile(srcFile, destFile);
                    numMoved = numMoved + 1;
                catch ME
                    warning('sortSctFiles:CopyError', ...
                        'Failed to copy file %s: %s', name, ME.message);
                end
            end
        end
    end

    fprintf('    SCT files: %d copied, %d already existed\n', numMoved, numSkipped);
end


function sortREGFiles(ctInfo, sourceDir, destDir) %#ok<INUSL>
%SORTREGFILES Copy all REG (Spatial/Image Registration) DICOM files into
%             the sct directory, with sequential REG_<n>.dcm filenames.

    if ~ismember('Modality', ctInfo.Properties.VariableNames)
        return;
    end

    regRows = strcmpi(ctInfo.Modality, 'REG');
    if ~any(regRows)
        fprintf('    No REG (image registration) files found.\n');
        return;
    end

    regTable = ctInfo(regRows, :);
    fprintf('    Found %d REG series\n', height(regTable));

    nCopied  = 0;
    nSkipped = 0;
    regIdx   = 0;

    for ri = 1:height(regTable)
        fileCell = regTable.Filenames{ri};
        if isempty(fileCell), continue; end
        for k = 1:numel(fileCell)
            srcFile = fileCell{k};
            if isempty(srcFile) || ~isfile(srcFile), continue; end

            regIdx   = regIdx + 1;
            destFile = fullfile(destDir, sprintf('REG_%03d.dcm', regIdx));

            if isfile(destFile)
                nSkipped = nSkipped + 1;
                continue;
            end

            try
                copyfile(srcFile, destFile);
                nCopied = nCopied + 1;
            catch ME
                warning('sortREGFiles:CopyError', ...
                    'Failed to copy REG file %s: %s', srcFile, ME.message);
            end
        end
    end

    fprintf('    REG files: %d copied, %d already existed\n', nCopied, nSkipped);
end


function sortRTFiles(ctInfo, sourceDir, destDir, sctSeriesUID) %#ok<INUSL>
%SORTRTFILES Classify RTPLAN files as REFERENCE or ADAPTED, then trace and
%            copy each plan's RTSTRUCT and RTDOSE.
%
%   sortRTFiles(ctInfo, sourceDir, destDir, sctSeriesUID)
%
%   For each RTPLAN, the field
%       metadata.ReferenceRTPlanSequence.Item_1.RTPlanRelationship
%   is read.  If it contains 'REFERENCE' the plan is classified as the
%   reference plan; if it contains 'ADAPTED' it is the adaptive plan.
%   The field name ReferencedRTPlanSequence is tried as a fallback.
%
%   For each classified plan the function:
%     1. Traces the RTSTRUCT via ReferencedStructureSetSequence
%     2. Copies it as RTSTRUCT_SCT.dcm (RTSTRUCT_SCT_2.dcm, ... if both
%        plans use different SCT structure sets) only when it references
%        the SCT. Non-SCT structure sets are handled by sortExtraStructFiles.
%     3. Traces the RTDOSE via its ReferencedRTPlanSequence
%     4. Copies files:  RTPLAN/RTDOSE_reference.dcm  or  ..._adapted.dcm

    % ---- collect all RTPLAN file paths --------------------------------
    if ~ismember('Modality', ctInfo.Properties.VariableNames)
        warning('sortRTFiles:NoModality', 'Modality column missing');
        return;
    end

    planRows = strcmp(ctInfo.Modality, 'RTPLAN');
    if ~any(planRows)
        warning('sortRTFiles:NoRTPLAN', 'No RTPLAN files found in collection');
        return;
    end
    planTable = ctInfo(planRows, :);
    fprintf('    Found %d RTPLAN series to inspect\n', height(planTable));

    % ---- collect all RTSTRUCT and RTDOSE paths ----------------------
    allStructPaths = collectModalityPaths(ctInfo, 'RTSTRUCT');
    allDosePaths   = collectModalityPaths(ctInfo, 'RTDOSE');
    fprintf('    Found %d RTSTRUCT, %d RTDOSE files\n', ...
        length(allStructPaths), length(allDosePaths));

    % ---- classify each RTPLAN by RTPlanRelationship -----------------
    refPlanPath = '';
    adpPlanPath = '';

    for pi = 1:height(planTable)
        fileCell = planTable.Filenames{pi};
        if isempty(fileCell) || isempty(fileCell{1}), continue; end
        planPath = fileCell{1};

        try
            meta = dicominfo(planPath);
        catch
            warning('sortRTFiles:DicomReadError', ...
                'Cannot read DICOM metadata from: %s', planPath);
            continue;
        end

        relationship = extractRTPlanRelationship(meta);

        if isempty(relationship)
            fprintf('    [SKIP] No RTPlanRelationship found in: %s\n', planPath);
            continue;
        end

        fprintf('    Plan: %s    RTPlanRelationship = "%s"\n', ...
            planPath, relationship);

        if contains(upper(relationship), 'REFERENCE')
            if ~isempty(refPlanPath)
                warning('sortRTFiles:DuplicatePlan', ...
                    'Multiple REFERENCE plans found; keeping first.');
            else
                refPlanPath = planPath;
                fprintf('       Classified as REFERENCE plan\n');
            end

        elseif contains(upper(relationship), 'ADAPTED')
            if ~isempty(adpPlanPath)
                warning('sortRTFiles:DuplicatePlan', ...
                    'Multiple ADAPTED plans found; keeping first.');
            else
                adpPlanPath = planPath;
                fprintf('       Classified as ADAPTED plan\n');
            end
        else
            fprintf('       Unrecognised relationship "%s" (skipped)\n', relationship);
        end
    end

    % ---- process each plan type ------------------------------------
    planTypes  = {'reference',  'adapted'};
    planPaths  = {refPlanPath,  adpPlanPath};

    % Count SCT RTSTRUCTs across both plan types for enumeration (_2, _3)
    nSctStructs = 0;

    for ti = 1:2
        label    = planTypes{ti};
        planPath = planPaths{ti};

        if isempty(planPath)
            warning('sortRTFiles:MissingPlan', ...
                'No %s plan found; skipping RTSTRUCT/RTPLAN/RTDOSE for this type.', upper(label));
            continue;
        end

        fprintf('\n  --- Processing %s plan ---\n', upper(label));

        % Copy RTPLAN
        destRP = fullfile(destDir, sprintf('RTPLAN_%s.dcm', label));
        copyFileAs(planPath, destRP, sprintf('RTPLAN_%s', label));

        planMeta = dicominfo(planPath);

        % ---- find and copy RTSTRUCT ---------------------------------
        structSOPUID = '';
        try
            refSSSeq = planMeta.ReferencedStructureSetSequence;
            if isfield(refSSSeq, 'Item_1') && ...
               isfield(refSSSeq.Item_1, 'ReferencedSOPInstanceUID')
                structSOPUID = refSSSeq.Item_1.ReferencedSOPInstanceUID;
            end
        catch
        end

        if isempty(structSOPUID)
            warning('sortRTFiles:NoStructRef', ...
                '%s plan has no ReferencedStructureSetSequence.', upper(label));
        else
            structPath = findBySopUID(allStructPaths, structSOPUID);
            if isempty(structPath)
                warning('sortRTFiles:StructNotFound', ...
                    'RTSTRUCT with SOPInstanceUID %s not found for %s plan.', ...
                    structSOPUID, upper(label));
            elseif isempty(sctSeriesUID) || ...
                   ~strcmp(extractReferencedCTUID(structPath), sctSeriesUID)
                fprintf('    [INFO] %s plan RTSTRUCT does not reference the SCT (goes to extras).\n', ...
                    upper(label));
            else
                nSctStructs = nSctStructs + 1;
                if nSctStructs == 1
                    ctLabel = 'SCT';
                else
                    ctLabel = sprintf('SCT_%d', nSctStructs);
                end
                destRS = fullfile(destDir, sprintf('RTSTRUCT_%s.dcm', ctLabel));
                copyFileAs(structPath, destRS, sprintf('RTSTRUCT_%s', ctLabel));
            end
        end

        % ---- find and copy RTDOSE -----------------------------------
        planSOPUID = planMeta.SOPInstanceUID;
        dosePath   = findDoseForPlan(allDosePaths, planSOPUID);

        if isempty(dosePath)
            warning('sortRTFiles:DoseNotFound', ...
                'No RTDOSE referencing %s plan (SOPInstanceUID %s).', ...
                upper(label), planSOPUID);
        else
            destRD = fullfile(destDir, sprintf('RTDOSE_%s.dcm', label));
            copyFileAs(dosePath, destRD, sprintf('RTDOSE_%s', label));
        end
    end
end


function relationship = extractRTPlanRelationship(meta)
%EXTRACTRTPLANRELATIONSHIP Read RTPlanRelationship from RTPLAN metadata.
%   Tries 'ReferenceRTPlanSequence' first, then 'ReferencedRTPlanSequence'.
%   Returns '' if not found.

    relationship = '';

    candidates = {'ReferenceRTPlanSequence', 'ReferencedRTPlanSequence'};

    for ci = 1:length(candidates)
        fieldName = candidates{ci};
        if isfield(meta, fieldName)
            seq = meta.(fieldName);
            if isstruct(seq) && isfield(seq, 'Item_1')
                item1 = seq.Item_1;
                if isfield(item1, 'RTPlanRelationship')
                    relationship = strtrim(item1.RTPlanRelationship);
                    return;
                end
            end
        end
    end
end


function paths = collectModalityPaths(ctInfo, modality)
%COLLECTMODALITYPATHS Return cell array of file paths for a given modality.

    paths = {};
    if ~ismember('Modality', ctInfo.Properties.VariableNames), return; end

    rows = strcmp(ctInfo.Modality, modality);
    if ~any(rows), return; end

    subTable = ctInfo(rows, :);
    for ri = 1:height(subTable)
        fileCell = subTable.Filenames{ri};
        if ~isempty(fileCell) && ~isempty(fileCell{1})
            paths{end+1} = fileCell{1}; %#ok<AGROW>
        end
    end
end


function matchPath = findBySopUID(filePaths, targetSOPUID)
%FINDBYSOPUID Return the first file whose SOPInstanceUID matches targetSOPUID.

    matchPath = '';
    for fi = 1:length(filePaths)
        try
            meta = dicominfo(filePaths{fi});
            if isfield(meta, 'SOPInstanceUID') && ...
               strcmp(meta.SOPInstanceUID, targetSOPUID)
                matchPath = filePaths{fi};
                return;
            end
        catch
        end
    end
end


function dosePath = findDoseForPlan(dosePaths, planSOPUID)
%FINDDOSEFORPLAN Return first RTDOSE whose ReferencedRTPlanSequence points
%               to planSOPUID.

    dosePath = '';
    for di = 1:length(dosePaths)
        try
            meta = dicominfo(dosePaths{di});
            if isfield(meta, 'ReferencedRTPlanSequence')
                seq = meta.ReferencedRTPlanSequence;
                if isstruct(seq) && isfield(seq, 'Item_1') && ...
                   isfield(seq.Item_1, 'ReferencedSOPInstanceUID')
                    if strcmp(seq.Item_1.ReferencedSOPInstanceUID, planSOPUID)
                        dosePath = dosePaths{di};
                        return;
                    end
                end
            end
        catch
        end
    end
end



function referencedUID = extractReferencedCTUID(structPath)
%EXTRACTREFERENCEDCTUID Extract the referenced CT SeriesInstanceUID from an RTSTRUCT.

    referencedUID = '';
    try
        meta = dicominfo(structPath);
        if isfield(meta, 'ReferencedFrameOfReferenceSequence')
            refFOR = meta.ReferencedFrameOfReferenceSequence;
            if isstruct(refFOR) && isfield(refFOR, 'Item_1')
                item1 = refFOR.Item_1;
                if isfield(item1, 'RTReferencedStudySequence')
                    studySeq = item1.RTReferencedStudySequence;
                    if isstruct(studySeq) && isfield(studySeq, 'Item_1')
                        studyItem = studySeq.Item_1;
                        if isfield(studyItem, 'RTReferencedSeriesSequence')
                            seriesSeq = studyItem.RTReferencedSeriesSequence;
                            if isstruct(seriesSeq) && isfield(seriesSeq, 'Item_1') && ...
                               isfield(seriesSeq.Item_1, 'SeriesInstanceUID')
                                referencedUID = seriesSeq.Item_1.SeriesInstanceUID;
                            end
                        end
                    end
                end
            end
        end
    catch ME
        warning('extractReferencedCTUID:ReadError', ...
            'Could not read CT reference from RTSTRUCT %s: %s', structPath, ME.message);
    end
end


function sortExtraStructFiles(ctInfo, sct_dir, sctSeriesUID)
%SORTEXTRASTRUCTFILES Copy every RTSTRUCT that does NOT reference the SCT
%   into sct_dir as RTSTRUCT_extra_<n>.dcm (n = order in the DICOM
%   collection). These are the CBCT structure sets; step06 moves them to
%   Raystation_Input/<pid>/<session>/extra_structs/.
%   Files already present in sct_dir are skipped.

    allStructPaths = collectModalityPaths(ctInfo, 'RTSTRUCT');
    if isempty(allStructPaths)
        fprintf('    No RTSTRUCT files found in collection.\n');
        return;
    end

    nExtra = 0;
    for fi = 1:numel(allStructPaths)
        structPath    = allStructPaths{fi};
        referencedUID = extractReferencedCTUID(structPath);

        if ~isempty(sctSeriesUID) && strcmp(referencedUID, sctSeriesUID)
            continue;   % SCT structure set, already copied by sortRTFiles
        end

        nExtra   = nExtra + 1;
        destFile = fullfile(sct_dir, sprintf('RTSTRUCT_extra_%d.dcm', nExtra));
        copyFileAs(structPath, destFile, sprintf('RTSTRUCT_extra_%d', nExtra));
    end

    fprintf('    %d non-SCT RTSTRUCT file(s) sorted as extras.\n', nExtra);
end


function copyFileAs(srcPath, destPath, label)
%COPYFILEAS Copy srcPath to destPath with logging.

    if isempty(srcPath) || ~isfile(srcPath)
        warning('copyFileAs:NotFound', 'Source not found for %s: %s', label, srcPath);
        return;
    end

    if isfile(destPath)
        fprintf('    %s already exists in destination (skipping)\n', label);
        return;
    end

    try
        copyfile(srcPath, destPath);
        fprintf('    Copied %s  %s\n', label, destPath);
    catch ME
        warning('copyFileAs:CopyError', 'Failed to copy %s: %s', label, ME.message);
    end
end



function appendResearchToPatientName(dirPath)
%APPENDRESEARCHTOPATIENTNAME Append " research" to PatientName in every .dcm
%   file in dirPath.
%
%   Reads the PatientName field (expected format: "FirstName LastName"),
%   and rewrites it as "FirstName LastName research" using dicomwrite with
%   CreateMode='copy' to preserve all other tags.  Files that already
%   contain 'research' in the patient name are skipped.

    if isempty(dirPath) || ~isfolder(dirPath)
        return;
    end

    dcmFiles = dir(fullfile(dirPath, '*.dcm'));
    if isempty(dcmFiles)
        return;
    end

    numPatched  = 0;
    numSkipped  = 0;
    numFailed   = 0;

    for k = 1:length(dcmFiles)
        fPath = fullfile(dcmFiles(k).folder, dcmFiles(k).name);

        try
            meta = dicominfo(fPath);
        catch
            numFailed = numFailed + 1;
            continue;
        end

        % Extract current patient name --------------------------------
        currentName = '';
        if isfield(meta, 'PatientName')
            pn = meta.PatientName;
            if isstruct(pn)
                % DICOM PN struct: FamilyName, GivenName, ...
                given  = '';
                family = '';
                if isfield(pn, 'GivenName'),  given  = strtrim(pn.GivenName);  end
                if isfield(pn, 'FamilyName'), family = strtrim(pn.FamilyName); end
                if ~isempty(given) && ~isempty(family)
                    currentName = sprintf('%s %s', given, family);
                elseif ~isempty(family)
                    currentName = family;
                elseif ~isempty(given)
                    currentName = given;
                end
            elseif ischar(pn) || isstring(pn)
                currentName = strtrim(char(pn));
            end
        end

        if isempty(currentName)
            numSkipped = numSkipped + 1;
            continue;
        end

        % Skip if already patched ------------------------------------
        if contains(lower(currentName), 'research')
            numSkipped = numSkipped + 1;
            continue;
        end

        newName = [currentName, ' research'];

        % Write patched file in-place --------------------------------
        try
            imgData = dicomread(fPath);
            meta.PatientName = newName;
            tmpPath = [fPath, '.tmp'];
            dicomwrite(imgData, tmpPath, meta, 'CreateMode', 'copy', ...
                'WritePrivate', true);
            movefile(tmpPath, fPath, 'f');
            numPatched = numPatched + 1;
        catch ME
            warning('appendResearchToPatientName:WriteError', ...
                'Failed to patch %s: %s', dcmFiles(k).name, ME.message);
            % Clean up temp file if it exists
            tmpPath = [fPath, '.tmp'];
            if isfile(tmpPath), delete(tmpPath); end
            numFailed = numFailed + 1;
        end
    end

    fprintf('    PatientName patched: %d updated, %d skipped, %d failed  [%s]\n', ...
        numPatched, numSkipped, numFailed, dirPath);
end


function printCollectionInfo(ctInfo)
%PRINTCOLLECTIONINFO Print diagnostic information about DICOM collection.

    fprintf('\n  --- DICOM Collection Summary ---\n');
    fprintf('  Total series: %d\n', height(ctInfo));

    if ismember('Modality', ctInfo.Properties.VariableNames)
        modalities = unique(ctInfo.Modality);
        for i = 1:length(modalities)
            count = sum(strcmp(ctInfo.Modality, modalities{i}));
            fprintf('    %s: %d series\n', modalities{i}, count);
        end
    end

    fprintf('\n  --- SCT Series ---\n');
    if ismember('SeriesDescription', ctInfo.Properties.VariableNames)
        sctRows = strcmpi(ctInfo.SeriesDescription, 'sct');
        if any(sctRows)
            sctInfo  = ctInfo(sctRows, :);
            fileCell = sctInfo.Filenames{1};
            if ~isempty(fileCell) && ~isempty(fileCell{1})
                try
                    metadata = dicominfo(fileCell{1});
                    fprintf('    Series     : %s\n', metadata.SeriesDescription);
                    fprintf('    Date/Time  : %s / %s\n', metadata.SeriesDate, metadata.SeriesTime);
                    fprintf('    SeriesUID  : %s\n', metadata.SeriesInstanceUID);
                    fprintf('    Image count: %d\n', length(fileCell));
                catch
                    fprintf('    (metadata unavailable)\n');
                end
            end
        else
            fprintf('    No SCT series found\n');
        end
    end

    fprintf('\n  --- RT File Summary ---\n');
    for mod = {'RTSTRUCT', 'RTPLAN', 'RTDOSE'}
        if ismember('Modality', ctInfo.Properties.VariableNames)
            n = sum(strcmp(ctInfo.Modality, mod{1}));
            if n > 0
                fprintf('    %s: %d found\n', mod{1}, n);
                % For RTPLANs, print the RTPlanRelationship for each
                if strcmp(mod{1}, 'RTPLAN')
                    planRows = strcmp(ctInfo.Modality, 'RTPLAN');
                    planTable = ctInfo(planRows, :);
                    for pi = 1:height(planTable)
                        fileCell = planTable.Filenames{pi};
                        if isempty(fileCell) || isempty(fileCell{1}), continue; end
                        try
                            meta = dicominfo(fileCell{1});
                            rel  = extractRTPlanRelationship(meta);
                            if isempty(rel), rel = '(none)'; end
                            fprintf('      Plan %d: RTPlanRelationship = "%s"\n', pi, rel);
                        catch
                            fprintf('      Plan %d: (unreadable)\n', pi);
                        end
                    end
                end
            end
        end
    end

    fprintf('  --------------------------------\n\n');
end