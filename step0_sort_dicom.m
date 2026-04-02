function [sct_dir, sim_ct_dir, RT_metadata] = step0_sort_dicom(patient_id, session, config)
%% STEP0_SORT_DICOM - Sort DICOM files for ETHOS pipeline
%
%   [sct_dir, sim_ct_dir, RT_metadata] = step0_sort_dicom(patient_id, session, config)
%
%   PURPOSE:
%   Organize raw ETHOS DICOM export by:
%     1. Identifying the SCT (synthetic CT) series via three-tier priority
%     2. Identifying the simulation CT (original planning CT)
%     3. Finding RTPLANs tagged as REFERENCE and ADAPTED via RTPlanRelationship
%     4. Confirming each plan's referenced image set matches the SCT
%     5. Sorting the RTSTRUCT referenced by each plan
%     6. Sorting the RTDOSE corresponding to each plan
%     7. Copying all matched files to a clean 'sct' subdirectory with
%        descriptive filenames
%     8. Collecting metadata for all RT files into RT_metadata
%
%   OUTPUTS:
%     sct_dir     - Path to 'sct' subdirectory containing sorted files:
%                     CT*.dcm          (SCT slices)
%                     RP_reference.dcm (REFERENCE plan)
%                     RP_adapted.dcm   (ADAPTED plan)
%                     RS_reference.dcm (RTSTRUCT for REFERENCE plan)
%                     RS_adapted.dcm   (RTSTRUCT for ADAPTED plan)
%                     RD_reference.dcm (RTDOSE for REFERENCE plan)
%                     RD_adapted.dcm   (RTDOSE for ADAPTED plan)
%     sim_ct_dir  - Path to 'sim_ct' subdirectory containing simulation CT slices
%     RT_metadata - Cell array of structs, one per RT file discovered across
%                   all modalities (RTPLAN, RTSTRUCT, RTDOSE). Indexed by an
%                   external counter that increments continuously across all
%                   modality loops. Each struct contains fields:
%                     .filepath           Full path to source file
%                     .filename           Filename only
%                     .modality           'RTPLAN' | 'RTSTRUCT' | 'RTDOSE'
%                     .sorted_label       Label assigned (e.g. 'REFERENCE',
%                                         'ADAPTED', 'unmatched')
%                     .sorted_filename    Destination filename in sct dir,
%                                         or '' if not sorted
%                     .SOPInstanceUID     Unique ID for this object
%                     .SOPClassUID        Class UID (if present)
%                     .StudyInstanceUID   Study UID
%                     .SeriesInstanceUID  Series UID
%                     .SeriesDate         DICOM series date string
%                     .SeriesTime         DICOM series time string
%                   RTPLAN-specific:
%                     .RTPlanLabel        Plan label
%                     .RTPlanName         Plan name
%                     .RTPlanRelationship Relationship string (REFERENCE/ADAPTED/etc)
%                     .ReferencedStructSetUID  SOPInstanceUID of linked RTSTRUCT
%                   RTSTRUCT-specific:
%                     .StructureSetLabel  Structure set label
%                     .StructureSetDate   Structure set date
%                     .ReferencedFrameUID FrameOfReferenceUID (if present)
%                   RTDOSE-specific:
%                     .DoseType           PHYSICAL | EFFECTIVE | ERROR
%                     .DoseUnits          GY | RELATIVE
%                     .DoseSummationType  PLAN | FRACTION | BEAM | etc
%                     .ReferencedPlanUID  SOPInstanceUID of linked RTPLAN
%                   Also saved as RT_metadata.mat in the sct directory.
%
%   SCT IDENTIFICATION (three-tier priority):
%     1. SeriesDescription contains 'sct' (case-insensitive)
%     2. Empty or unparseable SeriesDate/SeriesTime (ETHOS often exports sCT
%        without timestamps)
%     3. Oldest valid timestamp among remaining CT series
%
%   SIM CT IDENTIFICATION (three-tier priority):
%     1. SeriesDescription contains 'sim' (case-insensitive)
%     2. Empty or unparseable SeriesDate/SeriesTime (if not already used for SCT)
%     3. Oldest valid timestamp among remaining CT series

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('step0_sort_dicom:InvalidInput', ...
        'patient_id must be a char or string. Received: %s', class(patient_id));
end
patient_id = char(patient_id);

if ~ischar(session) && ~isstring(session)
    error('step0_sort_dicom:InvalidInput', ...
        'session must be a char or string. Received: %s', class(session));
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
sct_dir    = fullfile(rawwd, 'sct');
sim_ct_dir = '';
RT_metadata  = {};
rt_meta_idx  = 0;   % External counter — increments across ALL modality loops

fprintf('\n========== STEP 0: SORT DICOM ==========\n');
fprintf('  Patient:       %s\n', patient_id);
fprintf('  Session:       %s\n', session);
fprintf('  Raw directory: %s\n', rawwd);

%% ======================== VERIFY RAW DIRECTORY ========================

if ~isfolder(rawwd)
    warning('step0_sort_dicom:DirectoryNotFound', ...
        'Raw directory not found: %s', rawwd);
    sct_dir    = '';
    sim_ct_dir = '';
    RT_metadata = {};
    return;
end

%% ======================== SCAN DICOM COLLECTION ========================

fprintf('\n  [1/6] Scanning DICOM collection...\n');

try
    dcmCol = dicomCollection(rawwd);
catch ME
    error('step0_sort_dicom:DicomScanFailed', ...
        'Failed to scan DICOM directory: %s\nError: %s', rawwd, ME.message);
end

if isempty(dcmCol) || height(dcmCol) == 0
    warning('step0_sort_dicom:EmptyCollection', ...
        'No DICOM files found in: %s', rawwd);
    sct_dir    = '';
    sim_ct_dir = '';
    RT_metadata = {};
    return;
end

fprintf('    Found %d DICOM series\n', height(dcmCol));
printCollectionInfo(dcmCol);

%% ======================== CREATE OUTPUT DIRECTORIES ========================

if ~isfolder(sct_dir)
    mkdir(sct_dir);
    fprintf('\n  Created sct directory: %s\n', sct_dir);
else
    fprintf('\n  sct directory exists: %s\n', sct_dir);
end

sim_ct_dir = fullfile(rawwd, 'sim_ct');
if ~isfolder(sim_ct_dir)
    mkdir(sim_ct_dir);
end

%% ======================== SORT SCT FILES ========================

fprintf('\n  [2/6] Identifying and sorting SCT (synthetic CT)...\n');
[sctSeriesUID, sctStudyUID, sctFrameUID] = sortSctFiles(dcmCol, sct_dir);

if isempty(sctSeriesUID)
    warning('step0_sort_dicom:NoSCT', ...
        'No SCT series found for patient %s, %s', patient_id, session);
end

%% ======================== SORT SIM CT FILES ========================

fprintf('\n  [3/6] Identifying and sorting simulation CT...\n');
sim_ct_dir = sortSimCtFiles(dcmCol, sim_ct_dir, sctSeriesUID);

%% ======================== FIND ALL RTPLAN FILES ========================

fprintf('\n  [4/6] Scanning RTPLAN files for RTPlanRelationship...\n');
allRtplanFiles = findAllRtplanFiles(rawwd);

if isempty(allRtplanFiles)
    warning('step0_sort_dicom:NoRTPlans', 'No RTPLAN files found in %s', rawwd);
    return;
end

fprintf('    Found %d RTPLAN file(s) total\n', length(allRtplanFiles));

% Search each RTPLAN for RTPlanRelationship = contains REFERENCE or ADAPTED
% Also collect metadata for every RTPLAN encountered
referencePlanFile = '';
adaptedPlanFile   = '';
referencePlanUID  = '';
adaptedPlanUID    = '';
referenceStructUID = '';
adaptedStructUID   = '';

for i = 1:length(allRtplanFiles)
    fp = allRtplanFiles{i};
    try
        info = dicominfo(fp);
    catch
        fprintf('    [WARN] Could not read: %s\n', fp);
        continue;
    end

    relationship = getRTPlanRelationship(info);
    fprintf('    %s  -> RTPlanRelationship: "%s"\n', ...
        getFriendlyName(fp), relationship);

    % Determine sorted label for this plan
    sorted_label    = 'unmatched';
    sorted_filename = '';

    % Check for REFERENCE
    if contains(relationship, 'REFERENCE', 'IgnoreCase', true) && isempty(referencePlanFile)
        referencePlanFile  = fp;
        referencePlanUID   = info.SOPInstanceUID;
        referenceStructUID = getReferencedStructUID(info);
        sorted_label       = 'REFERENCE';
        sorted_filename    = 'RP_reference.dcm';
        fprintf('      ** Identified as REFERENCE plan\n');
        confirmImageSet(info, sctSeriesUID, sctStudyUID, sctFrameUID, 'REFERENCE');
    end

    % Check for ADAPTED
    if contains(relationship, 'ADAPTED', 'IgnoreCase', true) && isempty(adaptedPlanFile)
        adaptedPlanFile  = fp;
        adaptedPlanUID   = info.SOPInstanceUID;
        adaptedStructUID = getReferencedStructUID(info);
        sorted_label     = 'ADAPTED';
        sorted_filename  = 'RP_adapted.dcm';
        fprintf('      ** Identified as ADAPTED plan\n');
        confirmImageSet(info, sctSeriesUID, sctStudyUID, sctFrameUID, 'ADAPTED');
    end

    % Collect metadata entry
    rt_meta_idx = rt_meta_idx + 1;
    RT_metadata{rt_meta_idx} = collectRTMetadata(info, fp, 'RTPLAN', ...
        sorted_label, sorted_filename);
end

% Report
if isempty(referencePlanFile)
    warning('step0_sort_dicom:NoReferencePlan', ...
        'No RTPLAN with RTPlanRelationship containing "REFERENCE" found.');
end
if isempty(adaptedPlanFile)
    warning('step0_sort_dicom:NoAdaptedPlan', ...
        'No RTPLAN with RTPlanRelationship containing "ADAPTED" found.');
end

%% ======================== COPY RTPLAN FILES ========================

fprintf('\n  [5/6] Copying RTPLAN files to sct directory...\n');

if ~isempty(referencePlanFile)
    destRef = fullfile(sct_dir, 'RP_reference.dcm');
    copyDicomFile(referencePlanFile, destRef, 'RP_reference.dcm');
end

if ~isempty(adaptedPlanFile)
    destAdp = fullfile(sct_dir, 'RP_adapted.dcm');
    copyDicomFile(adaptedPlanFile, destAdp, 'RP_adapted.dcm');
end

%% ======================== SORT RTSTRUCT FILES ========================

fprintf('\n  [6a/6] Sorting RTSTRUCT files...\n');
allRtstructFiles = findAllRtFiles(rawwd, 'RS');

sortRtstructForPlan(allRtstructFiles, referenceStructUID, sct_dir, 'RS_reference.dcm', 'REFERENCE');
sortRtstructForPlan(allRtstructFiles, adaptedStructUID,   sct_dir, 'RS_adapted.dcm',   'ADAPTED');

% Collect metadata for all RTSTRUCT files
for i = 1:length(allRtstructFiles)
    fp = allRtstructFiles{i};
    try
        info = dicominfo(fp);
    catch
        fprintf('    [WARN] Could not read RTSTRUCT metadata: %s\n', getFriendlyName(fp));
        continue;
    end
    % Determine label
    sorted_label    = 'unmatched';
    sorted_filename = '';
    if isfield(info, 'SOPInstanceUID')
        uid = strtrim(char(info.SOPInstanceUID));
        if strcmp(uid, referenceStructUID)
            sorted_label    = 'REFERENCE';
            sorted_filename = 'RS_reference.dcm';
        elseif strcmp(uid, adaptedStructUID)
            sorted_label    = 'ADAPTED';
            sorted_filename = 'RS_adapted.dcm';
        end
    end
    rt_meta_idx = rt_meta_idx + 1;
    RT_metadata{rt_meta_idx} = collectRTMetadata(info, fp, 'RTSTRUCT', ...
        sorted_label, sorted_filename);
end

%% ======================== SORT RTDOSE FILES ========================

fprintf('\n  [6b/6] Sorting RTDOSE files...\n');
allRtdoseFiles = findAllRtFiles(rawwd, 'RD');

sortRtdoseForPlan(allRtdoseFiles, referencePlanUID, sct_dir, 'RD_reference.dcm', 'REFERENCE');
sortRtdoseForPlan(allRtdoseFiles, adaptedPlanUID,   sct_dir, 'RD_adapted.dcm',   'ADAPTED');

% Collect metadata for all RTDOSE files
for i = 1:length(allRtdoseFiles)
    fp = allRtdoseFiles{i};
    try
        info = dicominfo(fp);
    catch
        fprintf('    [WARN] Could not read RTDOSE metadata: %s\n', getFriendlyName(fp));
        continue;
    end
    % Determine label by matching referenced plan UID
    sorted_label    = 'unmatched';
    sorted_filename = '';
    refPlanUID = getReferencedPlanUID(info);
    if strcmp(refPlanUID, referencePlanUID) && ~isempty(referencePlanUID)
        sorted_label    = 'REFERENCE';
        sorted_filename = 'RD_reference.dcm';
    elseif strcmp(refPlanUID, adaptedPlanUID) && ~isempty(adaptedPlanUID)
        sorted_label    = 'ADAPTED';
        sorted_filename = 'RD_adapted.dcm';
    end
    rt_meta_idx = rt_meta_idx + 1;
    RT_metadata{rt_meta_idx} = collectRTMetadata(info, fp, 'RTDOSE', ...
        sorted_label, sorted_filename);
end

%% ======================== VERIFY OUTPUT ========================

fprintf('\n  === Sorting complete ===\n');
allOut = dir(fullfile(sct_dir, '*.dcm'));
ctCount  = sum(~cellfun(@isempty, regexp({allOut.name}, '^CT', 'once')));
rtCount  = length(allOut) - ctCount;
fprintf('  SCT directory: %d CT slices, %d RT files\n', ctCount, rtCount);
fprintf('  Expected RT files: RP_reference, RP_adapted, RS_reference, RS_adapted, RD_reference, RD_adapted\n');

% List RT files found
for k = 1:length(allOut)
    if ~startsWith(allOut(k).name, 'CT')
        fprintf('    %s\n', allOut(k).name);
    end
end

simOut = dir(fullfile(sim_ct_dir, '*.dcm'));
fprintf('  Sim CT directory: %d CT slices\n', length(simOut));

%% ======================== SAVE RT METADATA ========================

fprintf('\n  RT_metadata: %d RT file entries collected\n', rt_meta_idx);
fprintf('  %-12s  %-12s  %-30s  %s\n', 'Modality', 'Label', 'Sorted filename', 'Source file');
fprintf('  %s\n', repmat('-', 1, 90));
for k = 1:rt_meta_idx
    m = RT_metadata{k};
    fprintf('  %-12s  %-12s  %-30s  %s\n', ...
        m.modality, m.sorted_label, m.sorted_filename, m.filename);
end

% Save to sct directory for downstream reference
rt_meta_save_path = fullfile(sct_dir, 'RT_metadata.mat');
save(rt_meta_save_path, 'RT_metadata');
fprintf('\n  RT_metadata saved to: %s\n', rt_meta_save_path);

end % END MAIN FUNCTION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printCollectionInfo(dcmCol)
%PRINTCOLLECTIONINFO Display summary of DICOM collection
    fprintf('    %-40s %-15s %-20s %-20s\n', ...
        'SeriesDescription', 'Modality', 'SeriesDate', 'SeriesTime');
    fprintf('    %s\n', repmat('-', 1, 100));
    for i = 1:height(dcmCol)
        desc = '';
        mod  = '';
        dat  = '';
        tim  = '';
        if ismember('SeriesDescription', dcmCol.Properties.VariableNames)
            desc = char(dcmCol.SeriesDescription(i));
        end
        if ismember('Modality', dcmCol.Properties.VariableNames)
            mod = char(dcmCol.Modality(i));
        end
        if ismember('SeriesDate', dcmCol.Properties.VariableNames)
            dat = char(dcmCol.SeriesDate(i));
        end
        if ismember('SeriesTime', dcmCol.Properties.VariableNames)
            tim = char(dcmCol.SeriesTime(i));
        end
        fprintf('    %-40s %-15s %-20s %-20s\n', ...
            desc(1:min(end,40)), mod(1:min(end,15)), dat(1:min(end,20)), tim(1:min(end,20)));
    end
end


function [sctSeriesUID, sctStudyUID, sctFrameUID] = sortSctFiles(dcmCol, destDir)
%SORTSCTFILES Identify SCT series and copy CT slices to destDir
%   Three-tier priority:
%     1. SeriesDescription contains 'sct'
%     2. Empty/unparseable SeriesDate or SeriesTime
%     3. Oldest valid timestamp

    sctSeriesUID = '';
    sctStudyUID  = '';
    sctFrameUID  = '';

    % Restrict to CT modality rows
    if ismember('Modality', dcmCol.Properties.VariableNames)
        ctMask = strcmpi(dcmCol.Modality, 'CT');
        ctCol  = dcmCol(ctMask, :);
    else
        ctCol = dcmCol;
    end

    if height(ctCol) == 0
        warning('sortSctFiles:NoCT', 'No CT series found in collection');
        return;
    end

    % --- Tier 1: description contains 'sct' ---
    selectedRow = [];
    if ismember('SeriesDescription', ctCol.Properties.VariableNames)
        descMatch = find(contains(ctCol.SeriesDescription, 'sct', 'IgnoreCase', true), 1);
        if ~isempty(descMatch)
            selectedRow = descMatch;
            fprintf('    SCT identified by description match (Tier 1)\n');
        end
    end

    % --- Tier 2: empty datetime ---
    if isempty(selectedRow) && ismember('SeriesDate', ctCol.Properties.VariableNames)
        for i = 1:height(ctCol)
            d = char(ctCol.SeriesDate(i));
            t = char(ctCol.SeriesTime(i));
            if isempty(strtrim(d)) || isempty(strtrim(t))
                selectedRow = i;
                fprintf('    SCT identified by empty datetime (Tier 2)\n');
                break;
            end
        end
    end

    % --- Tier 3: oldest valid timestamp ---
    if isempty(selectedRow)
        oldestTime = Inf;
        for i = 1:height(ctCol)
            dt = parseDateTime(char(ctCol.SeriesDate(i)), char(ctCol.SeriesTime(i)));
            if ~isnan(dt) && dt < oldestTime
                oldestTime   = dt;
                selectedRow  = i;
            end
        end
        if ~isempty(selectedRow)
            fprintf('    SCT identified by oldest timestamp (Tier 3)\n');
        end
    end

    if isempty(selectedRow)
        warning('sortSctFiles:NotFound', 'Could not identify SCT series');
        return;
    end

    % Read first file for UIDs
    files = ctCol.Filenames{selectedRow};
    if isempty(files)
        warning('sortSctFiles:EmptySeries', 'SCT series has no files');
        return;
    end

    try
        firstInfo    = dicominfo(files{1});
        sctSeriesUID = firstInfo.SeriesInstanceUID;
        sctStudyUID  = firstInfo.StudyInstanceUID;
        if isfield(firstInfo, 'FrameOfReferenceUID')
            sctFrameUID = firstInfo.FrameOfReferenceUID;
        end
        fprintf('    SCT SeriesUID: %s\n', sctSeriesUID);
        fprintf('    SCT StudyUID:  %s\n', sctStudyUID);
        if ~isempty(sctFrameUID)
            fprintf('    SCT FrameUID:  %s\n', sctFrameUID);
        end
    catch ME
        warning('sortSctFiles:MetadataError', 'Failed to read SCT metadata: %s', ME.message);
        return;
    end

    % Copy CT slices
    numCopied  = 0;
    numSkipped = 0;
    for k = 1:length(files)
        src = files{k};
        [~, nm, ex] = fileparts(src);
        dst = fullfile(destDir, [nm ex]);
        if exist(dst, 'file')
            numSkipped = numSkipped + 1;
        else
            try
                copyfile(src, dst);
                numCopied = numCopied + 1;
            catch ME
                warning('sortSctFiles:CopyError', 'Failed to copy %s: %s', nm, ME.message);
            end
        end
    end
    fprintf('    SCT: %d slices copied, %d already existed\n', numCopied, numSkipped);
end


function sim_ct_dir = sortSimCtFiles(dcmCol, destDir, excludeSeriesUID)
%SORTSIMCTFILES Identify simulation CT and copy slices to destDir
%   Three-tier priority (excludes the SCT series):
%     1. SeriesDescription contains 'sim'
%     2. Empty/unparseable SeriesDate or SeriesTime
%     3. Oldest valid timestamp

    sim_ct_dir = destDir;

    if ismember('Modality', dcmCol.Properties.VariableNames)
        ctMask = strcmpi(dcmCol.Modality, 'CT');
        ctCol  = dcmCol(ctMask, :);
    else
        ctCol = dcmCol;
    end

    if height(ctCol) == 0
        warning('sortSimCtFiles:NoCT', 'No CT series found');
        return;
    end

    % Exclude the already-identified SCT series
    if ~isempty(excludeSeriesUID) && ismember('SeriesInstanceUID', ctCol.Properties.VariableNames)
        notSct = ~strcmp(ctCol.SeriesInstanceUID, excludeSeriesUID);
        ctCol  = ctCol(notSct, :);
    end

    if height(ctCol) == 0
        fprintf('    No additional CT series found for sim CT\n');
        return;
    end

    selectedRow = [];

    % Tier 1: description contains 'sim'
    if ismember('SeriesDescription', ctCol.Properties.VariableNames)
        descMatch = find(contains(ctCol.SeriesDescription, 'sim', 'IgnoreCase', true), 1);
        if ~isempty(descMatch)
            selectedRow = descMatch;
            fprintf('    Sim CT identified by description match (Tier 1)\n');
        end
    end

    % Tier 2: empty datetime
    if isempty(selectedRow) && ismember('SeriesDate', ctCol.Properties.VariableNames)
        for i = 1:height(ctCol)
            d = char(ctCol.SeriesDate(i));
            t = char(ctCol.SeriesTime(i));
            if isempty(strtrim(d)) || isempty(strtrim(t))
                selectedRow = i;
                fprintf('    Sim CT identified by empty datetime (Tier 2)\n');
                break;
            end
        end
    end

    % Tier 3: oldest valid timestamp
    if isempty(selectedRow)
        oldestTime = Inf;
        for i = 1:height(ctCol)
            dt = parseDateTime(char(ctCol.SeriesDate(i)), char(ctCol.SeriesTime(i)));
            if ~isnan(dt) && dt < oldestTime
                oldestTime  = dt;
                selectedRow = i;
            end
        end
        if ~isempty(selectedRow)
            fprintf('    Sim CT identified by oldest timestamp (Tier 3)\n');
        end
    end

    if isempty(selectedRow)
        fprintf('    [WARN] Could not identify simulation CT series\n');
        return;
    end

    files = ctCol.Filenames{selectedRow};
    numCopied  = 0;
    numSkipped = 0;
    for k = 1:length(files)
        src = files{k};
        [~, nm, ex] = fileparts(src);
        dst = fullfile(destDir, [nm ex]);
        if exist(dst, 'file')
            numSkipped = numSkipped + 1;
        else
            try
                copyfile(src, dst);
                numCopied = numCopied + 1;
            catch ME
                warning('sortSimCtFiles:CopyError', 'Failed to copy %s: %s', nm, ME.message);
            end
        end
    end
    fprintf('    Sim CT: %d slices copied, %d already existed\n', numCopied, numSkipped);
end


function rtplanFiles = findAllRtplanFiles(searchDir)
%FINDALLRTPLANFILES Recursively find all RTPLAN DICOM files
%   Identifies by Modality = RTPLAN (reads each .dcm header)

    rtplanFiles = {};
    allDcm = dir(fullfile(searchDir, '**', '*.dcm'));

    for i = 1:length(allDcm)
        fp = fullfile(allDcm(i).folder, allDcm(i).name);
        % Quick check via filename prefix
        if startsWith(allDcm(i).name, 'RP')
            rtplanFiles{end+1} = fp; %#ok<AGROW>
            continue;
        end
        % Otherwise check modality tag
        try
            info = dicominfo(fp, 'UseDictionaryVR', false);
            if isfield(info, 'Modality') && strcmpi(info.Modality, 'RTPLAN')
                rtplanFiles{end+1} = fp; %#ok<AGROW>
            end
        catch
            % Not a readable DICOM or not RTPLAN
        end
    end
end


function rtFiles = findAllRtFiles(searchDir, prefix)
%FINDALLRTFILES Find all RT files with given filename prefix (RS or RD)

    rtFiles = {};
    allDcm = dir(fullfile(searchDir, '**', '*.dcm'));
    for i = 1:length(allDcm)
        if startsWith(allDcm(i).name, prefix)
            rtFiles{end+1} = fullfile(allDcm(i).folder, allDcm(i).name); %#ok<AGROW>
        end
    end
    % Also catch files that match by modality if none found by prefix
    if isempty(rtFiles)
        modTarget = '';
        if strcmp(prefix, 'RS'), modTarget = 'RTSTRUCT'; end
        if strcmp(prefix, 'RD'), modTarget = 'RTDOSE';   end
        if ~isempty(modTarget)
            for i = 1:length(allDcm)
                fp = fullfile(allDcm(i).folder, allDcm(i).name);
                try
                    info = dicominfo(fp, 'UseDictionaryVR', false);
                    if isfield(info, 'Modality') && strcmpi(info.Modality, modTarget)
                        rtFiles{end+1} = fp; %#ok<AGROW>
                    end
                catch
                end
            end
        end
    end
end


function relationship = getRTPlanRelationship(planInfo)
%GETRTPLANRELATIONSHIP Extract RTPlanRelationship from plan DICOM info
%   Checks ReferenceRTPlanSequence.Item_1.RTPlanRelationship
%   Returns empty string if not found

    relationship = '';

    if ~isfield(planInfo, 'ReferenceRTPlanSequence')
        return;
    end

    seq = planInfo.ReferenceRTPlanSequence;

    % Try Item_1 first
    if isfield(seq, 'Item_1')
        item = seq.Item_1;
    elseif isstruct(seq) && ~isempty(fieldnames(seq))
        fields = fieldnames(seq);
        item = seq.(fields{1});
    else
        return;
    end

    if isfield(item, 'RTPlanRelationship')
        relationship = strtrim(char(item.RTPlanRelationship));
    end
end


function structUID = getReferencedStructUID(planInfo)
%GETREFERENCEDSTRUCTUID Get the SOPInstanceUID of the referenced RTSTRUCT
    structUID = '';
    if isfield(planInfo, 'ReferencedStructureSetSequence')
        seq = planInfo.ReferencedStructureSetSequence;
        if isfield(seq, 'Item_1') && isfield(seq.Item_1, 'ReferencedSOPInstanceUID')
            structUID = strtrim(char(seq.Item_1.ReferencedSOPInstanceUID));
        end
    end
end


function confirmImageSet(planInfo, sctSeriesUID, sctStudyUID, sctFrameUID, planLabel)
%CONFIRMIMAGESETMATCH Check that plan's referenced image set matches the SCT
%   Checks StudyInstanceUID and/or FrameOfReferenceUID linkage

    matched = false;
    details = '';

    % Check ReferencedStudySequence
    if isfield(planInfo, 'ReferencedStudySequence')
        seq = planInfo.ReferencedStudySequence;
        if isfield(seq, 'Item_1') && isfield(seq.Item_1, 'ReferencedSOPInstanceUID')
            planStudy = strtrim(char(seq.Item_1.ReferencedSOPInstanceUID));
            if strcmp(planStudy, sctStudyUID)
                matched = true;
                details = 'StudyInstanceUID matches SCT';
            end
        end
    end

    % Check StudyInstanceUID directly on plan
    if ~matched && isfield(planInfo, 'StudyInstanceUID')
        if strcmp(strtrim(char(planInfo.StudyInstanceUID)), sctStudyUID)
            matched = true;
            details = 'Plan StudyUID matches SCT StudyUID';
        end
    end

    % Check FrameOfReferenceUID linkage
    if ~matched && ~isempty(sctFrameUID) && isfield(planInfo, 'FrameOfReferenceUID')
        if strcmp(strtrim(char(planInfo.FrameOfReferenceUID)), sctFrameUID)
            matched = true;
            details = 'FrameOfReferenceUID matches SCT';
        end
    end

    if matched
        fprintf('      [OK] %s plan image set confirmed: %s\n', planLabel, details);
    else
        warning('step0_sort_dicom:ImageSetMismatch', ...
            '%s plan image set could NOT be confirmed to match SCT. ' ...
            'Verify manually that plan references the correct SCT.', planLabel);
    end
end


function sortRtstructForPlan(rtstructFiles, targetStructUID, destDir, outFilename, label)
%SORTRSTRUCTFORPLAN Find the RTSTRUCT with matching SOPInstanceUID and copy it

    if isempty(targetStructUID)
        fprintf('    [WARN] No referenced RTSTRUCT UID for %s plan\n', label);
        return;
    end

    for i = 1:length(rtstructFiles)
        fp = rtstructFiles{i};
        try
            info = dicominfo(fp);
        catch
            continue;
        end
        if isfield(info, 'SOPInstanceUID') && ...
                strcmp(strtrim(char(info.SOPInstanceUID)), targetStructUID)
            copyDicomFile(fp, fullfile(destDir, outFilename), outFilename);
            fprintf('    [OK] %s RTSTRUCT -> %s\n', label, outFilename);
            return;
        end
    end

    warning('step0_sort_dicom:RtstructNotFound', ...
        'RTSTRUCT with SOPInstanceUID %s (for %s plan) not found.', ...
        targetStructUID, label);
end


function sortRtdoseForPlan(rtdoseFiles, planSOPUID, destDir, outFilename, label)
%SORTRTSOSEFORPLAN Find the RTDOSE referencing the given plan SOPInstanceUID

    if isempty(planSOPUID)
        fprintf('    [WARN] No plan SOPInstanceUID for %s RTDOSE search\n', label);
        return;
    end

    for i = 1:length(rtdoseFiles)
        fp = rtdoseFiles{i};
        try
            info = dicominfo(fp);
        catch
            continue;
        end
        referencedUID = getReferencedPlanUID(info);
        if strcmp(referencedUID, planSOPUID)
            copyDicomFile(fp, fullfile(destDir, outFilename), outFilename);
            fprintf('    [OK] %s RTDOSE -> %s\n', label, outFilename);
            return;
        end
    end

    warning('step0_sort_dicom:RtdoseNotFound', ...
        'RTDOSE referencing plan %s (for %s) not found.', planSOPUID, label);
end


function planUID = getReferencedPlanUID(doseInfo)
%GETREFERENCEDPLANUID Extract the referenced RTPLAN SOPInstanceUID from RTDOSE
    planUID = '';
    if isfield(doseInfo, 'ReferencedRTPlanSequence')
        seq = doseInfo.ReferencedRTPlanSequence;
        if isfield(seq, 'Item_1') && isfield(seq.Item_1, 'ReferencedSOPInstanceUID')
            planUID = strtrim(char(seq.Item_1.ReferencedSOPInstanceUID));
        end
    end
end


function copyDicomFile(src, dst, label)
%COPYDICOMFILE Copy a single DICOM file with feedback
    if exist(dst, 'file')
        fprintf('    %s already exists, skipping copy\n', label);
        return;
    end
    try
        copyfile(src, dst);
        fprintf('    Copied: %s\n', label);
    catch ME
        warning('step0_sort_dicom:CopyError', 'Failed to copy %s: %s', label, ME.message);
    end
end


function name = getFriendlyName(filepath)
%GETFRIENDLYNAME Return just the filename from a full path
    [~, nm, ex] = fileparts(filepath);
    name = [nm ex];
end


function meta = collectRTMetadata(info, filepath, modality, sorted_label, sorted_filename)
%COLLECTRTMETADATA Extract key metadata fields from a DICOM info struct
%
%   Builds a flat struct with fields common to all RT modalities plus
%   modality-specific fields. Missing DICOM fields are silently set to ''.

    [~, nm, ex] = fileparts(filepath);
    meta.filepath        = filepath;
    meta.filename        = [nm ex];
    meta.modality        = modality;
    meta.sorted_label    = sorted_label;
    meta.sorted_filename = sorted_filename;

    % --- Common fields ---
    meta.SOPInstanceUID      = getField(info, 'SOPInstanceUID');
    meta.SOPClassUID         = getField(info, 'SOPClassUID');
    meta.StudyInstanceUID    = getField(info, 'StudyInstanceUID');
    meta.SeriesInstanceUID   = getField(info, 'SeriesInstanceUID');
    meta.SeriesDate          = getField(info, 'SeriesDate');
    meta.SeriesTime          = getField(info, 'SeriesTime');
    meta.SeriesDescription   = getField(info, 'SeriesDescription');
    meta.PatientID           = getField(info, 'PatientID');
    meta.FrameOfReferenceUID = getField(info, 'FrameOfReferenceUID');

    % --- Modality-specific fields ---
    switch upper(modality)

        case 'RTPLAN'
            meta.RTPlanLabel    = getField(info, 'RTPlanLabel');
            meta.RTPlanName     = getField(info, 'RTPlanName');
            meta.RTPlanDate     = getField(info, 'RTPlanDate');
            meta.RTPlanTime     = getField(info, 'RTPlanTime');
            meta.RTPlanGeometry = getField(info, 'RTPlanGeometry');
            meta.RTPlanRelationship     = getRTPlanRelationship(info);
            meta.ReferencedStructSetUID = '';
            if isfield(info, 'ReferencedStructureSetSequence')
                seq = info.ReferencedStructureSetSequence;
                if isfield(seq, 'Item_1') && isfield(seq.Item_1, 'ReferencedSOPInstanceUID')
                    meta.ReferencedStructSetUID = ...
                        strtrim(char(seq.Item_1.ReferencedSOPInstanceUID));
                end
            end
            meta.NumberOfBeams = 0;
            if isfield(info, 'BeamSequence')
                meta.NumberOfBeams = length(fieldnames(info.BeamSequence));
            end
            meta.NumberOfFractionsPlanned = '';
            if isfield(info, 'FractionGroupSequence') && ...
               isfield(info.FractionGroupSequence, 'Item_1') && ...
               isfield(info.FractionGroupSequence.Item_1, 'NumberOfFractionsPlanned')
                meta.NumberOfFractionsPlanned = ...
                    info.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
            end

        case 'RTSTRUCT'
            meta.StructureSetLabel = getField(info, 'StructureSetLabel');
            meta.StructureSetName  = getField(info, 'StructureSetName');
            meta.StructureSetDate  = getField(info, 'StructureSetDate');
            meta.StructureSetTime  = getField(info, 'StructureSetTime');
            meta.ReferencedSeriesUID = '';
            if isfield(info, 'ReferencedFrameOfReferenceSequence')
                seq = info.ReferencedFrameOfReferenceSequence;
                if isfield(seq, 'Item_1') && ...
                   isfield(seq.Item_1, 'RTReferencedStudySequence')
                    studySeq = seq.Item_1.RTReferencedStudySequence;
                    if isfield(studySeq, 'Item_1') && ...
                       isfield(studySeq.Item_1, 'RTReferencedSeriesSequence')
                        seriesSeq = studySeq.Item_1.RTReferencedSeriesSequence;
                        if isfield(seriesSeq, 'Item_1') && ...
                           isfield(seriesSeq.Item_1, 'SeriesInstanceUID')
                            meta.ReferencedSeriesUID = ...
                                strtrim(char(seriesSeq.Item_1.SeriesInstanceUID));
                        end
                    end
                end
            end
            meta.NumberOfROIs = 0;
            if isfield(info, 'StructureSetROISequence')
                meta.NumberOfROIs = length(fieldnames(info.StructureSetROISequence));
            end

        case 'RTDOSE'
            meta.DoseType          = getField(info, 'DoseType');
            meta.DoseUnits         = getField(info, 'DoseUnits');
            meta.DoseSummationType = getField(info, 'DoseSummationType');
            meta.DoseGridScaling   = '';
            if isfield(info, 'DoseGridScaling')
                meta.DoseGridScaling = info.DoseGridScaling;
            end
            meta.ReferencedPlanUID = getReferencedPlanUID(info);
            meta.Rows              = getField(info, 'Rows');
            meta.Columns           = getField(info, 'Columns');
            meta.NumberOfFrames    = getField(info, 'NumberOfFrames');
    end
end


function val = getField(info, fieldName)
%GETFIELD Safely extract a string or numeric field from a dicominfo struct
%   Returns '' if the field does not exist
    if isfield(info, fieldName)
        raw = info.(fieldName);
        if ischar(raw) || isstring(raw)
            val = strtrim(char(raw));
        else
            val = raw;   % Return numeric/logical as-is
        end
    else
        val = '';
    end
end


function dt = parseDateTime(dateStr, timeStr)
%PARSEDATETIME Convert DICOM date/time strings to numeric value for comparison
%   Returns NaN if either string is empty or unparseable

    dt = NaN;
    dateStr = strtrim(char(dateStr));
    timeStr = strtrim(char(timeStr));

    if isempty(dateStr) || isempty(timeStr)
        return;
    end

    % DICOM date: YYYYMMDD, time: HHMMSS or HHMMSS.FFFFFF
    try
        if length(dateStr) >= 8 && length(timeStr) >= 6
            yr  = str2double(dateStr(1:4));
            mo  = str2double(dateStr(5:6));
            dy  = str2double(dateStr(7:8));
            hr  = str2double(timeStr(1:2));
            mn  = str2double(timeStr(3:4));
            sc  = str2double(timeStr(5:6));
            if ~any(isnan([yr mo dy hr mn sc]))
                dt = datenum(yr, mo, dy, hr, mn, sc);
            end
        end
    catch
        % Leave as NaN
    end
end