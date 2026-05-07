% diagnostic_rtstruct.m
%
% PURPOSE:
%   Searches the EthosExports directory for a given patient and session,
%   finds all RTSTRUCT DICOM files (searching recursively), and prints
%   every ROI contour found in each file.
%
% USAGE:
%   Edit patient_id, session, and working_dir below, then run (F5 or Run).
%   Results are printed to the console.
%
% DEPENDENCIES:
%   MATLAB Image Processing Toolbox (dicominfo)
%
% See also: dicominfo, step0_sort_dicom, diagnostic_deformable_registration

% ------------------------------------------------------------------ %
%  CONFIGURATION — edit before running                                %
% ------------------------------------------------------------------ %
patient_id    = '1194203';
session       = 'Session_5';
working_dir   = 'C:/Users/80030361/ETHOS_Simulations';
treatment_site = 'Pancreas';

% ------------------------------------------------------------------ %
%  DERIVED PATHS                                                       %
% ------------------------------------------------------------------ %
search_root = fullfile(working_dir, 'EthosExports', patient_id, treatment_site, session);

% ------------------------------------------------------------------ %
%  HEADER                                                              %
% ------------------------------------------------------------------ %
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  RTSTRUCT ROI CONTOUR DIAGNOSTIC\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('Patient     : %s\n', patient_id);
fprintf('Session     : %s\n', session);
fprintf('Search root : %s\n', search_root);
fprintf('Timestamp   : %s\n\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

if ~isfolder(search_root)
    error('diagnostic_rtstruct:NotFound', ...
          'Directory not found: %s', search_root);
end

% ------------------------------------------------------------------ %
%  STEP 1: RECURSIVE FILE DISCOVERY                                   %
% ------------------------------------------------------------------ %
fprintf('[STEP 1] Scanning for DICOM files...\n');
allFiles = dir(fullfile(search_root, '**', '*'));
allFiles = allFiles(~[allFiles.isdir]);
fprintf('  Found %d files total.\n\n', numel(allFiles));

% ------------------------------------------------------------------ %
%  STEP 2: IDENTIFY RTSTRUCT FILES                                    %
% ------------------------------------------------------------------ %
fprintf('[STEP 2] Identifying RTSTRUCT files...\n');

rtstructPaths = {};

for fi = 1:numel(allFiles)
    fPath = fullfile(allFiles(fi).folder, allFiles(fi).name);

    try
        meta = dicominfo(fPath);
    catch
        continue;  % Not a DICOM file
    end

    if isfield(meta, 'Modality') && strcmpi(strtrim(meta.Modality), 'RTSTRUCT')
        rtstructPaths{end+1} = fPath; %#ok<AGROW>
    end
end

fprintf('  Found %d RTSTRUCT file(s).\n\n', numel(rtstructPaths));

if isempty(rtstructPaths)
    fprintf('  No RTSTRUCT files found. Exiting.\n');
    return;
end

% ------------------------------------------------------------------ %
%  STEP 3: PRINT ROI CONTOURS FOR EACH RTSTRUCT                      %
% ------------------------------------------------------------------ %
fprintf('[STEP 3] Printing ROI contours for each RTSTRUCT...\n');

for fi = 1:numel(rtstructPaths)
    fPath = rtstructPaths{fi};
    [~, fname, fext] = fileparts(fPath);

    fprintf('\n%s\n', repmat('-', 1, 70));
    fprintf('  File %d/%d: %s%s\n', fi, numel(rtstructPaths), fname, fext);
    fprintf('  Full path : %s\n', fPath);

    try
        meta = dicominfo(fPath);
    catch ME
        fprintf('  [ERROR] Could not read DICOM metadata: %s\n', ME.message);
        continue;
    end

    % -- Print StructureSet-level identifiers -------------------------
    fprintf('\n  --- Structure Set Info ---\n');
    printField(meta, 'SOPInstanceUID',      'SOP Instance UID');
    printField(meta, 'StructureSetLabel',   'Label           ');
    printField(meta, 'StructureSetName',    'Name            ');
    printField(meta, 'StructureSetDate',    'Date            ');
    printField(meta, 'StructureSetTime',    'Time            ');

    % -- Collect contour counts from ROIContourSequence ---------------
    contourCounts = buildContourCountMap(meta);

    % -- Print each ROI from StructureSetROISequence ------------------
    fprintf('\n  --- ROI List (%s) ---\n', getROICountStr(meta));

    if ~isfield(meta, 'StructureSetROISequence')
        fprintf('    [WARN] No StructureSetROISequence found in this file.\n');
        continue;
    end

    roiSeq   = meta.StructureSetROISequence;
    itemKeys = fieldnames(roiSeq);

    for ri = 1:numel(itemKeys)
        item = roiSeq.(itemKeys{ri});

        roiNum  = getItemField(item, 'ROINumber',   '?');
        roiName = getItemField(item, 'ROIName',     '(unnamed)');
        roiDesc = getItemField(item, 'ROIDescription', '');
        roiAlgo = getItemField(item, 'ROIGenerationAlgorithm', '');

        % Look up contour slice count for this ROI
        if isnumeric(roiNum)
            roiNumKey = roiNum;
        else
            roiNumKey = str2double(roiNum);
        end
        nSlices = getContourCount(contourCounts, roiNumKey);

        % Format output line
        descStr = '';
        if ~isempty(roiDesc),  descStr = sprintf('  [desc: %s]', roiDesc); end
        algoStr = '';
        if ~isempty(roiAlgo) && ~strcmpi(roiAlgo, 'NONE')
            algoStr = sprintf('  [algo: %s]', roiAlgo);
        end

        fprintf('    ROI %3s | %-40s | %3d slice(s)%s%s\n', ...
            num2str(roiNumKey), roiName, nSlices, descStr, algoStr);
    end
end

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  Diagnostic complete.\n');
fprintf('%s\n\n', repmat('=', 1, 70));


% ================================================================== %
%  LOCAL HELPER FUNCTIONS                                             %
% ================================================================== %

function printField(meta, fieldName, label)
%PRINTFIELD Print a scalar metadata field if it exists.
    if isfield(meta, fieldName)
        val = meta.(fieldName);
        if ischar(val) || isstring(val)
            fprintf('  %s : %s\n', label, strtrim(char(val)));
        elseif isnumeric(val)
            fprintf('  %s : %g\n', label, val);
        end
    end
end


function s = getItemField(item, fieldName, defaultVal)
%GETITEMFIELD Safely retrieve a field from a sequence item struct.
    if isfield(item, fieldName)
        val = item.(fieldName);
        if ischar(val) || isstring(val)
            s = strtrim(char(val));
        elseif isnumeric(val)
            s = val;
        else
            s = defaultVal;
        end
    else
        s = defaultVal;
    end
end


function countMap = buildContourCountMap(meta)
%BUILDCONTOURCOUNTMAP Build a map from ROINumber -> number of contour slices.
%   Iterates over ROIContourSequence items and counts ContourSequence items.
    countMap = containers.Map('KeyType', 'double', 'ValueType', 'double');

    if ~isfield(meta, 'ROIContourSequence')
        return;
    end

    roiContourSeq = meta.ROIContourSequence;
    itemKeys      = fieldnames(roiContourSeq);

    for ki = 1:numel(itemKeys)
        item = roiContourSeq.(itemKeys{ki});

        if ~isfield(item, 'ReferencedROINumber')
            continue;
        end

        refROINum = item.ReferencedROINumber;
        if ~isnumeric(refROINum)
            refROINum = str2double(refROINum);
        end
        if isnan(refROINum), continue; end

        nSlices = 0;
        if isfield(item, 'ContourSequence')
            contSeq = item.ContourSequence;
            nSlices = numel(fieldnames(contSeq));
        end

        countMap(refROINum) = nSlices;
    end
end


function n = getContourCount(countMap, roiNum)
%GETCONTOURCOUNT Return contour slice count for a given ROI number.
    if isnan(roiNum) || ~isKey(countMap, roiNum)
        n = 0;
    else
        n = countMap(roiNum);
    end
end


function s = getROICountStr(meta)
%GETROICOUNTSTR Return a string describing the number of ROIs.
    if ~isfield(meta, 'StructureSetROISequence')
        s = 'no ROIs';
        return;
    end
    n = numel(fieldnames(meta.StructureSetROISequence));
    s = sprintf('%d ROIs', n);
end
