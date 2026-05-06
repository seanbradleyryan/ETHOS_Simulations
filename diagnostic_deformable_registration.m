% diagnostic_deformable_registration.m
%
% PURPOSE:
%   Searches the ETHOS export directory tree for Deformable Spatial
%   Registration DICOM files (SOP class 1.2.840.10008.5.1.4.1.1.66.3) and
%   Spatial Registration files (SOP class 1.2.840.10008.5.1.4.1.1.66.1),
%   then reports key metadata and performs targeted checks against the
%   RayStation import error:
%
%     "Registration Type Code Sequence must contain exactly one item for
%      the registered item in Deformable Spatial Registration <UID>"
%
% USAGE:
%   Edit search_root below, then run this script directly (F5 or Run).
%   Results are printed to the console. The variable 'diag_regFiles' is
%   left in the workspace for interactive inspection afterward.
%
% DEPENDENCIES:
%   MATLAB Image Processing Toolbox (dicominfo)
%
% See also: dicominfo, step0_sort_dicom

% ------------------------------------------------------------------ %
%  CONFIGURATION — edit search_root before running                    %
% ------------------------------------------------------------------ %
search_root = 'C:/Users/seanr/ETHOS_Simulations/EthosExports';

SOP_DEFORMABLE = '1.2.840.10008.5.1.4.1.1.66.3';
SOP_SPATIAL    = '1.2.840.10008.5.1.4.1.1.66.1';

% ------------------------------------------------------------------ %
%  CONFIG — diagnostic options                                         %
% ------------------------------------------------------------------ %
CONFIG.remove_duplicates    = false;  % Remove DeformableRegistrationSequence items
    % that lack PreDeformationMatrixRegistrationSequence before inspection.

CONFIG.compare_registrations = true;  % Compare found deformable REG file(s) against
    % RayStationTestRegistration.dcm in the path below.
CONFIG.compare_ref_path = 'C:/Users/seanr/ETHOS_Simulations/Raystation_Input/1194203/Session_1/RayStationTestRegistration.dcm';

if ~isfolder(search_root)
    error('diagnostic_deformable_registration:NotFound', ...
          'Directory not found: %s', search_root);
end

% ------------------------------------------------------------------ %
%  STEP 1: RECURSIVE FILE DISCOVERY                                   %
% ------------------------------------------------------------------ %
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEFORMABLE REGISTRATION DICOM DIAGNOSTIC\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('Search root : %s\n', search_root);
fprintf('Timestamp   : %s\n\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

fprintf('[STEP 1] Scanning for DICOM files...\n');
allFiles = dir(fullfile(search_root, '**', '*'));
allFiles = allFiles(~[allFiles.isdir]);
fprintf('  Found %d files total.\n\n', numel(allFiles));

% ------------------------------------------------------------------ %
%  STEP 2: FILTER TO REGISTRATION DICOM FILES                         %
% ------------------------------------------------------------------ %
fprintf('[STEP 2] Filtering to DEFORMABLE registration DICOM files...\n');
fprintf('  (Files without DeformableRegistrationSequence are skipped.)\n');

regFiles     = struct('path', {}, 'info', {}, 'sop_class', {});
skippedCount = 0;

for k = 1:numel(allFiles)
    fpath = fullfile(allFiles(k).folder, allFiles(k).name);

    % Skip files with non-DICOM extensions (but allow files with no extension)
    [~, ~, ext] = fileparts(fpath);
    if ~isempty(ext) && ~strcmpi(ext, '.dcm')
        skippedCount = skippedCount + 1;
        continue;
    end

    try
        info = dicominfo(fpath);
    catch
        skippedCount = skippedCount + 1;
        continue;
    end

    % Only keep files that contain DeformableRegistrationSequence
    if ~isfield(info, 'DeformableRegistrationSequence')
        skippedCount = skippedCount + 1;
        continue;
    end

    sopClass = '';
    if isfield(info, 'SOPClassUID')
        sopClass = info.SOPClassUID;
    end

    regFiles(end+1).path      = fpath;      %#ok<SAGROW>
    regFiles(end).info        = info;
    regFiles(end).sop_class   = sopClass;
end

fprintf('  Skipped %d non-DICOM / non-deformable / unreadable files.\n', skippedCount);
fprintf('  Found %d deformable registration DICOM file(s).\n\n', numel(regFiles));

if isempty(regFiles)
    fprintf('[WARNING] No deformable or spatial registration DICOM files found.\n');
    fprintf('  Check that the ETHOS export includes REG files and that\n');
    fprintf('  search_root points to the correct export directory.\n\n');
    return;
end

% diag_regFiles is assigned after STEP 3 so it reflects any duplicate removal.

% ------------------------------------------------------------------ %
%  STEP 3: PER-FILE HEADER REPORT & DIAGNOSTICS                       %
% ------------------------------------------------------------------ %
totalErrors   = 0;
totalWarnings = 0;

for f = 1:numel(regFiles)
    info   = regFiles(f).info;
    fpath  = regFiles(f).path;
    sopUID = regFiles(f).sop_class;

    fprintf('%s\n', repmat('-', 1, 70));
    fprintf('FILE %d/%d: %s\n', f, numel(regFiles), fpath);
    fprintf('%s\n', repmat('-', 1, 70));

    % 3a. DICOM content summary
    fprintf('\n  [DICOM CONTENT SUMMARY]\n');
    summarize_dicom_contents(info);

    % 3b. Header metadata
    fprintf('\n  [HEADER METADATA]\n');
    print_field(info, 'SOPClassUID',           '  SOP Class UID      ');
    print_field(info, 'SOPInstanceUID',         '  SOP Instance UID   ');
    print_field(info, 'Modality',               '  Modality           ');
    print_field(info, 'Manufacturer',           '  Manufacturer       ');
    print_field(info, 'ManufacturerModelName',  '  Device Model       ');
    print_field(info, 'SoftwareVersions',       '  Software Version   ');
    print_field(info, 'StudyDate',              '  Study Date         ');
    print_field(info, 'ContentDate',            '  Content Date       ');
    print_field(info, 'FrameOfReferenceUID',    '  Frame of Ref UID   ');
    print_field(info, 'StudyInstanceUID',       '  Study UID          ');
    print_field(info, 'SeriesInstanceUID',      '  Series UID         ');
    print_field(info, 'ContentLabel',           '  Content Label      ');
    print_field(info, 'ContentDescription',     '  Content Desc       ');

    fprintf('\n  SOP Class type: ');
    if strcmp(sopUID, SOP_DEFORMABLE)
        fprintf('Deformable Spatial Registration (1.2.840.10008.5.1.4.1.1.66.3)\n');
    elseif strcmp(sopUID, SOP_SPATIAL)
        fprintf('Spatial Registration (Rigid) (1.2.840.10008.5.1.4.1.1.66.1)\n');
    else
        fprintf('Other/Unknown: %s\n', sopUID);
    end

    % 3c. Remove duplicate items (items without PreDeformationMatrix) if configured
    if CONFIG.remove_duplicates
        [info, nRemoved] = remove_non_predeform_items(info);
        if nRemoved > 0
            fprintf('\n  [REMOVE_DUPLICATES] Removed %d item(s) lacking PreDeformationMatrixRegistrationSequence.\n', nRemoved);
            regFiles(f).info = info;
        end

        % Save modified DICOM to same directory as source file
        outPath = fullfile(fileparts(fpath), 'Registration_Deformable_RSimport.dcm');
        try
            dicomwrite([], outPath, info, 'CreateMode', 'copy');
            fprintf('  [SAVED] %s\n', outPath);
        catch ME
            fprintf('  [ERROR] Could not save modified DICOM: %s\n', ME.message);
        end
    end

    % 3d. DeformableRegistrationSequence inspection
    fprintf('\n  [DEFORMABLE REGISTRATION SEQUENCE INSPECTION]\n');
    [nErr, nWarn] = inspect_deformable_sequence(info);
    totalErrors   = totalErrors   + nErr;
    totalWarnings = totalWarnings + nWarn;

    % 3e. Referenced image series
    fprintf('\n  [REFERENCED IMAGE SERIES]\n');
    inspect_referenced_series(info);

    fprintf('\n');
end

diag_regFiles = regFiles;
fprintf('  -> Variable ''diag_regFiles'' saved to workspace (%d entries).\n\n', numel(regFiles));

% ------------------------------------------------------------------ %
%  STEP 4: SUMMARY                                                     %
% ------------------------------------------------------------------ %
fprintf('%s\n', repmat('=', 1, 70));
fprintf('  DIAGNOSTIC SUMMARY\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('  REG files found : %d\n', numel(regFiles));
fprintf('  Errors flagged  : %d\n', totalErrors);
fprintf('  Warnings flagged: %d\n', totalWarnings);

if totalErrors == 0 && totalWarnings == 0
    fprintf('\n  [PASS] No conformance issues detected.\n');
    fprintf('  If RayStation still rejects the file, the issue may be in\n');
    fprintf('  RayStation version compatibility or a non-standard Elekta extension.\n');
else
    fprintf('\n  [ACTION REQUIRED] Review ERRORs above.\n');
    fprintf('  The most common fix for the RayStation import error is:\n');
    fprintf('    - RegistrationTypeCodeSequence must have exactly 1 item\n');
    fprintf('    - Check that each DeformableRegistrationSequence item has\n');
    fprintf('      one and only one RegistrationTypeCodeSequence entry.\n');
end
fprintf('%s\n\n', repmat('=', 1, 70));

% ------------------------------------------------------------------ %
%  STEP 5: COMPARE AGAINST RAYSTATION TEST REGISTRATION               %
% ------------------------------------------------------------------ %
if CONFIG.compare_registrations
    fprintf('%s\n', repmat('=', 1, 70));
    fprintf('  REGISTRATION COMPARISON\n');
    fprintf('%s\n', repmat('=', 1, 70));
    fprintf('  Reference file: %s\n\n', CONFIG.compare_ref_path);

    if ~isfile(CONFIG.compare_ref_path)
        fprintf('  [ERROR] Reference file not found: %s\n\n', CONFIG.compare_ref_path);
    else
        try
            refInfo = dicominfo(CONFIG.compare_ref_path);
        catch ME
            fprintf('  [ERROR] Could not read reference file: %s\n\n', ME.message);
            refInfo = [];
        end

        if ~isempty(refInfo)
            for f = 1:numel(regFiles)
                fprintf('%s\n', repmat('-', 1, 70));
                fprintf('  Comparing FILE %d: %s\n', f, regFiles(f).path);
                fprintf('%s\n\n', repmat('-', 1, 70));
                compare_dicom_files(regFiles(f).info, refInfo);
            end
        end
    end
    fprintf('%s\n\n', repmat('=', 1, 70));
end


% ======================================================================= %
%  LOCAL FUNCTIONS                                                         %
% ======================================================================= %

function compare_dicom_files(infoA, infoB)
    % compare_dicom_files  Recursively compare two dicominfo structs and report
    %   fields that are only in one, or present in both with different values.
    %   infoA = source (found) file; infoB = reference (RayStation test) file.
    fprintf('  Legend:  [A only]  field present only in source\n');
    fprintf('           [B only]  field present only in reference\n');
    fprintf('           [DIFF]    field present in both, values differ\n');
    fprintf('           [MATCH]   field present in both, values identical\n\n');

    nDiff   = 0;
    nAonly  = 0;
    nBonly  = 0;
    nMatch  = 0;
    compare_structs(infoA, infoB, '', nDiff, nAonly, nBonly, nMatch);
end


function compare_structs(A, B, prefix, ~, ~, ~, ~)
    % compare_structs  Recurse into two structs and print differences.
    %   Uses fprintf directly; counters are informational totals printed here.
    fieldsA = fieldnames(A);
    fieldsB = fieldnames(B);
    allFields = union(fieldsA, fieldsB);

    nDiff  = 0;
    nAonly = 0;
    nBonly = 0;
    nMatch = 0;

    for k = 1:numel(allFields)
        fname  = allFields{k};
        fqName = [prefix fname];
        inA    = isfield(A, fname);
        inB    = isfield(B, fname);

        if inA && ~inB
            nAonly = nAonly + 1;
            fprintf('  [A only]  %s\n', fqName);

        elseif ~inA && inB
            nBonly = nBonly + 1;
            fprintf('  [B only]  %s\n', fqName);

        else
            valA = A.(fname);
            valB = B.(fname);

            if isstruct(valA) && isstruct(valB)
                % Recurse one level, labelling sub-fields with dot notation
                compare_structs(valA, valB, [fqName '.'], 0, 0, 0, 0);

            elseif isstruct(valA) || isstruct(valB)
                % Type mismatch
                nDiff = nDiff + 1;
                fprintf('  [DIFF]    %s  (type mismatch: struct vs non-struct)\n', fqName);

            elseif isnumeric(valA) && isnumeric(valB)
                if isequal(size(valA), size(valB)) && all(abs(valA(:) - valB(:)) < 1e-9)
                    nMatch = nMatch + 1;
                    % Uncomment the line below to show matching numeric fields:
                    % fprintf('  [MATCH]   %s\n', fqName);
                else
                    nDiff = nDiff + 1;
                    szA = sprintf('%dx', size(valA)); szA(end) = [];
                    szB = sprintf('%dx', size(valB)); szB(end) = [];
                    if isscalar(valA) && isscalar(valB)
                        fprintf('  [DIFF]    %-55s  A=%g  B=%g\n', fqName, valA, valB);
                    elseif isequal(size(valA), size(valB))
                        maxDiff = max(abs(valA(:) - valB(:)));
                        fprintf('  [DIFF]    %-55s  [%s]  max|A-B|=%g\n', fqName, szA, maxDiff);
                    else
                        fprintf('  [DIFF]    %-55s  size A=[%s]  size B=[%s]\n', fqName, szA, szB);
                    end
                end

            elseif ischar(valA) && ischar(valB)
                if strcmp(valA, valB)
                    nMatch = nMatch + 1;
                    % Uncomment to show matching char fields:
                    % fprintf('  [MATCH]   %s = "%s"\n', fqName, valA);
                else
                    nDiff = nDiff + 1;
                    fprintf('  [DIFF]    %-55s  A="%s"  B="%s"\n', fqName, valA, valB);
                end

            elseif iscell(valA) && iscell(valB)
                strA = strjoin(cellfun(@num2str, valA(:)', 'UniformOutput', false), ',');
                strB = strjoin(cellfun(@num2str, valB(:)', 'UniformOutput', false), ',');
                if strcmp(strA, strB)
                    nMatch = nMatch + 1;
                else
                    nDiff = nDiff + 1;
                    fprintf('  [DIFF]    %-55s  A={%s}  B={%s}\n', fqName, strA, strB);
                end

            else
                % Fallback: convert to string and compare
                try
                    strA = num2str(valA);
                    strB = num2str(valB);
                    if strcmp(strA, strB)
                        nMatch = nMatch + 1;
                    else
                        nDiff = nDiff + 1;
                        fprintf('  [DIFF]    %-55s  A=%s  B=%s\n', fqName, strA, strB);
                    end
                catch
                    nDiff = nDiff + 1;
                    fprintf('  [DIFF]    %-55s  (could not compare type: %s)\n', fqName, class(valA));
                end
            end
        end
    end

    % Print totals only at the top level (no prefix)
    if isempty(prefix)
        fprintf('\n  Comparison complete: %d differ, %d A-only, %d B-only, %d match.\n\n', ...
            nDiff, nAonly, nBonly, nMatch);
    end
end


function summarize_dicom_contents(info)
    % summarize_dicom_contents  Print a concise one-line-per-field summary of
    %   all top-level fields in a dicominfo struct.
    fields = fieldnames(info);
    fprintf('  %d top-level fields:\n', numel(fields));
    for k = 1:numel(fields)
        fname = fields{k};
        val   = info.(fname);
        if isstruct(val)
            subFields = fieldnames(val);
            % Count Item_N fields as items
            nItems = sum(~cellfun(@isempty, regexp(subFields, '^Item_\d+$')));
            if nItems > 0
                fprintf('    %-45s: struct sequence (%d item(s))\n', fname, nItems);
            else
                fprintf('    %-45s: struct (%d field(s))\n', fname, numel(subFields));
            end
        elseif ischar(val)
            if numel(val) > 60
                fprintf('    %-45s: char  "%s..."\n', fname, val(1:57));
            else
                fprintf('    %-45s: char  "%s"\n', fname, val);
            end
        elseif isnumeric(val)
            sz = size(val);
            szStr = sprintf('%dx', sz);
            szStr(end) = [];  % remove trailing 'x'
            if isscalar(val)
                fprintf('    %-45s: %-10s = %g\n', fname, class(val), val);
            else
                fprintf('    %-45s: %-10s [%s]\n', fname, class(val), szStr);
            end
        elseif iscell(val)
            fprintf('    %-45s: cell  {%s}\n', fname, strjoin(val, ', '));
        else
            fprintf('    %-45s: %s\n', fname, class(val));
        end
    end
end


function [info, nRemoved] = remove_non_predeform_items(info)
    % remove_non_predeform_items  Remove items from DeformableRegistrationSequence
    %   that do not have a PreDeformationMatrixRegistrationSequence field.
    %
    % Returns the modified info struct and the count of removed items.
    nRemoved = 0;

    if ~isfield(info, 'DeformableRegistrationSequence')
        return;
    end

    seqData   = info.DeformableRegistrationSequence;
    itemNames = fieldnames(seqData);
    newSeq    = struct();
    newIdx    = 0;

    for i = 1:numel(itemNames)
        iName = itemNames{i};
        % Only process Item_N fields
        if isempty(regexp(iName, '^Item_\d+$', 'once'))
            newSeq.(iName) = seqData.(iName);
            continue;
        end
        item = seqData.(iName);
        if isfield(item, 'PreDeformationMatrixRegistrationSequence')
            newIdx = newIdx + 1;
            newSeq.(sprintf('Item_%d', newIdx)) = item;
        else
            nRemoved = nRemoved + 1;
            fprintf('    Removing %s (no PreDeformationMatrixRegistrationSequence)\n', iName);
        end
    end

    info.DeformableRegistrationSequence = newSeq;
end


function [nErrors, nWarnings] = inspect_deformable_sequence(info)
    % inspect_deformable_sequence  Validate DeformableRegistrationSequence
    %
    % Checks each item for:
    %   - RegistrationTypeCodeSequence with exactly 1 item (RayStation requirement)
    %   - SourceFrameOfReferenceUID presence
    %   - DeformableRegistrationGridSequence presence and content
    nErrors   = 0;
    nWarnings = 0;

    if ~isfield(info, 'DeformableRegistrationSequence')
        fprintf('  [ERROR] DeformableRegistrationSequence is MISSING.\n');
        nErrors = nErrors + 1;
        return;
    end

    seqData   = info.DeformableRegistrationSequence;
    itemNames = fieldnames(seqData);
    nItems    = numel(itemNames);
    fprintf('  DeformableRegistrationSequence: %d item(s)\n', nItems);

    for i = 1:nItems
        item = seqData.(itemNames{i});
        fprintf('\n    --- DeformableRegistrationSequence Item %d ---\n', i);

        % SourceFrameOfReferenceUID
        if isfield(item, 'SourceFrameOfReferenceUID')
            fprintf('    SourceFrameOfReferenceUID  : %s\n', item.SourceFrameOfReferenceUID);
        else
            fprintf('    [WARNING] SourceFrameOfReferenceUID : MISSING\n');
            nWarnings = nWarnings + 1;
        end

        % ---- RegistrationTypeCodeSequence (THE KEY CHECK) ----
        fprintf('\n    RegistrationTypeCodeSequence:\n');
        if ~isfield(item, 'RegistrationTypeCodeSequence')
            fprintf('      [ERROR] MISSING — RayStation requires exactly 1 item.\n');
            nErrors = nErrors + 1;
        else
            rtcSeq    = item.RegistrationTypeCodeSequence;
            rtcFields = fieldnames(rtcSeq);
            nRtc      = numel(rtcFields);
            fprintf('      Item count: %d  ', nRtc);

            if nRtc == 1
                fprintf('[OK]\n');
            elseif nRtc == 0
                fprintf('[ERROR — 0 items, must be exactly 1]\n');
                nErrors = nErrors + 1;
            else
                fprintf('[ERROR — %d items, must be exactly 1]\n', nRtc);
                nErrors = nErrors + 1;
            end

            for j = 1:nRtc
                rtcItem = rtcSeq.(rtcFields{j});
                fprintf('      Item %d:\n', j);
                print_seq_field(rtcItem, 'CodeValue',              '        CodeValue            ');
                print_seq_field(rtcItem, 'CodingSchemeDesignator', '        CodingSchemeDesignator');
                print_seq_field(rtcItem, 'CodingSchemeVersion',    '        CodingSchemeVersion   ');
                print_seq_field(rtcItem, 'CodeMeaning',            '        CodeMeaning          ');
                if ~isfield(rtcItem, 'CodeValue')
                    fprintf('        [WARNING] CodeValue missing.\n');
                    nWarnings = nWarnings + 1;
                end
                if ~isfield(rtcItem, 'CodingSchemeDesignator')
                    fprintf('        [WARNING] CodingSchemeDesignator missing.\n');
                    nWarnings = nWarnings + 1;
                end
                if ~isfield(rtcItem, 'CodeMeaning')
                    fprintf('        [WARNING] CodeMeaning missing.\n');
                    nWarnings = nWarnings + 1;
                end
            end
        end

        % ---- DeformableRegistrationGridSequence ----
        fprintf('\n    DeformableRegistrationGridSequence:\n');
        if ~isfield(item, 'DeformableRegistrationGridSequence')
            fprintf('      [ERROR] MISSING — no deformation grid in this item.\n');
            nErrors = nErrors + 1;
        else
            gridSeq   = item.DeformableRegistrationGridSequence;
            gridItems = fieldnames(gridSeq);
            nGrid     = numel(gridItems);
            fprintf('      Present — %d grid item(s).\n', nGrid);

            for g = 1:nGrid
                gItem = gridSeq.(gridItems{g});
                fprintf('      Grid item %d:\n', g);
                if isfield(gItem, 'Rows') && isfield(gItem, 'Columns')
                    nSlices = 0;
                    if isfield(gItem, 'GridFrameOffsetVector')
                        nSlices = numel(gItem.GridFrameOffsetVector);
                    end
                    fprintf('        Dimensions : %d x %d x %d (row x col x slice)\n', ...
                        gItem.Rows, gItem.Columns, nSlices);
                end
                if isfield(gItem, 'PixelSpacing')
                    fprintf('        PixelSpacing: [%.3f, %.3f] mm\n', ...
                        gItem.PixelSpacing(1), gItem.PixelSpacing(2));
                end
                if isfield(gItem, 'GridFrameOffsetVector')
                    ov = gItem.GridFrameOffsetVector;
                    fprintf('        Slice offsets: %d entries, [%.2f to %.2f] mm\n', ...
                        numel(ov), min(ov), max(ov));
                end
                if isfield(gItem, 'ImagePositionPatient')
                    fprintf('        Origin (mm): [%.2f, %.2f, %.2f]\n', ...
                        gItem.ImagePositionPatient(1), ...
                        gItem.ImagePositionPatient(2), ...
                        gItem.ImagePositionPatient(3));
                end
                if ~isfield(gItem, 'VectorGridData')
                    fprintf('        [WARNING] VectorGridData tag absent.\n');
                    nWarnings = nWarnings + 1;
                else
                    fprintf('        VectorGridData: present (%d bytes)\n', ...
                        numel(gItem.VectorGridData));
                end
            end
        end

        % ---- Optional pre/post rigid components ----
        if isfield(item, 'PreDeformationMatrixRegistrationSequence')
            fprintf('\n    PreDeformationMatrixRegistrationSequence: present\n');
            preSeq = item.PreDeformationMatrixRegistrationSequence;
            preF   = fieldnames(preSeq);
            for p = 1:numel(preF)
                pItem = preSeq.(preF{p});
                if isfield(pItem, 'FrameOfReferenceTransformationMatrix')
                    M = reshape(pItem.FrameOfReferenceTransformationMatrix, 4, 4)';
                    fprintf('      Matrix item %d (4x4):\n', p);
                    for row = 1:4
                        fprintf('        [%9.5f %9.5f %9.5f %9.5f ]\n', M(row,:));
                    end
                end
            end
        end

        if isfield(item, 'PostDeformationMatrixRegistrationSequence')
            fprintf('\n    PostDeformationMatrixRegistrationSequence: present\n');
        end
    end
end


function [nErrors, nWarnings] = inspect_spatial_sequence(info)
    % inspect_spatial_sequence  Validate RegistrationSequence in rigid REG files
    nErrors   = 0;
    nWarnings = 0;

    if ~isfield(info, 'RegistrationSequence')
        fprintf('  [ERROR] RegistrationSequence is MISSING.\n');
        nErrors = nErrors + 1;
        return;
    end

    seqData   = info.RegistrationSequence;
    itemNames = fieldnames(seqData);
    nItems    = numel(itemNames);
    fprintf('  RegistrationSequence: %d item(s)\n', nItems);

    for i = 1:nItems
        item = seqData.(itemNames{i});
        fprintf('\n    --- RegistrationSequence Item %d ---\n', i);

        if isfield(item, 'FrameOfReferenceUID')
            fprintf('    FrameOfReferenceUID: %s\n', item.FrameOfReferenceUID);
        end

        fprintf('\n    RegistrationTypeCodeSequence:\n');
        if ~isfield(item, 'RegistrationTypeCodeSequence')
            fprintf('      [WARNING] MISSING (may be optional for spatial reg).\n');
            nWarnings = nWarnings + 1;
        else
            rtcSeq = item.RegistrationTypeCodeSequence;
            rtcF   = fieldnames(rtcSeq);
            nRtc   = numel(rtcF);
            fprintf('      Item count: %d', nRtc);
            if nRtc == 1
                fprintf(' [OK]\n');
            else
                fprintf(' [WARNING — expected 1]\n');
                nWarnings = nWarnings + 1;
            end
            for j = 1:nRtc
                rtcItem = rtcSeq.(rtcF{j});
                print_seq_field(rtcItem, 'CodeValue',   '        CodeValue   ');
                print_seq_field(rtcItem, 'CodeMeaning', '        CodeMeaning ');
            end
        end

        if isfield(item, 'MatrixRegistrationSequence')
            matSeq = item.MatrixRegistrationSequence;
            matF   = fieldnames(matSeq);
            fprintf('\n    MatrixRegistrationSequence: %d item(s)\n', numel(matF));
            for m = 1:numel(matF)
                mItem = matSeq.(matF{m});
                if isfield(mItem, 'FrameOfReferenceTransformationMatrix')
                    M = reshape(mItem.FrameOfReferenceTransformationMatrix, 4, 4)';
                    fprintf('      Matrix item %d (4x4):\n', m);
                    for row = 1:4
                        fprintf('        [%9.5f %9.5f %9.5f %9.5f ]\n', M(row,:));
                    end
                end
                if isfield(mItem, 'FrameOfReferenceTransformationMatrixType')
                    fprintf('      Transform type: %s\n', ...
                        mItem.FrameOfReferenceTransformationMatrixType);
                end
            end
        end
    end
end


function inspect_referenced_series(info)
    % inspect_referenced_series  List referenced studies/series in the REG file
    if isfield(info, 'StudiesContainingOtherReferencedInstancesSequence')
        studySeq = info.StudiesContainingOtherReferencedInstancesSequence;
        studyF   = fieldnames(studySeq);
        fprintf('  StudiesContainingOtherReferencedInstancesSequence: %d study/studies\n', numel(studyF));
        for s = 1:numel(studyF)
            sItem = studySeq.(studyF{s});
            if isfield(sItem, 'StudyInstanceUID')
                fprintf('    Study %d UID: %s\n', s, sItem.StudyInstanceUID);
            end
            if isfield(sItem, 'ReferencedSeriesSequence')
                serSeq = sItem.ReferencedSeriesSequence;
                serF   = fieldnames(serSeq);
                for r = 1:numel(serF)
                    rItem = serSeq.(serF{r});
                    if isfield(rItem, 'SeriesInstanceUID')
                        fprintf('      Series UID: %s\n', rItem.SeriesInstanceUID);
                    end
                end
            end
        end
    else
        fprintf('  StudiesContainingOtherReferencedInstancesSequence: not present\n');
    end

    if isfield(info, 'ReferencedSeriesSequence')
        serSeq = info.ReferencedSeriesSequence;
        serF   = fieldnames(serSeq);
        fprintf('  ReferencedSeriesSequence (top-level): %d series\n', numel(serF));
        for r = 1:numel(serF)
            rItem = serSeq.(serF{r});
            if isfield(rItem, 'SeriesInstanceUID')
                fprintf('    Series UID: %s\n', rItem.SeriesInstanceUID);
            end
        end
    end
end


function print_field(info, tag, label)
    % print_field  Print a dicominfo top-level field if it exists
    if isfield(info, tag)
        val = info.(tag);
        if isnumeric(val)
            val = num2str(val);
        elseif iscell(val)
            val = strjoin(val, ', ');
        end
        fprintf('%s: %s\n', label, val);
    else
        fprintf('%s: [not present]\n', label);
    end
end


function print_seq_field(item, tag, label)
    % print_seq_field  Print a field from a sequence item struct if it exists
    if isfield(item, tag)
        val = item.(tag);
        if isnumeric(val)
            val = num2str(val);
        elseif iscell(val)
            val = strjoin(val, ', ');
        end
        fprintf('%s: %s\n', label, val);
    end
end
