function diagnostic_deformable_registration(search_root)
% diagnostic_deformable_registration  Inspect ETHOS deformable registration DICOMs
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
% INPUTS:
%   search_root - (char, optional) Root directory to search recursively.
%                 Defaults to 'C:/Users/seanr/ETHOS_Simulations/EthosExports'
%
% OUTPUTS:
%   Console report only. Variables are left in the base workspace via
%   assignin() for interactive inspection after the script runs.
%
% ALGORITHM:
%   1. Walk search_root recursively for all DICOM files
%   2. Filter to REG modality (Spatial / Deformable Spatial Registration)
%   3. For each REG file, print header metadata
%   4. Inspect DeformableRegistrationSequence items:
%      a. Count RegistrationTypeCodeSequence items (must be exactly 1)
%      b. Print code values found
%      c. Check SourceFrameOfReferenceUID
%      d. Check DeformableRegistrationGridSequence presence
%   5. Flag any non-conformances
%
% EXAMPLE:
%   diagnostic_deformable_registration()
%   diagnostic_deformable_registration('C:/MyData/Patient001/Session1')
%
% DEPENDENCIES:
%   MATLAB Image Processing Toolbox (dicominfo, dicomuid)
%
% See also: dicominfo, step0_sort_dicom

    % ------------------------------------------------------------------ %
    %  0. DEFAULTS & INPUT VALIDATION                                     %
    % ------------------------------------------------------------------ %
    SOP_DEFORMABLE = '1.2.840.10008.5.1.4.1.1.66.3';
    SOP_SPATIAL    = '1.2.840.10008.5.1.4.1.1.66.1';

    if nargin < 1 || isempty(search_root)
        search_root = 'C:/Users/seanr/ETHOS_Simulations/EthosExports';
    end

    if ~ischar(search_root) && ~isstring(search_root)
        error('diagnostic_deformable_registration:InvalidInput', ...
              'search_root must be a char or string path.');
    end
    search_root = char(search_root);

    if ~isfolder(search_root)
        error('diagnostic_deformable_registration:NotFound', ...
              'Directory not found: %s', search_root);
    end

    % ------------------------------------------------------------------ %
    %  1. RECURSIVE FILE DISCOVERY                                        %
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
    %  2. FILTER TO REGISTRATION DICOM FILES                              %
    % ------------------------------------------------------------------ %
    fprintf('[STEP 2] Filtering to REG / registration DICOM files...\n');

    regFiles    = struct('path', {}, 'info', {}, 'sop_class', {});
    skippedCount = 0;

    for k = 1:numel(allFiles)
        fpath = fullfile(allFiles(k).folder, allFiles(k).name);

        % Quick filter: skip obviously non-DICOM extensions
        [~, ~, ext] = fileparts(fpath);
        if ~isempty(ext) && ~strcmpi(ext, '.dcm') && ~isempty(ext)
            % Still attempt if no extension (many DICOM files lack .dcm)
            if ~isempty(ext)
                skippedCount = skippedCount + 1;
                continue;
            end
        end

        try
            info = dicominfo(fpath);
        catch
            skippedCount = skippedCount + 1;
            continue;
        end

        % Check SOP class or Modality
        sopClass = '';
        if isfield(info, 'SOPClassUID')
            sopClass = info.SOPClassUID;
        end

        modality = '';
        if isfield(info, 'Modality')
            modality = info.Modality;
        end

        isRegFile = strcmpi(modality, 'REG') || ...
                    strcmp(sopClass, SOP_DEFORMABLE) || ...
                    strcmp(sopClass, SOP_SPATIAL);

        if isRegFile
            regFiles(end+1).path      = fpath; %#ok<AGROW>
            regFiles(end).info         = info;
            regFiles(end).sop_class    = sopClass;
        end
    end

    fprintf('  Skipped %d non-DICOM / unreadable files.\n', skippedCount);
    fprintf('  Found %d REG/registration DICOM file(s).\n\n', numel(regFiles));

    if isempty(regFiles)
        fprintf('[WARNING] No deformable or spatial registration DICOM files found.\n');
        fprintf('  Check that the ETHOS export includes REG files and that\n');
        fprintf('  search_root points to the correct export directory.\n\n');
        return;
    end

    % Export discovery results to base workspace for interactive use
    assignin('base', 'diag_regFiles', regFiles);
    fprintf('  -> Variable ''diag_regFiles'' saved to workspace (%d entries).\n\n', numel(regFiles));

    % ------------------------------------------------------------------ %
    %  3. PER-FILE HEADER REPORT & DIAGNOSTICS                           %
    % ------------------------------------------------------------------ %
    totalErrors   = 0;
    totalWarnings = 0;

    for f = 1:numel(regFiles)
        info    = regFiles(f).info;
        fpath   = regFiles(f).path;
        sopUID  = regFiles(f).sop_class;

        fprintf('%s\n', repmat('-', 1, 70));
        fprintf('FILE %d/%d: %s\n', f, numel(regFiles), fpath);
        fprintf('%s\n', repmat('-', 1, 70));

        % --- 3a. Header metadata ---
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

        % Classify SOP class
        fprintf('\n  SOP Class type: ');
        if strcmp(sopUID, SOP_DEFORMABLE)
            fprintf('Deformable Spatial Registration (1.2.840.10008.5.1.4.1.1.66.3)\n');
        elseif strcmp(sopUID, SOP_SPATIAL)
            fprintf('Spatial Registration (Rigid) (1.2.840.10008.5.1.4.1.1.66.1)\n');
        else
            fprintf('Other/Unknown: %s\n', sopUID);
        end

        % --- 3b. DeformableRegistrationSequence inspection ---
        if strcmp(sopUID, SOP_DEFORMABLE)
            fprintf('\n  [DEFORMABLE REGISTRATION SEQUENCE INSPECTION]\n');
            [nErr, nWarn] = inspect_deformable_sequence(info);
            totalErrors   = totalErrors   + nErr;
            totalWarnings = totalWarnings + nWarn;
        end

        % --- 3c. Spatial (rigid) RegistrationSequence inspection ---
        if strcmp(sopUID, SOP_SPATIAL)
            fprintf('\n  [SPATIAL REGISTRATION SEQUENCE INSPECTION]\n');
            [nErr, nWarn] = inspect_spatial_sequence(info);
            totalErrors   = totalErrors   + nErr;
            totalWarnings = totalWarnings + nWarn;
        end

        % --- 3d. Referenced Image / Frame of Reference coverage ---
        fprintf('\n  [REFERENCED IMAGE SERIES]\n');
        inspect_referenced_series(info);

        fprintf('\n');
    end

    % ------------------------------------------------------------------ %
    %  4. SUMMARY                                                         %
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
end


% ======================================================================= %
%  LOCAL HELPER FUNCTIONS                                                  %
% ======================================================================= %

function inspect_deformable_sequence(info)
    % inspect_deformable_sequence  Check DeformableRegistrationSequence
    %
    % Inspects each item in DeformableRegistrationSequence and validates:
    %   - RegistrationTypeCodeSequence has exactly 1 item (RayStation requirement)
    %   - SourceFrameOfReferenceUID is present
    %   - DeformableRegistrationGridSequence is present and non-empty

    nErrors   = 0;
    nWarnings = 0;

    if ~isfield(info, 'DeformableRegistrationSequence')
        fprintf('  [ERROR] DeformableRegistrationSequence is MISSING from this file.\n');
        nErrors = nErrors + 1;
        return;
    end

    seqData = info.DeformableRegistrationSequence;

    % MATLAB dicominfo stores sequences as structs with Item_1, Item_2, ...
    itemNames = fieldnames(seqData);
    nItems    = numel(itemNames);
    fprintf('  DeformableRegistrationSequence: %d item(s)\n', nItems);

    for i = 1:nItems
        item = seqData.(itemNames{i});
        fprintf('\n    --- Sequence Item %d ---\n', i);

        % SourceFrameOfReferenceUID
        if isfield(item, 'SourceFrameOfReferenceUID')
            fprintf('    SourceFrameOfReferenceUID : %s\n', item.SourceFrameOfReferenceUID);
        else
            fprintf('    [WARNING] SourceFrameOfReferenceUID is missing.\n');
            nWarnings = nWarnings + 1;
        end

        % RegistrationTypeCodeSequence — THE KEY CHECK
        fprintf('\n    RegistrationTypeCodeSequence:\n');
        if ~isfield(item, 'RegistrationTypeCodeSequence')
            fprintf('      [ERROR] RegistrationTypeCodeSequence is MISSING.\n');
            fprintf('      RayStation requires exactly 1 item here.\n');
            nErrors = nErrors + 1;
        else
            rtcSeq    = item.RegistrationTypeCodeSequence;
            rtcItems  = fieldnames(rtcSeq);
            nRtcItems = numel(rtcItems);
            fprintf('      Item count: %d  ', nRtcItems);

            if nRtcItems == 1
                fprintf('[OK - exactly 1 item]\n');
            elseif nRtcItems == 0
                fprintf('[ERROR - 0 items, RayStation requires exactly 1]\n');
                nErrors = nErrors + 1;
            else
                fprintf('[ERROR - %d items found, RayStation requires exactly 1]\n', nRtcItems);
                nErrors = nErrors + 1;
            end

            % Print what's in each item of RegistrationTypeCodeSequence
            for j = 1:nRtcItems
                rtcItem = rtcSeq.(rtcItems{j});
                fprintf('      Item %d:\n', j);
                if isfield(rtcItem, 'CodeValue')
                    fprintf('        CodeValue            : %s\n', rtcItem.CodeValue);
                else
                    fprintf('        CodeValue            : [MISSING]\n');
                    nWarnings = nWarnings + 1;
                end
                if isfield(rtcItem, 'CodingSchemeDesignator')
                    fprintf('        CodingSchemeDesignator: %s\n', rtcItem.CodingSchemeDesignator);
                else
                    fprintf('        CodingSchemeDesignator: [MISSING]\n');
                    nWarnings = nWarnings + 1;
                end
                if isfield(rtcItem, 'CodeMeaning')
                    fprintf('        CodeMeaning          : %s\n', rtcItem.CodeMeaning);
                else
                    fprintf('        CodeMeaning          : [MISSING]\n');
                    nWarnings = nWarnings + 1;
                end
            end
        end

        % DeformableRegistrationGridSequence
        fprintf('\n    DeformableRegistrationGridSequence:\n');
        if ~isfield(item, 'DeformableRegistrationGridSequence')
            fprintf('      [ERROR] MISSING — no deformation grid found in this item.\n');
            nErrors = nErrors + 1;
        else
            gridSeq   = item.DeformableRegistrationGridSequence;
            gridItems = fieldnames(gridSeq);
            fprintf('      Present. %d grid item(s).\n', numel(gridItems));

            for g = 1:numel(gridItems)
                gItem = gridSeq.(gridItems{g});
                fprintf('      Grid Item %d:\n', g);
                if isfield(gItem, 'ImageOrientationPatient')
                    fprintf('        ImageOrientationPatient: [%s]\n', ...
                        num2str(gItem.ImageOrientationPatient(:)', '%.4f '));
                end
                if isfield(gItem, 'ImagePositionPatient')
                    fprintf('        ImagePositionPatient   : [%s] mm\n', ...
                        num2str(gItem.ImagePositionPatient(:)', '%.2f '));
                end
                if isfield(gItem, 'PixelSpacing')
                    fprintf('        PixelSpacing           : [%s] mm\n', ...
                        num2str(gItem.PixelSpacing(:)', '%.3f '));
                end
                if isfield(gItem, 'GridFrameOffsetVector')
                    offsets = gItem.GridFrameOffsetVector;
                    fprintf('        GridFrameOffsetVector  : %d elements, range [%.2f, %.2f] mm\n', ...
                        numel(offsets), min(offsets), max(offsets));
                end
                if isfield(gItem, 'Rows') && isfield(gItem, 'Columns')
                    nSlices = 0;
                    if isfield(gItem, 'GridFrameOffsetVector')
                        nSlices = numel(gItem.GridFrameOffsetVector);
                    end
                    fprintf('        Grid dimensions        : %d rows x %d cols x %d slices\n', ...
                        gItem.Rows, gItem.Columns, nSlices);
                end
                if isfield(gItem, 'VectorGridData')
                    vgd = gItem.VectorGridData;
                    fprintf('        VectorGridData size    : %s bytes raw\n', ...
                        num2str(numel(vgd)));
                else
                    fprintf('        [WARNING] VectorGridData tag not found.\n');
                    nWarnings = nWarnings + 1;
                end
            end
        end

        % PreDeformationMatrixRegistrationSequence (optional)
        if isfield(item, 'PreDeformationMatrixRegistrationSequence')
            fprintf('\n    PreDeformationMatrixRegistrationSequence: present\n');
            preSeq = item.PreDeformationMatrixRegistrationSequence;
            preItems = fieldnames(preSeq);
            for p = 1:numel(preItems)
                pItem = preSeq.(preItems{p});
                if isfield(pItem, 'FrameOfReferenceTransformationMatrix')
                    M = reshape(pItem.FrameOfReferenceTransformationMatrix, 4, 4)';
                    fprintf('      Pre-matrix item %d:\n', p);
                    for row = 1:4
                        fprintf('        [ %10.5f %10.5f %10.5f %10.5f ]\n', M(row,:));
                    end
                end
            end
        end
    end

    % Return counts (assign to caller via output args is not possible in
    % a void local function — counts are printed above)
end

function [nErrors, nWarnings] = inspect_deformable_sequence(info)
    % Wrapper that returns error/warning counts (overloaded via output args)
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
                pf = @(tag, label) print_seq_field(rtcItem, tag, label);
                pf('CodeValue',             '        CodeValue            ');
                pf('CodingSchemeDesignator','        CodingSchemeDesignator');
                pf('CodingSchemeVersion',   '        CodingSchemeVersion   ');
                pf('CodeMeaning',           '        CodeMeaning          ');
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

        % ---- Optional pre/post deformation rigid components ----
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
    % inspect_spatial_sequence  Check RegistrationSequence in rigid REG files
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

        % RegistrationTypeCodeSequence check (same rule applies here)
        fprintf('\n    RegistrationTypeCodeSequence:\n');
        if ~isfield(item, 'RegistrationTypeCodeSequence')
            fprintf('      [WARNING] MISSING (may be optional for spatial reg).\n');
            nWarnings = nWarnings + 1;
        else
            rtcSeq  = item.RegistrationTypeCodeSequence;
            rtcF    = fieldnames(rtcSeq);
            nRtc    = numel(rtcF);
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

        % Matrix
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
    % inspect_referenced_series  List referenced series/studies in REG file

    % Check ReferencedSeriesSequence at top level or inside study seq
    if isfield(info, 'StudiesContainingOtherReferencedInstancesSequence')
        studySeq  = info.StudiesContainingOtherReferencedInstancesSequence;
        studyF    = fieldnames(studySeq);
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
    % print_field  Print a dicominfo field if it exists
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
    % print_seq_field  Print a field from a sequence item struct
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
