function diagnostic_fix_registration_v2(inputFilePath)
% PURPOSE:
%   Diagnose and fix a DICOM Deformable Registration file where Item 2 of
%   DeformableRegistrationSequence has an empty RegistrationTypeCodeSequence.
%   RayStation requires every item in the sequence to have a fully populated
%   RegistrationTypeCodeSequence (CodeValue, CodingSchemeDesignator, CodeMeaning).
%   This script deep-copies the populated code sequence from Item 1 into Item 2.
%
% INPUTS:
%   inputFilePath (optional, char) - Full path to the .dcm file to fix.
%       Defaults to: C:\Users\80030361\ETHOS_Simulations\Raystation_Input\
%                    1194203\Session_5\REG_002.dcm
%
% OUTPUTS:
%   A fixed .dcm file written to the same directory as the input, with
%   '_fixed' appended before the extension (e.g. REG_002_fixed.dcm).
%   Console output summarises all findings and changes made.
%
% ALGORITHM:
%   1. Accept file path argument or fall back to default path
%   2. Read DICOM metadata via dicominfo
%   3. Check Modality / SOPClassUID confirms a Deformable Registration object
%   4. Verify DeformableRegistrationSequence exists with >= 2 items
%   5. Validate Item 1 RegistrationTypeCodeSequence is fully populated
%   6. Inspect Item 2 RegistrationTypeCodeSequence for empty / missing items
%   7. If defective: deep-copy the code item from Item 1 into Item 2
%   8. Write corrected file with '_fixed' suffix using dicomwrite CreateMode copy
%   9. Print a clear before/after summary
%
% EXAMPLE:
%   diagnostic_fix_registration_v2()
%   diagnostic_fix_registration_v2('C:\path\to\REG_002.dcm')
%
% DEPENDENCIES:
%   MATLAB Image Processing Toolbox (dicominfo, dicomwrite)
%
% See also: dicominfo, dicomwrite, diagnostic_fix_registration,
%           diagnostic_deformable_registration

    % Deformable Spatial Registration Storage SOP Class UID (DICOM PS3.4)
    DEFORMABLE_REG_SOP_CLASS = '1.2.840.10008.5.1.4.1.1.66.3';
    DEFAULT_INPUT = ['C:\Users\80030361\ETHOS_Simulations\' ...
                     'Raystation_Input\1194203\Session_5\REG_002.dcm'];

    % ------------------------------------------------------------------ %
    %  1. Resolve input path                                               %
    % ------------------------------------------------------------------ %
    if nargin < 1 || isempty(inputFilePath)
        inputFilePath = DEFAULT_INPUT;
        fprintf('[INFO] No argument supplied. Using default path:\n');
        fprintf('       %s\n\n', inputFilePath);
    end
    if ~ischar(inputFilePath) && ~isstring(inputFilePath)
        error('diagnostic_fix_registration_v2:InvalidInput', ...
              'inputFilePath must be a character array or string.');
    end
    inputFilePath = char(strtrim(inputFilePath));

    printSeparator('=');
    fprintf(' diagnostic_fix_registration_v2\n');
    printSeparator('=');

    % ------------------------------------------------------------------ %
    %  2. Verify file exists and read DICOM metadata                      %
    % ------------------------------------------------------------------ %
    fprintf('[STEP 1] Checking file and reading DICOM metadata...\n');

    if ~isfile(inputFilePath)
        error('diagnostic_fix_registration_v2:FileNotFound', ...
              'File not found:\n  %s\nCheck the path and try again.', inputFilePath);
    end
    fprintf('  File found: %s\n', inputFilePath);

    try
        info = dicominfo(inputFilePath, 'UseDictionaryVR', true);
    catch ME
        error('diagnostic_fix_registration_v2:DicomReadError', ...
              'dicominfo failed on file:\n  %s\nMATLAB error: %s', ...
              inputFilePath, ME.message);
    end
    fprintf('  [OK] DICOM metadata loaded.\n\n');

    % ------------------------------------------------------------------ %
    %  3. Confirm this is a Deformable Registration object                 %
    % ------------------------------------------------------------------ %
    fprintf('[STEP 2] Validating DICOM object type...\n');

    sopUID   = strtrim(getFieldSafe(info, 'SOPClassUID',  ''));
    modality = strtrim(getFieldSafe(info, 'Modality',     ''));

    isDeformable = strcmp(sopUID, DEFORMABLE_REG_SOP_CLASS) || ...
                   strcmpi(modality, 'REG');

    if isDeformable
        fprintf('  Modality:    %s\n', modality);
        fprintf('  SOPClassUID: %s\n', sopUID);
        fprintf('  [OK] Confirmed as Deformable Spatial Registration object.\n\n');
    else
        fprintf('  [WARNING] This file does not appear to be a Deformable Registration object.\n');
        fprintf('    Modality:    %s\n', modality);
        fprintf('    SOPClassUID: %s\n', sopUID);
        fprintf('    Expected Modality REG or SOPClassUID %s\n', DEFORMABLE_REG_SOP_CLASS);
        fprintf('  Proceeding anyway — results may be unreliable.\n\n');
    end

    % ------------------------------------------------------------------ %
    %  4. Verify DeformableRegistrationSequence exists with >= 2 items    %
    % ------------------------------------------------------------------ %
    fprintf('[STEP 3] Inspecting DeformableRegistrationSequence...\n');

    if ~isfield(info, 'DeformableRegistrationSequence')
        error('diagnostic_fix_registration_v2:MissingSequence', ...
              'DeformableRegistrationSequence tag is absent from this file.\n%s\n%s', ...
              'This is required for deformable registration objects.', ...
              'File may be corrupt or incorrectly exported.');
    end

    drs      = info.DeformableRegistrationSequence;
    allItems = fieldnames(drs);
    seqItems = allItems(startsWith(allItems, 'Item_'));
    numItems = numel(seqItems);

    fprintf('  DeformableRegistrationSequence item count: %d\n', numItems);

    if numItems < 2
        error('diagnostic_fix_registration_v2:InsufficientItems', ...
              'DeformableRegistrationSequence has only %d item(s); at least 2 are required.\n%s', ...
              numItems, ...
              'The file may already be incomplete or the wrong REG file was selected.');
    end
    fprintf('  [OK] At least 2 items found.\n\n');

    % ------------------------------------------------------------------ %
    %  5. Validate Item 1 RegistrationTypeCodeSequence                    %
    % ------------------------------------------------------------------ %
    fprintf('[STEP 4] Validating Item 1 RegistrationTypeCodeSequence (copy source)...\n');

    item1 = drs.Item_1;
    [item1Valid, codeFromItem1, item1Msg] = inspectCodeSequence(item1, 'Item 1');

    if ~item1Valid
        fprintf('  [FAIL] %s\n', item1Msg);
        error('diagnostic_fix_registration_v2:Item1Invalid', ...
              ['Item 1 RegistrationTypeCodeSequence is missing or empty.\n' ...
               'Cannot copy from Item 1. Manual inspection is required.\n' ...
               'Details: %s'], item1Msg);
    end

    fprintf('  CodeValue:              "%s"\n', getFieldSafe(codeFromItem1, 'CodeValue',              ''));
    fprintf('  CodingSchemeDesignator: "%s"\n', getFieldSafe(codeFromItem1, 'CodingSchemeDesignator', ''));
    fprintf('  CodeMeaning:            "%s"\n', getFieldSafe(codeFromItem1, 'CodeMeaning',            ''));
    fprintf('  [OK] Item 1 is fully populated and will serve as the copy source.\n\n');

    % ------------------------------------------------------------------ %
    %  6. Inspect Item 2 RegistrationTypeCodeSequence                     %
    % ------------------------------------------------------------------ %
    fprintf('[STEP 5] Inspecting Item 2 RegistrationTypeCodeSequence...\n');

    item2 = drs.Item_2;
    [item2Valid, ~, item2Msg] = inspectCodeSequence(item2, 'Item 2');

    needsFix = ~item2Valid;

    if needsFix
        fprintf('  [ISSUE] %s\n', item2Msg);
        fprintf('  --> This is the root cause of the RayStation import failure:\n');
        fprintf('      "Registration Type Code Sequence must contain exactly one item"\n\n');
    else
        fprintf('  [OK] Item 2 RegistrationTypeCodeSequence is already populated.\n');
        fprintf('  No structural fix is required. File will be re-saved as-is.\n\n');
    end

    % ------------------------------------------------------------------ %
    %  7. Apply fix in memory                                              %
    % ------------------------------------------------------------------ %
    if needsFix
        fprintf('[STEP 6] Applying fix: copying RegistrationTypeCodeSequence from Item 1 to Item 2...\n');

        % Build a fresh sequence struct with one item — a deep copy of Item 1's code
        copiedSeq.Item_1 = codeFromItem1;
        info.DeformableRegistrationSequence.Item_2.RegistrationTypeCodeSequence = copiedSeq;

        fprintf('  Injected into Item 2:\n');
        fprintf('    CodeValue:              "%s"\n', getFieldSafe(codeFromItem1, 'CodeValue',              ''));
        fprintf('    CodingSchemeDesignator: "%s"\n', getFieldSafe(codeFromItem1, 'CodingSchemeDesignator', ''));
        fprintf('    CodeMeaning:            "%s"\n', getFieldSafe(codeFromItem1, 'CodeMeaning',            ''));
        fprintf('  [OK] Fix applied in memory.\n\n');
    end

    % ------------------------------------------------------------------ %
    %  8. Write corrected file                                             %
    % ------------------------------------------------------------------ %
    fprintf('[STEP 7] Writing corrected DICOM file...\n');

    [dirPart, namePart, extPart] = fileparts(inputFilePath);
    outputFilePath = fullfile(dirPart, [namePart '_fixed' extPart]);
    fprintf('  Output path: %s\n', outputFilePath);

    try
        dicomwrite([], outputFilePath, info, 'CreateMode', 'copy', 'WritePrivate', true);
    catch ME
        error('diagnostic_fix_registration_v2:WriteError', ...
              'dicomwrite failed:\n  Output: %s\n  MATLAB error: %s', ...
              outputFilePath, ME.message);
    end
    fprintf('  [OK] File written successfully.\n\n');

    % ------------------------------------------------------------------ %
    %  9. Summary                                                          %
    % ------------------------------------------------------------------ %
    printSeparator('=');
    fprintf(' SUMMARY\n');
    printSeparator('=');
    fprintf('  Input file:  %s\n', inputFilePath);
    fprintf('  Output file: %s\n', outputFilePath);
    fprintf('\n');
    if needsFix
        fprintf('  STATUS: FIX APPLIED\n');
        fprintf('  Problem: Item 2 RegistrationTypeCodeSequence was empty / absent.\n');
        fprintf('  Action:  Copied RegistrationTypeCodeSequence.Item_1 from Item 1 to Item 2.\n');
        fprintf('\n');
        fprintf('  Copied values:\n');
        fprintf('    CodeValue              = "%s"\n', getFieldSafe(codeFromItem1, 'CodeValue',              ''));
        fprintf('    CodingSchemeDesignator = "%s"\n', getFieldSafe(codeFromItem1, 'CodingSchemeDesignator', ''));
        fprintf('    CodeMeaning            = "%s"\n', getFieldSafe(codeFromItem1, 'CodeMeaning',            ''));
    else
        fprintf('  STATUS: NO FIX REQUIRED\n');
        fprintf('  Both items already have valid RegistrationTypeCodeSequences.\n');
        fprintf('  File re-saved with no structural changes.\n');
    end
    fprintf('\n');
    fprintf('  Next step: Import "%s_fixed%s" into RayStation.\n', namePart, extPart);
    printSeparator('=');

end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function [isValid, codeItem, msg] = inspectCodeSequence(drsItem, label)
% Inspect a single DeformableRegistrationSequence item for a valid
% RegistrationTypeCodeSequence.  Returns:
%   isValid  - true if the sequence has >= 1 item with all 3 required tags
%   codeItem - the first code item struct (or [] if invalid)
%   msg      - human-readable description of what was found

    isValid  = false;
    codeItem = [];

    if ~isfield(drsItem, 'RegistrationTypeCodeSequence')
        msg = sprintf('%s: RegistrationTypeCodeSequence tag is ABSENT.', label);
        return;
    end

    rtcs = drsItem.RegistrationTypeCodeSequence;

    % Empty struct / not a struct at all
    if isempty(rtcs) || ~isstruct(rtcs)
        msg = sprintf('%s: RegistrationTypeCodeSequence is EMPTY (not a struct).', label);
        return;
    end

    % Look for Item_ sub-fields
    codeFields = fieldnames(rtcs);
    codeFields = codeFields(startsWith(codeFields, 'Item_'));

    if isempty(codeFields)
        msg = sprintf('%s: RegistrationTypeCodeSequence contains NO Item_ entries.', label);
        return;
    end

    firstCode = rtcs.(codeFields{1});

    hasCodeValue   = isfield(firstCode, 'CodeValue')              && ~isempty(firstCode.CodeValue);
    hasDesignator  = isfield(firstCode, 'CodingSchemeDesignator') && ~isempty(firstCode.CodingSchemeDesignator);
    hasMeaning     = isfield(firstCode, 'CodeMeaning')            && ~isempty(firstCode.CodeMeaning);

    if ~(hasCodeValue && hasDesignator && hasMeaning)
        msg = sprintf(['%s: RegistrationTypeCodeSequence item is INCOMPLETE.\n' ...
                       '    CodeValue present:              %s\n' ...
                       '    CodingSchemeDesignator present: %s\n' ...
                       '    CodeMeaning present:            %s'], ...
                       label, ...
                       tf2str(hasCodeValue), tf2str(hasDesignator), tf2str(hasMeaning));
        return;
    end

    isValid  = true;
    codeItem = firstCode;
    msg      = sprintf('%s: RegistrationTypeCodeSequence is fully populated.', label);
end


function val = getFieldSafe(s, fieldName, defaultVal)
% Return strtrim(char(s.fieldName)) if it exists and is non-empty,
% otherwise return defaultVal.
    if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
        val = strtrim(char(s.(fieldName)));
    else
        val = defaultVal;
    end
end


function printSeparator(ch)
% Print a 60-character separator line.
    fprintf('%s\n', repmat(ch, 1, 60));
end


function s = tf2str(val)
% Convert logical to 'YES' / 'NO' for display.
    if val
        s = 'YES';
    else
        s = 'NO';
    end
end
