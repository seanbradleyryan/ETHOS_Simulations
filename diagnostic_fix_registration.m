% diagnostic_fix_registration.m
%
% PURPOSE:
%   Loads a DICOM deformable registration file, swaps Item_1 and Item_2
%   inside DeformableRegistrationSequence, and writes the corrected file.
%
% USAGE:
%   Edit input_file below, then run (F5 or Run).
%   Output file is written to the same directory with "_fixed" appended.
%
% DEPENDENCIES:
%   MATLAB Image Processing Toolbox (dicominfo, dicomwrite)
%
% See also: dicominfo, dicomwrite, diagnostic_deformable_registration

% ------------------------------------------------------------------ %
%  CONFIGURATION                                                       %
% ------------------------------------------------------------------ %
input_file = 'C:\Users\80030361\ETHOS_Simulations\Raystation_Input\1194203\Session_5\REG_002.dcm';

% ------------------------------------------------------------------ %
%  BUILD OUTPUT PATH                                                   %
% ------------------------------------------------------------------ %
[dirPart, namePart, extPart] = fileparts(input_file);
output_file = fullfile(dirPart, [namePart '_fixed' extPart]);

% ------------------------------------------------------------------ %
%  LOAD                                                                %
% ------------------------------------------------------------------ %
fprintf('[FIX] Reading: %s\n', input_file);
info = dicominfo(input_file);

if ~isfield(info, 'DeformableRegistrationSequence')
    error('diagnostic_fix_registration:MissingField', ...
          'DeformableRegistrationSequence not found in %s', input_file);
end

seq = info.DeformableRegistrationSequence;

if ~isfield(seq, 'Item_1') || ~isfield(seq, 'Item_2')
    error('diagnostic_fix_registration:MissingItems', ...
          'DeformableRegistrationSequence must contain Item_1 and Item_2');
end

% ------------------------------------------------------------------ %
%  SWAP Item_1 and Item_2                                              %
% ------------------------------------------------------------------ %
fprintf('[FIX] Swapping Item_1 and Item_2 in DeformableRegistrationSequence...\n');

item1 = seq.Item_1;
item2 = seq.Item_2;

seq.Item_1 = item2;
seq.Item_2 = item1;

info.DeformableRegistrationSequence = seq;

% ------------------------------------------------------------------ %
%  WRITE                                                               %
% ------------------------------------------------------------------ %
fprintf('[FIX] Writing: %s\n', output_file);
dicomwrite([], output_file, info, 'CreateMode', 'copy');

fprintf('[FIX] Done.\n');
