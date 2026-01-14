function inspectStruct(s, name)
    % If no name is provided, start with 'data'
    if nargin < 2
        name = 'input';
    end

    % Get all field names at the current level
    fields = fieldnames(s);

    for i = 1:numel(fields)
        fieldName = fields{i};
        value = s.(fieldName);
        fullName = [name, '.', fieldName];

        if isstruct(value)
            % If it's a structure, print its info and dive deeper
            fprintf('%-25s : [Structure]\n', fullName);
            inspectStruct(value, fullName); % <--- Recursive call
        else
            % Get size and class information
            sz = size(value);
            szStr = strjoin(arrayfun(@numel2str, sz, 'UniformOutput', false), 'x');
            cl = class(value);
            
            % Print formatted output
            fprintf('%-25s : %-10s %s\n', fullName, szStr, cl);
        end
    end
end

function s = numel2str(n)
    s = num2str(n);
end