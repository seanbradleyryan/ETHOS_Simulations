function inspectStruct(s, name)
    % If no name is provided, start with 'input'
    if nargin < 2
        name = 'input';
    end

    % Handle struct arrays - check for redundancy
    if isstruct(s) && numel(s) > 1
        % Check if all elements have identical field structure
        allSame = true;
        firstFields = fieldnames(s(1));
        for idx = 2:numel(s)
            if ~isequal(fieldnames(s(idx)), firstFields)
                allSame = false;
                break;
            end
        end
        
        if allSame
            % Identical structures - only show first element
            sz = size(s);
            szStr = strjoin(arrayfun(@numel2str, sz, 'UniformOutput', false), 'x');
            fprintf('%-25s : [%s struct array - identical structure, showing first only]\n', name, szStr);
            inspectStruct(s(1), [name, '(1)']);
            return;
        else
            % Different structures - show each separately (rare case)
            for idx = 1:numel(s)
                inspectStruct(s(idx), sprintf('%s(%d)', name, idx));
            end
            return;
        end
    end

    % Handle non-struct input
    if ~isstruct(s)
        sz = size(s);
        szStr = strjoin(arrayfun(@numel2str, sz, 'UniformOutput', false), 'x');
        cl = class(s);
        fprintf('%-25s : %-10s %s\n', name, szStr, cl);
        return;
    end

    % Process single struct's fields
    fields = fieldnames(s);

    for i = 1:numel(fields)
        fieldName = fields{i};
        value = s.(fieldName);
        fullName = [name, '.', fieldName];

        if isstruct(value)
            % Recurse into structure
            inspectStruct(value, fullName);
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
```

**Key changes:**

1. **Detects struct arrays**: When the function encounters a struct array (multiple elements), it checks if all elements have identical field names
2. **Collapses redundancy**: If they're identical, it prints a summary line showing the array dimensions and only processes the first element
3. **Recursive collapse**: This applies at every level, so your 17 items with 166 sub-items each will show as two collapsed arrays instead of 2,822 individual entries

**Example output:**
```
input.items            : [17x1 struct array - identical structure, showing first only]
input.items(1).data    : [166x1 struct array - identical structure, showing first only]
input.items(1).data(1).field1 : 100x50    double
input.items(1).data(1).field2 : 1x1       char