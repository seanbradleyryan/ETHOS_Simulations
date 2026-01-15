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
    
    % Check for Item_N pattern (DICOM-style numbered fields)
    itemFields = regexp(fields, '^Item_(\d+)$', 'tokens');
    itemIndices = find(~cellfun(@isempty, itemFields));
    
    if ~isempty(itemIndices) && length(itemIndices) > 1
        % We have multiple Item_N fields - check if they're identical
        itemFieldNames = fields(itemIndices);
        firstItemFields = fieldnames(s.(itemFieldNames{1}));
        
        allSame = true;
        for idx = 2:length(itemFieldNames)
            if ~isequal(fieldnames(s.(itemFieldNames{idx})), firstItemFields)
                allSame = false;
                break;
            end
        end
        
        if allSame
            % Process non-Item fields first
            nonItemIndices = setdiff(1:length(fields), itemIndices);
            for i = nonItemIndices
                fieldName = fields{i};
                value = s.(fieldName);
                fullName = [name, '.', fieldName];
                
                if isstruct(value)
                    inspectStruct(value, fullName);
                else
                    sz = size(value);
                    szStr = strjoin(arrayfun(@numel2str, sz, 'UniformOutput', false), 'x');
                    cl = class(value);
                    fprintf('%-25s : %-10s %s\n', fullName, szStr, cl);
                end
            end
            
            % Show only the first Item_N and indicate how many exist
            firstItemName = itemFieldNames{1};
            fprintf('%-25s : [%d identical Item_N structures, showing first only]\n', ...
                    [name, '.', firstItemName], length(itemFieldNames));
            inspectStruct(s.(firstItemName), [name, '.', firstItemName]);
            return;
        end
    end
    
    % Normal processing - no Item_N pattern detected or they're not identical
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
