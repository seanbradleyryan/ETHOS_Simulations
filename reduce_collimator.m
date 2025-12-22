function pln = reduce_collimator(pln)
% REDUCE_COLLIMATOR Convert Halcyon dual-layer MLC to single-layer for matRad
%
% Description:
%   Varian Halcyon uses a dual-layer stacked MLC design where two banks of
%   leaves work together to create apertures. matRad expects single-layer
%   MLC data. This function combines the dual layers by selecting the most
%   restrictive leaf positions:
%     - For left (proximal) bank: max() of both layers
%     - For right (distal) bank: min() of both layers
%
% Usage:
%   pln = reduce_collimator(pln)
%
% Input:
%   pln - matRad planning structure containing beam shape information
%         in pln.propStf.beam(i).shape
%
% Output:
%   pln - Modified planning structure with single-layer MLC positions
%
% Notes:
%   - Call this function AFTER DICOM import but BEFORE STF generation
%   - Modifies pln.propStf.beam(i).shape in place
%   - Halcyon has 2x28 leaf pairs per layer = 56 leaves per bank per layer
%   - Combined single layer has 28 leaf pairs
%
% Author: Generated for ETHOS Simulation Analysis
% Date: 2024

fprintf('  [reduce_collimator] Starting dual-layer MLC reduction...\n');

% Check if pln structure has the required fields
if ~isfield(pln, 'propStf')
    warning('reduce_collimator: pln.propStf not found. Returning unchanged.');
    return;
end

if ~isfield(pln.propStf, 'beam')
    warning('reduce_collimator: pln.propStf.beam not found. Returning unchanged.');
    return;
end

numBeams = length(pln.propStf.beam);
fprintf('    Processing %d beams...\n', numBeams);

% Process each beam
for beamIdx = 1:numBeams
    
    beam = pln.propStf.beam(beamIdx);
    
    % Check if shape data exists
    if ~isfield(beam, 'shape')
        fprintf('    Beam %d: No shape data found, skipping\n', beamIdx);
        continue;
    end
    
    % Process each control point / segment
    numShapes = length(beam.shape);
    
    for shapeIdx = 1:numShapes
        
        shape = beam.shape(shapeIdx);
        
        % Check for MLC leaf positions
        % Halcyon dual-layer may have different field names depending on import
        
        % Try different possible field names for MLC data
        mlcFieldNames = {'MLCPositions', 'mlcPositions', 'leafPositions', ...
                        'LeafPositions', 'MLC', 'mlc'};
        
        mlcData = [];
        mlcFieldUsed = '';
        
        for fn = 1:length(mlcFieldNames)
            if isfield(shape, mlcFieldNames{fn})
                mlcData = shape.(mlcFieldNames{fn});
                mlcFieldUsed = mlcFieldNames{fn};
                break;
            end
        end
        
        if isempty(mlcData)
            % Try to find MLC data in nested structures
            if isfield(shape, 'collimator')
                if isfield(shape.collimator, 'leafPositions')
                    mlcData = shape.collimator.leafPositions;
                    mlcFieldUsed = 'collimator.leafPositions';
                elseif isfield(shape.collimator, 'MLC')
                    mlcData = shape.collimator.MLC;
                    mlcFieldUsed = 'collimator.MLC';
                end
            end
        end
        
        if isempty(mlcData)
            continue; % No MLC data found for this shape
        end
        
        % Determine MLC configuration
        [numRows, numCols] = size(mlcData);
        
        % Halcyon dual-layer configurations:
        % Option 1: [2*numLeaves x 2] - stacked leaves, 2 banks
        % Option 2: [numLeaves x 4] - 2 layers x 2 banks
        % Option 3: [4 x numLeaves] - transposed version
        
        % Standard Halcyon has 28 leaf pairs per layer = 56 per layer total
        % Dual layer = 112 leaves total for both banks
        
        if numRows == 56 && numCols == 2
            % Format: [56 leaves x 2 banks] - single layer 28 pairs
            % Each row: [left_position, right_position]
            % This might already be reduced, check if further reduction needed
            
            % Assume this is already dual-layer combined into 56 leaves
            % Need to reduce to 28 leaf pairs
            
            numLeafPairs = 28;
            reducedMLC = zeros(numLeafPairs, 2);
            
            % Layer 1: leaves 1-28, Layer 2: leaves 29-56
            for lp = 1:numLeafPairs
                layer1_left = mlcData(lp, 1);
                layer2_left = mlcData(lp + numLeafPairs, 1);
                layer1_right = mlcData(lp, 2);
                layer2_right = mlcData(lp + numLeafPairs, 2);
                
                % Most restrictive: max for left bank, min for right bank
                reducedMLC(lp, 1) = max(layer1_left, layer2_left);
                reducedMLC(lp, 2) = min(layer1_right, layer2_right);
            end
            
            mlcData = reducedMLC;
            
        elseif numRows == 112 && numCols == 2
            % Format: [112 leaves x 2 banks] - dual layer, 28 pairs each
            % Leaves 1-28: Layer 1 proximal, 29-56: Layer 1 distal
            % Leaves 57-84: Layer 2 proximal, 85-112: Layer 2 distal
            
            numLeafPairs = 28;
            reducedMLC = zeros(numLeafPairs, 2);
            
            for lp = 1:numLeafPairs
                % Get all layer positions for this leaf pair
                layer1_left = mlcData(lp, 1);
                layer1_right = mlcData(lp + numLeafPairs, 2);
                layer2_left = mlcData(lp + 2*numLeafPairs, 1);
                layer2_right = mlcData(lp + 3*numLeafPairs, 2);
                
                % Most restrictive combination
                reducedMLC(lp, 1) = max(layer1_left, layer2_left);
                reducedMLC(lp, 2) = min(layer1_right, layer2_right);
            end
            
            mlcData = reducedMLC;
            
        elseif numCols == 4 && numRows >= 28
            % Format: [numLeaves x 4] - [layer1_left, layer1_right, layer2_left, layer2_right]
            
            numLeafPairs = numRows;
            reducedMLC = zeros(numLeafPairs, 2);
            
            for lp = 1:numLeafPairs
                layer1_left = mlcData(lp, 1);
                layer1_right = mlcData(lp, 2);
                layer2_left = mlcData(lp, 3);
                layer2_right = mlcData(lp, 4);
                
                reducedMLC(lp, 1) = max(layer1_left, layer2_left);
                reducedMLC(lp, 2) = min(layer1_right, layer2_right);
            end
            
            mlcData = reducedMLC;
            
        elseif numRows == 4 && numCols >= 28
            % Transposed version of above
            mlcData = mlcData';
            [numRows, numCols] = size(mlcData);
            
            numLeafPairs = numRows;
            reducedMLC = zeros(numLeafPairs, 2);
            
            for lp = 1:numLeafPairs
                layer1_left = mlcData(lp, 1);
                layer1_right = mlcData(lp, 2);
                layer2_left = mlcData(lp, 3);
                layer2_right = mlcData(lp, 4);
                
                reducedMLC(lp, 1) = max(layer1_left, layer2_left);
                reducedMLC(lp, 2) = min(layer1_right, layer2_right);
            end
            
            mlcData = reducedMLC;
            
        elseif numRows == 28 && numCols == 2
            % Already in single-layer format, no reduction needed
            % Keep as is
            
        else
            % Unknown format - try generic reduction
            fprintf('    Beam %d, Shape %d: Unknown MLC format [%d x %d], attempting generic reduction\n', ...
                beamIdx, shapeIdx, numRows, numCols);
            
            % If more than 2 columns, assume dual-layer format
            if numCols > 2
                numLeafPairs = numRows;
                reducedMLC = zeros(numLeafPairs, 2);
                
                % Assume columns are: [left1, right1, left2, right2, ...]
                numLayers = floor(numCols / 2);
                
                for lp = 1:numLeafPairs
                    leftPositions = mlcData(lp, 1:2:end);
                    rightPositions = mlcData(lp, 2:2:end);
                    
                    reducedMLC(lp, 1) = max(leftPositions);
                    reducedMLC(lp, 2) = min(rightPositions);
                end
                
                mlcData = reducedMLC;
            end
        end
        
        % Update the shape with reduced MLC data
        if contains(mlcFieldUsed, '.')
            % Nested structure
            pln.propStf.beam(beamIdx).shape(shapeIdx).collimator.leafPositions = mlcData;
        else
            pln.propStf.beam(beamIdx).shape(shapeIdx).(mlcFieldUsed) = mlcData;
        end
        
    end % shape loop
    
    fprintf('    Beam %d: Processed %d control points\n', beamIdx, numShapes);
    
end % beam loop

fprintf('  [reduce_collimator] MLC reduction complete.\n');

end
