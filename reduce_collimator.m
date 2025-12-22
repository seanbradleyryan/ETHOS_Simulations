function pln = reduce_collimator(pln)
% REDUCE_COLLIMATOR - Reduce dual-layer MLC to single-layer for MATRAD compatibility
%
% Usage: pln = reduce_collimator(pln)
%
% Takes the most restrictive leaf positions from dual-layer MLC systems
% (like Halcyon) and creates a single-layer MLC that MATRAD can process.
%
% Input:
%   pln - MATRAD plan structure with pln.propStf.beam containing MLC data
%
% Output:
%   pln - Modified plan structure with single-layer MLC data
%
% The "most restrictive" position is:
%   - For left bank (negative positions): maximum value (most closed)
%   - For right bank (positive positions): minimum value (most closed)

fprintf('\n========================================\n');
fprintf('DUAL-LAYER MLC REDUCTION\n');
fprintf('========================================\n');

% Check if MLC data exists
if ~isfield(pln, 'propStf') || ~isfield(pln.propStf, 'beam')
    fprintf('No beam data found in pln.propStf\n');
    fprintf('No collimator reduction needed.\n');
    return;
end

numBeams = length(pln.propStf.beam);
beamsProcessed = 0;
beamsReduced = 0;

for beamIdx = 1:numBeams
    
    if ~isfield(pln.propStf.beam(beamIdx), 'shape')
        continue;
    end
    
    beamsProcessed = beamsProcessed + 1;
    
    % Process X-direction MLC (most common)
    if isfield(pln.propStf.beam(beamIdx).shape, 'x')
        mlcX = pln.propStf.beam(beamIdx).shape.x;
        
        % Determine if this is dual-layer by checking for duplicate patterns
        % or by the number of leaves (should be even for dual-layer)
        numLeaves = length(mlcX) / 2;
        
        fprintf('\nBeam %d - X MLC Analysis:\n', beamIdx);
        fprintf('  Total leaf positions: %d\n', length(mlcX));
        fprintf('  Apparent leaf pairs: %.1f\n', numLeaves);
        
        % Check if we have an even number (potential dual-layer)
        if mod(numLeaves, 2) == 0
            % Assume dual-layer format: [layer1_left, layer2_left, ...] then [layer1_right, layer2_right, ...]
            % OR: [layer1_all_leaves, layer2_all_leaves]
            
            halfLength = length(mlcX) / 2;
            
            % Strategy 1: Assume interleaved layers (left bank, right bank format)
            % MLC format: [leaf1_left, leaf2_left, ...] [leaf1_right, leaf2_right, ...]
            % For dual layer: every other leaf is the same physical leaf, different layer
            
            % Check if we have a dual-layer pattern (every other leaf similar)
            leftBank = mlcX(1:halfLength);
            rightBank = mlcX(halfLength+1:end);
            
            if mod(length(leftBank), 2) == 0
                % Test for dual-layer pattern in left bank
                layer1_left = leftBank(1:2:end);
                layer2_left = leftBank(2:2:end);
                
                layer1_right = rightBank(1:2:end);
                layer2_right = rightBank(2:2:end);
                
                % Check if this looks like dual-layer (positions are relatively close)
                if length(layer1_left) > 1
                    meanDiff_left = mean(abs(layer1_left - layer2_left));
                    meanDiff_right = mean(abs(layer1_right - layer2_right));
                    
                    fprintf('  Mean position difference between layers:\n');
                    fprintf('    Left bank:  %.2f mm\n', meanDiff_left);
                    fprintf('    Right bank: %.2f mm\n', meanDiff_right);
                    
                    % If differences are small-ish (< 50mm), assume dual-layer
                    if meanDiff_left < 50 && meanDiff_right < 50
                        fprintf('  ✓ Detected dual-layer MLC pattern\n');
                        fprintf('  Reducing to single layer (most restrictive)...\n');
                        
                        % Take most restrictive: max for left bank (negative), min for right bank (positive)
                        reduced_left = max([layer1_left; layer2_left], [], 1);
                        reduced_right = min([layer1_right; layer2_right], [], 1);
                        
                        % Combine into single-layer format
                        mlcX_reduced = [reduced_left(:); reduced_right(:)];
                        
                        fprintf('  Original: %d positions\n', length(mlcX));
                        fprintf('  Reduced:  %d positions\n', length(mlcX_reduced));
                        fprintf('  Left bank sample - Layer 1: [%.2f, %.2f, %.2f...]\n', ...
                            layer1_left(1), layer1_left(2), layer1_left(min(3,end)));
                        fprintf('  Left bank sample - Layer 2: [%.2f, %.2f, %.2f...]\n', ...
                            layer2_left(1), layer2_left(2), layer2_left(min(3,end)));
                        fprintf('  Left bank sample - Reduced: [%.2f, %.2f, %.2f...]\n', ...
                            reduced_left(1), reduced_left(2), reduced_left(min(3,end)));
                        
                        % Update the structure
                        pln.propStf.beam(beamIdx).shape.x = mlcX_reduced;
                        beamsReduced = beamsReduced + 1;
                        
                    else
                        fprintf('  - Pattern does not match dual-layer (differences too large)\n');
                        fprintf('  - Keeping original MLC configuration\n');
                    end
                end
            else
                fprintf('  - Cannot reduce: odd number of leaves in bank\n');
            end
        else
            fprintf('  - Odd total leaf count - likely already single-layer\n');
        end
    end
    
    % Process Y-direction MLC (if it exists - uncommon)
    if isfield(pln.propStf.beam(beamIdx).shape, 'y')
        mlcY = pln.propStf.beam(beamIdx).shape.y;
        
        fprintf('\nBeam %d - Y MLC Analysis:\n', beamIdx);
        fprintf('  Total leaf positions: %d\n', length(mlcY));
        
        % Apply same logic as X direction
        numLeaves = length(mlcY) / 2;
        
        if mod(numLeaves, 2) == 0
            halfLength = length(mlcY) / 2;
            leftBank = mlcY(1:halfLength);
            rightBank = mlcY(halfLength+1:end);
            
            if mod(length(leftBank), 2) == 0
                layer1_left = leftBank(1:2:end);
                layer2_left = leftBank(2:2:end);
                layer1_right = rightBank(1:2:end);
                layer2_right = rightBank(2:2:end);
                
                if length(layer1_left) > 1
                    meanDiff_left = mean(abs(layer1_left - layer2_left));
                    meanDiff_right = mean(abs(layer1_right - layer2_right));
                    
                    if meanDiff_left < 50 && meanDiff_right < 50
                        fprintf('  ✓ Detected dual-layer MLC pattern in Y\n');
                        
                        reduced_left = max([layer1_left; layer2_left], [], 1);
                        reduced_right = min([layer1_right; layer2_right], [], 1);
                        mlcY_reduced = [reduced_left(:); reduced_right(:)];
                        
                        pln.propStf.beam(beamIdx).shape.y = mlcY_reduced;
                        fprintf('  Reduced Y MLC from %d to %d positions\n', ...
                            length(mlcY), length(mlcY_reduced));
                    end
                end
            end
        end
    end
end

fprintf('\n========================================\n');
fprintf('REDUCTION SUMMARY\n');
fprintf('========================================\n');
fprintf('Beams with MLC data: %d\n', beamsProcessed);
fprintf('Beams reduced to single-layer: %d\n', beamsReduced);

if beamsReduced > 0
    fprintf('\n✓ Dual-layer MLC successfully reduced to single-layer\n');
    fprintf('  This should improve MATRAD compatibility\n');
else
    fprintf('\nNo dual-layer patterns detected or reduction not needed\n');
end

fprintf('========================================\n\n');

end