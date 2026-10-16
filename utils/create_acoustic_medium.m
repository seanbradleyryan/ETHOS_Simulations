function medium = create_acoustic_medium(sct_resampled, config)
%CREATE_ACOUSTIC_MEDIUM Assign acoustic properties from CT Hounsfield Units
%
%   medium = create_acoustic_medium(sct_resampled, config)
%
%   Segments the CT volume into tissue types based on HU thresholds and
%   assigns density, speed of sound, attenuation, and Gruneisen coefficient
%   for each voxel. Supports four methods: uniform, threshold_1 (9-tissue
%   HU anchors, linearly interpolated between anchors), threshold_2
%   (4-tissue simplified bins), and continuous_1 (Yu et al. 2019 Table I:
%   piecewise-linear HU -> density -> sound speed, attenuation, and
%   Gamma = alpha_V c^2 / Cp).
%
%   INPUTS:
%       sct_resampled - Struct from step15_process_doses:
%           .cubeHU       - 3D array of Hounsfield Units
%           .spacing      - [dx, dy, dz] in mm
%           .dimensions   - [nx, ny, nz]
%           .bodyMask     - 3D logical (optional). Needed so the region
%                           outside the body stays water for sensor coupling:
%                           threshold_1 looks up every outside-body voxel as
%                           water; continuous_1 sets every outside-body voxel
%                           to its .outside_body water values; the
%                           threshold_2 air row is applied ONLY to low-HU
%                           voxels inside the body.
%       config - Configuration struct:
%           .gruneisen_method  - 'uniform' | 'threshold_1' | 'threshold_2' |
%                                'continuous_1'
%           .tissue_tables     - Tissue property lookup tables (from
%                                define_tissue_tables() in master pipeline).
%                                A threshold table may carry an '.air' row and
%                                '.air_hu_threshold' scalar (threshold_2 does);
%                                if so, in-body air pockets are assigned real
%                                air properties (gruneisen 0 -> no PA signal).
%
%   OUTPUTS:
%       medium - Struct with fields:
%           .density       - 3D density array (kg/m^3)
%           .sound_speed   - 3D speed of sound array (m/s)
%           .alpha_coeff   - 3D attenuation coefficient array (dB/MHz^y/cm)
%           .alpha_power   - Scalar power law exponent (table value for
%                            threshold_1 / continuous_1; mean across tissues
%                            for threshold_2)
%           .gruneisen     - 3D Gruneisen parameter array (dimensionless)
%           .grid_size     - [nx, ny, nz] dimensions
%           .method        - String indicating method used
%
%   EXAMPLE:
%       config.gruneisen_method = 'threshold_2';
%       config.tissue_tables = define_tissue_tables();
%       medium = create_acoustic_medium(sct_resampled, config);
%
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: run_single_field_simulation, define_tissue_tables

    %% ======================== INPUT VALIDATION ========================

    if ~isstruct(sct_resampled)
        error('create_acoustic_medium:InvalidInput', ...
            'sct_resampled must be a struct.');
    end

    if ~isfield(sct_resampled, 'cubeHU')
        error('create_acoustic_medium:MissingField', ...
            'sct_resampled must contain field ''cubeHU''.');
    end

    if ~isfield(config, 'tissue_tables') || isempty(config.tissue_tables)
        error('create_acoustic_medium:MissingConfig', ...
            'config.tissue_tables is required. Run define_tissue_tables() first.');
    end

    %% ======================== EXTRACT CT DATA ========================

    ctCube   = sct_resampled.cubeHU;
    gridSize = size(ctCube);

    % Default Gruneisen method
    if ~isfield(config, 'gruneisen_method') || isempty(config.gruneisen_method)
        config.gruneisen_method = 'threshold_2';
    end

    fprintf('    Creating acoustic medium (%s method)...\n', config.gruneisen_method);
    fprintf('      Grid size: [%d x %d x %d]\n', gridSize(1), gridSize(2), gridSize(3));
    fprintf('      HU range: [%.0f, %.0f]\n', min(ctCube(:)), max(ctCube(:)));

    %% ======================== ASSIGN PROPERTIES ========================

    switch lower(config.gruneisen_method)

        case 'uniform'
            %% --- Uniform: single tissue properties everywhere ---
            tables = config.tissue_tables.uniform;

            density       = tables.density     * ones(gridSize);
            soundSpeed    = tables.sound_speed * ones(gridSize);
            alphaCoeff    = tables.alpha_coeff * ones(gridSize);
            alphaPower    = tables.alpha_power;
            gruneisenArr  = tables.gruneisen   * ones(gridSize);

            fprintf('      Uniform medium: rho=%.0f kg/m^3, c=%.0f m/s, Gamma=%.2f\n', ...
                tables.density, tables.sound_speed, tables.gruneisen);

        case 'threshold_1'
            %% --- HU anchor table, linear interpolation between anchors ---
            tables   = config.tissue_tables.threshold_1;
            huLookup = double(ctCube);

            % Water override: voxels outside the body contour are looked up
            % as water, so the sensor-coupling region is a water bath.
            if isfield(sct_resampled, 'bodyMask') && ~isempty(sct_resampled.bodyMask) ...
                    && isequal(size(sct_resampled.bodyMask), gridSize)
                outsideBody = ~logical(sct_resampled.bodyMask);
                huLookup(outsideBody) = tables.outside_body_hu;
                fprintf('      Outside body -> water (HU %g): %d voxels (%.1f%%)\n', ...
                    tables.outside_body_hu, sum(outsideBody(:)), 100 * mean(outsideBody(:)));
            else
                fprintf('      [Water override] No matching bodyMask; outside-body voxels keep their CT HU.\n');
            end

            % HU outside the anchor range takes the end-anchor values
            huLookup = min(max(huLookup, tables.hu_points(1)), tables.hu_points(end));

            density      = interp_anchors(tables.hu_points, tables.density,     huLookup);
            soundSpeed   = interp_anchors(tables.hu_points, tables.sound_speed, huLookup);
            alphaCoeff   = interp_anchors(tables.hu_points, tables.alpha_coeff, huLookup);
            gruneisenArr = interp_anchors(tables.hu_points, tables.gruneisen,   huLookup);
            alphaPower   = tables.alpha_power;

        case 'threshold_2'
            %% --- Threshold-based tissue segmentation ---
            tables = config.tissue_tables.threshold_2;

            hu_bounds  = tables.hu_boundaries;
            numTissues = length(tables.tissue_names);

            % Initialize arrays
            density       = zeros(gridSize);
            soundSpeed    = zeros(gridSize);
            alphaCoeff    = zeros(gridSize);
            alphaPowerArr = zeros(gridSize);
            gruneisenArr  = zeros(gridSize);

            % Assign properties for each tissue bin
            % hu_boundaries defines bin edges: [edge_1, edge_2, ..., edge_N+1]
            % Tissue i covers: hu_boundaries(i) <= HU < hu_boundaries(i+1)
            for t = 1:numTissues
                hu_lo = hu_bounds(t);
                hu_hi = hu_bounds(t + 1);

                % Last bin includes upper edge (>= hu_lo)
                if t < numTissues
                    mask = (ctCube >= hu_lo) & (ctCube < hu_hi);
                else
                    mask = (ctCube >= hu_lo);
                end

                numVoxels = sum(mask(:));

                if numVoxels > 0
                    density(mask)       = tables.density(t);
                    soundSpeed(mask)    = tables.sound_speed(t);
                    alphaCoeff(mask)    = tables.alpha_coeff(t);
                    alphaPowerArr(mask) = tables.alpha_power(t);
                    gruneisenArr(mask)  = tables.gruneisen(t);

                    fprintf('      %-12s (HU [%6.0f, %6.0f)): %8d voxels (%5.1f%%)\n', ...
                        tables.tissue_names{t}, hu_lo, hu_hi, ...
                        numVoxels, 100 * numVoxels / numel(ctCube));
                end
            end

            % Handle unassigned voxels (default to water-like properties)
            unassigned = (density == 0);
            if any(unassigned(:))
                fprintf('      Unassigned: %d voxels -> water defaults\n', sum(unassigned(:)));
                density(unassigned)       = 1000;
                soundSpeed(unassigned)    = 1480;
                alphaCoeff(unassigned)    = 0.0022;
                alphaPowerArr(unassigned) = 2.0;
                gruneisenArr(unassigned)  = 0.11;
            end

            % --- Air pockets (inside body only) ---
            % When the table carries an air row, assign real air acoustics to
            % low-HU voxels INSIDE the body contour. Outside the body, low-HU
            % voxels stay water so the sensor-coupling bath is not turned into
            % air. Air gruneisen = 0, so air generates no forward PA signal.
            if isfield(tables, 'air') && isfield(tables, 'air_hu_threshold')
                if isfield(sct_resampled, 'bodyMask') && ~isempty(sct_resampled.bodyMask) ...
                        && isequal(size(sct_resampled.bodyMask), gridSize)
                    airMask = (ctCube < tables.air_hu_threshold) & logical(sct_resampled.bodyMask);
                    nAir    = sum(airMask(:));
                    if nAir > 0
                        density(airMask)       = tables.air.density;
                        soundSpeed(airMask)    = tables.air.sound_speed;
                        alphaCoeff(airMask)    = tables.air.alpha_coeff;
                        alphaPowerArr(airMask) = tables.air.alpha_power;
                        gruneisenArr(airMask)  = tables.air.gruneisen;
                        fprintf('      %-12s (HU < %5.0f, in body): %8d voxels (%5.1f%%)\n', ...
                            'Air', tables.air_hu_threshold, nAir, 100 * nAir / numel(ctCube));
                    end
                else
                    fprintf('      [Air] No matching bodyMask; low-HU voxels left as water (coupling).\n');
                end
            end

            % k-Wave uses a scalar alpha_power; compute weighted mean
            % over assigned voxels
            assignedMask = (density > 0);
            if any(assignedMask(:))
                alphaPower = mean(alphaPowerArr(assignedMask));
            else
                alphaPower = 1.1;  % safe default
            end

        case 'continuous_1'
            %% --- Yu et al. 2019 Table I: HU -> density -> c, alpha_V, Cp, alpha_0 ---
            tables = config.tissue_tables.continuous_1;

            % HU below min_hu would give negative density
            hu  = max(double(ctCube), tables.min_hu);
            rho = piecewise_linear(hu, tables.hu_to_rho);                  % g/cm^3

            soundSpeed = piecewise_linear(rho, tables.rho_to_c);           % m/s
            alphaV     = piecewise_linear(rho, tables.rho_to_alphaV);      % 1/K
            heatCap    = piecewise_linear(rho, tables.rho_to_Cp);          % J/(kg K)
            % Yu's air line dips slightly below 0 near HU -1000
            alphaCoeff = max(piecewise_linear(rho, tables.rho_to_alpha0), 0);  % dB/(MHz^y cm)

            gruneisenArr = alphaV .* soundSpeed.^2 ./ heatCap;           % Gamma = alpha_V c^2 / Cp
            density      = 1000 * rho;                                   % kg/m^3
            alphaPower   = tables.alpha_power;

            % Water override: voxels outside the body contour become a
            % water bath for sensor coupling.
            if isfield(sct_resampled, 'bodyMask') && ~isempty(sct_resampled.bodyMask) ...
                    && isequal(size(sct_resampled.bodyMask), gridSize)
                outsideBody = ~logical(sct_resampled.bodyMask);
                density(outsideBody)      = tables.outside_body.density;
                soundSpeed(outsideBody)   = tables.outside_body.sound_speed;
                alphaCoeff(outsideBody)   = tables.outside_body.alpha_coeff;
                gruneisenArr(outsideBody) = tables.outside_body.gruneisen;
                fprintf('      Outside body -> water: %d voxels (%.1f%%)\n', ...
                    sum(outsideBody(:)), 100 * mean(outsideBody(:)));
            else
                fprintf('      [Water override] No matching bodyMask; outside-body voxels keep their CT values.\n');
            end

        otherwise
            error('create_acoustic_medium:UnknownMethod', ...
                'Unknown gruneisen_method: ''%s''.\nValid options: ''uniform'', ''threshold_1'', ''threshold_2'', ''continuous_1''.', ...
                config.gruneisen_method);
    end

    %% ======================== NUMERICAL STABILITY ========================

    % Ensure minimum values to prevent division-by-zero or NaN propagation
    density    = max(density, 30);       % minimum 1 kg/m^3
    soundSpeed = max(soundSpeed, 100);  % minimum 100 m/s

    %% ======================== BUILD OUTPUT STRUCT ========================

    medium = struct();
    medium.density     = density;
    medium.sound_speed = soundSpeed;
    medium.alpha_coeff = alphaCoeff;
    medium.alpha_power = alphaPower;
    medium.gruneisen   = gruneisenArr;
    medium.grid_size   = gridSize;
    medium.method      = config.gruneisen_method;

    fprintf('      Density range:     [%.0f, %.0f] kg/m^3\n', ...
        min(density(:)), max(density(:)));
    fprintf('      Sound speed range: [%.0f, %.0f] m/s\n', ...
        min(soundSpeed(:)), max(soundSpeed(:)));
    fprintf('      Alpha power:       %.2f\n', alphaPower);
    fprintf('    Acoustic medium created.\n');
end


function values = interp_anchors(huPoints, anchorValues, hu)
%INTERP_ANCHORS Linearly interpolate one tissue property over the HU anchors.
%   Anchors with NaN (no sourced value) are skipped, so that HU range is
%   interpolated from the neighbouring anchors.
    known  = ~isnan(anchorValues);
    values = interp1(huPoints(known), anchorValues(known), hu, 'linear');
end


function y = piecewise_linear(x, segments)
%PIECEWISE_LINEAR Evaluate a curve made of straight-line pieces.
%   Each row of segments is [from, to, slope, intercept]; that row gives
%   y = slope * x + intercept for from <= x < to.
    y = zeros(size(x));
    for s = 1:size(segments, 1)
        inRange    = (x >= segments(s, 1)) & (x < segments(s, 2));
        y(inRange) = segments(s, 3) * x(inRange) + segments(s, 4);
    end
end
