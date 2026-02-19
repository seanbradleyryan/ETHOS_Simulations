function [sensor_mask, sensor_info] = determine_sensor_mask(sct_resampled, field_dose, beam_metadata, config)
%DETERMINE_SENSOR_MASK Place a flat 2D ultrasound array on the anterior abdomen
%
%   [sensor_mask, sensor_info] = determine_sensor_mask(sct_resampled, field_dose, beam_metadata, config)
%
%   PURPOSE:
%   Compute the 3D binary mask for a rigid, flat 10x10 cm ultrasound sensor
%   pressed against the patient's anterior abdomen. The sensor must avoid all
%   beam field projections (jaw openings) on the anterior surface. Everything
%   outside the patient body is water (no coupling concerns).
%
%   ALGORITHM:
%   1. Compute anterior surface height map from body mask (excluding couch).
%   2. Project each beam's jaw opening from isocenter onto the anterior surface.
%      Exclude beams whose projections do not intersect the anterior surface
%      (e.g., lateral/posterior beams with gantry > ~60 deg from AP).
%   3. Find the largest contiguous anterior surface area outside the exclusion
%      zone, as close as possible to the dose centroid's X-Z projection.
%   4. Place a flat planar sensor at a fixed Y index (the most anterior body
%      surface point within the chosen region, minus standoff).
%   5. Validate: sensor outside body, outside exclusion zone, within grid bounds.
%   6. Optionally partition the sensor into element patches for signal averaging.
%
%   INPUTS:
%       sct_resampled - Struct with CT resampled to dose grid:
%           .bodyMask       - 3D logical (true = inside body)
%           .couchMask      - 3D logical (true = couch region)
%           .origin         - [x, y, z] mm (DICOM patient coordinates)
%           .spacing        - [dx, dy, dz] mm
%           .dimensions     - [ny, nx, nz] (rows=Y, cols=X, slices=Z)
%       field_dose - Struct with:
%           .dose_Gy        - 3D dose array (Gy)
%           .gantry_angle   - Gantry angle for this field (degrees)
%           .origin         - [x, y, z] mm
%           .spacing        - [dx, dy, dz] mm
%           .dimensions     - [ny, nx, nz]
%       beam_metadata - Struct array (ALL beams in plan) with:
%           .beam_number    - Beam number
%           .gantry_angle   - Gantry angle (degrees)
%           .isocenter      - [x, y, z] mm (DICOM patient coords)
%           .jaw_x          - [x1, x2] mm at isocenter
%           .jaw_y          - [y1, y2] mm at isocenter
%       config - Struct with sensor parameters:
%           .sensor_size_cm       - [X, Z] physical sensor dims in cm (default: [10, 10])
%           .sensor_standoff_mm   - Gap between body surface and sensor (default: 5)
%           .element_size_mm      - Element patch size for averaging (default: [])
%           .jaw_margin_mm        - Extra margin around jaw projection (default: 10)
%           .sensor_placement     - Placement side: 'anterior' (default)
%           .pml_size             - PML thickness in voxels (default: 10)
%
%   OUTPUTS:
%       sensor_mask - 3D logical array (same size as dose grid), true at sensor voxels.
%       sensor_info - Struct with diagnostic fields:
%           .sensor_y_index              - Y index of the sensor plane
%           .sensor_x_range              - [x_start, x_end] voxel indices
%           .sensor_z_range              - [z_start, z_end] voxel indices
%           .sensor_center_mm            - [x, y, z] physical position of sensor center
%           .exclusion_zone              - 2D logical on anterior surface (X x Z)
%           .element_map                 - 2D array mapping sensor voxels to element index
%           .num_elements                - Number of signal-averaged elements
%           .surface_to_dose_distance_mm - Distance from sensor plane to dose centroid
%           .sensor_size_voxels          - [nx_sensor, nz_sensor] sensor extent in voxels
%           .placement_valid             - Boolean, true if placement passed all checks
%
%   COORDINATE SYSTEM:
%       Dose grid uses DICOM patient coordinates:
%           X = patient left-right (columns, dim 2)
%           Y = patient anterior-posterior (rows, dim 1). Lower index = more anterior.
%           Z = patient superior-inferior (slices, dim 3)
%       Array indexing: array(Y, X, Z)
%
%   NOTES:
%       - Source-to-axis distance (SAD) for Halcyon/ETHOS is 100 cm.
%       - Sensor must be entirely outside body mask and all jaw projections.
%       - Uses warnings (not errors) for non-fatal issues to support batch processing.
%       - Logging style consistent with pipeline (indented, with step labels).
%
%   EXAMPLE:
%       config.sensor_size_cm = [10, 10];
%       config.sensor_standoff_mm = 5;
%       config.element_size_mm = 2;
%       [mask, info] = determine_sensor_mask(sct_resampled, field_dose, beam_metadata, config);
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (bwconncomp, regionprops)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: run_single_field_simulation, step15_process_doses

%% ======================== CONFIG DEFAULTS ========================

sensor_size_cm    = get_field(config, 'sensor_size_cm', [10, 10]);
standoff_mm       = get_field(config, 'sensor_standoff_mm', 5);
element_size_mm   = get_field(config, 'element_size_mm', []);
jaw_margin_mm     = get_field(config, 'jaw_margin_mm', 10);
placement_side    = get_field(config, 'sensor_placement', 'anterior');
pml_size          = get_field(config, 'pml_size', 10);

SAD_mm = 1000;  % Halcyon/ETHOS source-to-axis distance: 100 cm

fprintf('        [Sensor] Placing %s sensor (%.0fx%.0f cm, standoff %.0f mm)\n', ...
    placement_side, sensor_size_cm(1), sensor_size_cm(2), standoff_mm);

%% ======================== EXTRACT GRID INFO ========================

% Grid dimensions: array is (Y, X, Z) = (rows, cols, slices)
grid_dims = size(sct_resampled.bodyMask);
Ny = grid_dims(1);  % rows = Y (anterior-posterior)
Nx = grid_dims(2);  % cols = X (left-right)
Nz = grid_dims(3);  % slices = Z (superior-inferior)

% Physical spacing
dx = sct_resampled.spacing(1);  % mm, X direction
dy = sct_resampled.spacing(2);  % mm, Y direction
dz = sct_resampled.spacing(3);  % mm, Z direction

% Grid origin in DICOM patient coordinates
origin = sct_resampled.origin(:)';  % [x, y, z] mm

% Sensor size in voxels
sensor_nx = round(sensor_size_cm(1) * 10 / dx);  % X extent in voxels
sensor_nz = round(sensor_size_cm(2) * 10 / dz);  % Z extent in voxels

fprintf('        [Sensor] Grid: [%d(Y) x %d(X) x %d(Z)], spacing: [%.2f, %.2f, %.2f] mm\n', ...
    Ny, Nx, Nz, dx, dy, dz);
fprintf('        [Sensor] Sensor size in voxels: %d(X) x %d(Z)\n', sensor_nx, sensor_nz);

%% ======================== STEP 1: ANTERIOR SURFACE MAP ========================

% Body mask excluding couch
body = sct_resampled.bodyMask & ~sct_resampled.couchMask;

% Anterior surface height map: for each (X, Z) column, find minimum Y index
% where body is true. Lower Y = more anterior.
% anterior_surface(x, z) = min Y index with body==true, or NaN if no body
anterior_surface = NaN(Nx, Nz);

for ix = 1:Nx
    for iz = 1:Nz
        col = squeeze(body(:, ix, iz));  % Y-column
        y_indices = find(col);
        if ~isempty(y_indices)
            anterior_surface(ix, iz) = min(y_indices);
        end
    end
end

% Valid surface points (where body exists in column)
surface_valid = ~isnan(anterior_surface);
num_surface_pts = sum(surface_valid(:));
fprintf('        [Sensor] Anterior surface: %d valid columns\n', num_surface_pts);

if num_surface_pts == 0
    warning('determine_sensor_mask:NoSurface', ...
        'No anterior body surface found. Returning empty sensor mask.');
    [sensor_mask, sensor_info] = empty_result(grid_dims);
    return;
end

%% ======================== STEP 2: BEAM FIELD EXCLUSION ZONE ========================

% Initialize exclusion zone on the anterior surface (X x Z)
exclusion_zone = false(Nx, Nz);

% Compute physical coordinates for each X, Z index
x_coords = origin(1) + (0:Nx-1) * dx;  % X physical positions (mm)
z_coords = origin(3) + (0:Nz-1) * dz;  % Z physical positions (mm)

% For each beam, project jaw opening onto the anterior surface
if ~isempty(beam_metadata) && isstruct(beam_metadata)
    for b = 1:length(beam_metadata)
        ga = mod(beam_metadata(b).gantry_angle, 360);
        
        % Only exclude beams that project onto the anterior surface.
        % AP beam (gantry ~0°): source is anterior, beam enters anteriorly.
        % Beams with gantry > ~60° from anterior don't constrain anterior sensor.
        % Anterior-facing beams: gantry in [0, 60] or [300, 360].
        if ~((ga >= 0 && ga <= 60) || (ga >= 300 && ga <= 360))
            fprintf('        [Sensor] Beam %d (gantry %.1f°): lateral/posterior, no anterior exclusion\n', ...
                beam_metadata(b).beam_number, ga);
            continue;
        end
        
        % Check for required fields
        if ~isfield(beam_metadata(b), 'isocenter') || isempty(beam_metadata(b).isocenter)
            warning('determine_sensor_mask:NoIsocenter', ...
                'Beam %d missing isocenter. Using dose centroid as fallback.', ...
                beam_metadata(b).beam_number);
            % Fallback: use dose centroid
            iso = compute_dose_centroid_mm(field_dose, origin, dx, dy, dz);
        else
            iso = beam_metadata(b).isocenter(:)';  % [x, y, z] mm
        end
        
        if ~isfield(beam_metadata(b), 'jaw_x') || isempty(beam_metadata(b).jaw_x)
            % Default Halcyon 10x10 cm jaws
            warning('determine_sensor_mask:NoJaws', ...
                'Beam %d missing jaw data. Using default [-50, 50] mm.', ...
                beam_metadata(b).beam_number);
            jaw_x = [-50, 50];
            jaw_y = [-50, 50];
        else
            jaw_x = beam_metadata(b).jaw_x;
            jaw_y = beam_metadata(b).jaw_y;
        end
        
        % Project jaw opening from source through isocenter onto anterior surface.
        % For gantry ~0° (AP beam): source is above (anterior to) patient.
        % Source position: isocenter + SAD in the beam direction.
        %
        % For gantry 0°: beam travels in +Y direction (anterior to posterior).
        % Source is at Y = iso_y - SAD (more anterior).
        % Jaw X limits define left-right field extent at isocenter.
        % Jaw Y limits define sup-inf field extent at isocenter.
        %
        % Project to the anterior surface Y plane:
        % For a general anterior beam, we project the jaw rectangle at isocenter
        % onto the mean anterior surface Y. Divergence factor = SSD / SAD,
        % where SSD = distance from source to anterior surface.
        
        % Mean anterior surface Y position (physical)
        valid_surface_y = anterior_surface(surface_valid);
        mean_surface_y_idx = round(median(valid_surface_y));
        mean_surface_y_mm = origin(2) + (mean_surface_y_idx - 1) * dy;
        
        % For gantry ~0°: source Y = iso_y - SAD
        % SSD = source_y to surface_y distance
        ga_rad = deg2rad(ga);
        
        % Source position relative to isocenter (IEC gantry convention)
        % Gantry 0°: beam travels +Y (ant→post), source at -Y from iso
        % Using simplified projection for near-AP beams:
        source_y = iso(2) - SAD_mm * cosd(ga);
        
        % Distance from source to anterior surface
        SSD = mean_surface_y_mm - source_y;
        
        if SSD <= 0
            % Surface is behind the source - shouldn't happen for AP
            fprintf('        [Sensor] Beam %d: SSD <= 0, skipping exclusion\n', ...
                beam_metadata(b).beam_number);
            continue;
        end
        
        % Divergence magnification factor from isocenter to surface
        SAD_to_surface = mean_surface_y_mm - source_y;
        SAD_to_iso = iso(2) - source_y;
        
        if SAD_to_iso <= 0
            mag = 1.0;  % Fallback
        else
            mag = SAD_to_surface / SAD_to_iso;
        end
        
        % Projected field extent at anterior surface (mm, centered on isocenter X,Z)
        field_x_min = iso(1) + jaw_x(1) * mag;
        field_x_max = iso(1) + jaw_x(2) * mag;
        field_z_min = iso(3) + jaw_y(1) * mag;  % jaw_y maps to Z (sup-inf)
        field_z_max = iso(3) + jaw_y(2) * mag;
        
        % Add margin
        field_x_min = field_x_min - jaw_margin_mm;
        field_x_max = field_x_max + jaw_margin_mm;
        field_z_min = field_z_min - jaw_margin_mm;
        field_z_max = field_z_max + jaw_margin_mm;
        
        % Convert to voxel indices
        ix_min = max(1,  floor((field_x_min - origin(1)) / dx) + 1);
        ix_max = min(Nx, ceil( (field_x_max - origin(1)) / dx) + 1);
        iz_min = max(1,  floor((field_z_min - origin(3)) / dz) + 1);
        iz_max = min(Nz, ceil( (field_z_max - origin(3)) / dz) + 1);
        
        % Mark exclusion zone
        exclusion_zone(ix_min:ix_max, iz_min:iz_max) = true;
        
        fprintf('        [Sensor] Beam %d (gantry %.1f°): exclusion X=[%d,%d], Z=[%d,%d] (mag=%.2f)\n', ...
            beam_metadata(b).beam_number, ga, ix_min, ix_max, iz_min, iz_max, mag);
    end
else
    warning('determine_sensor_mask:NoBeamMetadata', ...
        'No beam metadata provided. Sensor placed without beam exclusion.');
end

fprintf('        [Sensor] Exclusion zone: %d voxels (%.1f%% of surface)\n', ...
    sum(exclusion_zone(:) & surface_valid(:)), ...
    100 * sum(exclusion_zone(:) & surface_valid(:)) / max(1, num_surface_pts));

%% ======================== STEP 3: FIND SENSOR PLACEMENT REGION ========================

% Available region: on body surface AND not in exclusion zone
available = surface_valid & ~exclusion_zone;

% Also exclude PML boundary regions
pml_margin_x = pml_size + 2;  % Extra safety margin beyond PML
pml_margin_z = pml_size + 2;
available(1:pml_margin_x, :) = false;
available(end-pml_margin_x+1:end, :) = false;
available(:, 1:pml_margin_z) = false;
available(:, end-pml_margin_z+1:end) = false;

if sum(available(:)) < sensor_nx * sensor_nz
    warning('determine_sensor_mask:InsufficientSpace', ...
        'Available anterior surface (%d voxels) may be too small for sensor (%d voxels).', ...
        sum(available(:)), sensor_nx * sensor_nz);
end

% Compute dose centroid in X-Z (for proximity targeting)
dose_centroid_mm = compute_dose_centroid_mm(field_dose, origin, dx, dy, dz);
dose_centroid_ix = round((dose_centroid_mm(1) - origin(1)) / dx) + 1;
dose_centroid_iz = round((dose_centroid_mm(3) - origin(3)) / dz) + 1;

fprintf('        [Sensor] Dose centroid (voxel): X=%d, Z=%d\n', dose_centroid_ix, dose_centroid_iz);

% Strategy: sweep candidate sensor rectangles across the available region,
% find the position closest to the dose centroid where the full sensor fits.
best_dist = Inf;
best_ix_start = [];
best_iz_start = [];

% Half-sensor extents for centering search
half_nx = floor(sensor_nx / 2);
half_nz = floor(sensor_nz / 2);

% Search grid: candidate top-left corners for the sensor rectangle
for ix_start = 1:(Nx - sensor_nx + 1)
    for iz_start = 1:(Nz - sensor_nz + 1)
        ix_end = ix_start + sensor_nx - 1;
        iz_end = iz_start + sensor_nz - 1;
        
        % Check if entire sensor rectangle is in available region
        patch = available(ix_start:ix_end, iz_start:iz_end);
        if all(patch(:))
            % Compute center of this candidate
            cx = ix_start + half_nx;
            cz = iz_start + half_nz;
            
            % Distance to dose centroid in X-Z plane (voxel units)
            dist = sqrt((cx - dose_centroid_ix)^2 + (cz - dose_centroid_iz)^2);
            
            if dist < best_dist
                best_dist = dist;
                best_ix_start = ix_start;
                best_iz_start = iz_start;
            end
        end
    end
end

% Fallback: if no full-sensor rectangle fits, try shrinking
if isempty(best_ix_start)
    warning('determine_sensor_mask:ShrinkSensor', ...
        'Full sensor does not fit in available region. Attempting to shrink.');
    
    % Try progressively smaller sizes (90%, 80%, ..., 50%)
    for shrink_pct = [0.9, 0.8, 0.7, 0.6, 0.5]
        trial_nx = round(sensor_nx * shrink_pct);
        trial_nz = round(sensor_nz * shrink_pct);
        trial_half_nx = floor(trial_nx / 2);
        trial_half_nz = floor(trial_nz / 2);
        
        for ix_start = 1:(Nx - trial_nx + 1)
            for iz_start = 1:(Nz - trial_nz + 1)
                ix_end = ix_start + trial_nx - 1;
                iz_end = iz_start + trial_nz - 1;
                
                patch = available(ix_start:ix_end, iz_start:iz_end);
                if all(patch(:))
                    cx = ix_start + trial_half_nx;
                    cz = iz_start + trial_half_nz;
                    dist = sqrt((cx - dose_centroid_ix)^2 + (cz - dose_centroid_iz)^2);
                    
                    if dist < best_dist
                        best_dist = dist;
                        best_ix_start = ix_start;
                        best_iz_start = iz_start;
                        sensor_nx = trial_nx;
                        sensor_nz = trial_nz;
                    end
                end
            end
        end
        
        if ~isempty(best_ix_start)
            fprintf('        [Sensor] Sensor shrunk to %.0f%% (now %d x %d voxels)\n', ...
                shrink_pct * 100, sensor_nx, sensor_nz);
            break;
        end
    end
end

if isempty(best_ix_start)
    warning('determine_sensor_mask:NoPlacement', ...
        'Could not find any valid sensor placement. Returning empty mask.');
    [sensor_mask, sensor_info] = empty_result(grid_dims);
    sensor_info.exclusion_zone = exclusion_zone;
    return;
end

% Final sensor X and Z ranges (voxel indices in the X and Z dimensions)
sensor_x_range = [best_ix_start, best_ix_start + sensor_nx - 1];
sensor_z_range = [best_iz_start, best_iz_start + sensor_nz - 1];

fprintf('        [Sensor] Placement: X=[%d,%d], Z=[%d,%d] (dist to dose: %.1f voxels)\n', ...
    sensor_x_range(1), sensor_x_range(2), sensor_z_range(1), sensor_z_range(2), best_dist);

%% ======================== STEP 4: DETERMINE SENSOR Y INDEX ========================

% Find the most anterior (minimum Y) body surface point within the sensor region
surface_patch = anterior_surface(sensor_x_range(1):sensor_x_range(2), ...
                                  sensor_z_range(1):sensor_z_range(2));

min_anterior_y = min(surface_patch(:));

% Place sensor at standoff distance anterior to the body surface
standoff_voxels = ceil(standoff_mm / dy);
sensor_y_index = max(1, min_anterior_y - standoff_voxels);

% Ensure sensor Y is within grid bounds (accounting for PML)
if sensor_y_index <= pml_size
    warning('determine_sensor_mask:SensorNearPML', ...
        'Sensor Y index (%d) is within PML region (size %d). Adjusting.', ...
        sensor_y_index, pml_size);
    sensor_y_index = pml_size + 1;
end

% Physical position of sensor center
sensor_center_x = origin(1) + (mean(sensor_x_range) - 1) * dx;
sensor_center_y = origin(2) + (sensor_y_index - 1) * dy;
sensor_center_z = origin(3) + (mean(sensor_z_range) - 1) * dz;

fprintf('        [Sensor] Y index: %d (surface min Y: %d, standoff: %d voxels)\n', ...
    sensor_y_index, min_anterior_y, standoff_voxels);
fprintf('        [Sensor] Center (mm): [%.1f, %.1f, %.1f]\n', ...
    sensor_center_x, sensor_center_y, sensor_center_z);

%% ======================== STEP 5: CREATE SENSOR MASK ========================

% Sensor is a flat rectangle at a single Y index spanning the X-Z range.
% Array indexing: mask(Y, X, Z)
sensor_mask = false(grid_dims);
sensor_mask(sensor_y_index, ...
            sensor_x_range(1):sensor_x_range(2), ...
            sensor_z_range(1):sensor_z_range(2)) = true;

num_sensor_voxels = sum(sensor_mask(:));
fprintf('        [Sensor] Mask created: %d voxels\n', num_sensor_voxels);

%% ======================== STEP 6: VALIDATION ========================

placement_valid = true;

% Check 1: No sensor voxels inside body
overlap_body = sensor_mask & body;
if any(overlap_body(:))
    num_overlap = sum(overlap_body(:));
    warning('determine_sensor_mask:BodyOverlap', ...
        'Sensor overlaps body mask at %d voxels. Removing overlapping voxels.', num_overlap);
    sensor_mask = sensor_mask & ~body;
    placement_valid = false;
end

% Check 2: No sensor voxels in exclusion zone
% Map exclusion zone (X, Z) to sensor Y plane
for ix = sensor_x_range(1):sensor_x_range(2)
    for iz = sensor_z_range(1):sensor_z_range(2)
        if exclusion_zone(ix, iz)
            if sensor_mask(sensor_y_index, ix, iz)
                warning('determine_sensor_mask:ExclusionOverlap', ...
                    'Sensor voxel at X=%d, Z=%d overlaps exclusion zone. Removing.', ix, iz);
                sensor_mask(sensor_y_index, ix, iz) = false;
                placement_valid = false;
            end
        end
    end
end

% Check 3: Sensor within grid bounds (accounting for PML)
if sensor_x_range(1) <= pml_size || sensor_x_range(2) > Nx - pml_size || ...
   sensor_z_range(1) <= pml_size || sensor_z_range(2) > Nz - pml_size || ...
   sensor_y_index <= pml_size
    warning('determine_sensor_mask:PMLOverlap', ...
        'Sensor extends into PML region (PML size: %d). Results may be affected.', pml_size);
    placement_valid = false;
end

num_sensor_voxels_final = sum(sensor_mask(:));
fprintf('        [Sensor] Final sensor: %d voxels (valid: %s)\n', ...
    num_sensor_voxels_final, mat2str(placement_valid));

%% ======================== STEP 7: ELEMENT GROUPINGS ========================

element_map = [];
num_elements = 0;

if ~isempty(element_size_mm) && element_size_mm > 0
    % Partition sensor into patches of element_size_mm x element_size_mm
    elem_nx = max(1, floor(element_size_mm / dx));  % Element size in X voxels
    elem_nz = max(1, floor(element_size_mm / dz));  % Element size in Z voxels
    
    % Create element map for the sensor rectangle (local coordinates)
    local_nx = sensor_x_range(2) - sensor_x_range(1) + 1;
    local_nz = sensor_z_range(2) - sensor_z_range(1) + 1;
    element_map = zeros(local_nx, local_nz);
    
    elem_idx = 0;
    for ex = 1:elem_nx:local_nx
        for ez = 1:elem_nz:local_nz
            elem_idx = elem_idx + 1;
            ex_end = min(ex + elem_nx - 1, local_nx);
            ez_end = min(ez + elem_nz - 1, local_nz);
            element_map(ex:ex_end, ez:ez_end) = elem_idx;
        end
    end
    
    num_elements = elem_idx;
    
    % Zero out element assignments for voxels that were removed during validation
    for local_ix = 1:local_nx
        for local_iz = 1:local_nz
            global_ix = sensor_x_range(1) + local_ix - 1;
            global_iz = sensor_z_range(1) + local_iz - 1;
            if ~sensor_mask(sensor_y_index, global_ix, global_iz)
                element_map(local_ix, local_iz) = 0;
            end
        end
    end
    
    fprintf('        [Sensor] Element grouping: %d elements (%.1f mm patches, %dx%d voxels each)\n', ...
        num_elements, element_size_mm, elem_nx, elem_nz);
end

%% ======================== COMPUTE DIAGNOSTICS ========================

% Distance from sensor plane to dose centroid
surface_to_dose_distance_mm = abs(dose_centroid_mm(2) - sensor_center_y);

%% ======================== PACK OUTPUT ========================

sensor_info = struct();
sensor_info.sensor_y_index              = sensor_y_index;
sensor_info.sensor_x_range              = sensor_x_range;
sensor_info.sensor_z_range              = sensor_z_range;
sensor_info.sensor_center_mm            = [sensor_center_x, sensor_center_y, sensor_center_z];
sensor_info.exclusion_zone              = exclusion_zone;
sensor_info.element_map                 = element_map;
sensor_info.num_elements                = num_elements;
sensor_info.surface_to_dose_distance_mm = surface_to_dose_distance_mm;
sensor_info.sensor_size_voxels          = [sensor_nx, sensor_nz];
sensor_info.placement_valid             = placement_valid;
sensor_info.num_sensor_voxels           = num_sensor_voxels_final;
sensor_info.anterior_surface            = anterior_surface;
sensor_info.standoff_voxels             = standoff_voxels;
sensor_info.gantry_angle                = field_dose.gantry_angle;

fprintf('        [Sensor] Surface-to-dose distance: %.1f mm\n', surface_to_dose_distance_mm);
fprintf('        [Sensor] Placement complete.\n');

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function centroid_mm = compute_dose_centroid_mm(field_dose, origin, dx, dy, dz)
%COMPUTE_DOSE_CENTROID_MM Compute the physical centroid of the dose distribution
%
%   Returns [x, y, z] in mm (DICOM patient coordinates).

    dose = field_dose.dose_Gy;
    dose(dose < 0.01 * max(dose(:))) = 0;  % Threshold at 1% of max
    
    total_dose = sum(dose(:));
    if total_dose == 0
        % Fallback: center of the grid
        dims = size(dose);
        centroid_mm = origin + ([dims(2), dims(1), dims(3)] - 1) / 2 .* [dx, dy, dz];
        return;
    end
    
    [Ny, Nx, Nz] = size(dose);
    
    % Create coordinate grids (physical mm)
    x_vec = origin(1) + (0:Nx-1) * dx;
    y_vec = origin(2) + (0:Ny-1) * dy;
    z_vec = origin(3) + (0:Nz-1) * dz;
    
    % Dose-weighted centroid
    % Sum over Y, Z to get X projection; etc.
    dose_x = squeeze(sum(sum(dose, 1), 3));  % [1 x Nx] -> [Nx x 1]
    dose_y = squeeze(sum(sum(dose, 2), 3));  % [Ny x 1]
    dose_z = squeeze(sum(sum(dose, 1), 2));  % [1 x 1 x Nz] -> [Nz x 1]
    
    cx = sum(x_vec(:) .* dose_x(:)) / total_dose;
    cy = sum(y_vec(:) .* dose_y(:)) / total_dose;
    cz = sum(z_vec(:) .* dose_z(:)) / total_dose;
    
    centroid_mm = [cx, cy, cz];
end


function [sensor_mask, sensor_info] = empty_result(grid_dims)
%EMPTY_RESULT Return empty/default sensor mask and info struct

    sensor_mask = false(grid_dims);
    sensor_info = struct();
    sensor_info.sensor_y_index              = 0;
    sensor_info.sensor_x_range              = [0, 0];
    sensor_info.sensor_z_range              = [0, 0];
    sensor_info.sensor_center_mm            = [0, 0, 0];
    sensor_info.exclusion_zone              = false(grid_dims(2), grid_dims(3));
    sensor_info.element_map                 = [];
    sensor_info.num_elements                = 0;
    sensor_info.surface_to_dose_distance_mm = 0;
    sensor_info.sensor_size_voxels          = [0, 0];
    sensor_info.placement_valid             = false;
    sensor_info.num_sensor_voxels           = 0;
end


function val = get_field(s, name, default)
%GET_FIELD Retrieve struct field with fallback to default
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default;
    end
end
