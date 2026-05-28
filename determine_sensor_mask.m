function [sensor_mask, sensor_info] = determine_sensor_mask(sct_resampled, field_dose, beam_metadata, config)
%DETERMINE_SENSOR_MASK Place a flat 2D ultrasound array on the anterior abdomen
%
%   [sensor_mask, sensor_info] = determine_sensor_mask(sct_resampled, field_dose, beam_metadata, config)
%
%   PURPOSE:
%   Compute the 3D binary mask for a rigid, flat 2D ultrasound array
%   pressed against the patient's anterior abdomen. The array is modeled as
%   a sparse grid of discrete piezoelectric elements (default 32x32) with a
%   physical pitch, active element size, and kerf (dead) gap between elements.
%   Only voxels falling within active element footprints become sensor points;
%   kerf regions are empty. The aperture footprint is derived from
%   elements_per_side * element_pitch_mm, not from a separate sensor_size_cm.
%   The sensor must avoid all beam field projections (jaw openings) on the
%   anterior surface. Everything outside the patient body is water
%   (no coupling concerns).
%
%   ALGORITHM:
%   1. Compute anterior surface height map from body mask (excluding couch).
%   2. Project each beam's jaw opening from isocenter onto the anterior surface.
%      Exclude beams whose projections do not intersect the anterior surface
%      (e.g., lateral/posterior beams with gantry > ~60 deg from AP).
%   3. Find the largest contiguous anterior surface area outside the exclusion
%      zone, as close as possible to the dose centroid's X-Z projection.
%      Aperture footprint = elements_per_side * element_pitch_mm.
%   4. Place a flat planar sensor at a fixed Y index (the most anterior body
%      surface point within the chosen region, minus standoff).
%   5. Build sparse mask: for each of N x N element centers, mark grid voxels
%      within the element's active footprint (element_size_mm x element_size_mm).
%      Kerf voxels remain unmarked. element_map encodes element index per voxel.
%   6. Validate: sensor outside body, outside exclusion zone, within grid bounds.
%   7. Aliasing check: if grid is too coarse to resolve all elements uniquely,
%      warn-and-degrade (continue with reduced effective element count).
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
%           .elements_per_side    - N for an NxN element array (default: 32)
%           .element_pitch_mm     - Center-to-center pitch in mm (default: 3.65)
%           .element_size_mm      - Active element width in mm (default: 2.43).
%                                   Kerf is derived: kerf = pitch - size.
%           .sensor_standoff_mm   - Gap between body surface and sensor (default: 5)
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
%           .element_map                 - 2D array (X x Z over aperture footprint)
%                                          mapping each voxel to element index;
%                                          0 marks kerf/inactive voxels.
%           .num_elements                - Effective number of resolved elements
%                                          (= effective_element_count, may be
%                                          < elements_per_side^2 if aliased).
%           .element_positions_mm        - [N x 3] physical centers of all elements.
%           .element_pitch_mm            - Echoed pitch (mm).
%           .element_size_mm             - Echoed active element width (mm).
%           .kerf_mm                     - Derived kerf = pitch - size (mm).
%           .elements_per_side           - Echoed N (default 32).
%           .effective_element_count     - Elements with >= 1 assigned voxel.
%           .aliased_element_count       - Elements that lost their voxel.
%           .voxels_per_element_stats    - [min, median, max] voxels per element.
%           .aperture_mm                 - Derived total aperture (mm).
%           .fill_factor_actual          - Actual mask fill within aperture
%                                          rectangle (compare to spec 0.44).
%           .surface_to_dose_distance_mm - Distance from sensor plane to dose centroid
%           .sensor_size_voxels          - [nx_sensor, nz_sensor] aperture footprint extent in voxels
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
%       config.elements_per_side = 32;
%       config.element_pitch_mm  = 3.65;
%       config.element_size_mm   = 2.43;
%       config.sensor_standoff_mm = 5;
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

elements_per_side = get_field(config, 'elements_per_side', 32);
element_pitch_mm  = get_field(config, 'element_pitch_mm', 3.65);
element_size_mm   = get_field(config, 'element_size_mm', 2.43);
standoff_mm       = get_field(config, 'sensor_standoff_mm', 5);
jaw_margin_mm     = get_field(config, 'jaw_margin_mm', 10);
placement_side    = get_field(config, 'sensor_placement', 'anterior');
pml_size          = get_field(config, 'pml_size', 10);
pml_size = 1; % Hardcoded because script logic assumes pml inside. 

% Kerf is derived; never accepted from config.
kerf_mm = element_pitch_mm - element_size_mm;
if kerf_mm < 0
    error('determine_sensor_mask:InvalidGeometry', ...
        'element_size_mm (%.3f) cannot exceed element_pitch_mm (%.3f).', ...
        element_size_mm, element_pitch_mm);
end

% Total aperture footprint (mm) derived from N x pitch.
aperture_mm = elements_per_side * element_pitch_mm;

SAD_mm = 1000;  % Halcyon/ETHOS source-to-axis distance: 100 cm

fprintf('        [Sensor] Placing %s array: %dx%d elements, pitch %.2f mm, size %.2f mm, kerf %.2f mm\n', ...
    placement_side, elements_per_side, elements_per_side, ...
    element_pitch_mm, element_size_mm, kerf_mm);
fprintf('        [Sensor] Aperture %.1f mm, standoff %.0f mm, fill factor (spec) %.0f%%\n', ...
    aperture_mm, standoff_mm, ...
    100 * (element_size_mm / element_pitch_mm)^2);

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

% Aperture footprint in voxels (for the placement search)
sensor_nx = ceil(aperture_mm / dx);  % X extent in voxels
sensor_nz = ceil(aperture_mm / dz);  % Z extent in voxels

fprintf('        [Sensor] Grid: [%d(Y) x %d(X) x %d(Z)], spacing: [%.2f, %.2f, %.2f] mm\n', ...
    Ny, Nx, Nz, dx, dy, dz);
fprintf('        [Sensor] Aperture footprint in voxels: %d(X) x %d(Z)\n', sensor_nx, sensor_nz);
if dx > element_pitch_mm || dz > element_pitch_mm
    warning('determine_sensor_mask:GridCoarserThanPitch', ...
        'Grid spacing (dx=%.2f, dz=%.2f mm) exceeds element pitch (%.2f mm). Severe element aliasing expected.', ...
        dx, dz, element_pitch_mm);
end

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
% beam_metadata = [];
if ~isempty(beam_metadata) && isstruct(beam_metadata)
    for b = 1:length(beam_metadata)
        ga = mod(beam_metadata(b).gantry_angle, 360);
        
        % Only exclude beams that project onto the anterior surface.
        % AP beam (gantry ~0): source is anterior, beam enters anteriorly.
        % Beams with gantry > ~60 from anterior don't constrain anterior sensor.
        % Anterior-facing beams: gantry in [0, 60] or [300, 360].
        if ~((ga >= 0 && ga <= 60) || (ga >= 300 && ga <= 360))
            fprintf('        [Sensor] Beam %d (gantry %.1f): lateral/posterior, no anterior exclusion\n', ...
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
        % For gantry ~0 (AP beam): source is above (anterior to) patient.
        % Source position: isocenter + SAD in the beam direction.
        %
        % For gantry 0: beam travels in +Y direction (anterior to posterior).
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
        
        % For gantry ~0: source Y = iso_y - SAD
        % SSD = source_y to surface_y distance
        ga_rad = deg2rad(ga);
        
        % Source position relative to isocenter (IEC gantry convention)
        % Gantry 0: beam travels +Y (antpost), source at -Y from iso
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
        
        fprintf('        [Sensor] Beam %d (gantry %.1f): exclusion X=[%d,%d], Z=[%d,%d] (mag=%.2f)\n', ...
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

% Track grid-expansion padding (set by the fallback path below).
% Caller must apply matching water-filled padding to medium/p0/etc.
grid_pad_x_pre  = 0; grid_pad_x_post = 0;
grid_pad_y_pre  = 0; grid_pad_y_post = 0;
grid_pad_z_pre  = 0; grid_pad_z_post = 0;
grid_was_expanded = false;

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

% Fallback: expand the computational grid (in X and/or Z) to fit the sensor
% outside the exclusion zone. New voxels are treated as water — body=false,
% exclusion=false, no anterior surface. The caller must apply matching
% water-filled padding to medium/p0 using sensor_info.grid_pad.
% Inferior Z placement is preferred over superior on ties (per user).
if isempty(best_ix_start)
    fprintf('        [Sensor] No placement in original grid; expanding grid past exclusion zone.\n');

    % --- X placement: center on dose centroid; pad X only if sensor wider than grid ---
    ix_start_target = dose_centroid_ix - half_nx;
    ix_start_min_orig = pml_margin_x + 1;
    ix_start_max_orig = Nx - pml_margin_x - sensor_nx + 1;

    if ix_start_max_orig >= ix_start_min_orig
        ix_start_orig = max(ix_start_min_orig, min(ix_start_max_orig, ix_start_target));
    else
        % Grid narrower than sensor+PML; pad X symmetrically around target.
        extra_x = (sensor_nx + 2 * pml_margin_x) - Nx;
        grid_pad_x_pre  = floor(extra_x / 2);
        grid_pad_x_post = extra_x - grid_pad_x_pre;
        % In pre-expansion coords, the sensor's left edge sits at
        % (pml_margin_x+1) - grid_pad_x_pre (may be <= 0; that range is water).
        ix_start_orig = (pml_margin_x + 1) - grid_pad_x_pre;
    end
    ix_end_orig = ix_start_orig + sensor_nx - 1;

    % --- Find exclusion Z extent within the sensor's X strip (clipped to original grid) ---
    ix_lo_clip = max(1, ix_start_orig);
    ix_hi_clip = min(Nx, ix_end_orig);
    if ix_lo_clip > ix_hi_clip
        excl_z_indices = [];
    else
        excl_strip = exclusion_zone(ix_lo_clip:ix_hi_clip, :);
        excl_z_indices = find(any(excl_strip, 1));
    end

    if isempty(excl_z_indices)
        % No exclusion under the sensor X strip — center on dose Z, no Z pad needed.
        iz_start_orig   = max(pml_margin_z + 1, ...
                          min(Nz - pml_margin_z - sensor_nz + 1, ...
                              dose_centroid_iz - half_nz));
        grid_pad_z_pre  = 0;
        grid_pad_z_post = 0;
        fprintf('        [Sensor] No exclusion in X strip; centered Z=%d (no Z pad)\n', iz_start_orig);
    else
        excl_z_max = max(excl_z_indices);
        excl_z_min = min(excl_z_indices);

        % Inferior candidate: butt up just past excl_z_max
        iz_start_inf  = excl_z_max + 1;
        iz_end_inf    = iz_start_inf + sensor_nz - 1;
        pad_z_post_inf = max(0, iz_end_inf - (Nz - pml_margin_z));

        % Superior candidate: butt up just before excl_z_min
        iz_end_sup    = excl_z_min - 1;
        iz_start_sup  = iz_end_sup - sensor_nz + 1;
        pad_z_pre_sup = max(0, (pml_margin_z + 1) - iz_start_sup);

        if pad_z_post_inf <= pad_z_pre_sup
            iz_start_orig   = iz_start_inf;
            grid_pad_z_post = pad_z_post_inf;
            fprintf('        [Sensor] Inferior to exclusion (Z=%d); padding Z+ by %d voxels (water)\n', ...
                iz_start_orig, grid_pad_z_post);
        else
            iz_start_orig   = iz_start_sup;
            grid_pad_z_pre  = pad_z_pre_sup;
            fprintf('        [Sensor] Superior to exclusion (Z=%d orig); padding Z- by %d voxels (water)\n', ...
                iz_start_orig, grid_pad_z_pre);
        end
    end

    % --- Pad anterior_surface / body / surface_valid / exclusion_zone with water defaults ---
    new_Nx = Nx + grid_pad_x_pre + grid_pad_x_post;
    new_Nz = Nz + grid_pad_z_pre + grid_pad_z_post;

    body_expanded = false(Ny, new_Nx, new_Nz);
    body_expanded(:, grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = body;
    body = body_expanded;

    surface_valid_expanded = false(new_Nx, new_Nz);
    surface_valid_expanded(grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = surface_valid;
    surface_valid = surface_valid_expanded;

    exclusion_zone_expanded = false(new_Nx, new_Nz);
    exclusion_zone_expanded(grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = exclusion_zone;
    exclusion_zone = exclusion_zone_expanded;

    anterior_surface_expanded = NaN(new_Nx, new_Nz);
    anterior_surface_expanded(grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = anterior_surface;
    anterior_surface = anterior_surface_expanded;

    % Update grid dims to expanded values; downstream mask building uses these.
    Nx_pre_expansion = Nx;
    Nz_pre_expansion = Nz;
    Nx = new_Nx;
    Nz = new_Nz;
    grid_dims = [Ny, Nx, Nz];

    % Convert pre-expansion placement coords to expanded coords.
    best_ix_start = ix_start_orig + grid_pad_x_pre;
    best_iz_start = iz_start_orig + grid_pad_z_pre;
    best_dist     = sqrt( (best_ix_start + half_nx - (dose_centroid_ix + grid_pad_x_pre))^2 + ...
                          (best_iz_start + half_nz - (dose_centroid_iz + grid_pad_z_pre))^2 );

    % Shift dose centroid voxel index into expanded coords (for any downstream use).
    dose_centroid_ix = dose_centroid_ix + grid_pad_x_pre;
    dose_centroid_iz = dose_centroid_iz + grid_pad_z_pre;

    grid_was_expanded = (grid_pad_x_pre + grid_pad_x_post + grid_pad_z_pre + grid_pad_z_post) > 0;

    fprintf('        [Sensor] Grid expanded: [Nx %d, Nz %d] -> [Nx %d, Nz %d]; sensor at X=[%d,%d], Z=[%d,%d]\n', ...
        Nx_pre_expansion, Nz_pre_expansion, Nx, Nz, ...
        best_ix_start, best_ix_start + sensor_nx - 1, ...
        best_iz_start, best_iz_start + sensor_nz - 1);
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

% Defensive: confirm the aperture footprint can host the configured element
% grid at the requested pitch. With grid-expansion, this should always pass;
% kept as a guard in case sensor_nx/sensor_nz were modified elsewhere.
actual_aperture_x_mm = sensor_nx * dx;
actual_aperture_z_mm = sensor_nz * dz;
N_fit_x = max(1, floor(actual_aperture_x_mm / element_pitch_mm));
N_fit_z = max(1, floor(actual_aperture_z_mm / element_pitch_mm));
N_fit   = min(N_fit_x, N_fit_z);
if N_fit < elements_per_side
    warning('determine_sensor_mask:ApertureShrunk', ...
        'Available aperture (%.1f x %.1f mm) fits %dx%d elements at pitch %.2f mm; reducing from %dx%d.', ...
        actual_aperture_x_mm, actual_aperture_z_mm, N_fit, N_fit, ...
        element_pitch_mm, elements_per_side, elements_per_side);
    elements_per_side = N_fit;
    aperture_mm = elements_per_side * element_pitch_mm;
end

%% ======================== STEP 4: DETERMINE SENSOR Y INDEX ========================

% Find the most anterior (minimum Y) body surface point within the sensor region
surface_patch = anterior_surface(sensor_x_range(1):sensor_x_range(2), ...
                                  sensor_z_range(1):sensor_z_range(2));

% Footprint may overhang the body (grid-expansion placement places the
% sensor in water past the body's Z extent), so NaN columns are expected.
% Take the min over real surface points only; fall back to the global
% median anterior surface Y if no body sits under the footprint at all.
min_anterior_y = min(surface_patch(:), [], 'omitnan');
if isnan(min_anterior_y)
    valid_surface_y_all = anterior_surface(surface_valid);
    min_anterior_y = round(median(valid_surface_y_all));
    fprintf('        [Sensor] No body under footprint; using median surface Y=%d\n', min_anterior_y);
end

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

%% ======================== STEP 5: BUILD SPARSE SENSOR MASK ========================

% The mask is the union of N x N active element footprints. Each element
% center sits at (ex - 0.5) * pitch from the aperture corner. Kerf voxels
% are never marked. element_map is built simultaneously so apply_element_averaging
% can group voxels into elements (0 = inactive/kerf/removed).

sensor_mask = false(grid_dims);

local_nx = sensor_x_range(2) - sensor_x_range(1) + 1;
local_nz = sensor_z_range(2) - sensor_z_range(1) + 1;
element_map = zeros(local_nx, local_nz);

% Physical coordinate of the aperture corner (top-left voxel center, in mm).
aperture_origin_x_mm = origin(1) + (sensor_x_range(1) - 1) * dx;
aperture_origin_z_mm = origin(3) + (sensor_z_range(1) - 1) * dz;

% Half-extent of an active element footprint, in voxel units.
half_nx_active = (element_size_mm / 2) / dx;
half_nz_active = (element_size_mm / 2) / dz;

% Pre-allocate element_positions_mm: [N_total x 3] in patient coords.
N_total_elements = elements_per_side * elements_per_side;
element_positions_mm = zeros(N_total_elements, 3);

% Track which elements actually got at least one voxel (for aliasing report).
element_voxel_counts = zeros(N_total_elements, 1);

for ex = 1:elements_per_side
    % Element center physical X (mm)
    cx_mm = aperture_origin_x_mm + (ex - 0.5) * element_pitch_mm;
    % Element center in local voxel coords (1-based, fractional)
    cx_local = (cx_mm - aperture_origin_x_mm) / dx + 1;

    % Active footprint X voxel range (clamped to local aperture)
    ix_lo = max(1,        ceil(cx_local - half_nx_active));
    ix_hi = min(local_nx, floor(cx_local + half_nx_active));
    % Sub-voxel element: snap to nearest voxel
    if ix_hi < ix_lo
        ix_lo = max(1, min(local_nx, round(cx_local)));
        ix_hi = ix_lo;
    end

    for ez = 1:elements_per_side
        cz_mm = aperture_origin_z_mm + (ez - 0.5) * element_pitch_mm;
        cz_local = (cz_mm - aperture_origin_z_mm) / dz + 1;

        iz_lo = max(1,        ceil(cz_local - half_nz_active));
        iz_hi = min(local_nz, floor(cz_local + half_nz_active));
        if iz_hi < iz_lo
            iz_lo = max(1, min(local_nz, round(cz_local)));
            iz_hi = iz_lo;
        end

        elem_idx = (ex - 1) * elements_per_side + ez;
        element_positions_mm(elem_idx, :) = ...
            [cx_mm, origin(2) + (sensor_y_index - 1) * dy, cz_mm];

        % Convert local (ix, iz) ranges to global grid indices and mark.
        for li = ix_lo:ix_hi
            for lj = iz_lo:iz_hi
                gi = sensor_x_range(1) + li - 1;
                gj = sensor_z_range(1) + lj - 1;
                if gi < 1 || gi > Nx || gj < 1 || gj > Nz
                    continue;
                end
                if element_map(li, lj) == 0
                    sensor_mask(sensor_y_index, gi, gj) = true;
                    element_map(li, lj) = elem_idx;
                    element_voxel_counts(elem_idx) = element_voxel_counts(elem_idx) + 1;
                else
                    % Another element already claimed this voxel (aliasing).
                    % Leave existing assignment; the colliding element gets no voxel.
                end
            end
        end
    end
end

num_sensor_voxels = sum(sensor_mask(:));
effective_element_count = sum(element_voxel_counts > 0);
aliased_element_count   = N_total_elements - effective_element_count;

fprintf('        [Sensor] Sparse mask: %d voxels across %d/%d elements\n', ...
    num_sensor_voxels, effective_element_count, N_total_elements);
if aliased_element_count > 0
    warning('determine_sensor_mask:ElementAliasing', ...
        '%d of %d elements collapsed (grid spacing too coarse for full resolution).', ...
        aliased_element_count, N_total_elements);
end

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
                element_map(ix - sensor_x_range(1) + 1, iz - sensor_z_range(1) + 1) = 0;
                placement_valid = false;
            end
        end
    end
end

% Also sync element_map for any body-overlap voxels removed in Check 1.
for ix = sensor_x_range(1):sensor_x_range(2)
    for iz = sensor_z_range(1):sensor_z_range(2)
        if ~sensor_mask(sensor_y_index, ix, iz)
            element_map(ix - sensor_x_range(1) + 1, iz - sensor_z_range(1) + 1) = 0;
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

%% ======================== STEP 7: ELEMENT DIAGNOSTICS ========================

% element_map and per-element voxel counts were built in STEP 5.
% Recompute counts after validation (voxels may have been removed for body
% overlap, exclusion zone, etc.) to get an accurate effective_element_count.
element_voxel_counts = zeros(N_total_elements, 1);
em_vec = element_map(:);
for v = 1:numel(em_vec)
    e = em_vec(v);
    if e > 0
        element_voxel_counts(e) = element_voxel_counts(e) + 1;
    end
end

effective_element_count = sum(element_voxel_counts > 0);
aliased_element_count   = N_total_elements - effective_element_count;
num_elements            = effective_element_count;

if effective_element_count > 0
    vpe_assigned = element_voxel_counts(element_voxel_counts > 0);
    voxels_per_element_stats = [min(vpe_assigned), median(vpe_assigned), max(vpe_assigned)];
else
    voxels_per_element_stats = [0, 0, 0];
end

% Actual fill factor within the aperture footprint rectangle
aperture_footprint_voxels = sensor_nx * sensor_nz;
if aperture_footprint_voxels > 0
    fill_factor_actual = sum(sensor_mask(:)) / aperture_footprint_voxels;
else
    fill_factor_actual = 0;
end

fprintf('        [Sensor] Element diagnostics: %d/%d effective, %d aliased/removed\n', ...
    effective_element_count, N_total_elements, aliased_element_count);
fprintf('        [Sensor] Voxels per element: min=%d, median=%d, max=%d\n', ...
    voxels_per_element_stats(1), voxels_per_element_stats(2), voxels_per_element_stats(3));
fprintf('        [Sensor] Fill factor (actual): %.1f%% (spec: %.1f%%)\n', ...
    100 * fill_factor_actual, 100 * (element_size_mm / element_pitch_mm)^2);

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
sensor_info.element_positions_mm        = element_positions_mm;
sensor_info.element_pitch_mm            = element_pitch_mm;
sensor_info.element_size_mm             = element_size_mm;
sensor_info.kerf_mm                     = kerf_mm;
sensor_info.elements_per_side           = elements_per_side;
sensor_info.effective_element_count     = effective_element_count;
sensor_info.aliased_element_count       = aliased_element_count;
sensor_info.voxels_per_element_stats    = voxels_per_element_stats;
sensor_info.aperture_mm                 = aperture_mm;
sensor_info.fill_factor_actual          = fill_factor_actual;
sensor_info.surface_to_dose_distance_mm = surface_to_dose_distance_mm;
sensor_info.sensor_size_voxels          = [sensor_nx, sensor_nz];
sensor_info.placement_valid             = placement_valid;
sensor_info.num_sensor_voxels           = num_sensor_voxels_final;
sensor_info.anterior_surface            = anterior_surface;
sensor_info.standoff_voxels             = standoff_voxels;
sensor_info.gantry_angle                = field_dose.gantry_angle;

% Grid expansion applied to fit sensor outside exclusion zone. Caller must
% apply matching water-filled padding to medium fields and p0 so the sensor
% mask's coordinate system matches the simulation grid. Order matches the
% native (Ny, Nx, Nz) array layout returned by size(sensor_mask).
sensor_info.grid_pad = struct( ...
    'y_pre',  grid_pad_y_pre,  'y_post', grid_pad_y_post, ...
    'x_pre',  grid_pad_x_pre,  'x_post', grid_pad_x_post, ...
    'z_pre',  grid_pad_z_pre,  'z_post', grid_pad_z_post, ...
    'expanded', grid_was_expanded);

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
    sensor_info.element_positions_mm        = zeros(0, 3);
    sensor_info.element_pitch_mm            = 0;
    sensor_info.element_size_mm             = 0;
    sensor_info.kerf_mm                     = 0;
    sensor_info.elements_per_side           = 0;
    sensor_info.effective_element_count     = 0;
    sensor_info.aliased_element_count       = 0;
    sensor_info.voxels_per_element_stats    = [0, 0, 0];
    sensor_info.aperture_mm                 = 0;
    sensor_info.fill_factor_actual          = 0;
    sensor_info.surface_to_dose_distance_mm = 0;
    sensor_info.sensor_size_voxels          = [0, 0];
    sensor_info.placement_valid             = false;
    sensor_info.num_sensor_voxels           = 0;
    sensor_info.grid_pad                    = struct( ...
        'y_pre', 0, 'y_post', 0, 'x_pre', 0, 'x_post', 0, ...
        'z_pre', 0, 'z_post', 0, 'expanded', false);
end


function val = get_field(s, name, default)
%GET_FIELD Retrieve struct field with fallback to default
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default;
    end
end
