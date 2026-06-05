function [sensor_mask, sensor_info] = determine_sensor_mask(sct_resampled, field_dose, beam_metadata, config)
%DETERMINE_SENSOR_MASK Place a 2D ultrasound array tilted toward beam isocenter
%
%   [sensor_mask, sensor_info] = determine_sensor_mask(sct_resampled, field_dose, beam_metadata, config)
%
%   PURPOSE:
%   Compute the 3D binary mask for a rigid 2D ultrasound array pressed
%   against the patient's anterior abdomen. The array is modeled as a sparse
%   grid of discrete piezoelectric elements (default 32x32) with a physical
%   pitch, active element size, and kerf (dead) gap between elements.
%
%   Placement is decoupled from tilt: the center is fixed first by a flat
%   sweep, then the rotation toward iso is scaled to the largest theta in
%   [0, 1] that keeps every element outside the body.
%
%   When config.aim_at_iso is false, the flat coronal sensor (single Y index
%   slab) is produced and no rotation is applied.
%
%   ALGORITHM:
%   1. Compute anterior surface height map from body mask (excluding couch).
%   2. Compute dose-based Z exclusion zone (>= 10% of peak per Z slice).
%      Caller should pass the SUMMED PLAN DOSE in field_dose.dose_Gy when
%      computing a session-level sensor.
%   3. Find the flat sensor center: sweep candidate rectangles on the
%      anterior surface outside the exclusion zone; pick the one closest
%      to the dose centroid in XZ. Placement is restricted to be at or
%      inferior to the dose centroid Z (toward the feet) - superior
%      placement is never permitted because the ribs above the abdomen
%      acoustically block the signal. Grid-expand inferior with water
%      padding if no in-grid candidate fits.
%   4. Place the sensor: set its Y to (anterior surface min) - standoff;
%      compute the sensor center in mm.
%   5. Tilt toward iso (only if aim_at_iso): build the full Rodrigues
%      rotation from coronal normal [0, -1, 0] to (c - iso)/|c - iso|,
%      then binary-search the largest theta in [0, 1] for which every
%      rotated element stays outside the body. Grid bounds and the dose
%      exclusion projection do NOT limit the search - PML is outside the
%      grid and an in-air element above the skin does not interact with
%      the dose volume.
%   5b. Expand the grid (water padding) if the chosen tilt pushes any
%      element past the original grid edge in X, Y, or Z, so the rotated
%      sensor fits entirely inside the simulation grid.
%   6. Voxelize: if theta > 0, round each rotated element center to its
%      nearest voxel (single-voxel-per-element). Otherwise voxelize the
%      flat sensor with multi-voxel-per-element element footprints.
%   7. Validate: drop sensor voxels that fell inside the body.
%   8. Element diagnostics: count voxels per element, fill factor, etc.
%
%   INPUTS:
%       sct_resampled - Struct with CT resampled to dose grid:
%           .bodyMask       - 3D logical (true = inside body)
%           .couchMask      - 3D logical (true = couch region)
%           .origin         - [x, y, z] mm (DICOM patient coordinates)
%           .spacing        - [dx, dy, dz] mm
%           .dimensions     - [ny, nx, nz] (rows=Y, cols=X, slices=Z)
%       field_dose - Struct with:
%           .dose_Gy        - 3D dose array (Gy). For session-level sensor
%                             placement pass the SUMMED PLAN DOSE.
%           .gantry_angle   - Gantry angle (degrees) - logged only.
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
%           .element_size_mm      - Active element width in mm (default: 2.43)
%           .sensor_standoff_mm   - Gap between body surface and sensor (default: 5)
%           .jaw_margin_mm        - Extra margin around jaw projection (default: 10)
%           .sensor_placement     - Placement side: 'anterior' (default)
%           .pml_size             - PML thickness in voxels (default: 10)
%           .aim_at_iso           - Enable iso-aimed tilt (default: true)
%           .force_turn_angle     - Forced tilt magnitude in degrees about
%                                   the iso-aim axis. 0 = normal feasibility-
%                                   limited binary search. Nonzero = apply
%                                   exactly this angle regardless of body
%                                   or grid; the grid is expanded to fit.
%                                   (default: 30)
%
%   OUTPUTS:
%       sensor_mask - 3D logical array (same size as dose grid), true at sensor voxels.
%       sensor_info - Struct with diagnostic fields:
%           .sensor_y_index              - Y index of the sensor CENTER voxel.
%                                          NOTE: when tilted, sensor voxels span
%                                          multiple Y indices; this is the centroid.
%           .sensor_x_range              - XZ bounding box of true voxels, X indices.
%           .sensor_z_range              - XZ bounding box of true voxels, Z indices.
%           .sensor_center_mm            - [x, y, z] physical position of sensor center.
%           .exclusion_zone              - 2D logical on anterior surface (X x Z)
%           .element_map                 - 2D (X x Z) projection of element index per
%                                          voxel. Best-effort for back-compat; not
%                                          authoritative for tilted sensors. Use
%                                          voxel_element_idx instead.
%           .voxel_element_idx           - Nvox x 1 element-index lookup keyed by
%                                          find(sensor_mask) linear-index ordering.
%                                          AUTHORITATIVE voxel->element mapping.
%           .num_elements                - Effective number of resolved elements.
%           .element_positions_mm        - [N x 3] physical centers of all elements
%                                          (reflects tilt).
%           .element_pitch_mm            - Echoed pitch (mm).
%           .element_size_mm             - Echoed active element width (mm).
%           .kerf_mm                     - Derived kerf = pitch - size (mm).
%           .elements_per_side           - Echoed N (default 32).
%           .effective_element_count     - Elements with >= 1 assigned voxel.
%           .aliased_element_count       - Elements that lost their voxel.
%           .voxels_per_element_stats    - [min, median, max] voxels per element.
%           .aperture_mm                 - Derived total aperture (mm).
%           .fill_factor_actual          - Actual mask fill within aperture rectangle.
%           .surface_to_dose_distance_mm - Distance from sensor center to dose centroid.
%           .sensor_size_voxels          - [nx_sensor, nz_sensor] aperture footprint in voxels.
%           .placement_valid             - Boolean, true if placement passed all checks.
%           .tilt_theta                  - Tilt scalar in [0, 1].
%           .tilt_angle_deg              - theta * alpha_full in degrees.
%           .tilt_alpha_full_deg         - Full aim-at-iso rotation magnitude in degrees.
%           .rotation_R                  - 3x3 rotation applied to the sensor.
%           .rotation_axis               - Unit rotation axis (3-vec).
%           .aim_target_mm               - Isocenter target, [x, y, z] mm.
%           .aim_normal                  - Unit normal vector used for aim.
%           .aim_enabled                 - True if aim-at-iso was applied.
%           .grid_pad                    - Grid-expansion fallback bookkeeping.
%
%   COORDINATE SYSTEM:
%       Dose grid uses DICOM patient coordinates (HFS convention):
%           X = patient left-right (columns, dim 2)
%           Y = patient anterior-posterior (rows, dim 1).
%               Lower index = more anterior.
%           Z = patient superior-inferior (slices, dim 3).
%               Higher index = SUPERIOR (cranial/head).
%               Lower  index = INFERIOR (caudal/feet).
%       Array indexing: array(Y, X, Z)
%
%   NOTES:
%       - Source-to-axis distance (SAD) for Halcyon/ETHOS is 100 cm.
%       - All beams in a plan share the same isocenter (DICOM guarantee);
%         the iso for tilt is read from beam_metadata(1).isocenter and
%         consistency is verified across remaining beams.
%       - Uses warnings (not errors) for non-fatal issues to support batch processing.
%
%   DEPENDENCIES:
%       - Image Processing Toolbox (bwconncomp, regionprops)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: May 2026
%   VERSION: 2.0 (tilted)
%
%   See also: run_single_field_simulation, step15_process_doses, apply_element_averaging

%% ======================== CONFIG DEFAULTS ========================

elements_per_side = get_field(config, 'elements_per_side', 32);
element_pitch_mm  = get_field(config, 'element_pitch_mm', 7.00);
element_size_mm   = get_field(config, 'element_size_mm', 2.43);
standoff_mm       = get_field(config, 'sensor_standoff_mm', 5);
jaw_margin_mm     = get_field(config, 'jaw_margin_mm', 10); %#ok<NASGU>
placement_side    = get_field(config, 'sensor_placement', 'anterior');
pml_size          = get_field(config, 'pml_size', 10);
pml_size = 1; % Hardcoded because script logic assumes pml inside.

aim_at_iso        = get_field(config, 'aim_at_iso', true);
force_turn_angle_deg = get_field(config, 'force_turn_angle', 300);

% Kerf is derived; never accepted from config.
kerf_mm = element_pitch_mm - element_size_mm;
if kerf_mm < 0
    error('determine_sensor_mask:InvalidGeometry', ...
        'element_size_mm (%.3f) cannot exceed element_pitch_mm (%.3f).', ...
        element_size_mm, element_pitch_mm);
end

% Total aperture footprint (mm) derived from N x pitch.
aperture_mm = elements_per_side * element_pitch_mm;

SAD_mm = 1000;  %#ok<NASGU>  % Halcyon/ETHOS source-to-axis distance: 100 cm

fprintf('        [Sensor] Placing %s array: %dx%d elements, pitch %.2f mm, size %.2f mm, kerf %.2f mm\n', ...
    placement_side, elements_per_side, elements_per_side, ...
    element_pitch_mm, element_size_mm, kerf_mm);
fprintf('        [Sensor] Aperture %.1f mm, standoff %.0f mm, fill factor (spec) %.0f%%\n', ...
    aperture_mm, standoff_mm, ...
    100 * (element_size_mm / element_pitch_mm)^2);
if aim_at_iso
    fprintf('        [Sensor] aim_at_iso=true (will tilt placed sensor toward iso)\n');
else
    fprintf('        [Sensor] aim_at_iso=false (flat coronal placement)\n');
end

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
anterior_surface = NaN(Nx, Nz);
for ix = 1:Nx
    for iz = 1:Nz
        col = squeeze(body(:, ix, iz));
        y_indices = find(col);
        if ~isempty(y_indices)
            anterior_surface(ix, iz) = min(y_indices);
        end
    end
end

surface_valid = ~isnan(anterior_surface);
num_surface_pts = sum(surface_valid(:));
fprintf('        [Sensor] Anterior surface: %d valid columns\n', num_surface_pts);

if num_surface_pts == 0
    warning('determine_sensor_mask:NoSurface', ...
        'No anterior body surface found. Returning empty sensor mask.');
    [sensor_mask, sensor_info] = empty_result(grid_dims);
    return;
end

%% ======================== STEP 2: DOSE-BASED EXCLUSION ZONE ========================

exclusion_zone = false(Nx, Nz);
dose = field_dose.dose_Gy;
dose_max = max(dose(:));

if ~isempty(dose) && dose_max > 0
    dose_thresh = 0.10 * dose_max;
    if ~isequal(size(dose), [Ny, Nx, Nz])
        warning('determine_sensor_mask:DoseSizeMismatch', ...
            'field_dose.dose_Gy size %s does not match sct grid [%d %d %d]. Skipping exclusion.', ...
            mat2str(size(dose)), Ny, Nx, Nz);
    else
        slice_max = squeeze(max(max(dose, [], 1), [], 2));
        excluded_z_mask = slice_max(:) >= dose_thresh;
        exclusion_zone(:, excluded_z_mask) = true;

        excluded_z_idx = find(excluded_z_mask);
        if isempty(excluded_z_idx)
            z_lo = 0; z_hi = 0;
        else
            z_lo = excluded_z_idx(1);
            z_hi = excluded_z_idx(end);
        end
        fprintf('        [Sensor] Dose-based exclusion: %d/%d Z-slices (Z=[%d,%d]) >= 10%% of peak dose (%.3f Gy)\n', ...
            sum(excluded_z_mask), Nz, z_lo, z_hi, dose_thresh);
    end
else
    warning('determine_sensor_mask:NoDose', ...
        'Dose array is empty or zero; exclusion zone left empty.');
end

fprintf('        [Sensor] Exclusion zone: %d voxels (%.1f%% of surface)\n', ...
    sum(exclusion_zone(:) & surface_valid(:)), ...
    100 * sum(exclusion_zone(:) & surface_valid(:)) / max(1, num_surface_pts));

%% ======================== RESOLVE AIM TARGET ========================

iso_mm = [];
if aim_at_iso
    if isempty(beam_metadata)
        warning('determine_sensor_mask:NoBeamMetadata', ...
            'aim_at_iso=true but beam_metadata is empty. Disabling tilt.');
        aim_at_iso = false;
    else
        iso_mm = resolve_iso(beam_metadata);
        if isempty(iso_mm)
            warning('determine_sensor_mask:NoIsocenter', ...
                'beam_metadata does not contain a usable isocenter. Disabling tilt.');
            aim_at_iso = false;
        else
            fprintf('        [Sensor] Aim target (isocenter): [%.1f, %.1f, %.1f] mm\n', ...
                iso_mm(1), iso_mm(2), iso_mm(3));
        end
    end
end

%% ======================== STEP 3: FIND FLAT SENSOR CENTER ========================

% Available region: on body surface AND not in exclusion zone
available = surface_valid & ~exclusion_zone;

% PML lives outside the grid: no edge buffer is required. If the placed (or
% later tilted) sensor extends past the grid, the grid is expanded to fit.
pml_margin_x = 0;
pml_margin_z = 0;

if sum(available(:)) < sensor_nx * sensor_nz
    warning('determine_sensor_mask:InsufficientSpace', ...
        'Available anterior surface (%d voxels) may be too small for sensor (%d voxels).', ...
        sum(available(:)), sensor_nx * sensor_nz);
end

% Compute dose centroid in X-Z (proximity targeting for the flat sweep)
dose_centroid_mm = compute_dose_centroid_mm(field_dose, origin, dx, dy, dz);
dose_centroid_ix = round((dose_centroid_mm(1) - origin(1)) / dx) + 1;
dose_centroid_iz = round((dose_centroid_mm(3) - origin(3)) / dz) + 1;
fprintf('        [Sensor] Dose centroid (voxel): X=%d, Z=%d\n', dose_centroid_ix, dose_centroid_iz);

% Standoff in Y voxels (places sensor anterior to body surface)
standoff_voxels = ceil(standoff_mm / dy);

% Track grid-expansion padding (set by fallback path)
grid_pad_x_pre  = 0; grid_pad_x_post = 0;
grid_pad_y_pre  = 0; grid_pad_y_post = 0;
grid_pad_z_pre  = 0; grid_pad_z_post = 0;
grid_was_expanded = false;

% Output allocation
sensor_mask = false(grid_dims);
sensor_element_label = zeros(grid_dims, 'uint32');  % 0 = no sensor; else element idx
voxel_element_idx = [];
element_positions_mm = [];
element_map = [];
tilt_R = eye(3);
tilt_theta = 0;
tilt_alpha = 0;
tilt_axis = [0; 0; 1];
aim_normal_used = [0, -1, 0];
N_total_elements = elements_per_side * elements_per_side;

% Sweep candidate flat-sensor rectangles, restricted to placements that sit
% at or inferior to (toward the FEET) the dose centroid Z. Superior slides
% (toward the head) would put the sensor under the ribs, blocking the
% acoustic signal, and are never considered.
%
% HFS convention (see CLAUDE.md "Coordinate System"): higher Z voxel index
% = superior (cranial/head); lower Z voxel index = inferior (caudal/feet).
best_dist = Inf;
best_ix_start = [];
best_iz_start = [];
half_nx = floor(sensor_nx / 2);
half_nz = floor(sensor_nz / 2);

% Only consider placements whose center is inferior to (or aligned with) the
% dose centroid Z.
for ix_start = 1:(Nx - sensor_nx + 1)
    for iz_start = 1:(Nz - sensor_nz + 1)
        cz = iz_start + half_nz;
        if cz > dose_centroid_iz
            continue;  % superior placement - never allowed
        end
        ix_end = ix_start + sensor_nx - 1;
        iz_end = iz_start + sensor_nz - 1;
        patch = available(ix_start:ix_end, iz_start:iz_end);
        if all(patch(:))
            cx = ix_start + half_nx;
            dist = sqrt((cx - dose_centroid_ix)^2 + (cz - dose_centroid_iz)^2);
            if dist < best_dist
                best_dist = dist;
                best_ix_start = ix_start;
                best_iz_start = iz_start;
            end
        end
    end
end

% Superior placement is NEVER permitted: a sensor that drifts to higher iz
% (toward the head) sits under the ribs, which acoustically block the signal.
% There is no superior fallback. If no inferior placement fits in the
% original grid, control falls through to the grid-expansion fallback below,
% which always pads inferior (lower iz).
% Grid-expansion fallback (water padding) if no inferior in-grid candidate.
% Always expands inferior (toward the feet); superior is never an option.
if isempty(best_ix_start)
    fprintf('        [Sensor] No inferior placement in grid; expanding grid inferior past exclusion zone.\n');

    ix_start_target = dose_centroid_ix - half_nx;
    ix_start_min_orig = pml_margin_x + 1;
    ix_start_max_orig = Nx - pml_margin_x - sensor_nx + 1;

    if ix_start_max_orig >= ix_start_min_orig
        ix_start_orig = max(ix_start_min_orig, min(ix_start_max_orig, ix_start_target));
    else
        extra_x = (sensor_nx + 2 * pml_margin_x) - Nx;
        grid_pad_x_pre  = floor(extra_x / 2);
        grid_pad_x_post = extra_x - grid_pad_x_pre;
        ix_start_orig = (pml_margin_x + 1) - grid_pad_x_pre;
    end
    ix_end_orig = ix_start_orig + sensor_nx - 1;

    ix_lo_clip = max(1, ix_start_orig);
    ix_hi_clip = min(Nx, ix_end_orig);
    if ix_lo_clip > ix_hi_clip
        excl_z_indices = [];
    else
        excl_strip = exclusion_zone(ix_lo_clip:ix_hi_clip, :);
        excl_z_indices = find(any(excl_strip, 1));
    end

    if isempty(excl_z_indices)
        iz_start_orig   = max(pml_margin_z + 1, ...
                          min(Nz - pml_margin_z - sensor_nz + 1, ...
                              dose_centroid_iz - half_nz));
        grid_pad_z_pre  = 0;
        grid_pad_z_post = 0;
        fprintf('        [Sensor] No exclusion in X strip; centered Z=%d (no Z pad)\n', iz_start_orig);
    else
        % HFS convention: higher iz = superior (cranial); lower iz = inferior
        % (caudal). The sensor must ALWAYS sit INFERIOR to the exclusion (lower
        % iz, toward feet) so the ribs above the abdomen do not block the
        % signal. Superior placement is never permitted; if the inferior
        % candidate runs off the lower grid edge we pad Z- (water) to fit it.
        excl_z_min = min(excl_z_indices);

        % Inferior candidate: butt up just BEFORE excl_z_min (lower Z = feet).
        iz_end_inf    = excl_z_min - 1;
        iz_start_inf  = iz_end_inf - sensor_nz + 1;
        pad_z_pre_inf = max(0, (pml_margin_z + 1) - iz_start_inf);

        iz_start_orig   = iz_start_inf;
        grid_pad_z_pre  = pad_z_pre_inf;
        fprintf('        [Sensor] Inferior to exclusion (Z=%d orig); padding Z- by %d voxels (water)\n', ...
            iz_start_orig, grid_pad_z_pre);
    end

    new_Nx = Nx + grid_pad_x_pre + grid_pad_x_post;
    new_Nz = Nz + grid_pad_z_pre + grid_pad_z_post;

    body_expanded = false(Ny, new_Nx, new_Nz);
    body_expanded(:, grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = body;
    body = body_expanded;

    surface_valid_expanded = false(new_Nx, new_Nz);
    surface_valid_expanded(grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = surface_valid;
    surface_valid = surface_valid_expanded; %#ok<NASGU>

    exclusion_zone_expanded = false(new_Nx, new_Nz);
    exclusion_zone_expanded(grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = exclusion_zone;
    exclusion_zone = exclusion_zone_expanded;

    anterior_surface_expanded = NaN(new_Nx, new_Nz);
    anterior_surface_expanded(grid_pad_x_pre + (1:Nx), grid_pad_z_pre + (1:Nz)) = anterior_surface;
    anterior_surface = anterior_surface_expanded;

    Nx_pre_expansion = Nx;
    Nz_pre_expansion = Nz;
    Nx = new_Nx;
    Nz = new_Nz;
    grid_dims = [Ny, Nx, Nz];

    sensor_mask = false(grid_dims);
    sensor_element_label = zeros(grid_dims, 'uint32');

    best_ix_start = ix_start_orig + grid_pad_x_pre;
    best_iz_start = iz_start_orig + grid_pad_z_pre;
    best_dist     = sqrt( (best_ix_start + half_nx - (dose_centroid_ix + grid_pad_x_pre))^2 + ...
                          (best_iz_start + half_nz - (dose_centroid_iz + grid_pad_z_pre))^2 );

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

sensor_x_range = [best_ix_start, best_ix_start + sensor_nx - 1];
sensor_z_range = [best_iz_start, best_iz_start + sensor_nz - 1];

fprintf('        [Sensor] Placement: X=[%d,%d], Z=[%d,%d] (dist to dose: %.1f voxels)\n', ...
    sensor_x_range(1), sensor_x_range(2), sensor_z_range(1), sensor_z_range(2), best_dist);

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
    N_total_elements = elements_per_side * elements_per_side;
end

%% ======================== STEP 4: SENSOR Y INDEX & CENTER ========================

surface_patch = anterior_surface(sensor_x_range(1):sensor_x_range(2), ...
                                  sensor_z_range(1):sensor_z_range(2));
min_anterior_y = min(surface_patch(:), [], 'omitnan');
if isnan(min_anterior_y)
    valid_surface_y_all = anterior_surface(~isnan(anterior_surface));
    min_anterior_y = round(median(valid_surface_y_all));
    fprintf('        [Sensor] No body under footprint; using median surface Y=%d\n', min_anterior_y);
end
sensor_y_index = max(1, min_anterior_y - standoff_voxels);
% PML is outside the grid; sensor may sit at vy = 1. If a tilted sensor
% later needs more anterior room, the grid is expanded below.

sensor_center_x = origin(1) + (mean(sensor_x_range) - 1) * dx;
sensor_center_y = origin(2) + (sensor_y_index - 1) * dy;
sensor_center_z = origin(3) + (mean(sensor_z_range) - 1) * dz;

fprintf('        [Sensor] Y index: %d (surface min Y: %d, standoff: %d voxels)\n', ...
    sensor_y_index, min_anterior_y, standoff_voxels);
fprintf('        [Sensor] Center (mm): [%.1f, %.1f, %.1f]\n', ...
    sensor_center_x, sensor_center_y, sensor_center_z);

%% ======================== STEP 4.5: TILT TOWARD ISOCENTER ========================
% Center is fixed by the flat sweep above. Compute the rotation that would
% point the sensor normal directly at iso, then binary-search the largest
% scaling theta in [0, 1] that keeps every element outside the body and
% outside the exclusion projection. Result populates tilt_R / tilt_theta.

if aim_at_iso
    c_mm = [sensor_center_x, sensor_center_y, sensor_center_z];
    [~, tilt_alpha, tilt_axis_candidate] = build_aim_rotation(c_mm, iso_mm);

    if force_turn_angle_deg ~= 0
        % --- FORCED TILT ---
        % Apply exactly the requested angle about the iso-aim axis; ignore
        % body collision and grid bounds. The grid is expanded later (Step
        % 4.6) to fit the rotated sensor.
        force_rad = force_turn_angle_deg * pi / 180;
        tilt_axis  = tilt_axis_candidate;
        tilt_alpha = force_rad;   % store applied angle in alpha
        tilt_theta = 1;           % theta=1 => tilt_angle_deg == force_turn_angle
        tilt_R     = rodrigues(tilt_axis, force_rad);
        d = c_mm - iso_mm;
        nd = norm(d);
        if nd > eps
            aim_normal_used = d / nd;
        end
        fprintf('        [Sensor] FORCED tilt: %.2f deg about iso-aim axis (config.force_turn_angle); feasibility checks bypassed\n', ...
            force_turn_angle_deg);
    elseif abs(tilt_alpha) < 1e-6
        fprintf('        [Sensor] Sensor center already aligned with iso direction; no tilt needed.\n');
    else
        theta_max = binary_search_theta_max(c_mm, tilt_axis_candidate, tilt_alpha, ...
            elements_per_side, element_pitch_mm, origin, dx, dy, dz, ...
            body, exclusion_zone, pml_size, Nx, Ny, Nz);
        if isnan(theta_max)
            theta_max = 0;
        end
        if theta_max > 0
            tilt_theta = theta_max;
            tilt_axis  = tilt_axis_candidate;
            tilt_R     = rodrigues(tilt_axis, tilt_theta * tilt_alpha);
            d = c_mm - iso_mm;
            nd = norm(d);
            if nd > eps
                aim_normal_used = d / nd;
            end
            fprintf('        [Sensor] Tilt: theta=%.3f (%.2f deg of full aim %.2f deg)\n', ...
                tilt_theta, tilt_theta * tilt_alpha * 180/pi, tilt_alpha * 180/pi);
        else
            fprintf('        [Sensor] Tilt search: max feasible theta=0; keeping flat placement.\n');
        end
    end
end

%% ======================== STEP 4.6: GRID EXPANSION TO FIT TILT ========================
% Rotation feasibility is body-only (PML is outside the grid). After the
% tilt is chosen, some elements may have rotated past the original grid
% edge in any of X, Y, or Z. Expand the grid (water padding) so every
% rotated element center voxelizes inside the grid.

if tilt_theta > 0
    c_mm_for_tilt = [sensor_center_x, sensor_center_y, sensor_center_z];
    half_span_pad = (elements_per_side - 1) / 2;
    elem_vx_all = zeros(N_total_elements, 1);
    elem_vy_all = zeros(N_total_elements, 1);
    elem_vz_all = zeros(N_total_elements, 1);
    ee = 0;
    for ex = 1:elements_per_side
        u_off = (ex - 1 - half_span_pad) * element_pitch_mm;
        for ez = 1:elements_per_side
            v_off = (ez - 1 - half_span_pad) * element_pitch_mm;
            ee = ee + 1;
            p_local = [u_off; v_off; 0];
            p_w = (tilt_R * p_local)' + c_mm_for_tilt;
            elem_vx_all(ee) = round((p_w(1) - origin(1)) / dx) + 1;
            elem_vy_all(ee) = round((p_w(2) - origin(2)) / dy) + 1;
            elem_vz_all(ee) = round((p_w(3) - origin(3)) / dz) + 1;
        end
    end

    pad_x_pre_tilt  = max(0, 1 - min(elem_vx_all));
    pad_x_post_tilt = max(0, max(elem_vx_all) - Nx);
    pad_y_pre_tilt  = max(0, 1 - min(elem_vy_all));
    pad_y_post_tilt = max(0, max(elem_vy_all) - Ny);
    pad_z_pre_tilt  = max(0, 1 - min(elem_vz_all));
    pad_z_post_tilt = max(0, max(elem_vz_all) - Nz);

    tilt_pad_total = pad_x_pre_tilt + pad_x_post_tilt + ...
                     pad_y_pre_tilt + pad_y_post_tilt + ...
                     pad_z_pre_tilt + pad_z_post_tilt;

    if tilt_pad_total > 0
        fprintf('        [Sensor] Tilt pushes elements past grid; expanding by [X -%d/+%d, Y -%d/+%d, Z -%d/+%d] voxels (water)\n', ...
            pad_x_pre_tilt, pad_x_post_tilt, ...
            pad_y_pre_tilt, pad_y_post_tilt, ...
            pad_z_pre_tilt, pad_z_post_tilt);

        new_Nx = Nx + pad_x_pre_tilt + pad_x_post_tilt;
        new_Ny = Ny + pad_y_pre_tilt + pad_y_post_tilt;
        new_Nz = Nz + pad_z_pre_tilt + pad_z_post_tilt;

        body_expanded = false(new_Ny, new_Nx, new_Nz);
        body_expanded(pad_y_pre_tilt + (1:Ny), ...
                      pad_x_pre_tilt + (1:Nx), ...
                      pad_z_pre_tilt + (1:Nz)) = body;
        body = body_expanded;

        exclusion_zone_expanded = false(new_Nx, new_Nz);
        exclusion_zone_expanded(pad_x_pre_tilt + (1:Nx), ...
                                pad_z_pre_tilt + (1:Nz)) = exclusion_zone;
        exclusion_zone = exclusion_zone_expanded;

        anterior_surface_expanded = NaN(new_Nx, new_Nz);
        anterior_surface_expanded(pad_x_pre_tilt + (1:Nx), ...
                                  pad_z_pre_tilt + (1:Nz)) = anterior_surface;
        anterior_surface = anterior_surface_expanded;

        % Shift voxel-space bookkeeping; mm coordinates are unchanged after
        % the origin shift below.
        sensor_x_range = sensor_x_range + pad_x_pre_tilt;
        sensor_z_range = sensor_z_range + pad_z_pre_tilt;
        sensor_y_index = sensor_y_index + pad_y_pre_tilt;
        dose_centroid_ix = dose_centroid_ix + pad_x_pre_tilt;
        dose_centroid_iz = dose_centroid_iz + pad_z_pre_tilt;

        % Origin = mm position of voxel (1,1,1); pre-padding moves it back.
        origin(1) = origin(1) - pad_x_pre_tilt * dx;
        origin(2) = origin(2) - pad_y_pre_tilt * dy;
        origin(3) = origin(3) - pad_z_pre_tilt * dz;

        Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
        grid_dims = [Ny, Nx, Nz];
        sensor_mask = false(grid_dims);
        sensor_element_label = zeros(grid_dims, 'uint32');

        grid_pad_x_pre  = grid_pad_x_pre  + pad_x_pre_tilt;
        grid_pad_x_post = grid_pad_x_post + pad_x_post_tilt;
        grid_pad_y_pre  = grid_pad_y_pre  + pad_y_pre_tilt;
        grid_pad_y_post = grid_pad_y_post + pad_y_post_tilt;
        grid_pad_z_pre  = grid_pad_z_pre  + pad_z_pre_tilt;
        grid_pad_z_post = grid_pad_z_post + pad_z_post_tilt;
        grid_was_expanded = true;
    end
end

%% ======================== STEP 5: BUILD SPARSE SENSOR MASK ========================

if tilt_theta > 0
    % --- TILTED voxelization: round each rotated element center to a voxel ---
    element_positions_mm = zeros(N_total_elements, 3);
    half_span = (elements_per_side - 1) / 2;
    c_mm = [sensor_center_x, sensor_center_y, sensor_center_z];

    e = 0;
    for ex = 1:elements_per_side
        u_off = (ex - 1 - half_span) * element_pitch_mm;
        for ez = 1:elements_per_side
            v_off = (ez - 1 - half_span) * element_pitch_mm;
            e = e + 1;
            p_local = [u_off; v_off; 0];
            p_w = (tilt_R * p_local)' + c_mm;
            element_positions_mm(e, :) = p_w;

            vx = round((p_w(1) - origin(1)) / dx) + 1;
            vy = round((p_w(2) - origin(2)) / dy) + 1;
            vz = round((p_w(3) - origin(3)) / dz) + 1;
            if vx < 1 || vx > Nx || vy < 1 || vy > Ny || vz < 1 || vz > Nz
                continue;
            end
            if sensor_element_label(vy, vx, vz) == 0
                sensor_mask(vy, vx, vz) = true;
                sensor_element_label(vy, vx, vz) = uint32(e);
            end
        end
    end

    % Update XZ bounding box / sensor_y centroid from the actual voxelization.
    sensor_lin_tmp = find(sensor_mask);
    if ~isempty(sensor_lin_tmp)
        [sv_y_tmp, sv_x_tmp, sv_z_tmp] = ind2sub(size(sensor_mask), sensor_lin_tmp);
        sensor_x_range = [min(sv_x_tmp), max(sv_x_tmp)];
        sensor_z_range = [min(sv_z_tmp), max(sv_z_tmp)];
        sensor_y_index = round(mean(sv_y_tmp));
    end

    % Best-effort 2D element_map projection (for back-compat plot tools)
    local_nx = sensor_x_range(2) - sensor_x_range(1) + 1;
    local_nz = sensor_z_range(2) - sensor_z_range(1) + 1;
    element_map = zeros(local_nx, local_nz);
    for ee = 1:N_total_elements
        p_w = element_positions_mm(ee, :);
        vx = round((p_w(1) - origin(1)) / dx) + 1;
        vz = round((p_w(3) - origin(3)) / dz) + 1;
        li = vx - sensor_x_range(1) + 1;
        lj = vz - sensor_z_range(1) + 1;
        if li >= 1 && li <= local_nx && lj >= 1 && lj <= local_nz
            if element_map(li, lj) == 0
                element_map(li, lj) = ee;
            end
        end
    end
else
    % --- FLAT voxelization (legacy multi-voxel-per-element) ---
    local_nx = sensor_x_range(2) - sensor_x_range(1) + 1;
    local_nz = sensor_z_range(2) - sensor_z_range(1) + 1;
    element_map = zeros(local_nx, local_nz);

    aperture_origin_x_mm = origin(1) + (sensor_x_range(1) - 1) * dx;
    aperture_origin_z_mm = origin(3) + (sensor_z_range(1) - 1) * dz;

    half_nx_active = (element_size_mm / 2) / dx;
    half_nz_active = (element_size_mm / 2) / dz;

    element_positions_mm = zeros(N_total_elements, 3);

    for ex = 1:elements_per_side
        cx_mm = aperture_origin_x_mm + (ex - 0.5) * element_pitch_mm;
        cx_local = (cx_mm - aperture_origin_x_mm) / dx + 1;
        ix_lo = max(1,        ceil(cx_local - half_nx_active));
        ix_hi = min(local_nx, floor(cx_local + half_nx_active));
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

            for li = ix_lo:ix_hi
                for lj = iz_lo:iz_hi
                    gi = sensor_x_range(1) + li - 1;
                    gj = sensor_z_range(1) + lj - 1;
                    if gi < 1 || gi > Nx || gj < 1 || gj > Nz
                        continue;
                    end
                    if element_map(li, lj) == 0
                        sensor_mask(sensor_y_index, gi, gj) = true;
                        sensor_element_label(sensor_y_index, gi, gj) = uint32(elem_idx);
                        element_map(li, lj) = elem_idx;
                    end
                end
            end
        end
    end
end

%% ======================== STEP 6: VALIDATION ========================

placement_valid = true;

% The only post-placement constraint is "no sensor voxel inside the body."
% PML is outside the grid and the dose-exclusion projection is intentionally
% allowed under the tilted sensor (elements in air above the skin do not
% interact with the dose).
overlap_body = sensor_mask & body;
if any(overlap_body(:))
    num_overlap = sum(overlap_body(:));
    warning('determine_sensor_mask:BodyOverlap', ...
        'Sensor overlaps body mask at %d voxels. Removing overlapping voxels.', num_overlap);
    sensor_mask = sensor_mask & ~body;
    sensor_element_label(overlap_body) = 0;
    placement_valid = false;
end

num_sensor_voxels_final = sum(sensor_mask(:));
fprintf('        [Sensor] Final sensor: %d voxels (valid: %s)\n', ...
    num_sensor_voxels_final, mat2str(placement_valid));

%% ======================== STEP 7: ELEMENT DIAGNOSTICS ========================

% Build authoritative voxel_element_idx from sensor_element_label.
sensor_lin_final = find(sensor_mask);
voxel_element_idx = double(sensor_element_label(sensor_lin_final));

% Count voxels per element using the authoritative label.
element_voxel_counts = zeros(N_total_elements, 1);
for v = 1:numel(voxel_element_idx)
    e = voxel_element_idx(v);
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

aperture_footprint_voxels = sensor_nx * sensor_nz;
if aperture_footprint_voxels > 0
    fill_factor_actual = num_sensor_voxels_final / aperture_footprint_voxels;
else
    fill_factor_actual = 0;
end

fprintf('        [Sensor] Element diagnostics: %d/%d effective, %d aliased/removed\n', ...
    effective_element_count, N_total_elements, aliased_element_count);
if ~isempty(vpe_assigned)
    fprintf('        [Sensor] Voxels per element: min=%d, median=%d, max=%d\n', ...
        voxels_per_element_stats(1), voxels_per_element_stats(2), voxels_per_element_stats(3));
end
fprintf('        [Sensor] Fill factor (actual): %.1f%% (spec: %.1f%%)\n', ...
    100 * fill_factor_actual, 100 * (element_size_mm / element_pitch_mm)^2);

if aliased_element_count > 0
    warning('determine_sensor_mask:ElementAliasing', ...
        '%d of %d elements collapsed (grid spacing too coarse or voxels removed by validation).', ...
        aliased_element_count, N_total_elements);
end

%% ======================== COMPUTE DIAGNOSTICS ========================

surface_to_dose_distance_mm = abs(dose_centroid_mm(2) - sensor_center_y);

%% ======================== PACK OUTPUT ========================

sensor_info = struct();
sensor_info.sensor_y_index              = sensor_y_index;
sensor_info.sensor_x_range              = sensor_x_range;
sensor_info.sensor_z_range              = sensor_z_range;
sensor_info.sensor_center_mm            = [sensor_center_x, sensor_center_y, sensor_center_z];
sensor_info.exclusion_zone              = exclusion_zone;
sensor_info.element_map                 = element_map;
sensor_info.voxel_element_idx           = voxel_element_idx;
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

% Tilt geometry
sensor_info.tilt_theta                  = tilt_theta;
sensor_info.tilt_angle_deg              = tilt_theta * tilt_alpha * 180 / pi;
sensor_info.tilt_alpha_full_deg         = tilt_alpha * 180 / pi;
sensor_info.rotation_R                  = tilt_R;
sensor_info.rotation_axis               = tilt_axis(:)';
sensor_info.aim_target_mm               = iso_mm;
sensor_info.aim_normal                  = aim_normal_used;
sensor_info.aim_enabled                 = aim_at_iso && (tilt_theta > 0);

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
    dose = field_dose.dose_Gy;
    dose(dose < 0.01 * max(dose(:))) = 0;
    total_dose = sum(dose(:));
    if total_dose == 0
        dims = size(dose);
        centroid_mm = origin + ([dims(2), dims(1), dims(3)] - 1) / 2 .* [dx, dy, dz];
        return;
    end
    [Ny, Nx, Nz] = size(dose);
    x_vec = origin(1) + (0:Nx-1) * dx;
    y_vec = origin(2) + (0:Ny-1) * dy;
    z_vec = origin(3) + (0:Nz-1) * dz;
    dose_x = squeeze(sum(sum(dose, 1), 3));
    dose_y = squeeze(sum(sum(dose, 2), 3));
    dose_z = squeeze(sum(sum(dose, 1), 2));
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
    sensor_info.voxel_element_idx           = [];
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
    sensor_info.tilt_theta                  = 0;
    sensor_info.tilt_angle_deg              = 0;
    sensor_info.tilt_alpha_full_deg         = 0;
    sensor_info.rotation_R                  = eye(3);
    sensor_info.rotation_axis               = [0, 0, 1];
    sensor_info.aim_target_mm               = [];
    sensor_info.aim_normal                  = [0, -1, 0];
    sensor_info.aim_enabled                 = false;
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


function iso_mm = resolve_iso(beam_metadata)
%RESOLVE_ISO Extract the plan isocenter from beam_metadata; warn if inconsistent.
%   Returns [x, y, z] in mm, or [] if no usable isocenter is present.

    iso_mm = [];
    isos = [];
    for i = 1:numel(beam_metadata)
        if isfield(beam_metadata(i), 'isocenter') && ...
                ~isempty(beam_metadata(i).isocenter) && ...
                numel(beam_metadata(i).isocenter) == 3
            isos(end+1, :) = beam_metadata(i).isocenter(:)'; %#ok<AGROW>
        end
    end
    if isempty(isos)
        return;
    end
    iso_mm = isos(1, :);
    if size(isos, 1) > 1
        max_drift = max(vecnorm(isos - iso_mm, 2, 2));
        if max_drift > 1.0  % mm
            warning('determine_sensor_mask:IsocenterDrift', ...
                'beam_metadata isocenters disagree by up to %.2f mm; using beam 1.', max_drift);
        end
    end
end


function R = rodrigues(k, ang)
%RODRIGUES Rotation matrix from unit axis k (3-vec) and angle ang (rad).
    k = k(:);
    nk = norm(k);
    if nk < eps || abs(ang) < eps
        R = eye(3);
        return;
    end
    k = k / nk;
    K = [   0   -k(3)  k(2);
          k(3)    0   -k(1);
         -k(2)  k(1)    0  ];
    R = eye(3) + sin(ang) * K + (1 - cos(ang)) * (K * K);
end


function [R_full, alpha, axis_k] = build_aim_rotation(c_mm, iso_mm)
%BUILD_AIM_ROTATION Rotation that rotates coronal normal [0,-1,0] to (c-iso)/|c-iso|.
%   c_mm    - 1x3 sensor center (mm)
%   iso_mm  - 1x3 isocenter (mm)
%   Returns R_full (3x3), full rotation angle alpha (rad), unit rotation axis axis_k.

    n0 = [0, -1, 0];
    d  = c_mm(:)' - iso_mm(:)';
    nd = norm(d);
    if nd < eps
        R_full = eye(3);
        alpha = 0;
        axis_k = [0; 0; 1];
        return;
    end
    n_aim = d / nd;

    cross_v = cross(n0, n_aim);
    sin_a = norm(cross_v);
    cos_a = dot(n0, n_aim);
    alpha = atan2(sin_a, cos_a);

    if sin_a < 1e-9
        % Either aligned (alpha~0) or anti-parallel (alpha~pi).
        if cos_a > 0
            R_full = eye(3);
            axis_k = [0; 0; 1];
            alpha = 0;
        else
            % 180 deg flip about any in-plane axis; choose X.
            axis_k = [1; 0; 0];
            R_full = rodrigues(axis_k, alpha);
        end
        return;
    end
    axis_k = cross_v(:) / sin_a;
    R_full = rodrigues(axis_k, alpha);
end


function [feasible, P_world, vox_idx] = check_placement(c_mm, R, ...
        elements_per_side, element_pitch_mm, origin, dx, dy, dz, ...
        body, ~, ~, Nx, Ny, Nz)
%CHECK_PLACEMENT Evaluate whether ALL elements at (c, R) stay outside the body.
%   PML lives outside the grid, so PML is NOT checked. Grid bounds are NOT
%   enforced either: any element that rotates past the grid edge is assumed
%   to be in air/water, and the caller will expand the grid to fit it. The
%   dose-exclusion XZ projection is also ignored - an element in air above
%   the skin does not interact with the dose volume even when its column
%   sits in an irradiated Z slab.
%   Returns:
%     feasible - logical (true iff no rotated element falls inside body)
%     P_world  - [N x 3] world positions (mm) for each element
%     vox_idx  - [N x 3] voxel indices (y, x, z) for each element

    N2 = elements_per_side * elements_per_side;
    P_world = zeros(N2, 3);
    vox_idx = nan(N2, 3);

    half_span = (elements_per_side - 1) / 2;

    feasible = true;
    e = 0;
    for ex = 1:elements_per_side
        u_off = (ex - 1 - half_span) * element_pitch_mm;
        for ez = 1:elements_per_side
            v_off = (ez - 1 - half_span) * element_pitch_mm;
            e = e + 1;
            p_local = [u_off; v_off; 0];
            p_w = (R * p_local)' + c_mm(:)';
            P_world(e, :) = p_w;

            vx = round((p_w(1) - origin(1)) / dx) + 1;
            vy = round((p_w(2) - origin(2)) / dy) + 1;
            vz = round((p_w(3) - origin(3)) / dz) + 1;
            vox_idx(e, :) = [vy, vx, vz];

            % Only body collision limits rotation. Out-of-grid elements are
            % handled by post-tilt grid expansion in the caller.
            if vx >= 1 && vx <= Nx && vy >= 1 && vy <= Ny && vz >= 1 && vz <= Nz
                if body(vy, vx, vz)
                    feasible = false;
                    return;
                end
            end
        end
    end
end


function theta_max = binary_search_theta_max(c_mm, axis_k, alpha, ...
        elements_per_side, element_pitch_mm, origin, dx, dy, dz, ...
        body, exclusion_zone, pml_size, Nx, Ny, Nz)
%BINARY_SEARCH_THETA_MAX Largest theta in [0, 1] for which all elements are
%   feasible at rotation rodrigues(axis_k, theta * alpha) about c_mm.
%   Returns NaN if theta = 0 is infeasible.

    R0 = eye(3);
    [ok0, ~, ~] = check_placement(c_mm, R0, ...
        elements_per_side, element_pitch_mm, origin, dx, dy, dz, ...
        body, exclusion_zone, pml_size, Nx, Ny, Nz);
    if ~ok0
        theta_max = NaN;
        return;
    end

    if abs(alpha) < 1e-6
        theta_max = 0;
        return;
    end

    R1 = rodrigues(axis_k, alpha);
    [ok1, ~, ~] = check_placement(c_mm, R1, ...
        elements_per_side, element_pitch_mm, origin, dx, dy, dz, ...
        body, exclusion_zone, pml_size, Nx, Ny, Nz);
    if ok1
        theta_max = 1;
        return;
    end

    lo = 0; hi = 1;
    for it = 1:12
        mid = 0.5 * (lo + hi);
        Rm = rodrigues(axis_k, mid * alpha);
        [okm, ~, ~] = check_placement(c_mm, Rm, ...
            elements_per_side, element_pitch_mm, origin, dx, dy, dz, ...
            body, exclusion_zone, pml_size, Nx, Ny, Nz);
        if okm
            lo = mid;
        else
            hi = mid;
        end
    end
    theta_max = lo;
end


