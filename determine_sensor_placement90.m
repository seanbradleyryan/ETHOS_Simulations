function [sensor_mask, sensor_info] = determine_sensor_placement90(sct_resampled, field_dose, beam_metadata, config)
%DETERMINE_SENSOR_PLACEMENT90  Flat sensor at 90° to beam axis in transverse plane
%
%   [sensor_mask, sensor_info] = determine_sensor_placement90(sct_resampled, field_dose, beam_metadata, config)
%
%   PURPOSE:
%   Computes the 3D binary mask for a rigid flat 10×10 cm ultrasound sensor
%   whose normal vector is rotated 90° CCW from the beam source direction in
%   the transverse plane. As the gantry angle changes across fields, the sensor
%   tracks around the outside of the body, always positioned on the lateral/
%   oblique side relative to the current beam.
%
%   SENSOR GEOMETRY (transverse plane view):
%   - The sensor appears as a STRAIGHT LINE in the transverse plane.
%   - The line's normal is at 90° to the beam axis (normal = 90° CCW of beam-
%     source direction in physical XY).
%   - The line runs PARALLEL to the beam axis direction in the transverse plane.
%   - The line is centered on the beam isocenter (equal extent on both sides).
%   - The line is displaced outward (in the normal direction) until outside
%     the body surface, plus a standoff gap.
%
%   COORDINATE SYSTEM (same as determine_sensor_mask):
%       DICOM patient coords: X = left, Y = posterior, Z = superior.
%       Array indexing:       array(iy, ix, iz) = array(Y_row, X_col, Z_slice).
%       spacing = [dx, dy, dz] mm  (X, Y, Z respectively).
%
%   GANTRY ANGLE CONVENTION (Varian IEC 61217, supine patient):
%       0°   = source ANTERIOR  (source at −Y, beam traveling +Y)
%      90°   = source PATIENT-LEFT  (source at +X)
%     180°   = source POSTERIOR (source at +Y)
%     270°   = source PATIENT-RIGHT (source at −X)
%
%   SENSOR NORMAL DIRECTION (90° CCW rotation of source direction in XY):
%       0°   gantry → normal toward PATIENT-LEFT   (+X, snorm = (1, 0))
%      90°   gantry → normal toward POSTERIOR       (+Y, snorm = (0, 1))
%     180°   gantry → normal toward PATIENT-RIGHT   (-X, snorm = (-1, 0))
%     270°   gantry → normal toward ANTERIOR        (-Y, snorm = (0, -1))
%
%   FLAT SENSOR GUARANTEE:
%   The algorithm finds the maximum body-surface offset along the full
%   sensor span (all positions along the beam direction). All sensor voxels
%   use that uniform offset, ensuring the sensor is coplanar and no voxel
%   is inside the body.
%
%   INPUTS:
%       sct_resampled - Struct (same interface as determine_sensor_mask):
%           .bodyMask       - 3D logical (true = inside body). Array(Y, X, Z).
%           .couchMask      - 3D logical (true = couch region).
%           .origin         - [x, y, z] mm — DICOM ImagePositionPatient
%           .spacing        - [dx, dy, dz] mm — voxel sizes in X, Y, Z
%           .dimensions     - [Ny, Nx, Nz] — grid size (rows, cols, slices)
%       field_dose - Struct:
%           .gantry_angle   - Gantry angle for this field (degrees)
%           .isocenter      - [x, y, z] mm (DICOM patient coords, from RTPLAN)
%           .dose_Gy        - 3D dose array (Gy) — used for centroid fallback
%       beam_metadata - Struct array (not used; kept for API compatibility with
%                       determine_sensor_mask)
%       config - Struct with optional fields:
%           .sensor_size_cm      - [in-plane_cm, Z_cm] (default: [10, 10])
%           .sensor_standoff_mm  - Gap between body surface and sensor (default: 5)
%           .element_size_mm     - Patch size for element averaging (default: [])
%           .pml_size            - PML thickness in voxels (default: 10)
%
%   OUTPUTS:
%       sensor_mask - 3D logical array, same size as sct_resampled.bodyMask.
%                     true = sensor voxel.
%       sensor_info - Struct with diagnostic fields:
%           .gantry_angle         - Input gantry angle (degrees)
%           .sensor_normal_deg    - Sensor normal direction angle (deg from +X in XY)
%           .t_body_exit_mm       - Max body-surface distance from isocenter along normal (mm)
%           .t_sensor_mm          - Total sensor offset from isocenter along normal (mm)
%           .sensor_center_mm     - [x, y, z] physical center of sensor (mm)
%           .iz_start, .iz_end    - Z slice range of sensor in array
%           .sensor_size_voxels   - [n_inplane_approx, n_z] approximate sensor extents
%           .num_sensor_voxels    - Total voxels placed
%           .element_map          - 2D element index array (or [])
%           .num_elements         - Number of signal-averaged elements
%           .placement_valid      - true if mask is non-empty and body-free
%           .snorm_physical       - [snorm_X, snorm_Y] unit normal in DICOM XY
%           .inplane_physical     - [inplane_X, inplane_Y] unit in-plane dir in DICOM XY
%
%   EXAMPLE:
%       config.sensor_size_cm     = [10, 10];
%       config.sensor_standoff_mm = 5;
%       [mask, info] = determine_sensor_placement90(sct_resampled, field_dose, [], config);
%
%   DEPENDENCIES:
%       None (no toolbox required)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE:   April 2026
%   VERSION: 1.0
%
%   See also: determine_sensor_mask, run_single_field_simulation

%% ======================== CONFIG DEFAULTS ========================

sensor_size_cm  = get_field(config, 'sensor_size_cm',     [10, 10]);
standoff_mm     = get_field(config, 'sensor_standoff_mm', 5);
element_size_mm = get_field(config, 'element_size_mm',    []);
pml_size        = get_field(config, 'pml_size',           10);

sensor_w_mm = sensor_size_cm(1) * 10;   % in-plane width  (beam direction), mm
sensor_z_mm = sensor_size_cm(2) * 10;   % height (Z / superior-inferior),   mm

gantry_deg = field_dose.gantry_angle;
ga_rad     = gantry_deg * pi / 180;

fprintf('        [Sensor90] Gantry %.1f°, sensor %.0fx%.0f cm, standoff %.0f mm\n', ...
    gantry_deg, sensor_size_cm(1), sensor_size_cm(2), standoff_mm);

%% ======================== GRID INFO ========================
% Array convention: array(iy, ix, iz) = array(Y_row, X_col, Z_slice)
% spacing = [dx, dy, dz] where dx=X spacing (mm), dy=Y spacing (mm), dz=Z spacing (mm)

Ny = sct_resampled.dimensions(1);
Nx = sct_resampled.dimensions(2);
Nz = sct_resampled.dimensions(3);

dx = sct_resampled.spacing(1);   % mm per voxel, X (cols)
dy = sct_resampled.spacing(2);   % mm per voxel, Y (rows)
dz = sct_resampled.spacing(3);   % mm per voxel, Z (slices)

ox = sct_resampled.origin(1);    % DICOM X origin (mm)
oy = sct_resampled.origin(2);    % DICOM Y origin (mm)
oz = sct_resampled.origin(3);    % DICOM Z origin (mm)

% Body mask excluding couch (same treatment as determine_sensor_mask)
body_mask = logical(sct_resampled.bodyMask) & ~logical(sct_resampled.couchMask);

fprintf('        [Sensor90] Grid [Ny=%d, Nx=%d, Nz=%d], spacing [%.2f, %.2f, %.2f] mm\n', ...
    Ny, Nx, Nz, dx, dy, dz);

%% ======================== BEAM GEOMETRY IN PHYSICAL XY ========================
%
% Varian IEC 61217 gantry convention (supine patient, from determine_sensor_mask):
%   source_y_mm = iso_y - SAD * cos(ga)
%   → Source direction from isocenter: src = (sin(ga), -cos(ga)) in (X, Y)
%     ga=0:   src=(0, -1) = anterior   (−Y)  ✓
%     ga=90:  src=(1,  0) = pat-left   (+X)  ✓
%     ga=180: src=(0, +1) = posterior  (+Y)  ✓
%     ga=270: src=(-1, 0) = pat-right  (−X)  ✓
%
% Sensor normal = 90° CCW rotation of src in XY: (a,b) → (−b, a)
%   snorm = (cos(ga), sin(ga))
%     ga=0:   snorm=(1,  0) = pat-left   (+X)
%     ga=90:  snorm=(0,  1) = posterior  (+Y)
%     ga=180: snorm=(-1, 0) = pat-right  (-X)
%     ga=270: snorm=(0, -1) = anterior   (-Y)
%
% Sensor in-plane direction = beam direction (from source toward iso) = −src:
%   inplane = (−sin(ga), cos(ga)) in (X, Y)
%     ga=0:   inplane=(0,  1) = +Y (posterior)   → sensor line runs A-P
%     ga=90:  inplane=(-1, 0) = -X (pat-right)   → sensor line runs L-R
%     ga=180: inplane=(0, -1) = -Y (anterior)    → sensor line runs A-P
%     ga=270: inplane=(1,  0) = +X (pat-left)    → sensor line runs L-R

snorm_X   =  cos(ga_rad);    % sensor normal, DICOM X component
snorm_Y   =  sin(ga_rad);    % sensor normal, DICOM Y component

inplane_X = -sin(ga_rad);    % in-plane direction, DICOM X component
inplane_Y =  cos(ga_rad);    % in-plane direction, DICOM Y component

sensor_normal_deg = atan2d(snorm_Y, snorm_X);

fprintf('        [Sensor90] Snorm: (X=%.3f, Y=%.3f) → %.1f° from +X\n', ...
    snorm_X, snorm_Y, sensor_normal_deg);
fprintf('        [Sensor90] In-plane: (X=%.3f, Y=%.3f)\n', inplane_X, inplane_Y);

%% ======================== ISOCENTER ========================

iso_mm = field_dose.isocenter(:)';
iso_X  = iso_mm(1);
iso_Y  = iso_mm(2);
iso_Z  = iso_mm(3);

% Isocenter in array indices (may be non-integer)
iso_ix = (iso_X - ox) / dx + 1;
iso_iy = (iso_Y - oy) / dy + 1;
iso_iz = (iso_Z - oz) / dz + 1;

% Clamp to grid
iso_ix = max(1, min(Nx, iso_ix));
iso_iy = max(1, min(Ny, iso_iy));
iso_iz = max(1, min(Nz, iso_iz));

fprintf('        [Sensor90] Isocenter: (X=%.1f, Y=%.1f, Z=%.1f) mm  →  grid (iy=%.1f, ix=%.1f, iz=%.1f)\n', ...
    iso_X, iso_Y, iso_Z, iso_iy, iso_ix, iso_iz);

%% ======================== SENSOR Z RANGE ========================

iz_center = round(iso_iz);
half_z_vox = floor((sensor_z_mm / dz) / 2);

iz_start = max(pml_size + 1, iz_center - half_z_vox);
iz_end   = min(Nz - pml_size, iz_center + half_z_vox);

fprintf('        [Sensor90] Z range: slices [%d, %d] = %.0f mm\n', ...
    iz_start, iz_end, (iz_end - iz_start) * dz);

%% ======================== FIND BODY SURFACE OFFSET ========================
%
% Walk outward from the isocenter in the snorm direction, starting at each
% (u) position along the in-plane (beam) direction.  Track the maximum body-
% exit distance; the flat sensor is placed at that uniform offset so every
% sensor voxel is guaranteed to be outside the body.
%
% Scan resolution: 1 mm steps along normal direction.
% u samples: enough to cover the full sensor width at sub-voxel density.

u_half_mm    = sensor_w_mm / 2;           % ±50 mm from isocenter
n_u_scan     = max(7, ceil(sensor_w_mm / min(dx, dy)));
u_scan_mm    = linspace(-u_half_mm, u_half_mm, n_u_scan);

max_search_mm = max(Nx * dx, Ny * dy) * 1.2;   % generous upper bound
step_scan_mm  = 1.0;                            % 1 mm walk steps

max_t_exit_mm = 0;

for ui = 1:numel(u_scan_mm)
    u = u_scan_mm(ui);

    % Base position in physical XY, displaced along in-plane direction from iso
    base_X = iso_X + u * inplane_X;
    base_Y = iso_Y + u * inplane_Y;

    in_body = false;
    t_exit  = 0;

    for t_mm = 0 : step_scan_mm : max_search_mm
        x_test = base_X + t_mm * snorm_X;
        y_test = base_Y + t_mm * snorm_Y;

        ix_t = round((x_test - ox) / dx) + 1;
        iy_t = round((y_test - oy) / dy) + 1;

        % Off the grid → treat as past body boundary
        if ix_t < 1 || ix_t > Nx || iy_t < 1 || iy_t > Ny
            t_exit = t_mm;
            break;
        end

        % Sample body over a few Z slices centred on isocenter
        z_lo   = max(1,  iz_center - 2);
        z_hi   = min(Nz, iz_center + 2);
        at_body = any(any(body_mask(iy_t, ix_t, z_lo:z_hi)));

        if at_body
            in_body = true;
        elseif in_body
            % First transition: body → outside
            t_exit = t_mm;
            break;
        end
        % If we never entered body, t_exit stays 0 (base is already outside)
    end

    if t_exit > max_t_exit_mm
        max_t_exit_mm = t_exit;
    end
end

t_sensor_mm = max_t_exit_mm + standoff_mm;

fprintf('        [Sensor90] Max body exit: %.1f mm | standoff: %.1f mm | total offset: %.1f mm\n', ...
    max_t_exit_mm, standoff_mm, t_sensor_mm);

%% ======================== BUILD SENSOR MASK ========================
%
% Step along the in-plane (beam) direction at sub-voxel resolution to avoid
% gaps in the oblique sensor plane.  For each u and each Z slice, compute the
% physical sensor position and mark the nearest voxel.
%
% Sub-voxel step: half the minimum voxel dimension ensures no gaps even at 45°.

step_u_mm  = min(dx, dy) * 0.5;
u_range_mm = -u_half_mm : step_u_mm : u_half_mm;

sensor_mask         = false(Ny, Nx, Nz);
placed_count        = 0;
body_conflict_count = 0;

for u = u_range_mm

    % Physical base position along sensor span (in-plane from isocenter)
    base_X = iso_X + u * inplane_X;
    base_Y = iso_Y + u * inplane_Y;

    % Physical sensor point: base + uniform normal offset
    x_sens = base_X + t_sensor_mm * snorm_X;
    y_sens = base_Y + t_sensor_mm * snorm_Y;

    % Array indices (round to nearest voxel)
    ix = round((x_sens - ox) / dx) + 1;
    iy = round((y_sens - oy) / dy) + 1;

    % Skip if outside grid or inside PML
    if ix < pml_size + 1 || ix > Nx - pml_size || ...
       iy < pml_size + 1 || iy > Ny - pml_size
        continue;
    end

    % Place sensor voxels across the full Z range
    for iz = iz_start:iz_end

        if iz < 1 || iz > Nz, continue; end

        % Body conflict guard (should not happen for a well-placed flat sensor,
        % but handles curved body surfaces and numerical edge cases)
        if body_mask(iy, ix, iz)
            body_conflict_count = body_conflict_count + 1;
            continue;
        end

        if ~sensor_mask(iy, ix, iz)
            sensor_mask(iy, ix, iz) = true;
            placed_count = placed_count + 1;
        end

    end
end

fprintf('        [Sensor90] Placed %d sensor voxels', placed_count);
if body_conflict_count > 0
    fprintf(' (%d body-conflict voxels skipped)', body_conflict_count);
end
fprintf('\n');

%% ======================== VALIDATION ========================

placement_valid = (placed_count > 0) && (body_conflict_count == 0);

% Final body-overlap check (defensive)
overlap = sensor_mask & body_mask;
if any(overlap(:))
    n_overlap = sum(overlap(:));
    warning('determine_sensor_placement90:BodyOverlap', ...
        'Sensor overlaps body at %d voxels. Removing.', n_overlap);
    sensor_mask    = sensor_mask & ~body_mask;
    placement_valid = false;
end

if placed_count == 0
    warning('determine_sensor_placement90:EmptySensor', ...
        'Sensor mask is empty (gantry %.1f°). Check isocenter, body mask, and grid bounds.', ...
        gantry_deg);
end

%% ======================== ELEMENT GROUPINGS ========================
%
% Partition the sensor into rectangular patches for signal averaging.
% Uses (u_bin, z_bin) coordinates.  u bins span the in-plane direction;
% z bins span the Z (S-I) direction.

element_map  = [];
num_elements = 0;

if ~isempty(element_size_mm) && element_size_mm > 0 && placed_count > 0

    n_u_bins = max(1, floor(sensor_w_mm / element_size_mm));
    n_z_bins = max(1, floor((iz_end - iz_start + 1) * dz / element_size_mm));
    num_elements = n_u_bins * n_z_bins;

    element_map = zeros(n_u_bins, n_z_bins, 'uint16');
    for eu = 1:n_u_bins
        for ez = 1:n_z_bins
            element_map(eu, ez) = uint16((eu - 1) * n_z_bins + ez);
        end
    end

    fprintf('        [Sensor90] Element grouping: %d elements (%d u-bins × %d z-bins, %.1f mm each)\n', ...
        num_elements, n_u_bins, n_z_bins, element_size_mm);
end

%% ======================== SENSOR CENTRE (physical) ========================

x_ctr = iso_X + t_sensor_mm * snorm_X;
y_ctr = iso_Y + t_sensor_mm * snorm_Y;
z_ctr = iso_Z;
sensor_center_mm = [x_ctr, y_ctr, z_ctr];

fprintf('        [Sensor90] Sensor centre: (X=%.1f, Y=%.1f, Z=%.1f) mm\n', ...
    sensor_center_mm(1), sensor_center_mm(2), sensor_center_mm(3));
fprintf('        [Sensor90] Placement valid: %s\n', mat2str(placement_valid));

%% ======================== PACK OUTPUT ========================

sensor_info = struct();
sensor_info.gantry_angle        = gantry_deg;
sensor_info.sensor_normal_deg   = sensor_normal_deg;
sensor_info.t_body_exit_mm      = max_t_exit_mm;
sensor_info.t_sensor_mm         = t_sensor_mm;
sensor_info.sensor_center_mm    = sensor_center_mm;
sensor_info.iz_start            = iz_start;
sensor_info.iz_end              = iz_end;
sensor_info.sensor_size_voxels  = [numel(u_range_mm), iz_end - iz_start + 1];
sensor_info.num_sensor_voxels   = sum(sensor_mask(:));
sensor_info.element_map         = element_map;
sensor_info.num_elements        = num_elements;
sensor_info.placement_valid     = placement_valid;
sensor_info.snorm_physical      = [snorm_X,   snorm_Y];
sensor_info.inplane_physical    = [inplane_X, inplane_Y];

end


%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function val = get_field(s, name, default)
%GET_FIELD  Retrieve struct field with fallback to default (same as determine_sensor_mask)
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default;
    end
end
