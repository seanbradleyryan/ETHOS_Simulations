function [sensor_mask, sensor_info] = determine_sensor_placement_fixed(config, sct_resampled, structures)
%DETERMINE_SENSOR_PLACEMENT_FIXED  Deterministic anterior sensor placement
%
%   [sensor_mask, sensor_info] = determine_sensor_placement_fixed(config, sct_resampled, structures)
%
%   PURPOSE:
%   Compute the 3D binary mask for a rigid, flat 10x10 cm ultrasound sensor
%   on the anterior abdomen using four deterministic placement rules:
%     1. Location determined using the total dose distribution.
%     2. Sensor must be on the anterior surface and fully outside the body.
%     3. Sensor superior-most edge must be anatomically inferior to the beam
%        field (beam field extent is defined by jaw_y in the transverse plane).
%     4. Sensor center is laterally aligned with the plan isocenter (X axis).
%
%   INPUTS:
%       config - Configuration struct:
%           .total_dose        - (Optional) 3D dose array (Gy). Overrides file.
%           .total_dose_file   - Path to .mat with total_rs_dose or dose_Gy var.
%           .sensor_size_cm    - [X, Z] physical dims in cm (default [10, 10])
%           .sensor_standoff_mm- Body-to-sensor gap, mm (default 5)
%           .jaw_margin_mm     - Safety margin below beam inferior edge, mm
%                                (default 10)
%           .pml_size          - PML thickness in voxels (default 10)
%           .element_size_mm   - Element patch size for averaging (default [])
%       sct_resampled - Struct with CT resampled to dose grid:
%           .bodyMask          - 3D logical (true = inside body)
%           .couchMask         - 3D logical (true = couch, optional)
%           .spacing           - [dx, dy, dz] mm
%           .origin            - [x, y, z] mm DICOM patient coords (optional)
%           .dimensions        - [ny, nx, nz] (optional, inferred from bodyMask)
%       structures - Struct with plan information:
%           .beam_metadata     - Struct array (per beam):
%               .gantry_angle  - degrees
%               .isocenter     - [x, y, z] mm
%               .jaw_x         - [x1, x2] mm at isocenter
%               .jaw_y         - [y1, y2] mm (IEC Y -> patient Z; negative =
%                                inferior jaw for HFS patients)
%           .total_dose        - (Optional) 3D dose array if not in config
%
%   OUTPUTS:
%       sensor_mask - 3D logical array (true at sensor voxels), same size as
%                     sct_resampled.bodyMask.
%       sensor_info - Diagnostic struct:
%           .sensor_y_index        - Y (AP) voxel index of the sensor plane
%           .sensor_x_range        - [ix_start, ix_end] voxel indices
%           .sensor_z_range        - [iz_start, iz_end] voxel indices
%           .sensor_center_mm      - [x, y, z] physical center of sensor
%           .isocenter_mm          - Isocenter used for lateral alignment
%           .beam_inferior_z_mm    - Inferior Z boundary of beam field (mm)
%           .z_limit_mm            - Z limit applied (beam_inferior - margin)
%           .iz_limit              - Z voxel limit for sensor superior edge
%           .iso_ix                - Isocenter X voxel index
%           .placement_valid       - true if all validation checks passed
%           .num_sensor_voxels     - Active sensor voxel count
%           .sensor_size_voxels    - [nx, nz] actual sensor extent in voxels
%           .element_map           - 2D element grouping array (if requested)
%           .num_elements          - Number of averaging elements
%           .anterior_surface      - Nx x Nz min-Y surface map (diagnostic)
%           .standoff_voxels       - Standoff in voxels
%
%   COORDINATE SYSTEM:
%       DICOM patient coords: X = left-right, Y = anterior-posterior (lower
%       index = more anterior), Z = superior-inferior.
%       Array indexing: array(Y, X, Z) = array(row, col, slice).
%       Assumes HFS convention: dz > 0, higher iz = more superior (cranial).
%       jaw_y(1) < 0 = inferior jaw, jaw_y(2) > 0 = superior jaw.
%
%   ALGORITHM:
%     1. Load total dose (config.total_dose > structures.total_dose >
%        config.total_dose_file, in that priority order).
%     2. Build anterior surface map from body mask (exclude couch).
%     3. Compute beam inferior Z boundary from beam_metadata jaw_y;
%        fall back to dose distribution if beam_metadata unavailable.
%     4. Determine isocenter X voxel for lateral centering.
%     5. Fix sensor X range centered on isocenter X (clamped to PML bounds).
%     6. Find the highest valid Z position such that sensor iz_end <= iz_limit
%        and the sensor footprint lies on the anterior surface.
%     7. Set sensor Y index at standoff from the most anterior surface point
%        within the chosen footprint.
%     8. Validate (no body overlap, Z constraint satisfied, PML clear).
%     9. Optionally partition sensor into element patches.
%
%   NOTES:
%     - Compatible output with run_standalone_simulation and
%       run_single_field_simulation (set sensor_placement_method =
%       'fixed_anterior' in CONFIG/config).
%     - When beam_metadata is unavailable, the total dose inferior Z boundary
%       (first Z slice with dose > 10% max) is used as a fallback.
%     - sct_resampled.origin defaults to [0,0,0] if not present.
%
%   EXAMPLE:
%       config.total_dose_file = '/path/to/total_rs_dose.mat';
%       config.sensor_size_cm  = [10, 10];
%       config.jaw_margin_mm   = 15;
%       [mask, info] = determine_sensor_placement_fixed(config, sct, structs);
%
%   DEPENDENCIES:
%       None beyond base MATLAB.
%
%   See also: determine_sensor_mask, run_single_field_simulation,
%             run_standalone_simulation

%% ======================== INPUT VALIDATION ========================

if ~isstruct(config)
    error('determine_sensor_placement_fixed:InvalidInput', ...
        'config must be a struct.');
end
if ~isstruct(sct_resampled) || ~isfield(sct_resampled, 'bodyMask')
    error('determine_sensor_placement_fixed:InvalidInput', ...
        'sct_resampled must be a struct with a bodyMask field.');
end
if ~isfield(sct_resampled, 'spacing') || isempty(sct_resampled.spacing)
    error('determine_sensor_placement_fixed:InvalidInput', ...
        'sct_resampled.spacing is required.');
end

%% ======================== CONFIG DEFAULTS ========================

sensor_size_cm  = get_field(config, 'sensor_size_cm',     [10, 10]);
standoff_mm     = get_field(config, 'sensor_standoff_mm', 5);
jaw_margin_mm   = get_field(config, 'jaw_margin_mm',      10);
pml_size        = get_field(config, 'pml_size',           10);
element_size_mm = get_field(config, 'element_size_mm',    []);

fprintf('[SensorFixed] Placing fixed-rule anterior sensor (%.0fx%.0f cm, standoff %.0f mm, margin %.0f mm)\n', ...
    sensor_size_cm(1), sensor_size_cm(2), standoff_mm, jaw_margin_mm);

%% ======================== STEP 1: LOAD TOTAL DOSE ========================

total_dose = load_total_dose(config, structures);

if isempty(total_dose)
    error('determine_sensor_placement_fixed:NoDose', ...
        ['Total dose not found. Supply config.total_dose, config.total_dose_file, ' ...
         'or structures.total_dose.']);
end

%% ======================== STEP 2: EXTRACT GRID INFO ========================

% Array layout: (Y, X, Z) = (rows, cols, slices)
grid_dims = size(sct_resampled.bodyMask);
Ny = grid_dims(1);  % rows  = Y (anterior-posterior)
Nx = grid_dims(2);  % cols  = X (lateral/medial, left-right)
Nz = grid_dims(3);  % slices = Z (superior-inferior)

dx = sct_resampled.spacing(1);  % mm, X
dy = sct_resampled.spacing(2);  % mm, Y
dz = sct_resampled.spacing(3);  % mm, Z

if isfield(sct_resampled, 'origin') && ~isempty(sct_resampled.origin)
    origin = sct_resampled.origin(:)';  % [x, y, z] mm
else
    origin = [0, 0, 0];
    fprintf('[SensorFixed] Warning: sct_resampled.origin not found; using [0,0,0].\n');
end

% Validate total dose matches grid
if ~isequal(size(total_dose), grid_dims)
    error('determine_sensor_placement_fixed:SizeMismatch', ...
        'Total dose [%d %d %d] does not match SCT grid [%d %d %d].', ...
        size(total_dose,1), size(total_dose,2), size(total_dose,3), ...
        Ny, Nx, Nz);
end

% Sensor size in voxels
sensor_nx = round(sensor_size_cm(1) * 10 / dx);
sensor_nz = round(sensor_size_cm(2) * 10 / dz);
sensor_nx = max(1, sensor_nx);
sensor_nz = max(1, sensor_nz);

fprintf('[SensorFixed] Grid: [%d(Y) x %d(X) x %d(Z)], spacing: [%.2f, %.2f, %.2f] mm\n', ...
    Ny, Nx, Nz, dx, dy, dz);
fprintf('[SensorFixed] Sensor size: %d(X) x %d(Z) voxels\n', sensor_nx, sensor_nz);

%% ======================== STEP 3: ANTERIOR SURFACE MAP ========================

body = sct_resampled.bodyMask;
if isfield(sct_resampled, 'couchMask') && ~isempty(sct_resampled.couchMask)
    body = body & ~sct_resampled.couchMask;
end

% anterior_surface(ix, iz) = minimum Y index where body is true
% Lower Y index = more anterior (smaller Y coordinate)
anterior_surface = NaN(Nx, Nz);
for ix = 1:Nx
    for iz = 1:Nz
        y_idx = find(body(:, ix, iz), 1, 'first');
        if ~isempty(y_idx)
            anterior_surface(ix, iz) = y_idx;
        end
    end
end

surface_valid = ~isnan(anterior_surface);
num_surface_pts = sum(surface_valid(:));

if num_surface_pts == 0
    error('determine_sensor_placement_fixed:NoSurface', ...
        'No anterior body surface found in body mask.');
end

fprintf('[SensorFixed] Anterior surface: %d valid columns\n', num_surface_pts);

%% ======================== STEP 4: BEAM INFERIOR Z BOUNDARY & ISOCENTER ========================

[beam_inferior_z_mm, isocenter_mm] = compute_beam_constraints( ...
    structures, total_dose, origin, dx, dy, dz, Ny, Nx, Nz);

fprintf('[SensorFixed] Isocenter: [%.1f, %.1f, %.1f] mm\n', isocenter_mm);
fprintf('[SensorFixed] Beam inferior Z boundary: %.1f mm\n', beam_inferior_z_mm);

% Z limit: sensor superior edge must be below this (with margin)
z_limit_mm = beam_inferior_z_mm - jaw_margin_mm;

% Convert to voxel index.
% Assumption: dz > 0 (Z increases superiorly in standard HFS DICOM).
% iz_limit = last voxel index such that z_coord <= z_limit_mm.
% z_coord(iz) = origin(3) + (iz - 1) * dz
% iz_limit = floor((z_limit_mm - origin(3)) / dz) + 1
iz_limit = floor((z_limit_mm - origin(3)) / dz) + 1;
iz_limit = min(Nz - pml_size, iz_limit);
iz_limit = max(pml_size + 1 + sensor_nz, iz_limit);  % Must fit at least one sensor

fprintf('[SensorFixed] Z limit: %.1f mm (beam_inf=%.1f, margin=%.1f) -> iz_limit=%d\n', ...
    z_limit_mm, beam_inferior_z_mm, jaw_margin_mm, iz_limit);

%% ======================== STEP 5: SENSOR X RANGE (centered on isocenter X) ========================

% Convert isocenter X to voxel index
iso_ix = round((isocenter_mm(1) - origin(1)) / dx) + 1;
iso_ix = max(1, min(Nx, iso_ix));

half_nx = floor(sensor_nx / 2);
ix_start = iso_ix - half_nx;
ix_end   = ix_start + sensor_nx - 1;

% Clamp to grid bounds (stay outside PML)
pml_margin = pml_size + 1;
if ix_start < pml_margin
    ix_start = pml_margin;
    ix_end   = ix_start + sensor_nx - 1;
end
if ix_end > Nx - pml_size
    ix_end   = Nx - pml_size;
    ix_start = max(pml_margin, ix_end - sensor_nx + 1);
end
ix_start = max(1, ix_start);
ix_end   = min(Nx, ix_end);

fprintf('[SensorFixed] Sensor X range: [%d, %d] (isocenter ix=%d)\n', ix_start, ix_end, iso_ix);

if ix_start > ix_end
    warning('determine_sensor_placement_fixed:InvalidXRange', ...
        'Sensor X range invalid after clamping. Returning empty mask.');
    [sensor_mask, sensor_info] = empty_result(grid_dims, isocenter_mm, beam_inferior_z_mm);
    return;
end

%% ======================== STEP 6: FIND Z PLACEMENT ========================
%  Sweep from iz_limit downward to find the highest valid Z position.
%  "Valid" means the full sensor footprint (ix_start:ix_end x iz_s:iz_e)
%  lies on the anterior body surface.

best_iz_start = [];

iz_search_top = iz_limit - sensor_nz + 1;  % Highest iz_start giving iz_end <= iz_limit
iz_search_top = max(pml_margin, iz_search_top);

% Primary pass: require full surface coverage
for iz_end_cand = iz_limit : -1 : pml_margin + sensor_nz - 1
    iz_start_cand = iz_end_cand - sensor_nz + 1;
    if iz_start_cand < pml_margin
        continue;
    end
    patch_valid = surface_valid(ix_start:ix_end, iz_start_cand:iz_end_cand);
    if all(patch_valid(:))
        best_iz_start = iz_start_cand;
        break;
    end
end

% Fallback: allow >= 50% surface coverage
if isempty(best_iz_start)
    fprintf('[SensorFixed] No fully-valid Z placement; relaxing to >= 50%% surface coverage.\n');
    best_coverage = 0;
    for iz_end_cand = iz_limit : -1 : pml_margin + sensor_nz - 1
        iz_start_cand = iz_end_cand - sensor_nz + 1;
        if iz_start_cand < pml_margin
            continue;
        end
        patch_valid = surface_valid(ix_start:ix_end, iz_start_cand:iz_end_cand);
        coverage = mean(patch_valid(:));
        if coverage > best_coverage
            best_coverage = coverage;
            best_iz_start = iz_start_cand;
        end
        if best_coverage >= 0.5 && iz_end_cand < iz_limit - sensor_nz
            break;  % Good enough; stop searching far from limit
        end
    end
    if ~isempty(best_iz_start)
        fprintf('[SensorFixed] Fallback placement: %.0f%% surface coverage.\n', ...
            best_coverage * 100);
    end
end

if isempty(best_iz_start)
    warning('determine_sensor_placement_fixed:NoPlacement', ...
        'Cannot find valid Z placement inferior to beam field (iz_limit=%d). ' ...
        'Check jaw_margin_mm or beam_metadata.', iz_limit);
    [sensor_mask, sensor_info] = empty_result(grid_dims, isocenter_mm, beam_inferior_z_mm);
    sensor_info.iz_limit = iz_limit;
    sensor_info.z_limit_mm = z_limit_mm;
    return;
end

iz_start = best_iz_start;
iz_end   = min(iz_start + sensor_nz - 1, Nz);

sensor_x_range = [ix_start, ix_end];
sensor_z_range = [iz_start, iz_end];

sensor_z_max_mm = origin(3) + (iz_end - 1) * dz;
fprintf('[SensorFixed] Z placement: iz=[%d, %d] (z_max=%.1f mm, z_limit=%.1f mm)\n', ...
    iz_start, iz_end, sensor_z_max_mm, z_limit_mm);

%% ======================== STEP 7: SENSOR Y INDEX ========================

surface_patch = anterior_surface(sensor_x_range(1):sensor_x_range(2), ...
                                  sensor_z_range(1):sensor_z_range(2));
valid_in_patch = ~isnan(surface_patch);

if any(valid_in_patch(:))
    min_anterior_y = min(surface_patch(valid_in_patch));
else
    % Fallback to global minimum
    min_anterior_y = min(anterior_surface(surface_valid));
    fprintf('[SensorFixed] Warning: no valid surface in sensor footprint; using global minimum.\n');
end

standoff_voxels = ceil(standoff_mm / dy);
sensor_y_index  = max(pml_margin, round(min_anterior_y) - standoff_voxels);

sensor_center_x = origin(1) + (mean(sensor_x_range) - 1) * dx;
sensor_center_y = origin(2) + (sensor_y_index - 1) * dy;
sensor_center_z = origin(3) + (mean(sensor_z_range) - 1) * dz;

fprintf('[SensorFixed] Sensor Y index: %d (surface min Y: %d, standoff: %d voxels)\n', ...
    sensor_y_index, round(min_anterior_y), standoff_voxels);
fprintf('[SensorFixed] Sensor center (mm): [%.1f, %.1f, %.1f]\n', ...
    sensor_center_x, sensor_center_y, sensor_center_z);

%% ======================== STEP 8: CREATE SENSOR MASK ========================

sensor_mask = false(grid_dims);
sensor_mask(sensor_y_index, ...
            sensor_x_range(1):sensor_x_range(2), ...
            sensor_z_range(1):sensor_z_range(2)) = true;

fprintf('[SensorFixed] Raw sensor: %d voxels\n', sum(sensor_mask(:)));

%% ======================== STEP 9: VALIDATE ========================

placement_valid = true;

% Check 1: No sensor voxels inside body
overlap = sensor_mask & body;
if any(overlap(:))
    n_overlap = sum(overlap(:));
    warning('determine_sensor_placement_fixed:BodyOverlap', ...
        'Sensor overlaps body at %d voxels. Removing overlapping voxels.', n_overlap);
    sensor_mask = sensor_mask & ~body;
    placement_valid = false;
end

% Check 2: Z constraint satisfied (sensor superior edge < beam inferior edge)
if sensor_z_max_mm >= beam_inferior_z_mm
    warning('determine_sensor_placement_fixed:ZConstraintViolation', ...
        'Sensor Z max (%.1f mm) >= beam inferior Z (%.1f mm). Constraint violated.', ...
        sensor_z_max_mm, beam_inferior_z_mm);
    placement_valid = false;
end

% Check 3: PML clearance
if sensor_x_range(1) <= pml_size || sensor_x_range(2) > Nx - pml_size || ...
   sensor_z_range(1) <= pml_size || sensor_z_range(2) > Nz - pml_size || ...
   sensor_y_index    <= pml_size
    warning('determine_sensor_placement_fixed:PMLOverlap', ...
        'Sensor placement extends into PML boundary region (pml_size=%d).', pml_size);
    placement_valid = false;
end

num_sensor_voxels = sum(sensor_mask(:));
fprintf('[SensorFixed] Final sensor: %d voxels (valid: %s)\n', ...
    num_sensor_voxels, mat2str(placement_valid));

%% ======================== STEP 10: ELEMENT GROUPINGS ========================

element_map  = [];
num_elements = 0;

if ~isempty(element_size_mm) && element_size_mm > 0
    elem_nx = max(1, floor(element_size_mm / dx));
    elem_nz = max(1, floor(element_size_mm / dz));

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

    % Remove elements for voxels removed during body-overlap validation
    for lix = 1:local_nx
        for liz = 1:local_nz
            gix = sensor_x_range(1) + lix - 1;
            giz = sensor_z_range(1) + liz - 1;
            if ~sensor_mask(sensor_y_index, gix, giz)
                element_map(lix, liz) = 0;
            end
        end
    end

    fprintf('[SensorFixed] Element grouping: %d elements (%dx%d voxels each)\n', ...
        num_elements, elem_nx, elem_nz);
end

%% ======================== PACK OUTPUT ========================

sensor_info = struct();
sensor_info.sensor_y_index      = sensor_y_index;
sensor_info.sensor_x_range      = sensor_x_range;
sensor_info.sensor_z_range      = sensor_z_range;
sensor_info.sensor_center_mm    = [sensor_center_x, sensor_center_y, sensor_center_z];
sensor_info.isocenter_mm        = isocenter_mm;
sensor_info.beam_inferior_z_mm  = beam_inferior_z_mm;
sensor_info.z_limit_mm          = z_limit_mm;
sensor_info.iz_limit            = iz_limit;
sensor_info.iso_ix              = iso_ix;
sensor_info.element_map         = element_map;
sensor_info.num_elements        = num_elements;
sensor_info.placement_valid     = placement_valid;
sensor_info.num_sensor_voxels   = num_sensor_voxels;
sensor_info.sensor_size_voxels  = [sensor_x_range(2)-sensor_x_range(1)+1, ...
                                    sensor_z_range(2)-sensor_z_range(1)+1];
sensor_info.anterior_surface    = anterior_surface;
sensor_info.standoff_voxels     = standoff_voxels;

fprintf('[SensorFixed] Placement complete.\n');

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function total_dose = load_total_dose(config, structures)
%LOAD_TOTAL_DOSE  Retrieve total dose from config or structures (priority order).
%  Priority: config.total_dose > structures.total_dose > config.total_dose_file

    total_dose = [];

    % 1. Direct array in config
    if isfield(config, 'total_dose') && ~isempty(config.total_dose) && ...
            isnumeric(config.total_dose)
        total_dose = double(config.total_dose);
        fprintf('[SensorFixed] Total dose: loaded from config.total_dose (%d voxels).\n', ...
            numel(total_dose));
        return;
    end

    % 2. Direct array in structures
    if isstruct(structures) && isfield(structures, 'total_dose') && ...
            ~isempty(structures.total_dose) && isnumeric(structures.total_dose)
        total_dose = double(structures.total_dose);
        fprintf('[SensorFixed] Total dose: loaded from structures.total_dose (%d voxels).\n', ...
            numel(total_dose));
        return;
    end

    % 3. File path in config
    if isfield(config, 'total_dose_file') && ~isempty(config.total_dose_file)
        fpath = config.total_dose_file;
        if ~isfile(fpath)
            warning('determine_sensor_placement_fixed:DoseFileNotFound', ...
                'total_dose_file not found: %s', fpath);
            return;
        end
        try
            d = load(fpath);
            fnames = fieldnames(d);
            if isfield(d, 'total_rs_dose')
                total_dose = double(d.total_rs_dose);
            elseif isfield(d, 'dose_Gy')
                total_dose = double(d.dose_Gy);
            elseif length(fnames) == 1
                total_dose = double(d.(fnames{1}));
            else
                % Use first 3D numeric field
                for k = 1:length(fnames)
                    v = d.(fnames{k});
                    if isnumeric(v) && ndims(v) == 3
                        total_dose = double(v);
                        break;
                    end
                end
            end
            if ~isempty(total_dose)
                fprintf('[SensorFixed] Total dose: loaded from file %s.\n', fpath);
            else
                warning('determine_sensor_placement_fixed:DoseVarNotFound', ...
                    'Could not identify dose variable in %s.', fpath);
            end
        catch ME
            warning('determine_sensor_placement_fixed:DoseLoadFail', ...
                'Failed to load dose file %s: %s', fpath, ME.message);
        end
    end
end


function [beam_inferior_z_mm, isocenter_mm] = compute_beam_constraints( ...
        structures, total_dose, origin, dx, dy, dz, Ny, Nx, Nz)
%COMPUTE_BEAM_CONSTRAINTS  Extract beam inferior Z limit and isocenter from
%  beam_metadata. Falls back to total dose distribution if metadata is absent.
%
%  Assumes HFS DICOM convention: jaw_y(1) < 0 = inferior jaw (caudal),
%  jaw_y(2) > 0 = superior jaw (cranial). dz > 0 (higher iz = more cranial).

    beam_inferior_z_mm = NaN;
    isocenter_mm       = [];

    has_beam_meta = isstruct(structures) && ...
                    isfield(structures, 'beam_metadata') && ...
                    ~isempty(structures.beam_metadata);

    if has_beam_meta
        bm = structures.beam_metadata;

        % Isocenter: use the first beam with a valid isocenter
        for b = 1:length(bm)
            if isfield(bm(b), 'isocenter') && numel(bm(b).isocenter) >= 3
                isocenter_mm = bm(b).isocenter(:)';
                break;
            end
        end

        % Beam inferior Z: minimum of (iso_z + min(jaw_y)) across all beams
        min_inf_z = Inf;
        for b = 1:length(bm)
            if ~isfield(bm(b), 'isocenter') || isempty(bm(b).isocenter)
                continue;
            end
            if ~isfield(bm(b), 'jaw_y') || isempty(bm(b).jaw_y)
                continue;
            end
            iso_z    = bm(b).isocenter(3);
            inf_edge = iso_z + min(bm(b).jaw_y(:));  % most inferior jaw edge
            if inf_edge < min_inf_z
                min_inf_z = inf_edge;
            end
        end

        if ~isinf(min_inf_z)
            beam_inferior_z_mm = min_inf_z;
            fprintf('[SensorFixed] Beam inferior Z from jaw_y: %.1f mm\n', ...
                beam_inferior_z_mm);
        end
    end

    % Fallback: derive from total dose distribution
    if isnan(beam_inferior_z_mm) || isempty(isocenter_mm)
        fprintf('[SensorFixed] Using total dose fallback for beam constraints.\n');

        dose = total_dose;
        dose_max = max(dose(:));

        if dose_max > 0
            thresh = 0.10 * dose_max;

            % Dose-weighted centroid as isocenter estimate
            dose_thresh = dose;
            dose_thresh(dose_thresh < thresh) = 0;
            total_d = sum(dose_thresh(:));

            if total_d > 0
                [Ny_d, Nx_d, Nz_d] = size(dose_thresh);
                x_vec = origin(1) + (0:Nx_d-1) * dx;
                y_vec = origin(2) + (0:Ny_d-1) * dy;
                z_vec = origin(3) + (0:Nz_d-1) * dz;

                dose_x = squeeze(sum(sum(dose_thresh, 1), 3));
                dose_y = squeeze(sum(sum(dose_thresh, 2), 3));
                dose_z = squeeze(sum(sum(dose_thresh, 1), 2));

                cx = sum(x_vec(:) .* dose_x(:)) / total_d;
                cy = sum(y_vec(:) .* dose_y(:)) / total_d;
                cz = sum(z_vec(:) .* dose_z(:)) / total_d;

                if isempty(isocenter_mm)
                    isocenter_mm = [cx, cy, cz];
                    fprintf('[SensorFixed] Isocenter fallback (dose centroid): [%.1f, %.1f, %.1f] mm\n', ...
                        isocenter_mm);
                end

                if isnan(beam_inferior_z_mm)
                    % Inferior boundary: first Z slice (from inferior end)
                    % where summed dose exceeds threshold.
                    % With dz > 0: iz=1 is most inferior.
                    dose_z_profile = squeeze(sum(sum(dose_thresh, 1), 2));
                    iz_inf = find(dose_z_profile > 0, 1, 'first');
                    if ~isempty(iz_inf)
                        beam_inferior_z_mm = origin(3) + (iz_inf - 1) * dz;
                    else
                        beam_inferior_z_mm = cz - 50;  % 5 cm below centroid fallback
                    end
                    fprintf('[SensorFixed] Beam inferior Z fallback (dose): %.1f mm\n', ...
                        beam_inferior_z_mm);
                end
            end
        end

        % Last resort
        if isempty(isocenter_mm)
            isocenter_mm = origin + [(Nx/2)*dx, (Ny/2)*dy, (Nz/2)*dz];
            fprintf('[SensorFixed] Isocenter last-resort (grid center): [%.1f, %.1f, %.1f] mm\n', ...
                isocenter_mm);
        end
        if isnan(beam_inferior_z_mm)
            beam_inferior_z_mm = isocenter_mm(3) - 50;
            fprintf('[SensorFixed] Beam inferior Z last-resort (iso-50mm): %.1f mm\n', ...
                beam_inferior_z_mm);
        end
    end
end


function [sensor_mask, sensor_info] = empty_result(grid_dims, isocenter_mm, beam_inferior_z_mm)
%EMPTY_RESULT  Return empty sensor mask and default sensor_info struct.

    sensor_mask = false(grid_dims);
    sensor_info = struct();
    sensor_info.sensor_y_index      = 0;
    sensor_info.sensor_x_range      = [0, 0];
    sensor_info.sensor_z_range      = [0, 0];
    sensor_info.sensor_center_mm    = [0, 0, 0];
    sensor_info.placement_valid     = false;
    sensor_info.num_sensor_voxels   = 0;
    sensor_info.sensor_size_voxels  = [0, 0];
    sensor_info.element_map         = [];
    sensor_info.num_elements        = 0;
    sensor_info.anterior_surface    = [];
    sensor_info.standoff_voxels     = 0;
    sensor_info.iz_limit            = 0;
    sensor_info.z_limit_mm          = 0;
    sensor_info.iso_ix              = 0;

    if nargin >= 2 && ~isempty(isocenter_mm)
        sensor_info.isocenter_mm = isocenter_mm;
    else
        sensor_info.isocenter_mm = [0, 0, 0];
    end
    if nargin >= 3 && ~isnan(beam_inferior_z_mm)
        sensor_info.beam_inferior_z_mm = beam_inferior_z_mm;
    else
        sensor_info.beam_inferior_z_mm = 0;
    end
end


function val = get_field(s, name, default)
%GET_FIELD  Retrieve struct field with fallback default.
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default;
    end
end
