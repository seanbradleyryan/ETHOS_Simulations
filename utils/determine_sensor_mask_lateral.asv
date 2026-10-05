function [sensor_mask, sensor_info] = determine_sensor_mask_lateral(sct_resampled, field_dose, beam_metadata, config)
%DETERMINE_SENSOR_MASK_LATERAL  Place a 2D ultrasound array on the patient's side
%
%   [sensor_mask, sensor_info] = determine_sensor_mask_lateral(sct_resampled, field_dose, beam_metadata, config)
%
%   PURPOSE:
%   Lateral counterpart to DETERMINE_SENSOR_MASK. Instead of pressing the rigid
%   NxN element array against the ANTERIOR abdomen, this places it against the
%   patient's RIGHT (default) or LEFT flank, snaps it to the body surface at a
%   fixed standoff, and tilts it to point as directly as possible at the
%   TOTAL-plan dose centroid without entering the beam/dose zone. Same inputs
%   and outputs as determine_sensor_mask, so it is a drop-in alternative.
%
%   WHY A SEPARATE FILE:
%   determine_sensor_mask is calibrated end-to-end for an anterior sensor: its
%   surface map is the min-Y skin, its flat panel is a constant-Y slab, its
%   panel normal is coronal [0,-1,0], and its grid-expansion fallback slides
%   the sensor inferior. A lateral sensor changes all four (constant-X slab,
%   min/max-X skin, panel normal +/-X, expansion mostly in X), so the geometry
%   is re-derived here rather than bolted on with flags.
%
%   PRIMARY BEHAVIOUR (implemented, robust):
%     1. Side: 'right' (default) or 'left' from config.sensor_side. Right =
%        lower X index; left = higher X index (HFS convention below).
%     2. Snap to the lateral body surface at config.sensor_standoff_mm.
%     3. Aim the panel at the dose centroid of the SUMMED PLAN DOSE (all
%        fields/segments), not the single field passed in. Tilt magnitude is
%        the largest that keeps every element outside the body (binary search),
%        exactly like the anterior version's iso-aim.
%     4. Avoid the beam/dose zone: the sensor footprint stays off the (Y,Z)
%        shadow of the summed dose (columns reaching >=10% of peak along X).
%     5. Physical elements: discrete NxN grid with pitch/size/kerf, voxelised
%        with a multi-voxel active footprint when flat and one-voxel-per-element
%        when tilted (identical model to determine_sensor_mask).
%     6. Grid expansion: if the standoff or the tilt pushes any element past a
%        grid edge in X, Y, or Z, the grid is water-padded to fit, and the
%        padding is reported in sensor_info.grid_pad for the caller (same
%        contract as determine_sensor_mask).
%
%   SECONDARY BEHAVIOUR (best-effort, config-gated, needs sct_resampled.cubeHU):
%     A. Bone avoidance (config.avoid_bone, default true when cubeHU present):
%        among near-optimal placements, prefer the one whose sensor->centroid
%        ray passes through the least bone. This nudges the panel away from the
%        ribs (too superior) and the ilium (too inferior) toward a soft-tissue
%        acoustic window, without overriding the primary distance objective.
%     B. Side auto-selection (config.auto_select_side, default false): pick the
%        side whose straight centroid->surface ray crosses the least intra-body
%        gas. Off by default; config.sensor_side wins when auto-select is off.
%
%   INPUTS:
%       sct_resampled - Struct with CT resampled to the dose grid:
%           .bodyMask   - 3D logical (true = inside body). array(Y, X, Z).
%           .couchMask  - 3D logical (true = couch). Optional.
%           .cubeHU     - 3D HU volume (same grid). Optional; enables secondary
%                         bone/gas features when present.
%           .origin     - [x, y, z] mm (DICOM ImagePositionPatient)
%           .spacing    - [dx, dy, dz] mm
%           .dimensions - [Ny, Nx, Nz] (rows, cols, slices). Inferred if absent.
%       field_dose - Struct with:
%           .dose_Gy       - 3D per-field dose (Gy). Fallback placement dose.
%           .total_dose_Gy - (Preferred) summed plan dose on the same grid, so
%                            placement is deterministic across all fields.
%           .gantry_angle  - degrees (logged only).
%           .origin/.spacing/.dimensions - as above (optional).
%       beam_metadata - Struct array (accepted for API compatibility; the
%                       isocenter is logged if present but placement aims at the
%                       dose centroid, not the isocenter).
%       config - Struct. Relevant fields (all optional, defaults in parens):
%           .sensor_side          - 'right' (default) | 'left'
%           .elements_per_side    - N for an NxN array (32)
%           .element_pitch_mm     - centre-to-centre pitch (7.00)
%           .element_size_mm      - active element width (2.43)
%           .sensor_standoff_mm   - body-to-panel gap (5)
%           .aim_at_centroid      - enable tilt toward centroid (true)
%           .force_turn_angle     - forced tilt magnitude, deg. 0 = feasibility-
%                                   limited search (default 0). Nonzero applies
%                                   exactly this angle and expands the grid to
%                                   fit. (NOTE: unlike determine_sensor_mask,
%                                   which defaults to 290 to work around its
%                                   anterior rotation bug, this file builds the
%                                   panel frame correctly and defaults to 0.)
%           .auto_select_side     - gas-based side auto-pick (false)
%           .avoid_bone           - bone-aware placement refinement (true)
%           .bone_hu_threshold    - HU at/above which a voxel is bone (150)
%           .gas_hu_threshold     - HU at/below which an in-body voxel is gas (-300)
%           .bone_penalty_mm      - mm added to a candidate's score per unit
%                                   bone fraction along its ray (30)
%           .bone_search_margin_mm- only candidates within this distance of the
%                                   best-by-distance placement are bone-scored (30)
%
%   OUTPUTS:
%       sensor_mask - 3D logical (same size as the possibly-expanded grid).
%       sensor_info - Diagnostic struct. Same fields as determine_sensor_mask
%                     plus: .sensor_side, .panel_normal (outward, flat),
%                     .aim_target_mm (the dose centroid), .side_auto_selected,
%                     .bone_frac_on_ray, .gas_frac_on_ray.
%
%   COORDINATE SYSTEM (HFS, identical to determine_sensor_mask):
%       array(Y, X, Z) = array(row, col, slice).
%       X (cols, dim 2): lower index = patient RIGHT, higher = patient LEFT.
%       Y (rows, dim 1): lower index = ANTERIOR, higher = POSTERIOR.
%       Z (slices, dim 3): higher index = SUPERIOR (head), lower = INFERIOR (feet).
%       World mm: x=origin(1)+(ix-1)*dx, y=origin(2)+(iy-1)*dy, z=origin(3)+(iz-1)*dz.
%
%   DEPENDENCIES: none beyond base MATLAB.
%   See also: determine_sensor_mask, determine_sensor_placement_fixed,
%             run_single_field_simulation

%% ======================== CONFIG DEFAULTS ========================

elements_per_side = get_field(config, 'elements_per_side', 32);
element_pitch_mm  = get_field(config, 'element_pitch_mm',  7.00);
element_size_mm   = get_field(config, 'element_size_mm',   2.43);
standoff_mm       = get_field(config, 'sensor_standoff_mm', 5);
aim_at_centroid   = get_field(config, 'aim_at_centroid',   true);
force_turn_deg    = get_field(config, 'force_turn_angle',  0);
sensor_side       = lower(get_field(config, 'sensor_side', 'right'));
auto_select_side  = get_field(config, 'auto_select_side',  false);
avoid_bone        = get_field(config, 'avoid_bone',        true);
bone_hu           = get_field(config, 'bone_hu_threshold', 150);
gas_hu            = get_field(config, 'gas_hu_threshold', -300);
bone_penalty_mm   = get_field(config, 'bone_penalty_mm',   30);
bone_margin_mm    = get_field(config, 'bone_search_margin_mm', 30);

if ~ismember(sensor_side, {'right', 'left'})
    warning('determine_sensor_mask_lateral:BadSide', ...
        'config.sensor_side "%s" not recognised; defaulting to right.', sensor_side);
    sensor_side = 'right';
end

kerf_mm = element_pitch_mm - element_size_mm;
if kerf_mm < 0
    error('determine_sensor_mask_lateral:InvalidGeometry', ...
        'element_size_mm (%.3f) cannot exceed element_pitch_mm (%.3f).', ...
        element_size_mm, element_pitch_mm);
end
aperture_mm = elements_per_side * element_pitch_mm;

fprintf('        [SensorLat] Placing %s-side array: %dx%d elements, pitch %.2f mm, size %.2f mm, kerf %.2f mm\n', ...
    sensor_side, elements_per_side, elements_per_side, element_pitch_mm, element_size_mm, kerf_mm);
fprintf('        [SensorLat] Aperture %.1f mm, standoff %.0f mm\n', ...
    aperture_mm, standoff_mm);

%% ======================== GRID INFO ========================

grid_dims = size(sct_resampled.bodyMask);
Ny = grid_dims(1);  % rows   = Y (anterior-posterior)
Nx = grid_dims(2);  % cols   = X (right-left)
Nz = grid_dims(3);  % slices = Z (superior-inferior)

dx = sct_resampled.spacing(1);
dy = sct_resampled.spacing(2);
dz = sct_resampled.spacing(3);

if isfield(sct_resampled, 'origin') && ~isempty(sct_resampled.origin)
    origin = sct_resampled.origin(:)';
else
    origin = [0, 0, 0];
    fprintf('        [SensorLat] Warning: origin missing; using [0,0,0].\n');
end

% Body mask excluding couch.
if isfield(sct_resampled, 'couchMask') && ~isempty(sct_resampled.couchMask)
    body = sct_resampled.bodyMask & ~sct_resampled.couchMask;
else
    body = logical(sct_resampled.bodyMask);
end

% Optional HU cube -> bone / gas masks for the secondary features.
have_hu = isfield(sct_resampled, 'cubeHU') && ~isempty(sct_resampled.cubeHU) && ...
          isequal(size(sct_resampled.cubeHU), grid_dims);
if have_hu
    hu = double(sct_resampled.cubeHU);
    bone_mask = (hu >= bone_hu) & body;
    gas_mask  = (hu <= gas_hu)  & body;
else
    bone_mask = false(grid_dims);
    gas_mask  = false(grid_dims);
    if avoid_bone || auto_select_side
        fprintf('        [SensorLat] No cubeHU available; disabling bone/gas features.\n');
    end
    avoid_bone       = false;
    auto_select_side = false;
end

fprintf('        [SensorLat] Grid: [%d(Y) x %d(X) x %d(Z)], spacing [%.2f, %.2f, %.2f] mm\n', ...
    Ny, Nx, Nz, dx, dy, dz);

% Isocenter is logged for reference only (aim is the dose centroid, not iso).
if ~isempty(beam_metadata) && isfield(beam_metadata(1), 'isocenter') && ...
        numel(beam_metadata(1).isocenter) == 3
    fprintf('        [SensorLat] Isocenter (reference only): [%.1f, %.1f, %.1f] mm\n', ...
        beam_metadata(1).isocenter);
end

%% ======================== PLACEMENT DOSE (SUMMED PLAN DOSE) ========================
% Placement (exclusion zone, centroid, and tilt target) must be identical for
% every field/segment of a plan, so it uses the SUMMED plan dose when supplied.
placement_dose = resolve_placement_dose(field_dose, config, [Ny, Nx, Nz]);

dose_max = max(placement_dose(:));
if isempty(placement_dose) || dose_max <= 0
    warning('determine_sensor_mask_lateral:NoDose', ...
        'Placement dose empty or zero; returning empty sensor mask.');
    [sensor_mask, sensor_info] = empty_result(grid_dims, sensor_side);
    return;
end

% Dose centroid (aim target), in mm and voxels.
centroid_src = field_dose;
centroid_src.dose_Gy = placement_dose;
centroid_mm = compute_dose_centroid_mm(centroid_src, origin, dx, dy, dz);
centroid_iy = clampi(round((centroid_mm(2) - origin(2)) / dy) + 1, 1, Ny);
centroid_ix = clampi(round((centroid_mm(1) - origin(1)) / dx) + 1, 1, Nx);
centroid_iz = clampi(round((centroid_mm(3) - origin(3)) / dz) + 1, 1, Nz);
fprintf('        [SensorLat] Dose centroid: [%.1f, %.1f, %.1f] mm (voxel Y=%d, X=%d, Z=%d)\n', ...
    centroid_mm, centroid_iy, centroid_ix, centroid_iz);

%% ======================== BEAM/DOSE EXCLUSION (Y,Z SHADOW) ========================
% A lateral sensor "in the beam zone" is one sitting over the irradiated volume
% as seen from the side. Mark (Y,Z) columns whose summed dose reaches >=10% of
% peak anywhere along X; the sensor footprint must avoid them.
dose_thresh   = 0.10 * dose_max;
col_max_yz    = squeeze(max(placement_dose, [], 2));   % max along X -> (Y,Z)
exclusion_zone = col_max_yz >= dose_thresh;            % (Ny x Nz) logical
fprintf('        [SensorLat] Beam/dose exclusion: %d/%d (Y,Z) columns >= 10%% peak (%.3f Gy)\n', ...
    sum(exclusion_zone(:)), Ny * Nz, dose_thresh);

%% ======================== OPTIONAL SIDE AUTO-SELECT (gas) ========================
side_auto_selected = false;
if auto_select_side
    gas_right = gas_frac_lateral('right', centroid_iy, centroid_ix, centroid_iz, body, gas_mask, Nx);
    gas_left  = gas_frac_lateral('left',  centroid_iy, centroid_ix, centroid_iz, body, gas_mask, Nx);
    fprintf('        [SensorLat] Gas along centroid->surface ray: right=%.3f, left=%.3f\n', ...
        gas_right, gas_left);
    if gas_left < gas_right
        sensor_side = 'left';
    else
        sensor_side = 'right';
    end
    side_auto_selected = true;
    fprintf('        [SensorLat] Auto-selected side: %s (less gas)\n', sensor_side);
end

is_right = strcmp(sensor_side, 'right');

% Gas along the straight centroid->side-surface ray (diagnostic). Computed on
% the ORIGINAL grid before any expansion so gas_mask/body stay the same size.
if have_hu
    gas_frac_on_ray = gas_frac_lateral(sensor_side, centroid_iy, centroid_ix, centroid_iz, body, gas_mask, Nx);
else
    gas_frac_on_ray = NaN;
end

%% ======================== LATERAL SURFACE MAP ========================
% surf(iy,iz) = the skin X index on the chosen side for column (iy,iz):
%   right -> min X with body (most rightward, lowest index)
%   left  -> max X with body (most leftward,  highest index)
lateral_surface = NaN(Ny, Nz);
for iy = 1:Ny
    for iz = 1:Nz
        xcol = find(body(iy, :, iz));
        if ~isempty(xcol)
            if is_right
                lateral_surface(iy, iz) = xcol(1);     % min X
            else
                lateral_surface(iy, iz) = xcol(end);   % max X
            end
        end
    end
end
surface_valid = ~isnan(lateral_surface);
if ~any(surface_valid(:))
    warning('determine_sensor_mask_lateral:NoSurface', ...
        'No lateral body surface found; returning empty sensor mask.');
    [sensor_mask, sensor_info] = empty_result(grid_dims, sensor_side);
    sensor_info.exclusion_zone = exclusion_zone;
    return;
end

%% ======================== FLAT PLACEMENT SWEEP (Y,Z footprint) ========================
% Footprint size in voxels (aperture is square in mm).
sensor_ny = ceil(aperture_mm / dy);
sensor_nz = ceil(aperture_mm / dz);

% Shrink the element grid if the footprint cannot fit the in-plane grid extent.
fit_ny = max(1, floor((Ny * dy) / element_pitch_mm));
fit_nz = max(1, floor((Nz * dz) / element_pitch_mm));
N_fit  = min([elements_per_side, fit_ny, fit_nz]);
if N_fit < elements_per_side
    warning('determine_sensor_mask_lateral:ApertureShrunk', ...
        'Grid fits %dx%d elements at pitch %.2f mm; reducing from %dx%d.', ...
        N_fit, N_fit, element_pitch_mm, elements_per_side, elements_per_side);
    elements_per_side = N_fit;
    aperture_mm = elements_per_side * element_pitch_mm;
    sensor_ny = ceil(aperture_mm / dy);
    sensor_nz = ceil(aperture_mm / dz);
end
N_total_elements = elements_per_side * elements_per_side;

half_ny = floor(sensor_ny / 2);
half_nz = floor(sensor_nz / 2);

% Candidate = footprint whose surface is fully valid and fully outside the
% exclusion shadow. Objective: minimise (Y,Z) distance to the centroid; when
% avoid_bone is on, break near-ties toward the least-bone acoustic path.
best_iy_start = []; best_iz_start = []; best_dist = Inf;
cand_iy = zeros(0, 1); cand_iz = zeros(0, 1); cand_dist = zeros(0, 1);
collect_cands = avoid_bone;

for iy_start = 1:(Ny - sensor_ny + 1)
    iy_end = iy_start + sensor_ny - 1;
    for iz_start = 1:(Nz - sensor_nz + 1)
        iz_end = iz_start + sensor_nz - 1;
        surf_patch = surface_valid(iy_start:iy_end, iz_start:iz_end);
        if ~all(surf_patch(:))
            continue;
        end
        excl_patch = exclusion_zone(iy_start:iy_end, iz_start:iz_end);
        if any(excl_patch(:))
            continue;   % footprint overlaps the beam/dose shadow
        end
        cy_mm = origin(2) + (iy_start + half_ny - 1) * dy;
        cz_mm = origin(3) + (iz_start + half_nz - 1) * dz;
        dist  = hypot(cy_mm - centroid_mm(2), cz_mm - centroid_mm(3));
        if dist < best_dist
            best_dist = dist; best_iy_start = iy_start; best_iz_start = iz_start;
        end
        if collect_cands
            cand_iy(end+1, 1) = iy_start; %#ok<AGROW>
            cand_iz(end+1, 1) = iz_start; %#ok<AGROW>
            cand_dist(end+1, 1) = dist;   %#ok<AGROW>
        end
    end
end

if isempty(best_iy_start)
    % No footprint clears the exclusion shadow inside the grid. Fall back to the
    % centroid's (Y,Z) projection (clamped) and let grid expansion fit it.
    fprintf('        [SensorLat] No exclusion-free footprint in grid; centering on centroid projection.\n');
    best_iy_start = clampi(centroid_iy - half_ny, 1, max(1, Ny - sensor_ny + 1));
    best_iz_start = clampi(centroid_iz - half_nz, 1, max(1, Nz - sensor_nz + 1));
end

% Secondary: bone-aware refinement among placements close to the best distance.
bone_frac_on_ray = NaN;
if avoid_bone && ~isempty(cand_dist)
    keep = cand_dist <= (best_dist + bone_margin_mm);
    ci = cand_iy(keep); cz = cand_iz(keep); cd = cand_dist(keep);
    best_score = Inf; ref_iy = best_iy_start; ref_iz = best_iz_start; ref_bone = NaN;
    for k = 1:numel(cd)
        iy0 = ci(k); iz0 = cz(k);
        cyc = iy0 + half_ny; czc = iz0 + half_nz;                 % footprint centre (Y,Z voxel)
        sxc = round(lateral_surface(clampi(cyc,1,Ny), clampi(czc,1,Nz)));  % skin X there
        if isnan(sxc), continue; end
        bfrac = bone_frac_on_line(bone_mask, [cyc, sxc, czc], ...
                                  [centroid_iy, centroid_ix, centroid_iz]);
        score = cd(k) + bone_penalty_mm * bfrac;
        if score < best_score
            best_score = score; ref_iy = iy0; ref_iz = iz0; ref_bone = bfrac;
        end
    end
    best_iy_start = ref_iy; best_iz_start = ref_iz; bone_frac_on_ray = ref_bone;
    fprintf('        [SensorLat] Bone-aware pick: footprint Y=%d, Z=%d, bone frac on ray=%.3f\n', ...
        best_iy_start, best_iz_start, bone_frac_on_ray);
end

sensor_y_range = [best_iy_start, best_iy_start + sensor_ny - 1];
sensor_z_range = [best_iz_start, best_iz_start + sensor_nz - 1];

%% ======================== SENSOR X SLAB & CENTER ========================
% Flat slab X: most-lateral skin over the footprint, offset by the standoff so
% the whole flat panel clears the body. Right = lower X, left = higher X.
surf_fp = lateral_surface(sensor_y_range(1):sensor_y_range(2), ...
                          sensor_z_range(1):sensor_z_range(2));
surf_fp = surf_fp(~isnan(surf_fp));
if isempty(surf_fp)
    skin_x = round(lateral_surface(clampi(centroid_iy,1,Ny), clampi(centroid_iz,1,Nz)));
else
    if is_right, skin_x = min(surf_fp); else, skin_x = max(surf_fp); end
end
standoff_vox_x = ceil(standoff_mm / dx);
if is_right
    ix_sensor = round(skin_x) - standoff_vox_x;   % further right = lower X
    panel_normal = [-1, 0, 0];                    % outward = patient right
else
    ix_sensor = round(skin_x) + standoff_vox_x;   % further left = higher X
    panel_normal = [ 1, 0, 0];                    % outward = patient left
end

% Sensor centre in mm (fixed physical point; survives grid expansion).
c_mm = [origin(1) + (ix_sensor - 1) * dx, ...
        origin(2) + (mean(sensor_y_range) - 1) * dy, ...
        origin(3) + (mean(sensor_z_range) - 1) * dz];

fprintf('        [SensorLat] Flat slab: X index %d (skin X=%d, standoff %d vox), centre [%.1f, %.1f, %.1f] mm\n', ...
    ix_sensor, round(skin_x), standoff_vox_x, c_mm);

%% ======================== TILT TOWARD CENTROID ========================
% Rotate the flat panel so its outward normal points along (centre - centroid).
% Magnitude is the largest that keeps every element outside the body, unless a
% forced angle is requested.
tilt_R = eye(3); tilt_theta = 0; tilt_alpha = 0; tilt_axis = [0; 0; 1];
aim_normal_used = panel_normal;

if aim_at_centroid
    [~, tilt_alpha, tilt_axis_cand] = aim_rotation(panel_normal, c_mm, centroid_mm);

    if force_turn_deg ~= 0
        tilt_axis  = tilt_axis_cand;
        tilt_alpha = force_turn_deg * pi / 180;
        tilt_theta = 1;
        tilt_R     = rodrigues(tilt_axis, tilt_alpha);
        fprintf('        [SensorLat] FORCED tilt: %.2f deg (feasibility checks bypassed)\n', force_turn_deg);
    elseif abs(tilt_alpha) < 1e-6
        fprintf('        [SensorLat] Panel already aimed at centroid; no tilt needed.\n');
    else
        theta_max = binary_search_theta(c_mm, tilt_axis_cand, tilt_alpha, ...
            elements_per_side, element_pitch_mm, ...
            origin, dx, dy, dz, body, Nx, Ny, Nz);
        if theta_max > 0
            tilt_theta = theta_max;
            tilt_axis  = tilt_axis_cand;
            tilt_R     = rodrigues(tilt_axis, tilt_theta * tilt_alpha);
            fprintf('        [SensorLat] Tilt: theta=%.3f (%.2f of full aim %.2f deg)\n', ...
                tilt_theta, tilt_theta * tilt_alpha * 180/pi, tilt_alpha * 180/pi);
        else
            fprintf('        [SensorLat] Max feasible tilt = 0; keeping flat placement.\n');
        end
    end
    d = c_mm - centroid_mm; nd = norm(d);
    if nd > eps, aim_normal_used = d / nd; end
end

%% ======================== GRID EXPANSION TO FIT SENSOR ========================
% Element centres in mm (flat when tilt_R = I). If any element (allowing for its
% active footprint) lands past a grid edge in X/Y/Z, water-pad the grid to fit.
P_mm = element_centers_mm(c_mm, tilt_R, elements_per_side, element_pitch_mm);
vx = round((P_mm(:,1) - origin(1)) / dx) + 1;
vy = round((P_mm(:,2) - origin(2)) / dy) + 1;
vz = round((P_mm(:,3) - origin(3)) / dz) + 1;

margin_x = 1 + (tilt_theta == 0) * ceil((element_size_mm/2) / dx);
margin_y = 1 + (tilt_theta == 0) * ceil((element_size_mm/2) / dy);
margin_z = 1 + (tilt_theta == 0) * ceil((element_size_mm/2) / dz);

pad_x_pre  = max(0, 1  - (min(vx) - margin_x));
pad_x_post = max(0, (max(vx) + margin_x) - Nx);
pad_y_pre  = max(0, 1  - (min(vy) - margin_y));
pad_y_post = max(0, (max(vy) + margin_y) - Ny);
pad_z_pre  = max(0, 1  - (min(vz) - margin_z));
pad_z_post = max(0, (max(vz) + margin_z) - Nz);

grid_was_expanded = (pad_x_pre + pad_x_post + pad_y_pre + pad_y_post + pad_z_pre + pad_z_post) > 0;

if grid_was_expanded
    fprintf('        [SensorLat] Expanding grid by [X -%d/+%d, Y -%d/+%d, Z -%d/+%d] voxels (water)\n', ...
        pad_x_pre, pad_x_post, pad_y_pre, pad_y_post, pad_z_pre, pad_z_post);

    new_Ny = Ny + pad_y_pre + pad_y_post;
    new_Nx = Nx + pad_x_pre + pad_x_post;
    new_Nz = Nz + pad_z_pre + pad_z_post;

    body_new = false(new_Ny, new_Nx, new_Nz);
    body_new(pad_y_pre + (1:Ny), pad_x_pre + (1:Nx), pad_z_pre + (1:Nz)) = body;
    body = body_new;

    Ny = new_Ny; Nx = new_Nx; Nz = new_Nz;
    grid_dims = [Ny, Nx, Nz];

    % Origin = mm of voxel (1,1,1); pre-padding moves it back. mm coords fixed.
    origin(1) = origin(1) - pad_x_pre * dx;
    origin(2) = origin(2) - pad_y_pre * dy;
    origin(3) = origin(3) - pad_z_pre * dz;

    sensor_y_range = sensor_y_range + pad_y_pre;
    sensor_z_range = sensor_z_range + pad_z_pre;
    ix_sensor      = ix_sensor + pad_x_pre;
end

%% ======================== VOXELISE SENSOR ========================
% Recompute element centres against the (possibly shifted) origin, then voxelise:
%   tilted -> one voxel per element centre (sparse array of discrete elements)
%   flat   -> multi-voxel active footprint per element in the constant-X slab
sensor_mask          = false(grid_dims);
sensor_element_label = zeros(grid_dims, 'uint32');
element_positions_mm = element_centers_mm(c_mm, tilt_R, elements_per_side, element_pitch_mm);

if tilt_theta > 0
    for e = 1:N_total_elements
        pw = element_positions_mm(e, :);
        ex = round((pw(1) - origin(1)) / dx) + 1;
        ey = round((pw(2) - origin(2)) / dy) + 1;
        ez = round((pw(3) - origin(3)) / dz) + 1;
        if ex < 1 || ex > Nx || ey < 1 || ey > Ny || ez < 1 || ez > Nz, continue; end
        if sensor_element_label(ey, ex, ez) == 0
            sensor_mask(ey, ex, ez) = true;
            sensor_element_label(ey, ex, ez) = uint32(e);
        end
    end
else
    half_ay = (element_size_mm / 2) / dy;   % active half-width, Y voxels
    half_az = (element_size_mm / 2) / dz;   % active half-width, Z voxels
    for e = 1:N_total_elements
        pw = element_positions_mm(e, :);
        eyc = (pw(2) - origin(2)) / dy + 1;
        ezc = (pw(3) - origin(3)) / dz + 1;
        iy_lo = max(1,  ceil(eyc - half_ay));  iy_hi = min(Ny, floor(eyc + half_ay));
        iz_lo = max(1,  ceil(ezc - half_az));  iz_hi = min(Nz, floor(ezc + half_az));
        if iy_hi < iy_lo, iy_lo = clampi(round(eyc),1,Ny); iy_hi = iy_lo; end
        if iz_hi < iz_lo, iz_lo = clampi(round(ezc),1,Nz); iz_hi = iz_lo; end
        for iy = iy_lo:iy_hi
            for iz = iz_lo:iz_hi
                if sensor_element_label(iy, ix_sensor, iz) == 0
                    sensor_mask(iy, ix_sensor, iz) = true;
                    sensor_element_label(iy, ix_sensor, iz) = uint32(e);
                end
            end
        end
    end
end

%% ======================== VALIDATION ========================
placement_valid = true;
overlap_body = sensor_mask & body;
if any(overlap_body(:))
    n_over = sum(overlap_body(:));
    warning('determine_sensor_mask_lateral:BodyOverlap', ...
        'Sensor overlaps body at %d voxels. Removing overlapping voxels.', n_over);
    sensor_mask = sensor_mask & ~body;
    sensor_element_label(overlap_body) = 0;
    placement_valid = false;
end
num_sensor_voxels = sum(sensor_mask(:));
if num_sensor_voxels == 0
    warning('determine_sensor_mask_lateral:EmptySensor', ...
        'Sensor mask empty after placement (side %s). Check body mask / dose.', sensor_side);
end
fprintf('        [SensorLat] Final sensor: %d voxels (valid: %s)\n', ...
    num_sensor_voxels, mat2str(placement_valid));

%% ======================== ELEMENT DIAGNOSTICS ========================
sensor_lin        = find(sensor_mask);
voxel_element_idx = double(sensor_element_label(sensor_lin));

element_voxel_counts = zeros(N_total_elements, 1);
for v = 1:numel(voxel_element_idx)
    e = voxel_element_idx(v);
    if e > 0, element_voxel_counts(e) = element_voxel_counts(e) + 1; end
end
effective_element_count = sum(element_voxel_counts > 0);
aliased_element_count   = N_total_elements - effective_element_count;
if effective_element_count > 0
    vpe = element_voxel_counts(element_voxel_counts > 0);
    voxels_per_element_stats = [min(vpe), median(vpe), max(vpe)];
else
    voxels_per_element_stats = [0, 0, 0];
end

if ~isempty(sensor_lin)
    [svy, svx, svz] = ind2sub(grid_dims, sensor_lin);
    sensor_x_range = [min(svx), max(svx)];
    sensor_z_range = [min(svz), max(svz)];
    sensor_y_range = [min(svy), max(svy)];
    sensor_y_index = round(mean(svy));
else
    sensor_x_range = [0, 0]; sensor_y_range = [0, 0]; sensor_z_range = [0, 0];
    sensor_y_index = 0;
end

aperture_footprint_voxels = sensor_ny * sensor_nz;
if aperture_footprint_voxels > 0
    fill_factor_actual = num_sensor_voxels / aperture_footprint_voxels;
else
    fill_factor_actual = 0;
end

fprintf('        [SensorLat] Elements: %d/%d effective, %d aliased/removed; fill %.1f%%\n', ...
    effective_element_count, N_total_elements, aliased_element_count, 100 * fill_factor_actual);
if aliased_element_count > 0
    warning('determine_sensor_mask_lateral:ElementAliasing', ...
        '%d of %d elements collapsed (grid coarser than pitch or voxels removed).', ...
        aliased_element_count, N_total_elements);
end

%% ======================== PACK OUTPUT ========================
sensor_info = struct();
sensor_info.sensor_side                 = sensor_side;
sensor_info.side_auto_selected          = side_auto_selected;
sensor_info.panel_normal                = panel_normal;
sensor_info.sensor_x_index              = ix_sensor;
sensor_info.sensor_y_index              = sensor_y_index;
sensor_info.sensor_x_range              = sensor_x_range;
sensor_info.sensor_y_range              = sensor_y_range;
sensor_info.sensor_z_range              = sensor_z_range;
sensor_info.sensor_center_mm            = c_mm;
sensor_info.exclusion_zone              = exclusion_zone;
sensor_info.voxel_element_idx           = voxel_element_idx;
sensor_info.num_elements                = effective_element_count;
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
sensor_info.surface_to_dose_distance_mm = norm(c_mm - centroid_mm);
sensor_info.sensor_size_voxels          = [sensor_ny, sensor_nz];
sensor_info.placement_valid             = placement_valid;
sensor_info.num_sensor_voxels           = num_sensor_voxels;
sensor_info.standoff_voxels             = standoff_vox_x;
sensor_info.gantry_angle                = get_field(field_dose, 'gantry_angle', NaN);

% Tilt geometry
sensor_info.tilt_theta                  = tilt_theta;
sensor_info.tilt_angle_deg              = tilt_theta * tilt_alpha * 180 / pi;
sensor_info.tilt_alpha_full_deg         = tilt_alpha * 180 / pi;
sensor_info.rotation_R                  = tilt_R;
sensor_info.rotation_axis               = tilt_axis(:)';
sensor_info.aim_target_mm               = centroid_mm;   % dose centroid, not iso
sensor_info.aim_normal                  = aim_normal_used;
sensor_info.aim_enabled                 = aim_at_centroid && (tilt_theta > 0);

% Secondary-feature diagnostics
sensor_info.bone_frac_on_ray            = bone_frac_on_ray;
sensor_info.gas_frac_on_ray             = gas_frac_on_ray;

% Grid-expansion contract (same field layout the caller consumes).
sensor_info.grid_pad = struct( ...
    'y_pre',  pad_y_pre,  'y_post', pad_y_post, ...
    'x_pre',  pad_x_pre,  'x_post', pad_x_post, ...
    'z_pre',  pad_z_pre,  'z_post', pad_z_post, ...
    'expanded', grid_was_expanded);

fprintf('        [SensorLat] Placement complete (%s side).\n', sensor_side);

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function P = element_centers_mm(c_mm, R, N, pitch)
%ELEMENT_CENTERS_MM  World-mm centres of an NxN panel. In-plane axes are Y and
%   Z (offsets orthogonal to the panel normal); R rigidly rotates the panel.
%   Element order e = (bi-1)*N + ai, bi along Y, ai along Z.
    P = zeros(N * N, 3);
    half = (N - 1) / 2;
    e = 0;
    for bi = 1:N
        b_off = (bi - 1 - half) * pitch;    % along Y
        for ai = 1:N
            a_off = (ai - 1 - half) * pitch; % along Z
            e = e + 1;
            off = [0; b_off; a_off];
            P(e, :) = c_mm + (R * off)';
        end
    end
end


function feasible = check_placement(c_mm, R, N, pitch, ...
        origin, dx, dy, dz, body, Nx, Ny, Nz)
%CHECK_PLACEMENT  True iff no element centre at (c_mm, R) lands inside the body.
%   Out-of-grid elements are treated as air/water (feasible); the caller expands
%   the grid to fit them.
    P = element_centers_mm(c_mm, R, N, pitch);
    feasible = true;
    for e = 1:size(P, 1)
        ex = round((P(e,1) - origin(1)) / dx) + 1;
        ey = round((P(e,2) - origin(2)) / dy) + 1;
        ez = round((P(e,3) - origin(3)) / dz) + 1;
        if ex >= 1 && ex <= Nx && ey >= 1 && ey <= Ny && ez >= 1 && ez <= Nz
            if body(ey, ex, ez)
                feasible = false;
                return;
            end
        end
    end
end


function theta_max = binary_search_theta(c_mm, axis_k, alpha, N, pitch, ...
        origin, dx, dy, dz, body, Nx, Ny, Nz)
%BINARY_SEARCH_THETA  Largest theta in [0,1] keeping all elements outside body
%   at rotation rodrigues(axis_k, theta*alpha). Returns 0 if even flat collides.
    if ~check_placement(c_mm, eye(3), N, pitch, origin, dx, dy, dz, body, Nx, Ny, Nz)
        theta_max = 0; return;
    end
    if abs(alpha) < 1e-6, theta_max = 0; return; end
    if check_placement(c_mm, rodrigues(axis_k, alpha), N, pitch, origin, dx, dy, dz, body, Nx, Ny, Nz)
        theta_max = 1; return;
    end
    lo = 0; hi = 1;
    for it = 1:12 %#ok<NASGU>
        mid = 0.5 * (lo + hi);
        if check_placement(c_mm, rodrigues(axis_k, mid*alpha), N, pitch, origin, dx, dy, dz, body, Nx, Ny, Nz)
            lo = mid;
        else
            hi = mid;
        end
    end
    theta_max = lo;
end


function [R, alpha, axis_k] = aim_rotation(n0, c_mm, target_mm)
%AIM_ROTATION  Shortest-arc rotation taking the flat outward normal n0 to the
%   direction (c_mm - target_mm)/|.|, so the tilted panel faces the target.
    n0 = n0(:)' / norm(n0);
    d  = c_mm(:)' - target_mm(:)';
    nd = norm(d);
    if nd < eps
        R = eye(3); alpha = 0; axis_k = [0; 0; 1]; return;
    end
    n_aim = d / nd;
    cross_v = cross(n0, n_aim);
    sin_a = norm(cross_v);
    cos_a = dot(n0, n_aim);
    alpha = atan2(sin_a, cos_a);
    if sin_a < 1e-9
        if cos_a > 0
            R = eye(3); axis_k = [0; 0; 1]; alpha = 0;
        else
            axis_k = [0; 0; 1];              % 180 deg flip about Z (n0 is +/-X)
            R = rodrigues(axis_k, alpha);
        end
        return;
    end
    axis_k = cross_v(:) / sin_a;
    R = rodrigues(axis_k, alpha);
end


function R = rodrigues(k, ang)
%RODRIGUES  Rotation matrix from unit axis k and angle ang (rad).
    k = k(:); nk = norm(k);
    if nk < eps || abs(ang) < eps, R = eye(3); return; end
    k = k / nk;
    K = [   0   -k(3)  k(2);
          k(3)    0   -k(1);
         -k(2)  k(1)    0  ];
    R = eye(3) + sin(ang) * K + (1 - cos(ang)) * (K * K);
end


function frac = bone_frac_on_line(bone_mask, p0_vox, p1_vox)
%BONE_FRAC_ON_LINE  Fraction of samples between two voxels (inclusive) that are
%   bone. p0_vox, p1_vox are [iy, ix, iz]. Simple equal-step line sampling.
    d = p1_vox - p0_vox;
    nsteps = max(1, round(max(abs(d))));
    hits = 0;
    [Ny, Nx, Nz] = size(bone_mask);
    for s = 0:nsteps
        t = s / nsteps;
        iy = clampi(round(p0_vox(1) + t*d(1)), 1, Ny);
        ix = clampi(round(p0_vox(2) + t*d(2)), 1, Nx);
        iz = clampi(round(p0_vox(3) + t*d(3)), 1, Nz);
        if bone_mask(iy, ix, iz), hits = hits + 1; end
    end
    frac = hits / (nsteps + 1);
end


function frac = gas_frac_lateral(side, iy, ix_centroid, iz, body, gas_mask, Nx)
%GAS_FRAC_LATERAL  Fraction of the straight lateral ray from the centroid to the
%   chosen-side skin (constant Y,Z) that is intra-body gas. Used to pick a side.
    if strcmp(side, 'right')
        xs = ix_centroid:-1:1;      % march toward lower X (patient right)
    else
        xs = ix_centroid:1:Nx;      % march toward higher X (patient left)
    end
    in_body_seen = false; total = 0; hits = 0;
    for ix = xs
        if body(iy, ix, iz)
            in_body_seen = true;
            total = total + 1;
            if gas_mask(iy, ix, iz), hits = hits + 1; end
        elseif in_body_seen
            break;                  % left the body -> reached the skin
        end
    end
    if total == 0, frac = 0; else, frac = hits / total; end
end


function v = clampi(v, lo, hi)
%CLAMPI  Clamp scalar/array to [lo, hi].
    v = max(lo, min(hi, v));
end


function dose = resolve_placement_dose(field_dose, config, expected_dims)
%RESOLVE_PLACEMENT_DOSE  Dose used for placement (exclusion + centroid + aim).
%   Prefers the SUMMED plan dose so placement is identical across all fields:
%     1) field_dose.total_dose_Gy  2) config.placement_dose_Gy
%     3) config.total_dose_file    4) field_dose.dose_Gy (per-field fallback)
    dose = [];
    if isfield(field_dose, 'total_dose_Gy') && ~isempty(field_dose.total_dose_Gy)
        dose = accept_if_matches(field_dose.total_dose_Gy, expected_dims, 'field_dose.total_dose_Gy');
    end
    if isempty(dose) && isfield(config, 'placement_dose_Gy') && ~isempty(config.placement_dose_Gy)
        dose = accept_if_matches(config.placement_dose_Gy, expected_dims, 'config.placement_dose_Gy');
    end
    if isempty(dose) && isfield(config, 'total_dose_file') && ...
            ~isempty(config.total_dose_file) && isfile(config.total_dose_file)
        dose = accept_if_matches(load_total_dose_file(config.total_dose_file), ...
            expected_dims, config.total_dose_file);
    end
    if isempty(dose)
        dose = field_dose.dose_Gy;
        fprintf('        [SensorLat] Placement dose: per-field dose (no summed plan dose supplied).\n');
    else
        fprintf('        [SensorLat] Placement dose: SUMMED PLAN DOSE (deterministic across fields).\n');
    end
end


function out = accept_if_matches(arr, expected_dims, label)
%ACCEPT_IF_MATCHES  Return arr (double) if its size matches; else [] + warning.
    out = [];
    if isempty(arr), return; end
    if isequal(size(arr), expected_dims)
        out = double(arr);
    else
        warning('determine_sensor_mask_lateral:PlacementDoseSizeMismatch', ...
            '%s size %s ~= grid %s; ignoring for placement.', ...
            label, mat2str(size(arr)), mat2str(expected_dims));
    end
end


function dose = load_total_dose_file(fpath)
%LOAD_TOTAL_DOSE_FILE  Load a step15 total dose (full or sparse) into a 3D array.
    dose = [];
    S = load(fpath);
    if isfield(S, 'total_rs_dose')
        dose = S.total_rs_dose;
    elseif isfield(S, 'total_rs_dose_sparse') && isfield(S, 'total_rs_dose_dims')
        dose = reshape(full(S.total_rs_dose_sparse), S.total_rs_dose_dims(:)');
    elseif isfield(S, 'ct_total')
        dose = S.ct_total;
    elseif isfield(S, 'ct_total_sparse') && isfield(S, 'ct_total_dims')
        dose = reshape(full(S.ct_total_sparse), S.ct_total_dims(:)');
    else
        warning('determine_sensor_mask_lateral:UnknownTotalDoseFile', ...
            'No recognised total-dose variable in %s.', fpath);
    end
    dose = double(dose);
end


function centroid_mm = compute_dose_centroid_mm(field_dose, origin, dx, dy, dz)
%COMPUTE_DOSE_CENTROID_MM  Physical centroid of the dose distribution (mm).
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


function val = get_field(s, name, default)
%GET_FIELD  Struct field with fallback default.
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default;
    end
end


function [sensor_mask, sensor_info] = empty_result(grid_dims, sensor_side)
%EMPTY_RESULT  Empty mask + default info struct (grid_pad contract preserved).
    sensor_mask = false(grid_dims);
    sensor_info = struct();
    sensor_info.sensor_side                 = sensor_side;
    sensor_info.side_auto_selected          = false;
    sensor_info.panel_normal                = [0, 0, 0];
    sensor_info.sensor_x_index              = 0;
    sensor_info.sensor_y_index              = 0;
    sensor_info.sensor_x_range              = [0, 0];
    sensor_info.sensor_y_range              = [0, 0];
    sensor_info.sensor_z_range              = [0, 0];
    sensor_info.sensor_center_mm            = [0, 0, 0];
    sensor_info.exclusion_zone              = false(grid_dims(1), grid_dims(3));
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
    sensor_info.standoff_voxels             = 0;
    sensor_info.gantry_angle                = NaN;
    sensor_info.tilt_theta                  = 0;
    sensor_info.tilt_angle_deg              = 0;
    sensor_info.tilt_alpha_full_deg         = 0;
    sensor_info.rotation_R                  = eye(3);
    sensor_info.rotation_axis               = [0, 0, 1];
    sensor_info.aim_target_mm               = [];
    sensor_info.aim_normal                  = [0, 0, 0];
    sensor_info.aim_enabled                 = false;
    sensor_info.bone_frac_on_ray            = NaN;
    sensor_info.gas_frac_on_ray             = NaN;
    sensor_info.grid_pad = struct( ...
        'y_pre', 0, 'y_post', 0, 'x_pre', 0, 'x_post', 0, ...
        'z_pre', 0, 'z_post', 0, 'expanded', false);
end
