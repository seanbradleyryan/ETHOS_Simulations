function placement = place_sensor(patient_id, session, CONFIG)
%PLACE_SENSOR  Interactive sensor placement UI for ETHOS IRAI pipeline.
%
%   placement = place_sensor(patient_id, session, CONFIG)
%
%   PURPOSE:
%   Render a composite CT + summed-dose + beam-projection view, accept an
%   interactive sensor centre via drawrectangle, build a tilted ultrasound
%   sensor footprint snapped to the anterior body surface, validate it, and
%   persist `sensor_placement.mat` for downstream consumption by
%   `run_single_field_simulation.m`. Intended to be invoked ONCE per
%   patient/session when deciding sensor geometry.
%
%   INPUTS:
%       patient_id - char/string. Patient identifier (e.g., '1194203').
%       session    - char/string. Session label (e.g., 'Session_1').
%       CONFIG     - Struct of placement parameters:
%           .working_dir         - Base directory containing RayStationFiles
%           .sensor_size_cm      - [Lx, Lz] sensor physical size in cm
%                                  (default [10, 10])
%           .sensor_standoff_mm  - Body-to-sensor gap in mm (default 5)
%           .tilt_angle_deg      - Tilt about X-axis toward isocenter, deg
%                                  (default 15)
%           .pml_size            - PML thickness in voxels (default 10)
%
%   OUTPUT:
%       placement - Struct with fields:
%           .sensor_mask         - [Ny x Nx x Nz] logical sensor mask
%           .sensor_center_mm    - [x, y, z] mm centre in DICOM patient coords
%           .sensor_normal_hat   - 3x1 unit normal (DICOM patient coords)
%           .tilt_angle_deg      - Signed tilt angle used (deg)
%           .placement_valid     - true if all validation checks passed
%           .patient_id          - char
%           .session             - char
%           .timestamp           - datestr ISO 8601
%           .config              - Echo of CONFIG used
%
%   ALGORITHM:
%       1. Validate inputs and load sct_resampled + all field_dose_*.mat
%          (accumulating summed dose one file at a time).
%       2. Compute dose-weighted centroid slice index isoIz.
%       3. Draw composite figure (CT + dose overlay + body contour + beam
%          projection lines through isocenter) at slice isoIz.
%       4. Accept user sensor centre via drawrectangle (mm).
%       5. Build flat axial sensor mask at standoff from anterior surface.
%       6. Tilt sensor about X-axis toward isocenter.
%       7. Snap each tilted (X,Z) column anteriorly to the body surface.
%       8. Validate: warn on body overlap or PML proximity.
%       9. Visualise final tilted sensor (2D overlay + 3D scatter).
%      10. Save sensor_placement.mat in processed/ directory.
%
%   CONSTRAINTS:
%       - Uses figure + drawrectangle only (no App Designer / uifigure).
%       - Interactive; not parfor-safe by design.
%       - Issues warning(), not error(), on non-fatal failures.
%
%   EXAMPLE:
%       CONFIG.working_dir        = '/mnt/weka/home/80030361/ETHOS_Simulations';
%       CONFIG.sensor_size_cm     = [10, 10];
%       CONFIG.sensor_standoff_mm = 5;
%       CONFIG.tilt_angle_deg     = 15;
%       CONFIG.pml_size           = 10;
%       placement = place_sensor('1194203', 'Session_1', CONFIG);
%
%   DEPENDENCIES:
%       Image Processing Toolbox (drawrectangle, bwboundaries).
%
%   See also: determine_sensor_placement_fixed, run_single_field_simulation,
%             step15_process_doses

%% ======================== INPUT VALIDATION ========================

if ~ischar(patient_id) && ~isstring(patient_id)
    error('place_sensor:InvalidInput', ...
        'patient_id must be a char or string.');
end
if ~ischar(session) && ~isstring(session)
    error('place_sensor:InvalidInput', ...
        'session must be a char or string.');
end
if ~isstruct(CONFIG)
    error('place_sensor:InvalidInput', 'CONFIG must be a struct.');
end
if ~isfield(CONFIG, 'working_dir') || isempty(CONFIG.working_dir)
    error('place_sensor:InvalidInput', ...
        'CONFIG.working_dir is required.');
end

patient_id = char(patient_id);
session    = char(session);

% Defaults
sensor_size_cm     = get_field(CONFIG, 'sensor_size_cm',     [10, 10]);
sensor_standoff_mm = get_field(CONFIG, 'sensor_standoff_mm', 5);
tilt_angle_deg     = get_field(CONFIG, 'tilt_angle_deg',     15);
pml_size           = get_field(CONFIG, 'pml_size',           10);

fprintf('[STEP 01] Interactive sensor placement: %s / %s\n', patient_id, session);
fprintf('[STEP 01]   sensor_size = [%.1f, %.1f] cm, standoff = %.1f mm, tilt = %.1f deg\n', ...
    sensor_size_cm(1), sensor_size_cm(2), sensor_standoff_mm, tilt_angle_deg);

%% ======================== STEP 1: PATHS AND LOAD ========================

processedDir = fullfile(CONFIG.working_dir, 'RayStationFiles', patient_id, ...
                        session, 'processed');
if ~isfolder(processedDir)
    error('place_sensor:MissingDir', ...
        'processed/ directory not found: %s', processedDir);
end

sctFile = fullfile(processedDir, 'sct_resampled.mat');
if ~isfile(sctFile)
    error('place_sensor:MissingSCT', ...
        'sct_resampled.mat not found: %s', sctFile);
end

fprintf('[STEP 01] Loading sct_resampled.mat...\n');
sctData = load(sctFile);
if isfield(sctData, 'sct_resampled')
    sct = sctData.sct_resampled;
elseif isfield(sctData, 'sct')
    sct = sctData.sct;
else
    fnms = fieldnames(sctData);
    sct  = sctData.(fnms{1});
end

if ~isfield(sct, 'bodyMask') || ~isfield(sct, 'spacing')
    error('place_sensor:InvalidSCT', ...
        'sct_resampled must contain bodyMask and spacing.');
end

bodyMask = logical(sct.bodyMask);
gridDims = size(bodyMask);
Ny       = gridDims(1);
Nx       = gridDims(2);
Nz       = gridDims(3);
spacing  = sct.spacing(:)';
dx = spacing(1); dy = spacing(2); dz = spacing(3);

if isfield(sct, 'origin') && ~isempty(sct.origin)
    origin = sct.origin(:)';
else
    origin = [0, 0, 0];
    warning('place_sensor:NoOrigin', ...
        'sct.origin missing — defaulting to [0,0,0].');
end

% Optional CT (for visualisation)
if isfield(sct, 'ct') && ~isempty(sct.ct)
    ctVol = double(sct.ct);
elseif isfield(sct, 'image') && ~isempty(sct.image)
    ctVol = double(sct.image);
elseif isfield(sct, 'HU') && ~isempty(sct.HU)
    ctVol = double(sct.HU);
else
    ctVol = [];
    warning('place_sensor:NoCT', ...
        'CT volume not found in sct_resampled — using bodyMask for background.');
end

fprintf('[STEP 01]   Grid: [%d(Y) x %d(X) x %d(Z)], spacing [%.2f, %.2f, %.2f] mm\n', ...
    Ny, Nx, Nz, dx, dy, dz);

% Field doses — accumulate summed dose loading one at a time
doseFiles = dir(fullfile(processedDir, 'field_dose_*.mat'));
if isempty(doseFiles)
    error('place_sensor:NoDose', ...
        'No field_dose_*.mat files found in %s.', processedDir);
end

fprintf('[STEP 01] Accumulating summed dose from %d field files...\n', ...
    numel(doseFiles));

summedDose    = zeros(gridDims);
beamIsoList   = zeros(numel(doseFiles), 3);  % per-file isocenter
beamGantryDeg = NaN(numel(doseFiles), 1);
nValidBeams   = 0;
for i = 1:numel(doseFiles)
    dpath = fullfile(doseFiles(i).folder, doseFiles(i).name);
    fd    = load(dpath);
    if isfield(fd, 'field_dose')
        fdose = fd.field_dose;
    else
        fnms  = fieldnames(fd);
        fdose = fd.(fnms{1});
    end

    if ~isfield(fdose, 'dose_Gy') || isempty(fdose.dose_Gy)
        warning('place_sensor:EmptyDose', ...
            'Skipping %s (no dose_Gy).', doseFiles(i).name);
        continue;
    end

    d = double(fdose.dose_Gy);
    if ~isequal(size(d), gridDims)
        warning('place_sensor:DoseGridMismatch', ...
            'Skipping %s (size [%s] != grid [%s]).', doseFiles(i).name, ...
            num2str(size(d)), num2str(gridDims));
        clear d fd fdose;
        continue;
    end
    summedDose = summedDose + d;

    nValidBeams = nValidBeams + 1;
    if isfield(fdose, 'isocenter') && numel(fdose.isocenter) >= 3
        beamIsoList(nValidBeams, :) = fdose.isocenter(1:3);
    else
        beamIsoList(nValidBeams, :) = NaN;
    end
    if isfield(fdose, 'gantry_angle') && ~isempty(fdose.gantry_angle)
        beamGantryDeg(nValidBeams) = double(fdose.gantry_angle);
    end

    clear d fd fdose;
end
beamIsoList   = beamIsoList(1:nValidBeams, :);
beamGantryDeg = beamGantryDeg(1:nValidBeams);

if max(summedDose(:)) <= 0
    error('place_sensor:ZeroDose', 'Summed dose is empty/zero.');
end

%% ======================== STEP 2: ISO SLICE ========================

% Dose-weighted centroid Z index
zProfile = squeeze(sum(sum(summedDose, 1), 2));
zIdxVec  = (1:Nz)';
zMass    = sum(zProfile);
if zMass > 0
    isoIz = round(sum(zIdxVec .* zProfile) / zMass);
else
    isoIz = round(Nz/2);
end
isoIz = max(1, min(Nz, isoIz));

% Mean isocenter (mm), used for beam projections + tilt direction
validIso = all(~isnan(beamIsoList), 2);
if any(validIso)
    iso_mm = mean(beamIsoList(validIso, :), 1);
else
    iso_mm = origin + [(Nx/2)*dx, (Ny/2)*dy, (isoIz-1)*dz];
    warning('place_sensor:NoIso', ...
        'No beam isocenter found — using grid centroid: [%.1f %.1f %.1f] mm.', ...
        iso_mm);
end
isoX_mm = iso_mm(1);
isoY_mm = iso_mm(2);
isoZ_mm = iso_mm(3);

fprintf('[STEP 01]   isoIz = %d (%.1f mm), iso_mm = [%.1f, %.1f, %.1f]\n', ...
    isoIz, origin(3) + (isoIz-1)*dz, isoX_mm, isoY_mm, isoZ_mm);

%% ======================== STEP 3: COMPOSITE FIGURE ========================

xAxis_mm = origin(1) + (0:Nx-1) * dx;
yAxis_mm = origin(2) + (0:Ny-1) * dy;

% CT slice with soft-tissue window centred at 40 HU, width 400 HU
if ~isempty(ctVol)
    ctSlice = ctVol(:, :, isoIz);
else
    ctSlice = double(bodyMask(:, :, isoIz));
end
hu_centre = 40;
hu_width  = 400;
hu_lo     = hu_centre - hu_width/2;
hu_hi     = hu_centre + hu_width/2;
ctNorm    = (ctSlice - hu_lo) / (hu_hi - hu_lo);
ctNorm    = max(0, min(1, ctNorm));

% Dose overlay, 5% threshold of slice
doseSlice = summedDose(:, :, isoIz);
doseMax   = max(summedDose(:));
doseNorm  = doseSlice / doseMax;
doseAlpha = 0.45 * (doseNorm >= 0.05);

% Body contour
bodySlice = bodyMask(:, :, isoIz);
boundaries = bwboundaries(bodySlice);

hFig = figure('Name', sprintf('place_sensor — %s / %s', patient_id, session), ...
              'Color', 'w');
hAx  = axes('Parent', hFig);
hold(hAx, 'on');

% CT background
imagesc(hAx, xAxis_mm, yAxis_mm, ctNorm);
colormap(hAx, gray);
set(hAx, 'CLim', [0, 1]);

% Dose hot overlay using independent image with alpha
hotMap   = hot(256);
doseIdx  = uint8(round(doseNorm * 255)) + 1;
doseIdx  = min(max(doseIdx, 1), 256);
doseRGB  = ind2rgb(doseIdx, hotMap);
hDose    = image(hAx, xAxis_mm, yAxis_mm, doseRGB);
set(hDose, 'AlphaData', doseAlpha);

% Body contour
for b = 1:numel(boundaries)
    bnd = boundaries{b};
    % bwboundaries returns [row, col] = [yIdx, xIdx]
    bx_mm = origin(1) + (bnd(:, 2) - 1) * dx;
    by_mm = origin(2) + (bnd(:, 1) - 1) * dy;
    plot(hAx, bx_mm, by_mm, 'w-', 'LineWidth', 1.5);
end

% Beam projections through isocenter (SAD = 1000 mm, IEC convention)
SAD = 1000;
nBeams   = numel(beamGantryDeg);
beamCols = lines(max(nBeams, 1));
legendH  = gobjects(0);
legendL  = {};
for b = 1:nBeams
    ga = beamGantryDeg(b);
    if isnan(ga)
        continue;
    end
    gaRad = deg2rad(ga);
    srcX  = isoX_mm + SAD * sin(gaRad);
    srcY  = isoY_mm - SAD * cos(gaRad);
    tgtX  = isoX_mm - SAD * sin(gaRad);
    tgtY  = isoY_mm + SAD * cos(gaRad);
    h = plot(hAx, [srcX, tgtX], [srcY, tgtY], '-', ...
             'Color', beamCols(b, :), 'LineWidth', 1.5);
    legendH(end+1) = h;                                              %#ok<AGROW>
    legendL{end+1} = sprintf('Beam %d: %.0f^\\circ', b, ga);         %#ok<AGROW>
end

% Isocenter marker
plot(hAx, isoX_mm, isoY_mm, 'y+', 'MarkerSize', 12, 'LineWidth', 2);

axis(hAx, 'image');
set(hAx, 'YDir', 'normal');
xlabel(hAx, 'X (mm) — Lateral');
ylabel(hAx, 'Y (mm) — Anterior-Posterior');
title(hAx, sprintf('Sensor placement — slice iz=%d (z=%.1f mm) — draw rectangle for sensor centre', ...
    isoIz, origin(3) + (isoIz-1)*dz));
if ~isempty(legendH)
    legend(hAx, legendH, legendL, 'Location', 'best', 'TextColor', 'k');
end

%% ======================== STEP 4: INTERACTIVE PLACEMENT ========================

fprintf('[STEP 01] Draw a rectangle to choose sensor centre (size is taken from CONFIG)...\n');
roi = drawrectangle(hAx, 'Label', 'Sensor', 'Color', [0.2 0.8 0.2], ...
                    'LineWidth', 2);
wait(roi);

if ~isvalid(roi) || isempty(roi.Position)
    error('place_sensor:NoROI', 'No ROI selected.');
end

pos    = roi.Position;  % [x, y, w, h] in axis data units (mm)
cx_mm  = pos(1) + pos(3)/2;
cy_mm  = pos(2) + pos(4)/2;
cz_mm  = origin(3) + (isoIz - 1) * dz;

% Convert centre to voxel indices
cxIdx = round((cx_mm - origin(1)) / dx) + 1;
cyIdx = round((cy_mm - origin(2)) / dy) + 1;
czIdx = isoIz;

cxIdx = max(1, min(Nx, cxIdx));
cyIdx = max(1, min(Ny, cyIdx));

fprintf('[STEP 01]   Selected centre (mm): [%.1f, %.1f, %.1f] -> idx [Y=%d, X=%d, Z=%d]\n', ...
    cx_mm, cy_mm, cz_mm, cyIdx, cxIdx, czIdx);

% Snap anteriorly if centre is inside the body
if bodyMask(cyIdx, cxIdx, czIdx)
    warning('place_sensor:CentreInsideBody', ...
        'Selected centre is inside the body — snapping anteriorly to surface.');
    yIdxFirst = find(bodyMask(:, cxIdx, czIdx), 1, 'first');
    if ~isempty(yIdxFirst)
        cyIdx = max(pml_size + 1, yIdxFirst - 1);
        cy_mm = origin(2) + (cyIdx - 1) * dy;
    end
end

%% ======================== STEP 5: FLAT MASK ========================

% Anterior surface map (min Y where body==true per X,Z column)
anteriorSurface = build_anterior_surface(bodyMask);

halfNx = round(sensor_size_cm(1) * 10 / (2 * dx));
halfNz = round(sensor_size_cm(2) * 10 / (2 * dz));

xRange = (cxIdx - halfNx) : (cxIdx + halfNx);
zRange = (czIdx - halfNz) : (czIdx + halfNz);

xRange = xRange(xRange >= pml_size + 1 & xRange <= Nx - pml_size);
zRange = zRange(zRange >= pml_size + 1 & zRange <= Nz - pml_size);

if isempty(xRange) || isempty(zRange)
    warning('place_sensor:OutOfBounds', ...
        'Sensor footprint outside grid after PML clamping — returning invalid placement.');
    placement = empty_placement(gridDims, patient_id, session, CONFIG, tilt_angle_deg);
    return;
end

standoff_vox = max(1, round(sensor_standoff_mm / dy));

% Sensor Y plane derived from anterior surface at (cxIdx, czIdx)
asAtCentre = anteriorSurface(cxIdx, czIdx);
if isnan(asAtCentre)
    % Fallback: minimum surface within footprint
    patch_ant = anteriorSurface(xRange, zRange);
    asAtCentre = min(patch_ant(~isnan(patch_ant)));
    if isempty(asAtCentre)
        warning('place_sensor:NoSurface', ...
            'No anterior surface within footprint — returning invalid placement.');
        placement = empty_placement(gridDims, patient_id, session, CONFIG, tilt_angle_deg);
        return;
    end
end

sensorY = round(asAtCentre) - standoff_vox;
sensorY = max(pml_size + 1, min(Ny - pml_size, sensorY));

sensorMaskFlat = false(gridDims);
sensorMaskFlat(sensorY, xRange, zRange) = true;

fprintf('[STEP 01]   Flat mask: Y=%d, X=[%d,%d], Z=[%d,%d] (%d voxels)\n', ...
    sensorY, xRange(1), xRange(end), zRange(1), zRange(end), ...
    sum(sensorMaskFlat(:)));

%% ======================== STEP 6: TILT ABOUT X-AXIS ========================

[iyF, ixF, izF] = ind2sub(gridDims, find(sensorMaskFlat));

% Voxel -> mm (DICOM patient coords [X Y Z])
ptsMm = [origin(1) + (ixF - 1) * dx, ...
         origin(2) + (iyF - 1) * dy, ...
         origin(3) + (izF - 1) * dz];

cMm  = [origin(1) + (cxIdx   - 1) * dx, ...
        origin(2) + (sensorY - 1) * dy, ...
        origin(3) + (czIdx   - 1) * dz];

signZ = sign(isoZ_mm - cMm(3));
if signZ == 0
    signZ = 1;
end
th = signZ * deg2rad(tilt_angle_deg);

R = [1, 0,        0;
     0, cos(th), -sin(th);
     0, sin(th),  cos(th)];

ptsRot = (R * (ptsMm - cMm)')' + cMm;

% Round back to voxel indices, clamp
ixR = round((ptsRot(:, 1) - origin(1)) / dx) + 1;
iyR = round((ptsRot(:, 2) - origin(2)) / dy) + 1;
izR = round((ptsRot(:, 3) - origin(3)) / dz) + 1;

ixR = max(1, min(Nx, ixR));
iyR = max(1, min(Ny, iyR));
izR = max(1, min(Nz, izR));

sensorMaskTilt = false(gridDims);
linIdx = sub2ind(gridDims, iyR, ixR, izR);
sensorMaskTilt(linIdx) = true;

% Normal (rotation of -Y in patient coords)
nHat = R * [0; -1; 0];
nHat = nHat / norm(nHat);

fprintf('[STEP 01]   Tilted (theta=%.2f deg, sign=%+d): %d voxels\n', ...
    rad2deg(th), signZ, sum(sensorMaskTilt(:)));

%% ======================== STEP 7: SNAP TO CONTOUR ========================

sensorMaskSnap = false(gridDims);
for ix = unique(ixR(:))'
    for iz = unique(izR(izR > 0))'
        colMask = sensorMaskTilt(:, ix, iz);
        if ~any(colMask)
            continue;
        end
        yInCol = find(colMask);
        if isempty(yInCol)
            continue;
        end
        % Does the column overlap the body?
        bodyCol = bodyMask(:, ix, iz);
        as_ixiz = anteriorSurface(ix, iz);
        if any(colMask & bodyCol) && ~isnan(as_ixiz)
            shift = round(as_ixiz) - standoff_vox - min(yInCol);
            yNew  = yInCol + shift;
            yNew  = yNew(yNew >= 1 & yNew <= Ny);
            sensorMaskSnap(yNew, ix, iz) = true;
        else
            sensorMaskSnap(yInCol, ix, iz) = true;
        end
    end
end

fprintf('[STEP 01]   Snapped to contour: %d voxels\n', sum(sensorMaskSnap(:)));

%% ======================== STEP 8: VALIDATE ========================

placement_valid = true;

% Remove body overlap
overlap = sensorMaskSnap & bodyMask;
if any(overlap(:))
    warning('place_sensor:BodyOverlap', ...
        'Sensor overlaps body at %d voxels — removing.', sum(overlap(:)));
    sensorMaskSnap = sensorMaskSnap & ~bodyMask;
    placement_valid = false;
end

% PML clearance
[iySV, ixSV, izSV] = ind2sub(gridDims, find(sensorMaskSnap));
if isempty(iySV)
    warning('place_sensor:EmptyMask', ...
        'Final sensor mask is empty.');
    placement_valid = false;
elseif any(iySV <= pml_size) || any(iySV > Ny - pml_size) || ...
       any(ixSV <= pml_size) || any(ixSV > Nx - pml_size) || ...
       any(izSV <= pml_size) || any(izSV > Nz - pml_size)
    warning('place_sensor:PMLClearance', ...
        'Sensor voxels lie within PML buffer (pml_size=%d).', pml_size);
    placement_valid = false;
end

if ~isempty(iySV)
    sensor_center_mm = [origin(1) + (mean(ixSV) - 1) * dx, ...
                        origin(2) + (mean(iySV) - 1) * dy, ...
                        origin(3) + (mean(izSV) - 1) * dz];
else
    sensor_center_mm = cMm;
end

fprintf('[STEP 01]   Validation: %d voxels, valid=%s\n', ...
    sum(sensorMaskSnap(:)), mat2str(placement_valid));
fprintf('[STEP 01]   Sensor centre (mm): [%.1f, %.1f, %.1f]\n', sensor_center_mm);
fprintf('[STEP 01]   Sensor normal: [%+.3f, %+.3f, %+.3f]\n', nHat);

%% ======================== STEP 9: VISUALISE ========================

% 9a: overlay tilted sensor on the existing figure (project to slice isoIz)
slMask = squeeze(any(sensorMaskSnap(:, :, max(1, czIdx-2):min(Nz, czIdx+2)), 3));
[snapBnds] = bwboundaries(slMask);
for b = 1:numel(snapBnds)
    bnd  = snapBnds{b};
    bx_mm = origin(1) + (bnd(:, 2) - 1) * dx;
    by_mm = origin(2) + (bnd(:, 1) - 1) * dy;
    patch(hAx, bx_mm, by_mm, [0, 1, 1], ...
        'FaceAlpha', 0.35, 'EdgeColor', 'c', 'LineWidth', 1.5);
end

% 9b: 3D view
hFig3D = figure('Name', sprintf('place_sensor 3D — %s / %s', ...
    patient_id, session), 'Color', 'w');
hAx3   = axes('Parent', hFig3D);
hold(hAx3, 'on');

% Body isosurface (downsample for speed)
ds = max(1, round(min([Ny, Nx, Nz]) / 96));
[Xs, Ys, Zs] = meshgrid(origin(1) + (0:ds:Nx-1) * dx, ...
                        origin(2) + (0:ds:Ny-1) * dy, ...
                        origin(3) + (0:ds:Nz-1) * dz);
bodyDS = bodyMask(1:ds:end, 1:ds:end, 1:ds:end);
try
    pBody = patch(hAx3, isosurface(Xs, Ys, Zs, double(bodyDS), 0.5));
    set(pBody, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.15);
catch ME
    warning('place_sensor:IsosurfaceFail', ...
        'Could not render body isosurface: %s', ME.message);
end

% Sensor scatter
if ~isempty(iySV)
    sx_mm = origin(1) + (ixSV - 1) * dx;
    sy_mm = origin(2) + (iySV - 1) * dy;
    sz_mm = origin(3) + (izSV - 1) * dz;
    scatter3(hAx3, sx_mm, sy_mm, sz_mm, 18, 'b', 'filled');
end

% Beam axis lines
for b = 1:nBeams
    ga = beamGantryDeg(b);
    if isnan(ga)
        continue;
    end
    gaRad = deg2rad(ga);
    srcX  = isoX_mm + SAD * sin(gaRad);
    srcY  = isoY_mm - SAD * cos(gaRad);
    tgtX  = isoX_mm - SAD * sin(gaRad);
    tgtY  = isoY_mm + SAD * cos(gaRad);
    plot3(hAx3, [srcX, tgtX], [srcY, tgtY], [isoZ_mm, isoZ_mm], '-', ...
        'Color', beamCols(b, :), 'LineWidth', 1.2);
end

plot3(hAx3, isoX_mm, isoY_mm, isoZ_mm, 'y+', 'MarkerSize', 12, 'LineWidth', 2);

xlabel(hAx3, 'X (mm)'); ylabel(hAx3, 'Y (mm)'); zlabel(hAx3, 'Z (mm)');
title(hAx3, sprintf('Tilted sensor (%.0f deg) — %d voxels', ...
    rad2deg(th), sum(sensorMaskSnap(:))));
axis(hAx3, 'equal'); grid(hAx3, 'on'); view(hAx3, 35, 25);
lighting(hAx3, 'gouraud'); camlight(hAx3, 'headlight');

%% ======================== STEP 10: PACK + SAVE ========================

placement = struct();
placement.sensor_mask       = sensorMaskSnap;
placement.sensor_center_mm  = sensor_center_mm;
placement.sensor_normal_hat = nHat;
placement.tilt_angle_deg    = rad2deg(th);
placement.placement_valid   = placement_valid;
placement.patient_id        = patient_id;
placement.session           = session;
placement.timestamp         = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
placement.config            = CONFIG;

outFile = fullfile(processedDir, 'sensor_placement.mat');
fprintf('[STEP 01] Saving sensor_placement.mat...\n');
save(outFile, 'placement', '-v7.3');
fprintf('[STEP 01] Done: %s\n', outFile);

end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function val = get_field(s, name, default)
%GET_FIELD  Return struct field if present and non-empty, else default.
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default;
    end
end


function as = build_anterior_surface(bodyMask)
%BUILD_ANTERIOR_SURFACE  Min Y voxel index where body is true, per (X,Z) col.
    [Ny, Nx, Nz] = size(bodyMask);                                       %#ok<ASGLU>
    as = NaN(Nx, Nz);
    for ix = 1:Nx
        for iz = 1:Nz
            y = find(bodyMask(:, ix, iz), 1, 'first');
            if ~isempty(y)
                as(ix, iz) = y;
            end
        end
    end
end


function placement = empty_placement(gridDims, patient_id, session, CONFIG, tilt_deg)
%EMPTY_PLACEMENT  Build an invalid-placement output struct (for warning paths).
    placement = struct();
    placement.sensor_mask       = false(gridDims);
    placement.sensor_center_mm  = [0, 0, 0];
    placement.sensor_normal_hat = [0; -1; 0];
    placement.tilt_angle_deg    = tilt_deg;
    placement.placement_valid   = false;
    placement.patient_id        = patient_id;
    placement.session           = session;
    placement.timestamp         = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
    placement.config            = CONFIG;
end
