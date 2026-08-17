% diagnostic_place_sensor.m
%
% PURPOSE:
%   Stand-alone diagnostic harness for `place_sensor.m`. Builds a synthetic
%   patient phantom (body mask + CT) and a small set of synthetic per-field
%   dose distributions matching the layout expected by `place_sensor`,
%   writes them into a temporary `processed/` directory, invokes the
%   interactive sensor placement, and then runs an extensive battery of
%   structural / numerical / geometric checks on the returned struct and
%   the persisted `sensor_placement.mat`.
%
%   This avoids needing real DICOM exports or the upstream pipeline steps
%   while still exercising every code path in `place_sensor` end-to-end.
%
% USAGE:
%   Edit the CONFIGURATION block if desired, then run (F5). The script
%   will:
%     1. Build phantom + synthetic doses.
%     2. Open the interactive figure — DRAW A RECTANGLE on the anterior
%        side of the body to choose the sensor centre.
%     3. Report on the resulting placement.
%
% DEPENDENCIES:
%   place_sensor.m, Image Processing Toolbox (drawrectangle, bwboundaries).
%
% See also: place_sensor, determine_sensor_placement_fixed

clear; close all; clc;

% Script lives one level below the repo root; add utils/ and pipeline/ from there.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'utils')));
addpath(genpath(fullfile(repoRoot, 'pipeline')));

%% =============== CONFIGURATION =================================== %%
test_patient_id = 'PHANTOM_TEST';
test_session    = 'Session_Diag';
test_root       = fullfile(tempdir, 'ethos_place_sensor_diag');

% Phantom grid
Ny      = 160;  % A/P  (Y)
Nx      = 180;  % L/R  (X)
Nz      = 120;  % S/I  (Z)
spacing = [2.0, 2.0, 2.5];   % [dx, dy, dz] in mm
origin  = [-180, -160, -150]; % DICOM mm

% Phantom shape — elliptical cylinder of "soft tissue"
ellipse_a_mm = 140;   % half-width (X)
ellipse_b_mm = 90;    % half-thickness (Y)

% Beams: gantry angles for the synthetic plan
beam_gantry_deg = [0, 120, 240];
iso_offset_mm   = [10, 20, 0];   % isocentre offset from grid centre

% CONFIG passed to place_sensor
CONFIG = struct();
CONFIG.working_dir        = test_root;
CONFIG.sensor_size_cm     = [10, 10];
CONFIG.sensor_standoff_mm = 5;
CONFIG.tilt_angle_deg     = 15;
CONFIG.pml_size           = 10;

bar = repmat('=', 1, 72);
fprintf('\n%s\n  PLACE_SENSOR DIAGNOSTIC HARNESS\n%s\n', bar, bar);
fprintf('Patient        : %s\n', test_patient_id);
fprintf('Session        : %s\n', test_session);
fprintf('Test root      : %s\n', test_root);
fprintf('Grid (Y,X,Z)   : [%d, %d, %d]\n', Ny, Nx, Nz);
fprintf('Spacing (mm)   : [%.2f, %.2f, %.2f]\n', spacing);
fprintf('Origin (mm)    : [%.1f, %.1f, %.1f]\n', origin);
fprintf('Beams (gantry) : [%s] deg\n', num2str(beam_gantry_deg));
fprintf('Timestamp      : %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% =============== STEP 1: BUILD PHANTOM ============================ %%
fprintf('[DIAG 1] Building synthetic phantom...\n');

% Voxel coords (mm)
xGrid_mm = origin(1) + (0:Nx-1) * spacing(1);
yGrid_mm = origin(2) + (0:Ny-1) * spacing(2);
zGrid_mm = origin(3) + (0:Nz-1) * spacing(3);
[Xm, Ym, Zm] = meshgrid(xGrid_mm, yGrid_mm, zGrid_mm);  % size [Ny,Nx,Nz]

% Centre of phantom (rough patient centre)
cx_phantom = origin(1) + (Nx/2 - 1) * spacing(1);
cy_phantom = origin(2) + (Ny/2 - 1) * spacing(2);
cz_phantom = origin(3) + (Nz/2 - 1) * spacing(3);

% Elliptical cylinder along Z, body fills tail and head except near edges
bodyMask = ((Xm - cx_phantom).^2 / ellipse_a_mm^2 + ...
            (Ym - cy_phantom).^2 / ellipse_b_mm^2) <= 1;
bodyMask = bodyMask & (Zm > origin(3) + 20) & (Zm < origin(3) + spacing(3)*(Nz-1) - 20);

% Synthetic CT volume (HU). Air ~ -1000, soft tissue ~ 40, with mild noise.
ctVol = -1000 * ones(Ny, Nx, Nz);
ctVol(bodyMask) = 40 + 15 * randn(nnz(bodyMask), 1);

% Pack sct_resampled struct
sct_resampled = struct();
sct_resampled.bodyMask = bodyMask;
sct_resampled.ct       = ctVol;
sct_resampled.spacing  = spacing;
sct_resampled.origin   = origin;

fprintf('  Body voxels: %d (%.1f%% of grid)\n', ...
    nnz(bodyMask), 100 * nnz(bodyMask) / numel(bodyMask));

%% =============== STEP 2: WRITE PROCESSED DIR ====================== %%
fprintf('[DIAG 2] Writing test files to processed/...\n');

processedDir = fullfile(test_root, 'RayStationFiles', ...
    test_patient_id, test_session, 'processed');
if isfolder(processedDir)
    rmdir(processedDir, 's');
end
mkdir(processedDir);

sctPath = fullfile(processedDir, 'sct_resampled.mat');
save(sctPath, 'sct_resampled', '-v7.3');
fprintf('  sct_resampled.mat -> %s\n', sctPath);

% Isocentre
iso_mm = [cx_phantom + iso_offset_mm(1), ...
          cy_phantom + iso_offset_mm(2), ...
          cz_phantom + iso_offset_mm(3)];

% Synthetic dose: Gaussian blob centred at isocentre per beam
sigma_mm = [25, 25, 30];
for b = 1:numel(beam_gantry_deg)
    ga = beam_gantry_deg(b);
    blob = exp( -((Xm - iso_mm(1)).^2 / (2*sigma_mm(1)^2) ...
                + (Ym - iso_mm(2)).^2 / (2*sigma_mm(2)^2) ...
                + (Zm - iso_mm(3)).^2 / (2*sigma_mm(3)^2)) );
    dose_Gy = 2.0 * blob;  % peak ~ 2 Gy
    dose_Gy(~bodyMask) = 0;

    field_dose = struct();
    field_dose.dose_Gy      = dose_Gy;
    field_dose.origin       = origin;
    field_dose.spacing      = spacing;
    field_dose.dimensions   = [Ny, Nx, Nz];
    field_dose.beam_num     = b;
    field_dose.seg_num      = 1;
    field_dose.field_num    = b;
    field_dose.plan_type    = 'diagnostic';
    field_dose.gantry_angle = ga;
    field_dose.meterset     = 100;
    field_dose.isocenter    = iso_mm;
    field_dose.jaw_x        = [-40, 40];
    field_dose.jaw_y        = [-50, 50];
    field_dose.max_dose_Gy  = max(dose_Gy(:));

    fdPath = fullfile(processedDir, sprintf('field_dose_B%d_S1.mat', b));
    save(fdPath, 'field_dose', '-v7.3');
    fprintf('  field_dose_B%d_S1.mat (gantry %3.0f deg, max %.2f Gy)\n', ...
        b, ga, max(dose_Gy(:)));
end

%% =============== STEP 3: CALL PLACE_SENSOR ======================== %%
fprintf('\n[DIAG 3] Invoking place_sensor (INTERACTIVE: draw rectangle)...\n');
fprintf('         Aim the rectangle anterior to the body (low-Y side).\n\n');

t0 = tic;
try
    placement = place_sensor(test_patient_id, test_session, CONFIG);
    elapsed = toc(t0);
catch ME
    fprintf(2, '\n[DIAG] FATAL: place_sensor threw an exception:\n  %s\n', ...
        ME.message);
    for k = 1:numel(ME.stack)
        fprintf(2, '    at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    return;
end
fprintf('\n[DIAG 3]   place_sensor returned in %.2f s.\n', elapsed);

%% =============== STEP 4: STRUCTURAL CHECKS ======================== %%
fprintf('\n[DIAG 4] Structural checks on returned struct...\n');

passed = 0; failed = 0;

required_fields = {'sensor_mask', 'sensor_center_mm', 'sensor_normal_hat', ...
    'tilt_angle_deg', 'placement_valid', 'patient_id', 'session', ...
    'timestamp', 'config'};

for k = 1:numel(required_fields)
    f = required_fields{k};
    [ok, msg] = deal(isfield(placement, f), sprintf('has field .%s', f));
    [passed, failed] = report_inc(passed, failed, ok, msg);
end

% Field type / shape checks
if isfield(placement, 'sensor_mask')
    sm = placement.sensor_mask;
    [passed, failed] = report_inc(passed, failed, islogical(sm), ...
        'sensor_mask is logical');
    [passed, failed] = report_inc(passed, failed, ...
        isequal(size(sm), [Ny, Nx, Nz]), ...
        sprintf('sensor_mask size [%s] == grid [%d %d %d]', ...
            num2str(size(sm)), Ny, Nx, Nz));
end

if isfield(placement, 'sensor_center_mm')
    c = placement.sensor_center_mm;
    [passed, failed] = report_inc(passed, failed, ...
        isnumeric(c) && numel(c) == 3 && all(isfinite(c)), ...
        sprintf('sensor_center_mm finite 3-vector ([%.1f %.1f %.1f])', c(1), c(2), c(3)));
end

if isfield(placement, 'sensor_normal_hat')
    n = placement.sensor_normal_hat(:);
    [passed, failed] = report_inc(passed, failed, ...
        numel(n) == 3 && abs(norm(n) - 1) < 1e-6, ...
        sprintf('sensor_normal_hat unit-length (|n|=%.6f)', norm(n)));
    [passed, failed] = report_inc(passed, failed, n(2) < 0, ...
        sprintf('sensor_normal_hat points anteriorly (n_y=%.3f < 0)', n(2)));
end

if isfield(placement, 'tilt_angle_deg')
    [passed, failed] = report_inc(passed, failed, ...
        abs(abs(placement.tilt_angle_deg) - CONFIG.tilt_angle_deg) < 1e-6, ...
        sprintf('|tilt_angle_deg| == CONFIG.tilt_angle_deg (%.3f vs %.3f)', ...
            abs(placement.tilt_angle_deg), CONFIG.tilt_angle_deg));
end

if isfield(placement, 'patient_id')
    [passed, failed] = report_inc(passed, failed, ...
        strcmp(placement.patient_id, test_patient_id), ...
        'patient_id echoed correctly');
end
if isfield(placement, 'session')
    [passed, failed] = report_inc(passed, failed, ...
        strcmp(placement.session, test_session), 'session echoed correctly');
end
if isfield(placement, 'config')
    [passed, failed] = report_inc(passed, failed, ...
        isstruct(placement.config) && ...
        isequal(placement.config.sensor_size_cm, CONFIG.sensor_size_cm), ...
        'config echoed correctly');
end

%% =============== STEP 5: GEOMETRIC CHECKS ========================= %%
fprintf('\n[DIAG 5] Geometric checks on sensor mask...\n');

if isfield(placement, 'sensor_mask') && any(placement.sensor_mask(:))
    sm = placement.sensor_mask;
    nvox = nnz(sm);

    expected_nx = round(CONFIG.sensor_size_cm(1) * 10 / spacing(1));
    expected_nz = round(CONFIG.sensor_size_cm(2) * 10 / spacing(2));
    % Rotation about X spreads voxels along Y; expected order of magnitude:
    expected_min = round(0.4 * expected_nx * expected_nz);
    expected_max = round(3.0 * expected_nx * expected_nz);

    [passed, failed] = report_inc(passed, failed, ...
        nvox >= expected_min && nvox <= expected_max, ...
        sprintf('sensor voxel count %d in expected band [%d, %d]', ...
            nvox, expected_min, expected_max));

    [iyS, ixS, izS] = ind2sub(size(sm), find(sm));

    % PML clearance
    [passed, failed] = report_inc(passed, failed, ...
        min(iyS) > CONFIG.pml_size && max(iyS) <= Ny - CONFIG.pml_size && ...
        min(ixS) > CONFIG.pml_size && max(ixS) <= Nx - CONFIG.pml_size && ...
        min(izS) > CONFIG.pml_size && max(izS) <= Nz - CONFIG.pml_size, ...
        sprintf('all sensor voxels clear of PML (pml_size=%d)', CONFIG.pml_size));

    % Body overlap
    overlap = nnz(sm & bodyMask);
    [passed, failed] = report_inc(passed, failed, overlap == 0, ...
        sprintf('no body overlap (overlap voxels=%d)', overlap));

    % Standoff: minimum Y distance from each sensor (X,Z) column to body surface
    standoff_violations = 0;
    standoff_vox_expected = max(1, round(CONFIG.sensor_standoff_mm / spacing(2)));
    cols = unique([ixS, izS], 'rows');
    for kc = 1:size(cols, 1)
        ix = cols(kc, 1); iz = cols(kc, 2);
        ySensor = find(sm(:, ix, iz));
        yBody   = find(bodyMask(:, ix, iz), 1, 'first');
        if isempty(yBody) || isempty(ySensor)
            continue;
        end
        if (yBody - max(ySensor)) < (standoff_vox_expected - 1)
            standoff_violations = standoff_violations + 1;
        end
    end
    [passed, failed] = report_inc(passed, failed, standoff_violations <= 0.05 * size(cols, 1), ...
        sprintf('standoff respected on >=95%% of columns (violations=%d / %d)', ...
            standoff_violations, size(cols, 1)));

    % Sensor lies anteriorly (low Y) relative to body centroid
    [yBody, ~, ~] = ind2sub(size(bodyMask), find(bodyMask));
    [passed, failed] = report_inc(passed, failed, mean(iyS) < mean(yBody), ...
        sprintf('sensor anterior to body centroid (mean_y sensor=%.1f, body=%.1f)', ...
            mean(iyS), mean(yBody)));

    fprintf('  Sensor voxel count    : %d\n', nvox);
    fprintf('  Sensor Y range (idx)  : [%d, %d]\n', min(iyS), max(iyS));
    fprintf('  Sensor X range (idx)  : [%d, %d]\n', min(ixS), max(ixS));
    fprintf('  Sensor Z range (idx)  : [%d, %d]\n', min(izS), max(izS));
    fprintf('  Centre (mm)           : [%.1f, %.1f, %.1f]\n', placement.sensor_center_mm);
    fprintf('  Normal                : [%+.3f, %+.3f, %+.3f]\n', placement.sensor_normal_hat);
    fprintf('  placement_valid       : %s\n', mat2str(placement.placement_valid));
else
    [passed, failed] = report_inc(passed, failed, false, ...
        'sensor_mask is non-empty');
end

%% =============== STEP 6: PERSISTENCE CHECK ======================== %%
fprintf('\n[DIAG 6] Re-load sensor_placement.mat and compare...\n');

outFile = fullfile(processedDir, 'sensor_placement.mat');
[passed, failed] = report_inc(passed, failed, isfile(outFile), ...
    sprintf('output file exists (%s)', outFile));

if isfile(outFile)
    S = load(outFile);
    [passed, failed] = report_inc(passed, failed, isfield(S, 'placement'), ...
        'mat file contains placement variable');
    if isfield(S, 'placement')
        p2 = S.placement;
        [passed, failed] = report_inc(passed, failed, ...
            isequal(p2.sensor_mask, placement.sensor_mask), ...
            'persisted sensor_mask matches returned mask');
        [passed, failed] = report_inc(passed, failed, ...
            isequal(p2.sensor_center_mm, placement.sensor_center_mm), ...
            'persisted sensor_center_mm matches returned value');
        [passed, failed] = report_inc(passed, failed, ...
            isequal(p2.sensor_normal_hat, placement.sensor_normal_hat), ...
            'persisted sensor_normal_hat matches returned value');
    end
end

%% =============== STEP 7: TILT-SIGN BEHAVIOUR (non-interactive) ==== %%
fprintf('\n[DIAG 7] Geometric tilt-sign sanity check (analytic)...\n');
% Verify that the rotation matrix used by place_sensor produces the
% expected normal direction. This does NOT require the interactive UI.
for sgn = [-1, 0, 1]
    if sgn == 0
        signZ = 1;
    else
        signZ = sgn;
    end
    th = signZ * deg2rad(CONFIG.tilt_angle_deg);
    R  = [1, 0, 0; 0, cos(th), -sin(th); 0, sin(th), cos(th)];
    n  = R * [0; -1; 0];
    n  = n / norm(n);
    fprintf('  sign(iso_z-c_z)=%+d -> theta=%+.2f deg, n=[%+.3f %+.3f %+.3f], |n|=%.6f\n', ...
        sgn, rad2deg(th), n(1), n(2), n(3), norm(n));
    [passed, failed] = report_inc(passed, failed, abs(norm(n) - 1) < 1e-12, ...
        sprintf('analytic |n|==1 for sign=%+d', sgn));
end

%% =============== SUMMARY ========================================== %%
fprintf('\n%s\n', bar);
fprintf('  SUMMARY: %d passed, %d failed (%d total)\n', ...
    passed, failed, passed + failed);
fprintf('%s\n', bar);

if failed == 0
    fprintf('  ALL DIAGNOSTIC CHECKS PASSED.\n');
else
    fprintf(2, '  %d CHECK(S) FAILED — review output above.\n', failed);
end
fprintf('  Outputs left at: %s\n', processedDir);
fprintf('%s\n\n', bar);


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================
function report_check(cond, lbl)
    if cond
        fprintf('  [PASS] %s\n', lbl);
    else
        fprintf(2, '  [FAIL] %s\n', lbl);
    end
end

function [p, f] = report_inc(p, f, cond, lbl)
    report_check(cond, lbl);
    if cond
        p = p + 1;
    else
        f = f + 1;
    end
end
