%% =========================================================================
%  RUN_STANDALONE_SIMULATION.m
%  Standalone k-Wave Photoacoustic Forward + Time-Reversal Simulation
%  =========================================================================
%
%  PURPOSE:
%  Self-contained script for running a single-field (or total-dose)
%  photoacoustic simulation without the full ETHOS pipeline.  Loads a
%  .mat dose file and sct_resampled.mat from the current directory,
%  builds the acoustic medium, places a sensor, runs the k-Wave forward
%  simulation, and performs time-reversal reconstruction.
%
%  CONFIGURABLE OPTIONS (see CONFIGURATION section):
%    - Dose file selection (default: total_rs_dose.mat)
%    - Sensor mode: 'full_anterior_plane' or 'dose_based'
%    - Tissue heterogeneity: 'uniform', 'threshold_1' (9 tissues),
%                            'threshold_2' (4 tissues)
%    - Individual toggles to disable heterogeneity in specific
%      acoustic properties (density, sound speed, attenuation, Gruneisen)
%    - GPU/CPU selection, PML size, CFL number, TR iterations
%
%  REQUIRED FILES (in working directory):
%    - total_rs_dose.mat  (contains total_rs_dose: 3D Gy array)
%    - sct_resampled.mat  (contains sct_resampled struct with cubeHU,
%                          cubeDensity, bodyMask, couchMask, spacing, etc.)
%
%  PREREQUISITES:
%    - MATLAB R2022a or later
%    - k-Wave Toolbox (http://www.k-wave.org)
%    - Image Processing Toolbox (optional, for dose_based sensor mode)
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: February 2026
%  VERSION: 1.0
%
%  See also: run_single_field_simulation, determine_sensor_mask
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION ================================

% --- File Selection ---
% Path to the dose .mat file. Set to a specific path or leave as-is to
% load from the current working directory.
CONFIG.dose_file = 'total_rs_dose.mat';      % dose file to load
CONFIG.sct_file  = 'sct_resampled.mat';      % CT / tissue data file

% --- Sensor Placement Mode ---
%   'full_anterior_plane' : Entire anteriormost voxel plane is the sensor.
%                           Simple, no geometry dependence, good baseline.
%   'dose_based'          : Sensor placed on beam-exit side based on dose
%                           extent and gantry angle (original pipeline logic).
CONFIG.sensor_mode = 'full_anterior_plane';  % 'full_anterior_plane' | 'dose_based'

% --- Gantry Angle (used only for dose_based sensor mode) ---
% If using total dose (not a single field), set a representative angle.
% 0 deg = beam from anterior; 180 deg = beam from posterior.
CONFIG.gantry_angle = 180;  % degrees (only used when sensor_mode = 'dose_based')

% --- Tissue Heterogeneity ---
%   'uniform'       : Homogeneous water-like medium everywhere
%   'threshold_1'   : 9-tissue model (air, lung, fat, water, blood,
%                     muscle, soft tissue, bone, metal)
%   'threshold_2'   : 4-tissue model (water, fat, soft tissue, bone)
CONFIG.gruneisen_method = 'uniform';

% --- Per-Property Heterogeneity Overrides ---
% When gruneisen_method is NOT 'uniform', you can selectively force
% individual properties to be spatially constant while keeping others
% heterogeneous.  Set to true to DISABLE heterogeneity for that property.
%   e.g., force_uniform_density = true means density = 1000 kg/m^3
%         everywhere, even if threshold_2 is selected.
CONFIG.force_uniform_density     = false;
CONFIG.force_uniform_sound_speed = false;
CONFIG.force_uniform_attenuation = false;
CONFIG.force_uniform_gruneisen   = false;

% --- Uniform Property Values (used when uniform or force_uniform_*) ---
CONFIG.uniform_density      = 1000;    % kg/m^3  (water)
CONFIG.uniform_sound_speed  = 1540;    % m/s     (soft tissue average)
CONFIG.uniform_alpha_coeff  = 0.5;     % dB/MHz^y/cm
CONFIG.uniform_alpha_power  = 1.1;     % exponent
CONFIG.uniform_gruneisen    = 1.0;     % dimensionless

% --- Simulation Parameters ---
CONFIG.dose_per_pulse_cGy     = 0.16;   % cGy per LINAC pulse
CONFIG.meterset               = 100;    % MU (monitor units) for total dose
CONFIG.pml_size               = 10;     % PML thickness (voxels)
CONFIG.cfl_number             = 0.3;    % CFL stability number
CONFIG.use_gpu                = true;   % Use GPU acceleration
CONFIG.num_time_reversal_iter = 1;      % Time-reversal iterations

% --- Output ---
CONFIG.save_results = true;             % Save reconstruction to .mat
CONFIG.output_file  = 'standalone_recon_results.mat';
CONFIG.plot_results = true;             % Show comparison figures

%% ========================= LOAD DATA ====================================

fprintf('=========================================================\n');
fprintf('  Standalone k-Wave Photoacoustic Simulation\n');
fprintf('=========================================================\n');
fprintf('  Dose file:       %s\n', CONFIG.dose_file);
fprintf('  SCT file:        %s\n', CONFIG.sct_file);
fprintf('  Sensor mode:     %s\n', CONFIG.sensor_mode);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  GPU:             %s\n', mat2str(CONFIG.use_gpu));
fprintf('=========================================================\n\n');

% ---- Load dose ----
fprintf('[1/6] Loading dose data...\n');
if ~isfile(CONFIG.dose_file)
    error('Dose file not found: %s', CONFIG.dose_file);
end
dose_data = load(CONFIG.dose_file);

% Auto-detect the dose variable name
dose_fields = fieldnames(dose_data);
if isfield(dose_data, 'total_rs_dose')
    doseGrid = dose_data.total_rs_dose;
    fprintf('       Loaded variable: total_rs_dose\n');
elseif isfield(dose_data, 'dose_Gy')
    doseGrid = dose_data.dose_Gy;
    fprintf('       Loaded variable: dose_Gy\n');
elseif length(dose_fields) == 1
    doseGrid = dose_data.(dose_fields{1});
    fprintf('       Loaded variable: %s\n', dose_fields{1});
else
    fprintf('       Available variables: %s\n', strjoin(dose_fields, ', '));
    error('Cannot auto-detect dose variable. Expected total_rs_dose or dose_Gy.');
end

if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
    error('Dose data must be a 3D numeric array. Got %s with %d dimensions.', ...
        class(doseGrid), ndims(doseGrid));
end

gridSize = size(doseGrid);
Nx = gridSize(1);
Ny = gridSize(2);
Nz = gridSize(3);
fprintf('       Grid size: [%d x %d x %d]\n', Nx, Ny, Nz);
fprintf('       Dose range: [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));

% ---- Load SCT ----
fprintf('[2/6] Loading SCT data...\n');
if ~isfile(CONFIG.sct_file)
    error('SCT file not found: %s', CONFIG.sct_file);
end
sct_data = load(CONFIG.sct_file);

if isfield(sct_data, 'sct_resampled')
    sct = sct_data.sct_resampled;
else
    error('sct_resampled variable not found in %s', CONFIG.sct_file);
end

% Validate SCT has required fields
required_sct_fields = {'cubeHU', 'spacing'};
for i = 1:length(required_sct_fields)
    if ~isfield(sct, required_sct_fields{i})
        error('sct_resampled missing required field: %s', required_sct_fields{i});
    end
end

% Extract spacing (mm -> m for k-Wave)
spacing_mm = sct.spacing(:)';
dx = spacing_mm(1) / 1000;
dy = spacing_mm(2) / 1000;
dz = spacing_mm(3) / 1000;

% Validate grid dimensions match
sctSize = size(sct.cubeHU);
if ~isequal(gridSize, sctSize)
    error('Dose grid [%d %d %d] does not match SCT grid [%d %d %d].', ...
        Nx, Ny, Nz, sctSize(1), sctSize(2), sctSize(3));
end

fprintf('       Spacing: [%.2f, %.2f, %.2f] mm\n', spacing_mm);
fprintf('       HU range: [%.0f, %.0f]\n', min(sct.cubeHU(:)), max(sct.cubeHU(:)));

%% ========================= CREATE ACOUSTIC MEDIUM ========================

fprintf('[3/6] Creating acoustic medium (method: %s)...\n', CONFIG.gruneisen_method);

medium = create_medium(sct, CONFIG);

fprintf('       Density range:     [%.0f, %.0f] kg/m^3\n', ...
    min(medium.density(:)), max(medium.density(:)));
fprintf('       Sound speed range: [%.0f, %.0f] m/s\n', ...
    min(medium.sound_speed(:)), max(medium.sound_speed(:)));
fprintf('       Alpha coeff range: [%.4f, %.4f] dB/MHz^y/cm\n', ...
    min(medium.alpha_coeff(:)), max(medium.alpha_coeff(:)));
fprintf('       Gruneisen range:   [%.4f, %.4f]\n', ...
    min(medium.gruneisen(:)), max(medium.gruneisen(:)));

% Print override status
overrides = {};
if CONFIG.force_uniform_density,     overrides{end+1} = 'density';     end
if CONFIG.force_uniform_sound_speed, overrides{end+1} = 'sound_speed'; end
if CONFIG.force_uniform_attenuation, overrides{end+1} = 'attenuation'; end
if CONFIG.force_uniform_gruneisen,   overrides{end+1} = 'gruneisen';   end
if ~isempty(overrides)
    fprintf('       Forced uniform:    {%s}\n', strjoin(overrides, ', '));
end

%% ========================= INITIAL PRESSURE p0 ==========================

fprintf('[4/6] Computing initial pressure...\n');

% Pulse calculation
meterset = CONFIG.meterset;
num_pulses = ceil(meterset / CONFIG.dose_per_pulse_cGy);
dose_per_pulse = doseGrid / num_pulses;

% p0(r) = D_per_pulse(r) * Gamma(r) * rho(r)
p0 = dose_per_pulse .* medium.gruneisen .* medium.density;

fprintf('       Meterset: %.2f MU -> %d pulses\n', meterset, num_pulses);
fprintf('       Max dose per pulse: %.6f Gy\n', max(dose_per_pulse(:)));
fprintf('       p0 range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));

% Check for significant dose
doseThreshold = 0.01 * max(doseGrid(:));
doseMask = doseGrid > doseThreshold;

if ~any(doseMask(:)) || max(p0(:)) == 0
    warning('No significant dose or zero initial pressure. Aborting.');
    return;
end

%% ========================= SENSOR PLACEMENT ==============================

fprintf('[5/6] Placing sensor (mode: %s)...\n', CONFIG.sensor_mode);

sensor = struct();

switch lower(CONFIG.sensor_mode)
    case 'full_anterior_plane'
        % ---- Full anteriormost plane ----
        % The sensor is the entire first slice along Y (dim 1),
        % which corresponds to the most anterior plane in DICOM.
        % Place it at y=1 (anteriormost row).
        sensor.mask = zeros(Nx, Ny, Nz);
        sensor.mask(1, :, :) = 1;  % first X-plane (anterior)

        % If bodyMask is available, find the actual anterior surface
        if isfield(sct, 'bodyMask')
            body = sct.bodyMask;
            if isfield(sct, 'couchMask')
                body = body & ~sct.couchMask;
            end
            % Find the minimum Y index that has any body voxel
            % (most anterior body surface)
            y_has_body = squeeze(any(any(body, 1), 3));  % [1 x Ny] logical
            anterior_y = find(y_has_body, 1, 'first');
            if ~isempty(anterior_y)
                % Place sensor a few voxels anterior to the body surface
                sensor_y = max(1, anterior_y - 3);
            else
                sensor_y = 1;
            end
        else
            sensor_y = 1;
        end

        sensor.mask = zeros(Nx, Ny, Nz);
        sensor.mask(:, sensor_y, :) = 1;

        fprintf('       Full anterior plane at Y index %d\n', sensor_y);
        fprintf('       Sensor plane size: %d x %d = %d voxels\n', ...
            Nx, Nz, Nx * Nz);

    case 'dose_based'
        % ---- Dose-extent-based placement (original pipeline logic) ----
        sensor = place_sensor_for_field(doseMask, Nx, Ny, Nz, CONFIG.gantry_angle);
        fprintf('       Gantry angle: %.1f deg\n', CONFIG.gantry_angle);

    otherwise
        error('Unknown sensor_mode: %s. Use ''full_anterior_plane'' or ''dose_based''.', ...
            CONFIG.sensor_mode);
end

numSensorPts = sum(sensor.mask(:));
fprintf('       Total sensor voxels: %d\n', numSensorPts);

if numSensorPts == 0
    warning('Sensor mask is empty. Aborting.');
    return;
end

%% ========================= k-WAVE SIMULATION ============================

fprintf('[6/6] Running k-Wave simulation...\n');

% ---- Create k-Wave grid ----
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% CFL-stable time step
maxC = max(medium.sound_speed(:));
minC = min(medium.sound_speed(medium.sound_speed > 0));
dt   = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;

% Simulation time: 2.5x grid diagonal traversal at minimum speed
gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
simTime  = 2.5 * gridDiag / minC;
Nt       = ceil(simTime / dt);

kgrid.dt = dt;
kgrid.Nt = Nt;

fprintf('       dt = %.2e s, Nt = %d, T_sim = %.2e s\n', dt, Nt, simTime);

% ---- k-Wave medium struct ----
kmedium = struct();
kmedium.density     = medium.density;
kmedium.sound_speed = medium.sound_speed;
kmedium.alpha_coeff = medium.alpha_coeff;
kmedium.alpha_power = medium.alpha_power;

% ---- Data cast (GPU/CPU) ----
if CONFIG.use_gpu
    try
        gpuDevice;
        dataCast = 'gpuArray-single';
        fprintf('       Compute: GPU\n');
    catch
        dataCast = 'single';
        fprintf('       Compute: CPU (GPU unavailable)\n');
    end
else
    dataCast = 'single';
    fprintf('       Compute: CPU\n');
end

inputArgs = {'Smooth', false, ...
             'PMLInside', false, ...
             'PMLSize', CONFIG.pml_size, ...
             'DataCast', dataCast, ...
             'PlotSim', false};

% ==== FORWARD SIMULATION ====
fprintf('       Running forward simulation...\n');

source_fwd    = struct();
source_fwd.p0 = p0;

try
    fwd_tic = tic;
    sensorData = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
    fwd_time = toc(fwd_tic);
    fprintf('       Forward complete (%.1f s). Sensor data: [%d x %d]\n', ...
        fwd_time, size(sensorData, 1), size(sensorData, 2));
catch ME
    fprintf('[ERROR] Forward simulation failed: %s\n', ME.message);
    return;
end

% ==== TIME REVERSAL RECONSTRUCTION ====
num_tr_iter = CONFIG.num_time_reversal_iter;
fprintf('       Running time reversal (%d iteration(s))...\n', num_tr_iter);

reconPressure = zeros(gridSize);

try
    tr_tic = tic;

    for tr_iter = 1:num_tr_iter
        if num_tr_iter > 1
            fprintf('         TR iteration %d/%d...\n', tr_iter, num_tr_iter);
        end

        % Time-reversed source on sensor locations
        source_tr        = struct();
        source_tr.p_mask = sensor.mask;
        source_tr.p      = fliplr(sensorData);
        source_tr.p_mode = 'dirichlet';

        % Record final pressure over entire grid
        sensor_tr        = struct();
        sensor_tr.mask   = ones(Nx, Ny, Nz);
        sensor_tr.record = {'p_final'};

        p0_recon = kspaceFirstOrder3D(kgrid, kmedium, source_tr, sensor_tr, inputArgs{:});

        % Extract 3D pressure field
        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
        else
            reconPressure = p0_recon;
        end

        % Positivity constraint
        reconPressure = max(reconPressure, 0);

        % Iterative TR: compute residual for next iteration
        if tr_iter < num_tr_iter
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, ...
                source_resid, sensor, inputArgs{:});
            sensorData = sensorData + (sensorData - sensorDataRecon);
        end
    end

    tr_time = toc(tr_tic);
    fprintf('       Time reversal complete (%.1f s).\n', tr_time);
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] Time reversal failed: %s\n', ME.message);
    return;
end

%% ========================= PRESSURE -> DOSE ==============================

conversionFactor = medium.gruneisen .* medium.density;
conversionFactor(conversionFactor == 0) = 1;  % prevent div-by-zero

reconDosePerPulse = reconPressure ./ conversionFactor;
recon_dose = reconDosePerPulse * num_pulses;

fprintf('\n========= RESULTS =========\n');
fprintf('  Original dose:       [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));
fprintf('  Reconstructed dose:  [%.6f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));

% Quick error metrics within dose region
dose_region = doseGrid > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
    rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
    fprintf('  Mean abs error (>1%% dose): %.6f Gy\n', mean(abs_error));
    fprintf('  Mean rel error (>1%% dose): %.2f%%\n', mean(rel_error));
    fprintf('  Max  rel error (>1%% dose): %.2f%%\n', max(rel_error));
end

total_time = fwd_time + tr_time;
fprintf('  Total simulation time: %.1f s (fwd: %.1f, TR: %.1f)\n', ...
    total_time, fwd_time, tr_time);
fprintf('===========================\n');

%% ========================= SAVE RESULTS =================================

if CONFIG.save_results
    results = struct();
    results.recon_dose     = recon_dose;
    results.original_dose  = doseGrid;
    results.p0             = p0;
    results.reconPressure  = reconPressure;
    results.sensor_mask    = sensor.mask;
    results.config         = CONFIG;
    results.spacing_mm     = spacing_mm;
    results.grid_size      = gridSize;
    results.fwd_time_sec   = fwd_time;
    results.tr_time_sec    = tr_time;

    save(CONFIG.output_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', CONFIG.output_file);
end

%% ========================= VISUALIZATION =================================

if CONFIG.plot_results
    plot_dose_comparison(doseGrid, recon_dose, sensor.mask, spacing_mm);
end

fprintf('\nStandalone simulation complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function medium = create_medium(sct, config)
%CREATE_MEDIUM Build acoustic medium from SCT data and tissue model config
%
%   Returns struct with density, sound_speed, alpha_coeff, alpha_power,
%   gruneisen, all on the same 3D grid as sct.cubeHU.

    HU = double(sct.cubeHU);
    gridSize = size(HU);

    method = lower(config.gruneisen_method);

    % --- Get tissue property lookup tables ---
    tables = define_tissue_tables();

    switch method
        case 'uniform'
            % Spatially constant properties everywhere
            medium.density     = ones(gridSize) * config.uniform_density;
            medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
            medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
            medium.alpha_power = config.uniform_alpha_power;
            medium.gruneisen   = ones(gridSize) * config.uniform_gruneisen;

        case {'threshold_1', 'threshold_2'}
            T = tables.(method);
            nTissues   = length(T.tissue_names);
            boundaries = T.hu_boundaries;

            % Initialize with water defaults
            medium.density     = ones(gridSize) * 1000;
            medium.sound_speed = ones(gridSize) * 1540;
            medium.alpha_coeff = ones(gridSize) * 0.5;
            medium.alpha_power = T.alpha_power(1);  % scalar (use first tissue)
            medium.gruneisen   = ones(gridSize) * 0.11;

            % Assign tissue properties by HU thresholding
            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                medium.density(mask)     = T.density(t);
                medium.sound_speed(mask) = T.sound_speed(t);
                medium.alpha_coeff(mask) = T.alpha_coeff(t);
                medium.gruneisen(mask)   = T.gruneisen(t);
            end

            % Use dominant alpha_power (most common tissue)
            % For simplicity, use the soft tissue value
            if isfield(T, 'alpha_power')
                % Find the soft tissue index
                st_idx = find(contains(lower(T.tissue_names), 'soft'), 1);
                if ~isempty(st_idx)
                    medium.alpha_power = T.alpha_power(st_idx);
                else
                    medium.alpha_power = T.alpha_power(1);
                end
            end

            fprintf('       Tissue model: %s (%d tissues)\n', method, nTissues);
            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                fprintf('         %-12s: %7d voxels (%.1f%%)\n', ...
                    T.tissue_names{t}, sum(mask(:)), ...
                    100 * sum(mask(:)) / numel(HU));
            end

        otherwise
            error('Unknown gruneisen_method: %s', method);
    end

    % --- Apply per-property uniform overrides ---
    if config.force_uniform_density
        medium.density = ones(gridSize) * config.uniform_density;
        fprintf('       [OVERRIDE] density -> %.0f kg/m^3 (uniform)\n', ...
            config.uniform_density);
    end

    if config.force_uniform_sound_speed
        medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
        fprintf('       [OVERRIDE] sound_speed -> %.0f m/s (uniform)\n', ...
            config.uniform_sound_speed);
    end

    if config.force_uniform_attenuation
        medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
        medium.alpha_power = config.uniform_alpha_power;
        fprintf('       [OVERRIDE] attenuation -> %.4f dB/MHz^%.1f/cm (uniform)\n', ...
            config.uniform_alpha_coeff, config.uniform_alpha_power);
    end

    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gridSize) * config.uniform_gruneisen;
        fprintf('       [OVERRIDE] gruneisen -> %.4f (uniform)\n', ...
            config.uniform_gruneisen);
    end

    % --- Fill regions outside body with water (if bodyMask available) ---
    if isfield(sct, 'bodyMask')
        outside_body = ~sct.bodyMask;
        if isfield(sct, 'couchMask')
            outside_body = outside_body | sct.couchMask;
        end
        medium.density(outside_body)     = 1000;   % water
        medium.sound_speed(outside_body) = 1480;   % water
        medium.alpha_coeff(outside_body) = 0.0022; % water
        medium.gruneisen(outside_body)   = 0.11;   % water
    end

    medium.grid_size = gridSize;
end


function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding
%
%   Matches the tables in ethos_master_pipeline_pseudocode.m

    % THRESHOLD_1: 9-tissue model
    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air', 'Lung', 'Fat', 'Water', 'Blood', ...
                                        'Muscle', 'SoftTissue', 'Bone', 'Metal'};
    tables.threshold_1.density       = [1.2,   400,   920,   1000,  1060, 1050, 1040, 1900, 7800];
    tables.threshold_1.sound_speed   = [343,   600,   1450,  1480,  1575, 1580, 1540, 3200, 5900];
    tables.threshold_1.alpha_coeff   = [0,     0.5,   0.48,  0.0022, 0.2, 0.5,  0.5,  10,   0];
    tables.threshold_1.alpha_power   = [1.0,   1.5,   1.5,   2.0,   1.3,  1.0,  1.1,  1.0,  1.0];
    tables.threshold_1.gruneisen     = [0,     0.5,   0.7,   0.11,  0.15, 0.2,  1.0,  0,    0];

    % THRESHOLD_2: 4-tissue model
    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water', 'Fat', 'SoftTissue', 'Bone'};
    tables.threshold_2.density       = [1000,    920,  1040,          1900];
    tables.threshold_2.sound_speed   = [1480,    1450, 1540,          3200];
    tables.threshold_2.alpha_coeff   = [0.0022,  0.48, 0.5,           10];
    tables.threshold_2.alpha_power   = [2.0,     1.5,  1.1,           1.0];
    tables.threshold_2.gruneisen     = [0.11,    0.7,  1.0,           0];
end


function sensor = place_sensor_for_field(doseMask, Nx, Ny, Nz, gantry_angle)
%PLACE_SENSOR_FOR_FIELD Create planar sensor based on dose extent and gantry
%
%   Sensor is placed on the beam exit side, determined by gantry angle:
%     IEC convention (from isocenter):
%       0   deg -> beam from anterior (+Y), sensor on posterior (-Y)
%       90  deg -> beam from right   (-X), sensor on left      (+X)
%       180 deg -> beam from posterior(-Y), sensor on anterior  (+Y)
%       270 deg -> beam from left    (+X), sensor on right     (-X)

    MARGIN = 10;
    PAD_XZ = 5;

    [xIdx, yIdx, zIdx] = ind2sub([Nx, Ny, Nz], find(doseMask));

    sensor = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    if isempty(xIdx)
        return;
    end

    xMin = max(1,  min(xIdx) - PAD_XZ);
    xMax = min(Nx, max(xIdx) + PAD_XZ);
    yMin = max(1,  min(yIdx) - PAD_XZ);
    yMax = min(Ny, max(yIdx) + PAD_XZ);
    zMin = max(1,  min(zIdx) - PAD_XZ);
    zMax = min(Nz, max(zIdx) + PAD_XZ);

    ga = mod(gantry_angle, 360);

    if (ga >= 315 || ga < 45)
        sY = max(1, min(yIdx) - MARGIN);
        sensor.mask(xMin:xMax, sY, zMin:zMax) = 1;
    elseif (ga >= 45 && ga < 135)
        sX = min(Nx, max(xIdx) + MARGIN);
        sensor.mask(sX, yMin:yMax, zMin:zMax) = 1;
    elseif (ga >= 135 && ga < 225)
        sY = min(Ny, max(yIdx) + MARGIN);
        sensor.mask(xMin:xMax, sY, zMin:zMax) = 1;
    else
        sX = max(1, min(xIdx) - MARGIN);
        sensor.mask(sX, yMin:yMax, zMin:zMax) = 1;
    end

    fprintf('       Gantry %.0f deg -> sensor: %d voxels\n', ...
        gantry_angle, sum(sensor.mask(:)));
end


function plot_dose_comparison(original, reconstructed, sensor_mask, spacing_mm)
%PLOT_DOSE_COMPARISON Visualize original vs reconstructed dose
%
%   Shows transverse, coronal, and sagittal slices at the dose centroid.

    gridSize = size(original);

    % Find dose centroid slice indices
    dose_thresh = original > 0.01 * max(original(:));
    [ix, iy, iz] = ind2sub(gridSize, find(dose_thresh));

    if isempty(ix)
        fprintf('  [plot] No significant dose to display.\n');
        return;
    end

    % Dose-weighted centroid
    dose_vals = original(dose_thresh);
    cx = round(sum(ix .* dose_vals) / sum(dose_vals));
    cy = round(sum(iy .* dose_vals) / sum(dose_vals));
    cz = round(sum(iz .* dose_vals) / sum(dose_vals));

    cx = max(1, min(gridSize(1), cx));
    cy = max(1, min(gridSize(2), cy));
    cz = max(1, min(gridSize(3), cz));

    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    % ---- Figure 1: Side-by-side comparison ----
    figure('Name', 'Dose Comparison', 'Position', [100, 100, 1400, 900]);
    sgtitle(sprintf('Dose Comparison (centroid slice: X=%d, Y=%d, Z=%d)', cx, cy, cz));

    % Transverse (axial) slice at Z = cz
    subplot(2, 3, 1);
    imagesc(squeeze(original(:, :, cz))');
    axis image; colorbar; title('Original - Transverse');
    xlabel('X'); ylabel('Y'); caxis([0, max_dose]);

    subplot(2, 3, 4);
    imagesc(squeeze(reconstructed(:, :, cz))');
    axis image; colorbar; title('Reconstructed - Transverse');
    xlabel('X'); ylabel('Y'); caxis([0, max_dose]);

    % Coronal slice at Y = cy
    subplot(2, 3, 2);
    imagesc(squeeze(original(:, cy, :))');
    axis image; colorbar; title('Original - Coronal');
    xlabel('X'); ylabel('Z'); caxis([0, max_dose]);

    subplot(2, 3, 5);
    imagesc(squeeze(reconstructed(:, cy, :))');
    axis image; colorbar; title('Reconstructed - Coronal');
    xlabel('X'); ylabel('Z'); caxis([0, max_dose]);

    % Sagittal slice at X = cx
    subplot(2, 3, 3);
    imagesc(squeeze(original(cx, :, :))');
    axis image; colorbar; title('Original - Sagittal');
    xlabel('Y'); ylabel('Z'); caxis([0, max_dose]);

    subplot(2, 3, 6);
    imagesc(squeeze(reconstructed(cx, :, :))');
    axis image; colorbar; title('Reconstructed - Sagittal');
    xlabel('Y'); ylabel('Z'); caxis([0, max_dose]);

    colormap('jet');

    % ---- Figure 2: Difference and sensor overlay ----
    figure('Name', 'Error & Sensor', 'Position', [150, 50, 1200, 500]);
    sgtitle('Reconstruction Error & Sensor Placement');

    % Difference map (transverse)
    subplot(1, 3, 1);
    diff_slice = squeeze(reconstructed(:, :, cz))' - squeeze(original(:, :, cz))';
    imagesc(diff_slice);
    axis image; colorbar; title('Difference (Recon - Orig) - Transverse');
    xlabel('X'); ylabel('Y');
    caxis([-max_dose * 0.2, max_dose * 0.2]);
    colormap(gca, 'jet');

    % Relative error (dose region only)
    subplot(1, 3, 2);
    orig_slice = squeeze(original(:, :, cz))';
    recon_slice = squeeze(reconstructed(:, :, cz))';
    rel_err_slice = zeros(size(orig_slice));
    dose_mask_slice = orig_slice > 0.1 * max_dose;
    rel_err_slice(dose_mask_slice) = 100 * abs(recon_slice(dose_mask_slice) - ...
        orig_slice(dose_mask_slice)) ./ orig_slice(dose_mask_slice);
    imagesc(rel_err_slice);
    axis image; colorbar; title('Relative Error (%) - Transverse');
    xlabel('X'); ylabel('Y');
    caxis([0, 50]);

    % Sensor overlay on dose (coronal)
    subplot(1, 3, 3);
    coronal_dose = squeeze(original(:, cy, :))';
    coronal_sensor = squeeze(sensor_mask(:, cy, :))';
    imagesc(coronal_dose); hold on;
    % Overlay sensor as red contour
    if any(coronal_sensor(:))
        contour(coronal_sensor, [0.5, 0.5], 'r', 'LineWidth', 2);
    end
    axis image; colorbar; title('Dose + Sensor (red) - Coronal');
    xlabel('X'); ylabel('Z');
    colormap(gca, 'jet');
    hold off;

    drawnow;
end
