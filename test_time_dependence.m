%% =========================================================================
%  TEST_TIME_DEPENDENCE.m
%  Compare two methods of applying Gaussian pulse time-dependence:
%    Sim 1 - Time-dependent source: p(x,t) = p0(x) * g(t)
%    Sim 2 - Instantaneous source + post-convolution: conv(sensor_data, g)
%
%  By linearity these should be identical. This script verifies that.
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';

CONFIG.dose_filename = 'total_rs_dose.mat';
CONFIG.sct_filename  = 'sct_resampled.mat';

CONFIG.dose_file_override = 'field1dose.mat';
CONFIG.sct_file_override  = '';

CONFIG.sensor_placement_method = 'full_plane_lateral';
CONFIG.sensor_x_index = 20;
CONFIG.sensor_y_index = 40;

CONFIG.gruneisen_method = 'threshold_2';

CONFIG.force_uniform_density     = false;
CONFIG.force_uniform_sound_speed = false;
CONFIG.force_uniform_attenuation = false;
CONFIG.force_uniform_gruneisen   = false;

CONFIG.uniform_density      = 1000;
CONFIG.uniform_sound_speed  = 1540;
CONFIG.uniform_alpha_coeff  = 0;
CONFIG.uniform_alpha_power  = 1.1;
CONFIG.uniform_gruneisen    = 1.0;

CONFIG.dose_per_pulse_cGy = 0.16;
CONFIG.meterset           = 140;
CONFIG.pml_size           = 10;
CONFIG.cfl_number         = 0.3;
CONFIG.use_gpu            = true;

CONFIG.downscale_factor = 1;
CONFIG.use_grid_padding = true;

% Gaussian pulse width (sigma), in seconds
CONFIG.gaussian_sigma_s = 5e-6;   % 5 µs

% Source point threshold for Sim 1: only voxels with p0 > this fraction
% of max(p0) are used as active source points (memory management)
CONFIG.source_threshold_frac = 0.01;

%% ========================= RESOLVE FILE PATHS ============================

if ~isempty(CONFIG.dose_file_override)
    dose_filepath = CONFIG.dose_file_override;
else
    processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    dose_filepath = fullfile(processed_dir, CONFIG.dose_filename);
end

if ~isempty(CONFIG.sct_file_override)
    sct_filepath = CONFIG.sct_file_override;
else
    if ~exist('processed_dir', 'var')
        processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
            CONFIG.patient_id, CONFIG.session, 'processed');
    end
    sct_filepath = fullfile(processed_dir, CONFIG.sct_filename);
end

%% ========================= PRINT CONFIGURATION ===========================

fprintf('=========================================================\n');
fprintf('  Time-Dependence Test  (Sim1 vs Sim2)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', dose_filepath);
fprintf('  SCT file:        %s\n', sct_filepath);
fprintf('  Gaussian sigma:  %.1f µs\n', CONFIG.gaussian_sigma_s * 1e6);
fprintf('=========================================================\n\n');

%% ========================= LOAD DATA =====================================

fprintf('[1/5] Loading dose data...\n');
if ~isfile(dose_filepath)
    error('Dose file not found: %s', dose_filepath);
end
dose_data   = load(dose_filepath);
dose_fields = fieldnames(dose_data);
if isfield(dose_data, 'total_rs_dose')
    doseGrid = dose_data.total_rs_dose;
elseif isfield(dose_data, 'dose_Gy')
    doseGrid = dose_data.dose_Gy;
elseif length(dose_fields) == 1
    doseGrid = dose_data.(dose_fields{1});
else
    error('Cannot auto-detect dose variable.');
end
if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
    error('Dose data must be a 3D numeric array.');
end

gridSize = size(doseGrid);
Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);
fprintf('       Grid size: [%d x %d x %d]\n', Nx, Ny, Nz);

fprintf('[2/5] Loading SCT data...\n');
if ~isfile(sct_filepath)
    error('SCT file not found: %s', sct_filepath);
end
sct_data = load(sct_filepath);
if isfield(sct_data, 'sct_resampled')
    sct = sct_data.sct_resampled;
else
    error('sct_resampled variable not found.');
end

spacing_mm = sct.spacing(:)';
dx = spacing_mm(1) / 1000;
dy = spacing_mm(2) / 1000;
dz = spacing_mm(3) / 1000;

if isfield(sct, 'bodyMask')
    doseGrid = doseGrid .* double(sct.bodyMask);
end

%% ========================= GRID DOWNSCALING ==============================

if CONFIG.downscale_factor ~= 1
    df     = CONFIG.downscale_factor;
    new_Nx = max(1, round(Nx / df));
    new_Ny = max(1, round(Ny / df));
    new_Nz = max(1, round(Nz / df));

    fprintf('[DS]  Downscaling: [%d x %d x %d] -> [%d x %d x %d]\n', ...
        Nx, Ny, Nz, new_Nx, new_Ny, new_Nz);

    doseGrid   = max(imresize3(doseGrid, [new_Nx, new_Ny, new_Nz]), 0);
    sct.cubeHU = imresize3(sct.cubeHU, [new_Nx, new_Ny, new_Nz]);

    if isfield(sct, 'bodyMask')
        sct.bodyMask = imresize3(single(sct.bodyMask), ...
            [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
    end
    if isfield(sct, 'couchMask')
        sct.couchMask = imresize3(single(sct.couchMask), ...
            [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
    end
    if isfield(sct, 'cubeDensity')
        sct.cubeDensity = imresize3(sct.cubeDensity, [new_Nx, new_Ny, new_Nz]);
    end

    spacing_mm = spacing_mm .* ([Nx, Ny, Nz] ./ [new_Nx, new_Ny, new_Nz]);
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;
    sct.spacing = spacing_mm;

    Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
    gridSize = [Nx, Ny, Nz];
end

%% ========================= ACOUSTIC MEDIUM & p0 =========================

fprintf('[3/5] Creating acoustic medium and initial pressure...\n');
medium = create_medium(sct, CONFIG);

meterset       = CONFIG.meterset;
num_pulses     = ceil(meterset / CONFIG.dose_per_pulse_cGy);
dose_per_pulse = doseGrid / num_pulses;

p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
p0 = smooth(p0);

fprintf('       p0 range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));

if max(p0(:)) == 0
    error('Initial pressure is zero everywhere. Aborting.');
end

%% ========================= GRID PADDING ==================================

Nx_orig = Nx; Ny_orig = Ny; Nz_orig = Nz;
medium_orig = medium;

if CONFIG.use_grid_padding
    Nx_pad = find_optimal_kwave_size(Nx, CONFIG.pml_size);
    Ny_pad = find_optimal_kwave_size(Ny, CONFIG.pml_size);
    Nz_pad = find_optimal_kwave_size(Nz, CONFIG.pml_size);
else
    Nx_pad = Nx; Ny_pad = Ny; Nz_pad = Nz;
end

if ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz])
    fprintf('[PAD] Padding grid: [%d %d %d] -> [%d %d %d]\n', ...
        Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);

    density_pad    = ones(Nx_pad, Ny_pad, Nz_pad) * 1000;
    soundSpeed_pad = ones(Nx_pad, Ny_pad, Nz_pad) * 1540;
    alphaCoeff_pad = zeros(Nx_pad, Ny_pad, Nz_pad);

    density_pad(1:Nx, 1:Ny, 1:Nz)    = medium.density;
    soundSpeed_pad(1:Nx, 1:Ny, 1:Nz) = medium.sound_speed;
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_pad(1:Nx, 1:Ny, 1:Nz) = medium.alpha_coeff;
    else
        alphaCoeff_pad(:) = medium.alpha_coeff;
    end

    medium.density     = density_pad;
    medium.sound_speed = soundSpeed_pad;
    medium.alpha_coeff = alphaCoeff_pad;

    p0_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
    p0_pad(1:Nx, 1:Ny, 1:Nz) = p0;
    p0 = p0_pad;

    Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
    gridSize = [Nx, Ny, Nz];
end

%% ========================= SENSOR PLACEMENT ==============================

sensor      = struct();
sensor.mask = zeros(Nx, Ny, Nz);

switch CONFIG.sensor_placement_method
    case 'full_plane_anterior'
        sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
    case 'full_plane_lateral'
        sensor.mask(:, CONFIG.sensor_y_index, :) = 1;
    case 'spherical'
        sph_radius  = floor(min([Nx, Ny, Nz]) / 2) - CONFIG.pml_size;
        sensor.mask = makeSphere(Nx, Ny, Nz, sph_radius);
    otherwise
        error('Unknown sensor_placement_method: "%s"', CONFIG.sensor_placement_method);
end
fprintf('       Sensor: %d active points\n', sum(sensor.mask(:)));

%% ========================= k-WAVE GRID SETUP =============================

fprintf('[4/5] Setting up k-Wave grid and Gaussian pulse...\n');

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

maxC = max(medium.sound_speed(:));
minC = min(medium.sound_speed(medium.sound_speed > 0));
dt   = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;

gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
simTime  = 2.5 * gridDiag / minC;
Nt       = ceil(simTime / dt);

kgrid.dt = dt;
kgrid.Nt = Nt;

fprintf('       dt = %.2e s,  Nt = %d,  T_sim = %.2e s\n', dt, Nt, simTime);

kmedium             = struct();
kmedium.density     = medium.density;
kmedium.sound_speed = medium.sound_speed;
kmedium.alpha_coeff = 0 * medium.alpha_coeff;
kmedium.alpha_power = 0;

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

inputArgs = {'Smooth', false, 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
             'DataCast', dataCast, 'PlotSim', false};

%% ========================= BUILD GAUSSIAN PULSE ==========================

sigma    = CONFIG.gaussian_sigma_s;          % 5 µs in seconds
t_vec    = (0:Nt-1) * dt;                   % time axis
t_center = 4 * sigma;                       % center Gaussian 4-sigma in

gauss_signal = exp(-0.5 * ((t_vec - t_center) / sigma).^2);
% Normalise so that sum(gauss)*dt = 1 (matches energy of instantaneous p0 injection)
gauss_signal = gauss_signal / (sum(gauss_signal) * dt);

fprintf('       Gaussian: sigma=%.1f µs, peak at t=%.1f µs, integral=%.4f\n', ...
    sigma*1e6, t_center*1e6, sum(gauss_signal)*dt);

%% ========================= SIMULATION 1: TIME-DEPENDENT SOURCE ===========
%  Source: p(x,t) = p0(x) * g(t)
%  Each nonzero p0 voxel is an active source point emitting a Gaussian pulse
%  with amplitude proportional to its p0 value.
%
%  Memory note: source.p is [Nsource x Nt]. Significant-voxel thresholding
%  keeps this tractable. Linearity is preserved by the threshold fraction.

fprintf('\n[5/5] Simulation 1: Time-dependent source...\n');

p0_thresh   = CONFIG.source_threshold_frac * max(p0(:));
source_mask = p0 > p0_thresh;
Nsource     = sum(source_mask(:));

fprintf('       Source points: %d  (threshold = %.1f%% of max p0)\n', ...
    Nsource, CONFIG.source_threshold_frac * 100);
fprintf('       Source matrix size: %.1f MB (single)\n', ...
    Nsource * Nt * 4 / 1e6);

p0_vals = single(p0(source_mask));          % [Nsource x 1]
g_row   = single(gauss_signal(:)');         % [1 x Nt]

% Build [Nsource x Nt] source matrix: outer product
source_p_matrix = p0_vals * g_row;          % [Nsource x Nt]

source_sim1         = struct();
source_sim1.p_mask  = source_mask;
source_sim1.p       = source_p_matrix;
source_sim1.p_mode  = 'additive';

try
    tic_sim1     = tic;
    sensorData1  = kspaceFirstOrder3D(kgrid, kmedium, source_sim1, sensor, inputArgs{:});
    time_sim1    = toc(tic_sim1);
    if isa(sensorData1, 'gpuArray'), sensorData1 = gather(sensorData1); end
    fprintf('       Sim 1 complete (%.1f s). Data: [%d x %d]\n', ...
        time_sim1, size(sensorData1, 1), size(sensorData1, 2));
catch ME
    fprintf('[ERROR] Sim 1 failed: %s\n', ME.message);
    return;
end

clear source_p_matrix;   % free memory before Sim 2

%% ========================= SIMULATION 2: INSTANTANEOUS + POST-CONV =======
%  Source: standard source.p0, then convolve each sensor row with g(t)

fprintf('\n       Simulation 2: Instantaneous source + post-convolution...\n');

source_sim2    = struct();
source_sim2.p0 = p0;

try
    tic_sim2     = tic;
    sensorData2_raw = kspaceFirstOrder3D(kgrid, kmedium, source_sim2, sensor, inputArgs{:});
    time_sim2       = toc(tic_sim2);
    if isa(sensorData2_raw, 'gpuArray'), sensorData2_raw = gather(sensorData2_raw); end
    fprintf('       Sim 2 forward complete (%.1f s).\n', time_sim2);
catch ME
    fprintf('[ERROR] Sim 2 failed: %s\n', ME.message);
    return;
end

% Convolve each sensor row with the Gaussian pulse, keep same length
fprintf('       Convolving Sim 2 sensor data with Gaussian...\n');
[Nsen, Nt_data] = size(sensorData2_raw);
sensorData2     = zeros(Nsen, Nt_data, 'like', sensorData2_raw);

for i = 1:Nsen
    convolved = conv(double(sensorData2_raw(i,:)), gauss_signal * dt, 'full');
    sensorData2(i,:) = convolved(1:Nt_data);
end

%% ========================= RESULTS & PLOTTING ============================

fprintf('\n===== Results =====\n');
absDiff = abs(double(sensorData1) - double(sensorData2));
fprintf('  Max |Sim1 - Sim2|:  %.4e Pa\n', max(absDiff(:)));
fprintf('  Mean |Sim1 - Sim2|: %.4e Pa\n', mean(absDiff(:)));
fprintf('  Max |Sim1|:         %.4e Pa\n', max(abs(sensorData1(:))));
fprintf('  Relative error:     %.4f%%\n', ...
    100 * max(absDiff(:)) / (max(abs(sensorData1(:))) + eps));

% Pick a representative sensor point (closest to peak dose)
[~, dose_max_lin] = max(doseGrid(:));
[cx, ~, cz] = ind2sub([Nx_orig, Ny_orig, Nz_orig], dose_max_lin);

switch CONFIG.sensor_placement_method
    case 'full_plane_lateral'
        % sensor is XZ plane; index maps as (ix-1)*Nz_active + iz
        [sx, sz] = find(squeeze(sensor.mask(:, CONFIG.sensor_y_index, :)));
        dists    = sqrt((sx - cx).^2 + (sz - cz).^2);
        [~, idx_rep] = min(dists);
        rep_label = sprintf('Sensor at x=%d, z=%d', sx(idx_rep), sz(idx_rep));
    case 'full_plane_anterior'
        [sy, sz] = find(squeeze(sensor.mask(CONFIG.sensor_x_index, :, :)));
        dists    = sqrt((sz - cz).^2);
        [~, idx_rep] = min(dists);
        rep_label = sprintf('Sensor at y=%d, z=%d', sy(idx_rep), sz(idx_rep));
    otherwise
        idx_rep   = ceil(Nsen / 2);
        rep_label = sprintf('Sensor point %d', idx_rep);
end

t_us = t_vec * 1e6;   % time axis in µs

figure('Name', 'Time Dependence Test', 'Color', 'w', ...
    'NumberTitle', 'off', 'Position', [50, 100, 1400, 420]);

sgtitle(sprintf('Time-Dependence Test  |  %s  |  Gaussian σ = %.0f µs', ...
    rep_label, sigma*1e6), 'FontWeight', 'bold', 'FontSize', 11);

% --- Panel 1: Sim 1 (time-dependent source) ---
ax1 = subplot(1, 3, 1);
plot(ax1, t_us, double(sensorData1(idx_rep, :)), 'b-', 'LineWidth', 1.4);
xlabel(ax1, 'Time (µs)');
ylabel(ax1, 'Pressure (Pa)');
title(ax1, 'Sim 1: Time-Dependent Source', 'FontWeight', 'bold');
grid(ax1, 'on');
xlim(ax1, [0, t_us(end)]);

% --- Panel 2: Sim 2 (instantaneous + post-conv) ---
ax2 = subplot(1, 3, 2);
plot(ax2, t_us, double(sensorData2(idx_rep, :)), 'r-', 'LineWidth', 1.4);
xlabel(ax2, 'Time (µs)');
ylabel(ax2, 'Pressure (Pa)');
title(ax2, 'Sim 2: Instantaneous + Post-Conv', 'FontWeight', 'bold');
grid(ax2, 'on');
xlim(ax2, [0, t_us(end)]);

% --- Panel 3: Absolute difference ---
ax3 = subplot(1, 3, 3);
plot(ax3, t_us, absDiff(idx_rep, :), 'k-', 'LineWidth', 1.4);
xlabel(ax3, 'Time (µs)');
ylabel(ax3, '|Sim1 - Sim2| (Pa)');
title(ax3, 'Absolute Difference', 'FontWeight', 'bold');
grid(ax3, 'on');
xlim(ax3, [0, t_us(end)]);

drawnow;

fprintf('\nDone.\n');

%% ========================= LOCAL FUNCTIONS ===============================

function medium = create_medium(sct, config)
    HU       = double(sct.cubeHU);
    gridSize = size(HU);
    tables   = define_tissue_tables();

    switch lower(config.gruneisen_method)
        case 'uniform'
            medium.density     = ones(gridSize) * config.uniform_density;
            medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
            medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
            medium.alpha_power = config.uniform_alpha_power;
            medium.gruneisen   = ones(gridSize) * config.uniform_gruneisen;

        case {'threshold_1', 'threshold_2'}
            T        = tables.(config.gruneisen_method);
            nTissues = length(T.tissue_names);
            bounds   = T.hu_boundaries;

            medium.density     = ones(gridSize) * 1000;
            medium.sound_speed = ones(gridSize) * 1540;
            medium.alpha_coeff = ones(gridSize) * 0.5;
            medium.alpha_power = T.alpha_power(1);
            medium.gruneisen   = ones(gridSize) * 0.11;

            for t = 1:nTissues
                mask = (HU >= bounds(t)) & (HU < bounds(t+1));
                medium.density(mask)     = T.density(t);
                medium.sound_speed(mask) = T.sound_speed(t);
                medium.alpha_coeff(mask) = T.alpha_coeff(t);
                medium.gruneisen(mask)   = T.gruneisen(t);
            end
            st_idx = find(contains(lower(T.tissue_names), 'soft'), 1);
            if ~isempty(st_idx)
                medium.alpha_power = T.alpha_power(st_idx);
            end

        otherwise
            error('Unknown gruneisen_method: %s', config.gruneisen_method);
    end

    if config.force_uniform_density
        medium.density = ones(gridSize) * config.uniform_density;
    end
    if config.force_uniform_sound_speed
        medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
    end
    if config.force_uniform_attenuation
        medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
        medium.alpha_power = config.uniform_alpha_power;
    end
    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gridSize) * config.uniform_gruneisen;
    end

    if isfield(sct, 'bodyMask')
        outside = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask')
            outside = outside | logical(sct.couchMask);
        end
        medium.density(outside)     = config.uniform_density;
        medium.sound_speed(outside) = config.uniform_sound_speed;
        medium.alpha_coeff(outside) = config.uniform_alpha_coeff;
        medium.gruneisen(outside)   = config.uniform_gruneisen;
    end

    medium.grid_size = gridSize;
end


function tables = define_tissue_tables()
    tables.threshold_1 = struct();
    tables.threshold_1.hu_boundaries = [-1000, -900, -500, -200, -50, 13, 50, 80, 300, 3000, Inf];
    tables.threshold_1.tissue_names  = {'Air','Lung','Fat','Water','Blood','Muscle','SoftTissue','Bone','Metal'};
    tables.threshold_1.density       = [1.2,  400,  920, 1000, 1060, 1050, 1040, 1900, 7800];
    tables.threshold_1.sound_speed   = [343,  600, 1450, 1480, 1575, 1580, 1540, 3200, 5900];
    tables.threshold_1.alpha_coeff   = [0,   0.5, 0.48, 0.0022, 0.2, 0.5, 0.5,  10,   0];
    tables.threshold_1.alpha_power   = [1.0, 1.5, 1.5,  2.0,   1.3, 1.0, 1.1,  1.0,  1.0];
    tables.threshold_1.gruneisen     = [0,   0.5, 0.7,  0.11,  0.15, 0.2, 1.0,  0,    0];

    tables.threshold_2 = struct();
    tables.threshold_2.hu_boundaries = [-1000, -200, -50, 100, Inf];
    tables.threshold_2.tissue_names  = {'Water','Fat','SoftTissue','Bone'};
    tables.threshold_2.density       = [1000,   920, 1040,         1900];
    tables.threshold_2.sound_speed   = [1480,  1450, 1540,         3200];
    tables.threshold_2.alpha_coeff   = [0.0022, 0.48, 0.5,         10];
    tables.threshold_2.alpha_power   = [2.0,    1.5,  1.1,         1.0];
    tables.threshold_2.gruneisen     = [0.11,   0.7,  1.0,         1.0];
end
