%% =========================================================================
%  RUN_STANDALONE_SIMULATION.m
%  Standalone k-Wave Photoacoustic Forward + Time-Reversal Simulation
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';

CONFIG.dose_filename = 'dose_1194203_Session_1_reference_B9_96.mat';
CONFIG.sct_filename  = 'sct_resampled.mat';

CONFIG.dose_file_override = '';
CONFIG.sct_file_override  = '';

CONFIG.sensor_placement_method = 'determine_sensor_mask';
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

CONFIG.dose_per_pulse_cGy     = 0.16;
CONFIG.meterset               = 140;
CONFIG.pml_size               = 10;
CONFIG.cfl_number             = 0.3;
CONFIG.use_gpu                = true;
CONFIG.correction_factor           = 1.9;
CONFIG.use_pressure_scale_correction = true;   % divide max(p0) / max(recon_pressure) before dose conversion

% --- Reconstruction method ---
%   'tr'  : iterative time-reversal (k-Wave back-propagation)
%   'DAS' : Delay-And-Sum back-projection (homogeneous c, non-iterative)
CONFIG.reconstruction_method = 'tr';

CONFIG.num_time_reversal_iter = 30;
CONFIG.convergence_tol        = 1e-3;

% --- Pulse Convolution / Noise / Deconvolution ---
% Mimics a finite transducer impulse response applied to forward sensor data.
% Set convolution_kernel to 0 to disable the entire block.
CONFIG.convolution_kernel  = 4e-6;   % Gaussian sigma in seconds (4 us)
CONFIG.conv_noise_level    = 0.01;   % Noise amplitude as fraction of peak sensor signal
CONFIG.conv_deconv_lambda  = 1e-4;   % Wiener regularization for deconvolution

CONFIG.use_psf_correction      = false;
CONFIG.regularization_lambda   = 0.05;

CONFIG.downscale_factor = 1;
CONFIG.use_grid_padding = true;

CONFIG.save_results = true;
CONFIG.output_file  = 'standalone_recon_results.mat';
CONFIG.plot_results = true;

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
fprintf('  Standalone k-Wave Photoacoustic Simulation  (v4.1)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', dose_filepath);
fprintf('  SCT file:        %s\n', sct_filepath);
fprintf('  Sensor:          %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  Recon method:    %s\n', CONFIG.reconstruction_method);
fprintf('  TR iterations:   %d (tol: %.1e)\n', CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);
fprintf('  GPU:             %s\n', mat2str(CONFIG.use_gpu));
if CONFIG.downscale_factor ~= 1
    fprintf('  Downscale factor: %g\n', CONFIG.downscale_factor);
end
fprintf('=========================================================\n\n');

%% ========================= LOAD DATA ====================================

fprintf('[1/7] Loading dose data...\n');
if ~isfile(dose_filepath)
    error('Dose file not found: %s', dose_filepath);
end
dose_data = load(dose_filepath);

dose_fields = fieldnames(dose_data);

% --- Auto-detect step15_process_doses output formats ---
if isfield(dose_data, 'field_dose')
    % Individual field dose file from step15_process_doses.
    % dose_Gy may be stored as sparse 2D: reshape(full(dose_Gy), dose_dims)
    fd = dose_data.field_dose;
    if ~isfield(fd, 'dose_Gy')
        error('field_dose struct missing dose_Gy field.');
    end
    if (isfield(fd, 'is_sparse') && fd.is_sparse) || issparse(fd.dose_Gy)
        if ~isfield(fd, 'dose_dims')
            error('field_dose.dose_dims missing — cannot reconstruct sparse dose.');
        end
        doseGrid = reshape(full(fd.dose_Gy), fd.dose_dims);
        fprintf('       Loaded: field_dose.dose_Gy (sparse -> [%d x %d x %d])\n', fd.dose_dims);
    else
        doseGrid = double(fd.dose_Gy);
        fprintf('       Loaded: field_dose.dose_Gy (dense)\n');
    end
    % Pull embedded metadata: override CONFIG only when value is non-trivial
    if isfield(fd, 'spacing') && ~isempty(fd.spacing)
        step15_spacing_mm = fd.spacing(:)';
        fprintf('       Spacing from file:  [%.3f %.3f %.3f] mm\n', step15_spacing_mm);
    end
    if isfield(fd, 'meterset') && ~isempty(fd.meterset) && fd.meterset > 0
        if CONFIG.meterset ~= fd.meterset
            fprintf('       [INFO] Overriding CONFIG.meterset: %.2f -> %.2f MU\n', ...
                CONFIG.meterset, fd.meterset);
            CONFIG.meterset = fd.meterset;
        end
    end
    if isfield(fd, 'gantry_angle')
        fprintf('       Gantry angle: %.1f deg\n', fd.gantry_angle);
    end

elseif isfield(dose_data, 'total_rs_dose_sparse')
    % Total dose file from step15_process_doses (sparse format).
    if ~isfield(dose_data, 'total_rs_dose_dims')
        error('total_rs_dose_dims missing — cannot reconstruct sparse total dose.');
    end
    doseGrid = reshape(full(dose_data.total_rs_dose_sparse), dose_data.total_rs_dose_dims);
    fprintf('       Loaded: total_rs_dose_sparse (reconstructed [%d x %d x %d])\n', ...
        dose_data.total_rs_dose_dims);

elseif isfield(dose_data, 'total_rs_dose')
    doseGrid = dose_data.total_rs_dose;
    fprintf('       Loaded variable: total_rs_dose\n');
elseif isfield(dose_data, 'dose_Gy')
    doseGrid = dose_data.dose_Gy;
    fprintf('       Loaded variable: dose_Gy\n');
elseif length(dose_fields) == 1
    doseGrid = dose_data.(dose_fields{1});
    fprintf('       Loaded variable: %s\n', dose_fields{1});
else
    error('Cannot auto-detect dose variable. Fields found: %s', strjoin(dose_fields, ', '));
end
doseGrid = double(doseGrid);

if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
    error('Dose data must be a 3D numeric array.');
end

gridSize = size(doseGrid);
Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);
fprintf('       Grid size: [%d x %d x %d]\n', Nx, Ny, Nz);
fprintf('       Dose range: [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));

fprintf('[2/7] Loading SCT data...\n');
if ~isfile(sct_filepath)
    error('SCT file not found: %s', sct_filepath);
end
sct_data = load(sct_filepath);
if isfield(sct_data, 'sct_resampled')
    sct = sct_data.sct_resampled;
else
    error('sct_resampled variable not found in %s', sct_filepath);
end

if ~isfield(sct, 'cubeHU')
    error('sct_resampled missing required field: cubeHU');
end

% Spacing: prefer SCT field; fall back to spacing embedded in step15 field_dose
if isfield(sct, 'spacing') && ~isempty(sct.spacing)
    spacing_mm = sct.spacing(:)';
elseif exist('step15_spacing_mm', 'var')
    spacing_mm = step15_spacing_mm;
    fprintf('       [INFO] Using spacing from field_dose file: [%.3f %.3f %.3f] mm\n', spacing_mm);
else
    error('sct_resampled missing required field: spacing');
end
dx = spacing_mm(1) / 1000;
dy = spacing_mm(2) / 1000;
dz = spacing_mm(3) / 1000;

sctSize = size(sct.cubeHU);
if ~isequal(gridSize, sctSize)
    error(['Dose grid [%d %d %d] does not match SCT grid [%d %d %d].\n' ...
           'Ensure sct_resampled.mat was produced by the same step15 run as the dose file.'], ...
        Nx, Ny, Nz, sctSize(1), sctSize(2), sctSize(3));
end

fprintf('       Spacing: [%.2f, %.2f, %.2f] mm\n', spacing_mm);
fprintf('       HU range: [%.0f, %.0f]\n', min(sct.cubeHU(:)), max(sct.cubeHU(:)));

% Mask to body
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
        sct.bodyMask = imresize3(single(sct.bodyMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
    end
    if isfield(sct, 'couchMask')
        sct.couchMask = imresize3(single(sct.couchMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
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

%% ========================= CREATE ACOUSTIC MEDIUM ========================

fprintf('[3/7] Creating acoustic medium (method: %s)...\n', CONFIG.gruneisen_method);
medium = create_medium(sct, CONFIG);

fprintf('       Density range:     [%.0f, %.0f] kg/m^3\n', min(medium.density(:)), max(medium.density(:)));
fprintf('       Sound speed range: [%.0f, %.0f] m/s\n',    min(medium.sound_speed(:)), max(medium.sound_speed(:)));
fprintf('       Gruneisen range:   [%.4f, %.4f]\n',        min(medium.gruneisen(:)), max(medium.gruneisen(:)));

%% ========================= INITIAL PRESSURE p0 ==========================

% Apply body mask to dose before p0 calculation
if isfield(sct, 'bodyMask')
    doseGrid = doseGrid .* double(sct.bodyMask);
end

fprintf('[4/7] Computing initial pressure...\n');

meterset       = CONFIG.meterset;
num_pulses     = ceil(meterset / CONFIG.dose_per_pulse_cGy);
dose_per_pulse = doseGrid / num_pulses;

p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
p0 = smooth(p0);

fprintf('       Meterset: %.2f MU -> %d pulses\n', meterset, num_pulses);
fprintf('       p0 range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));

doseThreshold = 0.01 * max(doseGrid(:));
doseMask      = doseGrid > doseThreshold;

if ~any(doseMask(:)) || max(p0(:)) == 0
    warning('No significant dose or zero initial pressure. Aborting.');
    return;
end

%% ========================= OPTIMAL GRID PADDING ==========================

fprintf('[PAD] Computing FFT-optimal padded dimensions...\n');

Nx_orig = Nx; Ny_orig = Ny; Nz_orig = Nz;
gridSize_orig = gridSize;
medium_orig   = medium;

if CONFIG.use_grid_padding
    Nx_pad = find_optimal_kwave_size(Nx, CONFIG.pml_size);
    Ny_pad = find_optimal_kwave_size(Ny, CONFIG.pml_size);
    Nz_pad = find_optimal_kwave_size(Nz, CONFIG.pml_size);
else
    Nx_pad = Nx; Ny_pad = Ny; Nz_pad = Nz;
end

did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
if did_pad
    fprintf('[PAD] Padding grid: [%d %d %d] -> [%d %d %d]\n', Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);

    density_pad    = ones(Nx_pad, Ny_pad, Nz_pad) * 1000;
    soundSpeed_pad = ones(Nx_pad, Ny_pad, Nz_pad) * 1540;
    alphaCoeff_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
    gruneisen_pad  = zeros(Nx_pad, Ny_pad, Nz_pad);

    density_pad(1:Nx, 1:Ny, 1:Nz)    = medium.density;
    soundSpeed_pad(1:Nx, 1:Ny, 1:Nz) = medium.sound_speed;
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_pad(1:Nx, 1:Ny, 1:Nz) = medium.alpha_coeff;
    else
        alphaCoeff_pad(:) = medium.alpha_coeff;
    end
    gruneisen_pad(1:Nx, 1:Ny, 1:Nz) = medium.gruneisen;

    medium.density     = density_pad;
    medium.sound_speed = soundSpeed_pad;
    medium.alpha_coeff = alphaCoeff_pad;
    medium.gruneisen   = gruneisen_pad;

    p0_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
    p0_pad(1:Nx, 1:Ny, 1:Nz) = p0;
    p0 = p0_pad;

    Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
    gridSize = [Nx, Ny, Nz];
else
    fprintf('[PAD] Grid [%d %d %d] already FFT-optimal.\n', Nx, Ny, Nz);
end

%% ========================= SENSOR PLACEMENT ==============================

fprintf('[5/7] Placing sensor (method: %s)...\n', CONFIG.sensor_placement_method);

sensor      = struct();
sensor.mask = zeros(Nx, Ny, Nz);

switch CONFIG.sensor_placement_method
    case 'full_plane_anterior'
        sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
        fprintf('       Sensor: YZ plane at x = %d\n', CONFIG.sensor_x_index);
    case 'full_plane_lateral'
        sensor.mask(:, CONFIG.sensor_y_index, :) = 1;
        fprintf('       Sensor: XZ plane at y = %d\n', CONFIG.sensor_y_index);
    case 'spherical'
        sph_radius  = floor(min([Nx, Ny, Nz]) / 2) - CONFIG.pml_size;
        sensor.mask = makeSphere(Nx, Ny, Nz, sph_radius);
        fprintf('       Sensor: spherical, radius %d voxels\n', sph_radius);
    case 'determine_sensor_mask'
        % Automatic placement via determine_sensor_mask: places a flat anterior
        % sensor avoiding beam exclusion zones, closest to the dose centroid.
        sct_for_sensor = sct;
        if ~isfield(sct_for_sensor, 'couchMask')
            sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
        end
        if ~isfield(sct_for_sensor, 'origin')
            sct_for_sensor.origin = [0, 0, 0];
        end
        sct_for_sensor.spacing = spacing_mm;

        field_dose_for_sensor = struct();
        field_dose_for_sensor.dose_Gy     = doseGrid;
        field_dose_for_sensor.gantry_angle = 0;
        field_dose_for_sensor.origin      = sct_for_sensor.origin;
        field_dose_for_sensor.spacing     = spacing_mm;
        field_dose_for_sensor.dimensions  = [Nx_orig, Ny_orig, Nz_orig];

        beam_meta = [];
        if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
            beam_meta = CONFIG.beam_metadata;
        end

        [sensor_mask_orig, ~] = determine_sensor_mask( ...
            sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);

        % Embed into the current (possibly padded) grid
        m1 = min(Nx, size(sensor_mask_orig, 1));
        m2 = min(Ny, size(sensor_mask_orig, 2));
        m3 = min(Nz, size(sensor_mask_orig, 3));
        sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
        fprintf('       Sensor: determine_sensor_mask — %d active points\n', sum(sensor_mask_orig(:)));
    case 'fixed_anterior'
        % Deterministic placement: anterior, inferior to beam field,
        % laterally centered on isocenter X.
        % Requires sct.bodyMask and either CONFIG.beam_metadata or the
        % dose already loaded as CONFIG.total_dose / CONFIG.total_dose_file.
        fixed_struct = struct();
        if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
            fixed_struct.beam_metadata = CONFIG.beam_metadata;
        end
        % Pass the pre-loaded dose (pre-padding, original grid size)
        fixed_struct.total_dose = doseGrid;
        % sct is at the current (downscaled if applicable, unpadded) grid
        [sensor_mask_orig, ~] = determine_sensor_placement_fixed(CONFIG, sct, fixed_struct);
        % Embed into the current (possibly padded) grid
        m1 = min(Nx, size(sensor_mask_orig, 1));
        m2 = min(Ny, size(sensor_mask_orig, 2));
        m3 = min(Nz, size(sensor_mask_orig, 3));
        sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
        fprintf('       Sensor: fixed_anterior — %d active points\n', sum(sensor_mask_orig(:)));
    otherwise
        error('Unknown sensor_placement_method: "%s"', CONFIG.sensor_placement_method);
end

numSensorPts = sum(sensor.mask(:));
fprintf('       Sensor: %d active points\n', numSensorPts);

if numSensorPts == 0
    warning('Sensor mask is empty. Aborting.');
    return;
end

%% ========================= 3D SENSOR vs DOSE VISUALIZATION ==============
%  Displayed before simulation starts  rotate to inspect geometry.

if CONFIG.plot_results
    sensor_vis    = logical(sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig));
    dose_mask_vis = double(doseGrid) >= 0.10 * max(double(doseGrid(:)));
    plot_sensor_dose_planes(dose_mask_vis, sensor_vis, spacing_mm, medium_orig.density, CONFIG);
    fprintf('       [Sensor vs dose mask visualization displayed]\n');
    drawnow;
end

%% ========================= k-WAVE GRID & MEDIUM SETUP ===================

fprintf('[6/7] Setting up k-Wave grid...\n');

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

maxC = max(medium.sound_speed(:));
minC = min(medium.sound_speed(medium.sound_speed > 0));
dt   = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;

gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
simTime  = 2.5 * gridDiag / minC;
Nt       = ceil(simTime / dt);

kgrid.dt = dt;
kgrid.Nt = Nt;

fprintf('       dt = %.2e s, Nt = %d, T_sim = %.2e s\n', dt, Nt, simTime);

kmedium             = struct();
kmedium.density     = medium.density;
kmedium.sound_speed = medium.sound_speed;
kmedium.alpha_coeff = medium.alpha_coeff;
kmedium.alpha_power = 1.1;

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

%% ========================= PSF CORRECTION FILTER =========================

psf_filter = [];
if CONFIG.use_psf_correction
    fprintf('\n[PSF] Computing PSF correction filter...\n');
    try
        psf_filter = get_psf(doseGrid, sct, medium_orig, CONFIG);
        fprintf('[PSF] Filter ready.\n');
    catch ME
        warning('get_psf failed: %s. Proceeding without PSF correction.', ME.message);
        psf_filter = [];
    end
end

%% ========================= FORWARD SIMULATION ============================

fprintf('\n[7/7] Running k-Wave forward simulation...\n');

source_fwd    = struct();
source_fwd.p0 = p0;

try
    fwd_tic    = tic;
    sensorData = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
    fwd_time   = toc(fwd_tic);
    fprintf('       Forward complete (%.1f s). Sensor data: [%d x %d]\n', ...
        fwd_time, size(sensorData, 1), size(sensorData, 2));
catch ME
    fprintf('[ERROR] Forward simulation failed: %s\n', ME.message);
    return;
end

sensorData_measured = sensorData;

%% ========================= PULSE CONVOLUTION / NOISE / DECONVOLUTION =====
%  Mimics a finite transducer impulse response:
%    1. Convolve each sensor time series with a Gaussian pulse kernel
%    2. Add white Gaussian noise at the specified level
%    3. Wiener-deconvolve with the same kernel to recover the broadband signal
%  Time reversal then proceeds on the deconvolved data as normal.

if CONFIG.convolution_kernel > 0
    conv_kernel_sigma  = CONFIG.convolution_kernel;
    conv_noise_level   = CONFIG.conv_noise_level;
    conv_deconv_lambda = CONFIG.conv_deconv_lambda;

    fprintf('       Pulse model: sigma=%.1f us, noise=%.1f%%, lambda=%.1e\n', ...
        conv_kernel_sigma * 1e6, conv_noise_level * 100, conv_deconv_lambda);

    % Build normalized Gaussian kernel in time (truncated at ±4 sigma)
    sigma_samples = conv_kernel_sigma / dt;
    kernel_half   = ceil(4 * sigma_samples);
    t_kernel      = (-kernel_half : kernel_half)';
    gauss_kernel  = exp(-t_kernel.^2 / (2 * sigma_samples^2));
    gauss_kernel  = gauss_kernel / sum(gauss_kernel);   % unit-sum normalization

    % Move to CPU for FFT operations
    sensorData_cpu = double(gather(sensorData));
    Nt_data        = size(sensorData_cpu, 2);

    % Kernel transfer function (zero-padded to signal length)
    H       = fft(gauss_kernel, Nt_data).';   % row vector [1 x Nt_data]
    H_conj  = conj(H);
    H_power = abs(H).^2;

    % 1. Convolve
    sensorData_conv = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));

    % 2. Add noise
    noise_amp        = conv_noise_level * max(abs(sensorData_conv(:)));
    sensorData_noisy = sensorData_conv + noise_amp * randn(size(sensorData_conv));

    % 3. Wiener deconvolution
    sensorData_deconv = real(ifft( ...
        fft(sensorData_noisy, [], 2) .* H_conj ./ (H_power + conv_deconv_lambda), ...
        [], 2));

    % Replace both working and reference sensor data with processed result
    sensorData          = single(sensorData_deconv);
    sensorData_measured = single(sensorData_deconv);

    fprintf('       Pulse model complete. Noise amp: %.3e Pa\n', noise_amp);
else
    fprintf('       Pulse convolution disabled.\n');
end

%% ========================= RECONSTRUCTION ================================

switch lower(CONFIG.reconstruction_method)
case 'tr'

fprintf('       Running iterative time reversal (%d iterations, tol=%.1e)...\n', ...
    CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);

reconPressure      = zeros(gridSize);
reconPressure_prev = zeros(gridSize);

% Convergence tracking
conv_max_pressure = zeros(CONFIG.num_time_reversal_iter, 1);
conv_rel_change   = nan(CONFIG.num_time_reversal_iter, 1);
num_iters_done    = 0;

% ---- Set up live TR figure ----
% Axial slice through the max-dose voxel (in original, pre-pad coords).
[~, dose_max_idx] = max(doseGrid(:));
[cx_live, cy_live, cz_live] = ind2sub([Nx_orig, Ny_orig, Nz_orig], dose_max_idx);

if CONFIG.plot_results
    fig_live = figure('Name', 'Live TR Reconstruction', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1060, 440]);

    % Panel 1  initial p0 (axial slice, fixed reference)
    ax_p0 = subplot(1, 3, 1);
    p0_orig_slice = squeeze(p0(1:Nx_orig, 1:Ny_orig, cz_live))';
    imagesc(ax_p0, p0_orig_slice);
    axis(ax_p0, 'xy'); axis(ax_p0, 'image');
    colormap(ax_p0, 'hot'); colorbar(ax_p0);
    clim_p0 = [0, max(p0_orig_slice(:)) + eps];
    caxis(ax_p0, clim_p0);
    xlabel(ax_p0, 'X (voxel)'); ylabel(ax_p0, 'Y (voxel)');
    title(ax_p0, sprintf('Initial p_0   (Z=%d)', cz_live), 'FontWeight', 'bold');

    % Panel 2  current reconstructed p0 (updates each iteration)
    ax_recon = subplot(1, 3, 2);
    hImg_recon = imagesc(ax_recon, zeros(Ny_orig, Nx_orig));
    axis(ax_recon, 'xy'); axis(ax_recon, 'image');
    colormap(ax_recon, 'hot'); colorbar(ax_recon);
    xlabel(ax_recon, 'X (voxel)'); ylabel(ax_recon, 'Y (voxel)');
    title(ax_recon, 'Reconstructed p_0   (iter 0)', 'FontWeight', 'bold');

    % Panel 3  live max-pressure convergence
    ax_conv = subplot(1, 3, 3);
    hLine_max = plot(ax_conv, NaN, NaN, 'b-o', 'LineWidth', 1.6, ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.2, 0.4, 1.0]);
    xlabel(ax_conv, 'TR Iteration');
    ylabel(ax_conv, 'Max Reconstructed p_0 (Pa)');
    title(ax_conv, 'Convergence (live)', 'FontWeight', 'bold');
    grid(ax_conv, 'on');
    xlim(ax_conv, [0.5, CONFIG.num_time_reversal_iter + 0.5]);

    sgtitle(fig_live, sprintf( ...
        'Live TR Reconstruction  Axial Z=%d   |   Patient %s', ...
        cz_live, CONFIG.patient_id), 'FontWeight', 'bold', 'FontSize', 11);
    drawnow;
end

% ---- TR iteration loop ----
try
    tr_total_tic = tic;

    for tr_iter = 1:CONFIG.num_time_reversal_iter

        fprintf('       --- TR Iteration %d/%d ---\n', tr_iter, CONFIG.num_time_reversal_iter);

        source_tr        = struct();
        source_tr.p_mask = sensor.mask;
        source_tr.p      = fliplr(sensorData);
        source_tr.p_mode = 'dirichlet';

        sensor_tr        = struct();
        sensor_tr.mask   = ones(Nx, Ny, Nz);
        sensor_tr.record = {'p_final'};

        p0_recon = kspaceFirstOrder3D(kgrid, kmedium, source_tr, sensor_tr, inputArgs{:});

        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
        else
            reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
        end

        reconPressure = max(reconPressure, 0);

        % Record convergence metrics
        conv_max_pressure(tr_iter) = max(reconPressure(:));
        num_iters_done = tr_iter;

        fprintf('       Max pressure: %.4e Pa\n', conv_max_pressure(tr_iter));

        % Convergence check (from iteration 2 onward)
        converged = false;
        if tr_iter > 1
            norm_prev = norm(reconPressure_prev(:));
            if norm_prev > 0
                rel_change = norm(reconPressure(:) - reconPressure_prev(:)) / norm_prev;
            else
                rel_change = Inf;
            end
            conv_rel_change(tr_iter) = rel_change;
            fprintf('       Rel change: %.4e\n', rel_change);
            if rel_change < CONFIG.convergence_tol
                fprintf('       *** Converged at iteration %d ***\n', tr_iter);
                converged = true;
            end
        end

        reconPressure_prev = reconPressure;

        % ---- Update live figure ----
        if CONFIG.plot_results && ishandle(fig_live)
            recon_slice_crop = squeeze( ...
                reconPressure(1:Nx_orig, 1:Ny_orig, cz_live))';
            recon_slice_crop = gather(recon_slice_crop);
            set(hImg_recon, 'CData', recon_slice_crop);
            caxis(ax_recon, [0, max(recon_slice_crop(:)) + eps]);
            if converged
                title(ax_recon, ...
                    sprintf('Reconstructed p_0   (iter %d  CONVERGED)', tr_iter), ...
                    'FontWeight', 'bold', 'Color', [0, 0.55, 0]);
            else
                title(ax_recon, ...
                    sprintf('Reconstructed p_0   (iter %d / %d)', ...
                    tr_iter, CONFIG.num_time_reversal_iter), 'FontWeight', 'bold');
            end
            set(hLine_max, 'XData', 1:tr_iter, ...
                'YData', conv_max_pressure(1:tr_iter));
            drawnow;
        end

        if converged
            break;
        end

        % Residual correction for next iteration
        if tr_iter < CONFIG.num_time_reversal_iter
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, source_resid, sensor, inputArgs{:});
            sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
        end
    end

    reconPressure = gather(reconPressure) * CONFIG.correction_factor;

    tr_time = toc(tr_total_tic);
    fprintf('       Time reversal complete (%.1f s).\n', tr_time);
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] Time reversal failed: %s\n', ME.message);
    return;
end

case 'das'

fprintf('       Running Delay-And-Sum reconstruction...\n');

try
    das_tic = tic;

    % 1. Sensor positions in physical coords on the padded grid
    sensor_lin = find(sensor.mask);
    [sx_i, sy_i, sz_i] = ind2sub([Nx, Ny, Nz], sensor_lin);
    sx_pos = (double(sx_i) - 1) * dx;
    sy_pos = (double(sy_i) - 1) * dy;
    sz_pos = (double(sz_i) - 1) * dz;
    nSens  = numel(sensor_lin);

    % 2. Voxel grid coords (padded grid)
    [Xg, Yg, Zg] = ndgrid((0:Nx-1)*dx, (0:Ny-1)*dy, (0:Nz-1)*dz);
    Xg = single(Xg);  Yg = single(Yg);  Zg = single(Zg);

    % 3. Effective sound speed: mean of nonzero medium speeds (homogeneous assumption)
    nonzero_c = medium.sound_speed(medium.sound_speed > 0);
    c_das = double(mean(nonzero_c, 'all'));

    % 4. Move sensor data to CPU single (avoids GPU memory blowup with full grid)
    sd     = single(gather(sensorData));
    Nt_sd  = size(sd, 2);

    reconPressure = zeros(Nx, Ny, Nz, 'single');
    dx_min        = single(min([dx, dy, dz]));

    % 5. Loop over sensors (memory-bounded; Nx*Ny*Nz*Nsens would be huge)
    for s = 1:nSens
        d_s   = sqrt((Xg - single(sx_pos(s))).^2 + ...
                     (Yg - single(sy_pos(s))).^2 + ...
                     (Zg - single(sz_pos(s))).^2);
        idx_f = d_s / single(c_das * dt) + 1;     % fractional sample index
        idx0  = floor(idx_f);
        frac  = idx_f - idx0;
        valid = idx0 >= 1 & idx0 < Nt_sd;

        samp = zeros(size(idx0), 'single');
        i0   = idx0(valid);
        v0   = sd(s, i0);
        v1   = sd(s, i0 + 1);
        samp(valid) = (1 - frac(valid)) .* v0(:) + frac(valid) .* v1(:);

        % 1/r weighting compensates spherical spreading
        w = single(1) ./ max(d_s, dx_min);
        reconPressure = reconPressure + w .* samp;

        if mod(s, max(1, round(nSens/10))) == 0
            fprintf('         DAS sensor %d/%d\n', s, nSens);
        end
    end

    reconPressure = max(double(reconPressure), 0) * CONFIG.correction_factor;

    % Synthesize TR-style outputs so downstream plot/save code is reused
    conv_max_pressure = max(reconPressure(:));
    conv_rel_change   = NaN;
    num_iters_done    = 1;
    tr_time           = toc(das_tic);

    fprintf('       DAS complete (%.1f s).\n', tr_time);
    fprintf('       c_DAS = %.1f m/s, sensors used: %d\n', c_das, nSens);
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] DAS reconstruction failed: %s\n', ME.message);
    return;
end

otherwise
    error('Unknown reconstruction_method: "%s" (use ''tr'' or ''DAS'')', ...
        CONFIG.reconstruction_method);
end

%% ========================= CROP TO ORIGINAL SIZE =========================

if did_pad
    fprintf('\n[CROP] Restoring dimensions: [%d %d %d] -> [%d %d %d]\n', ...
        Nx, Ny, Nz, Nx_orig, Ny_orig, Nz_orig);
    reconPressure = reconPressure(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    Nx = Nx_orig; Ny = Ny_orig; Nz = Nz_orig;
    gridSize    = gridSize_orig;
    medium      = medium_orig;
    sensor.mask = sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
end

%% ========================= PSF CORRECTION ================================

psf_applied = false;
if ~isempty(psf_filter) && isstruct(psf_filter) && isfield(psf_filter, 'F') ...
        && ~isempty(psf_filter.F)
    fprintf('       Applying PSF correction...\n');
    P_field       = fftn(reconPressure);
    corrected     = real(ifftn(P_field .* psf_filter.F));
    reconPressure = max(corrected, 0);
    psf_applied   = true;
    fprintf('       Corrected pressure range: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));
end

%% ========================= PRESSURE SCALE CORRECTION ====================

if CONFIG.use_pressure_scale_correction
    p0_max_orig = max(p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig), [], 'all');
    recon_max   = max(reconPressure(:));
    if recon_max > 0
        pressure_scale_cf = p0_max_orig / recon_max;
        reconPressure     = reconPressure * pressure_scale_cf;
        fprintf('       Pressure scale correction: %.4f  (p0_max = %.3e  /  recon_max = %.3e)\n', ...
            pressure_scale_cf, p0_max_orig, recon_max);
    else
        pressure_scale_cf = 1;
        warning('Reconstructed pressure max is zero; skipping pressure scale correction.');
    end
else
    pressure_scale_cf = 1;
end

%% ========================= PRESSURE -> DOSE ==============================

fprintf('\n[Post] Converting pressure to dose...\n');

conversionFactor = medium.gruneisen .* medium.density;
conversionFactor(conversionFactor == 0) = 1;

reconDosePerPulse = reconPressure ./ conversionFactor;

body_mask_plot = ones(gridSize);
if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), gridSize)
    body_mask_plot = double(sct.bodyMask);
end
recon_dose = reconDosePerPulse * num_pulses .* double(doseMask) .* body_mask_plot;

fprintf('       Reconstructed dose: [%.4f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));

%% Crop p0 to original size (if padded)
if did_pad
    p0 = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
end

%% ========================= RESULTS SUMMARY ===============================

fprintf('\n========= RESULTS =========\n');
fprintf('  Original dose:      [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));
fprintf('  Reconstructed dose: [%.6f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));
fprintf('  PSF correction:     %s\n', mat2str(psf_applied));

dose_region = doseGrid > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
    rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
    fprintf('  Mean abs err: %.6f Gy\n', mean(abs_error));
    fprintf('  Mean rel err: %.2f%%\n',  mean(rel_error));
    fprintf('  Max  rel err: %.2f%%\n',  max(rel_error));
end
fprintf('===========================\n');

%% ========================= GAMMA ANALYSIS ================================

gamma_results = struct();

if exist('CalcGamma', 'file') == 2

    fprintf('\n[Gamma] Running gamma analysis...\n');

    ref_struct.start = [0, 0, 0];
    ref_struct.width = spacing_mm;
    ref_struct.data  = double(doseGrid);

    tgt_struct.start = [0, 0, 0];
    tgt_struct.width = spacing_mm;
    tgt_struct.data  = double(recon_dose);

    low_dose_cutoff = 0.10 * max(doseGrid(:));
    gamma_eval_mask = doseGrid >= low_dose_cutoff;

    gamma_criteria = {10, 10, '10%/10mm'; 5, 5, '5%/5mm'; 3, 3, '3%/3mm'};
    gamma_maps     = cell(size(gamma_criteria, 1), 1);
    pass_rates     = zeros(size(gamma_criteria, 1), 1);

    for gc = 1:size(gamma_criteria, 1)
        pct_crit = gamma_criteria{gc, 1};
        dta_crit = gamma_criteria{gc, 2};
        lbl      = gamma_criteria{gc, 3};

        fprintf('  [%s] Computing...', lbl);
        try
            gmap = CalcGamma(ref_struct, tgt_struct, pct_crit, dta_crit, ...
                'local', 0, 'limit', dta_crit * 2, 'restrict', 1);
            gamma_maps{gc} = gmap;
            eval_vals       = gmap(gamma_eval_mask);
            pass_rate       = 100 * mean(eval_vals <= 1);
            pass_rates(gc)  = pass_rate;
            fprintf('  Pass rate: %.2f%%\n', pass_rate);
        catch ME
            warning('Gamma [%s] failed: %s', lbl, ME.message);
            gamma_maps{gc} = [];
            pass_rates(gc) = NaN;
            fprintf('  FAILED\n');
        end
    end

    gamma_results.maps        = gamma_maps;
    gamma_results.pass_rates  = pass_rates;
    gamma_results.criteria    = gamma_criteria;
    gamma_results.cutoff_Gy   = low_dose_cutoff;
    gamma_results.eval_mask   = gamma_eval_mask;

    fprintf('\n  ------ Gamma Pass Rates (10%% low-dose cutoff) ------\n');
    for gc = 1:size(gamma_criteria, 1)
        if isnan(pass_rates(gc))
            fprintf('  %-12s  FAILED\n', gamma_criteria{gc, 3});
        else
            fprintf('  %-12s  %.2f%%\n', gamma_criteria{gc, 3}, pass_rates(gc));
        end
    end
else
    warning('CalcGamma not found. Skipping gamma analysis.');
    gamma_results = [];
end

%% ========================= SAVE RESULTS =================================

if CONFIG.save_results
    % Build filename from configuration parameters
    output_fname = sprintf('standalone_results_%s_%s_%s_%d.mat', ...
        CONFIG.reconstruction_method, ...
        CONFIG.sensor_placement_method, CONFIG.gruneisen_method, ...
        CONFIG.num_time_reversal_iter);

    results.recon_dose    = recon_dose;
    results.original_dose = doseGrid;
    results.p0            = p0;
    results.reconPressure = reconPressure;
    results.sensor_mask   = sensor.mask;
    results.config        = CONFIG;
    results.spacing_mm    = spacing_mm;
    results.grid_size     = gridSize;
    results.fwd_time_sec  = fwd_time;
    results.tr_time_sec   = tr_time;
    results.psf_applied   = psf_applied;
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    save(output_fname, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', output_fname);
end

%% ========================= POST-SIMULATION VISUALIZATION ================

if CONFIG.plot_results

    % Figure 1  2x3 dose comparison (original top, recon bottom)
    %            Coronal | Sagittal | Axial
    %            Isocenter = max-dose voxel; sensor contour in red
    plot_dose_panels(doseGrid, recon_dose, sensor.mask, medium_orig.density, spacing_mm, ...
        'Dose Comparison: Original vs Reconstructed');

    % Figure 2  p0 convergence (max pressure + relative change) — TR only
    if strcmpi(CONFIG.reconstruction_method, 'tr')
        p0_max_for_plot = max(p0(:));
        plot_convergence_history(conv_max_pressure, conv_rel_change, ...
            num_iters_done, CONFIG.convergence_tol, p0_max_for_plot);
    end

    % Figure 3  Axial gamma (3 criteria) + absolute error
    if ~isempty(gamma_results) && isfield(gamma_results, 'maps')
        [~, max_dose_idx] = max(doseGrid(:));
        [~, ~, cz_gamma]  = ind2sub(gridSize, max_dose_idx);
        plot_gamma_and_error_axial(gamma_results, doseGrid, recon_dose, ...
            sensor.mask, cz_gamma);
    end

end

fprintf('\nStandalone simulation complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, density, config)
%PLOT_SENSOR_DOSE_PLANES  1x3 anatomical view of sensor geometry vs dose mask.
%  Shows three orthogonal projections (coronal, sagittal, axial).
%  CT density is rendered as a grayscale background (mean-projection).
%  Dose mask (dose >= 10% max) drawn as a filled semi-transparent blue region.
%  Sensor drawn as a solid red line/region — computed via max-projection so it
%  always appears regardless of which slice the dose centroid falls on.
%  This replaces the interactive 3-D isosurface view.

    [Nx3, Ny3, Nz3] = size(dose_mask);

    x_ax = (1:Nx3) * spacing_mm(1);
    y_ax = (1:Ny3) * spacing_mm(2);
    z_ax = (1:Nz3) * spacing_mm(3);

    % Max-projections of dose mask and sensor (so features always appear)
    dose_cor  = squeeze(any(dose_mask,   2));   % XZ  (Nx x Nz)
    dose_sag  = squeeze(any(dose_mask,   1));   % YZ  (Ny x Nz)
    dose_axi  = squeeze(any(dose_mask,   3));   % XY  (Nx x Ny)

    sens_cor  = squeeze(any(sensor_mask, 2));   % XZ
    sens_sag  = squeeze(any(sensor_mask, 1));   % YZ
    sens_axi  = squeeze(any(sensor_mask, 3));   % XY

    % CT density background: mean-projection (DRR-like anatomical context).
    % Soft-tissue window/level (kg/m^3) matches plot_dose_panels.
    have_density = ~isempty(density) && isequal(size(density), size(dose_mask));
    wl_center = 1050; wl_width = 350;
    wl_min    = wl_center - wl_width / 2;
    if have_density
        ct_projs = {
            squeeze(mean(double(density), 2))',   % Coronal  XZ
            squeeze(mean(double(density), 1))',   % Sagittal YZ
            squeeze(mean(double(density), 3))'    % Axial    XY
        };
    end

    figure('Name', 'Sensor Placement vs Dose Mask', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 1300, 420]);
    sgtitle(sprintf('Sensor Placement vs Dose Mask  (\\geq10%% max)   |   Sensor: %s', ...
        config.sensor_placement_method), 'FontWeight', 'bold', 'FontSize', 11);

    view_data = {
        dose_cor', sens_cor', x_ax, z_ax, 'X (mm)', 'Z (mm)', 'Coronal  (mean-proj along Y)';
        dose_sag', sens_sag', y_ax, z_ax, 'Y (mm)', 'Z (mm)', 'Sagittal  (mean-proj along X)';
        dose_axi', sens_axi', x_ax, y_ax, 'X (mm)', 'Y (mm)', 'Axial  (mean-proj along Z)';
    };

    dose_color   = [0.20, 0.50, 0.90];   % blue
    sensor_color = [0.90, 0.10, 0.10];   % red

    for col = 1:3
        ax  = subplot(1, 3, col);
        d2d = double(view_data{col, 1});
        s2d = double(view_data{col, 2});
        xv  = view_data{col, 3};
        yv  = view_data{col, 4};

        % Background: CT density as grayscale, or white if unavailable
        if have_density
            dn     = (ct_projs{col} - wl_min) / wl_width;
            dn     = max(0, min(1, dn));
            bg_rgb = repmat(dn, [1, 1, 3]);
        else
            bg_rgb = ones([size(d2d), 3]);   % white fallback
        end
        image(ax, xv, yv, bg_rgb);
        hold(ax, 'on');

        % Dose mask: filled blue region
        if any(d2d(:))
            dose_rgb = repmat(reshape(dose_color, [1,1,3]), [size(d2d), 1]);
            h_dose   = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = d2d * 0.45;
        end

        % Sensor: filled red region (appears as line for planar sensors)
        if any(s2d(:))
            sens_rgb = repmat(reshape(sensor_color, [1,1,3]), [size(s2d), 1]);
            h_sens   = image(ax, xv, yv, sens_rgb);
            h_sens.AlphaData = s2d * 0.85;
        end

        hold(ax, 'off');
        axis(ax, 'xy'); axis(ax, 'image');
        xlabel(ax, view_data{col, 5}); ylabel(ax, view_data{col, 6});
        title(ax, view_data{col, 7}, 'FontSize', 10);

        % Manual legend patches
        hold(ax, 'on');
        p1 = patch(ax, NaN, NaN, dose_color,   'FaceAlpha', 0.45, 'EdgeColor', 'none');
        p2 = patch(ax, NaN, NaN, sensor_color,  'FaceAlpha', 0.85, 'EdgeColor', 'none');
        hold(ax, 'off');
        if col == 3
            legend(ax, [p1, p2], {'Dose mask (>=10%)', 'Sensor'}, ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
        end
    end
    drawnow;
end


function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr)
%PLOT_DOSE_PANELS 2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = original dose,  Row 2 = reconstructed dose.
%  Dose (jet, semi-transparent) is overlaid on the density map (grayscale).
%  Voxels with dose < 10% of max are fully transparent (masked out).
%  Both rows share an identical dose colour range [0, max(original)]
%  so magnitudes are directly comparable.
%  Isocenter at max-dose voxel.  Sensor contour in red on every panel.
%
%  Pass density=[] to show dose only with a black background.

    gridSize = size(original);
    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end

    have_density = ~isempty(density) && isequal(size(density), gridSize);

    [~, max_idx] = max(original(:));
    [cx, cy, cz] = ind2sub(gridSize, max_idx);

    % Shared colour scale anchored to original (reference) dose
    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    x_ax = (1:gridSize(1)) * spacing_mm(1);
    y_ax = (1:gridSize(2)) * spacing_mm(2);
    z_ax = (1:gridSize(3)) * spacing_mm(3);

    % Jet LUT for manual RGB blending (avoids per-axes colormap conflict)
    cmap_jet = jet(256);

    % Soft-tissue window/level for the density background (kg/m^3).
    % W=350 / L=1050 keeps soft tissue in the mid-grey band and clips
    % bone (~1900 kg/m^3) to white, preventing bone from dominating the display.
    wl_center = 1050;          % kg/m^3  (soft tissue ~1000-1080)
    wl_width  = 350;           % kg/m^3
    wl_min    = wl_center - wl_width / 2;   % 875  kg/m^3
    wl_max    = wl_center + wl_width / 2;   % 1225 kg/m^3

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [50, 50, 1380, 700]);
    sgtitle(sprintf('%s\nIsocenter (max dose): X=%d  Y=%d  Z=%d voxel  |  Dose clim [0, %.4f] Gy', ...
        titleStr, cx, cy, cz, max_dose), 'FontWeight', 'bold', 'FontSize', 11);

    row_labels = {'Original', 'Reconstructed'};
    doses      = {original, recon};

    for row = 1:2
        d   = doses{row};
        lbl = row_labels{row};

        % Slice data for the three views
        dose_slices   = { squeeze(d(:, cy, :))',  squeeze(d(cx, :, :))',  squeeze(d(:, :, cz))' };
        sensor_slices = { squeeze(sensor_mask(:, cy, :))', ...
                          squeeze(sensor_mask(cx, :, :))', ...
                          squeeze(sensor_mask(:, :, cz))' };
        xvecs  = { x_ax, y_ax, x_ax };
        yvecs  = { z_ax, z_ax, y_ax };
        xlbls  = { 'X (mm)', 'Y (mm)', 'X (mm)' };
        ylbls  = { 'Z (mm)', 'Z (mm)', 'Y (mm)' };
        tsuffs = { sprintf('Coronal (Y=%d)', cy), ...
                   sprintf('Sagittal (X=%d)', cx), ...
                   sprintf('Axial (Z=%d)', cz) };
        if have_density
            density_slices = { squeeze(density(:, cy, :))', ...
                               squeeze(density(cx, :, :))', ...
                               squeeze(density(:, :, cz))' };
        end

        for col = 1:3
            ax         = subplot(2, 3, (row-1)*3 + col);
            dose_slice = double(dose_slices{col});
            xv         = xvecs{col};
            yv         = yvecs{col};

            % --- Background: density as grayscale RGB ---
            if have_density
                dn     = (density_slices{col} - wl_min) / wl_width;
                dn     = max(0, min(1, dn));   % clip to [0,1]
                bg_rgb = repmat(dn, [1, 1, 3]);
            else
                bg_rgb = zeros([size(dose_slice), 3]);   % black
            end
            image(ax, xv, yv, bg_rgb);
            hold(ax, 'on');

            % --- Foreground: dose as jet RGB with alpha mask ---
            %   mask: dose >= 10% of max_dose  (shared threshold)
            norm_d   = max(0, min(1, dose_slice / max_dose));
            idx      = max(1, min(256, round(norm_d * 255) + 1));
            sz       = size(dose_slice);
            dose_rgb = reshape(cmap_jet(idx(:), :), [sz, 3]);

            above      = dose_slice >= 0.10 * max_dose;
            ramp       = 0.45 + 0.35 * min(1, ...
                (dose_slice - 0.10*max_dose) / max(0.40*max_dose, 1e-12));
            dose_alpha = above .* ramp;

            h_dose           = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = dose_alpha;

            % --- Sensor contour ---
            s = sensor_slices{col};
            if any(s(:))
                contour(ax, xv, yv, s, [0.5, 0.5], 'r-', 'LineWidth', 2);
            end
            hold(ax, 'off');

            axis(ax, 'xy'); axis(ax, 'image');

            % Colorbar: attach jet LUT with shared dose clim
            colormap(ax, cmap_jet);
            caxis(ax, [0, max_dose]);
            cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';

            xlabel(ax, xlbls{col}); ylabel(ax, ylbls{col});
            title(ax, sprintf('%s  %s', lbl, tsuffs{col}));
        end
    end
    drawnow;
end

function plot_convergence_history(conv_max_pressure, conv_rel_change, num_iters, tol, p0_max_orig)
%PLOT_CONVERGENCE_HISTORY p0 convergence over TR iterations.
%  Left axis:  max reconstructed p0 per iteration (blue).
%              Dashed black line = max of original p0 distribution (reference).
%  Right axis: relative change between iterations (red, log-scale).

    iters   = 1:num_iters;
    p_vals  = conv_max_pressure(iters);
    rc_vals = conv_rel_change(iters);

    figure('Name', 'p0 Convergence', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [150, 520, 720, 390]);

    yyaxis left;
    plot(iters, p_vals, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.2, 0.4, 1.0]);
    hold on;
    if nargin >= 5 && ~isempty(p0_max_orig) && p0_max_orig > 0
        yline(p0_max_orig, 'k--', 'LineWidth', 1.8, ...
            'Label', sprintf('p_{0}^{orig} max = %.2e Pa', p0_max_orig), ...
            'LabelHorizontalAlignment', 'left', ...
            'LabelVerticalAlignment', 'bottom', 'FontSize', 8);
    end
    hold off;
    ylabel('Max Reconstructed p_0 (Pa)', 'Color', [0.2, 0.4, 1.0]);
    top_val = max([max(p_vals), p0_max_orig]) * 1.20;
    if top_val <= 0, top_val = 1; end
    ylim([0, top_val]);

    yyaxis right;
    valid = ~isnan(rc_vals);
    if any(valid)
        semilogy(iters(valid), rc_vals(valid), 'r-s', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'MarkerFaceColor', [0.9, 0.1, 0.1]);
        hold on;
        yline(tol, 'k--', sprintf('tol = %.0e', tol), 'LineWidth', 1.2, ...
            'LabelHorizontalAlignment', 'right');
        hold off;
    end
    ylabel('Relative Change ||p_n - p_{n-1}|| / ||p_{n-1}||', ...
        'Color', [0.8, 0.1, 0.1]);

    xlabel('TR Iteration');
    title(sprintf('p_0 Convergence  (%d/%d iterations)', num_iters, numel(conv_max_pressure)), ...
        'FontWeight', 'bold');
    xlim([0.5, num_iters + 0.5]);
    grid on;
    drawnow;
end


function plot_gamma_and_error_axial(gamma_results, original, recon, sensor_mask, cz)
%PLOT_GAMMA_AND_ERROR_AXIAL 1x4 axial figure:
%  gamma 10/10 | gamma 5/5 | gamma 3/3 | absolute dose error.
%  Sensor contour in red.  Slice at cz (max-dose axial index).

    gridSize = size(original);
    cz = max(1, min(gridSize(3), cz));

    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end

    eval_mask  = gamma_results.eval_mask;
    maps       = gamma_results.maps;
    criteria   = gamma_results.criteria;
    pass_rates = gamma_results.pass_rates;
    nCrit      = size(criteria, 1);

    figure('Name', 'Gamma & Absolute Error  Axial', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [50, 300, 1400, 370]);
    sgtitle(sprintf('Axial Plane (Z = %d voxel)  Gamma Index & Absolute Error', cz), ...
        'FontWeight', 'bold', 'FontSize', 11);

    gamma_clim   = [0, 2];
    sensor_slice = squeeze(sensor_mask(:, :, cz))';

    % Gamma panels
    for g = 1:nCrit
        ax = subplot(1, 4, g);
        gmap = maps{g};
        if ~isempty(gmap)
            gmap_disp = double(gmap);
            gmap_disp(~eval_mask) = NaN;
            slice2d = squeeze(gmap_disp(:, :, cz))';
        else
            slice2d = nan(gridSize(2), gridSize(1));
        end

        imagesc(ax, slice2d, gamma_clim);
        axis(ax, 'xy'); axis(ax, 'image');
        colormap(ax, gamma_colormap_gyr());
        cb = colorbar(ax); cb.Label.String = '\gamma';
        caxis(ax, gamma_clim);

        hold(ax, 'on');
        if any(sensor_slice(:))
            contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
        end
        hold(ax, 'off');

        pr_str = 'FAILED';
        if ~isnan(pass_rates(g))
            pr_str = sprintf('%.1f%%', pass_rates(g));
        end
        xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
        title(ax, sprintf('%s\nPass rate: %s', criteria{g, 3}, pr_str));
    end

    % Absolute error panel
    ax = subplot(1, 4, 4);
    orig_slice  = squeeze(original(:, :, cz))';
    recon_slice = squeeze(recon(:, :, cz))';
    err_slice   = abs(recon_slice - orig_slice);
    max_err     = max(err_slice(:));
    if max_err == 0, max_err = 1; end

    imagesc(ax, err_slice, [0, max_err]);
    axis(ax, 'xy'); axis(ax, 'image');
    colormap(ax, 'hot');
    cb = colorbar(ax); cb.Label.String = '|Error| (Gy)';

    hold(ax, 'on');
    if any(sensor_slice(:))
        contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
    end
    hold(ax, 'off');

    xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
    title(ax, sprintf('|Recon - Original|\nMax: %.4f Gy', max_err));

    drawnow;
end


function cmap = gamma_colormap_gyr()
%GAMMA_COLORMAP_GYR Green-yellow-red colormap for gamma index.
%  Green = pass (gamma <= 1),  Red = fail (gamma > 1).
    n    = 256;
    half = round(n / 2);
    rest = n - half;
    r = [linspace(0, 1, half)'; ones(rest, 1)];
    g = [linspace(0.85, 1, half)'; linspace(1, 0, rest)'];
    b = zeros(n, 1);
    cmap = [r, g, b];
end


function medium = create_medium(sct, config)
%CREATE_MEDIUM Build acoustic medium from SCT data and tissue model config.

    HU       = double(sct.cubeHU);
    gridSize = size(HU);
    tables   = define_tissue_tables();

    switch lower(config.gruneisen_method)
        case 'uniform'
            medium.density     = ones(gridSize) * config.uniform_density;
            medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
            medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
            medium.alpha_power = 1.1;
            medium.gruneisen   = ones(gridSize) * config.uniform_gruneisen;

        case {'threshold_1', 'threshold_2'}
            T          = tables.(config.gruneisen_method);
            nTissues   = length(T.tissue_names);
            boundaries = T.hu_boundaries;

            medium.density     = ones(gridSize) * 1000;
            medium.sound_speed = ones(gridSize) * 1540;
            medium.alpha_coeff = ones(gridSize) * 0.5;
            medium.alpha_power = 1.1;
            medium.gruneisen   = ones(gridSize) * 0.11;

            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                medium.density(mask)     = T.density(t);
                medium.sound_speed(mask) = T.sound_speed(t);
                medium.alpha_coeff(mask) = T.alpha_coeff(t);
                medium.gruneisen(mask)   = T.gruneisen(t);
            end

            fprintf('       Tissue model: %s (%d tissues)\n', config.gruneisen_method, nTissues);
            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                fprintf('         %-12s: %7d voxels (%.1f%%)\n', ...
                    T.tissue_names{t}, sum(mask(:)), 100 * sum(mask(:)) / numel(HU));
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
        medium.alpha_power = 1.1;
    end
    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gridSize) * config.uniform_gruneisen;
    end

    if isfield(sct, 'bodyMask')
        outside_body = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask')
            outside_body = outside_body | logical(sct.couchMask);
        end
        medium.density(outside_body)     = config.uniform_density;
        medium.sound_speed(outside_body) = config.uniform_sound_speed;
        medium.alpha_coeff(outside_body) = config.uniform_alpha_coeff;
        medium.gruneisen(outside_body)   = config.uniform_gruneisen;
    end

    medium.grid_size = gridSize;
end


function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding.

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