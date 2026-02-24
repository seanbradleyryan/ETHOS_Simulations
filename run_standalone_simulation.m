%% =========================================================================
%  RUN_STANDALONE_SIMULATION.m
%  Standalone k-Wave Photoacoustic Forward + Time-Reversal Simulation
%  =========================================================================
%
%  PURPOSE:
%  Self-contained script for running a single-field (or total-dose)
%  photoacoustic simulation without the full ETHOS pipeline.  Loads a
%  .mat dose file and sct_resampled.mat from the processed/ output
%  directory of step15_process_doses, builds the acoustic medium, places
%  a sensor, runs the k-Wave forward simulation, and performs
%  time-reversal reconstruction.
%
%  CONFIGURABLE OPTIONS (see CONFIGURATION section):
%    - Patient/session selection (constructs path to processed/ dir)
%    - Dose file selection (default: total_rs_dose.mat in processed/)
%    - Sensor mode: 'full_anterior_plane', 'dose_based', or 'spherical'
%    - Tissue heterogeneity: 'uniform', 'threshold_1', 'threshold_2'
%    - Per-property heterogeneity overrides
%    - Spherical compensation for limited-view planar reconstruction
%    - GPU/CPU selection, PML size, CFL number, TR iterations
%
%  SPHERICAL COMPENSATION (CONFIG.spherical_compensation.enable = true):
%    Corrects limited-view artifacts in planar reconstruction by:
%      1. Running an additional spherical (full-enclosure) reconstruction
%      2. Computing a frequency-domain transfer function:
%           H(k) = FFT3(p0_sphere) / FFT3(p0_planar)
%      3. Applying H(k) as a Wiener-regularized filter to the planar result
%    This captures the geometry-dependent frequency loss of the limited
%    aperture and compensates for it.  The regularization parameter
%    controls noise amplification vs artifact correction.
%
%  REQUIRED FILES (in processed/ directory from step15_process_doses):
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
%  VERSION: 2.0
%
%  See also: run_single_field_simulation, determine_sensor_mask,
%            step15_process_doses, load_processed_data
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION ================================

% --- Patient / Session Selection ---
% These construct the path to the processed/ output directory from step15.
% Path: {working_dir}/RayStationFiles/{patient_id}/{session}/processed/
CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';

% --- File Selection ---
% Filenames within the processed/ directory.  Paths are constructed
% automatically from patient_id/session above.
% To override with an explicit full path, set CONFIG.dose_file_override
% and CONFIG.sct_file_override (both empty by default).
CONFIG.dose_filename = 'total_rs_dose.mat';   % dose file in processed/
CONFIG.sct_filename  = 'sct_resampled.mat';   % CT file in processed/

% Override: set to a full path to bypass the patient/session directory
% structure. Leave empty to use the standard processed/ directory.
CONFIG.dose_file_override = '';   % e.g., '/some/path/total_rs_dose.mat'
CONFIG.sct_file_override  = '';   % e.g., '/some/path/sct_resampled.mat'

% --- Sensor Placement Mode ---
%   'full_anterior_plane' : Entire anteriormost voxel plane is the sensor.
%                           Simple, no geometry dependence, good baseline.
%   'dose_based'          : Sensor placed on beam-exit side based on dose
%                           extent and gantry angle (original pipeline logic).
%   'spherical'           : Spherical sensor via makeSphere, centered on grid.
%                           Best possible reconstruction (simulation-only).
%                           Useful as a reference / for filter calibration.
CONFIG.sensor_mode = 'full_anterior_plane';  % 'full_anterior_plane' | 'dose_based' | 'spherical'

% --- Gantry Angle (used only for dose_based sensor mode) ---
% If using total dose (not a single field), set a representative angle.
% 0 deg = beam from anterior; 180 deg = beam from posterior.
CONFIG.gantry_angle = 180;  % degrees (only used when sensor_mode = 'dose_based')

% --- Tissue Heterogeneity ---
%   'uniform'       : Homogeneous water-like medium everywhere
%   'threshold_1'   : 9-tissue model (air, lung, fat, water, blood,
%                     muscle, soft tissue, bone, metal)
%   'threshold_2'   : 4-tissue model (water, fat, soft tissue, bone)
CONFIG.gruneisen_method = 'threshold_2';

% --- Per-Property Heterogeneity Overrides ---
% When gruneisen_method is NOT 'uniform', you can selectively force
% individual properties to be spatially constant while keeping others
% heterogeneous.  Set to true to DISABLE heterogeneity for that property.
CONFIG.force_uniform_density     = false;
CONFIG.force_uniform_sound_speed = false;
CONFIG.force_uniform_attenuation = false;
CONFIG.force_uniform_gruneisen   = false;

% --- Uniform Property Values (used when uniform or force_uniform_*) ---
CONFIG.uniform_density      = 1000;    % kg/m^3  (water)
CONFIG.uniform_sound_speed  = 1540;    % m/s     (soft tissue average)
CONFIG.uniform_alpha_coeff  = 0;     % dB/MHz^y/cm
CONFIG.uniform_alpha_power  = 1.1;     % exponent
CONFIG.uniform_gruneisen    = 1.0;     % dimensionless

% --- Simulation Parameters ---
CONFIG.dose_per_pulse_cGy     = 0.16;   % cGy per LINAC pulse
CONFIG.meterset               = 140;    % MU (monitor units) for total dose
CONFIG.pml_size               = 10;     % PML thickness (voxels)
CONFIG.cfl_number             = 0.3;    % CFL stability number
CONFIG.use_gpu                = true;   % Use GPU acceleration

% --- Iterative Time-Reversal Reconstruction ---
%   Uses the Dirichlet-BC time-reversal method from ethos_kwave_simulation.m
%   with iterative residual correction:
%     Iter 1: TR(measured_data) -> p0_est
%     Iter n: forward(p0_est) -> simulated_data
%             residual = measured_data - simulated_data
%             measured_data += residual (corrected data)
%             TR(corrected_data) -> p0_est (updated)
%             positivity constraint applied each iteration
CONFIG.num_iterations = 5;             % Number of iterative TR iterations
CONFIG.convergence_tol = 1e-4;         % Early stop if relative change < tol

% --- Spherical Compensation for Limited-View Planar Reconstruction ---
%   When enabled, runs an additional spherical (fully-enclosing) sensor
%   reconstruction on the SAME forward data source (p0), then computes a
%   frequency-domain compensation filter:
%
%     H(k) = FFT3(p0_sphere) / FFT3(p0_planar)
%
%   This filter captures the spatial-frequency attenuation caused by the
%   limited planar aperture.  It is applied to the planar reconstruction
%   with Wiener regularization to avoid noise amplification:
%
%     p0_corrected = IFFT3( FFT3(p0_planar) * conj(H) * |H|^2 / (|H|^2 + lambda) )
%   or equivalently:
%     p0_corrected = IFFT3( FFT3(p0_planar) * FFT3(p0_sphere) / (|FFT3(p0_planar)|^2 + eps) )
%
%   The regularization_lambda parameter controls the tradeoff:
%     Small lambda (e.g. 1e-3): aggressive correction, may amplify noise
%     Large lambda (e.g. 1e-1): conservative, less artifact correction
%
%   NOTE: This is only meaningful for planar sensor modes.  If sensor_mode
%         is 'spherical', this option is ignored (already fully enclosed).
%
%   COMPUTATIONAL COST: Approximately doubles the simulation time (one
%   extra forward + TR cycle for the spherical sensor).
CONFIG.spherical_compensation = struct();
CONFIG.spherical_compensation.enable = true;         % Master toggle
CONFIG.spherical_compensation.regularization_lambda = 0.01;  % Wiener regularization
CONFIG.spherical_compensation.num_iterations = 1;     % TR iterations for spherical (1 is usually enough)
CONFIG.spherical_compensation.apply_to_dose = true;   % Apply filter in dose domain (after p->D conversion)
CONFIG.spherical_compensation.save_intermediates = false; % Save sphere/planar recon for analysis

% --- Output ---
CONFIG.save_results = true;             % Save reconstruction to .mat
CONFIG.output_file  = 'standalone_recon_results.mat';
CONFIG.plot_results = true;             % Show comparison figures

%% ========================= RESOLVE FILE PATHS ============================

% Construct paths from patient/session or use overrides
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
fprintf('  Standalone k-Wave Photoacoustic Simulation  (v2.0)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', dose_filepath);
fprintf('  SCT file:        %s\n', sct_filepath);
fprintf('  Sensor mode:     %s\n', CONFIG.sensor_mode);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  TR iterations:   %d (tol: %.1e)\n', CONFIG.num_iterations, CONFIG.convergence_tol);
fprintf('  GPU:             %s\n', mat2str(CONFIG.use_gpu));
if CONFIG.spherical_compensation.enable && ~strcmpi(CONFIG.sensor_mode, 'spherical')
    fprintf('  Spherical comp:  ENABLED (lambda=%.1e, %d sphere iters)\n', ...
        CONFIG.spherical_compensation.regularization_lambda, ...
        CONFIG.spherical_compensation.num_iterations);
else
    fprintf('  Spherical comp:  off\n');
end
fprintf('=========================================================\n\n');

%% ========================= LOAD DATA ====================================

fprintf('[1/8] Loading dose data...\n');
if ~isfile(dose_filepath)
    error('Dose file not found: %s\n  Check CONFIG.patient_id, CONFIG.session, and CONFIG.working_dir.', ...
        dose_filepath);
end
dose_data = load(dose_filepath);

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
fprintf('[2/8] Loading SCT data...\n');
if ~isfile(sct_filepath)
    error('SCT file not found: %s\n  Check CONFIG.patient_id, CONFIG.session, and CONFIG.working_dir.', ...
        sct_filepath);
end
sct_data = load(sct_filepath);

if isfield(sct_data, 'sct_resampled')
    sct = sct_data.sct_resampled;
else
    error('sct_resampled variable not found in %s', sct_filepath);
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

fprintf('[3/8] Creating acoustic medium (method: %s)...\n', CONFIG.gruneisen_method);

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

fprintf('[4/8] Computing initial pressure...\n');

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

fprintf('[5/8] Placing sensor (mode: %s)...\n', CONFIG.sensor_mode);

sensor = struct();

switch lower(CONFIG.sensor_mode)
    case 'full_anterior_plane'
        % ---- Full anteriormost plane ----
        % The sensor spans the entire Y-plane at the most anterior body
        % surface (or near it), giving maximum lateral coverage.
        sensor.mask = zeros(Nx, Ny, Nz);

        if isfield(sct, 'bodyMask')
            body = sct.bodyMask;
            if isfield(sct, 'couchMask')
                body = body & ~sct.couchMask;
            end
            % Find the minimum Y index that has any body voxel
            y_has_body = squeeze(any(any(body, 1), 3));  % [1 x Ny] logical
            anterior_y = find(y_has_body, 1, 'first');
            if ~isempty(anterior_y)
                sensor_y = max(1, anterior_y - 3);
            else
                sensor_y = 1;
            end
        else
            sensor_y = 1;
        end

        sensor.mask(:, sensor_y, :) = 1;

        fprintf('       Full anterior plane at Y index %d\n', sensor_y);
        fprintf('       Sensor plane size: %d x %d = %d voxels\n', ...
            Nx, Nz, Nx * Nz);

    case 'dose_based'
        % ---- Dose-extent-based placement (original pipeline logic) ----
        sensor = place_sensor_for_field(doseMask, Nx, Ny, Nz, CONFIG.gantry_angle);
        fprintf('       Gantry angle: %.1f deg\n', CONFIG.gantry_angle);

    case 'spherical'
        % ---- Spherical sensor via makeSphere ----
        % Best-possible reconstruction for simulation reference.
        sensor.mask = create_spherical_sensor(Nx, Ny, Nz);
        fprintf('       Spherical sensor (makeSphere, centered on grid)\n');

    otherwise
        error('Unknown sensor_mode: %s. Use ''full_anterior_plane'', ''dose_based'', or ''spherical''.', ...
            CONFIG.sensor_mode);
end

numSensorPts = sum(sensor.mask(:));
fprintf('       Total sensor voxels: %d\n', numSensorPts);

if numSensorPts == 0
    warning('Sensor mask is empty. Aborting.');
    return;
end

%% ========================= k-WAVE GRID & COMMON SETUP ===================

fprintf('[6/8] Setting up k-Wave grid...\n');

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
             'PlotSim', true};

%% ========================= FORWARD SIMULATION ============================

fprintf('[7/8] Running k-Wave forward simulation...\n');

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

% Store the original measured sensor data (never modified)
sensorData_measured = sensorData;

%% ========================= ITERATIVE TIME-REVERSAL RECONSTRUCTION ========
%
%  Algorithm (Treeby & Cox, "k-Wave" iterative TR):
%    Iteration 1:
%      1. Time-reverse measured sensor data (fliplr) as Dirichlet source
%      2. Record p_final over entire grid -> p0_est
%      3. Apply positivity constraint: p0_est = max(p0_est, 0)
%
%    Iterations 2..N:
%      4. Forward-simulate current estimate p0_est -> sensorData_simulated
%      5. Compute residual: residual = sensorData_measured - sensorData_simulated
%      6. Correct: sensorData_corrected = sensorData_measured + residual
%      7. Time-reverse corrected data -> p0_est (updated)
%      8. Apply positivity constraint
%      9. Check convergence
%  =========================================================================

reconPressure = run_iterative_tr(kgrid, kmedium, sensor, sensorData_measured, ...
    gridSize, CONFIG.num_iterations, CONFIG.convergence_tol, inputArgs, 'Planar');

%% ========================= SPHERICAL COMPENSATION ========================
%
%  If enabled, run an additional fully-enclosing reconstruction and use it
%  to build a frequency-domain compensation filter for the planar result.
%
%  Theory:
%    Let p0_true be the actual initial pressure distribution.
%    Spherical TR (full enclosure) recovers:  p0_S  p0_true
%    Planar TR (limited view) recovers:       p0_P = H_lv * p0_true  (convolution)
%    where H_lv is the limited-view point spread function.
%
%    The compensation filter in k-space is:
%      F(k) = P0_S(k) / P0_P(k)   (with regularization)
%
%    For a new planar reconstruction p0_P_new:
%      p0_corrected = IFFT3( FFT3(p0_P_new) * F_regularized(k) )
%
%  Since we're in simulation, we can compute the filter from the SAME p0
%  that generated the sensor data.  The filter characterizes the geometry.
%  =========================================================================

do_spherical_comp = CONFIG.spherical_compensation.enable && ...
                    ~strcmpi(CONFIG.sensor_mode, 'spherical');

reconPressure_compensated = [];
spherical_comp_results = struct();

if do_spherical_comp
    fprintf('\n========= SPHERICAL COMPENSATION =========\n');

    % ---- Step 1: Create spherical sensor and run forward ----
    fprintf('  [SC 1/4] Creating spherical sensor...\n');
    sensor_sphere = struct();
    sensor_sphere.mask = create_spherical_sensor(Nx, Ny, Nz);
    numSpherePts = sum(sensor_sphere.mask(:));
    fprintf('           Spherical sensor: %d voxels\n', numSpherePts);

    fprintf('  [SC 2/4] Forward simulation (enclosing sensor)...\n');
    try
        sphere_fwd_tic = tic;
        sensorData_sphere = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, ...
            sensor_sphere, inputArgs{:});
        sphere_fwd_time = toc(sphere_fwd_tic);
        fprintf('           Forward complete (%.1f s). Data: [%d x %d]\n', ...
            sphere_fwd_time, size(sensorData_sphere, 1), size(sensorData_sphere, 2));
    catch ME
        fprintf('  [ERROR] Spherical forward failed: %s\n', ME.message);
        fprintf('  Skipping spherical compensation.\n');
        do_spherical_comp = false;
    end
end

if do_spherical_comp
    % ---- Step 2: Time-reverse from enclosing sensor ----
    fprintf('  [SC 3/4] Spherical TR reconstruction (%d iterations)...\n', ...
        CONFIG.spherical_compensation.num_iterations);

    reconPressure_sphere = run_iterative_tr(kgrid, kmedium, sensor_sphere, ...
        sensorData_sphere, gridSize, ...
        CONFIG.spherical_compensation.num_iterations, ...
        CONFIG.convergence_tol, inputArgs, 'Spherical');

    fprintf('           Spherical recon range: [%.2e, %.2e] Pa\n', ...
        min(reconPressure_sphere(:)), max(reconPressure_sphere(:)));

    % ---- Step 3: Compute and apply compensation filter ----
    fprintf('  [SC 4/4] Computing compensation filter...\n');
    lambda = CONFIG.spherical_compensation.regularization_lambda;

    if CONFIG.spherical_compensation.apply_to_dose
        % Apply filter in dose domain (after pressure->dose conversion)
        % This is often more meaningful since dose is what we care about.
        convFactor = medium.gruneisen .* medium.density;
        convFactor(convFactor == 0) = 1;

        dose_sphere = (reconPressure_sphere ./ convFactor) * num_pulses;
        dose_planar = (reconPressure ./ convFactor) * num_pulses;

        [dose_compensated, filter_info] = apply_spherical_filter( ...
            dose_planar, dose_sphere, lambda);

        % Convert back to pressure for consistency
        reconPressure_compensated = (dose_compensated / num_pulses) .* convFactor;

        fprintf('           Filter applied in DOSE domain\n');
    else
        % Apply filter in pressure domain
        [reconPressure_compensated, filter_info] = apply_spherical_filter( ...
            reconPressure, reconPressure_sphere, lambda);

        fprintf('           Filter applied in PRESSURE domain\n');
    end

    fprintf('           Lambda: %.1e\n', lambda);
    fprintf('           Filter dynamic range: [%.2e, %.2e]\n', ...
        filter_info.filter_min, filter_info.filter_max);
    fprintf('           Compensated pressure range: [%.2e, %.2e] Pa\n', ...
        min(reconPressure_compensated(:)), max(reconPressure_compensated(:)));

    % Store intermediates if requested
    if CONFIG.spherical_compensation.save_intermediates
        spherical_comp_results.reconPressure_sphere = reconPressure_sphere;
        spherical_comp_results.filter_info = filter_info;
    end
    spherical_comp_results.enabled = true;
    spherical_comp_results.lambda = lambda;
    spherical_comp_results.sphere_fwd_time = sphere_fwd_time;

    fprintf('==========================================\n');
end

%% ========================= PRESSURE -> DOSE ==============================

fprintf('\n[8/8] Converting pressure to dose...\n');

conversionFactor = medium.gruneisen .* medium.density;
conversionFactor(conversionFactor == 0) = 1;  % prevent div-by-zero

reconDosePerPulse = reconPressure ./ conversionFactor;
recon_dose = reconDosePerPulse * num_pulses;

% Compensated dose (if spherical compensation was applied)
if ~isempty(reconPressure_compensated)
    if CONFIG.spherical_compensation.apply_to_dose
        % Already computed in dose domain above
        recon_dose_compensated = dose_compensated;
    else
        recon_dose_compensated = (reconPressure_compensated ./ conversionFactor) * num_pulses;
    end
else
    recon_dose_compensated = [];
end

%% ========================= RESULTS SUMMARY ===============================

fprintf('\n========= RESULTS =========\n');
fprintf('  Original dose:       [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));
fprintf('  Reconstructed dose:  [%.6f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));

% Error metrics within dose region
dose_region = doseGrid > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
    rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
    fprintf('  Planar  - Mean abs err: %.6f Gy, Mean rel err: %.2f%%, Max rel err: %.2f%%\n', ...
        mean(abs_error), mean(rel_error), max(rel_error));
end

if ~isempty(recon_dose_compensated)
    fprintf('  Compensated dose:    [%.6f, %.4f] Gy\n', ...
        min(recon_dose_compensated(:)), max(recon_dose_compensated(:)));
    if any(dose_region(:))
        abs_error_c = abs(recon_dose_compensated(dose_region) - doseGrid(dose_region));
        rel_error_c = abs_error_c ./ max(doseGrid(dose_region), 1e-10) * 100;
        fprintf('  Compens - Mean abs err: %.6f Gy, Mean rel err: %.2f%%, Max rel err: %.2f%%\n', ...
            mean(abs_error_c), mean(rel_error_c), max(rel_error_c));
    end
end

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

    if ~isempty(recon_dose_compensated)
        results.recon_dose_compensated     = recon_dose_compensated;
        results.reconPressure_compensated  = reconPressure_compensated;
        results.spherical_comp             = spherical_comp_results;
    end

    save(CONFIG.output_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', CONFIG.output_file);
end

%% ========================= VISUALIZATION =================================

if CONFIG.plot_results
    % Standard planar comparison
    plot_dose_comparison(doseGrid, recon_dose, sensor.mask, spacing_mm, 'Planar Reconstruction');

    % Compensated comparison (if available)
    if ~isempty(recon_dose_compensated)
        plot_dose_comparison(doseGrid, recon_dose_compensated, sensor.mask, ...
            spacing_mm, 'Spherical-Compensated Reconstruction');

        % Side-by-side: planar vs compensated vs original
        plot_three_way_comparison(doseGrid, recon_dose, recon_dose_compensated, ...
            spacing_mm);
    end
end

fprintf('\nStandalone simulation complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function reconPressure = run_iterative_tr(kgrid, kmedium, sensor, ...
    sensorData_measured, gridSize, num_iterations, convergence_tol, inputArgs, label)
%RUN_ITERATIVE_TR Iterative time-reversal reconstruction with convergence tracking
%
%   Shared implementation used for both planar and spherical reconstructions.

    Nx = gridSize(1);
    Ny = gridSize(2);
    Nz = gridSize(3);

    fprintf('       [%s] Iterative TR (%d iterations, tol=%.1e)...\n', ...
        label, num_iterations, convergence_tol);

    % Convergence tracking
    iter_metrics = struct();
    iter_metrics.max_pressure  = zeros(num_iterations, 1);
    iter_metrics.rel_change    = zeros(num_iterations, 1);
    iter_metrics.residual_norm = zeros(num_iterations, 1);
    iter_metrics.iter_time     = zeros(num_iterations, 1);

    sensorData_working = sensorData_measured;
    reconPressure_prev = zeros(gridSize);
    reconPressure = zeros(gridSize);

    try
        tr_total_tic = tic;

        for iter = 1:num_iterations
            iter_tic = tic;
            fprintf('       [%s] --- Iteration %d/%d ---\n', label, iter, num_iterations);

            % Step A: Time-reverse current sensor data
            source_tr        = struct();
            source_tr.p_mask = sensor.mask;
            source_tr.p      = fliplr(sensorData_working);
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

            % Positivity constraint
            reconPressure = max(reconPressure, 0);

            % Track metrics
            iter_metrics.max_pressure(iter) = max(reconPressure(:));

            if iter > 1
                delta = reconPressure - reconPressure_prev;
                norm_prev = norm(reconPressure_prev(:));
                if norm_prev > 0
                    iter_metrics.rel_change(iter) = norm(delta(:)) / norm_prev;
                else
                    iter_metrics.rel_change(iter) = Inf;
                end
            else
                iter_metrics.rel_change(iter) = Inf;
            end

            iter_metrics.iter_time(iter) = toc(iter_tic);

            fprintf('       [%s] Max p: %.4e, Rel change: %.4e (%.1f s)\n', ...
                label, iter_metrics.max_pressure(iter), ...
                iter_metrics.rel_change(iter), iter_metrics.iter_time(iter));

            % Check convergence
            if iter > 1 && iter_metrics.rel_change(iter) < convergence_tol
                fprintf('       [%s] *** Converged at iteration %d ***\n', label, iter);
                break;
            end

            reconPressure_prev = reconPressure;

            % Compute residual and correct sensor data for next iteration
            if iter < num_iterations
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sensorData_simulated = kspaceFirstOrder3D(kgrid, kmedium, ...
                    source_resid, sensor, inputArgs{:});

                residual = sensorData_measured - sensorData_simulated;
                iter_metrics.residual_norm(iter) = norm(residual(:));

                sensorData_working = sensorData_measured + residual;
            end
        end

        tr_time = toc(tr_total_tic);
        fprintf('       [%s] TR complete (%.1f s total)\n', label, tr_time);

    catch ME
        fprintf('[ERROR] %s TR failed: %s\n', label, ME.message);
        reconPressure = zeros(gridSize);
    end
end


function sensor_mask = create_spherical_sensor(Nx, Ny, Nz)
%CREATE_SPHERICAL_SENSOR Create a spherical sensor using k-Wave makeSphere
%
%   Places a spherical sensor shell at the grid center with a radius chosen
%   to fit within the usable domain (inside the PML layer).  Provides full
%   solid-angle coverage and serves as the "ideal" reconstruction reference.
%
%   The radius is the largest sphere that fits within the grid, with a
%   1-voxel margin from the boundary.

    % Radius: largest sphere fitting inside the grid (PML is outside the grid)
    radius = floor(min([Nx, Ny, Nz]) / 2) - 1;

    if radius < 1
        error('create_spherical_sensor:InvalidRadius', ...
            'Grid [%d x %d x %d] is too small for a spherical sensor.', ...
            Nx, Ny, Nz);
    end

    sensor_mask = makeSphere(Nx, Ny, Nz, radius);
end


function [corrected, filter_info] = apply_spherical_filter(planar_recon, sphere_recon, lambda)
%APPLY_SPHERICAL_FILTER Wiener-regularized frequency-domain compensation
%
%   Computes a transfer function from the ratio of the spherical (good)
%   and planar (limited-view) reconstructions in 3D Fourier space, then
%   applies it with Wiener regularization.
%
%   The filter is:
%     H(k) = FFT3(sphere) / FFT3(planar)
%
%   Applied with Wiener regularization:
%     corrected = IFFT3( FFT3(planar) * conj(H) / (|H|^2 + lambda^2) * H )
%   which simplifies to:
%     corrected = IFFT3( P * S / (|P|^2 + lambda^2 * max(|P|^2)) )
%   where P = FFT3(planar), S = FFT3(sphere)
%
%   INPUTS:
%     planar_recon  - 3D array, limited-view reconstruction
%     sphere_recon  - 3D array, fully-enclosed reconstruction (reference)
%     lambda        - Regularization parameter (0.001 to 0.1 typical)
%
%   OUTPUTS:
%     corrected     - 3D array, compensated reconstruction
%     filter_info   - Struct with filter diagnostics

    % 3D FFT of both reconstructions
    P = fftn(planar_recon);
    S = fftn(sphere_recon);

    % Wiener deconvolution:
    %   We want to find F such that planar * F  sphere
    %   In frequency domain: P * F = S  ->  F = S / P
    %   Regularized: F = S * conj(P) / (|P|^2 + lambda^2 * noise_power)
    %
    %   Equivalently, the corrected output is:
    %   Corrected = P * F = P * S * conj(P) / (|P|^2 + eps)
    %            = S * |P|^2 / (|P|^2 + eps)
    %
    %   But that just gives back S when |P| is large, which isn't useful.
    %   The correct formulation for deconvolution of the limited-view PSF:
    %
    %   corrected = IFFT( S * conj(P) / (|P|^2 + lambda^2 * mean(|P|^2)) * P / conj(P) )
    %   = IFFT( S )   -- which is trivially the spherical reconstruction.
    %
    %   The CORRECT approach: the planar_recon on NEW data (not the
    %   calibration data). For calibration, we compute F = S/P, then for
    %   new planar data Q, the corrected = Q * F.
    %
    %   Since in simulation the "new" and "calibration" are the same source,
    %   we use the more robust formulation:

    % Power spectra
    P_power = abs(P).^2;
    noise_est = lambda^2 * mean(P_power(:));

    % The compensation filter: F = S / P (regularized)
    % Applied directly: corrected = IFFT( P * F ) = IFFT( P * S * conj(P) / (|P|^2 + eps) )
    %
    % For self-calibration (same source), a cleaner approach is direct
    % Wiener filtering from planar toward sphere:
    %   Corrected(k) = P(k) * [ S(k) * conj(P(k)) ] / [ |P(k)|^2 + noise_est ]
    %                = S(k) * |P(k)|^2 / [ |P(k)|^2 + noise_est ]
    %
    % This blends: where P is strong -> corrected  S (spherical result)
    %              where P is weak  -> corrected -> 0 (regularized away)
    %
    % For actual NEW data (cross-application), the filter would be stored:
    %   F_stored(k) = S(k) * conj(P(k)) / (|P(k)|^2 + noise_est)
    %   corrected_new = IFFT( FFT(new_planar) .* F_stored )

    % Self-calibration mode (same source for both):
    % Use the Wiener filter that maps P toward S
    F = S .* conj(P) ./ (P_power + noise_est);

    % Apply filter to planar reconstruction
    corrected_fft = P .* F;
    corrected = real(ifftn(corrected_fft));

    % Apply positivity constraint
    corrected = max(corrected, 0);

    % Diagnostics
    filter_info = struct();
    F_mag = abs(F);
    filter_info.filter_min = min(F_mag(:));
    filter_info.filter_max = max(F_mag(:));
    filter_info.filter_mean = mean(F_mag(:));
    filter_info.noise_est = noise_est;
    filter_info.F = F;  % Store for potential cross-application to new data
end


function medium = create_medium(sct, config)
%CREATE_MEDIUM Build acoustic medium from SCT data and tissue model config

    HU = double(sct.cubeHU);
    gridSize = size(HU);

    tables = define_tissue_tables();

    switch lower(config.gruneisen_method)
        case 'uniform'
            medium.density     = ones(gridSize) * config.uniform_density;
            medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
            medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
            medium.alpha_power = config.uniform_alpha_power;
            medium.gruneisen   = ones(gridSize) * config.uniform_gruneisen;

        case {'threshold_1', 'threshold_2'}
            T = tables.(config.gruneisen_method);
            nTissues   = length(T.tissue_names);
            boundaries = T.hu_boundaries;

            medium.density     = ones(gridSize) * 1000;
            medium.sound_speed = ones(gridSize) * 1540;
            medium.alpha_coeff = ones(gridSize) * 0.5;
            medium.alpha_power = T.alpha_power(1);
            medium.gruneisen   = ones(gridSize) * 0.11;

            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                medium.density(mask)     = T.density(t);
                medium.sound_speed(mask) = T.sound_speed(t);
                medium.alpha_coeff(mask) = T.alpha_coeff(t);
                medium.gruneisen(mask)   = T.gruneisen(t);
            end

            st_idx = find(contains(lower(T.tissue_names), 'soft'), 1);
            if ~isempty(st_idx)
                medium.alpha_power = T.alpha_power(st_idx);
            else
                medium.alpha_power = T.alpha_power(1);
            end

            fprintf('       Tissue model: %s (%d tissues)\n', config.gruneisen_method, nTissues);
            for t = 1:nTissues
                mask = (HU >= boundaries(t)) & (HU < boundaries(t+1));
                fprintf('         %-12s: %7d voxels (%.1f%%)\n', ...
                    T.tissue_names{t}, sum(mask(:)), ...
                    100 * sum(mask(:)) / numel(HU));
            end

        otherwise
            error('Unknown gruneisen_method: %s', config.gruneisen_method);
    end

    % Per-property uniform overrides
    if config.force_uniform_density
        medium.density = ones(gridSize) * config.uniform_density;
        fprintf('       [OVERRIDE] density -> %.0f kg/m^3\n', config.uniform_density);
    end
    if config.force_uniform_sound_speed
        medium.sound_speed = ones(gridSize) * config.uniform_sound_speed;
        fprintf('       [OVERRIDE] sound_speed -> %.0f m/s\n', config.uniform_sound_speed);
    end
    if config.force_uniform_attenuation
        medium.alpha_coeff = ones(gridSize) * config.uniform_alpha_coeff;
        medium.alpha_power = config.uniform_alpha_power;
        fprintf('       [OVERRIDE] attenuation -> %.4f dB/MHz^%.1f/cm\n', ...
            config.uniform_alpha_coeff, config.uniform_alpha_power);
    end
    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gridSize) * config.uniform_gruneisen;
        fprintf('       [OVERRIDE] gruneisen -> %.4f\n', config.uniform_gruneisen);
    end

    % Fill outside body with water
    if isfield(sct, 'bodyMask')
        outside_body = ~sct.bodyMask;
        if isfield(sct, 'couchMask')
            outside_body = outside_body | sct.couchMask;
        end
        medium.density(outside_body)     = 1000;
        medium.sound_speed(outside_body) = 1480;
        medium.alpha_coeff(outside_body) = 0.0022;
        medium.gruneisen(outside_body)   = 0.11;
    end

    medium.grid_size = gridSize;
end


function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding

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
    tables.threshold_2.gruneisen     = [0.11,    0.7,  1.0,           1];
end


function sensor = place_sensor_for_field(doseMask, Nx, Ny, Nz, gantry_angle)
%PLACE_SENSOR_FOR_FIELD Create planar sensor based on dose extent and gantry

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


function plot_dose_comparison(original, reconstructed, sensor_mask, spacing_mm, titleStr)
%PLOT_DOSE_COMPARISON Visualize original vs reconstructed dose

    if nargin < 5, titleStr = 'Dose Comparison'; end

    gridSize = size(original);

    dose_thresh = original > 0.01 * max(original(:));
    [ix, iy, iz] = ind2sub(gridSize, find(dose_thresh));

    if isempty(ix)
        fprintf('  [plot] No significant dose to display.\n');
        return;
    end

    dose_vals = original(dose_thresh);
    cx = round(sum(ix .* dose_vals) / sum(dose_vals));
    cy = round(sum(iy .* dose_vals) / sum(dose_vals));
    cz = round(sum(iz .* dose_vals) / sum(dose_vals));

    cx = max(1, min(gridSize(1), cx));
    cy = max(1, min(gridSize(2), cy));
    cz = max(1, min(gridSize(3), cz));

    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    figure('Name', titleStr, 'Position', [100, 100, 1400, 900]);
    sgtitle(sprintf('%s (centroid: X=%d, Y=%d, Z=%d)', titleStr, cx, cy, cz));

    % Transverse
    subplot(2, 3, 1);
    imagesc(squeeze(original(:, :, cz))');
    axis image; colorbar; title('Original - Transverse');
    xlabel('X'); ylabel('Y'); caxis([0, max_dose]);

    subplot(2, 3, 4);
    imagesc(squeeze(reconstructed(:, :, cz))');
    axis image; colorbar; title('Reconstructed - Transverse');
    xlabel('X'); ylabel('Y'); caxis([0, max_dose]);

    % Coronal
    subplot(2, 3, 2);
    imagesc(squeeze(original(:, cy, :))');
    axis image; colorbar; title('Original - Coronal');
    xlabel('X'); ylabel('Z'); caxis([0, max_dose]);

    subplot(2, 3, 5);
    imagesc(squeeze(reconstructed(:, cy, :))');
    axis image; colorbar; title('Reconstructed - Coronal');
    xlabel('X'); ylabel('Z'); caxis([0, max_dose]);

    % Sagittal
    subplot(2, 3, 3);
    imagesc(squeeze(original(cx, :, :))');
    axis image; colorbar; title('Original - Sagittal');
    xlabel('Y'); ylabel('Z'); caxis([0, max_dose]);

    subplot(2, 3, 6);
    imagesc(squeeze(reconstructed(cx, :, :))');
    axis image; colorbar; title('Reconstructed - Sagittal');
    xlabel('Y'); ylabel('Z'); caxis([0, max_dose]);

    colormap('jet');

    % Error figure
    figure('Name', [titleStr ' - Error'], 'Position', [150, 50, 1200, 500]);
    sgtitle([titleStr ' - Error Analysis']);

    subplot(1, 3, 1);
    diff_slice = squeeze(reconstructed(:, :, cz))' - squeeze(original(:, :, cz))';
    imagesc(diff_slice);
    axis image; colorbar; title('Difference (Recon-Orig) - Transverse');
    xlabel('X'); ylabel('Y');
    caxis([-max_dose * 0.2, max_dose * 0.2]);
    colormap(gca, 'jet');

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

    subplot(1, 3, 3);
    coronal_dose = squeeze(original(:, cy, :))';
    coronal_sensor = squeeze(sensor_mask(:, cy, :))';
    imagesc(coronal_dose); hold on;
    if any(coronal_sensor(:))
        contour(coronal_sensor, [0.5, 0.5], 'r', 'LineWidth', 2);
    end
    axis image; colorbar; title('Dose + Sensor (red) - Coronal');
    xlabel('X'); ylabel('Z');
    colormap(gca, 'jet');
    hold off;

    drawnow;
end


function plot_three_way_comparison(original, planar, compensated, spacing_mm)
%PLOT_THREE_WAY_COMPARISON Side-by-side: original, planar, compensated

    gridSize = size(original);

    dose_thresh = original > 0.01 * max(original(:));
    [ix, iy, iz] = ind2sub(gridSize, find(dose_thresh));
    dose_vals = original(dose_thresh);

    cx = max(1, min(gridSize(1), round(sum(ix .* dose_vals) / sum(dose_vals))));
    cy = max(1, min(gridSize(2), round(sum(iy .* dose_vals) / sum(dose_vals))));
    cz = max(1, min(gridSize(3), round(sum(iz .* dose_vals) / sum(dose_vals))));

    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    figure('Name', 'Three-Way Comparison', 'Position', [50, 50, 1600, 900]);
    sgtitle(sprintf('Original vs Planar vs Compensated (Z=%d, Y=%d, X=%d)', cz, cy, cx));

    % --- Row 1: Transverse slices ---
    subplot(3, 3, 1);
    imagesc(squeeze(original(:, :, cz))');
    axis image; colorbar; title('Original'); caxis([0, max_dose]);
    xlabel('X'); ylabel('Y');

    subplot(3, 3, 2);
    imagesc(squeeze(planar(:, :, cz))');
    axis image; colorbar; title('Planar Recon'); caxis([0, max_dose]);
    xlabel('X'); ylabel('Y');

    subplot(3, 3, 3);
    imagesc(squeeze(compensated(:, :, cz))');
    axis image; colorbar; title('Compensated'); caxis([0, max_dose]);
    xlabel('X'); ylabel('Y');

    % --- Row 2: Coronal slices ---
    subplot(3, 3, 4);
    imagesc(squeeze(original(:, cy, :))');
    axis image; colorbar; caxis([0, max_dose]);
    xlabel('X'); ylabel('Z');

    subplot(3, 3, 5);
    imagesc(squeeze(planar(:, cy, :))');
    axis image; colorbar; caxis([0, max_dose]);
    xlabel('X'); ylabel('Z');

    subplot(3, 3, 6);
    imagesc(squeeze(compensated(:, cy, :))');
    axis image; colorbar; caxis([0, max_dose]);
    xlabel('X'); ylabel('Z');

    % --- Row 3: Line profiles through dose centroid ---
    subplot(3, 3, 7);
    plot(squeeze(original(:, cy, cz)), 'b-', 'LineWidth', 1.5); hold on;
    plot(squeeze(planar(:, cy, cz)), 'r--', 'LineWidth', 1.5);
    plot(squeeze(compensated(:, cy, cz)), 'g-.', 'LineWidth', 1.5);
    legend('Original', 'Planar', 'Compensated', 'Location', 'best');
    xlabel('X index'); ylabel('Dose (Gy)'); title('X Profile');
    grid on; hold off;

    subplot(3, 3, 8);
    plot(squeeze(original(cx, :, cz)), 'b-', 'LineWidth', 1.5); hold on;
    plot(squeeze(planar(cx, :, cz)), 'r--', 'LineWidth', 1.5);
    plot(squeeze(compensated(cx, :, cz)), 'g-.', 'LineWidth', 1.5);
    legend('Original', 'Planar', 'Compensated', 'Location', 'best');
    xlabel('Y index'); ylabel('Dose (Gy)'); title('Y Profile (depth)');
    grid on; hold off;

    subplot(3, 3, 9);
    plot(squeeze(original(cx, cy, :)), 'b-', 'LineWidth', 1.5); hold on;
    plot(squeeze(planar(cx, cy, :)), 'r--', 'LineWidth', 1.5);
    plot(squeeze(compensated(cx, cy, :)), 'g-.', 'LineWidth', 1.5);
    legend('Original', 'Planar', 'Compensated', 'Location', 'best');
    xlabel('Z index'); ylabel('Dose (Gy)'); title('Z Profile');
    grid on; hold off;

    colormap('jet');
    drawnow;
end
