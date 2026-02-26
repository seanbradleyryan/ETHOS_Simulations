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
%  Structural pattern mirrors run_single_field_simulation.m:
%    forward simulation -> iterative TR -> PSF correction -> dose recovery
%
%  CONFIGURABLE OPTIONS (see CONFIGURATION section):
%    - Patient/session selection (constructs path to processed/ dir)
%    - Dose file selection (default: total_rs_dose.mat in processed/)
%    - Sensor x-index for planar sensor placement
%    - Tissue heterogeneity: 'uniform', 'threshold_1', 'threshold_2'
%    - Per-property heterogeneity overrides
%    - PSF correction via get_psf (Wiener-regularised frequency-domain filter)
%    - GPU/CPU selection, PML size, CFL number, TR iterations
%    - Movie recording for first forward simulation and first TR iteration
%
%  PSF CORRECTION (CONFIG.use_psf_correction = true):
%    Calls get_psf() once to compute a Wiener-regularised frequency-domain
%    filter that compensates for limited-angle planar-sensor artifacts:
%      F(k) = FFT3(p0_sphere) * conj(P(k)) / (|P(k)|^2 + lambda^2 * mean(|P|^2))
%    Applied to the planar reconstruction:
%      p0_corrected = IFFT3( FFT3(p0_planar) * F(k) )
%
%  REQUIRED FILES (in processed/ directory from step15_process_doses):
%    - total_rs_dose.mat  (contains total_rs_dose: 3D Gy array)
%    - sct_resampled.mat  (contains sct_resampled struct with cubeHU,
%                          cubeDensity, bodyMask, couchMask, spacing, etc.)
%
%  PREREQUISITES:
%    - MATLAB R2022a or later
%    - k-Wave Toolbox (http://www.k-wave.org)
%    - Image Processing Toolbox (optional)
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: February 2026
%  VERSION: 3.0 (PSF correction via get_psf, conforms to run_single_field_simulation)
%
%  See also: run_single_field_simulation, get_psf, determine_sensor_mask,
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

% --- Sensor Placement ---
% YZ plane sensor at a fixed X index, matching run_single_field_simulation.
% X = 1 corresponds to the face at x-min (lateral edge).
CONFIG.sensor_x_index = 1;

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
CONFIG.uniform_alpha_coeff  = 0;       % dB/MHz^y/cm
CONFIG.uniform_alpha_power  = 1.1;     % exponent
CONFIG.uniform_gruneisen    = 1.0;     % dimensionless

% --- Simulation Parameters ---
CONFIG.dose_per_pulse_cGy     = 0.16;   % cGy per LINAC pulse
CONFIG.meterset               = 140;    % MU (monitor units) for total dose
CONFIG.pml_size               = 10;     % PML thickness (voxels)
CONFIG.cfl_number             = 0.3;    % CFL stability number
CONFIG.use_gpu                = true;   % Use GPU acceleration

% --- Iterative Time-Reversal Reconstruction ---
%   Uses Dirichlet-BC time-reversal with iterative residual correction:
%     Iter 1: TR(measured_data) -> p0_est, positivity constraint
%     Iter n: forward(p0_est) -> simulated_data
%             residual = measured_data - simulated_data
%             measured_data += residual   (corrected data)
%             TR(corrected_data) -> p0_est (updated), positivity constraint
%             check convergence
CONFIG.num_time_reversal_iter = 5;       % Maximum TR iterations
CONFIG.convergence_tol        = 1e-4;   % Early stop if relative change < tol

% --- PSF Correction ---
%   Calls get_psf() once (using the total dose as calibration source) to
%   compute a Wiener-regularised frequency-domain filter that compensates
%   for limited-angle artifacts in the planar reconstruction.
%   The filter is then applied after the iterative TR loop.
CONFIG.use_psf_correction      = true;  % Master toggle
CONFIG.regularization_lambda   = 0.01;  % Wiener regularisation (get_psf)

% --- Movie Recording ---
%   k-Wave can record the visualised simulation as a movie file.
%   Enabled for the FIRST forward simulation and the FIRST TR iteration.
%   Subsequent calls (residual passes, iterations 2+) use no recording.
CONFIG.record_movie = true;   % Record movies for first fwd + first TR

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
fprintf('  Standalone k-Wave Photoacoustic Simulation  (v3.0)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', dose_filepath);
fprintf('  SCT file:        %s\n', sct_filepath);
fprintf('  Sensor:          YZ plane at x = %d\n', CONFIG.sensor_x_index);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  TR iterations:   %d (tol: %.1e)\n', CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);
fprintf('  GPU:             %s\n', mat2str(CONFIG.use_gpu));
if CONFIG.use_psf_correction
    fprintf('  PSF correction:  ENABLED (lambda=%.1e)\n', CONFIG.regularization_lambda);
else
    fprintf('  PSF correction:  off\n');
end
fprintf('  Movie recording: %s\n', mat2str(CONFIG.record_movie));
fprintf('=========================================================\n\n');

%% ========================= LOAD DATA ====================================

fprintf('[1/7] Loading dose data...\n');
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
fprintf('[2/7] Loading SCT data...\n');
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

fprintf('[3/7] Creating acoustic medium (method: %s)...\n', CONFIG.gruneisen_method);

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

fprintf('[4/7] Computing initial pressure...\n');

% Pulse calculation
meterset   = CONFIG.meterset;
num_pulses = ceil(meterset / CONFIG.dose_per_pulse_cGy);
dose_per_pulse = doseGrid / num_pulses;

% p0(r) = D_per_pulse(r) * Gamma(r) * rho(r)
p0 = dose_per_pulse .* medium.gruneisen .* medium.density;

fprintf('       Meterset: %.2f MU -> %d pulses\n', meterset, num_pulses);
fprintf('       Max dose per pulse: %.6f Gy\n', max(dose_per_pulse(:)));
fprintf('       p0 range: [%.2e, %.2e] Pa\n', min(p0(:)), max(p0(:)));

% Check for significant dose
doseThreshold = 0.01 * max(doseGrid(:));
doseMask      = doseGrid > doseThreshold;

if ~any(doseMask(:)) || max(p0(:)) == 0
    warning('No significant dose or zero initial pressure. Aborting.');
    return;
end

%% ========================= SENSOR PLACEMENT ==============================

fprintf('[5/7] Placing sensor (YZ plane at x = %d)...\n', CONFIG.sensor_x_index);

sensor        = struct();
sensor.mask   = zeros(Nx, Ny, Nz);
sensor.mask(CONFIG.sensor_x_index, :, :) = 1;

numSensorPts = sum(sensor.mask(:));
fprintf('       Sensor: %d active points (%d x %d YZ plane)\n', ...
    numSensorPts, Ny, Nz);

if numSensorPts == 0
    warning('Sensor mask is empty. Aborting.');
    return;
end

% --- COMMENTED OUT: multi-mode sensor placement (for reference) ----------
% CONFIG.sensor_mode  — 'full_anterior_plane' | 'dose_based' | 'spherical'
% CONFIG.gantry_angle — degrees (only for dose_based)
%
% switch lower(CONFIG.sensor_mode)
%     case 'full_anterior_plane'
%         sensor.mask = zeros(Nx, Ny, Nz);
%         if isfield(sct, 'bodyMask')
%             body = sct.bodyMask;
%             if isfield(sct, 'couchMask')
%                 body = body & ~sct.couchMask;
%             end
%             y_has_body = squeeze(any(any(body, 1), 3));
%             anterior_y = find(y_has_body, 1, 'first');
%             if ~isempty(anterior_y)
%                 sensor_y = max(1, anterior_y - 3);
%             else
%                 sensor_y = 1;
%             end
%         else
%             sensor_y = 1;
%         end
%         sensor.mask(:, sensor_y, :) = 1;
%     case 'dose_based'
%         sensor = place_sensor_for_field(doseMask, Nx, Ny, Nz, CONFIG.gantry_angle);
%     case 'spherical'
%         radius = floor(min([Nx, Ny, Nz]) / 2) - 1;
%         sensor.mask = makeSphere(Nx, Ny, Nz, radius);
% end
% -------------------------------------------------------------------------

%% ========================= OPTIMAL GRID PADDING ==========================
%  Pad grid to FFT-friendly dimensions for k-Wave performance.
%  Original data sits at indices 1:N_orig; padding at N_orig+1:N_pad.
%  Padding region filled with water medium properties.

fprintf('[PAD] Computing FFT-optimal padded dimensions...\n');

Nx_orig = Nx;  Ny_orig = Ny;  Nz_orig = Nz;
gridSize_orig    = gridSize;
medium_orig      = medium;
sensor_mask_orig = sensor.mask;

Nx_pad = find_optimal_kwave_size(Nx, CONFIG.pml_size);
Ny_pad = find_optimal_kwave_size(Ny, CONFIG.pml_size);
Nz_pad = find_optimal_kwave_size(Nz, CONFIG.pml_size);

did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
if did_pad
    fprintf('[PAD] Padding grid: [%d %d %d] -> [%d %d %d] (FFT-optimal)\n', ...
        Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);

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
    gruneisen_pad(1:Nx, 1:Ny, 1:Nz)  = medium.gruneisen;

    medium.density     = density_pad;
    medium.sound_speed = soundSpeed_pad;
    medium.alpha_coeff = alphaCoeff_pad;
    medium.gruneisen   = gruneisen_pad;

    p0_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
    p0_pad(1:Nx, 1:Ny, 1:Nz) = p0;
    p0 = p0_pad;

    sensor_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
    sensor_pad(1:Nx, 1:Ny, 1:Nz) = sensor.mask;
    sensor.mask = sensor_pad;

    Nx = Nx_pad;  Ny = Ny_pad;  Nz = Nz_pad;
    gridSize = [Nx, Ny, Nz];
else
    fprintf('[PAD] Grid [%d %d %d] already FFT-optimal, no padding needed.\n', Nx, Ny, Nz);
end

%% ========================= k-WAVE GRID & MEDIUM SETUP ===================

fprintf('[6/7] Setting up k-Wave grid...\n');

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
kmedium             = struct();
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

% ---- inputArgs ----
% Standard: no movie, PlotSim off (used for residual passes and iterations 2+)
inputArgs = {'Smooth', false, ...
             'PMLInside', false, ...
             'PMLSize', CONFIG.pml_size, ...
             'DataCast', dataCast, ...
             'PlotSim', false};

% Movie-enabled: PlotSim on, movie recorded.
% Used for FIRST forward simulation and FIRST TR iteration only.
if CONFIG.record_movie
    inputArgs_fwd_movie = {'Smooth', false, ...
                           'PMLInside', false, ...
                           'PMLSize', CONFIG.pml_size, ...
                           'DataCast', dataCast, ...
                           'PlotSim', true, ...
                           'RecordMovie', true, ...
                           'MovieName', 'forward_sim'};
    inputArgs_tr_movie  = {'Smooth', false, ...
                           'PMLInside', false, ...
                           'PMLSize', CONFIG.pml_size, ...
                           'DataCast', dataCast, ...
                           'PlotSim', true, ...
                           'RecordMovie', true, ...
                           'MovieName', 'time_reversal'};
else
    inputArgs_fwd_movie = inputArgs;
    inputArgs_tr_movie  = inputArgs;
end

%% ========================= PSF CORRECTION FILTER =========================

psf_filter = [];
if CONFIG.use_psf_correction
    fprintf('\n[PSF] Computing PSF correction filter (get_psf)...\n');
    try
        psf_filter = get_psf(doseGrid, sct, medium_orig, CONFIG);
        fprintf('[PSF] Filter ready. Grid: [%d x %d x %d]\n', psf_filter.grid_size);
    catch ME
        warning('run_standalone_simulation:PSFFail', ...
            'get_psf failed: %s. Proceeding without PSF correction.', ME.message);
        psf_filter = [];
    end
end

%% ========================= FORWARD SIMULATION ============================

fprintf('\n[7/7] Running k-Wave forward simulation...\n');

source_fwd    = struct();
source_fwd.p0 = p0;

try
    fwd_tic = tic;
    sensorData = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, ...
        inputArgs_fwd_movie{:});
    fwd_time = toc(fwd_tic);
    fprintf('       Forward complete (%.1f s). Sensor data: [%d x %d]\n', ...
        fwd_time, size(sensorData, 1), size(sensorData, 2));
catch ME
    fprintf('[ERROR] Forward simulation failed: %s\n', ME.message);
    return;
end

% Store the original measured sensor data (used in residual correction)
sensorData_measured = sensorData;

%% ========================= ITERATIVE TIME-REVERSAL RECONSTRUCTION ========
%
%  Algorithm (Dirichlet-BC iterative TR):
%    Iteration 1:
%      1. Time-reverse measured sensor data (fliplr) as Dirichlet source
%      2. Record p_final over entire grid -> p0_est
%      3. Apply positivity constraint
%
%    Iterations 2..N:
%      4. Forward-simulate current estimate -> simulated sensor data
%      5. Compute residual: measured - simulated
%      6. Correct: working_data = measured + residual
%      7. Time-reverse corrected data -> p0_est (updated)
%      8. Apply positivity constraint
%      9. Check convergence
%  =========================================================================

fprintf('       Running iterative time reversal (%d iterations, tol=%.1e)...\n', ...
    CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);

reconPressure      = zeros(gridSize);
reconPressure_prev = zeros(gridSize);

try
    tr_total_tic = tic;

    for tr_iter = 1:CONFIG.num_time_reversal_iter

        fprintf('       --- TR Iteration %d/%d ---\n', tr_iter, CONFIG.num_time_reversal_iter);

        % Use movie inputArgs for first iteration only
        if tr_iter == 1
            iter_inputArgs = inputArgs_tr_movie;
        else
            iter_inputArgs = inputArgs;
        end

        % Step A: Time-reverse current sensor data
        source_tr        = struct();
        source_tr.p_mask = sensor.mask;
        source_tr.p      = fliplr(sensorData);  % time-reversed
        source_tr.p_mode = 'dirichlet';

        sensor_tr        = struct();
        sensor_tr.mask   = ones(Nx, Ny, Nz);
        sensor_tr.record = {'p_final'};

        p0_recon = kspaceFirstOrder3D(kgrid, kmedium, source_tr, sensor_tr, iter_inputArgs{:});

        % Extract 3D pressure field
        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
        else
            reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
        end

        % Positivity constraint (dose and pressure are non-negative)
        reconPressure = max(reconPressure, 0);

        fprintf('       Max pressure: %.4e Pa\n', max(reconPressure(:)));

        % Convergence check (from iteration 2 onward)
        if tr_iter > 1
            norm_prev = norm(reconPressure_prev(:));
            if norm_prev > 0
                rel_change = norm(reconPressure(:) - reconPressure_prev(:)) / norm_prev;
            else
                rel_change = Inf;
            end
            fprintf('       Rel change: %.4e\n', rel_change);
            if rel_change < CONFIG.convergence_tol
                fprintf('       *** Converged at iteration %d ***\n', tr_iter);
                break;
            end
        end

        reconPressure_prev = reconPressure;

        % Step B: Residual correction for next iteration
        if tr_iter < CONFIG.num_time_reversal_iter
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, ...
                source_resid, sensor, inputArgs{:});

            % Update working sensor data with residual
            sensorData = sensorData_measured + (sensorData_measured - sensorDataRecon);
        end
    end

    tr_time = toc(tr_total_tic);
    fprintf('       Time reversal complete (%.1f s).\n', tr_time);
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] Time reversal failed: %s\n', ME.message);
    return;
end

%% ========================= CROP TO ORIGINAL SIZE =========================
%  Remove padding before PSF correction and dose conversion.
%  psf_filter.F is computed at the original grid size; reconPressure must
%  match.  medium.gruneisen/density used for dose conversion must also be
%  the original tissue-property arrays, not the water-padded versions.

if did_pad
    fprintf('\n[CROP] Restoring original dimensions: [%d %d %d] -> [%d %d %d]\n', ...
        Nx, Ny, Nz, Nx_orig, Ny_orig, Nz_orig);
    reconPressure = reconPressure(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    Nx = Nx_orig;  Ny = Ny_orig;  Nz = Nz_orig;
    gridSize    = gridSize_orig;
    medium      = medium_orig;
    sensor.mask = sensor_mask_orig;
end

%% ========================= PSF CORRECTION ================================
%
%  If a pre-computed PSF filter (from get_psf) is available, apply it to
%  the planar reconstruction to compensate for limited-view artifacts.
%  Matches the PSF application block in run_single_field_simulation.m.
%  =========================================================================

psf_applied = false;
if ~isempty(psf_filter) && isstruct(psf_filter) && isfield(psf_filter, 'F') ...
        && ~isempty(psf_filter.F)
    fprintf('       Applying pre-computed PSF correction...\n');
    P_field   = fftn(reconPressure);
    corrected = real(ifftn(P_field .* psf_filter.F));
    reconPressure = max(corrected, 0);
    psf_applied = true;
    fprintf('       Corrected pressure range: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));
end

%% ========================= PRESSURE -> DOSE ==============================

fprintf('\n[Post] Converting pressure to dose...\n');

conversionFactor = medium.gruneisen .* medium.density;
conversionFactor(conversionFactor == 0) = 1;  % prevent div-by-zero

reconDosePerPulse = reconPressure ./ conversionFactor;
recon_dose        = reconDosePerPulse * num_pulses;

fprintf('       Reconstructed dose: [%.4f, %.4f] Gy\n', ...
    min(recon_dose(:)), max(recon_dose(:)));

%% ========================= CROP p0 TO ORIGINAL SIZE ======================
%  reconPressure and recon_dose are already original size (cropped above).
%  Crop p0 here so the saved results struct is self-consistent.

if did_pad
    p0 = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
end

%% ========================= RESULTS SUMMARY ===============================

fprintf('\n========= RESULTS =========\n');
fprintf('  Original dose:      [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));
fprintf('  Reconstructed dose: [%.6f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));
fprintf('  PSF correction:     %s\n', mat2str(psf_applied));

% Error metrics within dose region
dose_region = doseGrid > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
    rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
    fprintf('  Mean abs err: %.6f Gy\n', mean(abs_error));
    fprintf('  Mean rel err: %.2f%%\n',  mean(rel_error));
    fprintf('  Max  rel err: %.2f%%\n',  max(rel_error));
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
    results.tr_time_sec    = tr_time;
    results.psf_applied    = psf_applied;
    if psf_applied
        results.psf_filter = psf_filter;
    end

    save(CONFIG.output_file, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', CONFIG.output_file);
end

%% ========================= VISUALIZATION =================================

if CONFIG.plot_results
    plot_dose_comparison(doseGrid, recon_dose, sensor.mask, spacing_mm, ...
        'Standalone Reconstruction');
end

fprintf('\nStandalone simulation complete.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

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
    orig_slice  = squeeze(original(:, :, cz))';
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
    coronal_dose   = squeeze(original(:, cy, :))';
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
