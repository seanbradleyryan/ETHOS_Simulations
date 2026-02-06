function recon_dose = run_single_field_simulation(field_dose, sct_resampled, medium, config)
%RUN_SINGLE_FIELD_SIMULATION k-Wave forward + time-reversal for one field
%
%   recon_dose = run_single_field_simulation(field_dose, sct_resampled, medium, config)
%
%   Converts a single radiation field dose to initial acoustic pressure
%   (p0 = D * Gamma * rho), runs the k-Wave forward simulation to generate
%   synthetic sensor data, then applies time-reversal reconstruction to
%   recover the initial pressure distribution and converts back to dose.
%
%   INPUTS:
%       field_dose - Struct from step15_process_doses:
%           .dose_Gy       - 3D dose array (Gy)
%           .gantry_angle  - Gantry angle (degrees)
%           .meterset      - Monitor units (MU)
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%       sct_resampled - Struct with:
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%       medium - Struct from create_acoustic_medium():
%           .density       - 3D array (kg/m^3)
%           .sound_speed   - 3D array (m/s)
%           .alpha_coeff   - 3D array (dB/MHz^y/cm)
%           .alpha_power   - Scalar exponent
%           .gruneisen     - 3D array (dimensionless)
%           .grid_size     - [nx, ny, nz]
%       config - Configuration struct:
%           .dose_per_pulse_cGy      - Dose per LINAC pulse (default: 0.16)
%           .pml_size                - PML thickness, voxels (default: 10)
%           .cfl_number              - CFL stability number (default: 0.3)
%           .use_gpu                 - Boolean (default: true)
%           .num_time_reversal_iter  - TR iterations (default: 1)
%
%   OUTPUTS:
%       recon_dose - 3D reconstructed dose array (Gy), same size as input
%
%   NOTES:
%       - Stateless: safe for parfor execution
%       - PlotSim disabled for batch mode
%       - Returns zeros if no significant dose or on simulation failure
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: create_acoustic_medium, kspaceFirstOrder3D, kWaveGrid

    %% ======================== CONFIG DEFAULTS ========================

    dose_per_pulse_cGy = safe_config(config, 'dose_per_pulse_cGy', 0.16);
    pml_size           = safe_config(config, 'pml_size', 10);
    cfl                = safe_config(config, 'cfl_number', 0.3);
    use_gpu            = safe_config(config, 'use_gpu', true);
    num_tr_iter        = safe_config(config, 'num_time_reversal_iter', 1);

    %% ======================== EXTRACT DATA ========================

    doseGrid = field_dose.dose_Gy;
    gridSize = size(doseGrid);
    Nx = gridSize(1);
    Ny = gridSize(2);
    Nz = gridSize(3);

    % Grid spacing: mm -> m
    spacing_mm = sct_resampled.spacing(:)';
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;

    % Local copies of acoustic property arrays (parfor safety)
    density    = medium.density;
    soundSpeed = medium.sound_speed;
    gruneisen  = medium.gruneisen;

    %% ======================== VALIDATE DIMENSIONS ========================

    if ~isequal(gridSize, medium.grid_size)
        error('run_single_field_simulation:SizeMismatch', ...
            'Field dose [%d %d %d] does not match medium [%d %d %d].', ...
            Nx, Ny, Nz, ...
            medium.grid_size(1), medium.grid_size(2), medium.grid_size(3));
    end

    %% ======================== PULSE CALCULATION ========================

    meterset = field_dose.meterset;
    if isempty(meterset) || meterset <= 0
        warning('run_single_field_simulation:NoMeterset', ...
            'Invalid meterset (%.2f). Falling back to 100 MU.', meterset);
        meterset = 100;
    end
    num_pulses = ceil(meterset / dose_per_pulse_cGy);

    fprintf('        Meterset: %.2f MU, Pulses: %d\n', meterset, num_pulses);

    %% ======================== INITIAL PRESSURE p0 ========================

    % p0(r) = D(r)/N_pulses * Gamma(r) * rho(r)
    dose_per_pulse = doseGrid / num_pulses;
    p0 = dose_per_pulse .* gruneisen .* density;

    fprintf('        Max dose: %.4f Gy, Per-pulse max: %.6f Gy\n', ...
        max(doseGrid(:)), max(dose_per_pulse(:)));
    fprintf('        Initial pressure range: [%.2e, %.2e] Pa\n', ...
        min(p0(:)), max(p0(:)));

    %% ======================== CHECK FOR SIGNIFICANT DOSE ========================

    doseThreshold = 0.01 * max(doseGrid(:));  % 1% of max
    doseMask = doseGrid > doseThreshold;

    if ~any(doseMask(:)) || max(p0(:)) == 0
        warning('run_single_field_simulation:NoDose', ...
            'No significant dose or zero initial pressure. Returning zeros.');
        recon_dose = zeros(gridSize);
        return;
    end

    %% ======================== k-WAVE GRID ========================

    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    % CFL-stable time step
    maxC = max(soundSpeed(:));
    minC = min(soundSpeed(soundSpeed > 0));
    dt   = cfl * min([dx, dy, dz]) / maxC;

    % Simulation time: 2.5x grid diagonal traversal at minimum speed
    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);

    kgrid.dt = dt;
    kgrid.Nt = Nt;

    fprintf('        Grid: [%d x %d x %d], spacing: [%.3f, %.3f, %.3f] mm\n', ...
        Nx, Ny, Nz, dx*1000, dy*1000, dz*1000);
    fprintf('        dt = %.2e s, Nt = %d, T_sim = %.2e s\n', dt, Nt, simTime);

    %% ======================== k-WAVE MEDIUM ========================

    kmedium = struct();
    kmedium.density     = density;
    kmedium.sound_speed = soundSpeed;
    kmedium.alpha_coeff = medium.alpha_coeff;
    kmedium.alpha_power = medium.alpha_power;  % scalar

    %% ======================== SENSOR PLACEMENT ========================

    gantry_angle = field_dose.gantry_angle;
    sensor = place_sensor_for_field(doseMask, Nx, Ny, Nz, gantry_angle);

    numSensorPts = sum(sensor.mask(:));
    fprintf('        Sensor: %d points (gantry = %.1f deg)\n', numSensorPts, gantry_angle);

    if numSensorPts == 0
        warning('run_single_field_simulation:EmptySensor', ...
            'Sensor mask is empty. Returning zeros.');
        recon_dose = zeros(gridSize);
        return;
    end

    %% ======================== DATA CAST (GPU/CPU) ========================

    if use_gpu
        try
            gpuDevice;
            dataCast = 'gpuArray-single';
            fprintf('        Compute: GPU\n');
        catch
            dataCast = 'single';
            fprintf('        Compute: CPU (GPU unavailable)\n');
        end
    else
        dataCast = 'single';
        fprintf('        Compute: CPU\n');
    end

    inputArgs = {'Smooth', false, ...
                 'PMLInside', false, ...
                 'PMLSize', pml_size, ...
                 'DataCast', dataCast, ...
                 'PlotSim', false};

    %% ======================== FORWARD SIMULATION ========================

    fprintf('        Running forward simulation...\n');

    source_fwd    = struct();
    source_fwd.p0 = p0;

    try
        fwd_tic = tic;
        sensorData = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
        fwd_time = toc(fwd_tic);
        fprintf('        Forward complete (%.1f s). Sensor data: [%d x %d]\n', ...
            fwd_time, size(sensorData, 1), size(sensorData, 2));
    catch ME
        warning('run_single_field_simulation:ForwardFail', ...
            'Forward simulation failed: %s', ME.message);
        recon_dose = zeros(gridSize);
        return;
    end

    %% ======================== TIME REVERSAL RECONSTRUCTION ========================

    fprintf('        Running time reversal (%d iteration(s))...\n', num_tr_iter);

    reconPressure = zeros(gridSize);

    try
        tr_tic = tic;

        for tr_iter = 1:num_tr_iter
            if num_tr_iter > 1
                fprintf('          TR iteration %d/%d...\n', tr_iter, num_tr_iter);
            end

            % Time-reversed source on sensor locations
            source_tr        = struct();
            source_tr.p_mask = sensor.mask;
            source_tr.p      = fliplr(sensorData);  % time-reversed
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

            % Positivity constraint (dose and pressure are non-negative)
            reconPressure = max(reconPressure, 0);

            % Iterative TR: compute residual for next iteration
            if tr_iter < num_tr_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, ...
                    source_resid, sensor, inputArgs{:});

                % Residual correction
                sensorData = sensorData + (sensorData - sensorDataRecon);
            end
        end

        tr_time = toc(tr_tic);
        fprintf('        Time reversal complete (%.1f s).\n', tr_time);
        fprintf('        Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
            min(reconPressure(:)), max(reconPressure(:)));

    catch ME
        warning('run_single_field_simulation:TRFail', ...
            'Time reversal failed: %s', ME.message);
        recon_dose = zeros(gridSize);
        return;
    end

    %% ======================== PRESSURE -> DOSE CONVERSION ========================

    % D_recon = p0_recon / (Gamma * rho) * num_pulses
    conversionFactor = gruneisen .* density;
    conversionFactor(conversionFactor == 0) = 1;  % prevent div-by-zero

    reconDosePerPulse = reconPressure ./ conversionFactor;
    recon_dose = reconDosePerPulse * num_pulses;

    fprintf('        Reconstructed dose: [%.4f, %.4f] Gy\n', ...
        min(recon_dose(:)), max(recon_dose(:)));
    fprintf('        Field simulation complete.\n');
end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function sensor = place_sensor_for_field(doseMask, Nx, Ny, Nz, gantry_angle)
%PLACE_SENSOR_FOR_FIELD Create planar sensor based on dose extent and gantry
%
%   Sensor is placed on the beam exit side, determined by gantry angle:
%     IEC gantry angle convention (looking from isocenter):
%       0   deg -> beam from anterior (+Y),  sensor on posterior (-Y side)
%       90  deg -> beam from right   (-X),   sensor on left      (+X side)
%       180 deg -> beam from posterior(-Y),  sensor on anterior  (+Y side)
%       270 deg -> beam from left    (+X),   sensor on right     (-X side)
%
%   The sensor plane covers the X-Z (or Y-Z) extent of the dose with a
%   5-voxel pad, offset 10 voxels beyond the dose boundary on the exit side.

    MARGIN = 10;   % voxels beyond dose extent
    PAD_XZ = 5;    % padding around dose extent in plane

    % Find dose voxel indices
    [xIdx, yIdx, zIdx] = ind2sub([Nx, Ny, Nz], find(doseMask));

    sensor = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    if isempty(xIdx)
        return;
    end

    % Dose extent with padding
    xMin = max(1,  min(xIdx) - PAD_XZ);
    xMax = min(Nx, max(xIdx) + PAD_XZ);
    yMin = max(1,  min(yIdx) - PAD_XZ);
    yMax = min(Ny, max(yIdx) + PAD_XZ);
    zMin = max(1,  min(zIdx) - PAD_XZ);
    zMax = min(Nz, max(zIdx) + PAD_XZ);

    % Determine exit side from gantry angle (quantized to quadrants)
    ga = mod(gantry_angle, 360);

    if (ga >= 315 || ga < 45)
        % Beam from anterior -> sensor on posterior (low Y)
        sY = max(1, min(yIdx) - MARGIN);
        sensor.mask(xMin:xMax, sY, zMin:zMax) = 1;

    elseif (ga >= 45 && ga < 135)
        % Beam from right -> sensor on left (high X)
        sX = min(Nx, max(xIdx) + MARGIN);
        sensor.mask(sX, yMin:yMax, zMin:zMax) = 1;

    elseif (ga >= 135 && ga < 225)
        % Beam from posterior -> sensor on anterior (high Y)
        sY = min(Ny, max(yIdx) + MARGIN);
        sensor.mask(xMin:xMax, sY, zMin:zMax) = 1;

    else  % 225 <= ga < 315
        % Beam from left -> sensor on right (low X)
        sX = max(1, min(xIdx) - MARGIN);
        sensor.mask(sX, yMin:yMax, zMin:zMax) = 1;
    end
end


function val = safe_config(config, field_name, default_val)
%SAFE_CONFIG Retrieve config field with fallback to default
    if isfield(config, field_name) && ~isempty(config.(field_name))
        val = config.(field_name);
    else
        val = default_val;
    end
end
