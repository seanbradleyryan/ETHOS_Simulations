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
%           .enable_spherical_correction - Boolean (default: true)
%           .regularization_lambda   - Wiener filter lambda (default: 0.01)
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
    enable_sph_corr    = safe_config(config, 'enable_spherical_correction', true);
    reg_lambda         = safe_config(config, 'regularization_lambda', 0.01);

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

    %% ============== SPHERICAL COMPENSATION CORRECTION ==============
    %  Limited-angle planar sensors introduce aperture artifacts.
    %  A Wiener-regularized correction is computed by comparing the planar
    %  reconstruction (p0_planar) with a reference full-angle reconstruction
    %  (p0_sphere) obtained from an enclosing 6-face sensor mask.
    %
    %  In Fourier space the compensation filter is:
    %     F(k) = S(k) * conj(P(k)) / ( |P(k)|^2 + lambda^2 * mean(|P|^2) )
    %  and the corrected image is:
    %     corrected = max( real( IFFT3( P(k) .* F(k) ) ), 0 )

    p0_planar = reconPressure;  % save planar result before overwriting

    if enable_sph_corr
        fprintf('        Running spherical compensation correction...\n');

        try
            sph_tic = tic;

            % --- 1. Build enclosing sensor mask (all 6 faces, 1 voxel inside PML) ---
            enclosing_sensor = build_enclosing_sensor(Nx, Ny, Nz, pml_size);
            numEncPts = sum(enclosing_sensor.mask(:));
            fprintf('          Enclosing sensor: %d points on 6 faces\n', numEncPts);

            % --- 2. Forward simulation with enclosing sensor ---
            fprintf('          Forward simulation (enclosing sensor)...\n');
            source_enc    = struct();
            source_enc.p0 = p0;   % same initial pressure as original forward

            enc_fwd_tic = tic;
            sensorData_enc = kspaceFirstOrder3D(kgrid, kmedium, source_enc, ...
                enclosing_sensor, inputArgs{:});
            enc_fwd_time = toc(enc_fwd_tic);
            fprintf('          Enclosing forward complete (%.1f s). Data: [%d x %d]\n', ...
                enc_fwd_time, size(sensorData_enc, 1), size(sensorData_enc, 2));

            % --- 3. Time-reversal from enclosing sensor -> p0_sphere ---
            fprintf('          Time reversal (enclosing sensor, 1 iter)...\n');
            source_enc_tr        = struct();
            source_enc_tr.p_mask = enclosing_sensor.mask;
            source_enc_tr.p      = fliplr(sensorData_enc);
            source_enc_tr.p_mode = 'dirichlet';

            sensor_enc_tr        = struct();
            sensor_enc_tr.mask   = ones(Nx, Ny, Nz);
            sensor_enc_tr.record = {'p_final'};

            enc_tr_tic = tic;
            p0_enc_raw = kspaceFirstOrder3D(kgrid, kmedium, source_enc_tr, ...
                sensor_enc_tr, inputArgs{:});
            enc_tr_time = toc(enc_tr_tic);

            if isstruct(p0_enc_raw) && isfield(p0_enc_raw, 'p_final')
                p0_sphere = reshape(p0_enc_raw.p_final, [Nx, Ny, Nz]);
            else
                p0_sphere = p0_enc_raw;
            end
            p0_sphere = max(p0_sphere, 0);

            fprintf('          Enclosing TR complete (%.1f s).\n', enc_tr_time);
            fprintf('          p0_sphere range: [%.2e, %.2e] Pa\n', ...
                min(p0_sphere(:)), max(p0_sphere(:)));

            % --- 4. Wiener-regularized compensation in Fourier space ---
            fprintf('          Computing Wiener compensation filter (lambda = %.4f)...\n', ...
                reg_lambda);

            P = fftn(p0_planar);
            S = fftn(p0_sphere);

            P_abs_sq  = abs(P).^2;
            noise_reg = reg_lambda^2 * mean(P_abs_sq(:));

            F = S .* conj(P) ./ (P_abs_sq + noise_reg);

            corrected = real(ifftn(P .* F));
            corrected = max(corrected, 0);  % positivity constraint

            % --- Report improvement ---
            fprintf('          Planar  pressure range: [%.2e, %.2e] Pa\n', ...
                min(p0_planar(:)), max(p0_planar(:)));
            fprintf('          Corrected pressure range: [%.2e, %.2e] Pa\n', ...
                min(corrected(:)), max(corrected(:)));

            reconPressure = corrected;

            sph_time = toc(sph_tic);
            fprintf('        Spherical compensation complete (%.1f s total).\n', sph_time);

        catch ME
            warning('run_single_field_simulation:SphCorrFail', ...
                'Spherical compensation failed: %s. Using planar result.', ME.message);
            reconPressure = p0_planar;
        end
    else
        fprintf('        Spherical compensation: disabled.\n');
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


function sensor = build_enclosing_sensor(Nx, Ny, Nz, pml_size)
%BUILD_ENCLOSING_SENSOR Create sensor mask on all 6 faces of the 3D grid
%
%   Builds a binary mask with active sensor points on each face of the
%   computational grid, inset from the PML boundary by 1 voxel. This
%   provides near-complete angular coverage for reference reconstruction.
%
%   INPUTS:
%       Nx, Ny, Nz - Grid dimensions
%       pml_size   - PML layer thickness in voxels
%
%   OUTPUTS:
%       sensor - Struct with .mask field (Nx x Ny x Nz binary array)

    sensor = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    % Face positions: 1 voxel inside the PML boundary
    face_lo = pml_size + 1;
    face_hi_x = Nx - pml_size;
    face_hi_y = Ny - pml_size;
    face_hi_z = Nz - pml_size;

    % Interior range (excluding edge overlaps for clean face assignment)
    ix = face_lo:face_hi_x;
    iy = face_lo:face_hi_y;
    iz = face_lo:face_hi_z;

    % -X face (x = face_lo) and +X face (x = face_hi_x)
    sensor.mask(face_lo, iy, iz) = 1;
    sensor.mask(face_hi_x, iy, iz) = 1;

    % -Y face (y = face_lo) and +Y face (y = face_hi_y)
    sensor.mask(ix, face_lo, iz) = 1;
    sensor.mask(ix, face_hi_y, iz) = 1;

    % -Z face (z = face_lo) and +Z face (z = face_hi_z)
    sensor.mask(ix, iy, face_lo) = 1;
    sensor.mask(ix, iy, face_hi_z) = 1;
end


function val = safe_config(config, field_name, default_val)
%SAFE_CONFIG Retrieve config field with fallback to default
    if isfield(config, field_name) && ~isempty(config.(field_name))
        val = config.(field_name);
    else
        val = default_val;
    end
end
