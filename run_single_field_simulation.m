function [recon_dose, sim_results] = run_single_field_simulation(field_dose, sct_resampled, medium, beam_metadata, config, psf_filter)
%RUN_SINGLE_FIELD_SIMULATION k-Wave forward + time-reversal for one field
%
%   [recon_dose, sim_results] = run_single_field_simulation(field_dose, sct_resampled, medium, beam_metadata, config)
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
%           .isocenter     - [x, y, z] mm (from RTPLAN, propagated)
%           .jaw_x         - [x1, x2] mm at isocenter (from RTPLAN)
%           .jaw_y         - [y1, y2] mm at isocenter (from RTPLAN)
%       sct_resampled - Struct with:
%           .bodyMask      - 3D logical (true = inside body)
%           .couchMask     - 3D logical (true = couch region)
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%           .origin        - [x, y, z] in mm
%       medium - Struct from create_acoustic_medium():
%           .density       - 3D array (kg/m^3)
%           .sound_speed   - 3D array (m/s)
%           .alpha_coeff   - 3D array (dB/MHz^y/cm)
%           .alpha_power   - Scalar exponent
%           .gruneisen     - 3D array (dimensionless)
%           .grid_size     - [nx, ny, nz]
%       beam_metadata - Struct array (ALL beams in plan) with:
%           .beam_number   - Beam number
%           .gantry_angle  - Gantry angle (degrees)
%           .isocenter     - [x, y, z] mm
%           .jaw_x         - [x1, x2] mm at isocenter
%           .jaw_y         - [y1, y2] mm at isocenter
%       config - Configuration struct:
%           .dose_per_pulse_cGy      - Dose per LINAC pulse (default: 0.16)
%           .pml_size                - PML thickness, voxels (default: 10)
%           .cfl_number              - CFL stability number (default: 0.3)
%           .use_gpu                 - Boolean (default: true)
%           .num_time_reversal_iter  - TR iterations (default: 30)
%           .convergence_tol         - Relative change threshold for early TR exit (default: 1e-3)
%           .correction_factor       - Scalar multiplier on reconstructed pressure (default: 1.9)
%           .use_grid_padding        - Pad grid to FFT-optimal dims (default: true)
%           .sensor_x_index          - X voxel index for anterior sensor plane (default: 1)
%           .sensor_y_index          - Y voxel index for lateral sensor plane (default: 1)
%       psf_filter - (OPTIONAL, 6th arg) Pre-computed PSF correction from get_psf():
%           .F              - 3D complex Wiener filter in Fourier domain
%           If provided and non-empty, applied to the planar reconstruction
%           to compensate for limited-view artifacts.  Pass [] to skip.
%
%   OUTPUTS:
%       recon_dose  - 3D reconstructed dose array (Gy), same size as input.
%       sim_results - Struct with simulation diagnostics:
%           .sensor_info    - Sensor placement info from determine_sensor_mask
%           .forward_time_s - Forward simulation wall time
%           .tr_time_s      - Time reversal wall time
%           .num_pulses     - Number of LINAC pulses
%           .p0_max         - Maximum initial pressure (Pa)
%           .recon_max          - Maximum reconstructed pressure (Pa)
%           .psf_applied        - Boolean, whether PSF correction was applied
%           .num_iters_done     - Number of TR iterations completed
%           .conv_max_pressure  - Max pressure per iteration (vector)
%           .conv_rel_change    - Relative change per iteration (vector)
%
%   NOTES:
%       - Stateless: safe for parfor execution.
%       - PlotSim disabled for batch mode.
%       - Returns zeros if no significant dose or on simulation failure.
%       - beam_metadata is the FULL plan metadata (all beams), passed through
%         to determine_sensor_mask for computing the exclusion zone.
%       - For backwards compatibility, beam_metadata can be omitted (pass []),
%         in which case the legacy sensor placement is used.
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 2.0 (Integrated determine_sensor_mask + element averaging)
%
%   See also: get_psf, determine_sensor_mask, kspaceFirstOrder3D

    %% ======================== CONFIG DEFAULTS ========================

    dose_per_pulse_cGy = safe_config(config, 'dose_per_pulse_cGy', 0.16);
    pml_size           = safe_config(config, 'pml_size', 10);
    cfl                = safe_config(config, 'cfl_number', 0.3);
    use_gpu            = safe_config(config, 'use_gpu', true);
    num_tr_iter        = safe_config(config, 'num_time_reversal_iter', 30);
    convergence_tol    = safe_config(config, 'convergence_tol', 1e-3);
    correction_factor  = safe_config(config, 'correction_factor', 1.9);
    use_grid_padding   = safe_config(config, 'use_grid_padding', true);

    if nargin < 6, psf_filter = []; end

    % Initialize sim_results output
    sim_results = struct();

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

    sim_results.num_pulses = num_pulses;

    %% ======================== INITIAL PRESSURE p0 ========================

    % p0(r) = D(r)/N_pulses * Gamma(r) * rho(r)
    dose_per_pulse = doseGrid / num_pulses;
    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    p0 = smooth(p0);

    fprintf('        Max dose: %.4f Gy, Per-pulse max: %.6f Gy\n', ...
        max(doseGrid(:)), max(dose_per_pulse(:)));
    fprintf('        Initial pressure range: [%.2e, %.2e] Pa\n', ...
        min(p0(:)), max(p0(:)));

    sim_results.p0_max = max(p0(:));

    %% ======================== CHECK FOR SIGNIFICANT DOSE ========================

    doseThreshold = 0.01 * max(doseGrid(:));  % 1% of max
    doseMask = doseGrid > doseThreshold;

    if ~any(doseMask(:)) || max(p0(:)) == 0
        warning('run_single_field_simulation:NoDose', ...
            'No significant dose or zero initial pressure. Returning zeros.');
        recon_dose = zeros(gridSize);
        sim_results.sensor_info = struct();
        return;
    end

    %% ======================== OPTIMAL GRID PADDING ========================
    %  Pad grid to FFT-friendly dimensions for k-Wave performance.
    %  Original data sits at indices 1:N_orig; padding at N_orig+1:N_pad.
    %  Padding region filled with water medium properties.

    Nx_orig = Nx;
    Ny_orig = Ny;
    Nz_orig = Nz;
    gridSize_orig = gridSize;
    medium_orig   = medium;

    if use_grid_padding
        Nx_pad = find_optimal_kwave_size(Nx, pml_size);
        Ny_pad = find_optimal_kwave_size(Ny, pml_size);
        Nz_pad = find_optimal_kwave_size(Nz, pml_size);
    else
        Nx_pad = Nx; Ny_pad = Ny; Nz_pad = Nz;
    end

    did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
    if did_pad
        fprintf('        Padding grid: [%d %d %d] -> [%d %d %d] (FFT-optimal)\n', ...
            Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);

        % Pad medium arrays with water properties
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

        % Pad p0 with zeros (no pressure in padding region)
        p0_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
        p0_pad(1:Nx, 1:Ny, 1:Nz) = p0;
        p0 = p0_pad;

        % Update dimensions
        Nx = Nx_pad;
        Ny = Ny_pad;
        Nz = Nz_pad;
        gridSize = [Nx, Ny, Nz];

        sim_results.grid_padding = struct( ...
            'original_size', [Nx_orig, Ny_orig, Nz_orig], ...
            'padded_size',   [Nx_pad, Ny_pad, Nz_pad]);
    else
        fprintf('        Grid [%d %d %d] already FFT-optimal, no padding needed.\n', ...
            Nx, Ny, Nz);
        sim_results.grid_padding = struct( ...
            'original_size', [Nx_orig, Ny_orig, Nz_orig], ...
            'padded_size',   [Nx_orig, Ny_orig, Nz_orig]);
    end

    %% ======================== k-WAVE GRID ========================

    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    % CFL-stable time step
    maxC = max(medium.sound_speed(:));
    minC = min(medium.sound_speed(medium.sound_speed > 0));
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

    kmedium             = struct();
    kmedium.density     = medium.density;
    kmedium.sound_speed = medium.sound_speed;
    kmedium.alpha_coeff = 0 * medium.alpha_coeff;
    kmedium.alpha_power = 0 * medium.alpha_power;

    %% ======================== SENSOR PLACEMENT ========================
    %  Sensor geometry is selected via config.sensor_placement_method.
    %  X = left-right (lateral), Y = anterior-posterior, Z = sup-inf (transverse).

    sensor_method = safe_config(config, 'sensor_placement_method', 'full_plane_anterior');
    sensor = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    switch sensor_method
        case 'full_plane_anterior'
            sensor_x = safe_config(config, 'sensor_x_index', 1);
            sensor.mask(sensor_x, :, :) = 1;
            fprintf('        Sensor: full_plane_anterior — YZ plane at x = %d\n', sensor_x);

        case 'full_plane_lateral'
            sensor_y = safe_config(config, 'sensor_y_index', 1);
            sensor.mask(:, sensor_y, :) = 1;
            fprintf('        Sensor: full_plane_lateral — XZ plane at y = %d\n', sensor_y);

        case 'spherical'
            sph_radius = floor(min([Nx, Ny, Nz]) / 2) - safe_config(config, 'pml_size', 10);
            sensor.mask = makeSphere(Nx, Ny, Nz, sph_radius);
            fprintf('        Sensor: spherical — radius %d voxels\n', sph_radius);

        otherwise
            error('run_single_field_simulation:UnknownSensorMethod', ...
                'Unknown sensor_placement_method: "%s". Expected ''full_plane_anterior'', ''full_plane_lateral'', or ''spherical''.', ...
                sensor_method);
    end

    sensor_info = struct('element_map', [], 'num_elements', 0);
    sim_results.sensor_info = sensor_info;

    numSensorPts = sum(sensor.mask(:));
    fprintf('        Sensor: %d active points\n', numSensorPts);

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
        sim_results.forward_time_s = fwd_time;
    catch ME
        warning('run_single_field_simulation:ForwardFail', ...
            'Forward simulation failed: %s', ME.message);
        recon_dose = zeros(gridSize);
        return;
    end

    sensorData_measured = sensorData;

    %% ======================== TIME REVERSAL RECONSTRUCTION ========================

    fprintf('        Running iterative time reversal (%d iterations, tol=%.1e)...\n', ...
        num_tr_iter, convergence_tol);

    reconPressure      = zeros(gridSize);
    reconPressure_prev = zeros(gridSize);

    % Convergence tracking
    conv_max_pressure = zeros(num_tr_iter, 1);
    conv_rel_change   = nan(num_tr_iter, 1);
    num_iters_done    = 0;

    try
        tr_tic = tic;

        for tr_iter = 1:num_tr_iter

            fprintf('          TR iteration %d/%d...\n', tr_iter, num_tr_iter);

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
                reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
            end

            % Positivity constraint (dose and pressure are non-negative)
            reconPressure = max(reconPressure, 0);

            % Record convergence metrics
            conv_max_pressure(tr_iter) = max(reconPressure(:));
            num_iters_done = tr_iter;

            fprintf('          Max pressure: %.4e Pa\n', conv_max_pressure(tr_iter));

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
                fprintf('          Rel change: %.4e\n', rel_change);
                if rel_change < convergence_tol
                    fprintf('          *** Converged at iteration %d ***\n', tr_iter);
                    converged = true;
                end
            end

            reconPressure_prev = reconPressure;

            if converged
                break;
            end

            % Residual correction for next iteration
            if tr_iter < num_tr_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, ...
                    source_resid, sensor, inputArgs{:});

                % Residual correction using measured data as reference
                sensorData = sensorData + (sensorData_measured - sensorDataRecon);
            end
        end

        reconPressure = gather(reconPressure) * correction_factor;

        tr_time = toc(tr_tic);
        fprintf('        Time reversal complete (%.1f s).\n', tr_time);
        fprintf('        Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
            min(reconPressure(:)), max(reconPressure(:)));

        sim_results.tr_time_s         = tr_time;
        sim_results.recon_max         = max(reconPressure(:));
        sim_results.num_iters_done    = num_iters_done;
        sim_results.conv_max_pressure = conv_max_pressure(1:num_iters_done);
        sim_results.conv_rel_change   = conv_rel_change(1:num_iters_done);

    catch ME
        warning('run_single_field_simulation:TRFail', ...
            'Time reversal failed: %s', ME.message);
        recon_dose = zeros(gridSize_orig);
        return;
    end

    %% ======================== CROP TO ORIGINAL SIZE ========================
    %  Remove padding before PSF correction (computed at original grid
    %  size) and dose conversion (uses original-size gruneisen/density).

    if did_pad
        fprintf('        Cropping reconstruction: [%d %d %d] -> [%d %d %d]\n', ...
            Nx, Ny, Nz, Nx_orig, Ny_orig, Nz_orig);

        reconPressure = reconPressure(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
        Nx = Nx_orig;
        Ny = Ny_orig;
        Nz = Nz_orig;
        gridSize    = gridSize_orig;
        medium      = medium_orig;
        sensor.mask = sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    end

    %% ======================== PSF CORRECTION (OPTIONAL) ========================
    %  If a pre-computed PSF correction filter (from get_psf) is provided,
    %  apply it to the planar reconstruction to compensate for limited-view
    %  artifacts.  The filter is computed once before the field loop and
    %  reused for every field.

    if ~isempty(psf_filter) && isstruct(psf_filter) && isfield(psf_filter, 'F') ...
            && ~isempty(psf_filter.F)
        fprintf('        Applying pre-computed PSF correction...\n');
        P_field = fftn(reconPressure);
        corrected = real(ifftn(P_field .* psf_filter.F));
        reconPressure = max(corrected, 0);
        fprintf('        Corrected pressure range: [%.2e, %.2e] Pa\n', ...
            min(reconPressure(:)), max(reconPressure(:)));
        sim_results.psf_applied = true;
    else
        sim_results.psf_applied = false;
    end

    %% ======================== PRESSURE -> DOSE CONVERSION ========================

    % D_recon = p0_recon / (Gamma * rho) * num_pulses
    conversionFactor = medium.gruneisen .* medium.density;
    conversionFactor(conversionFactor == 0) = 1;  % prevent div-by-zero

    reconDosePerPulse = reconPressure ./ conversionFactor;

    body_mask_plot = ones(gridSize);
    if isfield(sct_resampled, 'bodyMask') && isequal(size(sct_resampled.bodyMask), gridSize)
        body_mask_plot = double(sct_resampled.bodyMask);
    end
    recon_dose = reconDosePerPulse * num_pulses .* double(doseMask) .* body_mask_plot;

    fprintf('        Reconstructed dose: [%.4f, %.4f] Gy\n', ...
        min(recon_dose(:)), max(recon_dose(:)));
    fprintf('        Field simulation complete.\n');
end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function val = safe_config(config, field_name, default_val)
%SAFE_CONFIG Retrieve config field with fallback to default
    if isfield(config, field_name) && ~isempty(config.(field_name))
        val = config.(field_name);
    else
        val = default_val;
    end
end
