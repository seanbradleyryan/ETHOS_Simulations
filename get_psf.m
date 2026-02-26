function psf_filter = get_psf(total_rs_dose, sct_resampled, medium, config)
%GET_PSF Compute limited-view PSF correction filter via spherical reference
%
%   psf_filter = get_psf(total_rs_dose, sct_resampled, medium, config)
%
%   PURPOSE:
%   Computes a Wiener-regularized frequency-domain correction filter that
%   compensates for limited-angle artifacts in planar sensor reconstruction.
%   The filter is computed ONCE using the total dose as a calibration source,
%   then passed to every per-field run_single_field_simulation call.
%
%   ALGORITHM:
%     1. Compute dose centroid from total_rs_dose (weighted centre of mass)
%     2. Place a 5-voxel-radius ball at the centroid via makeBall, then
%        smooth it with the k-Wave smooth() function (Hann window) to give
%        a compact, well-conditioned calibration source.  Using a simple
%        ball rather than the full dose avoids the large amplitude /
%        complex-spatial-structure issues that cause simulation divergence.
%     3. Forward simulate with planar sensor, time-reverse -> p0_planar
%     4. Forward simulate with spherical sensor (makeSphere), time-reverse
%        -> p0_sphere
%     5. Compute Wiener filter in Fourier space:
%          F(k) = S(k) * conj(P(k)) / (|P(k)|^2 + lambda^2 * mean(|P|^2))
%        where P = FFT3(p0_planar), S = FFT3(p0_sphere)
%     6. Return F for application to per-field planar reconstructions
%
%   The filter characterises the geometry-dependent frequency loss of the
%   limited planar aperture and is scale-invariant (the calibration source
%   magnitude does not affect the result).
%
%   INPUTS:
%       total_rs_dose   - 3D total dose array (Gy) from step15_process_doses.
%                         Used ONLY to locate the dose centroid; the actual
%                         dose values are not used as the calibration source.
%       sct_resampled   - Struct with .spacing, .cubeHU, .bodyMask, etc.
%       medium          - Acoustic medium struct from create_acoustic_medium:
%                         .density, .sound_speed, .alpha_coeff, .alpha_power,
%                         .gruneisen, .grid_size
%       config          - Configuration struct:
%           .pml_size               - PML thickness in voxels (default: 10)
%           .cfl_number             - CFL stability number (default: 0.3)
%           .use_gpu                - Boolean (default: true)
%           .regularization_lambda  - Wiener regularization (default: 0.01)
%           .sensor_x_index         - X index for planar sensor (default: 1)
%
%   OUTPUTS:
%       psf_filter - Struct containing:
%           .F                      - 3D complex Wiener filter (Fourier domain)
%           .grid_size              - [Nx, Ny, Nz]
%           .regularization_lambda  - Lambda used
%           .computation_time_s     - Total wall time (seconds)
%           .planar_sensor_pts      - Number of planar sensor voxels
%           .sphere_sensor_pts      - Number of spherical sensor voxels
%           .sphere_radius          - Radius used for makeSphere
%           .ball_centroid          - [cx, cy, cz] voxel indices of ball centre
%           .ball_radius            - Radius of calibration ball (voxels)
%
%   EXAMPLE:
%       psf = get_psf(total_rs_dose, sct_resampled, medium, CONFIG);
%       [recon, results] = run_single_field_simulation(...
%           field_dose, sct_resampled, medium, beam_meta, CONFIG, psf);
%
%   DEPENDENCIES:
%       - k-Wave Toolbox (kspaceFirstOrder3D, kWaveGrid, makeSphere)
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: February 2026
%   VERSION: 1.0
%
%   See also: run_single_field_simulation, makeSphere

    total_tic = tic;

    %% ======================== CONFIG DEFAULTS ========================

    pml_size   = safe_config(config, 'pml_size', 10);
    cfl        = safe_config(config, 'cfl_number', 0.3);
    use_gpu    = safe_config(config, 'use_gpu', true);
    reg_lambda = safe_config(config, 'regularization_lambda', 0.01);
    sensor_x   = safe_config(config, 'sensor_x_index', 1);

    %% ======================== GRID SETUP ========================

    gridSize = size(total_rs_dose);
    Nx = gridSize(1);
    Ny = gridSize(2);
    Nz = gridSize(3);

    spacing_mm = sct_resampled.spacing(:)';
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;

    fprintf('[PSF] Grid: [%d x %d x %d], spacing: [%.3f, %.3f, %.3f] mm\n', ...
        Nx, Ny, Nz, dx*1000, dy*1000, dz*1000);

    %% ======================== DOSE CENTROID ========================
    % Use total_rs_dose only to locate the centre of mass of the dose
    % distribution.  The dose values themselves are not used as p0 because
    % the large amplitudes and complex spatial structure cause the forward
    % simulation to diverge.

    doseSum = sum(total_rs_dose(:));
    if doseSum == 0
        warning('get_psf:ZeroDose', ...
            'total_rs_dose is all zero. Returning empty filter.');
        psf_filter = struct('F', [], 'grid_size', gridSize);
        return;
    end

    [xg, yg, zg] = ndgrid(1:Nx, 1:Ny, 1:Nz);
    cx = round(sum(xg(:) .* total_rs_dose(:)) / doseSum);
    cy = round(sum(yg(:) .* total_rs_dose(:)) / doseSum);
    cz = round(sum(zg(:) .* total_rs_dose(:)) / doseSum);

    % Clamp centroid to valid interior (keep away from PML boundary)
    margin = pml_size + 6;  % ball radius (5) + 1
    cx = max(margin, min(Nx - margin + 1, cx));
    cy = max(margin, min(Ny - margin + 1, cy));
    cz = max(margin, min(Nz - margin + 1, cz));

    fprintf('[PSF] Dose centroid (voxel): [%d, %d, %d]\n', cx, cy, cz);

    %% ======================== BALL CALIBRATION SOURCE ========================
    % Place a compact 5-voxel-radius ball at the centroid and smooth it.
    % This gives a well-conditioned, amplitude-controlled calibration source
    % whose PSF is purely geometry-dependent (independent of dose magnitude).

    ball_radius = 5;
    p0_ball = makeBall(Nx, Ny, Nz, cx, cy, cz, ball_radius);
    p0      = smooth(p0_ball, true);          % k-Wave Hann-window smoothing
    p0      = double(p0);                     % ensure double for kspaceFirstOrder3D

    fprintf('[PSF] Ball source: centre [%d,%d,%d], radius %d voxels\n', ...
        cx, cy, cz, ball_radius);
    fprintf('[PSF] p0 range after smoothing: [%.2e, %.2e]\n', ...
        min(p0(:)), max(p0(:)));
    %% ======================== OPTIMAL GRID PADDING ========================
    %  Pad grid to FFT-friendly dimensions for k-Wave performance.
    %  Original data sits at indices 1:N_orig; padding at N_orig+1:N_pad.
    %  Padding region filled with water medium properties.

    Nx_orig = Nx;  Ny_orig = Ny;  Nz_orig = Nz;
    gridSize_orig = gridSize;

    Nx_pad = find_optimal_kwave_size(Nx, pml_size);
    Ny_pad = find_optimal_kwave_size(Ny, pml_size);
    Nz_pad = find_optimal_kwave_size(Nz, pml_size);

    did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
    if did_pad
        fprintf('[PSF] Padding grid: [%d %d %d] -> [%d %d %d] (FFT-optimal)\n', ...
            Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);

        density_pad    = ones(Nx_pad, Ny_pad, Nz_pad) * 1000;
        soundSpeed_pad = ones(Nx_pad, Ny_pad, Nz_pad) * 1540;
        density_pad(1:Nx, 1:Ny, 1:Nz)    = medium.density;
        soundSpeed_pad(1:Nx, 1:Ny, 1:Nz) = medium.sound_speed;

        if numel(medium.alpha_coeff) > 1
            alphaCoeff_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
            alphaCoeff_pad(1:Nx, 1:Ny, 1:Nz) = medium.alpha_coeff;
        else
            alphaCoeff_pad = medium.alpha_coeff;   % scalar — no padding needed
        end

        p0_pad = zeros(Nx_pad, Ny_pad, Nz_pad);
        p0_pad(1:Nx, 1:Ny, 1:Nz) = p0;
        p0 = p0_pad;

        Nx = Nx_pad;  Ny = Ny_pad;  Nz = Nz_pad;
        gridSize = [Nx, Ny, Nz];
    else
        fprintf('[PSF] Grid [%d %d %d] already FFT-optimal, no padding needed.\n', Nx, Ny, Nz);
        density_pad    = medium.density;
        soundSpeed_pad = medium.sound_speed;
        alphaCoeff_pad = medium.alpha_coeff;
    end

    %% ======================== k-WAVE GRID & TIMING ========================

    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    maxC = max(medium.sound_speed(:));
    minC = min(medium.sound_speed(medium.sound_speed > 0));
    dt   = cfl * min([dx, dy, dz]) / maxC;

    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);

    kgrid.dt = dt;
    kgrid.Nt = Nt;

    fprintf('[PSF] dt = %.2e s, Nt = %d, T_sim = %.2e s\n', dt, Nt, simTime);

    %% ======================== MEDIUM ========================

    kmedium = struct();
    kmedium.density     = density_pad;
    kmedium.sound_speed = soundSpeed_pad;
    kmedium.alpha_coeff = alphaCoeff_pad;
    kmedium.alpha_power = medium.alpha_power;   % scalar, unchanged

    %% ======================== DATA CAST ========================

    if use_gpu
        try
            gpuDevice;
            dataCast = 'gpuArray-single';
            fprintf('[PSF] Compute: GPU\n');
        catch
            dataCast = 'single';
            fprintf('[PSF] Compute: CPU (GPU unavailable)\n');
        end
    else
        dataCast = 'single';
        fprintf('[PSF] Compute: CPU\n');
    end

    inputArgs = {'Smooth', false, ...
                 'PMLInside', false, ...
                 'PMLSize', pml_size, ...
                 'DataCast', dataCast, ...
                 'PlotSim', false};

    source_fwd    = struct();
    source_fwd.p0 = p0;

    %% ======================== PLANAR RECONSTRUCTION ========================

    fprintf('[PSF] Step 1/4: Forward simulation (planar sensor at x=%d)...\n', sensor_x);

    sensor_planar = struct();
    sensor_planar.mask = zeros(Nx, Ny, Nz);
    sensor_planar.mask(sensor_x, :, :) = 1;
    planarPts = sum(sensor_planar.mask(:));

    planar_fwd_tic = tic;
    sensorData_planar = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, ...
        sensor_planar, inputArgs{:});
    planar_fwd_time = toc(planar_fwd_tic);
    fprintf('[PSF]   Planar forward complete (%.1f s). Data: [%d x %d]\n', ...
        planar_fwd_time, size(sensorData_planar, 1), size(sensorData_planar, 2));

    fprintf('[PSF] Step 2/4: Time reversal (planar)...\n');
    planar_tr_tic = tic;
    p0_planar = run_single_tr(kgrid, kmedium, sensor_planar, ...
        sensorData_planar, gridSize, inputArgs);
    planar_tr_time = toc(planar_tr_tic);
    fprintf('[PSF]   Planar TR complete (%.1f s). Range: [%.2e, %.2e] Pa\n', ...
        planar_tr_time, min(p0_planar(:)), max(p0_planar(:)));

    %% ======================== SPHERICAL RECONSTRUCTION ========================

    fprintf('[PSF] Step 3/4: Forward simulation (spherical sensor)...\n');

    radius = floor(min([Nx, Ny, Nz]) / 2) - 1;
    sensor_sphere = struct();
    sensor_sphere.mask = makeSphere(Nx, Ny, Nz, radius);
    spherePts = sum(sensor_sphere.mask(:));

    fprintf('[PSF]   Spherical sensor: radius = %d, %d points\n', radius, spherePts);

    sphere_fwd_tic = tic;
    sensorData_sphere = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, ...
        sensor_sphere, inputArgs{:});
    sphere_fwd_time = toc(sphere_fwd_tic);
    fprintf('[PSF]   Sphere forward complete (%.1f s). Data: [%d x %d]\n', ...
        sphere_fwd_time, size(sensorData_sphere, 1), size(sensorData_sphere, 2));

    fprintf('[PSF] Step 4/4: Time reversal (spherical)...\n');
    sphere_tr_tic = tic;
    p0_sphere = run_single_tr(kgrid, kmedium, sensor_sphere, ...
        sensorData_sphere, gridSize, inputArgs);
    sphere_tr_time = toc(sphere_tr_tic);
    fprintf('[PSF]   Sphere TR complete (%.1f s). Range: [%.2e, %.2e] Pa\n', ...
        sphere_tr_time, min(p0_sphere(:)), max(p0_sphere(:)));

    %% ======================== CROP TO ORIGINAL SIZE ========================
    %  Remove padding so the Wiener filter is computed at the original grid
    %  size, matching the size expected by run_single_field_simulation.

    if did_pad
        fprintf('[PSF] Cropping TR results: [%d %d %d] -> [%d %d %d]\n', ...
            Nx, Ny, Nz, Nx_orig, Ny_orig, Nz_orig);
        p0_planar = p0_planar(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
        p0_sphere = p0_sphere(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
        gridSize  = gridSize_orig;
    end

    %% ======================== WIENER FILTER ========================

    fprintf('[PSF] Computing Wiener compensation filter (lambda = %.4f)...\n', reg_lambda);

    P = fftn(p0_planar);
    S = fftn(p0_sphere);

    P_abs_sq  = abs(P).^2;
    noise_reg = reg_lambda^2 * mean(P_abs_sq(:));

    F = S .* conj(P) ./ (P_abs_sq + noise_reg);

    %% ======================== PACK OUTPUT ========================

    total_time = toc(total_tic);

    psf_filter = struct();
    psf_filter.F                      = F;
    psf_filter.grid_size              = gridSize;
    psf_filter.regularization_lambda  = reg_lambda;
    psf_filter.computation_time_s     = total_time;
    psf_filter.planar_sensor_pts      = planarPts;
    psf_filter.sphere_sensor_pts      = spherePts;
    psf_filter.sphere_radius          = radius;
    psf_filter.ball_centroid          = [cx, cy, cz];
    psf_filter.ball_radius            = ball_radius;

    fprintf('[PSF] Complete (%.1f s). Filter dynamic range: [%.2e, %.2e]\n', ...
        total_time, min(abs(F(:))), max(abs(F(:))));

end


%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function p0_recon = run_single_tr(kgrid, kmedium, sensor, sensorData, gridSize, inputArgs)
%RUN_SINGLE_TR Single-iteration Dirichlet time reversal
%
%   Performs one pass of time reversal with Dirichlet boundary conditions.
%   Used for PSF calibration where a single iteration is sufficient.

    Nx = gridSize(1);
    Ny = gridSize(2);
    Nz = gridSize(3);

    source_tr        = struct();
    source_tr.p_mask = sensor.mask;
    source_tr.p      = fliplr(sensorData);
    source_tr.p_mode = 'dirichlet';

    sensor_tr        = struct();
    sensor_tr.mask   = ones(Nx, Ny, Nz);
    sensor_tr.record = {'p_final'};

    p0_raw = kspaceFirstOrder3D(kgrid, kmedium, source_tr, sensor_tr, inputArgs{:});

    if isstruct(p0_raw) && isfield(p0_raw, 'p_final')
        p0_recon = reshape(p0_raw.p_final, [Nx, Ny, Nz]);
    else
        p0_recon = reshape(p0_raw, [Nx, Ny, Nz]);
    end

    p0_recon = max(p0_recon, 0);  % positivity constraint
end


function val = safe_config(config, field_name, default_val)
%SAFE_CONFIG Retrieve config field with fallback to default
    if isfield(config, field_name) && ~isempty(config.(field_name))
        val = config.(field_name);
    else
        val = default_val;
    end
end
