function [recon_dose, sim_results] = run_single_field_simulation(field_dose, cbct_resampled, medium, beam_metadata, config, precomputed_sensor, medium_recon)
%RUN_SINGLE_FIELD_SIMULATION k-Wave forward + time-reversal for one field
%
%   [recon_dose, sim_results] = run_single_field_simulation(field_dose, cbct_resampled, medium, beam_metadata, config, precomputed_sensor, medium_recon)
%
%   Converts a single radiation field dose to initial acoustic pressure
%   (p0 = D * Gamma * rho), runs the k-Wave forward simulation to generate
%   synthetic sensor data, then applies time-reversal (or DAS) reconstruction
%   to recover the initial pressure distribution and converts back to dose.
%
%   TWO RECONSTRUCTION MODES (selected by CONFIG.blind_recon):
%     * FULL ACCESS (CONFIG.blind_recon = false, the default): the forward
%       simulation AND the reconstruction both run on `medium`. Identical to the
%       historical behavior of this function; `medium_recon` is ignored.
%     * BLIND (CONFIG.blind_recon = true AND a non-empty `medium_recon`): models
%       the live IRAI case where the acoustic waves are generated in and travel
%       through the patient's TRUE (adapted, e.g. CT_3) anatomy `medium`, but the
%       ONLY geometry we can invert on is the reference (e.g. CT_1) anatomy
%       `medium_recon`. Forward propagation (and p0) use `medium`; iterative time
%       reversal, the residual re-forward inside the TR loop, DAS delays, and the
%       final pressure->dose conversion all use `medium_recon`. Both media live on
%       the same base grid, so every grid transform (force-uniform overrides,
%       downscale, FFT padding, determine_sensor_mask grid expansion / re-padding)
%       is applied to `medium_recon` in lockstep to keep the two voxel-aligned.
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
%           .origin        - [x, y, z] mm (optional; used by determine_sensor_mask)
%       cbct_resampled - Per-CBCT geometry struct (CBCT1_resampled or
%           CBCT3_resampled from step15_process_doses). Fields used here:
%           .bodyMask      - 3D logical (true = inside body)
%           .couchMask     - 3D logical (true = couch region) [optional]
%           .spacing       - [dx, dy, dz] in mm
%           .dimensions    - [nx, ny, nz]
%           .origin        - [x, y, z] in mm
%           .cubeHU        - 3D HU array [optional; used by determine_sensor_mask]
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
%       config - Configuration struct (defaults shown in []):
%           .dose_per_pulse_cGy         [0.16]   - Dose per LINAC pulse (cGy)
%           .pml_size                   [10]     - PML thickness (voxels)
%           .cfl_number                 [0.3]    - CFL stability number
%           .use_gpu                    [true]   - GPU acceleration
%           .num_time_reversal_iter     [30]     - Max TR iterations
%           .convergence_tol            [1e-3]   - Relative-change TR early exit
%           .correction_factor          [1.9]    - Scalar multiplier on recon pressure
%           .use_grid_padding           [true]   - Pad grid to FFT-optimal dims
%           .sensor_placement_method    ['full_plane_anterior']
%                                                  Options: 'full_plane_anterior',
%                                                  'full_plane_lateral', 'spherical',
%                                                  'box', 'determine_sensor_mask',
%                                                  'fixed_anterior'.
%           .sensor_x_index             [1]      - YZ plane x-index
%           .sensor_y_index             [1]      - XZ plane y-index
%           .convolution_kernel         [4e-6]   - Gaussian pulse sigma (s); 0 disables
%           .conv_noise_level           [0.01]   - Noise fraction of peak sensor
%           .conv_deconv_lambda         [1e-4]   - Wiener deconvolution regularisation
%           .noise_only                 [false]  - Null-signal control: run the
%                                                  forward sim to fix the noise
%                                                  amplitude, then discard the true
%                                                  signal so ONLY electronic noise
%                                                  is reconstructed. Whatever dose
%                                                  comes out is pure noise artifact,
%                                                  a baseline for "is the recon
%                                                  actually recovering signal?".
%           .reconstruction_method      ['tr']   - 'tr' or 'das'
%           .blind_recon                [false]  - When true (and medium_recon is
%                                                  supplied), forward on `medium`
%                                                  but reconstruct on medium_recon.
%           .use_attenuation            [true]   - Use medium.alpha_coeff in kspaceFirstOrder3D
%           .alpha_power                [1.1]    - Acoustic attenuation power exponent
%           .use_pressure_scale_correction [false]
%                                                  Rescale reconPressure by max(p0)/max(recon)
%           .normalize                  [false]  - Divide doseGrid and recon_dose by their max
%           .mask_recon_to_dose_region  [true]   - Zero recon dose outside the
%                                                  >1% dose mask. False keeps the
%                                                  full reconstruction.
%           .downscale_factor           [1]      - imresize3 grids by 1/df before sim
%           .force_uniform_density      [false]
%           .force_uniform_sound_speed  [false]
%           .force_uniform_attenuation  [false]
%           .force_uniform_gruneisen    [false]
%           .uniform_density            [1000]
%           .uniform_sound_speed        [1540]
%           .uniform_alpha_coeff        [0]
%           .uniform_gruneisen          [1.0]
%           .plot_results               [false]  - Opt-in: pre-sim sensor view, live TR
%                                                  figure, post-sim dose panels +
%                                                  convergence plot. Must remain false
%                                                  when called under parfor.
%       precomputed_sensor - (Optional) Struct holding a plan-level sensor mask
%           computed ONCE by the caller (pipeline_simulate) and reused for every
%           field, so determine_sensor_mask is not re-run per field. Fields:
%           .sensor_mask  - Sensor mask returned by determine_sensor_mask, on
%                           the original (possibly grid-expanded) dose grid.
%           .sensor_info  - Companion diagnostics struct (carries .grid_pad,
%                           which drives the grid-expansion handling below).
%           Honored ONLY for sensor_placement_method='determine_sensor_mask',
%           and only when its base grid matches this field's grid (e.g. it is
%           ignored under a non-unity downscale_factor). Pass [] or omit to
%           compute the mask inline per field as before.
%       medium_recon - (Optional) RECONSTRUCTION (reference-geometry) acoustic
%           medium from create_acoustic_medium() for the reference CBCT (e.g.
%           CBCT1). Used for iterative time reversal, the residual re-forward
%           inside the TR loop, DAS delays, and the pressure->dose conversion.
%           Same fields as `medium` and the SAME base grid. Honored only when
%           CONFIG.blind_recon is true; pass [] or omit for full-access
%           reconstruction (both passes on `medium`).
%
%   OUTPUTS:
%       recon_dose  - 3D reconstructed dose array (Gy), cropped back to the
%                     input field_dose.dose_Gy size (or, with downscale_factor,
%                     to the downscaled size).
%       sim_results - Struct with diagnostics:
%           .sensor_info        - Sensor placement info
%           .forward_time_s     - Forward sim wall time
%           .tr_time_s          - Reconstruction wall time (TR or DAS)
%           .num_pulses         - Number of LINAC pulses
%           .p0_max             - Maximum initial pressure (Pa)
%           .recon_max          - Maximum reconstructed pressure (Pa)
%           .num_iters_done     - Iterations completed
%           .conv_max_pressure  - Per-iteration max pressure
%           .conv_rel_change    - Per-iteration relative change
%           .grid_padding       - Original / padded / expanded sizes
%           .recon_method       - Which branch ran ('tr' | 'das')
%           .blind_recon        - Logical: true when forward and reconstruction
%                                 media differed (blind geometry)
%
%   NOTES:
%       - Stateless: safe for parfor execution. Diagnostic plotting must be
%         disabled under parfor (figures from workers do not render).
%       - PlotSim disabled inside kspaceFirstOrder3D in all cases.
%       - Returns zeros if no significant dose or on simulation failure.
%       - beam_metadata is the FULL plan metadata (all beams), passed through
%         to determine_sensor_mask for computing the exclusion zone.
%       - For backwards compatibility, beam_metadata can be omitted (pass []),
%         in which case the legacy sensor placement is used.
%       - precomputed_sensor lets the caller compute the determine_sensor_mask
%         placement once for a whole plan and reuse it for every field. When it
%         is absent or its grid does not match, the mask is computed inline.
%
%   AUTHOR: ETHOS Pipeline Team
%   DATE: May 2026
%   VERSION: 3.0 (DAS, determine_sensor_mask + grid expansion, scale correction,
%                 dose normalisation, downscale_factor, attenuation toggle,
%                 force-uniform medium overrides, opt-in plotting)
%
%   See also: determine_sensor_mask, determine_sensor_placement_fixed,
%             kspaceFirstOrder3D, run_standalone_simulation

    %% ======================== CONFIG DEFAULTS ========================

    dose_per_pulse_cGy  = safe_config(config, 'dose_per_pulse_cGy', 0.16);
    pml_size            = safe_config(config, 'pml_size', 10);
    cfl                 = safe_config(config, 'cfl_number', 0.3);
    use_gpu             = safe_config(config, 'use_gpu', true);
    num_tr_iter         = safe_config(config, 'num_time_reversal_iter', 30);
    convergence_tol     = safe_config(config, 'convergence_tol', 1e-3);
    correction_factor   = safe_config(config, 'correction_factor', 1.9);
    use_grid_padding    = safe_config(config, 'use_grid_padding', true);
    conv_kernel_sigma   = safe_config(config, 'convolution_kernel', 4e-6);
    conv_noise_level    = safe_config(config, 'conv_noise_level', 0.01);
    conv_deconv_lambda  = safe_config(config, 'conv_deconv_lambda', 1e-4);
    noise_only          = safe_config(config, 'noise_only', false);

    recon_method            = lower(safe_config(config, 'reconstruction_method', 'tr'));
    use_pressure_scale_corr = safe_config(config, 'use_pressure_scale_correction', false);
    normalize_doses         = safe_config(config, 'normalize', false);
    mask_recon_to_dose      = safe_config(config, 'mask_recon_to_dose_region', true);
    use_attenuation         = safe_config(config, 'use_attenuation', true);
    alpha_power_value       = safe_config(config, 'alpha_power', 1.1);
    downscale_factor        = safe_config(config, 'downscale_factor', 1);

    force_uniform_density   = safe_config(config, 'force_uniform_density', false);
    force_uniform_speed     = safe_config(config, 'force_uniform_sound_speed', false);
    force_uniform_atten     = safe_config(config, 'force_uniform_attenuation', false);
    force_uniform_gruneisen = safe_config(config, 'force_uniform_gruneisen', false);
    uniform_density_val     = safe_config(config, 'uniform_density', 1000);
    uniform_speed_val       = safe_config(config, 'uniform_sound_speed', 1540);
    uniform_alpha_coeff_val = safe_config(config, 'uniform_alpha_coeff', 0);
    uniform_gruneisen_val   = safe_config(config, 'uniform_gruneisen', 1.0);

    plot_results = safe_config(config, 'plot_results', false);
    return_diag  = safe_config(config, 'return_diagnostics', false);

    sensor_method = safe_config(config, 'sensor_placement_method', 'full_plane_anterior');

    if nargin < 4, beam_metadata = []; end
    if nargin < 6, precomputed_sensor = []; end
    if nargin < 7, medium_recon = []; end

    % Blind reconstruction: forward-propagate on `medium` (true / CT_3 geometry)
    % but time-reverse on `medium_recon` (reference / CT_1 geometry). Enabled by
    % CONFIG.blind_recon AND a non-empty medium_recon. Otherwise it is full
    % access -- every reconstruction stage runs on `medium` exactly as before.
    blind_recon = safe_config(config, 'blind_recon', false) && ~isempty(medium_recon);
    if safe_config(config, 'blind_recon', false) && isempty(medium_recon)
        warning('run_single_field_simulation:BlindNoReconMedium', ...
            ['CONFIG.blind_recon is true but no reconstruction medium was ' ...
             'supplied; running full-access reconstruction.']);
    end

    % Initialize sim_results output
    sim_results              = struct();
    sim_results.recon_method = recon_method;
    sim_results.blind_recon  = blind_recon;
    sim_results.noise_only   = noise_only;

    %% ======================== EXTRACT DATA ========================

    doseGrid = field_dose.dose_Gy;
    gridSize = size(doseGrid);
    Nx = gridSize(1);
    Ny = gridSize(2);
    Nz = gridSize(3);

    % Grid spacing: mm -> m
    spacing_mm = cbct_resampled.spacing(:)';
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
    if blind_recon && isfield(medium_recon, 'grid_size') && ...
            ~isequal(gridSize, medium_recon.grid_size)
        error('run_single_field_simulation:ReconSizeMismatch', ...
            'Field dose [%d %d %d] does not match reconstruction medium [%d %d %d].', ...
            Nx, Ny, Nz, ...
            medium_recon.grid_size(1), medium_recon.grid_size(2), medium_recon.grid_size(3));
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
    if blind_recon
        fprintf('        [BLIND] Forward on true medium, reconstruction on reference medium.\n');
    end
    sim_results.num_pulses = num_pulses;

    %% ======================== FORCE-UNIFORM MEDIUM OVERRIDES ========================
    %  Optional sensitivity-test toggles: override one or more medium properties
    %  with uniform values. Applied to the passed-in medium in place.

    if force_uniform_density
        medium.density = ones(gridSize) * uniform_density_val;
        if blind_recon, medium_recon.density = ones(gridSize) * uniform_density_val; end
        fprintf('        [Override] density forced uniform = %g kg/m^3\n', uniform_density_val);
    end
    if force_uniform_speed
        medium.sound_speed = ones(gridSize) * uniform_speed_val;
        if blind_recon, medium_recon.sound_speed = ones(gridSize) * uniform_speed_val; end
        fprintf('        [Override] sound_speed forced uniform = %g m/s\n', uniform_speed_val);
    end
    if force_uniform_atten
        medium.alpha_coeff = ones(gridSize) * uniform_alpha_coeff_val;
        if blind_recon, medium_recon.alpha_coeff = ones(gridSize) * uniform_alpha_coeff_val; end
        fprintf('        [Override] alpha_coeff forced uniform = %g\n', uniform_alpha_coeff_val);
    end
    if force_uniform_gruneisen
        medium.gruneisen = ones(gridSize) * uniform_gruneisen_val;
        if blind_recon, medium_recon.gruneisen = ones(gridSize) * uniform_gruneisen_val; end
        fprintf('        [Override] gruneisen forced uniform = %g\n', uniform_gruneisen_val);
    end

    %% ======================== DOWNSCALE FACTOR ========================
    %  Optional pre-sim downscale of the whole grid (dose, medium, body/couch
    %  masks). Useful for rapid testing. Note that recon_dose returned from
    %  this function will then be at the downscaled size; the caller must
    %  accommodate that if it sums across fields.

    if downscale_factor ~= 1
        df     = downscale_factor;
        new_Nx = max(1, round(Nx / df));
        new_Ny = max(1, round(Ny / df));
        new_Nz = max(1, round(Nz / df));

        fprintf('        [DS] Downscaling: [%d %d %d] -> [%d %d %d]\n', ...
            Nx, Ny, Nz, new_Nx, new_Ny, new_Nz);

        doseGrid           = max(imresize3(doseGrid, [new_Nx, new_Ny, new_Nz]), 0);
        medium.density     = imresize3(medium.density,     [new_Nx, new_Ny, new_Nz]);
        medium.sound_speed = imresize3(medium.sound_speed, [new_Nx, new_Ny, new_Nz]);
        if numel(medium.alpha_coeff) > 1
            medium.alpha_coeff = imresize3(medium.alpha_coeff, [new_Nx, new_Ny, new_Nz]);
        end
        medium.gruneisen   = imresize3(medium.gruneisen,   [new_Nx, new_Ny, new_Nz]);
        medium.grid_size   = [new_Nx, new_Ny, new_Nz];

        % Blind: downscale the reconstruction medium in lockstep.
        if blind_recon
            medium_recon.density     = imresize3(medium_recon.density,     [new_Nx, new_Ny, new_Nz]);
            medium_recon.sound_speed = imresize3(medium_recon.sound_speed, [new_Nx, new_Ny, new_Nz]);
            if numel(medium_recon.alpha_coeff) > 1
                medium_recon.alpha_coeff = imresize3(medium_recon.alpha_coeff, [new_Nx, new_Ny, new_Nz]);
            end
            medium_recon.gruneisen   = imresize3(medium_recon.gruneisen,   [new_Nx, new_Ny, new_Nz]);
            medium_recon.grid_size   = [new_Nx, new_Ny, new_Nz];
        end

        if isfield(cbct_resampled, 'bodyMask') && ~isempty(cbct_resampled.bodyMask)
            cbct_resampled.bodyMask = imresize3(single(cbct_resampled.bodyMask), ...
                [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        if isfield(cbct_resampled, 'couchMask') && ~isempty(cbct_resampled.couchMask)
            cbct_resampled.couchMask = imresize3(single(cbct_resampled.couchMask), ...
                [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        if isfield(cbct_resampled, 'cubeHU') && ~isempty(cbct_resampled.cubeHU)
            cbct_resampled.cubeHU = imresize3(cbct_resampled.cubeHU, [new_Nx, new_Ny, new_Nz]);
        end

        spacing_mm = spacing_mm .* ([Nx, Ny, Nz] ./ [new_Nx, new_Ny, new_Nz]);
        dx = spacing_mm(1) / 1000;
        dy = spacing_mm(2) / 1000;
        dz = spacing_mm(3) / 1000;
        if isfield(cbct_resampled, 'spacing')
            cbct_resampled.spacing = spacing_mm;
        end

        Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
        gridSize = [Nx, Ny, Nz];
    end

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

    doseThreshold = 0.01 * max(doseGrid(:));
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

    Nx_orig       = Nx;
    Ny_orig       = Ny;
    Nz_orig       = Nz;
    gridSize_orig = gridSize;
    medium_orig   = medium;

    % Reconstruction medium at original size (pre-FFT-pad). Set only on the blind
    % path; drives the pressure->dose conversion below. Left [] for full access.
    medium_recon_orig = [];
    if blind_recon
        medium_recon_orig = medium_recon;
    end

    % Track the original (post-downscale, pre-expansion) input size and any
    % expansion offsets so the final recon_dose can be cropped back to match
    % what the caller expects from field_dose.dose_Gy.
    input_dims_for_crop = gridSize_orig;
    expansion_offsets   = [0, 0, 0];

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

        [medium, p0] = pad_medium_and_p0(medium, p0, Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);
        if blind_recon
            medium_recon = pad_medium_to(medium_recon, [Nx_pad, Ny_pad, Nz_pad]);
        end

        Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
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

    %% ======================== SENSOR PLACEMENT ========================
    %  Sensor placement happens BEFORE the k-Wave grid is built because the
    %  'determine_sensor_mask' method may expand the grid (water padding) and
    %  re-run FFT-optimal padding, which changes Nx/Ny/Nz.
    %  X = left-right (lateral), Y = anterior-posterior, Z = sup-inf.

    sensor      = struct();
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor_info = struct('element_map', [], 'num_elements', 0);

    switch sensor_method
        case 'full_plane_anterior'
            sensor_x = safe_config(config, 'sensor_x_index', 1);
            sensor.mask(sensor_x, :, :) = 1;
            fprintf('        Sensor: full_plane_anterior - YZ plane at x = %d\n', sensor_x);

        case 'full_plane_lateral'
            sensor_y = safe_config(config, 'sensor_y_index', 1);
            sensor.mask(:, sensor_y, :) = 1;
            fprintf('        Sensor: full_plane_lateral - XZ plane at y = %d\n', sensor_y);

        case 'spherical'
            sph_radius  = floor(min([Nx, Ny, Nz]) / 2) - pml_size;
            sensor.mask = makeSphere(Nx, Ny, Nz, sph_radius);
            % Zero p0 outside the enclosing ball: anything outside is
            % unobservable by this geometry and would otherwise pollute the
            % forward simulation and downstream pressure scaling.
            sph_cx = floor(Nx/2) + 1;
            sph_cy = floor(Ny/2) + 1;
            sph_cz = floor(Nz/2) + 1;
            [Xg_sph, Yg_sph, Zg_sph] = ndgrid(1:Nx, 1:Ny, 1:Nz);
            ball_mask = (Xg_sph - sph_cx).^2 + (Yg_sph - sph_cy).^2 + ...
                        (Zg_sph - sph_cz).^2 <= sph_radius^2;
            n_zeroed = nnz(p0 ~= 0 & ~ball_mask);
            p0 = p0 .* ball_mask;
            clear Xg_sph Yg_sph Zg_sph
            fprintf('        Sensor: spherical, radius %d voxels (zeroed %d p0 voxels outside)\n', ...
                sph_radius, n_zeroed);

        case 'box'
            % Six-face bounding box: planes at index 3 and (N-3) on each axis.
            bx_lo   = 3;
            bx_hi_x = Nx - 3;
            bx_hi_y = Ny - 3;
            bx_hi_z = Nz - 3;
            if bx_hi_x <= bx_lo || bx_hi_y <= bx_lo || bx_hi_z <= bx_lo
                error('run_single_field_simulation:BoxTooSmall', ...
                    'Grid [%d %d %d] too small for box sensor (need each dim > 6).', Nx, Ny, Nz);
            end
            sensor.mask(bx_lo,   bx_lo:bx_hi_y, bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_hi_x, bx_lo:bx_hi_y, bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_lo,   bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_hi_y, bx_lo:bx_hi_z) = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_lo:bx_hi_y, bx_lo)   = 1;
            sensor.mask(bx_lo:bx_hi_x, bx_lo:bx_hi_y, bx_hi_z) = 1;
            fprintf('        Sensor: box faces at x=[%d,%d], y=[%d,%d], z=[%d,%d]\n', ...
                bx_lo, bx_hi_x, bx_lo, bx_hi_y, bx_lo, bx_hi_z);

        case 'determine_sensor_mask'
            % Automatic placement via determine_sensor_mask: tilts a 2D array
            % toward the beam isocenter avoiding beam exclusion zones. May
            % expand the simulation grid in X/Y/Z with water padding to place
            % the sensor outside the exclusion zone.
            sct_for_sensor = struct();
            if isfield(cbct_resampled, 'cubeHU')
                sct_for_sensor.cubeHU = cbct_resampled.cubeHU;
            end
            if isfield(cbct_resampled, 'bodyMask')
                sct_for_sensor.bodyMask = cbct_resampled.bodyMask;
            end
            if isfield(cbct_resampled, 'couchMask') && ~isempty(cbct_resampled.couchMask)
                sct_for_sensor.couchMask = cbct_resampled.couchMask;
            else
                sct_for_sensor.couchMask = false(size(cbct_resampled.bodyMask));
            end
            if isfield(cbct_resampled, 'origin')
                sct_for_sensor.origin = cbct_resampled.origin;
            else
                sct_for_sensor.origin = [0, 0, 0];
            end
            sct_for_sensor.spacing = spacing_mm;

            field_dose_for_sensor = struct();
            field_dose_for_sensor.dose_Gy = doseGrid;
            if isfield(field_dose, 'gantry_angle')
                field_dose_for_sensor.gantry_angle = field_dose.gantry_angle;
            else
                field_dose_for_sensor.gantry_angle = 0;
            end
            if isfield(field_dose, 'origin')
                field_dose_for_sensor.origin = field_dose.origin;
            else
                field_dose_for_sensor.origin = sct_for_sensor.origin;
            end
            field_dose_for_sensor.spacing    = spacing_mm;
            field_dose_for_sensor.dimensions = [Nx_orig, Ny_orig, Nz_orig];

            if precomputed_sensor_fits(precomputed_sensor, [Nx_orig, Ny_orig, Nz_orig])
                % Reuse the plan-level sensor mask computed once by the caller
                % (deterministic placement from the summed plan dose), instead
                % of re-running determine_sensor_mask for every field/segment.
                sensor_mask_orig = precomputed_sensor.sensor_mask;
                sensor_info_orig = precomputed_sensor.sensor_info;
                fprintf('        Sensor: using precomputed plan sensor mask (%d active points)\n', ...
                    sum(sensor_mask_orig(:)));
            else
                [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask( ...
                    sct_for_sensor, field_dose_for_sensor, beam_metadata, config);
            end

            sensor_info = sensor_info_orig;

            % ---- Grid expansion ----
            % determine_sensor_mask labels dim 1=Y, dim 2=X, dim 3=Z but
            % preserves the caller's actual dim order. Since cbct_resampled
            % is passed with dim 1=Nx, dim 2=Ny, dim 3=Nz, gp.y_* maps to
            % dim 1 (Nx), gp.x_* to dim 2 (Ny), gp.z_* to dim 3 (Nz).
            if isfield(sensor_info_orig, 'grid_pad') && sensor_info_orig.grid_pad.expanded
                gp = sensor_info_orig.grid_pad;
                fprintf(['        [Sensor] Grid expansion: ' ...
                         'Y(+%d/+%d), X(+%d/+%d), Z(+%d/+%d). Re-padding with water.\n'], ...
                        gp.y_pre, gp.y_post, gp.x_pre, gp.x_post, gp.z_pre, gp.z_post);

                % Strip current FFT padding back to (Nx_orig, Ny_orig, Nz_orig)
                density_unp    = medium.density(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                soundSpeed_unp = medium.sound_speed(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                if numel(medium.alpha_coeff) > 1
                    alphaCoeff_unp = medium.alpha_coeff(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                else
                    alphaCoeff_unp = medium.alpha_coeff;
                end
                gruneisen_unp = medium.gruneisen(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                p0_unp        = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

                Nx_exp = Nx_orig + gp.y_pre + gp.y_post;
                Ny_exp = Ny_orig + gp.x_pre + gp.x_post;
                Nz_exp = Nz_orig + gp.z_pre + gp.z_post;

                density_exp    = ones(Nx_exp, Ny_exp, Nz_exp) * 1000;
                soundSpeed_exp = ones(Nx_exp, Ny_exp, Nz_exp) * 1540;
                alphaCoeff_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
                gruneisen_exp  = zeros(Nx_exp, Ny_exp, Nz_exp);
                p0_exp         = zeros(Nx_exp, Ny_exp, Nz_exp);

                xr = gp.y_pre + (1:Nx_orig);
                yr = gp.x_pre + (1:Ny_orig);
                zr = gp.z_pre + (1:Nz_orig);

                density_exp(xr, yr, zr)    = density_unp;
                soundSpeed_exp(xr, yr, zr) = soundSpeed_unp;
                if numel(alphaCoeff_unp) > 1
                    alphaCoeff_exp(xr, yr, zr) = alphaCoeff_unp;
                else
                    alphaCoeff_exp(:) = alphaCoeff_unp;
                end
                gruneisen_exp(xr, yr, zr)  = gruneisen_unp;
                p0_exp(xr, yr, zr)         = p0_unp;

                medium.density     = density_exp;
                medium.sound_speed = soundSpeed_exp;
                medium.alpha_coeff = alphaCoeff_exp;
                medium.gruneisen   = gruneisen_exp;
                p0 = p0_exp;

                % medium_orig drives the post-recon pressure->dose conversion in
                % full access; update it to the expanded-but-unpadded medium so
                % its size matches the cropped reconPressure.
                medium_orig = struct( ...
                    'density',     density_exp, ...
                    'sound_speed', soundSpeed_exp, ...
                    'alpha_coeff', alphaCoeff_exp, ...
                    'gruneisen',   gruneisen_exp, ...
                    'alpha_power', medium.alpha_power, ...
                    'grid_size',   [Nx_exp, Ny_exp, Nz_exp]);

                % Blind: mirror the strip->water-expand onto the reconstruction
                % medium so it occupies the identical expanded grid. Nx_orig/
                % Ny_orig/Nz_orig are still the PRE-expansion sizes here (bumped
                % below). medium_recon_orig drives the pressure->dose conversion,
                % so it must track the expanded-but-unpadded reconstruction medium.
                if blind_recon
                    medium_recon = crop_medium(medium_recon, [Nx_orig, Ny_orig, Nz_orig]);
                    medium_recon = expand_medium(medium_recon, [Nx_exp, Ny_exp, Nz_exp], xr, yr, zr);
                    medium_recon.alpha_power = alpha_power_value;
                    medium_recon.grid_size   = [Nx_exp, Ny_exp, Nz_exp];
                    medium_recon_orig = medium_recon;
                end

                % Expand doseGrid / doseMask / bodyMask / couchMask
                doseGrid_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
                doseGrid_exp(xr, yr, zr) = doseGrid;
                doseGrid = doseGrid_exp;

                doseMask_exp = false(Nx_exp, Ny_exp, Nz_exp);
                doseMask_exp(xr, yr, zr) = doseMask;
                doseMask = doseMask_exp;

                if isfield(cbct_resampled, 'bodyMask') && ...
                        isequal(size(cbct_resampled.bodyMask), [Nx_orig, Ny_orig, Nz_orig])
                    body_exp = false(Nx_exp, Ny_exp, Nz_exp);
                    body_exp(xr, yr, zr) = cbct_resampled.bodyMask;
                    cbct_resampled.bodyMask = body_exp;
                end
                if isfield(cbct_resampled, 'couchMask') && ~isempty(cbct_resampled.couchMask) && ...
                        isequal(size(cbct_resampled.couchMask), [Nx_orig, Ny_orig, Nz_orig])
                    couch_exp = false(Nx_exp, Ny_exp, Nz_exp);
                    couch_exp(xr, yr, zr) = cbct_resampled.couchMask;
                    cbct_resampled.couchMask = couch_exp;
                end

                % Save crop offsets BEFORE bumping Nx_orig, so the final
                % recon_dose can be sliced back to the input dose size.
                expansion_offsets = [gp.y_pre, gp.x_pre, gp.z_pre];

                Nx_orig = Nx_exp; Ny_orig = Ny_exp; Nz_orig = Nz_exp;
                gridSize_orig = [Nx_orig, Ny_orig, Nz_orig];

                % Re-run FFT-optimal padding on the expanded grid
                if use_grid_padding
                    Nx_pad2 = find_optimal_kwave_size(Nx_orig, pml_size);
                    Ny_pad2 = find_optimal_kwave_size(Ny_orig, pml_size);
                    Nz_pad2 = find_optimal_kwave_size(Nz_orig, pml_size);
                else
                    Nx_pad2 = Nx_orig; Ny_pad2 = Ny_orig; Nz_pad2 = Nz_orig;
                end

                if ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], [Nx_orig, Ny_orig, Nz_orig])
                    fprintf('        [Sensor] Re-pad to FFT-optimal: [%d %d %d] -> [%d %d %d]\n', ...
                        Nx_orig, Ny_orig, Nz_orig, Nx_pad2, Ny_pad2, Nz_pad2);
                    [medium, p0] = pad_medium_and_p0(medium, p0, ...
                        Nx_orig, Ny_orig, Nz_orig, Nx_pad2, Ny_pad2, Nz_pad2);
                    if blind_recon
                        medium_recon = pad_medium_to(medium_recon, [Nx_pad2, Ny_pad2, Nz_pad2]);
                    end
                    did_pad = true;
                else
                    did_pad = ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], gridSize_orig);
                end

                Nx = Nx_pad2; Ny = Ny_pad2; Nz = Nz_pad2;
                gridSize    = [Nx, Ny, Nz];
                sensor.mask = zeros(Nx, Ny, Nz);

                sim_results.grid_padding = struct( ...
                    'original_size', [Nx_exp, Ny_exp, Nz_exp], ...
                    'padded_size',   [Nx, Ny, Nz], ...
                    'expanded',      true, ...
                    'grid_pad',      gp);
            end

            % Embed sensor mask in current (possibly re-padded) grid
            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
            fprintf('        Sensor: determine_sensor_mask - %d active points\n', ...
                sum(sensor_mask_orig(:)));

        case 'fixed_anterior'
            % Deterministic placement: anterior, inferior to beam field,
            % centered on isocenter X.
            fixed_struct = struct();
            if ~isempty(beam_metadata)
                fixed_struct.beam_metadata = beam_metadata;
            end
            if isfield(config, 'total_dose') && ~isempty(config.total_dose)
                fixed_struct.total_dose = config.total_dose;
            else
                fixed_struct.total_dose = doseGrid;
            end
            [sensor_mask_orig, sensor_info] = determine_sensor_placement_fixed( ...
                config, cbct_resampled, fixed_struct);
            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
            fprintf('        Sensor: fixed_anterior - %d voxels (valid: %s)\n', ...
                sum(sensor_mask_orig(:)), mat2str(sensor_info.placement_valid));

        otherwise
            error('run_single_field_simulation:UnknownSensorMethod', ...
                ['Unknown sensor_placement_method: "%s". Expected ' ...
                 '''full_plane_anterior'', ''full_plane_lateral'', ''spherical'', ' ...
                 '''box'', ''determine_sensor_mask'', or ''fixed_anterior''.'], ...
                sensor_method);
    end

    sim_results.sensor_info = sensor_info;

    numSensorPts = sum(sensor.mask(:));
    fprintf('        Sensor: %d active points\n', numSensorPts);

    if numSensorPts == 0
        warning('run_single_field_simulation:EmptySensor', ...
            'Sensor mask is empty. Returning zeros.');
        recon_dose = zeros(input_dims_for_crop);
        return;
    end

    %% ======================== PRE-SIM SENSOR/DOSE PLOT (opt-in) ========================

    if plot_results
        try
            sensor_vis    = logical(sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig));
            dose_mask_vis = double(doseGrid) >= 0.10 * max(double(doseGrid(:)));
            plot_sensor_dose_planes(dose_mask_vis, sensor_vis, spacing_mm, ...
                medium_orig.density, config);
            drawnow;
        catch ME
            warning('run_single_field_simulation:PreSimPlotFail', ...
                'Pre-sim sensor/dose plot failed: %s', ME.message);
        end
    end

    %% ======================== k-WAVE GRID ========================

    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    % CFL-stable time step. The SAME kgrid.dt / Nt drive both the forward
    % propagation and (when blind) the reconstruction propagation on the
    % reference medium, so bound the step by the fastest tissue across BOTH media
    % and set the recording length by the slowest across both. For full access
    % the two media are the same, so this reduces to the historical single-medium
    % calculation.
    if blind_recon
        maxC = max(max(medium.sound_speed(:)), max(medium_recon.sound_speed(:)));
        minC = min(min(medium.sound_speed(medium.sound_speed > 0)), ...
                   min(medium_recon.sound_speed(medium_recon.sound_speed > 0)));
    else
        maxC = max(medium.sound_speed(:));
        minC = min(medium.sound_speed(medium.sound_speed > 0));
    end
    dt   = cfl * min([dx, dy, dz]) / maxC;

    % Simulation time: 2.5x grid diagonal traversal at minimum speed
    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);

    % Nt_scaling: air (~343 m/s) is the slowest tissue, so it sets minC and
    % inflates simTime/Nt. When air drives minC, shorten the recording length
    % by config.Nt_scaling (0 = disabled).
    Nt_scaling = safe_config(config, 'Nt_scaling', 0);
    if Nt_scaling ~= 0
        air_c = 343;
        if isfield(config, 'tissue_tables') && isfield(config.tissue_tables, 'threshold_2') ...
                && isfield(config.tissue_tables.threshold_2, 'air')
            air_c = config.tissue_tables.threshold_2.air.sound_speed;
        end
        if minC <= air_c
            Nt = max(1, ceil(Nt / Nt_scaling));
        end
    end

    kgrid.dt = dt;
    kgrid.Nt = Nt;

    fprintf('        Grid: [%d x %d x %d], spacing: [%.3f, %.3f, %.3f] mm\n', ...
        Nx, Ny, Nz, dx*1000, dy*1000, dz*1000);
    fprintf('        dt = %.2e s, Nt = %d, T_sim = %.2e s\n', dt, Nt, simTime);

    %% ======================== k-WAVE MEDIUM ========================

    %  kmedium     : FORWARD medium -- generates the sensor data.
    %  kmedium_rec : RECONSTRUCTION medium -- time reversal inverts on this. For
    %                full access it is the SAME struct as kmedium.

    kmedium             = struct();
    kmedium.density     = medium.density;
    kmedium.sound_speed = medium.sound_speed;
    if use_attenuation
        kmedium.alpha_coeff = medium.alpha_coeff;
        kmedium.alpha_power = alpha_power_value;
    else
        kmedium.alpha_coeff = 0 * medium.alpha_coeff;
        kmedium.alpha_power = 0;
    end

    if blind_recon
        kmedium_rec             = struct();
        kmedium_rec.density     = medium_recon.density;
        kmedium_rec.sound_speed = medium_recon.sound_speed;
        if use_attenuation
            kmedium_rec.alpha_coeff = medium_recon.alpha_coeff;
            kmedium_rec.alpha_power = alpha_power_value;
        else
            kmedium_rec.alpha_coeff = 0 * medium_recon.alpha_coeff;
            kmedium_rec.alpha_power = 0;
        end
    else
        kmedium_rec = kmedium;   % full access: reconstruct on the forward medium
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
        recon_dose = zeros(input_dims_for_crop);
        return;
    end

    % Sampling frequency for the sensor frequency-response filter, which is now
    % applied below AFTER the pulse-profile convolution (see physical chain).
    FS         = 1 / kgrid.dt;

    sensorData_measured = sensorData;

    %% ================= PULSE CONVOLUTION / FREQUENCY RESPONSE / NOISE / DECONV =================
    %  Physical measurement chain in acquisition order:
    %    1. Convolve with the pulse profile (finite radiation pulse)
    %    2. Apply the sensor frequency response (band-limit)
    %    3. Add electronic noise (downstream of the transducer, after band-limit)
    %    4. Wiener-deconvolve the pulse kernel to recover the broadband signal

    if conv_kernel_sigma > 0
        fprintf('        Pulse model: sigma=%.1f us, noise=%.1f%%, lambda=%.1e\n', ...
            conv_kernel_sigma * 1e6, conv_noise_level * 100, conv_deconv_lambda);

        sigma_samples = conv_kernel_sigma / dt;
        kernel_half   = ceil(4 * sigma_samples);
        t_kernel      = (-kernel_half : kernel_half)';
        gauss_kernel  = exp(-t_kernel.^2 / (2 * sigma_samples^2));
        gauss_kernel  = gauss_kernel / sum(gauss_kernel);

        sensorData_cpu = double(gather(sensorData));
        Nt_data        = size(sensorData_cpu, 2);

        H       = fft(gauss_kernel, Nt_data).';
        H_conj  = conj(H);
        H_power = abs(H).^2;

        % 1. Convolve with the pulse profile
        sensorData_conv = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));

        % 2. Sensor frequency response (band-limit AFTER pulse convolution)
        sensorData_resp = gaussianFilter(sensorData_conv, FS, 0.35e6, 100, true);

        % 3. Add electronic noise (after frequency-response filtering).
        %    noise_amp is set from the TRUE signal peak so a noise_only run adds
        %    the exact noise the real run would have; then, when noise_only, the
        %    true signal is discarded so ONLY noise reaches the reconstruction.
        noise_amp = conv_noise_level * max(abs(sensorData_resp(:)));
        if noise_only
            sensorData_resp(:) = 0;
            fprintf('        [NOISE ONLY] True signal nulled; reconstructing noise alone.\n');
        end
        sensorData_noisy = sensorData_resp + noise_amp * randn(size(sensorData_resp));

        % 4. Wiener deconvolution of the pulse kernel
        sensorData_deconv = real(ifft( ...
            fft(sensorData_noisy, [], 2) .* H_conj ./ (H_power + conv_deconv_lambda), ...
            [], 2));

        sensorData          = single(sensorData_deconv);
        sensorData_measured = single(sensorData_deconv);

        fprintf('        Pulse model complete. Noise amp: %.3e Pa\n', noise_amp);
    else
        % No pulse model: apply only the sensor frequency response.
        sensorData = gaussianFilter(sensorData, FS, 0.35e6, 100, true);
        if noise_only
            % Null-signal test: keep only electronic noise at the amplitude the
            % true signal would have set (see noise_only note above).
            noise_amp     = conv_noise_level * max(abs(sensorData(:)));
            sensorData(:) = 0;
            sensorData    = sensorData + noise_amp * randn(size(sensorData));
            fprintf('        [NOISE ONLY] True signal nulled; reconstructing noise alone (amp %.3e Pa).\n', ...
                noise_amp);
        end
        sensorData_measured = sensorData;
        fprintf('        Pulse convolution disabled; frequency response applied.\n');
    end

    %% ======================== RECONSTRUCTION ========================

    reconPressure     = zeros(gridSize);
    conv_max_pressure = zeros(num_tr_iter, 1);
    conv_rel_change   = nan(num_tr_iter, 1);
    num_iters_done    = 0;
    tr_time           = 0;

    switch recon_method
        case 'tr'
            % ---- Set up live TR figure (opt-in) ----
            fig_live   = [];
            ax_recon   = [];
            hImg_recon = [];
            hLine_max  = [];
            cz_live    = 1;
            if plot_results
                try
                    [~, dose_max_idx]               = max(doseGrid(:));
                    [~, ~, cz_live]                 = ind2sub([Nx_orig, Ny_orig, Nz_orig], dose_max_idx);
                    [fig_live, ax_recon, hImg_recon, hLine_max] = ...
                        setup_live_recon_figure(p0, Nx_orig, Ny_orig, cz_live, num_tr_iter, config);
                catch ME
                    warning('run_single_field_simulation:LiveFigFail', ...
                        'Live TR figure setup failed: %s', ME.message);
                    fig_live = [];
                end
            end

            fprintf('        Running iterative time reversal (%d iterations, tol=%.1e)...\n', ...
                num_tr_iter, convergence_tol);

            reconPressure_prev = zeros(gridSize);

            try
                tr_tic = tic;

                for tr_iter = 1:num_tr_iter

                    fprintf('          TR iteration %d/%d...\n', tr_iter, num_tr_iter);

                    source_tr        = struct();
                    source_tr.p_mask = sensor.mask;
                    source_tr.p      = fliplr(sensorData);
                    source_tr.p_mode = 'dirichlet';

                    sensor_tr        = struct();
                    sensor_tr.mask   = ones(Nx, Ny, Nz);
                    sensor_tr.record = {'p_final'};

                    p0_recon = kspaceFirstOrder3D(kgrid, kmedium_rec, source_tr, sensor_tr, inputArgs{:});

                    if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
                        reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
                    else
                        reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
                    end

                    reconPressure = max(reconPressure, 0);

                    conv_max_pressure(tr_iter) = max(reconPressure(:));
                    num_iters_done = tr_iter;

                    fprintf('          Max pressure: %.4e Pa\n', conv_max_pressure(tr_iter));

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

                    % Live figure update
                    if plot_results && ~isempty(fig_live) && ishandle(fig_live)
                        try
                            recon_slice_crop = squeeze( ...
                                reconPressure(1:Nx_orig, 1:Ny_orig, cz_live))';
                            recon_slice_crop = gather(recon_slice_crop);
                            set(hImg_recon, 'CData', recon_slice_crop);
                            caxis(ax_recon, [0, max(recon_slice_crop(:)) + eps]);
                            if converged
                                title(ax_recon, ...
                                    sprintf('Reconstructed p_0 (iter %d CONVERGED)', tr_iter), ...
                                    'FontWeight', 'bold', 'Color', [0, 0.55, 0]);
                            else
                                title(ax_recon, ...
                                    sprintf('Reconstructed p_0 (iter %d / %d)', ...
                                        tr_iter, num_tr_iter), 'FontWeight', 'bold');
                            end
                            set(hLine_max, 'XData', 1:tr_iter, ...
                                'YData', conv_max_pressure(1:tr_iter));
                            drawnow;
                        catch
                            % Ignore plot errors; don't break the recon
                        end
                    end

                    if converged
                        break;
                    end

                    % Residual correction for next iteration
                    if tr_iter < num_tr_iter
                        source_resid    = struct();
                        source_resid.p0 = reconPressure;
                        sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium_rec, ...
                            source_resid, sensor, inputArgs{:});
                        sensorData = sensorData + (sensorData_measured - sensorDataRecon);
                    end
                end

                reconPressure = gather(reconPressure) * correction_factor;
                tr_time = toc(tr_tic);
                fprintf('        Time reversal complete (%.1f s).\n', tr_time);
                fprintf('        Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
                    min(reconPressure(:)), max(reconPressure(:)));

            catch ME
                warning('run_single_field_simulation:TRFail', ...
                    'Time reversal failed: %s', ME.message);
                recon_dose = zeros(input_dims_for_crop);
                return;
            end

        case 'das'
            fprintf('        Running Delay-And-Sum reconstruction...\n');
            try
                das_tic = tic;
                das_opts = struct( ...
                    'use_elements', safe_config(config, 'das_use_elements', true), ...
                    'envelope',     safe_config(config, 'das_envelope',     true), ...
                    'use_aperture', safe_config(config, 'das_use_aperture', true), ...
                    'aperture_cos', safe_config(config, 'das_aperture_cos', 0.997), ...
                    'depth_weight', safe_config(config, 'das_depth_weight', true), ...
                    'interp',       safe_config(config, 'das_interp',       'nearest'));
                if blind_recon, das_medium = medium_recon; else, das_medium = medium; end
                reconPressure = das_reconstruct(sensorData, sensor, sensor_info, das_medium, ...
                                                Nx, Ny, Nz, dx, dy, dz, dt, das_opts);
                reconPressure = reconPressure * correction_factor;
                conv_max_pressure(1) = max(reconPressure(:));
                num_iters_done       = 1;
                tr_time              = toc(das_tic);
                fprintf('        DAS complete (%.1f s).\n', tr_time);
                fprintf('        Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
                    min(reconPressure(:)), max(reconPressure(:)));
            catch ME
                warning('run_single_field_simulation:DASFail', ...
                    'DAS reconstruction failed: %s', ME.message);
                recon_dose = zeros(input_dims_for_crop);
                return;
            end

        otherwise
            error('run_single_field_simulation:UnknownReconMethod', ...
                'Unknown reconstruction_method: "%s" (use ''tr'' or ''das'')', recon_method);
    end

    sim_results.tr_time_s         = tr_time;
    sim_results.recon_max         = max(reconPressure(:));
    sim_results.num_iters_done    = num_iters_done;
    sim_results.conv_max_pressure = conv_max_pressure(1:max(num_iters_done, 1));
    sim_results.conv_rel_change   = conv_rel_change(1:max(num_iters_done, 1));

    %% ======================== CROP TO ORIGINAL SIZE ========================
    %  Remove FFT padding before dose conversion (which uses the original-size
    %  gruneisen/density).

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

    %% ======================== PRESSURE SCALE CORRECTION ========================

    if use_pressure_scale_corr
        p0_max_orig = max(p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig), [], 'all');
        recon_max   = max(reconPressure(:));
        if recon_max > 0
            pressure_scale_cf = p0_max_orig / recon_max;
            reconPressure     = reconPressure * pressure_scale_cf;
            fprintf('        Pressure scale correction: %.4f (p0_max=%.3e / recon_max=%.3e)\n', ...
                pressure_scale_cf, p0_max_orig, recon_max);
            sim_results.pressure_scale_correction = pressure_scale_cf;
        else
            warning('run_single_field_simulation:RecoMaxZero', ...
                'Reconstructed pressure max is zero; skipping pressure scale correction.');
            sim_results.pressure_scale_correction = 1;
        end
    end

    %% ======================== PRESSURE -> DOSE CONVERSION ========================

    % Blind: the recovered pressure lives on the RECONSTRUCTION (reference)
    % geometry, so convert pressure->dose with the reconstruction medium's
    % gruneisen/density (medium_recon_orig). Full access: convert with the forward
    % medium (medium, restored to medium_orig by the crop-to-original step above),
    % exactly as before. The body mask still comes from cbct_resampled.
    if blind_recon
        conv_medium = medium_recon_orig;
    else
        conv_medium = medium;
    end

    conversionFactor = conv_medium.gruneisen .* conv_medium.density;
    conversionFactor(conversionFactor == 0) = 1;

    reconDosePerPulse = reconPressure ./ conversionFactor;

    body_mask_plot = ones(gridSize);
    if isfield(cbct_resampled, 'bodyMask') && isequal(size(cbct_resampled.bodyMask), gridSize)
        body_mask_plot = double(cbct_resampled.bodyMask);
    end
    if mask_recon_to_dose
        recon_dose = reconDosePerPulse * num_pulses .* double(doseMask) .* body_mask_plot;
    else
        recon_dose = reconDosePerPulse * num_pulses .* body_mask_plot;
    end

    % Mask air pockets: air-classified voxels (density ~1.2 kg/m^3, assigned
    % only inside the body for threshold_2) carry no real PA dose signal and
    % would otherwise reconstruct as hotspots. Use the reconstruction medium's
    % density (the geometry the recovered dose lives on) to find them.
    recon_dose(conv_medium.density < 100) = 0;

    fprintf('        Reconstructed dose: [%.4f, %.4f] Gy\n', ...
        min(recon_dose(:)), max(recon_dose(:)));

    %% ======================== DOSE NORMALIZATION ========================
    %  Optional: divide doseGrid and recon_dose by their respective maxima.
    %  Changes units from Gy to dimensionless.

    if normalize_doses
        dg_max = max(doseGrid(:));
        rd_max = max(recon_dose(:));
        fprintf('        [Norm] Original max: %.4f Gy, recon max: %.4f Gy -> normalising both to 1.\n', ...
            dg_max, rd_max);
        if dg_max > 0
            doseGrid = doseGrid / dg_max;
        end
        if rd_max > 0
            recon_dose = recon_dose / rd_max;
        end
        doseThreshold = 0.01 * max(doseGrid(:));
    end

    %% ======================== CROP RECON BACK TO INPUT SIZE ========================
    %  If grid expansion ran, the recon_dose currently includes water padding
    %  around the original dose region. Crop it out so the returned array
    %  matches the caller's expected field_dose.dose_Gy size (or the downscaled
    %  size, if downscale_factor != 1 was used).

    if ~isequal(size(recon_dose), input_dims_for_crop)
        eo  = expansion_offsets;
        idc = input_dims_for_crop;
        try
            recon_dose = recon_dose(eo(1) + (1:idc(1)), ...
                                    eo(2) + (1:idc(2)), ...
                                    eo(3) + (1:idc(3)));
            doseGrid   = doseGrid  (eo(1) + (1:idc(1)), ...
                                    eo(2) + (1:idc(2)), ...
                                    eo(3) + (1:idc(3)));
            fprintf('        Cropped recon_dose back to input size [%d %d %d]\n', idc);
        catch ME
            warning('run_single_field_simulation:CropFail', ...
                ['Could not crop recon_dose back to input size [%d %d %d] ' ...
                 '(current size [%d %d %d]): %s'], idc(1), idc(2), idc(3), ...
                size(recon_dose, 1), size(recon_dose, 2), size(recon_dose, 3), ...
                ME.message);
        end
    end

    %% ======================== DIAGNOSTIC ARRAYS (opt-in) ========================
    %  When CONFIG.return_diagnostics is true, hand back the pressure/geometry
    %  arrays the standalone/comparison analyses need, each windowed to the same
    %  input-dose region as recon_dose. Default-off so pipeline_simulate (parfor,
    %  memory-sensitive) is unaffected. reconPressure / sensor.mask / conv_medium
    %  are at gridSize_orig here; p0 is still at the padded working size, so take
    %  its leading original block before windowing.
    if return_diag
        eo  = expansion_offsets;
        idc = input_dims_for_crop;
        sim_results.recon_dose = recon_dose;     % input dims (convenience copy)
        sim_results.doseGrid   = doseGrid;       % input dims, masked (+ normalized if enabled)
        sim_results.spacing    = spacing_mm;
        try
            sim_results.reconPressure = window_to_input(reconPressure,       eo, idc);
            sim_results.sensor_mask   = logical(window_to_input(sensor.mask, eo, idc));
            sim_results.density       = window_to_input(conv_medium.density, eo, idc);
            p0_orig                   = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
            sim_results.p0            = window_to_input(p0_orig,             eo, idc);
        catch ME
            warning('run_single_field_simulation:DiagCropFail', ...
                'Could not window diagnostic arrays to input size: %s', ME.message);
        end
    end

    %% ======================== RESULTS SUMMARY ========================

    fprintf('\n        ========= RESULTS =========\n');
    fprintf('          Original dose:      [%.6f, %.4f]\n', min(doseGrid(:)), max(doseGrid(:)));
    fprintf('          Reconstructed dose: [%.6f, %.4f]\n', min(recon_dose(:)), max(recon_dose(:)));
    dose_region = doseGrid > doseThreshold;
    if any(dose_region(:))
        abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
        rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
        fprintf('          Mean abs err: %.6f\n', mean(abs_error));
        fprintf('          Mean rel err: %.2f%%\n', mean(rel_error));
        fprintf('          Max  rel err: %.2f%%\n', max(rel_error));
    end
    fprintf('        ===========================\n\n');

    %% ======================== POST-SIM PLOTTING (opt-in) ========================

    if plot_results
        % Crop sensor.mask and density to the same window as doseGrid /
        % recon_dose so the panels line up after grid expansion / FFT padding.
        eo  = expansion_offsets;
        idc = input_dims_for_crop;
        try
            sensor_for_plot = sensor.mask(eo(1) + (1:idc(1)), ...
                                          eo(2) + (1:idc(2)), ...
                                          eo(3) + (1:idc(3)));
            % Density background = geometry the recovered dose lives on:
            % reconstruction medium when blind, forward medium for full access.
            if blind_recon
                plot_density_full = medium_recon_orig.density;
                panel_title = 'Dose Comparison: True Dose vs Blind Reconstruction';
            else
                plot_density_full = medium_orig.density;
                panel_title = 'Dose Comparison: Original vs Reconstructed';
            end
            density_for_plot = plot_density_full(eo(1) + (1:idc(1)), ...
                                                 eo(2) + (1:idc(2)), ...
                                                 eo(3) + (1:idc(3)));
            plot_dose_panels(doseGrid, recon_dose, sensor_for_plot, ...
                density_for_plot, spacing_mm, panel_title);
        catch ME
            warning('run_single_field_simulation:DosePanelsFail', ...
                'plot_dose_panels failed: %s', ME.message);
        end
        if strcmpi(recon_method, 'tr')
            try
                p0_max_for_plot = max(p0(:));
                plot_convergence_history(conv_max_pressure, conv_rel_change, ...
                    num_iters_done, convergence_tol, p0_max_for_plot);
            catch ME
                warning('run_single_field_simulation:ConvPlotFail', ...
                    'plot_convergence_history failed: %s', ME.message);
            end
        end
    end

    fprintf('        Field simulation complete.\n');
end


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function val = safe_config(config, field_name, default_val)
%SAFE_CONFIG Retrieve config field with fallback to default.
    if isfield(config, field_name) && ~isempty(config.(field_name))
        val = config.(field_name);
    else
        val = default_val;
    end
end


function B = window_to_input(A, eo, idc)
%WINDOW_TO_INPUT Slice the input-dose window (offset eo, size idc) out of a
%  gridSize_orig-sized array so a diagnostic array aligns with recon_dose. With
%  no grid expansion eo = [0 0 0] and idc == size(A), so this is the identity.
    B = A(eo(1) + (1:idc(1)), eo(2) + (1:idc(2)), eo(3) + (1:idc(3)));
end


function tf = precomputed_sensor_fits(precomputed_sensor, base_dims)
%PRECOMPUTED_SENSOR_FITS True when a caller-supplied sensor mask can be reused.
%  Requires a struct with .sensor_mask and .sensor_info whose un-expanded base
%  grid (mask size minus any determine_sensor_mask grid_pad) equals this
%  field's grid base_dims = [Nx_orig, Ny_orig, Nz_orig]. This guards against
%  a mismatch (e.g. a non-unity downscale_factor shrank this field's grid),
%  in which case the caller's mask is ignored and the mask is computed inline.
    tf = false;
    if isempty(precomputed_sensor) || ~isstruct(precomputed_sensor) || ...
            ~isfield(precomputed_sensor, 'sensor_mask') || ...
            ~isfield(precomputed_sensor, 'sensor_info')
        return;
    end
    sm = precomputed_sensor.sensor_mask;
    info = precomputed_sensor.sensor_info;
    if isstruct(info) && isfield(info, 'grid_pad') && isstruct(info.grid_pad)
        gp = info.grid_pad;
        b1 = size(sm, 1) - (gp.y_pre + gp.y_post);   % grid_pad y_* -> dim 1
        b2 = size(sm, 2) - (gp.x_pre + gp.x_post);   % grid_pad x_* -> dim 2
        b3 = size(sm, 3) - (gp.z_pre + gp.z_post);   % grid_pad z_* -> dim 3
    else
        b1 = size(sm, 1); b2 = size(sm, 2); b3 = size(sm, 3);
    end
    tf = isequal([b1, b2, b3], base_dims(:)');
end


function [medium, p0] = pad_medium_and_p0(medium, p0, Nx_src, Ny_src, Nz_src, Nx_dst, Ny_dst, Nz_dst)
%PAD_MEDIUM_AND_P0 Embed an Nx_src x Ny_src x Nz_src medium + p0 into a
%  larger Nx_dst x Ny_dst x Nz_dst grid, filling the padding region with
%  water properties (density 1000, sound_speed 1540, alpha 0, gruneisen 0)
%  and zero pressure.

    density_pad    = ones(Nx_dst, Ny_dst, Nz_dst) * 1000;
    soundSpeed_pad = ones(Nx_dst, Ny_dst, Nz_dst) * 1540;
    alphaCoeff_pad = zeros(Nx_dst, Ny_dst, Nz_dst);
    gruneisen_pad  = zeros(Nx_dst, Ny_dst, Nz_dst);
    p0_pad         = zeros(Nx_dst, Ny_dst, Nz_dst);

    density_pad(1:Nx_src, 1:Ny_src, 1:Nz_src)    = medium.density(1:Nx_src, 1:Ny_src, 1:Nz_src);
    soundSpeed_pad(1:Nx_src, 1:Ny_src, 1:Nz_src) = medium.sound_speed(1:Nx_src, 1:Ny_src, 1:Nz_src);
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_pad(1:Nx_src, 1:Ny_src, 1:Nz_src) = ...
            medium.alpha_coeff(1:Nx_src, 1:Ny_src, 1:Nz_src);
    else
        alphaCoeff_pad(:) = medium.alpha_coeff;
    end
    gruneisen_pad(1:Nx_src, 1:Ny_src, 1:Nz_src) = ...
        medium.gruneisen(1:Nx_src, 1:Ny_src, 1:Nz_src);
    p0_pad(1:Nx_src, 1:Ny_src, 1:Nz_src) = p0(1:Nx_src, 1:Ny_src, 1:Nz_src);

    medium.density     = density_pad;
    medium.sound_speed = soundSpeed_pad;
    medium.alpha_coeff = alphaCoeff_pad;
    medium.gruneisen   = gruneisen_pad;
    p0 = p0_pad;
end


function m = pad_medium_to(m, sz)
%PAD_MEDIUM_TO Water-pad the four medium fields into an sz grid.
%  Copies each existing field into the leading sub-block of a fresh sz array
%  filled with the background (water: density 1000, c 1540, alpha/gruneisen 0),
%  mirroring the forward-medium FFT/grid padding so the reconstruction medium
%  stays voxel-aligned with it. Scalar .alpha_power / .grid_size are preserved.
    d = ones(sz) * 1000;
    d(1:size(m.density,1),     1:size(m.density,2),     1:size(m.density,3))     = m.density;
    c = ones(sz) * 1540;
    c(1:size(m.sound_speed,1), 1:size(m.sound_speed,2), 1:size(m.sound_speed,3)) = m.sound_speed;
    a = zeros(sz);
    if numel(m.alpha_coeff) > 1
        a(1:size(m.alpha_coeff,1), 1:size(m.alpha_coeff,2), 1:size(m.alpha_coeff,3)) = m.alpha_coeff;
    else
        a(:) = m.alpha_coeff;
    end
    g = zeros(sz);
    g(1:size(m.gruneisen,1),   1:size(m.gruneisen,2),   1:size(m.gruneisen,3))   = m.gruneisen;
    m.density = d; m.sound_speed = c; m.alpha_coeff = a; m.gruneisen = g;
end


function m = expand_medium(m_unp, sz, xr, yr, zr)
%EXPAND_MEDIUM Place cropped medium fields into a water-filled sz grid at the
%  index ranges (xr,yr,zr), matching the sensor-driven grid expansion applied
%  to the forward medium.
    d = ones(sz) * 1000;  d(xr, yr, zr) = m_unp.density;
    c = ones(sz) * 1540;  c(xr, yr, zr) = m_unp.sound_speed;
    a = zeros(sz);
    if numel(m_unp.alpha_coeff) > 1
        a(xr, yr, zr) = m_unp.alpha_coeff;
    else
        a(:) = m_unp.alpha_coeff;
    end
    g = zeros(sz);        g(xr, yr, zr) = m_unp.gruneisen;
    m = struct('density', d, 'sound_speed', c, 'alpha_coeff', a, 'gruneisen', g);
end


function m = crop_medium(m, sz)
%CROP_MEDIUM Crop the four medium fields to the leading sz sub-block.
    m.density     = m.density(1:sz(1),     1:sz(2),     1:sz(3));
    m.sound_speed = m.sound_speed(1:sz(1), 1:sz(2),     1:sz(3));
    if numel(m.alpha_coeff) > 1
        m.alpha_coeff = m.alpha_coeff(1:sz(1), 1:sz(2), 1:sz(3));
    end
    m.gruneisen   = m.gruneisen(1:sz(1),   1:sz(2),     1:sz(3));
end


function [fig_live, ax_recon, hImg_recon, hLine_max] = ...
        setup_live_recon_figure(p0, Nx_orig, Ny_orig, cz_live, N_iter, config)
%SETUP_LIVE_RECON_FIGURE Build the 1x3 live reconstruction figure for the
%  'tr' branch. Panels: (1) initial p0 slice, (2) current recon, (3) live
%  max-pressure convergence.

    fig_live = figure('Name', 'Live Reconstruction', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1060, 440]);

    ax_p0 = subplot(1, 3, 1);
    p0_orig_slice = squeeze(p0(1:Nx_orig, 1:Ny_orig, cz_live))';
    imagesc(ax_p0, p0_orig_slice);
    axis(ax_p0, 'xy'); axis(ax_p0, 'image');
    colormap(ax_p0, 'hot'); colorbar(ax_p0);
    clim_p0 = [0, max(p0_orig_slice(:)) + eps];
    caxis(ax_p0, clim_p0);
    xlabel(ax_p0, 'X (voxel)'); ylabel(ax_p0, 'Y (voxel)');
    title(ax_p0, sprintf('Initial p_0 (Z=%d)', cz_live), 'FontWeight', 'bold');

    ax_recon = subplot(1, 3, 2);
    hImg_recon = imagesc(ax_recon, zeros(Ny_orig, Nx_orig));
    axis(ax_recon, 'xy'); axis(ax_recon, 'image');
    colormap(ax_recon, 'hot'); colorbar(ax_recon);
    xlabel(ax_recon, 'X (voxel)'); ylabel(ax_recon, 'Y (voxel)');
    title(ax_recon, 'Reconstructed p_0 (iter 0)', 'FontWeight', 'bold');

    ax_conv = subplot(1, 3, 3);
    hLine_max = plot(ax_conv, NaN, NaN, 'b-o', 'LineWidth', 1.6, ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.2, 0.4, 1.0]);
    xlabel(ax_conv, 'Iteration');
    ylabel(ax_conv, 'Max Reconstructed p_0 (Pa)');
    title(ax_conv, 'Convergence (live)', 'FontWeight', 'bold');
    grid(ax_conv, 'on');
    xlim(ax_conv, [0.5, N_iter + 0.5]);

    recon_method_label = safe_config(config, 'reconstruction_method', 'tr');
    sgtitle(fig_live, sprintf('Live Reconstruction (%s)   Axial Z=%d', ...
        recon_method_label, cz_live), 'FontWeight', 'bold', 'FontSize', 11);
    drawnow;
end


function plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, density, config) %#ok<INUSD>
%PLOT_SENSOR_DOSE_PLANES Three anatomical views of dose mask + sensor.
%  Lightweight 1x3 figure: axial / sagittal / coronal slices through the
%  dose centroid. Sensor mask drawn as red outline; dose mask as semi-
%  transparent yellow fill; density background in grayscale.

    if ~any(dose_mask(:))
        return;
    end

    sz = size(dose_mask);
    [Xc, Yc, Zc] = ind2sub(sz, find(dose_mask));
    cx = round(mean(Xc));
    cy = round(mean(Yc));
    cz = round(mean(Zc));

    fig = figure('Name', 'Sensor vs Dose Mask', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [120, 120, 1140, 360]);

    views = { ...
        struct('name', 'Axial (Z = const)',    'slice_dim', 3, 'idx', cz), ...
        struct('name', 'Coronal (Y = const)',  'slice_dim', 2, 'idx', cy), ...
        struct('name', 'Sagittal (X = const)', 'slice_dim', 1, 'idx', cx)};

    for v = 1:3
        ax = subplot(1, 3, v);
        switch views{v}.slice_dim
            case 1
                d_slice = squeeze(double(density(views{v}.idx, :, :)))';
                m_slice = squeeze(dose_mask(views{v}.idx, :, :))';
                s_slice = squeeze(sensor_mask(views{v}.idx, :, :))';
            case 2
                d_slice = squeeze(double(density(:, views{v}.idx, :)))';
                m_slice = squeeze(dose_mask(:, views{v}.idx, :))';
                s_slice = squeeze(sensor_mask(:, views{v}.idx, :))';
            case 3
                d_slice = squeeze(double(density(:, :, views{v}.idx)))';
                m_slice = squeeze(dose_mask(:, :, views{v}.idx))';
                s_slice = squeeze(sensor_mask(:, :, views{v}.idx))';
        end

        imagesc(ax, d_slice);
        axis(ax, 'xy'); axis(ax, 'image');
        colormap(ax, 'gray');
        hold(ax, 'on');

        if any(m_slice(:))
            [Bd, ~] = bwboundaries(m_slice);
            for k = 1:numel(Bd)
                plot(ax, Bd{k}(:, 2), Bd{k}(:, 1), 'y-', 'LineWidth', 1.2);
            end
        end
        if any(s_slice(:))
            [Bs, ~] = bwboundaries(s_slice);
            for k = 1:numel(Bs)
                plot(ax, Bs{k}(:, 2), Bs{k}(:, 1), 'r-', 'LineWidth', 1.6);
            end
        end

        title(ax, views{v}.name, 'FontWeight', 'bold');
        hold(ax, 'off');
    end

    sgtitle(fig, sprintf('Sensor (red) vs >10%% Dose (yellow)   |   spacing %.2f mm', ...
        spacing_mm(1)), 'FontWeight', 'bold');
    drawnow;
end


function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr) %#ok<INUSD>
%PLOT_DOSE_PANELS 2x3 panel of original (top) vs recon (bottom).
%  Columns: coronal / sagittal / axial through the max-dose voxel.
%  Density shown as grayscale background; dose overlaid hot+alpha;
%  sensor mask as red contour.

    if ~any(original(:) > 0)
        return;
    end

    sz = size(original);
    [~, idx] = max(original(:));
    [cx, cy, cz] = ind2sub(sz, idx);

    fig = figure('Name', 'Dose Comparison', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 1280, 720]);

    dose_max = max([max(original(:)), max(recon(:)), eps]);

    cols = { ...
        struct('name', 'Coronal (Y)',  'slice_dim', 2, 'idx', cy), ...
        struct('name', 'Sagittal (X)', 'slice_dim', 1, 'idx', cx), ...
        struct('name', 'Axial (Z)',    'slice_dim', 3, 'idx', cz)};

    rows = { ...
        struct('label', 'Original dose',      'data', original), ...
        struct('label', 'Reconstructed dose', 'data', recon)};

    for r = 1:2
        for c = 1:3
            ax  = subplot(2, 3, (r - 1) * 3 + c);
            sd  = cols{c}.slice_dim;
            idx_s = cols{c}.idx;

            switch sd
                case 1
                    d_slice = squeeze(double(density(idx_s, :, :)))';
                    z_slice = squeeze(rows{r}.data(idx_s, :, :))';
                    s_slice = squeeze(sensor_mask(idx_s, :, :))';
                case 2
                    d_slice = squeeze(double(density(:, idx_s, :)))';
                    z_slice = squeeze(rows{r}.data(:, idx_s, :))';
                    s_slice = squeeze(sensor_mask(:, idx_s, :))';
                case 3
                    d_slice = squeeze(double(density(:, :, idx_s)))';
                    z_slice = squeeze(rows{r}.data(:, :, idx_s))';
                    s_slice = squeeze(sensor_mask(:, :, idx_s))';
            end

            imagesc(ax, d_slice);
            colormap(ax, 'gray');
            axis(ax, 'xy'); axis(ax, 'image');
            hold(ax, 'on');

            hot_layer = z_slice;
            hImg = imagesc(ax, hot_layer);
            colormap(ax, 'hot');
            set(hImg, 'AlphaData', 0.45 * (hot_layer > 0));
            caxis(ax, [0, dose_max]);

            if any(s_slice(:))
                [Bs, ~] = bwboundaries(s_slice);
                for k = 1:numel(Bs)
                    plot(ax, Bs{k}(:, 2), Bs{k}(:, 1), 'r-', 'LineWidth', 1.4);
                end
            end

            if r == 1
                title(ax, cols{c}.name, 'FontWeight', 'bold');
            end
            if c == 1
                ylabel(ax, rows{r}.label, 'FontWeight', 'bold');
            end
            hold(ax, 'off');
        end
    end

    sgtitle(fig, titleStr, 'FontWeight', 'bold');
    drawnow;
end


function plot_convergence_history(conv_max_pressure, conv_rel_change, ...
                                  num_iters, tol, p0_max_orig)
%PLOT_CONVERGENCE_HISTORY 1x2 plot of max pressure and relative change.
%  Left: max recon pressure per iter vs initial p0 reference (dashed).
%  Right: relative change per iter vs convergence tolerance (dashed).

    iters = 1:num_iters;
    cm    = conv_max_pressure(1:num_iters);
    rc    = conv_rel_change(1:num_iters);

    fig = figure('Name', 'TR Convergence', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [200, 200, 980, 420]);

    ax1 = subplot(1, 2, 1);
    plot(ax1, iters, cm, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', [0.2, 0.4, 1]);
    hold(ax1, 'on');
    if ~isempty(p0_max_orig) && p0_max_orig > 0
        yline(ax1, p0_max_orig, 'k--', 'LineWidth', 1.2, ...
            'Label', sprintf('max p_0 = %.2e', p0_max_orig));
    end
    xlabel(ax1, 'Iteration');
    ylabel(ax1, 'Max p (Pa)');
    title(ax1, 'Max reconstructed pressure', 'FontWeight', 'bold');
    grid(ax1, 'on');
    hold(ax1, 'off');

    ax2 = subplot(1, 2, 2);
    semilogy(ax2, iters, max(rc, eps), 'r-s', 'LineWidth', 1.5, ...
        'MarkerFaceColor', [1, 0.4, 0.4]);
    hold(ax2, 'on');
    if ~isempty(tol) && tol > 0
        yline(ax2, tol, 'k--', 'LineWidth', 1.2, ...
            'Label', sprintf('tol = %.1e', tol));
    end
    xlabel(ax2, 'Iteration');
    ylabel(ax2, 'Relative change');
    title(ax2, 'Iteration relative change', 'FontWeight', 'bold');
    grid(ax2, 'on');
    hold(ax2, 'off');

    sgtitle(fig, 'Time-Reversal Convergence', 'FontWeight', 'bold');
    drawnow;
end
