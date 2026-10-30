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
%                                                  'determine_sensor_mask_lateral',
%                                                  'fixed_anterior'.
%           .sensor_x_index             [1]      - YZ plane x-index
%           .sensor_y_index             [1]      - XZ plane y-index
%           .pulse_shape                ['gaussian'] - 'gaussian'|'rectangular' pulse profile
%           .convolution_kernel         [4e-6]   - Gaussian sigma / rectangular width (s); 0 disables
%           .conv_noise_level           [0.01]   - Legacy noise fraction of THIS
%                                                  field's peak (used only when
%                                                  noise_amp_Pa is unset)
%           .noise_amp_Pa               []       - Fixed absolute electronic-noise
%                                                  amplitude (Pa). When >0 it is
%                                                  used for every field (SNR then
%                                                  floats field to field); when
%                                                  empty, conv_noise_level is used.
%                                                  Set it via calibrate_noise_amp.
%           .forward_peak_only          [false]  - Calibration fast-path: return
%                                                  after computing the pre-noise
%                                                  sensor peak (sim_results
%                                                  .sensor_signal_peak); skips
%                                                  noise, deconv, and reconstruction.
%           .sensor_stages_only         [false]  - Diagnostic fast-path (used by
%                                                  verify_pipeline_simulate): record
%                                                  the sensor-averaged |p(t)| after
%                                                  each processing stage in
%                                                  sim_results.sensor_stages and
%                                                  return before reconstruction.
%           .conv_deconv_lambda         [1e-4]   - Wiener deconvolution regularisation
%           .conv_variable_lambda       [false]  - Deconvolve pulse x transducer with
%                                                  lambda = noise power / signal
%                                                  power of THIS field, instead of
%                                                  the fixed conv_deconv_lambda on
%                                                  the pulse only. Assumes the signal
%                                                  spectrum has the sensor response's
%                                                  shape (warns). See
%                                                  deconv_pulse_transducer.
%           .conv_use_deconv            [true]   - false skips the deconvolution: the
%                                                  data keep the pulse smear and the
%                                                  pulse kernel's time delay
%                                                  (build_pulse_kernel is not
%                                                  zero-lag centred).
%           .noise_only                 [false]  - Null-signal control: run the
%                                                  forward sim to fix the noise
%                                                  amplitude, then discard the true
%                                                  signal so ONLY electronic noise
%                                                  is reconstructed. Whatever dose
%                                                  comes out is pure noise artifact,
%                                                  a baseline for "is the recon
%                                                  actually recovering signal?".
%           .reconstruction_method      ['tr']   - 'tr', 'das', 'ubp', 'adjoint', or
%                                                  'adjoint_tv'
%                                                  ('ubp' = universal back-
%                                                  projection via ubp3d; honors
%                                                  ubp_use_elements, ubp_band_limit,
%                                                  ubp_aperture_cos, ubp_area_weight,
%                                                  ubp_interp. 'adjoint' = Arridge
%                                                  et al. 2016 gradient descent
%                                                  with the k-Wave adjoint.
%                                                  'adjoint_tv' = the same operators
%                                                  with a TV penalty, solved by
%                                                  FISTA, Huang et al. 2013)
%           .adjoint_num_iter           [5]      - 'adjoint' / 'adjoint_tv'
%                                                  iterations ('adjoint': up to 3
%                                                  k-Wave runs each; 'adjoint_tv':
%                                                  2 each); stops early on
%                                                  convergence_tol
%           .adjoint_tv_lambda          [0.01]   - 'adjoint_tv' TV weight,
%                                                  dimensionless: the TV-denoise
%                                                  strength is this fraction of the
%                                                  peak of the first unregularised
%                                                  step (0 = no TV, positivity only)
%           .adjoint_tv_force_power_recalc [false] - 'adjoint_tv': true ignores the
%                                                  cached power iteration
%                                                  (processed/simulated_doses/
%                                                  adjoint_power_cache_<hash>.mat),
%                                                  recomputes it and
%                                                  overwrites the cache
%           .adjoint_tv_power_iter      [5]      - 'adjoint_tv' power iterations for
%                                                  the step size (2 k-Wave runs each,
%                                                  skipped when cached)
%           .adjoint_tv_inner_iter      [10]     - 'adjoint_tv' TV-denoise iterations
%                                                  per outer iteration (no k-Wave)
%           .blind_recon                [false]  - When true (and medium_recon is
%                                                  supplied), forward on `medium`
%                                                  but reconstruct on medium_recon.
%           .use_attenuation            [true]   - Use medium.alpha_coeff in kspaceFirstOrder3DG
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
%           .plot_results               [false]  - Opt-in: pre-sim sensor view, live TR/
%                                                  adjoint figure, post-sim dose panels +
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
%       recon_dose  - 3D reconstructed dose array (Gy), always returned at the
%                     input field_dose.dose_Gy size. With a non-unity
%                     downscale_factor the recon is resampled back up from the
%                     coarse simulation grid to these native dimensions.
%       sim_results - Struct with diagnostics:
%           .sensor_info        - Sensor placement info
%           .forward_time_s     - Forward sim wall time
%           .tr_time_s          - Reconstruction wall time (TR or DAS)
%           .num_pulses         - Number of LINAC pulses
%           .p0_max             - Maximum initial pressure (Pa)
%           .sensor_signal_peak - Pre-noise peak sensor pressure (Pa); the
%                                 signal the SNR/noise amplitude is measured
%                                 against (pulse-model path: after the pulse
%                                 convolution, before the band-limit)
%           .noise_amp          - Electronic-noise amplitude applied (Pa)
%           .snr                - This field's SNR (sensor_signal_peak/noise_amp)
%           .recon_max          - Maximum reconstructed pressure (Pa)
%           .num_iters_done     - Iterations completed
%           .conv_max_pressure  - Per-iteration max pressure
%           .conv_rel_change    - Per-iteration relative change
%           .grid_padding       - Original / padded / expanded sizes
%           .recon_method       - Which branch ran ('tr' | 'das' | 'ubp' | 'adjoint' |
%                                 'adjoint_tv')
%           .conv_cost          - 'adjoint_tv' only: objective ||A p - d||^2 +
%                                 lambda*TV(p) per iteration
%           .blind_recon        - Logical: true when forward and reconstruction
%                                 media differed (blind geometry)
%           .sensor_stages      - Only with CONFIG.sensor_stages_only: .names
%                                 (stage labels), .mean_abs ([stage x Nt] mean
%                                 |p| over all sensor points, Pa), .dt (s)
%
%   NOTES:
%       - Stateless: safe for parfor execution. Diagnostic plotting must be
%         disabled under parfor (figures from workers do not render).
%       - PlotSim disabled inside kspaceFirstOrder3DG in all cases.
%       - Returns zeros if no significant dose or on simulation failure.
%       - beam_metadata is the FULL plan metadata (all beams), passed through
%         to determine_sensor_mask for computing the exclusion zone.
%       - For backwards compatibility, beam_metadata can be omitted (pass []),
%         in which case the legacy sensor placement is used.
%       - precomputed_sensor lets the caller compute the determine_sensor_mask
%         placement once for a whole plan and reuse it for every field. When it
%         is absent or its grid does not match, the mask is computed inline.
%       - Element averaging (apply_element_averaging) runs on the raw forward
%         data, before the pulse/noise/deconvolution chain; electronic noise
%         is then one independent trace per element (element_noise). Every
%         reconstruction method gets this element-averaged data; the TR
%         residual re-forward and the adjoint forward model are averaged the
%         same way. No-op for sensors without an element map.
%
%   DATE: May 2026
%   VERSION: 3.0 (DAS, determine_sensor_mask + grid expansion, scale correction,
%                 dose normalisation, downscale_factor, attenuation toggle,
%                 force-uniform medium overrides, opt-in plotting)
%
%   See also: determine_sensor_mask, determine_sensor_placement_fixed,
%             kspaceFirstOrder3DG, run_standalone_simulation

    %% ======================== CONFIG DEFAULTS ========================

    dose_per_pulse_cGy  = safe_config(config, 'dose_per_pulse_cGy', 0.16);
    pml_size            = safe_config(config, 'pml_size', 10);
    cfl                 = safe_config(config, 'cfl_number', 0.3);
    use_gpu             = safe_config(config, 'use_gpu', true);
    num_tr_iter         = safe_config(config, 'num_time_reversal_iter', 30);
    num_adj_iter        = safe_config(config, 'adjoint_num_iter', 5);
    tv_lambda           = safe_config(config, 'adjoint_tv_lambda', 0.01);
    tv_force_power_recalc = safe_config(config, 'adjoint_tv_force_power_recalc', false);
    tv_power_iter       = safe_config(config, 'adjoint_tv_power_iter', 5);
    tv_inner_iter       = safe_config(config, 'adjoint_tv_inner_iter', 10);
    convergence_tol     = safe_config(config, 'convergence_tol', 1e-3);
    correction_factor   = safe_config(config, 'correction_factor', 1.9);
    use_grid_padding    = safe_config(config, 'use_grid_padding', true);
    conv_kernel_sigma   = safe_config(config, 'convolution_kernel', 4e-6);
    pulse_shape         = safe_config(config, 'pulse_shape', 'gaussian');
    conv_noise_level    = safe_config(config, 'conv_noise_level', 0.01);
    conv_deconv_lambda  = safe_config(config, 'conv_deconv_lambda', 1e-4);
    conv_variable_lambda = safe_config(config, 'conv_variable_lambda', false);
    conv_use_deconv     = safe_config(config, 'conv_use_deconv', true);
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
    %  Debug sensitivity-test toggles: override one or more medium properties
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
    %  Debug pre-sim downscale of the whole grid (dose, medium, body/couch
    %  masks). Useful for rapid testing / large speedups (coarser dx also grows
    %  the CFL time step, so cost falls faster than voxel count alone). The
    %  reconstructed dose is resampled back to the native dimensions at the end
    %  of the function, so the caller always receives recon_dose at the input
    %  field_dose.dose_Gy size.

    % Native dose dims/spacing (pre-downscale) for the resample-back at the end.
    native_dims       = gridSize;
    native_spacing_mm = spacing_mm;

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
        recon_dose = zeros(native_dims);
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

        case {'determine_sensor_mask', 'determine_sensor_mask_lateral'}
            % Automatic placement: tilts a 2D array toward the dose centroid
            % avoiding beam exclusion zones. 'determine_sensor_mask' presses the
            % array against the ANTERIOR abdomen; 'determine_sensor_mask_lateral'
            % presses it against the RIGHT/LEFT flank (config.sensor_side). Both
            % return the same [sensor_mask, sensor_info] contract (including
            % sensor_info.grid_pad) and may expand the simulation grid in X/Y/Z
            % with water padding to place the sensor outside the exclusion zone.
            lateral_sensor = strcmp(sensor_method, 'determine_sensor_mask_lateral');
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
            elseif lateral_sensor
                warning("Precomputed sensor did not fit. Recalculating.")
                [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask_lateral( ...
                    sct_for_sensor, field_dose_for_sensor, beam_metadata, config);
            else
                warning("Precomputed sensor did not fit. Recalculating.")

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
            fprintf('        Sensor: %s - %d active points\n', ...
                sensor_method, sum(sensor_mask_orig(:)));

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
                 '''box'', ''determine_sensor_mask'', ''determine_sensor_mask_lateral'', ' ...
                 'or ''fixed_anterior''.'], ...
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

    gpuIndex = [];   % MATLAB's GPU for this process (1-based); [] = none
    if use_gpu
        try
            gpu = gpuDevice;
            gpuIndex = gpu.Index;
            dataCast = 'gpuArray-single';
            fprintf('        Compute: GPU %d\n', gpuIndex);
        catch
            dataCast = 'single';
            fprintf('        Compute: CPU (GPU unavailable)\n');
        end
    else
        dataCast = 'single';
        fprintf('        Compute: CPU\n');
    end

    % kspaceFirstOrder3DG runs a separate CUDA binary that reads and writes
    % .h5 files in tempdir. By default the files are named by the time to the
    % second, so parfor workers starting in the same second overwrite each
    % other's files (forward sim fails, or simulates another field's dose).
    % A unique DataName per field avoids that. The binary also ignores MATLAB's
    % GPU selection, so DeviceNum passes this worker's GPU (binary counts from 0).
    [~, kwaveDataName] = fileparts(tempname);
    inputArgs = {'Smooth', false, ...
                 'PMLInside', false, ...
                 'PMLSize', pml_size, ...
                 'DataCast', dataCast, ...
                 'PlotSim', false, ...
                 'DataName', kwaveDataName};
    if ~isempty(gpuIndex)
        inputArgs = [inputArgs, {'DeviceNum', gpuIndex - 1}];
    end

    %% ======================== FORWARD SIMULATION ========================

    fprintf('        Running forward simulation...\n');

    source_fwd    = struct();
    source_fwd.p0 = p0;

    try
        fwd_tic = tic;
        setenv("HDF5_DISABLE_VERSION_CHECK","1")
        sensorData = kspaceFirstOrder3DG(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
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

    %% ============ DIAGNOSTIC EXPORT: RAW FORWARD SINOGRAM PAIR (opt-in) ============
    %  For diagnostic_blind_fidelity. When config.export_forward_sinograms is set,
    %  hand back the RAW (pre-pulse, pre-noise) sensor data for the SAME p0 and
    %  SAME sensor propagated through BOTH media: the forward medium (kmedium) and,
    %  on the blind path, the reconstruction medium (kmedium_rec). This is the
    %  clean CT1-vs-CT3 channel-by-channel comparison (medium isolated: p0 and
    %  sensor held fixed). Guarded and inert for every normal caller.
    if safe_config(config, 'export_forward_sinograms', false)
        sim_results.sinogram_fwd = gather(double(sensorData));
        if blind_recon
            sd_alt = kspaceFirstOrder3DG(kgrid, kmedium_rec, source_fwd, sensor, inputArgs{:});
            sim_results.sinogram_alt = gather(double(sd_alt));
        else
            % Full access: no second medium, so the pair is degenerate.
            sim_results.sinogram_alt = sim_results.sinogram_fwd;
        end
        sim_results.sinogram_dt   = kgrid.dt;
        sim_results.sensor_info   = sensor_info;   % element_positions_mm, voxel_element_idx
        % Skip the recon entirely when only the sinogram pair is wanted (fast,
        % and avoids the diverging blind time reversal the diagnostic is probing).
        if safe_config(config, 'forward_sinograms_only', false)
            recon_dose = zeros(input_dims_for_crop);
            return;
        end
    end

    %% ======================== ELEMENT AVERAGING ========================
    %  The physical array reports one signal per element, not one per voxel:
    %  each element integrates the pressure over its face. Average the forward
    %  data over each element's voxels (every voxel of an element then carries
    %  its element's signal), so every later stage -- including the noise and
    %  every reconstruction method -- works on what the real array measures.
    %  The pulse and band-limit filters act on each channel separately, so
    %  averaging before them gives the same result as averaging after.
    %  apply_element_averaging passes the data through unchanged when the
    %  sensor has no element map (plane, sphere, box sensors). DAS/UBP average
    %  again internally, which changes nothing.
    [~, sensorData]     = apply_element_averaging(double(gather(sensorData)), sensor_info);
    sensorData          = single(sensorData);
    sensorData_measured = sensorData;

    %% ================= PULSE CONVOLUTION / NOISE / FREQUENCY RESPONSE / DECONV =================
    %  Physical measurement chain in acquisition order (on element-averaged data):
    %    1. Convolve with the pulse profile (finite radiation pulse)
    %    2. Add broadband electronic noise, one independent trace per element
    %       channel (element_noise), BEFORE the band-limit so the transducer
    %       filters noise and signal alike
    %    3. Apply the sensor frequency response (band-limit)
    %    4. Deconvolve (deconv_pulse_transducer): Wiener on the pulse kernel, or
    %       with conv_variable_lambda the noise-weighted Wiener on pulse x transducer

    if conv_kernel_sigma > 0
        fprintf('        Pulse model: %s, width=%.1f us, lambda=%.1e\n', ...
            pulse_shape, conv_kernel_sigma * 1e6, conv_deconv_lambda);

        pulse_kernel = build_pulse_kernel(pulse_shape, conv_kernel_sigma, dt);

        sensorData_cpu = double(gather(sensorData));
        Nt_data        = size(sensorData_cpu, 2);

        H = fft(pulse_kernel, Nt_data).';

        % 1. Convolve with the pulse profile
        sensorData_conv = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));

        % Pre-noise peak of THIS field's sensor signal (Pa), taken where the noise
        % is added (after the pulse convolution, before the band-limit). This is
        % the quantity the noise amplitude is measured against, and what
        % calibrate_noise_amp samples to fix a single blanket noise amplitude for
        % a target mean SNR.
        sensor_signal_peak             = max(abs(sensorData_conv(:)));
        sim_results.sensor_signal_peak = sensor_signal_peak;

        % Calibration fast-path: calibrate_noise_amp only needs the pre-noise
        % signal peak, so skip noise, deconvolution, and the expensive
        % time-reversal reconstruction entirely.
        if safe_config(config, 'forward_peak_only', false)
            recon_dose = zeros(input_dims_for_crop);
            fprintf('        [FORWARD PEAK ONLY] Signal peak %.3e Pa; skipping recon.\n', ...
                sensor_signal_peak);
            return;
        end

        % 2. Add broadband electronic noise (after the pulse convolution, before
        %    the band-limit).
        %    noise_amp is the FIXED absolute amplitude (config.noise_amp_Pa) once
        %    the pipeline has been calibrated, else the legacy per-field fraction
        %    of this field's own peak (config.conv_noise_level). A noise_only run
        %    keeps that same amplitude but discards the true signal so ONLY noise
        %    reaches the reconstruction (see resolve_noise_amp below).
        noise_amp = resolve_noise_amp(config, sensor_signal_peak, conv_noise_level);
        if noise_only
            sensorData_noisy = element_noise(noise_amp, size(sensorData_conv), sensor_info);
            fprintf('        [NOISE ONLY] True signal nulled; reconstructing noise alone.\n');
        else
            sensorData_noisy = sensorData_conv + element_noise(noise_amp, size(sensorData_conv), sensor_info);
        end

        % Record the SNR this field actually saw (peak signal / noise amplitude).
        % With a fixed absolute noise this floats field to field; recorded so it
        % can later be correlated with the per-segment gamma index.
        sim_results.noise_amp = noise_amp;
        if noise_amp > 0
            sim_results.snr = sensor_signal_peak / noise_amp;
        else
            sim_results.snr = Inf;
        end

        % 3. Sensor frequency response (band-limit), on signal + noise together
        sensorData_resp = gaussianFilter(sensorData_noisy, FS, 0.35e6, 100, false);

        % 4. Deconvolution. pulse_deconv_gain is the net gain of the pulse
        %    convolution + deconvolution, reused by the 'adjoint' sensor model.
        %    sensorData_conv is noise-free (and not nulled under noise_only), so
        %    the variable lambda matches the signal run's. With conv_use_deconv
        %    false the data keep the pulse smear, so the net gain is the pulse H.
        if conv_use_deconv
            [sensorData_deconv, pulse_deconv_gain] = deconv_pulse_transducer(sensorData_resp, ...
                sensorData_conv, H, FS, noise_amp, conv_deconv_lambda, conv_variable_lambda);
        else
            fprintf('        Deconvolution disabled (conv_use_deconv = false).\n');
            sensorData_deconv = sensorData_resp;
            pulse_deconv_gain = H;
        end

        sensorData          = single(sensorData_deconv);
        sensorData_measured = single(sensorData_deconv);

        fprintf('        Pulse model complete. Signal peak %.3e Pa, noise amp %.3e Pa, SNR %.2f\n', ...
            sensor_signal_peak, noise_amp, sim_results.snr);

        % Diagnostic fast-path (verify_pipeline_simulate): keep the sensor signal
        % after each stage above, then skip the reconstruction.
        if safe_config(config, 'sensor_stages_only', false)
            final_stage = 'Wiener deconvolved';
            if ~conv_use_deconv
                final_stage = 'no deconvolution';
            end
            sim_results.sensor_stages = sensor_stage_traces( ...
                {'raw (forward sim)', 'pulse convolved', 'noise added', ...
                 'sensor band-limit', final_stage}, ...
                {sensorData_cpu, sensorData_conv, sensorData_noisy, ...
                 sensorData_resp, sensorData_deconv}, dt);
            recon_dose = zeros(input_dims_for_crop);
            return;
        end
    else
        % No pulse model: apply only the sensor frequency response.
        sensorData = gaussianFilter(sensorData, FS, 0.35e6, 100, false);

        sensor_signal_peak             = max(abs(sensorData(:)));
        sim_results.sensor_signal_peak = sensor_signal_peak;

        % Calibration fast-path (see the pulse branch above).
        if safe_config(config, 'forward_peak_only', false)
            recon_dose = zeros(input_dims_for_crop);
            fprintf('        [FORWARD PEAK ONLY] Signal peak %.3e Pa; skipping recon.\n', ...
                sensor_signal_peak);
            return;
        end
        if noise_only
            % Null-signal test: keep only electronic noise at the amplitude the
            % real run would have (fixed absolute amp when calibrated, else the
            % legacy per-field fraction). Electronic noise otherwise enters only
            % through the pulse model, so a normal no-pulse run stays noise-free.
            noise_amp     = resolve_noise_amp(config, sensor_signal_peak, conv_noise_level);
            sensorData(:) = 0;
            sensorData    = sensorData + element_noise(noise_amp, size(sensorData), sensor_info);
            sim_results.noise_amp = noise_amp;
            if noise_amp > 0
                sim_results.snr = sensor_signal_peak / noise_amp;
            else
                sim_results.snr = Inf;
            end
            fprintf('        [NOISE ONLY] True signal nulled; reconstructing noise alone (amp %.3e Pa).\n', ...
                noise_amp);
        end

        % Diagnostic fast-path (see the pulse branch above). sensorData_measured
        % still holds the raw forward data here; it is overwritten just below.
        if safe_config(config, 'sensor_stages_only', false)
            sim_results.sensor_stages = sensor_stage_traces( ...
                {'raw (forward sim)', 'sensor band-limit'}, ...
                {sensorData_measured, sensorData}, dt);
            recon_dose = zeros(input_dims_for_crop);
            return;
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

    % ---- Set up live reconstruction figure (opt-in, 'tr', 'adjoint', 'adjoint_tv') ----
    fig_live   = [];
    ax_recon   = [];
    hImg_recon = [];
    hLine_max  = [];
    cz_live    = 1;
    if plot_results && any(strcmp(recon_method, {'tr', 'adjoint', 'adjoint_tv'}))
        if strcmp(recon_method, 'tr')
            num_live_iter = num_tr_iter;
        else
            num_live_iter = num_adj_iter;
        end
        try
            [~, dose_max_idx]               = max(doseGrid(:));
            [~, ~, cz_live]                 = ind2sub([Nx_orig, Ny_orig, Nz_orig], dose_max_idx);
            [fig_live, ax_recon, hImg_recon, hLine_max] = ...
                setup_live_recon_figure(p0, Nx_orig, Ny_orig, cz_live, num_live_iter, config);
        catch ME
            warning('run_single_field_simulation:LiveFigFail', ...
                'Live reconstruction figure setup failed: %s', ME.message);
            fig_live = [];
        end
    end

    switch recon_method
        case 'tr'
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

                    p0_recon = kspaceFirstOrder3DG(kgrid, kmedium_rec, source_tr, sensor_tr, inputArgs{:});

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
                    update_live_recon_figure(fig_live, ax_recon, hImg_recon, hLine_max, ...
                        reconPressure(1:Nx_orig, 1:Ny_orig, cz_live), tr_iter, num_tr_iter, ...
                        converged, conv_max_pressure);

                    if converged
                        break;
                    end

                    % Residual correction for next iteration
                    if tr_iter < num_tr_iter
                        source_resid    = struct();
                        source_resid.p0 = reconPressure;
                        sensorDataRecon = kspaceFirstOrder3DG(kgrid, kmedium_rec, ...
                            source_resid, sensor, inputArgs{:});
                        % Element-average the re-forward like the measured data
                        [~, sensorDataRecon] = apply_element_averaging( ...
                            double(gather(sensorDataRecon)), sensor_info);
                        sensorData = sensorData + (sensorData_measured - single(sensorDataRecon));
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

        case 'ubp'
            fprintf('        Running Universal Back-Projection reconstruction...\n');
            try
                ubp_tic = tic;
                ubp_opts = struct( ...
                    'use_elements', safe_config(config, 'ubp_use_elements', true), ...
                    'band_limit',   safe_config(config, 'ubp_band_limit',   'auto'), ...
                    'aperture_cos', safe_config(config, 'ubp_aperture_cos', 0), ...
                    'area_weight',  safe_config(config, 'ubp_area_weight',  true), ...
                    'interp',       safe_config(config, 'ubp_interp',       'linear'));
                % Non-ideal sensor (determine_sensor_mask) corrections happen
                % inside ubp3d: element averaging, band-limited derivative, and
                % the directivity/obliquity weight are driven by sensor_info.
                if blind_recon, ubp_medium = medium_recon; else, ubp_medium = medium; end
                reconPressure = ubp3d(sensorData, sensor, sensor_info, ubp_medium, ...
                                      Nx, Ny, Nz, dx, dy, dz, dt, ubp_opts);
                reconPressure = reconPressure * correction_factor;
                conv_max_pressure(1) = max(reconPressure(:));
                num_iters_done       = 1;
                tr_time              = toc(ubp_tic);
                fprintf('        UBP complete (%.1f s).\n', tr_time);
                fprintf('        Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
                    min(reconPressure(:)), max(reconPressure(:)));
            catch ME
                warning('run_single_field_simulation:UBPFail', ...
                    'UBP reconstruction failed: %s', ME.message);
                recon_dose = zeros(input_dims_for_crop);
                return;
            end

        case 'adjoint'
            % Gradient descent with the k-Wave adjoint (Arridge et al., Inverse
            % Problems 32, 115012 (2016); k-Wave adjoint reconstruction example).
            % Each iteration:
            %   1. difference = modelled - measured sensor data
            %   2. adjoint run: the time-reversed difference is ADDED at the
            %      sensor points ('additive' source; TR overwrites them). The
            %      final pressure field is the image update (the gradient).
            %   3. step length nu by exact line search: the nu that minimises
            %      ||difference - nu * A(image_update)||. The fixed nu = 0.25 of
            %      the 2D k-Wave example does not carry over to this 3D grid,
            %      because k-Wave rescales additive sources internally.
            %   4. p0_estimate = max(p0_estimate - nu * image_update, 0)
            %   5. forward run on the new estimate -> new modelled data
            % Sensor model: modelled signals get the same frequency response as
            % the measured data (band-limit, pulse convolution + Wiener
            % deconvolution) and the same element averaging, so
            % A is the measurement model and step 2 is its adjoint (up to a
            % constant the line search absorbs). Runs on kmedium_rec, so blind
            % reconstruction works as for TR. No correction_factor: the result
            % is already a least-squares fit to the data, in Pa.
            fprintf('        Running adjoint reconstruction (%d iterations, tol=%.1e)...\n', ...
                num_adj_iter, convergence_tol);

            conv_max_pressure = zeros(num_adj_iter, 1);
            conv_rel_change   = nan(num_adj_iter, 1);

            try
                adj_tic = tic;

                % Net frequency response of the pulse convolution +
                % deconvolution applied to the measured data above.
                if conv_kernel_sigma > 0
                    wiener_gain = pulse_deconv_gain;   % set in the pulse block above
                else
                    wiener_gain = ones(1, size(sensorData_measured, 2));
                end

                % Measured data (already element-averaged above)
                measured = double(sensorData_measured);

                p0_estimate = zeros(gridSize);
                modelled    = zeros(size(measured));   % a zero image predicts zero signals

                for adj_iter = 1:num_adj_iter

                    % 1. Difference between the modelled and measured data
                    difference = modelled - measured;
                    fprintf('          Adjoint iteration %d/%d: data mismatch %.4e\n', ...
                        adj_iter, num_adj_iter, norm(difference(:)) / norm(measured(:)));

                    % 2. Adjoint run (adjoint_backward_model). The final
                    %    pressure field is the image update (the gradient).
                    image_update = adjoint_backward_model(difference, kgrid, kmedium_rec, ...
                        sensor, inputArgs, FS, wiener_gain);

                    % 3. Step length by exact line search (one forward run)
                    update_signals = adjoint_forward_model(image_update, kgrid, kmedium_rec, ...
                        sensor, inputArgs, sensor_info, FS, wiener_gain);
                    update_energy  = sum(update_signals(:) .^ 2);
                    if update_energy == 0
                        fprintf('          Image update predicts no signal; stopping.\n');
                        break;
                    end
                    nu = sum(update_signals(:) .* difference(:)) / update_energy;

                    % 4. Update the image, then apply the positivity condition
                    p0_prev     = p0_estimate;
                    p0_estimate = p0_estimate - nu * image_update;
                    p0_estimate = p0_estimate .* (p0_estimate >= 0);

                    conv_max_pressure(adj_iter) = max(p0_estimate(:));
                    num_iters_done = adj_iter;
                    fprintf('          Step length %.4e, max pressure %.4e Pa\n', ...
                        nu, conv_max_pressure(adj_iter));

                    converged = false;
                    norm_prev = norm(p0_prev(:));
                    if norm_prev > 0
                        conv_rel_change(adj_iter) = norm(p0_estimate(:) - p0_prev(:)) / norm_prev;
                        fprintf('          Rel change: %.4e\n', conv_rel_change(adj_iter));
                        if conv_rel_change(adj_iter) < convergence_tol
                            fprintf('          *** Converged at iteration %d ***\n', adj_iter);
                            converged = true;
                        end
                    end

                    % Live figure update
                    update_live_recon_figure(fig_live, ax_recon, hImg_recon, hLine_max, ...
                        p0_estimate(1:Nx_orig, 1:Ny_orig, cz_live), adj_iter, num_adj_iter, ...
                        converged, conv_max_pressure);

                    if converged
                        break;
                    end

                    % 5. Modelled data of the new estimate, for the next iteration
                    if adj_iter < num_adj_iter
                        modelled = adjoint_forward_model(p0_estimate, kgrid, kmedium_rec, ...
                            sensor, inputArgs, sensor_info, FS, wiener_gain);
                    end
                end

                reconPressure = p0_estimate;
                tr_time       = toc(adj_tic);
                fprintf('        Adjoint reconstruction complete (%.1f s).\n', tr_time);
                fprintf('        Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
                    min(reconPressure(:)), max(reconPressure(:)));
            catch ME
                warning('run_single_field_simulation:AdjointFail', ...
                    'Adjoint reconstruction failed: %s', ME.message);
                recon_dose = zeros(input_dims_for_crop);
                return;
            end

        case 'adjoint_tv'
            % 'adjoint' with a total-variation penalty (Huang et al., IEEE TMI
            % 32, 1097 (2013), eq. 43 and Appendix B). Minimises
            %     J(p) = ||A p - d||^2 + lambda * TV(p),   p >= 0
            % A = adjoint_forward_model and its adjoint = adjoint_backward_model,
            % the same operators as 'adjoint'. TV = eq. 44 (calculate_tv).
            % Setup:
            %   a. Power iteration, started from the back-projected data. The
            %      k-Wave adjoint is the true adjoint times a constant (k-Wave
            %      rescales additive sources); a dot-product test measures it
            %      and it is divided out. The largest eigenvalue of A'A gives
            %      the Lipschitz constant Lip = 2 * max eigenvalue (step 1/Lip).
            %      Both depend on the operator A, not on this field (the field
            %      is only the starting image), so they are cached in
            %      processed/simulated_doses/adjoint_power_cache_<hash>.mat (adjoint_power_hash:
            %      reconstruction medium, grid and dt/Nt, sensor, PML, sensor
            %      response, power_iter) and reused on the next run.
            %   b. lambda = adjoint_tv_lambda * Lip * peak of the first
            %      unregularised step, so adjoint_tv_lambda is dimensionless and
            %      one value works for weak and strong fields.
            % Each iteration (FISTA, Algorithm 1):
            %   1. gradient of the data term at the look-ahead point y:
            %      2 A'(A y - d)
            %   2. gradient step y - gradient / Lip, then TV denoising with
            %      positivity (tv_denoise = F_Dnoise, Algorithm 2)
            %   3. forward run on the new estimate -> modelled data and J
            %   4. momentum: y = p + (t-1)/t_next * (p - p_prev). A is linear,
            %      so A y is the same combination of modelled data (no k-Wave run)
            % k-Wave runs: 1 + 2 * adjoint_tv_power_iter for setup (1 when
            % cached), 2 per iteration. Runs on kmedium_rec (blind works as for TR). No
            % correction_factor (a least-squares fit in Pa, as 'adjoint').
            fprintf('        Running adjoint_tv reconstruction (%d iterations, tol=%.1e)...\n', ...
                num_adj_iter, convergence_tol);

            conv_max_pressure = zeros(num_adj_iter, 1);
            conv_rel_change   = nan(num_adj_iter, 1);
            conv_cost         = nan(num_adj_iter, 1);

            try
                adj_tic = tic;

                % Net frequency response of the pulse convolution +
                % deconvolution applied to the measured data above.
                if conv_kernel_sigma > 0
                    wiener_gain = pulse_deconv_gain;   % set in the pulse block above
                else
                    wiener_gain = ones(1, size(sensorData_measured, 2));
                end

                % Measured data (already element-averaged above)
                measured = double(sensorData_measured);

                % a. Power iteration: adjoint scale and largest eigenvalue of A'A
                back_projection = adjoint_backward_model(measured, kgrid, kmedium_rec, ...
                    sensor, inputArgs, FS, wiener_gain);

                % Cache file in the simulation folder (no caching without a path)
                power_cache_file = '';
                if isfield(config, 'working_dir') && isfield(cbct_resampled, 'patient_id') ...
                        && isfield(cbct_resampled, 'session')
                    power_hash = adjoint_power_hash(kmedium_rec, kgrid, sensor, sensor_info, ...
                        pml_size, wiener_gain, tv_power_iter);
                    power_cache_file = fullfile(get_simulation_dir(config.working_dir, ...
                        cbct_resampled.patient_id, cbct_resampled.session), ...
                        sprintf('adjoint_power_cache_%s.mat', power_hash));
                end

                if ~tv_force_power_recalc && ~isempty(power_cache_file) && isfile(power_cache_file)
                    cached         = load(power_cache_file, 'adjoint_scale', 'max_eigenvalue');
                    adjoint_scale  = cached.adjoint_scale;
                    max_eigenvalue = cached.max_eigenvalue;
                    fprintf('          Power iteration loaded from %s: adjoint scale %.4e, max eigenvalue %.4e\n', ...
                        power_cache_file, adjoint_scale, max_eigenvalue);
                else
                    power_image = back_projection / norm(back_projection(:));
                    for power_iter = 1:tv_power_iter
                        power_signals = adjoint_forward_model(power_image, kgrid, kmedium_rec, ...
                            sensor, inputArgs, sensor_info, FS, wiener_gain);
                        round_trip    = adjoint_backward_model(power_signals, kgrid, kmedium_rec, ...
                            sensor, inputArgs, FS, wiener_gain);
                        % True adjoint: sum(x .* A'(A x)) == sum((A x).^2)
                        adjoint_scale  = sum(power_image(:) .* round_trip(:)) / sum(power_signals(:) .^ 2);
                        max_eigenvalue = norm(round_trip(:)) / abs(adjoint_scale);   % power_image has unit length
                        power_image    = round_trip / norm(round_trip(:));
                        fprintf('          Power iteration %d/%d: adjoint scale %.4e, max eigenvalue %.4e\n', ...
                            power_iter, tv_power_iter, adjoint_scale, max_eigenvalue);
                    end
                    if ~isempty(power_cache_file)
                        try
                            if ~isfolder(fileparts(power_cache_file))
                                mkdir(fileparts(power_cache_file));
                            end
                            save(power_cache_file, 'adjoint_scale', 'max_eigenvalue');
                            fprintf('          Power iteration saved to %s\n', power_cache_file);
                        catch ME
                            warning('run_single_field_simulation:PowerCacheSaveFail', ...
                                'Could not save power iteration cache %s: %s', power_cache_file, ME.message);
                        end
                    end
                end
                % Power iteration approaches the eigenvalue from below: 10% margin
                lipschitz = 2 * 1.1 * max_eigenvalue;

                % b. TV weight. denoise_weight = lambda / Lip is what tv_denoise uses.
                first_step     = 2 * (back_projection / adjoint_scale) / lipschitz;
                denoise_weight = tv_lambda * max(first_step(:));
                lambda_abs     = denoise_weight * lipschitz;
                fprintf('          Lipschitz constant %.4e, TV lambda %.4e (denoise weight %.4e Pa)\n', ...
                    lipschitz, lambda_abs, denoise_weight);
                clear back_projection first_step power_image power_signals round_trip

                p0_estimate    = zeros(gridSize);
                p0_prev        = zeros(gridSize);
                p0_lookahead   = zeros(gridSize);           % y
                modelled_prev  = zeros(size(measured));
                modelled_ahead = zeros(size(measured));     % A y (a zero image predicts zero signals)
                fista_t        = 1;

                for adj_iter = 1:num_adj_iter

                    % 1. Gradient of the data term at y
                    difference    = modelled_ahead - measured;
                    data_gradient = 2 * adjoint_backward_model(difference, kgrid, kmedium_rec, ...
                        sensor, inputArgs, FS, wiener_gain) / adjoint_scale;

                    % 2. Gradient step, then TV denoising with positivity
                    p0_estimate = tv_denoise(p0_lookahead - data_gradient / lipschitz, ...
                        denoise_weight, tv_inner_iter);

                    % 3. Modelled data of the new estimate and the objective J
                    modelled  = adjoint_forward_model(p0_estimate, kgrid, kmedium_rec, ...
                        sensor, inputArgs, sensor_info, FS, wiener_gain);
                    data_term = sum((modelled(:) - measured(:)) .^ 2);
                    tv_term   = lambda_abs * calculate_tv(p0_estimate);
                    conv_cost(adj_iter) = data_term + tv_term;

                    conv_max_pressure(adj_iter) = max(p0_estimate(:));
                    num_iters_done = adj_iter;
                    fprintf(['          Adjoint-TV iteration %d/%d: data mismatch %.4e, ' ...
                             'cost %.4e (TV %.1f%%), max pressure %.4e Pa\n'], ...
                        adj_iter, num_adj_iter, sqrt(data_term) / norm(measured(:)), ...
                        conv_cost(adj_iter), 100 * tv_term / conv_cost(adj_iter), ...
                        conv_max_pressure(adj_iter));

                    converged = false;
                    norm_prev = norm(p0_prev(:));
                    if norm_prev > 0
                        conv_rel_change(adj_iter) = norm(p0_estimate(:) - p0_prev(:)) / norm_prev;
                        fprintf('          Rel change: %.4e\n', conv_rel_change(adj_iter));
                        if conv_rel_change(adj_iter) < convergence_tol
                            fprintf('          *** Converged at iteration %d ***\n', adj_iter);
                            converged = true;
                        end
                    end

                    % Live figure update
                    update_live_recon_figure(fig_live, ax_recon, hImg_recon, hLine_max, ...
                        p0_estimate(1:Nx_orig, 1:Ny_orig, cz_live), adj_iter, num_adj_iter, ...
                        converged, conv_max_pressure);

                    if converged
                        break;
                    end

                    % 4. Momentum step (y and A y)
                    fista_t_next   = (1 + sqrt(1 + 4 * fista_t^2)) / 2;
                    momentum       = (fista_t - 1) / fista_t_next;
                    p0_lookahead   = p0_estimate + momentum * (p0_estimate - p0_prev);
                    modelled_ahead = modelled + momentum * (modelled - modelled_prev);
                    p0_prev        = p0_estimate;
                    modelled_prev  = modelled;
                    fista_t        = fista_t_next;
                end

                reconPressure         = p0_estimate;
                sim_results.conv_cost = conv_cost(1:max(num_iters_done, 1));
                tr_time               = toc(adj_tic);
                fprintf('        Adjoint-TV reconstruction complete (%.1f s).\n', tr_time);
                fprintf('        Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
                    min(reconPressure(:)), max(reconPressure(:)));
            catch ME
                warning('run_single_field_simulation:AdjointTVFail', ...
                    'Adjoint-TV reconstruction failed: %s', ME.message);
                recon_dose = zeros(input_dims_for_crop);
                return;
            end

        otherwise
            error('run_single_field_simulation:UnknownReconMethod', ...
                'Unknown reconstruction_method: "%s" (use ''tr'', ''das'', ''ubp'', ''adjoint'', or ''adjoint_tv'')', recon_method);
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

    %% ============ RESAMPLE RECON BACK TO NATIVE DOSE DIMENSIONS ============
    %  With a non-unity downscale_factor the recon currently sits on the coarse
    %  simulation grid. Resample it back up to the native dose dimensions so the
    %  caller can compare it against the full-resolution truth (gamma/SSIM) with
    %  no further resampling. Clamp negatives from the cubic interpolation.
    %  Diagnostic/pressure arrays above stay on the simulation grid.
    if downscale_factor ~= 1 && ~isequal(size(recon_dose), native_dims)
        recon_dose = max(imresize3(recon_dose, native_dims), 0);
        fprintf('        Resampled recon_dose back to native dims [%d %d %d]\n', native_dims);
        if return_diag
            sim_results.recon_dose = recon_dose;
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


function amp = resolve_noise_amp(config, signal_peak, conv_noise_level)
%RESOLVE_NOISE_AMP  Electronic-noise amplitude (Pa) added to the sensor data.
%   Returns the FIXED absolute amplitude config.noise_amp_Pa when it is set to a
%   positive scalar -- the calibrated blanket noise used across every
%   beam/segment, so the SNR floats field to field. Otherwise falls back to the
%   legacy behavior: a fraction (conv_noise_level) of THIS field's own signal
%   peak, which pins every field to the same SNR = 1/conv_noise_level.
%   calibrate_noise_amp.m computes noise_amp_Pa for a target mean SNR.
    if isfield(config, 'noise_amp_Pa') && ~isempty(config.noise_amp_Pa) ...
            && isnumeric(config.noise_amp_Pa) && isscalar(config.noise_amp_Pa) ...
            && config.noise_amp_Pa > 0
        amp = double(config.noise_amp_Pa);
    else
        amp = conv_noise_level * signal_peak;
    end
end


function stages = sensor_stage_traces(names, arrays, dt)
%SENSOR_STAGE_TRACES Mean |p(t)| over all sensor points for each processing
%  stage (arrays{k} is [num_sensor_points x Nt]). The absolute value is taken
%  before averaging so bipolar signals from different sensors do not cancel.
    stages = struct('names', {names}, 'mean_abs', [], 'dt', dt);
    for k = 1:numel(arrays)
        stages.mean_abs(k, :) = mean(abs(double(gather(arrays{k}))), 1);
    end
end


function noise = element_noise(noise_amp, sz, sensor_info)
%ELEMENT_NOISE Electronic noise (std noise_amp, Pa) for sensor data of size sz.
%  Real electronic noise belongs to an element's channel, not to each voxel, so
%  every voxel of an element gets the SAME noise trace and the data stay one
%  signal per element. Uses sensor_info.voxel_element_idx (keyed by
%  find(sensor.mask), like apply_element_averaging). Sensors without it (plane,
%  sphere, box, legacy element_map) get one independent trace per voxel.
    elem = [];
    if isfield(sensor_info, 'voxel_element_idx') && ~isempty(sensor_info.voxel_element_idx) ...
            && isfield(sensor_info, 'num_elements') && sensor_info.num_elements > 0
        elem = double(sensor_info.voxel_element_idx(:));
        elem = elem(elem > 0);
    end
    if numel(elem) == sz(1)
        noise_per_element = noise_amp * randn(sensor_info.num_elements, sz(2));
        noise = noise_per_element(elem, :);
    else
        noise = noise_amp * randn(sz);
    end
end


function [deconvolved, deconv_gain] = deconv_pulse_transducer(data, clean_conv, H, FS, ...
        noise_amp, conv_deconv_lambda, use_variable_lambda)
%DECONV_PULSE_TRANSDUCER Wiener deconvolution of the band-limited sensor data.
%  data        - [num_sensor_points x Nt] noisy, band-limited sensor data
%  clean_conv  - [num_sensor_points x Nt] noise-free pulse-convolved data; only
%                its energy is used (the signal level for the variable lambda)
%  H           - [1 x Nt] FFT of the pulse kernel
%
%  use_variable_lambda = false: Wiener filter on the pulse only, fixed lambda
%      G = conj(H) / (|H|^2 + conv_deconv_lambda)
%  The transducer band-limit stays in the data.
%
%  use_variable_lambda = true: Wiener filter on the whole system Hs = H * T
%  (pulse times transducer response T), with lambda(f) = N(f) / S(f):
%      N(f) = Nt * noise_amp^2 * T(f)^2   noise power (added before the transducer)
%      S(f) = S0 * T(f)^2                 ASSUMED signal power, shaped like T
%  T^2 cancels, so lambda = Nt * noise_amp^2 / S0: flat in frequency, set by
%  this field's SNR. S0 comes from the clean signal energy (Parseval), giving
%      lambda = noise_amp^2 * sum(|Hs|^2) / mean(sum(clean_conv.^2, 2))
%  Each frequency is boosted only as far as its SNR, |Hs|^2 / lambda, allows.
%
%  deconv_gain - [1 x Nt] net gain H * G of the pulse convolution and this
%                filter: what the 'adjoint' sensor model applies after the
%                band-limit (|H|^2 / (|H|^2 + lambda) for the fixed-lambda case).
    lambda = conv_deconv_lambda;
    if use_variable_lambda
        warning('deconv_pulse_transducer:AssumedSignalSpectrum', ...
            ['conv_variable_lambda assumes the signal power spectrum has the shape of ' ...
             'the sensor frequency response; the true spectrum differs, so the noise ' ...
             'weighting is approximate.']);

        % Transducer response T(f): the same gaussianFilter call as the data,
        % applied to a unit impulse. The filter is zero-phase, so T is real.
        impulse  = [1, zeros(1, size(data, 2) - 1)];
        T        = real(fft(gaussianFilter(impulse, FS, 0.35e6, 100, false)));
        H_system = H .* T;

        % With no noise the noise-based lambda is 0 (an unstable pure inverse),
        % so conv_deconv_lambda is kept.
        if noise_amp > 0
            lambda = noise_amp^2 * sum(abs(H_system).^2) / mean(sum(clean_conv.^2, 2));
        end
        fprintf('        Variable-lambda Wiener (pulse x transducer): lambda=%.2e\n', lambda);
    else
        H_system = H;
    end

    G           = conj(H_system) ./ (abs(H_system).^2 + lambda);
    deconvolved = real(ifft(fft(data, [], 2) .* G, [], 2));
    deconv_gain = real(H .* G);
end


function signals = apply_sensor_response(signals, FS, wiener_gain)
%APPLY_SENSOR_RESPONSE Frequency response the measured sensor data went through:
%  the transducer band-limit (same gaussianFilter call as the forward chain),
%  then the pulse convolution + deconvolution, whose net gain wiener_gain comes
%  from deconv_pulse_transducer (all ones without the pulse model; the pulse H
%  itself when conv_use_deconv is false). The band-limit is real and zero-phase,
%  so the adjoint of this operator is the same call with conj(wiener_gain),
%  which the 'adjoint' branch passes on the adjoint side.
    signals = gaussianFilter(double(gather(signals)), FS, 0.35e6, 100, false);
    signals = real(ifft(fft(signals, [], 2) .* wiener_gain, [], 2));
end


function signals = adjoint_forward_model(p0, kgrid, kmedium, sensor, inputArgs, sensor_info, FS, wiener_gain)
%ADJOINT_FORWARD_MODEL Model of the measured data for the 'adjoint' branch:
%  k-Wave forward run from p0, the sensor frequency response, then averaging
%  over each array element (apply_element_averaging passes the data through
%  unchanged when the sensor has no element map). Returns
%  [num_sensor_voxels x Nt], each voxel carrying its element's signal.
    source       = struct('p0', p0);
    signals      = kspaceFirstOrder3DG(kgrid, kmedium, source, sensor, inputArgs{:});
    signals      = apply_sensor_response(signals, FS, wiener_gain);
    [~, signals] = apply_element_averaging(signals, sensor_info);
end


function p_final = adjoint_backward_model(signals, kgrid, kmedium, sensor, inputArgs, FS, wiener_gain)
%ADJOINT_BACKWARD_MODEL Adjoint of adjoint_forward_model, up to a constant
%  (k-Wave rescales additive sources), for the 'adjoint' and 'adjoint_tv'
%  branches. The signals go through the adjoint of the sensor response (conj
%  of the gain; unchanged when the gain is real, i.e. with deconvolution),
%  then are time-reversed and ADDED at the sensor points ('additive' source;
%  TR overwrites them) as the adjoint source of Appendix B in Arridge et al.
%  2016. Returns the final pressure field [Nx x Ny x Nz]. No element
%  averaging: the inputs are already element-averaged, and averaging is its
%  own adjoint, so it would leave them unchanged.
    time_reversed = fliplr(apply_sensor_response(signals, FS, conj(wiener_gain)));
    source        = struct();
    source.p_mask = sensor.mask;
    source.p      = [time_reversed(:, 1), time_reversed(:, 1), time_reversed(:, 1:end-1)] + ...
                    [zeros(size(time_reversed(:, 1))), time_reversed(:, 1:end-1), 2 * time_reversed(:, end)];
    source.p      = source.p(:, 1:kgrid.Nt);   % k-Wave only uses the first Nt points
    source.p_mode = 'additive';

    recorder        = struct();
    recorder.mask   = ones(size(sensor.mask));
    recorder.record = {'p_final'};

    out     = kspaceFirstOrder3DG(kgrid, kmedium, source, recorder, inputArgs{:});
    p_final = reshape(double(gather(out.p_final)), size(sensor.mask));
end


function hash_hex = adjoint_power_hash(kmedium, kgrid, sensor, sensor_info, pml_size, wiener_gain, num_power_iter)
%ADJOINT_POWER_HASH 8-character MD5 key of the cached 'adjoint_tv' power
%  iteration: everything the operator A (adjoint_forward_model) depends on --
%  reconstruction medium, grid and time step, sensor mask and element map, PML
%  size, sensor frequency response -- plus the number of power iterations.
%  The field dose is left out: it is only the power iteration's starting image.
    parts = {kmedium.density, kmedium.sound_speed, kmedium.alpha_coeff, kmedium.alpha_power, ...
             [kgrid.Nx, kgrid.Ny, kgrid.Nz, kgrid.dx, kgrid.dy, kgrid.dz, kgrid.dt, kgrid.Nt], ...
             find(sensor.mask), pml_size, real(wiener_gain), imag(wiener_gain), num_power_iter};
    if isfield(sensor_info, 'voxel_element_idx'), parts{end+1} = sensor_info.voxel_element_idx; end
    if isfield(sensor_info, 'element_map'),       parts{end+1} = sensor_info.element_map;       end

    md = java.security.MessageDigest.getInstance('MD5');
    for k = 1:numel(parts)
        if ~isempty(parts{k})
            md.update(typecast(double(gather(parts{k}(:))), 'uint8'));
        end
    end
    raw      = typecast(md.digest(), 'uint8');
    hash_hex = sprintf('%02x', raw(1:4));
end


function tv = calculate_tv(p)
%CALCULATE_TV Total variation of a 3D image, Huang et al. 2013 eq. 44:
%  the sum over voxels n of sqrt(sum over the 3 dimensions of
%  (p(n) - p(neighbour before n))^2). The first voxel along a dimension has
%  no neighbour before it, so that term is 0.
    [d1, d2, d3] = tv_backward_diff(p);
    tv = sum(sqrt(d1 .^ 2 + d2 .^ 2 + d3 .^ 2), 'all');
end


function p = tv_denoise(b, weight, num_iter)
%TV_DENOISE F_Dnoise of Huang et al. 2013 (Appendix B, Algorithm 2; the FGP
%  algorithm of Beck & Teboulle 2009, extended to 3D). Solves
%      min over p >= 0 of  ||p - b||^2 + 2 * weight * calculate_tv(p)
%  through its dual: one vector (q1, q2, q3) per voxel, length <= 1. Each
%  iteration takes a projected gradient step on the dual (step 1/(12*weight);
%  12 bounds ||D||^2 in 3D, 4 per dimension, D = tv_backward_diff), with
%  momentum. The image is p = max(b - weight * D'(q), 0).
%  weight <= 0 means no TV: only the positivity is applied.
    if weight <= 0
        p = max(b, 0);
        return;
    end

    q1 = zeros(size(b));  q2 = zeros(size(b));  q3 = zeros(size(b));   % dual
    r1 = zeros(size(b));  r2 = zeros(size(b));  r3 = zeros(size(b));   % dual look-ahead
    t  = 1;
    for k = 1:num_iter
        p = max(b - weight * tv_backward_diff_transpose(r1, r2, r3), 0);
        [d1, d2, d3] = tv_backward_diff(p);
        s1 = r1 + d1 / (12 * weight);
        s2 = r2 + d2 / (12 * weight);
        s3 = r3 + d3 / (12 * weight);

        % Project each voxel's vector onto the unit ball
        vec_len = max(1, sqrt(s1 .^ 2 + s2 .^ 2 + s3 .^ 2));
        s1 = s1 ./ vec_len;
        s2 = s2 ./ vec_len;
        s3 = s3 ./ vec_len;

        t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
        r1 = s1 + ((t - 1) / t_next) * (s1 - q1);
        r2 = s2 + ((t - 1) / t_next) * (s2 - q2);
        r3 = s3 + ((t - 1) / t_next) * (s3 - q3);
        q1 = s1;  q2 = s2;  q3 = s3;
        t  = t_next;
    end
    p = max(b - weight * tv_backward_diff_transpose(q1, q2, q3), 0);
end


function [d1, d2, d3] = tv_backward_diff(p)
%TV_BACKWARD_DIFF Difference to the neighbour before each voxel (eq. 44),
%  p(n) - p(n-1) along array dims 1, 2 and 3; 0 at the first voxel.
    d1 = zeros(size(p));  d2 = zeros(size(p));  d3 = zeros(size(p));
    d1(2:end, :, :) = p(2:end, :, :) - p(1:end-1, :, :);
    d2(:, 2:end, :) = p(:, 2:end, :) - p(:, 1:end-1, :);
    d3(:, :, 2:end) = p(:, :, 2:end) - p(:, :, 1:end-1);
end


function p = tv_backward_diff_transpose(d1, d2, d3)
%TV_BACKWARD_DIFF_TRANSPOSE Transpose of tv_backward_diff: for any image x,
%  sum(D(x) .* (d1, d2, d3)) == sum(x .* tv_backward_diff_transpose(d1, d2, d3)).
%  Each difference x(n) - x(n-1) adds its weight to voxel n and subtracts it
%  from voxel n-1.
    p = zeros(size(d1));
    p(2:end, :, :)   = p(2:end, :, :)   + d1(2:end, :, :);
    p(1:end-1, :, :) = p(1:end-1, :, :) - d1(2:end, :, :);
    p(:, 2:end, :)   = p(:, 2:end, :)   + d2(:, 2:end, :);
    p(:, 1:end-1, :) = p(:, 1:end-1, :) - d2(:, 2:end, :);
    p(:, :, 2:end)   = p(:, :, 2:end)   + d3(:, :, 2:end);
    p(:, :, 1:end-1) = p(:, :, 1:end-1) - d3(:, :, 2:end);
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
%  'tr' and 'adjoint' branches. Panels: (1) initial p0 slice, (2) current
%  recon, (3) live max-pressure convergence.

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


function update_live_recon_figure(fig_live, ax_recon, hImg_recon, hLine_max, ...
        recon_slice, iter, N_iter, converged, conv_max_pressure)
%UPDATE_LIVE_RECON_FIGURE Redraw the live figure after one 'tr' or 'adjoint'
%  iteration: current recon slice (X-Y at cz_live), panel title, and the
%  max-pressure convergence line. Does nothing when the figure is not open.

    if isempty(fig_live) || ~ishandle(fig_live)
        return;
    end
    try
        recon_slice_crop = gather(squeeze(recon_slice)');
        set(hImg_recon, 'CData', recon_slice_crop);
        caxis(ax_recon, [0, max(recon_slice_crop(:)) + eps]);
        if converged
            title(ax_recon, ...
                sprintf('Reconstructed p_0 (iter %d CONVERGED)', iter), ...
                'FontWeight', 'bold', 'Color', [0, 0.55, 0]);
        else
            title(ax_recon, ...
                sprintf('Reconstructed p_0 (iter %d / %d)', ...
                    iter, N_iter), 'FontWeight', 'bold');
        end
        set(hLine_max, 'XData', 1:iter, ...
            'YData', conv_max_pressure(1:iter));
        drawnow;
    catch
        % Ignore plot errors; don't break the recon
    end
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
