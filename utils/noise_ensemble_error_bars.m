function ens = noise_ensemble_error_bars(config, ref_truth, spacing_mm, sct, gantry_angle, beam_meta, varargin)
%NOISE_ENSEMBLE_ERROR_BARS  Mean/std gamma pass rate of the noise-only null.
%
%   ens = noise_ensemble_error_bars(config, ref_truth, spacing_mm, sct, ...
%                                   gantry_angle, beam_meta, Name, Value, ...)
%
%   PURPOSE
%     A single noise-only reconstruction (true acoustic signal nulled, only one
%     random electronic-noise draw reconstructed) has a gamma pass rate that
%     jitters run to run (empirically ~18-20%). This helper builds the null band
%     by running an ENSEMBLE of such reconstructions and returning the MEAN and
%     STANDARD DEVIATION of the pass rate, so the null hypothesis can be drawn as
%     an error bar instead of a single noisy point.
%
%   FORWARD RUNS ONCE  (the efficiency point)
%     The forward k-Wave simulation is run a SINGLE time (build_forward_bundle),
%     capturing everything needed to re-reconstruct. Every realization then only
%     redraws electronic noise and runs ONE reconstruction (redraw_noisy_deconv +
%     reconstruct_recon_dose). At ~8 s per k-Wave stage that roughly doubles the
%     samples achievable in a fixed time budget versus re-running the forward.
%
%     This is faithful because for the noise-only null the forward's true signal
%     is discarded; its ONLY dose-dependent contribution is the noise amplitude
%     (noise_amp = conv_noise_level * peak true signal), and that overall scale is
%     divided straight back out by the least-squares normalization before gamma.
%     So only the RECONSTRUCTION geometry/medium matters, which is fixed. The
%     bundle is therefore built on the reconstruction (reference) geometry passed
%     in via sct -- a single-medium bundle, no blind forward/recon split needed.
%
%   CACHING  (the reason this exists)
%     The cache key is compute_sim_config_hash(config), which depends on the
%     sensor placement + noise + TR + reconstruction settings but NOT on the
%     dose/CBCT file. So the SAME cached ensemble is reused for every beam and
%     segment in a session -- you only pay the compute once.
%
%     >>> NOTE -- REEXAMINE IF THE SENSOR POSITION CHANGES <<<
%     Reuse across beams is only valid because the ultrasound array sits in the
%     SAME place for every beam/segment: determine_sensor_mask builds the mask
%     from the session-level beam_metadata (avoiding ALL beams' jaw projections),
%     not from the single field being reconstructed. If sensor placement ever
%     becomes per-beam -- or determine_sensor_mask is changed so the mask depends
%     on the specific dose -- this cache is no longer valid across beams. Then add
%     the sensor mask into the cache key (sensor_mask_nnz is stored below as a
%     starting point) and drop the "compute once" assumption.
%
%   INPUTS
%     config       - simulation CONFIG (the same one used for the real recon). Its
%                    sim fields set the cache key; .noise_only/.plot_results are
%                    forced here. Needs .working_dir/.patient_id/.session, and the
%                    engine knobs (sensor, recon, pulse, noise, downscale).
%     ref_truth    - 3D reference-CT truth dose the noise recon is gamma-compared
%                    against (supplies the 10% gamma eval mask). On the INPUT grid.
%     spacing_mm   - [dx dy dz] mm (kept for the record; gamma uses the bundle's
%                    own spacing in case of downscaling).
%     sct          - reconstruction-geometry CBCT struct (cubeHU/bodyMask/spacing/
%                    couchMask/origin) -- e.g. run_standalone_field's out.sct for
%                    the reference field.
%     gantry_angle - gantry angle of the reference field (for sensor placement).
%     beam_meta    - FULL plan beam_metadata (all beams) for determine_sensor_mask.
%
%   NAME/VALUE OPTIONS (defaults)
%     'TimeBudgetMin'  30    - stop once the projected time for one more realization
%                              would exceed this many minutes (bundle build included).
%     'MinSamples'     20    - always run at least this many (for a usable std).
%     'MaxSamples'     1000  - hard cap regardless of the time budget.
%     'ForceRecompute' false - ignore any cache and recompute from scratch.
%     'Normalize'      config.normalize (else true) - lsq-scale each noise recon to
%                              ref_truth before gamma (matches the pipeline; also
%                              what makes the single forward faithful, see above).
%     'GammaCriteria'  {3,3,'3%/3mm'} - single {pct, dta_mm, label} row.
%     'CacheDir'       AnalysisResults/<patient>/<session>
%     'BaseSeed'       10000 - realization r uses rng seed BaseSeed + r (reproducible).
%
%   OUTPUT (struct ens)
%     .mean_pass_rate .std_pass_rate .pass_rates .num_samples .criteria
%     .config_hash .from_cache .cache_file .conv_noise_level .num_tr_iter
%     .reconstruction_method .sensor_placement_method .sensor_mask_nnz .timestamp
%     (No 3D volume is returned or cached -- the deliverable is the null band.)
%
%   See also: run_standalone_comparison, run_standalone_field, compute_gamma,
%             compute_sim_config_hash, study_pass_rates_individual (noise ensemble)

    %% ---- Options ----
    p = struct('TimeBudgetMin', 30, 'MinSamples', 20, 'MaxSamples', 1000, ...
               'ForceRecompute', false, 'GammaCriteria', {{3, 3, '3%/3mm'}}, ...
               'CacheDir', '', 'Normalize', [], 'BaseSeed', 10000);
    for i = 1:2:numel(varargin)
        p.(varargin{i}) = varargin{i+1};
    end
    if isempty(p.Normalize)
        p.Normalize = ~isfield(config, 'normalize') || config.normalize;
    end

    config.noise_only   = true;    % this helper only makes sense for the null
    config.plot_results = false;    % never pop figures inside the ensemble loop
    crit = p.GammaCriteria;

    %% ---- Cache location + key ----
    config_hash = compute_sim_config_hash(config);
    if isempty(p.CacheDir)
        p.CacheDir = fullfile(config.working_dir, 'AnalysisResults', ...
            config.patient_id, config.session);
    end
    if ~isfolder(p.CacheDir)
        mkdir(p.CacheDir);
    end
    cache_file = fullfile(p.CacheDir, sprintf('noise_ensemble_%s_%s_%s.mat', ...
        config.patient_id, config.session, config_hash));

    %% ---- Load cache (default behaviour) ----
    if ~p.ForceRecompute && isfile(cache_file)
        S = load(cache_file);
        ens = S.ens;
        ens.from_cache = true;
        ens.cache_file = cache_file;
        fprintf(['[noise ensemble] Loaded cache %s\n' ...
                 '                 Null gamma %.2f +/- %.2f %% over %d samples.\n'], ...
                cache_file, ens.mean_pass_rate, ens.std_pass_rate, ens.num_samples);
        return;
    end

    %% ---- Forward simulation ONCE ----
    fprintf(['[noise ensemble] No cache -- computing null ensemble.\n' ...
             '                 Budget %d min, %d-%d samples, key %s\n'], ...
            p.TimeBudgetMin, p.MinSamples, p.MaxSamples, config_hash);

    t_start = tic;
    B = build_forward_bundle(double(ref_truth), sct, gantry_angle, beam_meta, config, 'noise-null');
    B.sensorData_resp(:) = 0;   % NOISE-ONLY null: discard true signal, keep noise_amp
    ref_embed  = embed_on_grid(double(ref_truth), B.gridSize_orig, B.embed_offset, 0);
    sensor_nnz = nnz(B.sensor.mask);
    bundle_time = toc(t_start);
    fprintf('[noise ensemble] forward bundle built (%.1f s). Redrawing noise...\n', bundle_time);

    %% ---- Redraw noise + reconstruct per realization ----
    budget_s   = p.TimeBudgetMin * 60;
    pass_rates = [];
    n_fail     = 0;

    while true
        r = numel(pass_rates) + 1;
        try
            sd        = redraw_noisy_deconv(B, p.BaseSeed + r);   % pure-noise trace
            recon_raw = reconstruct_recon_dose(B, sd);            % ONE reconstruction
            recon_g   = recon_raw;
            if p.Normalize
                recon_g = recon_g * least_squares_gain(ref_embed, recon_g);
            end
            gr = compute_gamma(ref_embed, recon_g, B.spacing_mm, ...
                'Criteria', {crit{1}, crit{2}, crit{3}});
            pass_rates(end+1) = gr.pass_rates(1);   %#ok<AGROW>
        catch ME
            n_fail = n_fail + 1;
            warning('noise_ensemble_error_bars:RunFail', ...
                'Realization %d failed (%s); skipping.', r, ME.message);
            if n_fail > 5
                error('noise_ensemble_error_bars:TooManyFailures', ...
                    'Aborting: more than 5 noise-only realizations failed.');
            end
            continue;
        end

        elapsed = toc(t_start);
        n_done  = numel(pass_rates);
        per_run = (elapsed - bundle_time) / n_done;   % avg per realization
        fprintf('[noise ensemble] sample %d: pass %.2f%%  (elapsed %.1f min)\n', ...
                n_done, pass_rates(end), elapsed / 60);

        if n_done >= p.MaxSamples
            break;
        end
        if n_done >= p.MinSamples && (elapsed + per_run) > budget_s
            break;
        end
    end

    if isempty(pass_rates)
        error('noise_ensemble_error_bars:NoSamples', ...
            'No noise-only realization succeeded; cannot form an ensemble.');
    end

    %% ---- Assemble + cache the stats (not the volumes) ----
    ens = struct();
    ens.mean_pass_rate          = mean(pass_rates);
    ens.std_pass_rate           = std(pass_rates);
    ens.pass_rates              = pass_rates;
    ens.num_samples             = numel(pass_rates);
    ens.criteria                = crit;
    ens.config_hash             = config_hash;
    ens.conv_noise_level        = getf(config, 'conv_noise_level', NaN);
    ens.num_tr_iter             = getf(config, 'num_time_reversal_iter', NaN);
    ens.reconstruction_method   = getf(config, 'reconstruction_method', '');
    ens.sensor_placement_method = getf(config, 'sensor_placement_method', '');
    ens.sensor_mask_nnz         = sensor_nnz;
    ens.timestamp               = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    save(cache_file, 'ens', '-v7');   % stats only -- small file, no 3D volumes

    ens.from_cache = false;
    ens.cache_file = cache_file;

    fprintf(['[noise ensemble] Done: null gamma %.2f +/- %.2f %% over %d samples ' ...
             '(%.1f min).\n                 Cached to %s\n'], ...
            ens.mean_pass_rate, ens.std_pass_rate, ens.num_samples, ...
            toc(t_start) / 60, cache_file);
end


function v = getf(s, f, d)
%GETF Struct field with default.
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end


%% ========================================================================
%  FORWARD-BUNDLE MACHINERY
%
%  The functions below (build_forward_bundle, pad_medium_p0, apply_grid_expansion,
%  embed_on_grid, redraw_noisy_deconv, reconstruct_recon_dose) are copied from the
%  in-production noise ensemble in study_pass_rates_individual.m (the accepted
%  null-hypothesis methodology; see CLAUDE-SIMULATION_CONTEXT.md). They are
%  duplicated here rather than shared because MATLAB local functions are not
%  callable across files, and extracting them would mean editing the working study
%  driver. Keep them in sync if the study file's versions change.
%
%  IMPORTANT: the medium is built by build_medium_with_bath -> create_acoustic_medium
%  (the single canonical builder, with the in-body air row), NOT the study driver's
%  local create_medium (which omits air and caused the null to disagree with the
%  real recon). reconstruct_recon_dose also masks air pockets exactly as the engine
%  does. This keeps the null on identical physics to run_standalone_field.
%% ========================================================================

function v = embed_on_grid(vol, dims, off, fillval)
%EMBED_ON_GRID  Place an original-grid volume into a (larger) display grid.
    if nargin < 4 || isempty(fillval), fillval = 0; end
    if isempty(vol) || (isequal(size(vol), dims) && all(off == 0))
        v = vol;
        return;
    end
    v = fillval * ones(dims, 'like', double(vol));
    s = size(vol);
    v(off(1) + (1:s(1)), off(2) + (1:s(2)), off(3) + (1:s(3))) = vol;
end


function B = build_forward_bundle(truthDose, sct, gantry_angle, beam_meta, CONFIG, label)
%BUILD_FORWARD_BUNDLE  Dose -> medium -> p0 -> sensor -> k-Wave forward -> pulse
%  model, returning the reconstruction bundle B (see reconstruct_recon_dose).

    fprintf('  [Bundle %s] forward simulation...\n', label);

    doseGrid = double(truthDose);
    if ~isfield(sct, 'cubeHU')
        error('build_forward_bundle:NoHU', 'CBCT struct missing cubeHU.');
    end
    spacing_mm = sct.spacing(:)';
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;

    gridSize = size(doseGrid);
    Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);

    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    % --- Optional downscaling ---
    if CONFIG.downscale_factor ~= 1
        df     = CONFIG.downscale_factor;
        new_Nx = max(1, round(Nx / df));
        new_Ny = max(1, round(Ny / df));
        new_Nz = max(1, round(Nz / df));
        doseGrid   = max(imresize3(doseGrid, [new_Nx, new_Ny, new_Nz]), 0);
        sct.cubeHU = imresize3(sct.cubeHU, [new_Nx, new_Ny, new_Nz]);
        if isfield(sct, 'bodyMask')
            sct.bodyMask = imresize3(single(sct.bodyMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        if isfield(sct, 'couchMask') && ~isempty(sct.couchMask)
            sct.couchMask = imresize3(single(sct.couchMask), [new_Nx, new_Ny, new_Nz], 'nearest') > 0.5;
        end
        spacing_mm = spacing_mm .* ([Nx, Ny, Nz] ./ [new_Nx, new_Ny, new_Nz]);
        dx = spacing_mm(1) / 1000; dy = spacing_mm(2) / 1000; dz = spacing_mm(3) / 1000;
        sct.spacing = spacing_mm;
        Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
        gridSize = [Nx, Ny, Nz];
    end

    % --- Acoustic medium ---
    % Single source of truth: create_acoustic_medium (+ coupling-bath override),
    % identical to run_standalone_field. Do NOT reimplement medium construction --
    % create_acoustic_medium is the only builder that applies the in-body air row.
    medium = build_medium_with_bath(sct, CONFIG);

    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    % --- Initial pressure p0 ---
    meterset       = CONFIG.meterset;
    num_pulses     = ceil(meterset / CONFIG.dose_per_pulse_cGy);
    dose_per_pulse = doseGrid / num_pulses;
    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    p0 = smooth(p0);

    doseThreshold = 0.01 * max(doseGrid(:));
    doseMask      = doseGrid > doseThreshold;
    if ~any(doseMask(:)) || max(p0(:)) == 0
        error('build_forward_bundle:NoDose', 'No significant dose / zero p0 for %s.', label);
    end

    % --- FFT-optimal grid padding ---
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
        [medium, p0] = pad_medium_p0(medium, p0, Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);
        Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
        gridSize = [Nx, Ny, Nz];
    end

    % Track the pre-expansion recon-grid dims and the offset that sensor
    % grid-expansion introduces.
    dims_pre_expand = [Nx_orig, Ny_orig, Nz_orig];
    embed_offset    = [0, 0, 0];

    % --- Sensor placement ---
    sensor      = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    switch CONFIG.sensor_placement_method
        case 'full_plane_anterior'
            sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
            sensor_info_orig = struct('num_elements', 0);

        case 'determine_sensor_mask'
            sct_for_sensor = sct;
            if ~isfield(sct_for_sensor, 'couchMask') || isempty(sct_for_sensor.couchMask)
                sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
            end
            if ~isfield(sct_for_sensor, 'origin') || isempty(sct_for_sensor.origin)
                sct_for_sensor.origin = [0, 0, 0];
            end
            sct_for_sensor.spacing = spacing_mm;

            field_dose_for_sensor              = struct();
            field_dose_for_sensor.dose_Gy      = doseGrid;
            field_dose_for_sensor.gantry_angle = gantry_angle;
            field_dose_for_sensor.origin       = sct_for_sensor.origin;
            field_dose_for_sensor.spacing      = spacing_mm;
            field_dose_for_sensor.dimensions   = [Nx_orig, Ny_orig, Nz_orig];

            [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask( ...
                sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);

            % --- Grid expansion handling (water padding to clear the exclusion
            % zone), then re-run FFT-optimal padding.
            gp = sensor_info_orig.grid_pad;
            if gp.expanded
                [medium, p0, doseGrid, doseMask, sct, medium_orig, ...
                 Nx_orig, Ny_orig, Nz_orig, Nx, Ny, Nz] = ...
                    apply_grid_expansion(medium, p0, doseGrid, doseMask, sct, ...
                        gp, Nx_orig, Ny_orig, Nz_orig, CONFIG);
                gridSize_orig = [Nx_orig, Ny_orig, Nz_orig];
                gridSize      = [Nx, Ny, Nz];
                sensor.mask   = zeros(Nx, Ny, Nz);
                embed_offset  = [gp.y_pre, gp.x_pre, gp.z_pre];
            end

            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));

        otherwise
            error('build_forward_bundle:Sensor', ...
                'Unsupported sensor_placement_method: %s', CONFIG.sensor_placement_method);
    end

    if sum(sensor.mask(:)) == 0
        error('build_forward_bundle:EmptySensor', 'Sensor mask is empty for %s.', label);
    end

    % --- k-Wave grid & medium ---
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    maxC  = max(medium.sound_speed(:));
    minC  = min(medium.sound_speed(medium.sound_speed > 0));
    dt    = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;
    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);

    % Nt_scaling: air (~343 m/s) is the slowest tissue, so it sets minC and
    % inflates simTime/Nt. When air drives minC, shorten the recording length
    % by CONFIG.Nt_scaling (0 = disabled).
    if isfield(CONFIG, 'Nt_scaling') && CONFIG.Nt_scaling ~= 0
        air_c = 343;
        if isfield(CONFIG, 'tissue_tables') && isfield(CONFIG.tissue_tables, 'threshold_2') ...
                && isfield(CONFIG.tissue_tables.threshold_2, 'air')
            air_c = CONFIG.tissue_tables.threshold_2.air.sound_speed;
        end
        if minC <= air_c
            Nt = max(1, ceil(Nt / CONFIG.Nt_scaling));
        end
    end

    kgrid.dt = dt;
    kgrid.Nt = Nt;

    kmedium             = struct();
    kmedium.density     = medium.density;
    kmedium.sound_speed = medium.sound_speed;
    kmedium.alpha_coeff = medium.alpha_coeff;
    kmedium.alpha_power = 1.1;

    if CONFIG.use_gpu
        try
            gpuDevice;
            dataCast = 'gpuArray-single';
        catch
            dataCast = 'single';
        end
    else
        dataCast = 'single';
    end
    inputArgs = {'Smooth', false, 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
                 'DataCast', dataCast, 'PlotSim', false};

    % --- Forward simulation ---
    source_fwd    = struct();
    source_fwd.p0 = p0;
    sensorData    = kspaceFirstOrder3D(kgrid, kmedium, source_fwd, sensor, inputArgs{:});
    sensorData    = smooth(sensorData);
    FS            = 1 / kgrid.dt;

    % --- Pulse model: convolve, frequency response, noise, Wiener deconv ---
    conv_kernel_sigma  = CONFIG.convolution_kernel;
    conv_deconv_lambda = CONFIG.conv_deconv_lambda;

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

    sensorData_conv = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));
    sensorData_resp = gaussianFilter(sensorData_conv, FS, 0.35e6, 100, true);
    noise_amp       = CONFIG.conv_noise_level * max(abs(sensorData_resp(:)));

    fprintf('  [Bundle %s] forward done. Sensor [%d x %d], noise amp %.3e Pa\n', ...
        label, size(sensorData,1), size(sensorData,2), noise_amp);

    % --- Capture the bundle ---
    bodyMask = [];
    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), gridSize_orig)
        bodyMask = double(sct.bodyMask);
    end

    B = struct();
    B.kgrid                = kgrid;
    B.kmedium              = kmedium;
    B.sensor               = sensor;
    B.inputArgs            = inputArgs;
    B.p0                   = p0;
    B.gridSize             = gridSize;
    B.gridSize_orig        = gridSize_orig;
    B.did_pad              = did_pad;
    B.Nx_orig              = Nx_orig;
    B.Ny_orig              = Ny_orig;
    B.Nz_orig              = Nz_orig;
    B.dx = dx; B.dy = dy; B.dz = dz; B.dt = dt;
    B.medium_pad           = medium;
    B.medium_orig          = medium_orig;
    B.doseMask             = doseMask;
    B.num_pulses           = num_pulses;
    B.sensor_info_orig     = sensor_info_orig;
    B.das_opts             = struct();
    B.bodyMask             = bodyMask;
    B.reconstruction_method         = CONFIG.reconstruction_method;
    B.num_time_reversal_iter        = CONFIG.num_time_reversal_iter;
    B.convergence_tol               = CONFIG.convergence_tol;
    B.correction_factor             = CONFIG.correction_factor;
    B.use_pressure_scale_correction = CONFIG.use_pressure_scale_correction;
    B.mask_recon_to_dose_region     = CONFIG.mask_recon_to_dose_region;
    B.sensorData_resp      = double(gather(sensorData_resp));
    B.H_conj               = H_conj;
    B.H_power              = H_power;
    B.conv_deconv_lambda   = conv_deconv_lambda;
    B.noise_amp            = noise_amp;
    B.truth_dose           = doseGrid;
    B.spacing_mm           = spacing_mm;
    B.gantry_angle         = gantry_angle;
    B.pre_expand_dims      = dims_pre_expand;
    B.embed_offset         = embed_offset;
    B.label                = label;
end

function [medium, p0] = pad_medium_p0(medium, p0, Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad)
%PAD_MEDIUM_P0 Zero/water-pad medium + p0 from [Nx Ny Nz] to the padded size.
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
end

function [medium, p0, doseGrid, doseMask, sct, medium_orig, ...
          Nx_orig, Ny_orig, Nz_orig, Nx, Ny, Nz] = ...
        apply_grid_expansion(medium, p0, doseGrid, doseMask, sct, gp, ...
            Nx_orig, Ny_orig, Nz_orig, CONFIG)
%APPLY_GRID_EXPANSION Water-pad the medium/p0/dose to the sensor-expanded grid,
%  then re-run FFT-optimal padding.

    density_unp    = medium.density(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    soundSpeed_unp = medium.sound_speed(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_unp = medium.alpha_coeff(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    else
        alphaCoeff_unp = medium.alpha_coeff;
    end
    gruneisen_unp  = medium.gruneisen(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    p0_unp         = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

    % Caller dim order: dim1=Nx, dim2=Ny, dim3=Nz -> gp.y_*->dim1, x_*->dim2, z_*->dim3.
    Nx_exp = Nx_orig + gp.y_pre + gp.y_post;
    Ny_exp = Ny_orig + gp.x_pre + gp.x_post;
    Nz_exp = Nz_orig + gp.z_pre + gp.z_post;

    density_exp    = ones(Nx_exp, Ny_exp, Nz_exp)  * 1000;
    soundSpeed_exp = ones(Nx_exp, Ny_exp, Nz_exp)  * 1540;
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

    medium_orig = struct( ...
        'density',     density_exp, ...
        'sound_speed', soundSpeed_exp, ...
        'alpha_coeff', alphaCoeff_exp, ...
        'gruneisen',   gruneisen_exp);

    doseGrid_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
    doseGrid_exp(xr, yr, zr) = doseGrid;
    doseGrid = doseGrid_exp;

    doseMask_exp = false(Nx_exp, Ny_exp, Nz_exp);
    doseMask_exp(xr, yr, zr) = doseMask;
    doseMask = doseMask_exp;

    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), [Nx_orig, Ny_orig, Nz_orig])
        body_exp = false(Nx_exp, Ny_exp, Nz_exp);
        body_exp(xr, yr, zr) = sct.bodyMask;
        sct.bodyMask = body_exp;
    end
    if isfield(sct, 'couchMask') && ~isempty(sct.couchMask) && ...
            isequal(size(sct.couchMask), [Nx_orig, Ny_orig, Nz_orig])
        couch_exp = false(Nx_exp, Ny_exp, Nz_exp);
        couch_exp(xr, yr, zr) = sct.couchMask;
        sct.couchMask = couch_exp;
    end

    Nx_orig = Nx_exp; Ny_orig = Ny_exp; Nz_orig = Nz_exp;

    if CONFIG.use_grid_padding
        Nx_pad2 = find_optimal_kwave_size(Nx_orig, CONFIG.pml_size);
        Ny_pad2 = find_optimal_kwave_size(Ny_orig, CONFIG.pml_size);
        Nz_pad2 = find_optimal_kwave_size(Nz_orig, CONFIG.pml_size);
    else
        Nx_pad2 = Nx_orig; Ny_pad2 = Ny_orig; Nz_pad2 = Nz_orig;
    end

    if ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], [Nx_orig, Ny_orig, Nz_orig])
        [medium, p0] = pad_medium_p0(medium, p0, Nx_orig, Ny_orig, Nz_orig, ...
            Nx_pad2, Ny_pad2, Nz_pad2);
    end

    Nx = Nx_pad2; Ny = Ny_pad2; Nz = Nz_pad2;
end

function sd = redraw_noisy_deconv(B, seed)
%REDRAW_NOISY_DECONV  Re-apply a fresh electronic-noise realization to the
%  cached pre-noise sensor response, then Wiener-deconvolve the pulse kernel.
%  For the noise-only null, B.sensorData_resp has been zeroed by the caller, so
%  this returns a purely-noise, deconvolved sensor trace.
    rng(seed);
    noisy = B.sensorData_resp + B.noise_amp * randn(size(B.sensorData_resp));
    deconv = real(ifft( ...
        fft(noisy, [], 2) .* B.H_conj ./ (B.H_power + B.conv_deconv_lambda), [], 2));
    sd = single(deconv);
end

function recon_dose = reconstruct_recon_dose(B, sensorData_measured)
%RECONSTRUCT_RECON_DOSE  Reconstruct a dose from a (noisy, deconvolved) sensor
%  trace using the cached forward bundle B. parfor/GPU-safe (no plotting).

    gridSize = B.gridSize;
    Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);
    inputArgs = B.inputArgs;
    N_iter = B.num_time_reversal_iter;

    switch lower(B.reconstruction_method)
    case 'tr'
        reconPressure      = zeros(gridSize);
        reconPressure_prev = zeros(gridSize);
        sd      = sensorData_measured;
        sd_meas = sensorData_measured;
        for it = 1:N_iter
            source_tr        = struct();
            source_tr.p_mask = B.sensor.mask;
            source_tr.p      = fliplr(sd);
            source_tr.p_mode = 'dirichlet';

            sensor_tr        = struct();
            sensor_tr.mask   = ones(Nx, Ny, Nz);
            sensor_tr.record = {'p_final'};

            pr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_tr, sensor_tr, inputArgs{:});
            if isstruct(pr) && isfield(pr, 'p_final')
                reconPressure = reshape(pr.p_final, [Nx, Ny, Nz]);
            else
                reconPressure = reshape(pr, [Nx, Ny, Nz]);
            end
            reconPressure = max(reconPressure, 0);

            if it > 1
                np = norm(reconPressure_prev(:));
                if np > 0
                    rc = norm(reconPressure(:) - reconPressure_prev(:)) / np;
                else
                    rc = Inf;
                end
                if rc < B.convergence_tol
                    break;
                end
            end
            reconPressure_prev = reconPressure;

            if it < N_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sdr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_resid, B.sensor, inputArgs{:});
                sd  = sd + (sd_meas - sdr);
            end
        end
        reconPressure = gather(reconPressure) * B.correction_factor;

    case 'das'
        reconPressure = das_reconstruct(sensorData_measured, B.sensor, B.sensor_info_orig, ...
            B.medium_pad, Nx, Ny, Nz, B.dx, B.dy, B.dz, B.dt, B.das_opts);
        reconPressure = reconPressure * B.correction_factor;

    case 'hybrid'
        reconPressure = das_reconstruct(sensorData_measured, B.sensor, B.sensor_info_orig, ...
            B.medium_pad, Nx, Ny, Nz, B.dx, B.dy, B.dz, B.dt, B.das_opts);
        reconPressure_prev = reconPressure;
        sd      = sensorData_measured;
        sd_meas = sensorData_measured;
        if N_iter > 1
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sdr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_resid, B.sensor, inputArgs{:});
            sd  = sd + (sd_meas - sdr);
        end
        for it = 2:N_iter
            source_tr        = struct();
            source_tr.p_mask = B.sensor.mask;
            source_tr.p      = fliplr(sd);
            source_tr.p_mode = 'dirichlet';

            sensor_tr        = struct();
            sensor_tr.mask   = ones(Nx, Ny, Nz);
            sensor_tr.record = {'p_final'};

            pr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_tr, sensor_tr, inputArgs{:});
            if isstruct(pr) && isfield(pr, 'p_final')
                reconPressure = reshape(pr.p_final, [Nx, Ny, Nz]);
            else
                reconPressure = reshape(pr, [Nx, Ny, Nz]);
            end
            reconPressure = max(reconPressure, 0);

            np = norm(reconPressure_prev(:));
            if np > 0
                rc = norm(reconPressure(:) - reconPressure_prev(:)) / np;
            else
                rc = Inf;
            end
            if rc < B.convergence_tol
                break;
            end
            reconPressure_prev = reconPressure;

            if it < N_iter
                source_resid    = struct();
                source_resid.p0 = reconPressure;
                sdr = kspaceFirstOrder3D(B.kgrid, B.kmedium, source_resid, B.sensor, inputArgs{:});
                sd  = sd + (sd_meas - sdr);
            end
        end
        reconPressure = gather(reconPressure) * B.correction_factor;

    otherwise
        error('reconstruct_recon_dose:UnknownMethod', ...
            'Unknown reconstruction_method: "%s"', B.reconstruction_method);
    end

    % --- Crop to original size ---
    if B.did_pad
        reconPressure = reconPressure(1:B.Nx_orig, 1:B.Ny_orig, 1:B.Nz_orig);
    end

    % --- Pressure scale correction (optional) ---
    if B.use_pressure_scale_correction
        p0_max_orig = max(B.p0(1:B.Nx_orig, 1:B.Ny_orig, 1:B.Nz_orig), [], 'all');
        recon_max   = max(reconPressure(:));
        if recon_max > 0
            reconPressure = reconPressure * (p0_max_orig / recon_max);
        end
    end

    % --- Pressure -> dose ---
    conversionFactor = B.medium_orig.gruneisen .* B.medium_orig.density;
    conversionFactor(conversionFactor == 0) = 1;
    reconDosePerPulse = reconPressure ./ conversionFactor;

    cropSize = B.gridSize_orig;
    if ~isempty(B.bodyMask) && isequal(size(B.bodyMask), cropSize)
        body_mask_plot = B.bodyMask;
    else
        body_mask_plot = ones(cropSize);
    end

    if ~B.mask_recon_to_dose_region
        recon_dose = reconDosePerPulse * B.num_pulses .* body_mask_plot;
    else
        recon_dose = reconDosePerPulse * B.num_pulses .* double(B.doseMask) .* body_mask_plot;
    end

    % Mask air pockets (identical to run_single_field_simulation): air-classified
    % voxels (density ~1.2 kg/m^3) carry no real PA dose and would otherwise
    % reconstruct as hotspots. Use the conversion medium's density.
    recon_dose(B.medium_orig.density < 100) = 0;
end

function medium = build_medium_with_bath(sct, config)
%BUILD_MEDIUM_WITH_BATH  create_acoustic_medium (the single, canonical medium
%  builder -- the only one that applies the in-body air row) + force the coupling
%  bath (outside body / couch) to the uniform medium. Identical to the local of
%  the same name in run_standalone_field, so the noise-only null reconstructs on
%  exactly the medium the real recon uses.
    medium = create_acoustic_medium(sct, config);
    ud = getf(config, 'uniform_density',     1000);
    uc = getf(config, 'uniform_sound_speed', 1540);
    ua = getf(config, 'uniform_alpha_coeff', 0);
    ug = getf(config, 'uniform_gruneisen',   1.0);
    if isfield(sct, 'bodyMask')
        outside = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask') && ~isempty(sct.couchMask)
            outside = outside | logical(sct.couchMask);
        end
        medium.density(outside)     = ud;
        medium.sound_speed(outside) = uc;
        if numel(medium.alpha_coeff) > 1
            medium.alpha_coeff(outside) = ua;
        end
        medium.gruneisen(outside)   = ug;
    end
end
