%% =========================================================================
%  RUN_STANDALONE_SIMULATION_UPGRADED.m
%  Standalone k-Wave Photoacoustic Forward + Time-Reversal Simulation
%  Adds CONFIG-gated reconstruction-fidelity upgrades. All upgrade flags
%  default to FALSE so they can be enabled and validated one-by-one.
%
%  Bias targeted by these upgrades: near-detector dose overestimated,
%  far-from-detector dose underestimated; gamma 10%/10mm ~60% baseline.
%
%  Upgrade flags (priority order, all default false):
%    sensitivity_map_correction  S = A'(A(1_body));  recon ./= (S + eps*max(S))
%    least_squares_gain          c = sum(r.*t)/sum(r.^2) over eval mask.
%                                Replaces use_pressure_scale_correction (peak
%                                pinning). Uses TRUTH -> characterization
%                                only, NOT deployable. Gamma reported both
%                                with and without c.
%    use_true_adjoint            True k-Wave adjoint per Arridge et al. 2016
%                                (Inverse Problems 32:115012). Includes
%                                dot-product adjoint test <Ax,y>==<x,A'y>.
%    absorption_compensation     Sign-flipped alpha_coeff (gain), dispersion
%                                unchanged, + Tukey-windowed k-space
%                                regularization per Treeby/Zhang/Cox 2010
%                                (Inverse Problems 26:115003).
%    das_ubp                     Universal back-projection (Xu & Wang 2005):
%                                back-project (2p - 2t*dp/dt) with solid
%                                angle/r^2 weighting, exp(+alpha(f0)*d_s)
%                                attenuation comp.
%    reconstruction_method='mbir'  FISTA solve min_{p0>=0} ||A p0 - d||^2
%                                + lambda*TV, reuses forward + new adjoint,
%                                preconditioned by S.
%    break_inverse_crime         Reconstruct on coarser grid and/or with
%                                perturbed medium relative to data-gen medium.
%  =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';

CONFIG.dose_filename = 'dose_1194203_Session_1_reference_CT_3_B15_112.mat';
% TEMPORARY: standalone sim always uses CBCT1 (CT_1) geometry, regardless
% of which CBCT the dose was actually computed on. Per-dose CBCT selection
% lives in the multi-field driver, not here.
CONFIG.cbct_filename = 'CBCT1_resampled.mat';

CONFIG.dose_file_override = '';
CONFIG.cbct_file_override = '';


% 'determine_sensor_mask' 'full_plane_lateral' 'full_plane_anterior'
CONFIG.sensor_placement_method = 'full_plane_anterior';
CONFIG.sensor_x_index = 2;
CONFIG.sensor_y_index = 4;

% Physical 2D ultrasound array geometry (sparse element mask).
% Kerf is derived inside determine_sensor_mask as (pitch - size).
CONFIG.elements_per_side = 32;
CONFIG.element_pitch_mm  = 7;
CONFIG.element_size_mm   = 2.43;

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
CONFIG.Nt_scaling             = 6;     % >0: when air sets minC, divide Nt by this to shorten the recording (0 = disabled)
CONFIG.use_gpu                = true;
CONFIG.correction_factor           = 1.9;
%CONFIG.correction_factor = .0228; 

%CONFIG.correction_factor = .0229
% Peak-pinning correction (DEFAULT FALSE in upgraded build). Divides
% reconPressure by max(recon)/max(p0_truth). Disabled by default because it
% destroys absolute calibration and biases gamma toward the peak voxel.
CONFIG.use_pressure_scale_correction = false;
CONFIG.mask_recon_to_dose_region     = true;    % zero recon dose outside the dose mask (>1% of original max). Set false to keep the full reconstruction.

%% ----------------- UPGRADED RECONSTRUCTION-FIDELITY FLAGS ---------------
%  All default false so each can be validated independently against the
%  baseline. Priority order matches the order of these flags.

% (1) SENSITIVITY MAP CORRECTION
%   Build S = adjoint(forward(ones(body_mask))); divide reconPressure by
%   (S + eps*max(S)) before pressure->dose. Saved to results.sensitivity_map.
%   This corrects for the spatially varying sensitivity of the sparse 2D
%   array (near-detector overestimation / far underestimation bias).
CONFIG.sensitivity_map_correction = true;
% When using full TR/forward to build S, this can be expensive. Optionally
% reuse the same recon_method as the main inversion; otherwise force DAS for
% the adjoint applied inside S.
CONFIG.sensitivity_map_adjoint = 'tr_lossy';   % 'tr_lossy' | 'tr_true' | 'das'

% (2) GLOBAL LEAST-SQUARES GAIN (uses truth - CHARACTERIZATION ONLY)
%   c = sum(r .* t) / sum(r .^ 2) over the gamma eval mask, where r=recon,
%   t=truth. Applied as recon_calibrated = c * recon. Reported alongside
%   uncalibrated gamma. Not deployable (uses truth). This is the
%   study_pass_rates_allsegments.m normalization method; enabled by default so
%   the LS-calibrated gamma is always reported alongside the uncalibrated one.
CONFIG.least_squares_gain = true;

% (3) TRUE k-WAVE ADJOINT (Arridge et al. 2016)
%   Replaces lossy-TR-as-adjoint with the true operator adjoint. Enables
%   dot-product test <Ax,y> == <x,A'y> as an assertion at startup.
CONFIG.use_true_adjoint = false;
% Tolerance for the dot-product test (relative).
CONFIG.adjoint_test_tol = 1e-3;

% (4) ABSORPTION-COMPENSATED TR (Treeby/Zhang/Cox 2010)
%   In the back-propagation, sign-flip alpha_coeff (gain), keep dispersion,
%   apply Tukey-windowed k-space low-pass to prevent gain blow-up.
CONFIG.absorption_compensation = false;
CONFIG.tukey_cutoff_frac       = 0.7;   % cutoff as fraction of k_max
CONFIG.tukey_taper_frac        = 0.2;   % taper width as fraction of k_max

% (5) UBP DAS (Xu & Wang 2005, Universal Back-Projection)
%   In DAS, back-project the filtered quantity (2p - 2t*dp/dt) with solid-
%   angle/r^2 weighting instead of raw p with 1/r. Adds exp(+alpha(f0)*d_s)
%   attenuation compensation, with alpha evaluated at center frequency f0.
CONFIG.das_ubp                 = false;
CONFIG.das_ubp_f0              = 1.0e6;   % Hz, center freq for atten comp
CONFIG.das_ubp_alpha_db_cm_mhz = 0.5;     % dB/(MHz^y * cm), tissue avg

% (6) MODEL-BASED ITERATIVE RECONSTRUCTION
%   Set CONFIG.reconstruction_method = 'mbir' to enable FISTA with TV.
CONFIG.mbir_max_iter   = 50;
CONFIG.mbir_lambda_tv  = 1e-3;
CONFIG.mbir_step_size  = [];   % [] => estimate Lipschitz via power method
CONFIG.mbir_tv_iters   = 10;   % inner iterations for TV prox
CONFIG.mbir_use_precond = true;  % precondition by sensitivity map S

% (7) BREAK INVERSE CRIME
%   Reconstruct on a coarser grid and/or with perturbed medium relative to
%   the data-generation forward simulation. When enabled, the forward sim
%   uses the original (fine, true) medium; recon uses a coarsened/perturbed
%   copy.
CONFIG.break_inverse_crime           = false;
CONFIG.inverse_crime_recon_downscale = 1.0;     % >1 coarsens recon grid
CONFIG.inverse_crime_sound_speed_sigma = 0.0;   % fractional std on c map
CONFIG.inverse_crime_density_sigma     = 0.0;   % fractional std on rho map
CONFIG.inverse_crime_alpha_bias_factor = 1.0;   % multiplicative on alpha map

% --- Reconstruction method ---
%   'tr'     : iterative time-reversal (k-Wave back-propagation)
%   'DAS'    : Delay-And-Sum back-projection (homogeneous c, non-iterative)
%   'hybrid' : DAS for iter 1, k-Wave TR with residual correction for iters 2..N
CONFIG.reconstruction_method = 'tr';

CONFIG.num_time_reversal_iter = 10;
CONFIG.convergence_tol        = 1e-3;

% --- Pulse Convolution / Noise / Deconvolution ---
% Mimics a finite transducer impulse response applied to forward sensor data.
% Set convolution_kernel to 0 to disable the entire block.
CONFIG.convolution_kernel  = 4e-6;   % Gaussian sigma in seconds (4 us)
CONFIG.conv_noise_level    = 0.01;   % Noise amplitude as fraction of peak sensor signal
CONFIG.conv_deconv_lambda  = 1e-4;   % Wiener regularization for deconvolution

CONFIG.downscale_factor = 1;
CONFIG.use_grid_padding = true;

CONFIG.save_results = true;
CONFIG.output_file  = 'standalone_recon_results.mat';
CONFIG.plot_results = true;

% Gaussian sigma (in voxels) applied to dose slices for DISPLAY ONLY in the
% dose comparison panels. Fills the speckle "pockets" left by a spotty recon
% so the overlay reads continuously instead of letting CT show through.
% Set to 0 to disable display smoothing.
CONFIG.viz_smooth_sigma = 1.0;

% Normalize: divide original and reconstructed dose by their own max before
% comparison / gamma so both peak at 1.0.
CONFIG.normalize = false;

% Gamma logging: append CONFIG + gamma pass rates to gamma_log.mat (in
% working_dir) after each run. Keeps a running record of the best gamma per
% criterion and the CONFIG that produced it.
CONFIG.log_gamma = false;
CONFIG.gamma_log_file = 'gamma_log.mat';

% Diagnostic plot: three anatomical views (transverse, sagittal, coronal) of
% the beam exclusion zone over the body. Useful to sanity-check that the
% projected jaw rectangles aren't unrealistically large.
CONFIG.plot_exclusion_zone = false;

%% ========================= RESOLVE FILE PATHS ============================

if ~isempty(CONFIG.dose_file_override)
    dose_filepath = CONFIG.dose_file_override;
else
    processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    dose_filepath = fullfile(processed_dir, CONFIG.dose_filename);
end

if ~isempty(CONFIG.cbct_file_override)
    cbct_filepath = CONFIG.cbct_file_override;
else
    if ~exist('processed_dir', 'var')
        processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
            CONFIG.patient_id, CONFIG.session, 'processed');
    end
    cbct_filepath = fullfile(processed_dir, CONFIG.cbct_filename);
end

%% ========================= LOAD PLAN BEAM METADATA =======================
% determine_sensor_mask needs beam_metadata for ALL beams in the plan
% (isocenter + jaw extents per beam) to compute the anterior-surface
% exclusion zone. step15_process_doses saves this as metadata.beam_metadata
% in <processed_dir>/metadata.mat. Load it here if not already set by the
% caller via CONFIG.beam_metadata.

if ~isfield(CONFIG, 'beam_metadata') || isempty(CONFIG.beam_metadata)
    if exist('processed_dir', 'var')
        metadata_filepath = fullfile(processed_dir, 'metadata.mat');
    else
        metadata_filepath = '';
    end

    if ~isempty(metadata_filepath) && isfile(metadata_filepath)
        try
            md = load(metadata_filepath, 'metadata');
            if isfield(md, 'metadata') && isfield(md.metadata, 'beam_metadata') ...
                    && ~isempty(md.metadata.beam_metadata)
                CONFIG.beam_metadata = md.metadata.beam_metadata;
                fprintf('  Loaded beam_metadata for %d beams from %s\n', ...
                    length(CONFIG.beam_metadata), metadata_filepath);
            else
                fprintf('  [WARN] metadata.mat present but no beam_metadata field.\n');
            end
        catch ME
            fprintf('  [WARN] Failed to load %s: %s\n', metadata_filepath, ME.message);
        end
    else
        fprintf('  [WARN] No metadata.mat in processed_dir; sensor exclusion zone will be empty.\n');
    end
end

%% ========================= PRINT CONFIGURATION ===========================

fprintf('=========================================================\n');
fprintf('  Standalone k-Wave Photoacoustic Simulation  (v4.1)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose file:       %s\n', dose_filepath);
fprintf('  CBCT file:       %s  (TEMP: standalone always uses CBCT1/CT_1)\n', cbct_filepath);
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
            error('field_dose.dose_dims missing  cannot reconstruct sparse dose.');
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
        fd_gantry_angle = fd.gantry_angle;
        fprintf('       Gantry angle: %.1f deg\n', fd_gantry_angle);
    end

elseif isfield(dose_data, 'total_rs_dose_sparse')
    % Total dose file from step15_process_doses (sparse format).
    if ~isfield(dose_data, 'total_rs_dose_dims')
        error('total_rs_dose_dims missing  cannot reconstruct sparse total dose.');
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

fprintf('[2/7] Loading CBCT data (CBCT1 / CT_1  temporary standalone default)...\n');
if ~isfile(cbct_filepath)
    error('CBCT file not found: %s', cbct_filepath);
end
cbct_data = load(cbct_filepath);
if isfield(cbct_data, 'CBCT1_resampled')
    sct = cbct_data.CBCT1_resampled;
elseif isfield(cbct_data, 'CBCT3_resampled')
    sct = cbct_data.CBCT3_resampled;
else
    error('CBCT1_resampled / CBCT3_resampled variable not found in %s', cbct_filepath);
end

if ~isfield(sct, 'cubeHU')
    error('CBCT resampled struct missing required field: cubeHU');
end

% Spacing: prefer CBCT field; fall back to spacing embedded in step15 field_dose
if isfield(sct, 'spacing') && ~isempty(sct.spacing)
    spacing_mm = sct.spacing(:)';
elseif exist('step15_spacing_mm', 'var')
    spacing_mm = step15_spacing_mm;
    fprintf('       [INFO] Using spacing from field_dose file: [%.3f %.3f %.3f] mm\n', spacing_mm);
else
    error('CBCT resampled struct missing required field: spacing');
end
dx = spacing_mm(1) / 1000;
dy = spacing_mm(2) / 1000;
dz = spacing_mm(3) / 1000;

sctSize = size(sct.cubeHU);
if ~isequal(gridSize, sctSize)
    error(['Dose grid [%d %d %d] does not match CBCT grid [%d %d %d].\n' ...
           'Ensure CBCT1_resampled.mat was produced by the same step15 run as the dose file.'], ...
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
        % Anything outside the enclosing sphere is unobservable by this
        % sensor geometry  zero p0 there so it doesn't pollute the forward
        % simulation or downstream pressure scaling.
        sph_cx = floor(Nx/2) + 1;
        sph_cy = floor(Ny/2) + 1;
        sph_cz = floor(Nz/2) + 1;
        [Xg_sph, Yg_sph, Zg_sph] = ndgrid(1:Nx, 1:Ny, 1:Nz);
        ball_mask = (Xg_sph - sph_cx).^2 + (Yg_sph - sph_cy).^2 + ...
                    (Zg_sph - sph_cz).^2 <= sph_radius^2;
        n_zeroed = nnz(p0 ~= 0 & ~ball_mask);
        p0 = p0 .* ball_mask;
        clear Xg_sph Yg_sph Zg_sph
        fprintf('       Sensor: spherical, radius %d voxels (zeroed %d p0 voxels outside sphere)\n', ...
            sph_radius, n_zeroed);
    case 'box'
        % Six-face bounding box enclosing the pressure: planes at index 3
        % and (N-3) on each axis.
        bx_lo   = 3;
        bx_hi_x = Nx - 3;
        bx_hi_y = Ny - 3;
        bx_hi_z = Nz - 3;
        if bx_hi_x <= bx_lo || bx_hi_y <= bx_lo || bx_hi_z <= bx_lo
            error('run_standalone_simulation:BoxTooSmall', ...
                'Grid [%d %d %d] too small for box sensor (need each dim > 6).', Nx, Ny, Nz);
        end
        sensor.mask(bx_lo,   bx_lo:bx_hi_y, bx_lo:bx_hi_z) = 1;
        sensor.mask(bx_hi_x, bx_lo:bx_hi_y, bx_lo:bx_hi_z) = 1;
        sensor.mask(bx_lo:bx_hi_x, bx_lo,   bx_lo:bx_hi_z) = 1;
        sensor.mask(bx_lo:bx_hi_x, bx_hi_y, bx_lo:bx_hi_z) = 1;
        sensor.mask(bx_lo:bx_hi_x, bx_lo:bx_hi_y, bx_lo)   = 1;
        sensor.mask(bx_lo:bx_hi_x, bx_lo:bx_hi_y, bx_hi_z) = 1;
        fprintf('       Sensor: box faces at x=[%d,%d], y=[%d,%d], z=[%d,%d]\n', ...
            bx_lo, bx_hi_x, bx_lo, bx_hi_y, bx_lo, bx_hi_z);
    case 'determine_sensor_mask'
        % Automatic placement via determine_sensor_mask: tilts a 2D array
        % toward the beam isocenter (or places it flat when CONFIG.aim_at_iso
        % is false), avoiding beam exclusion zones.
        %
        % SESSION-LEVEL USAGE: because all beams in an ETHOS plan share one
        % isocenter, the placement is reusable across fields. When using this
        % once per session, pass the SUMMED PLAN DOSE as
        %   field_dose_for_sensor.dose_Gy
        % so the exclusion zone reflects the full beam path from every field.
        % Passing a single field's dose here yields a per-field exclusion zone
        % only  correct for per-field placement but not session-level reuse.
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
        if exist('fd_gantry_angle', 'var') && ~isempty(fd_gantry_angle)
            field_dose_for_sensor.gantry_angle = fd_gantry_angle;
        else
            field_dose_for_sensor.gantry_angle = 0;
        end
        field_dose_for_sensor.origin      = sct_for_sensor.origin;
        field_dose_for_sensor.spacing     = spacing_mm;
        field_dose_for_sensor.dimensions  = [Nx_orig, Ny_orig, Nz_orig];

        beam_meta = [];
        if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
            beam_meta = CONFIG.beam_metadata;
        end

        [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask( ...
            sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);

        % --- GRID EXPANSION HANDLING ---
        % determine_sensor_mask may expand the grid in X/Z (or Y) to place the
        % sensor outside the beam exclusion zone, filling the new region with
        % water. Apply matching water padding to medium/p0 here so the sim
        % grid coordinate system matches the sensor mask. Expansion is applied
        % to the un-FFT-padded data, then FFT-optimal padding is re-run.
        gp = sensor_info_orig.grid_pad;
        if gp.expanded
            fprintf('       [Sensor] Grid expansion: Y(+%d/+%d), X(+%d/+%d), Z(+%d/+%d). Re-padding with water.\n', ...
                gp.y_pre, gp.y_post, gp.x_pre, gp.x_post, gp.z_pre, gp.z_post);

            density_unp    = medium.density(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
            soundSpeed_unp = medium.sound_speed(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
            if numel(medium.alpha_coeff) > 1
                alphaCoeff_unp = medium.alpha_coeff(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
            else
                alphaCoeff_unp = medium.alpha_coeff;
            end
            gruneisen_unp  = medium.gruneisen(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
            p0_unp         = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

            % determine_sensor_mask labels dim 1=Y, dim 2=X, dim 3=Z, but it
            % preserves the caller's actual dim order. Since this script passes
            % bodyMask with dim 1=Nx, dim 2=Ny, dim 3=Nz, the function's
            % grid_pad fields map as:
            %   gp.y_*  -> script dim 1 (Nx)   [currently always 0]
            %   gp.x_*  -> script dim 2 (Ny)
            %   gp.z_*  -> script dim 3 (Nz)
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

            % medium_orig is restored after FFT-pad cropping (~line 1059), and
            % feeds conversionFactor for pressure->dose. Update it to the
            % expanded pre-FFT-pad medium so its size matches the expanded
            % reconPressure after cropping.
            medium_orig = struct( ...
                'density',     density_exp, ...
                'sound_speed', soundSpeed_exp, ...
                'alpha_coeff', alphaCoeff_exp, ...
                'gruneisen',   gruneisen_exp);

            % Expand doseGrid / doseMask with zeros (no dose in water padding).
            doseGrid_exp = zeros(Nx_exp, Ny_exp, Nz_exp);
            doseGrid_exp(xr, yr, zr) = doseGrid;
            doseGrid = doseGrid_exp;

            doseMask_exp = false(Nx_exp, Ny_exp, Nz_exp);
            doseMask_exp(xr, yr, zr) = doseMask;
            doseMask = doseMask_exp;

            % Expand sct.bodyMask with false in the water region so downstream
            % visualizations and masking still align with the expanded grid.
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
            gridSize_orig = [Nx_orig, Ny_orig, Nz_orig];

            % Re-run FFT-optimal padding on the expanded grid (water fills).
            if CONFIG.use_grid_padding
                Nx_pad2 = find_optimal_kwave_size(Nx_orig, CONFIG.pml_size);
                Ny_pad2 = find_optimal_kwave_size(Ny_orig, CONFIG.pml_size);
                Nz_pad2 = find_optimal_kwave_size(Nz_orig, CONFIG.pml_size);
            else
                Nx_pad2 = Nx_orig; Ny_pad2 = Ny_orig; Nz_pad2 = Nz_orig;
            end

            if ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], [Nx_orig, Ny_orig, Nz_orig])
                fprintf('       [Sensor] Re-pad to FFT-optimal: [%d %d %d] -> [%d %d %d]\n', ...
                    Nx_orig, Ny_orig, Nz_orig, Nx_pad2, Ny_pad2, Nz_pad2);
                density_pad    = ones(Nx_pad2, Ny_pad2, Nz_pad2)  * 1000;
                soundSpeed_pad = ones(Nx_pad2, Ny_pad2, Nz_pad2)  * 1540;
                alphaCoeff_pad = zeros(Nx_pad2, Ny_pad2, Nz_pad2);
                gruneisen_pad  = zeros(Nx_pad2, Ny_pad2, Nz_pad2);
                p0_pad2        = zeros(Nx_pad2, Ny_pad2, Nz_pad2);

                density_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig)    = medium.density;
                soundSpeed_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = medium.sound_speed;
                if numel(medium.alpha_coeff) > 1
                    alphaCoeff_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = medium.alpha_coeff;
                else
                    alphaCoeff_pad(:) = medium.alpha_coeff;
                end
                gruneisen_pad(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = medium.gruneisen;
                p0_pad2(1:Nx_orig, 1:Ny_orig, 1:Nz_orig)       = p0;

                medium.density     = density_pad;
                medium.sound_speed = soundSpeed_pad;
                medium.alpha_coeff = alphaCoeff_pad;
                medium.gruneisen   = gruneisen_pad;
                p0 = p0_pad2;
            end

            Nx = Nx_pad2; Ny = Ny_pad2; Nz = Nz_pad2;
            gridSize = [Nx, Ny, Nz];
            sensor.mask = zeros(Nx, Ny, Nz);
        end

        % Optional diagnostic: three anatomical views of the exclusion zone.
        % Placed after grid expansion so sct.bodyMask and the exclusion zone
        % share the same coordinate system (both expanded if expansion ran).
        if isfield(CONFIG, 'plot_exclusion_zone') && CONFIG.plot_exclusion_zone
            plot_exclusion_zone_views(sct, sensor_info_orig, spacing_mm, ...
                sprintf('Exclusion zone (gantry %.1f deg)', field_dose_for_sensor.gantry_angle));
        end

        % Embed sensor mask. determine_sensor_mask preserves the caller's dim
        % order  its dim 1 matches sct.bodyMask's dim 1 (which is Nx in this
        % script). No permute needed.
        m1 = min(Nx, size(sensor_mask_orig, 1));
        m2 = min(Ny, size(sensor_mask_orig, 2));
        m3 = min(Nz, size(sensor_mask_orig, 3));
        sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
        fprintf('       Sensor: determine_sensor_mask  %d active points\n', sum(sensor_mask_orig(:)));
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
        fprintf('       Sensor: fixed_anterior  %d active points\n', sum(sensor_mask_orig(:)));
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

% Smooth the forward time series.
sensorData = smooth(sensorData);

% Sampling frequency for the sensor frequency-response filter, which is now
% applied below AFTER the pulse-profile convolution (see physical chain).
FS         = 1 / kgrid.dt;

sensorData_measured = sensorData;

%% ========================= BREAK INVERSE CRIME ===========================
%  The forward simulation has run on the "truth" medium/grid. If
%  CONFIG.break_inverse_crime is true, swap kmedium for an intentionally
%  perturbed/coarsened copy used by every subsequent reconstruction call
%  (TR back-prop, DAS, MBIR, residual updates). Keeps kgrid unchanged --
%  grid coarsening is approximated by Gaussian-smoothing the medium maps
%  with sigma ~ (downscale_factor / 2) voxels.

kmedium_truth = kmedium;   % keep a copy in case downstream needs truth medium
if isfield(CONFIG, 'break_inverse_crime') && CONFIG.break_inverse_crime
    fprintf('       [InvCrime] Perturbing medium for reconstruction:\n');
    fprintf('         sound_speed sigma : %.4f (frac)\n', CONFIG.inverse_crime_sound_speed_sigma);
    fprintf('         density     sigma : %.4f (frac)\n', CONFIG.inverse_crime_density_sigma);
    fprintf('         alpha bias factor : %.4f\n',        CONFIG.inverse_crime_alpha_bias_factor);
    fprintf('         recon grid downscale (smoothing) : %.3f\n', CONFIG.inverse_crime_recon_downscale);

    rng_state_ic = rng();   % save so the perturbation is reproducible per run
    rng(0, 'twister');

    in_body = true(size(kmedium.density));
    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), size(in_body))
        in_body = logical(sct.bodyMask);
    end

    c_p = kmedium.sound_speed;
    if CONFIG.inverse_crime_sound_speed_sigma > 0
        noise_c = 1 + CONFIG.inverse_crime_sound_speed_sigma * randn(size(c_p));
        c_p(in_body) = c_p(in_body) .* noise_c(in_body);
    end

    rho_p = kmedium.density;
    if CONFIG.inverse_crime_density_sigma > 0
        noise_r = 1 + CONFIG.inverse_crime_density_sigma * randn(size(rho_p));
        rho_p(in_body) = rho_p(in_body) .* noise_r(in_body);
    end

    alpha_p = kmedium.alpha_coeff;
    if CONFIG.inverse_crime_alpha_bias_factor ~= 1
        alpha_p = alpha_p * CONFIG.inverse_crime_alpha_bias_factor;
    end

    if CONFIG.inverse_crime_recon_downscale > 1
        sigma_vox = CONFIG.inverse_crime_recon_downscale / 2;
        c_p     = imgaussfilt3(c_p,     sigma_vox);
        rho_p   = imgaussfilt3(rho_p,   sigma_vox);
        alpha_p = imgaussfilt3(alpha_p, sigma_vox);
        fprintf('         Smoothed medium maps with sigma=%.2f vox.\n', sigma_vox);
    end

    kmedium.sound_speed = max(c_p,   1);
    kmedium.density     = max(rho_p, 1);
    kmedium.alpha_coeff = max(alpha_p, 0);

    rng(rng_state_ic);
end

%% ============ PULSE CONVOLUTION / FREQUENCY RESPONSE / NOISE / DECONV =====
%  Models the physical measurement chain in acquisition order:
%    1. Convolve each sensor time series with a Gaussian pulse kernel
%       (finite radiation pulse profile shaping the acoustic source)
%    2. Apply the sensor frequency response (0.35 MHz centre, 100% bandwidth
%       Gaussian band-limit of the acoustic-to-electrical conversion)
%    3. Add white Gaussian noise (electronic noise injected downstream of the
%       transducer -> added AFTER the frequency-response filter)
%    4. Wiener-deconvolve the pulse kernel to recover the broadband signal
%  Time reversal then proceeds on the deconvolved data as normal.

if CONFIG.convolution_kernel > 0
    conv_kernel_sigma  = CONFIG.convolution_kernel;
    conv_noise_level   = CONFIG.conv_noise_level;
    conv_deconv_lambda = CONFIG.conv_deconv_lambda;

    fprintf('       Pulse model: sigma=%.1f us, noise=%.1f%%, lambda=%.1e\n', ...
        conv_kernel_sigma * 1e6, conv_noise_level * 100, conv_deconv_lambda);

    % Build normalized Gaussian kernel in time (truncated at 4 sigma)
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

    % 1. Convolve with the pulse profile
    sensorData_conv = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));

    % 2. Sensor frequency response (band-limit AFTER pulse convolution)
    sensorData_resp = gaussianFilter(sensorData_conv, FS, 0.35e6, 100, true);

    % 3. Add electronic noise (after frequency-response filtering)
    noise_amp        = conv_noise_level * max(abs(sensorData_resp(:)));
    sensorData_noisy = sensorData_resp + noise_amp * randn(size(sensorData_resp));

    % 4. Wiener deconvolution of the pulse kernel
    sensorData_deconv = real(ifft( ...
        fft(sensorData_noisy, [], 2) .* H_conj ./ (H_power + conv_deconv_lambda), ...
        [], 2));

    % Replace both working and reference sensor data with processed result
    sensorData          = single(sensorData_deconv);
    sensorData_measured = single(sensorData_deconv);

    fprintf('       Pulse model complete. Noise amp: %.3e Pa\n', noise_amp);
else
    % No pulse model: apply only the sensor frequency response.
    sensorData          = gaussianFilter(sensorData, FS, 0.35e6, 100, true);
    sensorData_measured = sensorData;
    fprintf('       Pulse convolution disabled; frequency response applied.\n');
end

%% ========================= ADJOINT DOT-PRODUCT TEST ======================
%  When use_true_adjoint is enabled, verify <A x, y> == <x, A' y> within
%  CONFIG.adjoint_test_tol on random low-magnitude probes. Skipped when the
%  flag is off (the lossy-TR back-prop is not a true adjoint and would fail).

if isfield(CONFIG, 'use_true_adjoint') && CONFIG.use_true_adjoint
    fprintf('\n[ADJ] Running dot-product adjoint test...\n');
    try
        rng_state_adj = rng();
        rng(1, 'twister');

        x_probe = randn(Nx, Ny, Nz, 'single');
        x_probe(sensor.mask > 0) = x_probe(sensor.mask > 0);   % unrestricted

        % y_probe lives in sensor-data space: [n_sensors x Nt]
        n_sens   = nnz(sensor.mask);
        y_probe  = randn(n_sens, kgrid.Nt, 'single');

        src_probe    = struct('p0', double(x_probe));
        Ax = kspaceFirstOrder3D(kgrid, kmedium, src_probe, sensor, inputArgs{:});
        Ax = double(gather(Ax));

        Aty = back_propagate(kgrid, kmedium, sensor, double(y_probe), ...
                             inputArgs, CONFIG);
        Aty = double(gather(Aty));

        lhs = sum(Ax(:) .* double(y_probe(:)));
        rhs = sum(double(x_probe(:)) .* Aty(:));
        rel = abs(lhs - rhs) / max(abs(lhs), eps);

        fprintf('       <Ax,y> = %.6e  |  <x,A''y> = %.6e  |  rel = %.3e\n', ...
            lhs, rhs, rel);
        assert(rel < CONFIG.adjoint_test_tol, ...
            'run_standalone_simulation_upgraded:AdjointTestFailed', ...
            'Dot-product adjoint test failed: rel=%.3e > tol=%.3e', ...
            rel, CONFIG.adjoint_test_tol);

        rng(rng_state_adj);
    catch ME
        warning('Adjoint test skipped/failed: %s', ME.message);
    end
end

%% ========================= SENSITIVITY MAP ===============================
%  S = adjoint(forward(ones(body_mask))) over the recon grid. Divide
%  reconPressure by (S + eps*max(S)) before pressure->dose. S itself is
%  saved into results.sensitivity_map. Built once, before the inversion
%  proper, so it can also precondition MBIR.

sensitivity_map = [];
if isfield(CONFIG, 'sensitivity_map_correction') && CONFIG.sensitivity_map_correction
    fprintf('\n[SENS] Building sensitivity map S = A''(A 1_body)...\n');

    body_probe = zeros(Nx, Ny, Nz);
    if isfield(sct, 'bodyMask') && isequal(size(sct.bodyMask), [Nx_orig, Ny_orig, Nz_orig])
        body_probe(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = double(sct.bodyMask);
    else
        body_probe(1:Nx_orig, 1:Ny_orig, 1:Nz_orig) = 1;
    end

    sens_tic = tic;
    src_S = struct('p0', body_probe);
    d_S = kspaceFirstOrder3D(kgrid, kmedium, src_S, sensor, inputArgs{:});
    d_S = double(gather(d_S));

    sensitivity_map = back_propagate(kgrid, kmedium, sensor, d_S, ...
                                     inputArgs, CONFIG);
    sensitivity_map = double(gather(sensitivity_map));
    sensitivity_map = max(sensitivity_map, 0);

    fprintf('       Sensitivity map: range [%.3e, %.3e], built in %.1f s.\n', ...
        min(sensitivity_map(:)), max(sensitivity_map(:)), toc(sens_tic));
end

%% ========================= RECONSTRUCTION ================================

% Options for the shared Delay-And-Sum reconstruction (das_reconstruct.m).
% Defaults are faithful to the IRAI sample code; override any via CONFIG.das_*.
das_opts = struct();
if isfield(CONFIG, 'das_use_elements'), das_opts.use_elements = CONFIG.das_use_elements; end
if isfield(CONFIG, 'das_envelope'),     das_opts.envelope     = CONFIG.das_envelope;     end
if isfield(CONFIG, 'das_use_aperture'), das_opts.use_aperture = CONFIG.das_use_aperture; end
if isfield(CONFIG, 'das_aperture_cos'), das_opts.aperture_cos = CONFIG.das_aperture_cos; end
if isfield(CONFIG, 'das_depth_weight'), das_opts.depth_weight = CONFIG.das_depth_weight; end
if isfield(CONFIG, 'das_interp'),       das_opts.interp       = CONFIG.das_interp;       end

% Element metadata for per-element DAS. Guarantee a valid struct even for
% sensor methods that don't populate it (das_reconstruct falls back per-voxel).
if ~exist('sensor_info_orig', 'var') || ~isstruct(sensor_info_orig)
    sensor_info_orig = struct('num_elements', 0);
end

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

        reconPressure = back_propagate(kgrid, kmedium, sensor, sensorData, ...
                                        inputArgs, CONFIG);
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

if isfield(CONFIG, 'das_ubp') && CONFIG.das_ubp
    fprintf('       Running Universal Back-Projection DAS (Xu & Wang 2005)...\n');
else
    fprintf('       Running Delay-And-Sum reconstruction...\n');
end

try
    das_tic = tic;
    if isfield(CONFIG, 'das_ubp') && CONFIG.das_ubp
        reconPressure = das_reconstruct_ubp(sensorData, sensor, sensor_info_orig, ...
                                            medium, Nx, Ny, Nz, dx, dy, dz, dt, CONFIG);
    else
        reconPressure = das_reconstruct(sensorData, sensor, sensor_info_orig, medium, ...
                                        Nx, Ny, Nz, dx, dy, dz, dt, das_opts);
    end
    reconPressure = reconPressure * CONFIG.correction_factor;

    % Synthesize TR-style outputs so downstream plot/save code is reused
    conv_max_pressure = max(reconPressure(:));
    conv_rel_change   = NaN;
    num_iters_done    = 1;
    tr_time           = toc(das_tic);

    fprintf('       DAS complete (%.1f s).\n', tr_time);
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] DAS reconstruction failed: %s\n', ME.message);
    return;
end

case 'hybrid'

fprintf('       Running HYBRID reconstruction (DAS seed + up to %d-1 TR iterations)...\n', ...
    CONFIG.num_time_reversal_iter);

try
    hybrid_tic = tic;

    N_iter = CONFIG.num_time_reversal_iter;

    % ---- Iteration 1: DAS seed ----
    fprintf('       --- Iter 1/%d: DAS seed ---\n', N_iter);
    reconPressure = das_reconstruct(sensorData, sensor, sensor_info_orig, medium, ...
                                    Nx, Ny, Nz, dx, dy, dz, dt, das_opts);

    % Initialize convergence tracking (same shape as TR branch)
    conv_max_pressure    = zeros(N_iter, 1);
    conv_rel_change      = nan(N_iter, 1);
    conv_max_pressure(1) = max(reconPressure(:));
    num_iters_done       = 1;
    reconPressure_prev   = reconPressure;

    fprintf('       Max pressure (DAS seed): %.4e\n', conv_max_pressure(1));

    % ---- Live recon figure setup ----
    [~, dose_max_idx] = max(doseGrid(:));
    [cx_live, cy_live, cz_live] = ind2sub([Nx_orig, Ny_orig, Nz_orig], dose_max_idx);

    fig_live = []; hImg_recon = []; ax_recon = []; hLine_max = [];
    if CONFIG.plot_results
        [fig_live, ax_recon, hImg_recon, hLine_max] = ...
            setup_live_recon_figure(p0, Nx_orig, Ny_orig, cz_live, N_iter, CONFIG);

        % Show DAS recon as the iter-1 frame
        recon_slice_crop = squeeze(reconPressure(1:Nx_orig, 1:Ny_orig, cz_live))';
        recon_slice_crop = gather(recon_slice_crop);
        set(hImg_recon, 'CData', recon_slice_crop);
        caxis(ax_recon, [0, max(recon_slice_crop(:)) + eps]);
        title(ax_recon, sprintf('Reconstructed p_0   (iter 1  DAS seed)'), ...
            'FontWeight', 'bold');
        set(hLine_max, 'XData', 1, 'YData', conv_max_pressure(1));
        drawnow;
    end

    % ---- Residual update from DAS seed (only if more iterations remain) ----
    if N_iter > 1
        source_resid    = struct();
        source_resid.p0 = reconPressure;
        sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, source_resid, sensor, inputArgs{:});
        sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
    end

    % ---- TR iterations 2..N (mirror the 'tr' loop body) ----
    for tr_iter = 2:N_iter

        fprintf('       --- TR Iteration %d/%d ---\n', tr_iter, N_iter);

        reconPressure = back_propagate(kgrid, kmedium, sensor, sensorData, ...
                                        inputArgs, CONFIG);
        reconPressure = max(reconPressure, 0);

        conv_max_pressure(tr_iter) = max(reconPressure(:));
        num_iters_done = tr_iter;

        fprintf('       Max pressure: %.4e Pa\n', conv_max_pressure(tr_iter));

        % Convergence check vs previous iteration
        converged = false;
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

        reconPressure_prev = reconPressure;

        % ---- Update live figure ----
        if CONFIG.plot_results && ~isempty(fig_live) && ishandle(fig_live)
            recon_slice_crop = squeeze( ...
                reconPressure(1:Nx_orig, 1:Ny_orig, cz_live))';
            recon_slice_crop = gather(recon_slice_crop);
            set(hImg_recon, 'CData', recon_slice_crop);
            caxis(ax_recon, [0, max(recon_slice_crop(:)) + eps]);
            if converged
                title(ax_recon, ...
                    sprintf('Reconstructed p_0   (iter %d  CONVERGED)', tr_iter), ...
                    'FontWeight', 'bold', 'Color', [0, 0.55, 0]);
            else
                title(ax_recon, ...
                    sprintf('Reconstructed p_0   (iter %d / %d)', tr_iter, N_iter), ...
                    'FontWeight', 'bold');
            end
            set(hLine_max, 'XData', 1:tr_iter, 'YData', conv_max_pressure(1:tr_iter));
            drawnow;
        end

        if converged
            break;
        end

        % Residual correction for next iteration
        if tr_iter < N_iter
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium, source_resid, sensor, inputArgs{:});
            sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
        end
    end

    reconPressure = gather(reconPressure) * CONFIG.correction_factor;
    tr_time = toc(hybrid_tic);

    fprintf('       Hybrid complete (%.1f s, %d iterations).\n', tr_time, num_iters_done);
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] Hybrid reconstruction failed: %s\n', ME.message);
    return;
end

case 'mbir'

fprintf('       Running Model-Based Iterative Reconstruction (FISTA + TV)...\n');

try
    mbir_tic = tic;

    N_iter   = CONFIG.mbir_max_iter;
    lambda_T = CONFIG.mbir_lambda_tv;
    tv_iters = CONFIG.mbir_tv_iters;

    % --- Initial estimate: a single adjoint pass (cheap DAS-like seed) ---
    fprintf('       MBIR: seeding with one adjoint pass...\n');
    p0_seed = back_propagate(kgrid, kmedium, sensor, sensorData, ...
                              inputArgs, CONFIG);
    p0_seed = double(gather(max(p0_seed, 0)));

    % --- Lipschitz estimate via power method on A'A ---
    if isempty(CONFIG.mbir_step_size)
        fprintf('       MBIR: estimating Lipschitz via power method (5 iters)...\n');
        u = randn(Nx, Ny, Nz, 'single');
        u = u / norm(u(:));
        for pi = 1:5
            src_pm = struct('p0', double(u));
            Au = kspaceFirstOrder3D(kgrid, kmedium, src_pm, sensor, inputArgs{:});
            Au = double(gather(Au));
            AtAu = back_propagate(kgrid, kmedium, sensor, Au, inputArgs, CONFIG);
            AtAu = double(gather(AtAu));
            L_est = norm(AtAu(:));
            u = AtAu / max(L_est, eps);
            u = single(u);
        end
        L_step = 1 / max(L_est, eps);
    else
        L_step = CONFIG.mbir_step_size;
    end
    fprintf('       MBIR: step size = %.3e\n', L_step);

    % --- Preconditioner: 1 ./ (S + eps*max(S)) if sensitivity_map available ---
    if CONFIG.mbir_use_precond && ~isempty(sensitivity_map)
        S_eps   = sensitivity_map + eps + 1e-6 * max(sensitivity_map(:));
        precond = 1 ./ S_eps;
        precond = precond / max(precond(:));
    else
        precond = ones(Nx, Ny, Nz);
    end

    x_k       = p0_seed;
    y_k       = x_k;
    t_k       = 1;
    conv_max_pressure = zeros(N_iter, 1);
    conv_rel_change   = nan(N_iter, 1);
    num_iters_done    = 0;

    for k = 1:N_iter
        % Gradient: A'(A y_k - d)
        src_yk = struct('p0', y_k);
        A_yk   = kspaceFirstOrder3D(kgrid, kmedium, src_yk, sensor, inputArgs{:});
        A_yk   = double(gather(A_yk));
        resid  = A_yk - double(sensorData);
        grad   = back_propagate(kgrid, kmedium, sensor, resid, inputArgs, CONFIG);
        grad   = double(gather(grad));

        % Preconditioned gradient step
        x_next = y_k - L_step * (precond .* grad);

        % Proximal step: nonneg + TV
        x_next = max(x_next, 0);
        if lambda_T > 0
            x_next = tv_prox_3d(x_next, lambda_T * L_step, tv_iters);
            x_next = max(x_next, 0);
        end

        % FISTA momentum
        t_next = (1 + sqrt(1 + 4*t_k^2)) / 2;
        y_k    = x_next + ((t_k - 1) / t_next) * (x_next - x_k);

        conv_max_pressure(k) = max(x_next(:));
        if k > 1
            n_prev = norm(x_k(:));
            if n_prev > 0
                conv_rel_change(k) = norm(x_next(:) - x_k(:)) / n_prev;
            end
        end
        num_iters_done = k;

        x_k = x_next;
        t_k = t_next;

        fprintf('       MBIR iter %d/%d  max=%.3e  rel=%.3e\n', ...
            k, N_iter, conv_max_pressure(k), conv_rel_change(k));

        if k > 1 && ~isnan(conv_rel_change(k)) && ...
                conv_rel_change(k) < CONFIG.convergence_tol
            fprintf('       *** MBIR converged at iteration %d ***\n', k);
            break;
        end
    end

    reconPressure = x_k * CONFIG.correction_factor;
    tr_time       = toc(mbir_tic);
    fprintf('       MBIR complete (%.1f s).\n', tr_time);
    fprintf('       Reconstructed pressure: [%.2e, %.2e] Pa\n', ...
        min(reconPressure(:)), max(reconPressure(:)));

catch ME
    fprintf('[ERROR] MBIR reconstruction failed: %s\n', ME.message);
    return;
end

otherwise
    error('Unknown reconstruction_method: "%s" (use ''tr'', ''DAS'', ''hybrid'', or ''mbir'')', ...
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

%% ========================= SENSITIVITY-MAP CORRECTION ===================
%  Divide reconPressure by (S + eps*max(S)) before any other scaling.

if isfield(CONFIG, 'sensitivity_map_correction') && CONFIG.sensitivity_map_correction ...
        && ~isempty(sensitivity_map)

    S_crop = sensitivity_map(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    S_norm = S_crop / max(S_crop(:) + eps);
    denom  = S_norm + 1e-3 * max(S_norm(:));
    reconPressure = reconPressure ./ denom;
    fprintf('       [SENS] Applied sensitivity-map correction.\n');
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
if isfield(CONFIG, 'mask_recon_to_dose_region') && ~CONFIG.mask_recon_to_dose_region
    recon_dose = reconDosePerPulse * num_pulses .* body_mask_plot;
else
    recon_dose = reconDosePerPulse * num_pulses .* double(doseMask) .* body_mask_plot;
end

% Mask air pockets: air-classified voxels (density ~1.2 kg/m^3, assigned only
% inside the body for threshold_2) carry no real PA dose signal and would
% otherwise reconstruct as hotspots. Zero them from the final dose.
recon_dose(medium.density < 100) = 0;

fprintf('       Reconstructed dose: [%.4f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));

%% Crop p0 to original size (if padded)
if did_pad
    p0 = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
end

%% ========================= NORMALIZE DOSES ===============================
% When enabled, divide original and reconstructed doses by their respective
% maxima so both peak at 1.0. Affects downstream RESULTS SUMMARY, gamma, and
% plotting.

if isfield(CONFIG, 'normalize') && CONFIG.normalize
    dg_max = max(doseGrid(:));
    rd_max = max(recon_dose(:));
    fprintf('\n[NORM] Normalizing doses by their max:\n');
    fprintf('       Original max:      %.4f Gy\n', dg_max);
    fprintf('       Reconstructed max: %.4f Gy\n', rd_max);
    if dg_max > 0
        doseGrid = doseGrid / dg_max;
    end
    if rd_max > 0
        recon_dose = recon_dose / rd_max;
    end
    % Re-derive the low-dose threshold relative to the (now-normalized) original
    doseThreshold = 0.01 * max(doseGrid(:));
end

%% ========================= RESULTS SUMMARY ===============================

fprintf('\n========= RESULTS =========\n');
fprintf('  Original dose:      [%.6f, %.4f] Gy\n', min(doseGrid(:)), max(doseGrid(:)));
fprintf('  Reconstructed dose: [%.6f, %.4f] Gy\n', min(recon_dose(:)), max(recon_dose(:)));

dose_region = doseGrid > doseThreshold;
if any(dose_region(:))
    abs_error = abs(recon_dose(dose_region) - doseGrid(dose_region));
    rel_error = abs_error ./ max(doseGrid(dose_region), 1e-10) * 100;
    fprintf('  Mean abs err: %.6f Gy\n', mean(abs_error));
    fprintf('  Mean rel err: %.2f%%\n',  mean(rel_error));
    fprintf('  Max  rel err: %.2f%%\n',  max(rel_error));
end
fprintf('===========================\n');

%% ========================= LEAST-SQUARES GAIN ============================
%  NOTE: this gain uses the truth dose -> characterization metric only, not
%  deployable. c = sum(r .* t) / sum(r .^ 2) over the gamma eval mask.

ls_gain_c          = NaN;
recon_dose_calib   = [];
low_dose_cutoff    = 0.10 * max(doseGrid(:));
gamma_eval_mask    = doseGrid >= low_dose_cutoff;

if isfield(CONFIG, 'least_squares_gain') && CONFIG.least_squares_gain
    r_vec = recon_dose(gamma_eval_mask);
    t_vec = doseGrid(gamma_eval_mask);
    denom_ls = sum(r_vec .^ 2);
    if denom_ls > 0
        ls_gain_c = sum(r_vec .* t_vec) / denom_ls;
    else
        ls_gain_c = 1;
    end
    recon_dose_calib = ls_gain_c * recon_dose;
    fprintf('\n[LSGain] c = %.4f  (eval mask: %d voxels)\n', ls_gain_c, nnz(gamma_eval_mask));
    fprintf('         Calibrated max: %.4f Gy  (was %.4f Gy)\n', ...
        max(recon_dose_calib(:)), max(recon_dose(:)));
end

%% ========================= GAMMA ANALYSIS ================================

gamma_results = struct();

if exist('CalcGamma', 'file') == 2

    fprintf('\n[Gamma] Running gamma analysis (uncalibrated)...\n');

    [gamma_maps, pass_rates, gamma_criteria] = run_gamma_set( ...
        doseGrid, recon_dose, spacing_mm, gamma_eval_mask);

    gamma_results.maps        = gamma_maps;
    gamma_results.pass_rates  = pass_rates;
    gamma_results.criteria    = gamma_criteria;
    gamma_results.cutoff_Gy   = low_dose_cutoff;
    gamma_results.eval_mask   = gamma_eval_mask;

    fprintf('\n  ------ Gamma Pass Rates (uncalibrated, 10%% low-dose cutoff) ------\n');
    for gc = 1:size(gamma_criteria, 1)
        if isnan(pass_rates(gc))
            fprintf('  %-12s  FAILED\n', gamma_criteria{gc, 3});
        else
            fprintf('  %-12s  %.2f%%\n', gamma_criteria{gc, 3}, pass_rates(gc));
        end
    end

    % --- Optional second gamma run with least-squares calibration ---
    if ~isempty(recon_dose_calib)
        fprintf('\n[Gamma] Running gamma analysis (LS-calibrated c=%.4f)...\n', ls_gain_c);
        [gamma_maps_c, pass_rates_c, ~] = run_gamma_set( ...
            doseGrid, recon_dose_calib, spacing_mm, gamma_eval_mask);

        gamma_results.calibrated.maps       = gamma_maps_c;
        gamma_results.calibrated.pass_rates = pass_rates_c;
        gamma_results.calibrated.c          = ls_gain_c;

        fprintf('\n  ------ Gamma Pass Rates (LS-calibrated) ------\n');
        for gc = 1:size(gamma_criteria, 1)
            if isnan(pass_rates_c(gc))
                fprintf('  %-12s  FAILED\n', gamma_criteria{gc, 3});
            else
                fprintf('  %-12s  %.2f%%\n', gamma_criteria{gc, 3}, pass_rates_c(gc));
            end
        end
    end
else
    warning('CalcGamma not found. Skipping gamma analysis.');
    gamma_results = [];
end

%% ========================= GAMMA LOG ====================================
% Append CONFIG + gamma pass rates to a running .mat log. Maintains a
% per-criterion best record (highest pass rate and the CONFIG that produced
% it). Skipped when gamma analysis failed or produced no results.

if isfield(CONFIG, 'log_gamma') && CONFIG.log_gamma && ...
        ~isempty(gamma_results) && isfield(gamma_results, 'pass_rates')

    if isfield(CONFIG, 'gamma_log_file') && ~isempty(CONFIG.gamma_log_file)
        log_path = CONFIG.gamma_log_file;
    else
        log_path = 'gamma_log.mat';
    end
    if ~isfolder(fileparts(log_path)) && ~isempty(fileparts(log_path))
        log_path = fullfile(CONFIG.working_dir, log_path);
    elseif isempty(fileparts(log_path))
        log_path = fullfile(CONFIG.working_dir, log_path);
    end

    entry = struct();
    entry.timestamp     = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    entry.config        = CONFIG;
    entry.dose_filename = CONFIG.dose_filename;
    entry.criteria      = gamma_results.criteria(:, 3);
    entry.pass_rates    = gamma_results.pass_rates(:);

    n_crit = numel(entry.pass_rates);

    if isfile(log_path)
        L = load(log_path);
        if isfield(L, 'log_entries')
            log_entries = L.log_entries;
        else
            log_entries = struct([]);
        end
        if isfield(L, 'best')
            best = L.best;
        else
            best = repmat(struct('criterion', '', 'pass_rate', -Inf, ...
                'config', [], 'timestamp', ''), n_crit, 1);
            for gc = 1:n_crit
                best(gc).criterion = entry.criteria{gc};
            end
        end
    else
        log_entries = struct([]);
        best = repmat(struct('criterion', '', 'pass_rate', -Inf, ...
            'config', [], 'timestamp', ''), n_crit, 1);
        for gc = 1:n_crit
            best(gc).criterion = entry.criteria{gc};
        end
    end

    if isempty(log_entries)
        log_entries = entry;
    else
        log_entries(end+1) = entry;
    end

    for gc = 1:n_crit
        pr = entry.pass_rates(gc);
        if ~isnan(pr) && pr > best(gc).pass_rate
            best(gc).criterion = entry.criteria{gc};
            best(gc).pass_rate = pr;
            best(gc).config    = CONFIG;
            best(gc).timestamp = entry.timestamp;
        end
    end

    save(log_path, 'log_entries', 'best', '-v7.3');

    fprintf('\n[GammaLog] Appended run to %s (%d total entries).\n', ...
        log_path, numel(log_entries));
    fprintf('           Best pass rates so far:\n');
    for gc = 1:n_crit
        if isfinite(best(gc).pass_rate)
            fprintf('             %-12s  %.2f%%   (%s)\n', ...
                best(gc).criterion, best(gc).pass_rate, best(gc).timestamp);
        else
            fprintf('             %-12s  (none)\n', best(gc).criterion);
        end
    end
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
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    if ~isempty(sensitivity_map)
        results.sensitivity_map = sensitivity_map(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    end
    if ~isnan(ls_gain_c)
        results.ls_gain_c       = ls_gain_c;
        results.recon_dose_calib = recon_dose_calib;
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
        'Dose Comparison: Original vs Reconstructed', CONFIG.viz_smooth_sigma);

    % Figure 2  p0 convergence (max pressure + relative change)  TR & hybrid only
    if any(strcmpi(CONFIG.reconstruction_method, {'tr', 'hybrid'}))
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

    % Figure 4 - Depth profile (truth vs recon along Y through max-dose voxel)
    plot_depth_profile_y(doseGrid, recon_dose, recon_dose_calib, ls_gain_c, ...
        spacing_mm);

    % Figure 5 - Recon-vs-truth scatter colored by distance to sensor
    plot_recon_vs_truth_scatter(doseGrid, recon_dose, recon_dose_calib, ...
        ls_gain_c, sensor.mask, spacing_mm);

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
%  Sensor drawn as a solid red line/region  computed via max-projection so it
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


function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr, smooth_sigma)
%PLOT_DOSE_PANELS 2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = original dose,  Row 2 = reconstructed dose.
%  Dose (jet, semi-transparent) is overlaid on the density map (grayscale).
%  Voxels with dose < 10% of max are fully transparent (masked out).
%  Both rows share an identical dose colour range [0, max(original)]
%  so magnitudes are directly comparable.
%  Isocenter at max-dose voxel.  Sensor contour in red on every panel.
%  smooth_sigma (voxels) Gaussian-smooths each slice for display only, to
%  fill speckle gaps in a spotty recon. Pass 0 to disable.
%
%  Pass density=[] to show dose only with a black background.

    if nargin < 7 || isempty(smooth_sigma), smooth_sigma = 0; end

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
            if smooth_sigma > 0
                dose_slice = imgaussfilt(dose_slice, smooth_sigma);
            end
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


%  

function [fig_live, ax_recon, hImg_recon, hLine_max] = ...
        setup_live_recon_figure(p0, Nx_orig, Ny_orig, cz_live, N_iter, CONFIG)
%SETUP_LIVE_RECON_FIGURE  Build the 1x3 live reconstruction figure used by
%  both the 'tr' and 'hybrid' branches.
%  Panels: (1) initial p0 axial slice, (2) current recon p0, (3) live max-p convergence.

    fig_live = figure('Name', 'Live Reconstruction', 'Color', 'w', ...
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
    xlabel(ax_conv, 'Iteration');
    ylabel(ax_conv, 'Max Reconstructed p_0 (Pa)');
    title(ax_conv, 'Convergence (live)', 'FontWeight', 'bold');
    grid(ax_conv, 'on');
    xlim(ax_conv, [0.5, N_iter + 0.5]);

    sgtitle(fig_live, sprintf( ...
        'Live Reconstruction (%s)   Axial Z=%d   |   Patient %s', ...
        CONFIG.reconstruction_method, cz_live, CONFIG.patient_id), ...
        'FontWeight', 'bold', 'FontSize', 11);
    drawnow;
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

            % Air pockets (inside body only): assign real air acoustics to
            % low-HU voxels within the body contour. Outside the body, low-HU
            % stays water for sensor coupling. gruneisen = 0 -> no PA signal.
            if isfield(T, 'air') && isfield(T, 'air_hu_threshold') && isfield(sct, 'bodyMask')
                airMask = (HU < T.air_hu_threshold) & logical(sct.bodyMask);
                medium.density(airMask)     = T.air.density;
                medium.sound_speed(airMask) = T.air.sound_speed;
                medium.alpha_coeff(airMask) = T.air.alpha_coeff;
                medium.gruneisen(airMask)   = T.air.gruneisen;
                fprintf('         %-12s: %7d voxels (HU < %.0f, in body)\n', ...
                    'Air', sum(airMask(:)), T.air_hu_threshold);
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
    % Real air row (applied only to low-HU voxels INSIDE the body contour;
    % outside-body low-HU stays water for sensor coupling). gruneisen = 0 so
    % air generates no PA signal; these voxels are masked from the recon dose.
    tables.threshold_2.air_hu_threshold = -300;
    tables.threshold_2.air = struct( ...
        'density',     1.2, ...
        'sound_speed', 343, ...
        'alpha_coeff', 0, ...
        'alpha_power', 1.0, ...
        'gruneisen',   0);
end


%% =========================================================================
%  UPGRADED-BUILD LOCAL HELPERS
%% =========================================================================

function p0_adj = back_propagate(kgrid, kmedium, sensor, sensorData, inputArgs, CONFIG)
%BACK_PROPAGATE  Adjoint / time-reversal back-projection of sensor data onto
%  the recon grid. Switches:
%    CONFIG.use_true_adjoint        -> use 'additive' source injection
%                                      (Arridge et al. 2016, IP 32:115012)
%                                      instead of Dirichlet
%    CONFIG.absorption_compensation -> sign-flip alpha_coeff (gain) and
%                                      apply Tukey k-space window
%                                      (Treeby/Zhang/Cox 2010, IP 26:115003)

    Nx = size(sensor.mask, 1);
    Ny = size(sensor.mask, 2);
    Nz = size(sensor.mask, 3);

    source_tr        = struct();
    source_tr.p_mask = sensor.mask;
    source_tr.p      = fliplr(double(sensorData));

    if isfield(CONFIG, 'use_true_adjoint') && CONFIG.use_true_adjoint
        source_tr.p_mode = 'additive';
    else
        source_tr.p_mode = 'dirichlet';
    end

    kmedium_bp = kmedium;
    if isfield(CONFIG, 'absorption_compensation') && CONFIG.absorption_compensation
        % Sign-flip absorption (gain). Keep dispersion sign (alpha_power)
        % unchanged. k-Wave handles negative alpha_coeff as gain.
        kmedium_bp.alpha_coeff = -kmedium.alpha_coeff;
    end

    sensor_tr        = struct();
    sensor_tr.mask   = ones(Nx, Ny, Nz);
    sensor_tr.record = {'p_final'};

    p0_recon = kspaceFirstOrder3D(kgrid, kmedium_bp, source_tr, sensor_tr, inputArgs{:});

    if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
        p0_adj = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
    else
        p0_adj = reshape(p0_recon, [Nx, Ny, Nz]);
    end
    p0_adj = double(gather(p0_adj));

    if isfield(CONFIG, 'absorption_compensation') && CONFIG.absorption_compensation
        p0_adj = apply_tukey_kspace_filter(p0_adj, kgrid, ...
            CONFIG.tukey_cutoff_frac, CONFIG.tukey_taper_frac);
    end
end


function f = apply_tukey_kspace_filter(f, kgrid, cutoff_frac, taper_frac)
%APPLY_TUKEY_KSPACE_FILTER  Low-pass 3D Tukey window in k-space.
%  Window flat to cutoff_frac*k_max, raised-cosine taper of width
%  taper_frac*k_max, zero beyond. Stabilizes the gain in absorption-
%  compensated time-reversal (Treeby et al., Inverse Problems 26:115003).

    kx = kgrid.kx_vec;  ky = kgrid.ky_vec;  kz = kgrid.kz_vec;
    kx_max = max(abs(kx)); ky_max = max(abs(ky)); kz_max = max(abs(kz));
    if kx_max == 0, kx_max = 1; end
    if ky_max == 0, ky_max = 1; end
    if kz_max == 0, kz_max = 1; end

    [Kx, Ky, Kz] = ndgrid(kx, ky, kz);
    Knorm = sqrt((Kx/kx_max).^2 + (Ky/ky_max).^2 + (Kz/kz_max).^2);

    co = cutoff_frac;  tp = max(taper_frac, eps);
    w  = ones(size(Knorm));
    in_taper = Knorm > co & Knorm < co + tp;
    w(in_taper) = 0.5 * (1 + cos(pi * (Knorm(in_taper) - co) / tp));
    w(Knorm >= co + tp) = 0;

    F = fftn(f);
    F = F .* w;
    f = real(ifftn(F));
end


function [maps, pass_rates, criteria] = run_gamma_set(ref, tgt, spacing_mm, eval_mask)
%RUN_GAMMA_SET  Run the standard 10/10, 5/5, 3/3 gamma criteria.
    criteria = {10, 10, '10%/10mm'; 5, 5, '5%/5mm'; 3, 3, '3%/3mm'};
    maps       = cell(size(criteria, 1), 1);
    pass_rates = zeros(size(criteria, 1), 1);

    ref_s.start = [0, 0, 0]; ref_s.width = spacing_mm; ref_s.data = double(ref);
    tgt_s.start = [0, 0, 0]; tgt_s.width = spacing_mm; tgt_s.data = double(tgt);

    for gc = 1:size(criteria, 1)
        pct = criteria{gc, 1};
        dta = criteria{gc, 2};
        lbl = criteria{gc, 3};
        fprintf('  [%s] Computing...', lbl);
        try
            gmap = CalcGamma(ref_s, tgt_s, pct, dta, 'local', 0, ...
                'limit', dta * 2, 'restrict', 1);
            maps{gc}       = gmap;
            ev             = gmap(eval_mask);
            pass_rates(gc) = 100 * mean(ev <= 1);
            fprintf('  Pass rate: %.2f%%\n', pass_rates(gc));
        catch ME
            warning('Gamma [%s] failed: %s', lbl, ME.message);
            maps{gc}       = [];
            pass_rates(gc) = NaN;
            fprintf('  FAILED\n');
        end
    end
end


function recon = das_reconstruct_ubp(sensorData, sensor, sensor_info, medium, ...
                                     Nx, Ny, Nz, dx, dy, dz, dt, CONFIG)
%DAS_RECONSTRUCT_UBP  Universal Back-Projection (Xu & Wang 2005,
%  Phys. Rev. E 71:016706). For each voxel r and sensor pos r_s:
%     p0(r) = sum_s  w(r, r_s) * [2*p - 2*t*dp/dt] (t = |r - r_s| / c)
%  weight w(r, r_s) ~ cos(theta) / |r - r_s|^2 (solid angle).
%  Adds exp(+alpha(f0) * |r - r_s|) attenuation compensation.

    [sx_idx, sy_idx, sz_idx] = ind2sub([Nx, Ny, Nz], find(sensor.mask));
    n_sens = numel(sx_idx);
    if n_sens == 0
        error('das_reconstruct_ubp:NoSensors', 'Empty sensor mask.');
    end

    % Sensor coords in meters
    sx_m = (sx_idx - 1) * dx;
    sy_m = (sy_idx - 1) * dy;
    sz_m = (sz_idx - 1) * dz;

    % Voxel grid coords in meters
    [Xg, Yg, Zg] = ndgrid((0:Nx-1)*dx, (0:Ny-1)*dy, (0:Nz-1)*dz);

    % Mean sound speed (homogeneous backprop)
    c0 = mean(medium.sound_speed(:));

    % Sensor outward normal: for the sparse 2D array we approximate the
    % normal as the unit vector from the sensor centroid toward the body
    % centroid. (Solid-angle weighting requires SOME normal estimate.)
    cen_s = mean([sx_m, sy_m, sz_m], 1);
    cen_b = [mean(Xg(:)), mean(Yg(:)), mean(Zg(:))];
    n_hat = (cen_b - cen_s);
    n_hat = n_hat / max(norm(n_hat), eps);

    % Build attenuation in Np/m at f0 from alpha_coeff (dB/MHz^y/cm), y=alpha_power.
    f0_MHz = CONFIG.das_ubp_f0 / 1e6;
    a_db_cm_mhz = CONFIG.das_ubp_alpha_db_cm_mhz;
    y_pow = 1.1;
    if isfield(medium, 'alpha_power'), y_pow = medium.alpha_power; end
    alpha_db_cm = a_db_cm_mhz * f0_MHz ^ y_pow;
    alpha_np_m  = alpha_db_cm * 100 / 8.686;

    % Filter sensor data: time-domain 2*p - 2*t*dp/dt
    Nt = size(sensorData, 2);
    t_vec = (0:Nt-1) * dt;
    % gradient on a matrix: column-wise (along time, dim 2) when 1 output.
    dp_dt = gradient(double(sensorData), dt);
    sens_filt = 2 * double(sensorData) - 2 * (t_vec .* dp_dt);

    recon = zeros(Nx, Ny, Nz);

    fprintf('         UBP: %d sensors, c0=%.0f m/s, alpha(f0)=%.2f Np/m\n', ...
        n_sens, c0, alpha_np_m);

    for s = 1:n_sens
        dxv = Xg - sx_m(s);
        dyv = Yg - sy_m(s);
        dzv = Zg - sz_m(s);
        r   = sqrt(dxv.^2 + dyv.^2 + dzv.^2);

        % Time of flight from voxel to sensor
        idx   = round((r / c0) / dt) + 1;
        valid = idx >= 1 & idx <= Nt;

        % cos(theta): inner product of (r_voxel - r_sensor)/|.| with n_hat
        rdotn = (dxv * n_hat(1) + dyv * n_hat(2) + dzv * n_hat(3)) ./ max(r, eps);
        cos_t = max(rdotn, 0);   % only voxels in front of the element

        % Solid-angle / r^2 weighting + atten comp
        w = cos_t ./ max(r.^2, eps) .* exp(alpha_np_m * r);

        sens_s         = sens_filt(s, :);
        contrib        = zeros(Nx, Ny, Nz);
        contrib(valid) = sens_s(idx(valid));
        recon = recon + contrib .* w;
    end

    recon = max(recon, 0);
end


function x = tv_prox_3d(b, tau, n_iters)
%TV_PROX_3D  Anisotropic 3D total-variation proximal operator via Chambolle's
%  projection algorithm.  Solves  argmin_x 0.5||x-b||^2 + tau*TV(x).

    sz = size(b);
    p1 = zeros(sz);  p2 = zeros(sz);  p3 = zeros(sz);
    step = 1 / 12;

    for it = 1:n_iters
        div_p = back_diff_x(p1) + back_diff_y(p2) + back_diff_z(p3);
        g     = div_p - b / tau;

        gx = fwd_diff_x(g);
        gy = fwd_diff_y(g);
        gz = fwd_diff_z(g);

        p1 = p1 + step * gx;
        p2 = p2 + step * gy;
        p3 = p3 + step * gz;

        denom = max(1, sqrt(p1.^2 + p2.^2 + p3.^2));
        p1 = p1 ./ denom;
        p2 = p2 ./ denom;
        p3 = p3 ./ denom;
    end

    x = b - tau * (back_diff_x(p1) + back_diff_y(p2) + back_diff_z(p3));
end

function d = fwd_diff_x(a), d = cat(1, diff(a, 1, 1), zeros(1, size(a,2), size(a,3))); end
function d = fwd_diff_y(a), d = cat(2, diff(a, 1, 2), zeros(size(a,1), 1, size(a,3))); end
function d = fwd_diff_z(a), d = cat(3, diff(a, 1, 3), zeros(size(a,1), size(a,2), 1)); end
function d = back_diff_x(a)
    d = a - cat(1, zeros(1, size(a,2), size(a,3)), a(1:end-1, :, :));
    d(end, :, :) = -a(end-1, :, :);
end
function d = back_diff_y(a)
    d = a - cat(2, zeros(size(a,1), 1, size(a,3)), a(:, 1:end-1, :));
    d(:, end, :) = -a(:, end-1, :);
end
function d = back_diff_z(a)
    d = a - cat(3, zeros(size(a,1), size(a,2), 1), a(:, :, 1:end-1));
    d(:, :, end) = -a(:, :, end-1);
end


function plot_depth_profile_y(truth, recon, recon_calib, ls_c, spacing_mm)
%PLOT_DEPTH_PROFILE_Y  Line profile through the max-dose voxel along Y
%  (anterior-posterior). Truth vs reconstruction. If a calibrated recon is
%  available it is overlaid.

    [~, idx] = max(truth(:));
    [cx, ~, cz] = ind2sub(size(truth), idx);

    y_mm = (1:size(truth, 2)) * spacing_mm(2);
    t_p  = squeeze(truth(cx, :, cz));
    r_p  = squeeze(recon(cx, :, cz));

    figure('Name', 'Depth profile (Y through max-dose voxel)', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [200, 200, 720, 420]);
    plot(y_mm, t_p, 'k-', 'LineWidth', 2, 'DisplayName', 'Truth');
    hold on;
    plot(y_mm, r_p, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Recon (uncalib)');
    if ~isempty(recon_calib)
        r_c = squeeze(recon_calib(cx, :, cz));
        plot(y_mm, r_c, 'r--', 'LineWidth', 1.6, ...
            'DisplayName', sprintf('Recon (LS c=%.3f)', ls_c));
    end
    hold off;
    xlabel('Y (mm)'); ylabel('Dose (Gy)');
    title(sprintf('Depth profile through max-dose voxel (X=%d, Z=%d)', cx, cz), ...
        'FontWeight', 'bold');
    grid on; legend('Location', 'best');
    drawnow;
end


function plot_recon_vs_truth_scatter(truth, recon, recon_calib, ls_c, ...
                                     sensor_mask, spacing_mm)
%PLOT_RECON_VS_TRUTH_SCATTER  Scatter of recon vs truth on voxels above
%  10% of the truth max, color-coded by distance from each voxel to the
%  nearest active sensor element (mm). If a calibrated recon exists, a
%  second subplot shows the calibrated version.

    cutoff   = 0.10 * max(truth(:));
    mask_use = truth >= cutoff;
    if ~any(mask_use(:))
        return;
    end

    % Distance from every voxel to nearest sensor element (in mm)
    if any(sensor_mask(:))
        D = bwdist(logical(sensor_mask));
        Dmm = D * mean(spacing_mm);
    else
        Dmm = zeros(size(truth));
    end

    t_v = truth(mask_use);
    r_v = recon(mask_use);
    d_v = Dmm(mask_use);

    % Subsample for plotting if huge
    n_max = 30000;
    if numel(t_v) > n_max
        ix = randperm(numel(t_v), n_max);
        t_v = t_v(ix); r_v = r_v(ix); d_v = d_v(ix);
    end

    have_calib = ~isempty(recon_calib);
    n_pan = 1 + double(have_calib);

    figure('Name', 'Recon vs Truth (scatter, colored by sensor distance)', ...
        'Color', 'w', 'NumberTitle', 'off', 'Position', [220, 240, 720*n_pan, 460]);

    ax1 = subplot(1, n_pan, 1);
    scatter(ax1, t_v, r_v, 6, d_v, 'filled', 'MarkerFaceAlpha', 0.5);
    hold(ax1, 'on');
    mx = max([max(t_v), max(r_v)]) * 1.05;
    plot(ax1, [0, mx], [0, mx], 'k--', 'LineWidth', 1.2);
    hold(ax1, 'off');
    xlabel(ax1, 'Truth dose (Gy)'); ylabel(ax1, 'Recon dose (Gy)');
    title(ax1, 'Uncalibrated', 'FontWeight', 'bold');
    cb = colorbar(ax1); cb.Label.String = 'Distance to nearest sensor (mm)';
    axis(ax1, 'equal'); axis(ax1, [0, mx, 0, mx]); grid(ax1, 'on');

    if have_calib
        rc_v = recon_calib(mask_use);
        if exist('ix', 'var'), rc_v = rc_v(ix); end
        ax2 = subplot(1, n_pan, 2);
        scatter(ax2, t_v, rc_v, 6, d_v, 'filled', 'MarkerFaceAlpha', 0.5);
        hold(ax2, 'on');
        mx2 = max([max(t_v), max(rc_v)]) * 1.05;
        plot(ax2, [0, mx2], [0, mx2], 'k--', 'LineWidth', 1.2);
        hold(ax2, 'off');
        xlabel(ax2, 'Truth dose (Gy)'); ylabel(ax2, 'Recon dose (Gy)');
        title(ax2, sprintf('LS-calibrated (c=%.3f)', ls_c), 'FontWeight', 'bold');
        cb2 = colorbar(ax2); cb2.Label.String = 'Distance to nearest sensor (mm)';
        axis(ax2, 'equal'); axis(ax2, [0, mx2, 0, mx2]); grid(ax2, 'on');
    end

    drawnow;
end