%% =========================================================================
%  RUN_STANDALONE_COMPARISON.m
%  Two-dose k-Wave photoacoustic comparison.
%  Runs the standalone forward + reconstruction pipeline (identical to
%  run_standalone_simulation.m) on the listed dose file AND on its counterpart
%  on the other CT image, then performs three gamma analyses (10%/10mm,
%  5%/5mm, 3%/3mm) between the two reconstructed doses.
%  =========================================================================

clear; clc; close all;

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';

CONFIG.dose_filename = 'dose_1194203_Session_1_reference_CT_3_B15_112.mat';

% The two CT image indices to compare. The listed dose file above carries a
% _CT_k token selecting one of these; this script automatically locates that
% dose file's counterpart on the OTHER CT image, runs the identical forward +
% reconstruction pipeline on both, and gamma-compares the two reconstructions.
% Each dose is FORWARD-propagated on its OWN CBCT geometry (CBCT<k>_resampled.mat).
CONFIG.ct_pair = [1, 3];

% Reference geometry we experimentally have access to and RECONSTRUCT on.
% A dose whose own CT index differs from this is reconstructed "blind": the
% forward simulation (waves reaching the detector) runs on the dose's TRUE CT
% geometry, but iterative time reversal switches to this reference CT, because
% in a live IRAI system we never have the adapted (CT_3) geometry to invert on.
CONFIG.reference_ct = 1;

CONFIG.dose_file_override = '';
CONFIG.cbct_file_override = '';

CONFIG.sensor_placement_method = 'determine_sensor_mask';
CONFIG.sensor_x_index = 2;
CONFIG.sensor_y_index = 4;

% Physical 2D ultrasound array geometry (sparse element mask).
% Kerf is derived inside determine_sensor_mask as (pitch - size).
CONFIG.elements_per_side = 32;
CONFIG.element_pitch_mm  = 4.35; % Kerf is .7mm ?
CONFIG.element_size_mm   = 3.65;

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
CONFIG.correction_factor = 20;
%CONFIG.correction_factor = .0229;
%CONFIG.correction_factor = 0;
CONFIG.use_pressure_scale_correction = false;   % divide max(p0) / max(recon_pressure) before dose conversion
CONFIG.mask_recon_to_dose_region     = true;    % zero recon dose outside the dose mask (>1% of original max). Set false to keep the full reconstruction.

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
CONFIG.conv_noise_level    = 0.125;   % Noise amplitude as fraction of peak sensor signal
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

%% ===================== RESOLVE DOSE PAIR & CBCT PATHS ====================
%  Resolve the listed dose file (A), then derive its counterpart (B) on the
%  other CT image by swapping the _CT_k token. Each dose is paired with its
%  own CBCT geometry: CBCT<k>_resampled.mat.

if ~isempty(CONFIG.dose_file_override)
    dose_filepath_A = CONFIG.dose_file_override;
    [processed_dir, baseA, extA] = fileparts(dose_filepath_A);
    dose_basename_A = [baseA, extA];
else
    processed_dir   = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
    dose_basename_A = CONFIG.dose_filename;
    dose_filepath_A = fullfile(processed_dir, dose_basename_A);
end

% Identify this dose's CT image index from the _CT_k token.
ct_tok = regexp(dose_basename_A, '_CT_(\d+)_', 'tokens', 'once');
if isempty(ct_tok)
    error('run_standalone_comparison:NoCTtoken', ...
        ['Dose filename "%s" has no _CT_k token; cannot locate its ' ...
         'counterpart on the other CT image.'], dose_basename_A);
end
ct_a = str2double(ct_tok{1});
if ~ismember(ct_a, CONFIG.ct_pair)
    error('run_standalone_comparison:CTnotInPair', ...
        'Dose CT index %d is not in CONFIG.ct_pair = [%s].', ...
        ct_a, num2str(CONFIG.ct_pair));
end
ct_b_all = CONFIG.ct_pair(CONFIG.ct_pair ~= ct_a);
ct_b     = ct_b_all(1);

% Per-iteration CT index (di=1 -> ct_a listed, di=2 -> ct_b counterpart).
ct_list = [ct_a, ct_b];

% Validate the reconstruction reference CT and resolve its CBCT geometry. Any
% dose whose own CT differs from reference_ct is reconstructed on this CBCT.
if ~isfield(CONFIG, 'reference_ct') || isempty(CONFIG.reference_ct)
    CONFIG.reference_ct = 1;
end
if ~ismember(CONFIG.reference_ct, CONFIG.ct_pair)
    error('run_standalone_comparison:BadReferenceCT', ...
        'CONFIG.reference_ct = %d must be one of CONFIG.ct_pair = [%s].', ...
        CONFIG.reference_ct, num2str(CONFIG.ct_pair));
end
recon_cbct_filepath = fullfile(processed_dir, ...
    sprintf('CBCT%d_resampled.mat', CONFIG.reference_ct));

% Counterpart dose file: same name with the CT token swapped.
dose_basename_B = regexprep(dose_basename_A, '_CT_\d+_', sprintf('_CT_%d_', ct_b));
dose_filepath_B = fullfile(processed_dir, dose_basename_B);

% CBCT geometry for each dose (per-dose, not the legacy CBCT1-always default).
if ~isempty(CONFIG.cbct_file_override)
    cbct_filepath_A = CONFIG.cbct_file_override;
else
    cbct_filepath_A = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_a));
end
cbct_filepath_B = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_b));

dose_filepath_list = {dose_filepath_A, dose_filepath_B};
cbct_filepath_list = {cbct_filepath_A, cbct_filepath_B};
label_list         = {sprintf('CT_%d (listed)', ct_a), ...
                      sprintf('CT_%d (counterpart)', ct_b)};

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
fprintf('  Standalone k-Wave Photoacoustic Comparison  (v4.1)\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose A:          %s\n', dose_filepath_list{1});
fprintf('    on CBCT:       %s\n', cbct_filepath_list{1});
fprintf('  Dose B:          %s\n', dose_filepath_list{2});
fprintf('    on CBCT:       %s\n', cbct_filepath_list{2});
fprintf('  Sensor:          %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  Recon method:    %s\n', CONFIG.reconstruction_method);
fprintf('  TR iterations:   %d (tol: %.1e)\n', CONFIG.num_time_reversal_iter, CONFIG.convergence_tol);
fprintf('  GPU:             %s\n', mat2str(CONFIG.use_gpu));
if CONFIG.downscale_factor ~= 1
    fprintf('  Downscale factor: %g\n', CONFIG.downscale_factor);
end
fprintf('=========================================================\n\n');

%% ===================== TISSUE TABLES (shared medium builder) =============
% Tissue assignment uses the shared create_acoustic_medium (the same function
% the pipeline uses), which carries the threshold_2 air row. Build the lookup
% tables once and map the standalone uniform-medium config into the uniform
% table so the 'uniform' method reproduces the previous create_medium behavior.
CONFIG.tissue_tables = define_tissue_tables();
CONFIG.tissue_tables.uniform = struct( ...
    'density',     CONFIG.uniform_density, ...
    'sound_speed', CONFIG.uniform_sound_speed, ...
    'alpha_coeff', CONFIG.uniform_alpha_coeff, ...
    'alpha_power', 1.1, ...
    'gruneisen',   CONFIG.uniform_gruneisen);

%% ===================== PER-DOSE PIPELINE (RUN TWICE) =====================
%  Iterate over the two dose files. Each iteration runs the full, identical
%  forward + reconstruction pipeline and stores its reconstructed dose in
%  RES(di). After the loop the two reconstructions are gamma-compared.

for di = 1:2

clearvars fd_gantry_angle step15_spacing_mm sensor_info_orig

dose_filepath = dose_filepath_list{di};
cbct_filepath = cbct_filepath_list{di};
label         = label_list{di};

% Blind reconstruction: TRUE (forward) geometry is this dose's own CT, but if
% that differs from the reference CT we time-reverse on the reference geometry.
ct_this     = ct_list(di);
blind_recon = (ct_this ~= CONFIG.reference_ct);

fprintf('\n##################### PROCESSING DOSE %d/2: %s #####################\n', di, label);
if blind_recon
    fprintf('   [BLIND] Forward on CT_%d, reconstruction on reference CT_%d.\n', ...
        ct_this, CONFIG.reference_ct);
else
    fprintf('   [FULL ACCESS] Forward and reconstruction both on CT_%d.\n', ct_this);
end

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

fprintf('[2/7] Loading CBCT data for %s...\n', label); %  CT_1  temporary standalone default)...\n');
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
medium = create_acoustic_medium(sct, CONFIG);
% Standalone-only overrides (force-uniform toggles + coupling bath) applied
% after the shared builder, matching run_standalone_simulation.
medium = apply_standalone_medium_overrides(medium, sct, CONFIG);

fprintf('       Density range:     [%.0f, %.0f] kg/m^3\n', min(medium.density(:)), max(medium.density(:)));
fprintf('       Sound speed range: [%.0f, %.0f] m/s\n',    min(medium.sound_speed(:)), max(medium.sound_speed(:)));
fprintf('       Gruneisen range:   [%.4f, %.4f]\n',        min(medium.gruneisen(:)), max(medium.gruneisen(:)));

% --- Reconstruction medium (BLIND GEOMETRY HANDLING) ---
% `medium` above is the TRUE geometry, used to generate the sensor data the
% detector would actually measure (forward simulation). Iterative time reversal
% only has access to the reference CT, so build a SEPARATE reconstruction
% medium from the reference CBCT when this dose is on a different CT. Both media
% share the same grid, so every downstream padding/expansion/crop is applied to
% both in lockstep (see pad_medium_to / expand_medium / crop_medium helpers).
if blind_recon
    fprintf('       [BLIND] Loading reference-CT medium from %s\n', recon_cbct_filepath);
    if ~isfile(recon_cbct_filepath)
        error('Reference CBCT for reconstruction not found: %s', recon_cbct_filepath);
    end
    rec_cbct = load(recon_cbct_filepath);
    if isfield(rec_cbct, 'CBCT1_resampled')
        sct_recon = rec_cbct.CBCT1_resampled;
    elseif isfield(rec_cbct, 'CBCT3_resampled')
        sct_recon = rec_cbct.CBCT3_resampled;
    else
        error('CBCT1_resampled / CBCT3_resampled variable not found in %s', recon_cbct_filepath);
    end
    if ~isfield(sct_recon, 'cubeHU')
        error('Reference CBCT resampled struct missing required field: cubeHU');
    end
    % Match the forward sct grid (which may already be downscaled at this point).
    if CONFIG.downscale_factor ~= 1
        sct_recon.cubeHU = imresize3(sct_recon.cubeHU, size(sct.cubeHU));
        if isfield(sct_recon, 'bodyMask')
            sct_recon.bodyMask = imresize3(single(sct_recon.bodyMask), size(sct.cubeHU), 'nearest') > 0.5;
        end
        if isfield(sct_recon, 'couchMask')
            sct_recon.couchMask = imresize3(single(sct_recon.couchMask), size(sct.cubeHU), 'nearest') > 0.5;
        end
    elseif ~isequal(size(sct_recon.cubeHU), size(sct.cubeHU))
        error(['Reference CBCT grid [%s] does not match this dose''s CBCT grid [%s].\n' ...
               'Both CBCTs must be resampled onto the same dose grid.'], ...
            num2str(size(sct_recon.cubeHU)), num2str(size(sct.cubeHU)));
    end
    medium_recon = create_acoustic_medium(sct_recon, CONFIG);
    medium_recon = apply_standalone_medium_overrides(medium_recon, sct_recon, CONFIG);
    fprintf('       [BLIND] Recon medium density range: [%.0f, %.0f] kg/m^3\n', ...
        min(medium_recon.density(:)), max(medium_recon.density(:)));
else
    medium_recon = medium;   % full access: reconstruction geometry == forward geometry
end

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

    if blind_recon
        medium_recon = pad_medium_to(medium_recon, [Nx_pad, Ny_pad, Nz_pad]);
    else
        medium_recon = medium;
    end

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

            % Mirror the crop->water-expand onto the reference-CT medium so it
            % occupies the identical expanded grid. Nx_orig/Ny_orig/Nz_orig are
            % still the pre-expansion sizes here (updated below).
            if blind_recon
                medium_recon = crop_medium(medium_recon, [Nx_orig, Ny_orig, Nz_orig]);
                medium_recon = expand_medium(medium_recon, [Nx_exp, Ny_exp, Nz_exp], xr, yr, zr);
            end

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

                if blind_recon
                    medium_recon = pad_medium_to(medium_recon, [Nx_pad2, Ny_pad2, Nz_pad2]);
                end
            end

            if ~blind_recon
                medium_recon = medium;
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

% Forward-propagation medium = TRUE geometry (generates the measured signal).
kmedium_fwd             = struct();
kmedium_fwd.density     = medium.density;
kmedium_fwd.sound_speed = medium.sound_speed;
kmedium_fwd.alpha_coeff = medium.alpha_coeff;
kmedium_fwd.alpha_power = 1.1;

% Reconstruction medium = reference CT (what time reversal actually inverts on).
% Identical to the forward medium for full-access doses (medium_recon == medium).
kmedium_recon             = struct();
kmedium_recon.density     = medium_recon.density;
kmedium_recon.sound_speed = medium_recon.sound_speed;
kmedium_recon.alpha_coeff = medium_recon.alpha_coeff;
kmedium_recon.alpha_power = 1.1;

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
    sensorData = kspaceFirstOrder3D(kgrid, kmedium_fwd, source_fwd, sensor, inputArgs{:});
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
    %sensorData_noisy = wdenoise(sensorData_noisy)

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

        source_tr        = struct();
        source_tr.p_mask = sensor.mask;
        source_tr.p      = fliplr(sensorData);
        source_tr.p_mode = 'dirichlet';

        sensor_tr        = struct();
        sensor_tr.mask   = ones(Nx, Ny, Nz);
        sensor_tr.record = {'p_final'};

        p0_recon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_tr, sensor_tr, inputArgs{:});

        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
        else
            reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
        end

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
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_resid, sensor, inputArgs{:});
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

fprintf('       Running Delay-And-Sum reconstruction...\n');

try
    das_tic = tic;
    reconPressure = das_reconstruct(sensorData, sensor, sensor_info_orig, medium_recon, ...
                                    Nx, Ny, Nz, dx, dy, dz, dt, das_opts);
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
    reconPressure = das_reconstruct(sensorData, sensor, sensor_info_orig, medium_recon, ...
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
        sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_resid, sensor, inputArgs{:});
        sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
    end

    % ---- TR iterations 2..N (mirror the 'tr' loop body) ----
    for tr_iter = 2:N_iter

        fprintf('       --- TR Iteration %d/%d ---\n', tr_iter, N_iter);

        source_tr        = struct();
        source_tr.p_mask = sensor.mask;
        source_tr.p      = fliplr(sensorData);
        source_tr.p_mode = 'dirichlet';

        sensor_tr        = struct();
        sensor_tr.mask   = ones(Nx, Ny, Nz);
        sensor_tr.record = {'p_final'};

        p0_recon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_tr, sensor_tr, inputArgs{:});

        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = reshape(p0_recon.p_final, [Nx, Ny, Nz]);
        else
            reconPressure = reshape(p0_recon, [Nx, Ny, Nz]);
        end

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
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_resid, sensor, inputArgs{:});
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

otherwise
    error('Unknown reconstruction_method: "%s" (use ''tr'', ''DAS'', or ''hybrid'')', ...
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

% Pressure->dose conversion below uses `medium`. For a blind reconstruction the
% recovered pressure lives on the REFERENCE-CT geometry, so convert with the
% reference-CT medium (gruneisen/density), cropped to the reconstruction grid.
if blind_recon
    medium = crop_medium(medium_recon, gridSize);
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

%% ===================== CAPTURE PER-DOSE RESULTS ==========================
RES(di).recon_dose        = recon_dose;
RES(di).doseGrid          = doseGrid;
RES(di).spacing_mm        = spacing_mm;
RES(di).sensor_mask       = sensor.mask;
RES(di).density           = medium.density;
RES(di).p0                = p0;
RES(di).reconPressure     = reconPressure;
RES(di).gridSize          = gridSize;
RES(di).doseMask          = doseMask;
RES(di).num_pulses        = num_pulses;
RES(di).fwd_time          = fwd_time;
RES(di).tr_time           = tr_time;
RES(di).conv_max_pressure = conv_max_pressure;
RES(di).conv_rel_change   = conv_rel_change;
RES(di).num_iters_done    = num_iters_done;
RES(di).label             = label;

end   % di loop over the two dose files

%% ===================== EXTRACT DOSES & RECONSTRUCTIONS ===================

recon_A    = RES(1).recon_dose;   % listed dose reconstruction (for panels/save)
recon_B    = RES(2).recon_dose;   % counterpart dose reconstruction
spacing_mm = RES(1).spacing_mm;

if ~isequal(RES(1).gridSize, RES(2).gridSize)
    error('run_standalone_comparison:GridMismatch', ...
        ['The two doses are on different grids ([%s] vs [%s]).\n' ...
         'A common grid is required for the gamma comparisons.'], ...
        num2str(RES(1).gridSize), num2str(RES(2).gridSize));
end

% --- Map RES entries to CT-indexed doses used in the gamma comparisons ---
%   dose1 / recon1 : reference-CT dose (full-access reconstruction)
%   dose2 / recon2 : other-CT dose (blind reconstruction on the reference CT)
idx1 = find(ct_list == CONFIG.reference_ct, 1);
idx2 = find(ct_list ~= CONFIG.reference_ct, 1);
if isempty(idx2)
    tmp  = setdiff([1, 2], idx1);
    idx2 = tmp(1);
end

dose1  = double(RES(idx1).doseGrid);      % CT_ref truth (input dose)
recon1 = double(RES(idx1).recon_dose);    % CT_ref reconstruction
dose2  = double(RES(idx2).doseGrid);      % other-CT truth (input dose)
recon2 = double(RES(idx2).recon_dose);    % other-CT blind reconstruction
label1 = RES(idx1).label;
label2 = RES(idx2).label;
sensor1 = RES(idx1).sensor_mask;
sensor2 = RES(idx2).sensor_mask;

%% ========================= NORMALIZE DOSES ===============================
% When enabled, divide each volume by its own max so all peak at 1.0 before
% the gamma comparisons. Default (CONFIG.normalize=false) keeps absolute Gy so
% the recon-vs-truth comparisons remain physically meaningful.

if isfield(CONFIG, 'normalize') && CONFIG.normalize
    fprintf('\n[NORM] Normalizing each dose / reconstruction by its own max.\n');
    nrm = @(v) v / max([max(v(:)), eps]);
    dose1 = nrm(dose1); recon1 = nrm(recon1);
    dose2 = nrm(dose2); recon2 = nrm(recon2);
    recon_A = nrm(recon_A); recon_B = nrm(recon_B);
end

%% ========================= RESULTS SUMMARY ===============================

fprintf('\n========= COMPARISON RESULTS =========\n');
fprintf('  %-26s [%.6f, %.4f] Gy\n', [label1 ' truth (dose1):'],  min(dose1(:)),  max(dose1(:)));
fprintf('  %-26s [%.6f, %.4f] Gy\n', [label1 ' recon (recon1):'], min(recon1(:)), max(recon1(:)));
fprintf('  %-26s [%.6f, %.4f] Gy\n', [label2 ' truth (dose2):'],  min(dose2(:)),  max(dose2(:)));
fprintf('  %-26s [%.6f, %.4f] Gy\n', [label2 ' recon (recon2):'], min(recon2(:)), max(recon2(:)));
fprintf('=====================================\n');

%% ========================= GAMMA ANALYSIS ================================
% Four gamma comparisons (each at 10%/10mm, 5%/5mm, 3%/3mm):
%   1. dose1  vs dose2   truth change between the two plans/geometries
%   2. recon1 vs dose1   reference-CT detector accuracy
%   3. recon2 vs dose2   blind reconstruction accuracy (vs its own truth)
%   4. recon2 vs dose1   blind recon of the other-CT dose vs the reference truth
% Each pair is stored as {name, reference, target, title, sensor_mask}. The
% reference (ground truth) defines the 10% low-dose evaluation mask.

gamma_criteria = {3, 3, '3%/3mm'};

gamma_pairs = {
    'dose1_vs_dose2',  dose1,  dose2,  sprintf('%s truth  vs  %s truth', label1, label2), sensor1;
    'recon1_vs_dose1', dose1,  recon1, sprintf('%s recon  vs  %s truth', label1, label1), sensor1;
    'recon2_vs_dose2', dose2,  recon2, sprintf('%s recon  vs  %s truth', label2, label2), sensor2;
    'recon2_vs_dose1', dose1,  recon2, sprintf('%s recon  vs  %s truth', label2, label1), sensor2;
};

gamma_results = struct();

if exist('CalcGamma', 'file') == 2

    for pp = 1:size(gamma_pairs, 1)
        pair_name  = gamma_pairs{pp, 1};
        ref_vol    = gamma_pairs{pp, 2};
        tgt_vol    = gamma_pairs{pp, 3};
        pair_title = gamma_pairs{pp, 4};

        fprintf('\n[Gamma] %s\n', pair_title);
        gr = run_gamma_pair(ref_vol, tgt_vol, spacing_mm, gamma_criteria);
        gr.name  = pair_name;
        gr.title = pair_title;
        gr.ref   = ref_vol;
        gr.tgt   = tgt_vol;
        gr.sensor_mask = gamma_pairs{pp, 5};

        for gc = 1:size(gamma_criteria, 1)
            if isnan(gr.pass_rates(gc))
                fprintf('    %-12s  FAILED\n', gamma_criteria{gc, 3});
            else
                fprintf('    %-12s  %.2f%%\n', gamma_criteria{gc, 3}, gr.pass_rates(gc));
            end
        end

        gamma_results.(pair_name) = gr;
    end

    % Combined pass-rate table across all four comparisons.
    fprintf('\n  ------ Gamma Pass Rates (10%% low-dose cutoff on reference) ------\n');
    fprintf('  %-34s %10s %10s %10s\n', 'Comparison', gamma_criteria{1,3}, ...
        gamma_criteria{2,3}, gamma_criteria{3,3});
    pair_fn = fieldnames(gamma_results);
    for pp = 1:numel(pair_fn)
        gr = gamma_results.(pair_fn{pp});
        pr = gr.pass_rates;
        fprintf('  %-34s %9.2f%% %9.2f%% %9.2f%%\n', gr.title, pr(1), pr(2), pr(3));
    end
else
    warning('CalcGamma not found. Skipping gamma analysis.');
    gamma_results = [];
end

%% ========================= GAMMA LOG ====================================
% Append CONFIG + the four-comparison gamma pass-rate table to a running .mat
% log. Skipped when gamma analysis failed or produced no results.

if isfield(CONFIG, 'log_gamma') && CONFIG.log_gamma && ...
        ~isempty(gamma_results) && isstruct(gamma_results) && ~isempty(fieldnames(gamma_results))

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

    pair_fn = fieldnames(gamma_results);
    crit_labels = gamma_results.(pair_fn{1}).criteria(:, 3);

    % pass_rates: [nComparisons x nCriteria]; rows named by comparisons.
    PR = nan(numel(pair_fn), numel(crit_labels));
    for pp = 1:numel(pair_fn)
        PR(pp, :) = gamma_results.(pair_fn{pp}).pass_rates(:)';
    end

    entry = struct();
    entry.timestamp     = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    entry.config        = CONFIG;
    entry.dose_filename = CONFIG.dose_filename;
    entry.comparisons   = pair_fn;
    entry.criteria      = crit_labels;
    entry.pass_rates    = PR;

    if isfile(log_path)
        L = load(log_path);
        if isfield(L, 'log_entries')
            log_entries = L.log_entries;
        else
            log_entries = struct([]);
        end
    else
        log_entries = struct([]);
    end

    if isempty(log_entries)
        log_entries = entry;
    else
        log_entries(end+1) = entry;
    end

    save(log_path, 'log_entries', '-v7.3');

    fprintf('\n[GammaLog] Appended run to %s (%d total entries).\n', ...
        log_path, numel(log_entries));
end

%% ========================= SAVE RESULTS =================================

if CONFIG.save_results
    % Build filename from configuration parameters
    output_fname = sprintf('standalone_comparison_%s_%s_%s_%d.mat', ...
        CONFIG.reconstruction_method, ...
        CONFIG.sensor_placement_method, CONFIG.gruneisen_method, ...
        CONFIG.num_time_reversal_iter);

    results.recon_dose_A    = recon_A;
    results.recon_dose_B    = recon_B;
    results.label_A         = RES(1).label;
    results.label_B         = RES(2).label;
    results.original_dose_A = RES(1).doseGrid;
    results.original_dose_B = RES(2).doseGrid;
    results.p0_A            = RES(1).p0;
    results.p0_B            = RES(2).p0;
    results.reconPressure_A = RES(1).reconPressure;
    results.reconPressure_B = RES(2).reconPressure;
    results.sensor_mask_A   = RES(1).sensor_mask;
    results.sensor_mask_B   = RES(2).sensor_mask;
    results.config          = CONFIG;
    results.spacing_mm      = spacing_mm;
    results.grid_size       = RES(1).gridSize;
    results.fwd_time_sec    = [RES(1).fwd_time, RES(2).fwd_time];
    results.tr_time_sec     = [RES(1).tr_time,  RES(2).tr_time];
    if ~isempty(gamma_results)
        results.gamma = gamma_results;
    end
    save(output_fname, '-struct', 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', output_fname);
end

%% ========================= POST-SIMULATION VISUALIZATION ================

if CONFIG.plot_results

    % Figure 1  2x3 dose comparison (original top, recon bottom)
    %            Coronal | Sagittal | Axial
    %            Isocenter = max-dose voxel; sensor contour in red
    plot_dose_panels(recon_A, recon_B, RES(1).sensor_mask, RES(1).density, spacing_mm, ...
        'Dose Comparison: Two Reconstructed Doses', {RES(1).label, RES(2).label}, ...
        CONFIG.viz_smooth_sigma);

    % Figure 2  p0 convergence (max pressure + relative change)  TR & hybrid only
    if any(strcmpi(CONFIG.reconstruction_method, {'tr', 'hybrid'}))
        p0_max_for_plot = max(RES(1).p0(:));
        plot_convergence_history(RES(1).conv_max_pressure, RES(1).conv_rel_change, ...
            RES(1).num_iters_done, CONFIG.convergence_tol, p0_max_for_plot);
    end

    % Figure 3  Axial gamma (3 criteria) + absolute error
    if ~isempty(gamma_results) && isstruct(gamma_results)
        pair_fn = fieldnames(gamma_results);
        for k = 1:numel(pair_fn)
            gr = gamma_results.(pair_fn{k});
            [~, max_dose_idx] = max(gr.ref(:));
            [~, ~, cz_gamma]  = ind2sub(size(gr.ref), max_dose_idx);
            plot_gamma_and_error_axial(gr, gr.ref, gr.tgt, ...
                gr.sensor_mask, cz_gamma, gr.title);
        end
    end

end

fprintf('\nStandalone comparison complete.\n');


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


function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr, rowLabels, smooth_sigma)
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

    if nargin < 8 || isempty(smooth_sigma), smooth_sigma = 0; end

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

    if nargin < 7 || isempty(rowLabels)
        row_labels = {'Original', 'Reconstructed'};
    else
        row_labels = rowLabels;
    end
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


function plot_gamma_and_error_axial(gamma_results, original, recon, sensor_mask, cz, titleStr)
%PLOT_GAMMA_AND_ERROR_AXIAL 1x4 axial figure:
%  gamma 10/10 | gamma 5/5 | gamma 3/3 | absolute dose error.
%  Sensor contour in red.  Slice at cz (max-dose axial index).
%  original = gamma reference volume, recon = gamma target volume.
%  titleStr (optional) names the comparison in the figure title.

    if nargin < 6 || isempty(titleStr), titleStr = 'Gamma Index & Absolute Error'; end

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

    figure('Name', ['Gamma & Absolute Error - ' titleStr], 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [50, 300, 1400, 370]);
    sgtitle(sprintf('%s\nAxial Plane (Z = %d voxel)  Gamma Index & |Target - Reference|', titleStr, cz), ...
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


function medium = apply_standalone_medium_overrides(medium, sct, config)
%APPLY_STANDALONE_MEDIUM_OVERRIDES  Standalone-only medium adjustments applied
%  after the shared create_acoustic_medium: the force-uniform sensitivity
%  toggles, then force the coupling bath (outside the body / couch) to the
%  uniform medium. Mirrors run_standalone_simulation exactly.

    gs = medium.grid_size;
    if config.force_uniform_density
        medium.density = ones(gs) * config.uniform_density;
    end
    if config.force_uniform_sound_speed
        medium.sound_speed = ones(gs) * config.uniform_sound_speed;
    end
    if config.force_uniform_attenuation
        medium.alpha_coeff = ones(gs) * config.uniform_alpha_coeff;
        medium.alpha_power = 1.1;
    end
    if config.force_uniform_gruneisen
        medium.gruneisen = ones(gs) * config.uniform_gruneisen;
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


function m = pad_medium_to(m, sz)
%PAD_MEDIUM_TO Water-pad the four medium fields into an sz grid.
%  Copies each existing field into the leading sub-block of a fresh sz array
%  filled with the background (water: density 1000, c 1540, alpha/gruneisen 0),
%  mirroring the forward-medium FFT/grid padding so the reference-CT medium
%  stays voxel-aligned with it.
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


function gr = run_gamma_pair(ref_vol, tgt_vol, spacing_mm, criteria)
%RUN_GAMMA_PAIR Gamma index of tgt_vol against ref_vol at each criterion row.
%  ref_vol is the reference (ground truth): it sets the 10% low-dose evaluation
%  mask. criteria is an {pct, dta, label} cell array. Returns a struct with
%  .maps (per-criterion gamma maps), .pass_rates (% gamma<=1 within the mask),
%  .criteria, .cutoff_Gy, and .eval_mask.
    nCrit = size(criteria, 1);
    gr = struct('maps', {cell(nCrit, 1)}, 'pass_rates', nan(nCrit, 1), ...
                'criteria', {criteria}, 'cutoff_Gy', 0, 'eval_mask', []);

    ref_struct.start = [0, 0, 0]; ref_struct.width = spacing_mm; ref_struct.data = double(ref_vol);
    tgt_struct.start = [0, 0, 0]; tgt_struct.width = spacing_mm; tgt_struct.data = double(tgt_vol);

    cutoff    = 0.10 * max(ref_vol(:));
    eval_mask = ref_vol >= cutoff;
    gr.cutoff_Gy = cutoff;
    gr.eval_mask = eval_mask;

    for gc = 1:nCrit
        pct = criteria{gc, 1};
        dta = criteria{gc, 2};
        lbl = criteria{gc, 3};
        try
            gmap = CalcGamma(ref_struct, tgt_struct, pct, dta, ...
                'local', 0, 'limit', dta * 2, 'restrict', 1);
            gr.maps{gc}       = gmap;
            gr.pass_rates(gc) = 100 * mean(gmap(eval_mask) <= 1);
        catch ME
            warning('Gamma [%s] failed: %s', lbl, ME.message);
            gr.maps{gc}       = [];
            gr.pass_rates(gc) = NaN;
        end
    end
end