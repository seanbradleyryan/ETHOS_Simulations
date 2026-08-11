%% =========================================================================
%  STUDY_GAMMA_INDEX_CONVERGENCE.m
%  Gamma pass-rate vs number of time-reversal iterations.
%
%  Modeled on run_standalone_simulation.m / run_nt_convergence_sweep.m. Instead
%  of reconstructing once with a fixed iteration count, this study records the
%  3%/3 mm gamma pass rate (truth = reference, recon = target) AFTER EVERY time-
%  reversal iteration, for CONFIG.num_time_reversal_iter iterations, so the
%  reconstruction's convergence toward the truth can be read directly.
%
%  TWO STUDIES (identical sensor geometry + k-Wave grid for a fair comparison):
%
%    STUDY A - matched CT_1:
%       Forward-propagate the CT_1 field dose (D1) through the CT_1 acoustic
%       medium, then reconstruct on the CT_1 medium (matched). After each TR
%       iteration, gamma vs D1.  -> per-iteration detector accuracy.
%
%    STUDY B - cross-CT (CT_3 delivery, CT_1 model):
%       Forward-propagate the CT_3 counterpart field dose (D3) through the CT_3
%       acoustic medium (CT_3 geometry for the initial propagation), but
%       reconstruct on the CT_1 medium (CT_1 geometry for recon). After each TR
%       iteration, gamma of the recon vs D1 (the CT_1 truth). -> mimics adapted-
%       vs-reference change detection where the model anatomy (CT_1) differs from
%       the delivered anatomy (CT_3).
%
%  NORMALIZATION (from study_pass_rates_allsegments.m): before each gamma, the
%  recon is rescaled by the least-squares gain g = sum(truth.*recon)/sum(recon^2)
%  over the truth's >=10% region -- the g minimizing ||truth - g*recon||^2. The
%  truth (D1) is never rescaled. This cancels the stored absolute pressure scale
%  (CONFIG.correction_factor) so only the spatial agreement drives gamma.
%
%  DELIVERABLES:
%    Figure 1 - Study A: gamma pass rate vs TR iteration.
%    Figure 2 - Study B: gamma pass rate vs TR iteration.
%    Figure 3 - Both studies overlaid vs TR iteration.
%
%  NOTE: HIPAA / remote-execution - this file is WRITTEN here but must be RUN on
%  the remote device. Do not execute locally.
%% =========================================================================

function RESULTS = study_gamma_index_convergence(dose_file, CONFIG)
%STUDY_GAMMA_INDEX_CONVERGENCE  Gamma pass-rate convergence vs TR iteration for a
%  single field dose, plus truth-vs-recon dose / gamma comparison panels.
%
%  For one beam/segment (a CT_1 field dose and its derived CT_3 counterpart) this
%  runs both Study A (matched CT_1) and Study B (cross-CT: CT_3 delivery, CT_1
%  model) and, after EVERY time-reversal iteration, records the n%/n mm gamma
%  pass rate of the reconstructed dose against the CT_1 truth (D1). It also
%  snapshots the reconstructed dose and its gamma-index map at selected
%  iterations (default 1 and 10) for a truth-vs-recon comparison panel.
%
%  USAGE
%    study_gamma_index_convergence()                       % built-in defaults
%    study_gamma_index_convergence(dose_file)              % one CT_1 dose file
%    RESULTS = study_gamma_index_convergence(dose_file, CONFIG)
%
%  INPUTS
%    dose_file - [] / '' to keep CONFIG.dose_filename_ct1; a bare CT_1 field-dose
%                filename (looked up under the patient/session processed dir); or
%                a FULL PATH to a CT_1 field-dose .mat (its folder overrides the
%                processed dir). The CT_3 counterpart is derived by swapping
%                _CT_1_ -> _CT_3_ unless CONFIG.dose_filename_ct3 is set.
%    CONFIG    - optional struct; any field overrides the defaults returned by
%                default_config() (see the bottom of this file). Notable fields:
%                  .snapshot_iters  TR iterations to snapshot for the dose/gamma
%                                   panels (default [1 10], clamped to N_iter).
%                  .plot_results    draw figures (default true).
%                  .save_figures    save each figure as PNG to CONFIG.figure_dir.
%                  .close_figures   close each figure after saving (batch mode).
%                  .save_results    save the RESULTS struct to CONFIG.output_file.
%
%  OUTPUT
%    RESULTS - per-iteration pass rates (A/B), the snapshot recon doses + gamma
%              maps for both studies, and grid/metadata. Also written to disk when
%              CONFIG.save_results is true.
%
%  DELIVERABLES (figures):
%    - Study A / Study B truth-vs-recon dose + gamma comparison panels (iters 1,10)
%    - Gamma pass rate vs TR iteration: Study A, Study B, and both overlaid.
%
%  DEPENDENCIES: create_acoustic_medium, determine_sensor_mask,
%    find_optimal_kwave_size, kspaceFirstOrder3D, gaussianFilter, CalcGamma.
%  See also: run_standalone_simulation, study_pass_rates_individual,
%    run_gamma_convergence_batch.
%
%  NOTE: HIPAA / remote-execution - this file is WRITTEN here but must be RUN on
%  the remote device. Do not execute locally.
%% =========================================================================

study_timer = tic;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% ========================= CONFIGURATION ================================
% Start from the shared defaults, then overlay any caller-supplied overrides and
% the dose-file argument.

if nargin < 2 || isempty(CONFIG)
    CONFIG = struct();
end
CONFIG = merge_config(default_config(), CONFIG);

% Resolve the dose-file argument. A bare name is looked up under the processed
% dir; a full path additionally overrides that dir (used by the batch driver).
processed_dir_override = '';
if nargin >= 1 && ~isempty(dose_file)
    [dose_arg_dir, dose_arg_base, dose_arg_ext] = fileparts(char(dose_file));
    CONFIG.dose_filename_ct1 = [dose_arg_base, dose_arg_ext];
    if ~isempty(dose_arg_dir)
        processed_dir_override = dose_arg_dir;
    end
end

%% ========================= RESOLVE FILE PATHS ============================

if ~isempty(processed_dir_override)
    processed_dir = processed_dir_override;
else
    processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
        CONFIG.patient_id, CONFIG.session, 'processed');
end

% Derive the CT_3 dose filename from the CT_1 one unless explicitly given.
if isempty(CONFIG.dose_filename_ct3)
    CONFIG.dose_filename_ct3 = strrep(CONFIG.dose_filename_ct1, '_CT_1_', '_CT_3_');
    if strcmp(CONFIG.dose_filename_ct3, CONFIG.dose_filename_ct1)
        error('study_gamma_index_convergence:NoCT3Name', ...
            ['Could not derive a CT_3 dose filename from "%s" (no _CT_1_ token). ' ...
             'Set CONFIG.dose_filename_ct3 explicitly.'], CONFIG.dose_filename_ct1);
    end
end

dose_ct1_path = fullfile(processed_dir, CONFIG.dose_filename_ct1);
dose_ct3_path = fullfile(processed_dir, CONFIG.dose_filename_ct3);
cbct_ct1_path = fullfile(processed_dir, CONFIG.cbct_filename_ct1);
cbct_ct3_path = fullfile(processed_dir, CONFIG.cbct_filename_ct3);

% A short beam/segment tag (e.g. B15_112) parsed from the dose filename, used to
% name figures/results per field. Falls back to the base filename.
file_tag = dose_file_tag(CONFIG.dose_filename_ct1);

% Fill in figure / results output locations if the caller left them empty.
if isempty(CONFIG.figure_dir)
    CONFIG.figure_dir = fullfile(CONFIG.working_dir, 'AnalysisResults', ...
        CONFIG.patient_id, CONFIG.session, 'gamma_convergence');
end
if isempty(CONFIG.results_dir)
    CONFIG.results_dir = CONFIG.figure_dir;
end
if isempty(CONFIG.output_file)
    CONFIG.output_file = fullfile(CONFIG.results_dir, ...
        ['gamma_index_convergence_' file_tag '.mat']);
end

%% ========================= LOAD PLAN BEAM METADATA =======================
% determine_sensor_mask needs beam_metadata (isocenter + jaw extents per beam)
% to build the anterior-surface exclusion zone. Load it from metadata.mat.

if ~isfield(CONFIG, 'beam_metadata') || isempty(CONFIG.beam_metadata)
    metadata_filepath = fullfile(processed_dir, 'metadata.mat');
    if isfile(metadata_filepath)
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
fprintf('  Gamma Index TR-Convergence Study\n');
fprintf('=========================================================\n');
fprintf('  Patient:         %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Dose CT_1 (D1):  %s\n', CONFIG.dose_filename_ct1);
fprintf('  Dose CT_3 (D3):  %s\n', CONFIG.dose_filename_ct3);
fprintf('  Sensor:          %s\n', CONFIG.sensor_placement_method);
fprintf('  Tissue model:    %s\n', CONFIG.gruneisen_method);
fprintf('  TR iterations:   %d  (gamma %g%%/%g mm after each)\n', ...
    CONFIG.num_time_reversal_iter, CONFIG.gamma_n, CONFIG.gamma_n);
fprintf('  Normalization:   LSQ recon->truth  (from study_pass_rates_allsegments)\n');
fprintf('=========================================================\n\n');

if exist('CalcGamma', 'file') ~= 2
    error('study_gamma_index_convergence:NoCalcGamma', ...
        'CalcGamma not found on the path; cannot compute gamma pass rates.');
end

% Tissue tables (shared by both media). Mirrors run_standalone_simulation.
CONFIG.tissue_tables = define_tissue_tables();
CONFIG.tissue_tables.uniform = struct( ...
    'density',     CONFIG.uniform_density, ...
    'sound_speed', CONFIG.uniform_sound_speed, ...
    'alpha_coeff', CONFIG.uniform_alpha_coeff, ...
    'alpha_power', 1.1, ...
    'gruneisen',   CONFIG.uniform_gruneisen);

%% ========================= LOAD DATA ====================================

fprintf('[LOAD] CT_1 field dose (D1)...\n');
[D1, spacing1, meterset1, gantry1] = load_field_dose(dose_ct1_path, CONFIG);
fprintf('[LOAD] CT_3 field dose (D3)...\n');
[D3, ~, ~, ~] = load_field_dose(dose_ct3_path, CONFIG);

fprintf('[LOAD] CBCT_1 (CT_1 geometry)...\n');
sct1 = load_cbct(cbct_ct1_path);
fprintf('[LOAD] CBCT_3 (CT_3 geometry)...\n');
sct3 = load_cbct(cbct_ct3_path);

% Spacing: prefer CBCT, fall back to the dose file's embedded spacing.
if isfield(sct1, 'spacing') && ~isempty(sct1.spacing)
    spacing_mm = sct1.spacing(:)';
elseif ~isempty(spacing1)
    spacing_mm = spacing1(:)';
else
    error('study_gamma_index_convergence:NoSpacing', 'No spacing available from CBCT_1 or D1.');
end

% Meterset: the same beam/segment on either CT delivers the same MU, so use D1's.
if ~isempty(meterset1) && meterset1 > 0 && CONFIG.meterset ~= meterset1
    fprintf('  [INFO] Overriding CONFIG.meterset: %.2f -> %.2f MU\n', CONFIG.meterset, meterset1);
    CONFIG.meterset = meterset1;
end

% --- Grid-consistency checks: all four volumes must share one dose grid. ---
gridSize = size(D1);
assert_same_grid('D1',    gridSize, size(D1));
assert_same_grid('D3',    gridSize, size(D3));
assert_same_grid('CBCT1', gridSize, size(sct1.cubeHU));
assert_same_grid('CBCT3', gridSize, size(sct3.cubeHU));

Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);
fprintf('  Grid size: [%d x %d x %d], spacing [%.2f %.2f %.2f] mm\n', Nx, Ny, Nz, spacing_mm);

% Mask each dose to its own body contour.
if isfield(sct1, 'bodyMask'), D1 = D1 .* double(sct1.bodyMask); end
if isfield(sct3, 'bodyMask'), D3 = D3 .* double(sct3.bodyMask); end

dx = spacing_mm(1) / 1000;
dy = spacing_mm(2) / 1000;
dz = spacing_mm(3) / 1000;

%% ========================= BUILD BOTH MEDIA =============================

fprintf('\n[MEDIUM] Building CT_1 and CT_3 acoustic media (method: %s)...\n', ...
    CONFIG.gruneisen_method);
medium1 = build_medium(sct1, CONFIG);   % CT_1 (matched recon; Study A fwd+recon)
medium3 = build_medium(sct3, CONFIG);   % CT_3 (Study B forward propagation)

fprintf('  CT_1 density [%.0f, %.0f], c [%.0f, %.0f], Gr [%.3f, %.3f]\n', ...
    min(medium1.density(:)), max(medium1.density(:)), ...
    min(medium1.sound_speed(:)), max(medium1.sound_speed(:)), ...
    min(medium1.gruneisen(:)), max(medium1.gruneisen(:)));
fprintf('  CT_3 density [%.0f, %.0f], c [%.0f, %.0f], Gr [%.3f, %.3f]\n', ...
    min(medium3.density(:)), max(medium3.density(:)), ...
    min(medium3.sound_speed(:)), max(medium3.sound_speed(:)), ...
    min(medium3.gruneisen(:)), max(medium3.gruneisen(:)));

%% ========================= INITIAL PRESSURE p0 =========================
% Each study's source is built on its OWN delivery anatomy:
%   Study A: p0_A = D1_per_pulse .* Gr_CT1 .* rho_CT1
%   Study B: p0_B = D3_per_pulse .* Gr_CT3 .* rho_CT3

num_pulses = ceil(CONFIG.meterset / CONFIG.dose_per_pulse_cGy);
fprintf('\n[p0] Meterset %.2f MU -> %d pulses.\n', CONFIG.meterset, num_pulses);

p0_A_orig = compute_p0(D1, medium1, num_pulses);
p0_B_orig = compute_p0(D3, medium3, num_pulses);

% Low-dose masks (recon is only meaningful where the forward source had signal).
doseMaskA = D1 > 0.01 * max(D1(:));   % Study A forward = D1
doseMaskB = D3 > 0.01 * max(D3(:));   % Study B forward = D3

if max(p0_A_orig(:)) == 0 || max(p0_B_orig(:)) == 0
    error('study_gamma_index_convergence:NoSignal', 'One of the initial-pressure fields is all zero.');
end

%% ========================= SENSOR PLACEMENT (shared) ===================
% Placed ONCE from D1 + CT_1 so both studies share an identical sensor and grid
% (the same physical array on the patient). determine_sensor_mask may expand the
% grid to clear the beam exclusion zone; the expansion parameters are reused for
% both media / sources via the embed helpers below.

fprintf('\n[SENSOR] Placing sensor (method: %s)...\n', CONFIG.sensor_placement_method);

gp          = struct('expanded', false);
sensor_info = struct('num_elements', 0);
sensor_mask_orig = [];

switch CONFIG.sensor_placement_method
    case 'determine_sensor_mask'
        sct_for_sensor = sct1;
        if ~isfield(sct_for_sensor, 'couchMask')
            sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
        end
        if ~isfield(sct_for_sensor, 'origin')
            sct_for_sensor.origin = [0, 0, 0];
        end
        sct_for_sensor.spacing = spacing_mm;

        field_dose_for_sensor = struct();
        field_dose_for_sensor.dose_Gy = D1;
        if ~isempty(gantry1)
            field_dose_for_sensor.gantry_angle = gantry1;
        else
            field_dose_for_sensor.gantry_angle = 0;
        end
        field_dose_for_sensor.origin     = sct_for_sensor.origin;
        field_dose_for_sensor.spacing    = spacing_mm;
        field_dose_for_sensor.dimensions = [Nx, Ny, Nz];

        beam_meta = [];
        if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
            beam_meta = CONFIG.beam_metadata;
        end

        [sensor_mask_orig, sensor_info] = determine_sensor_mask( ...
            sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);
        gp = sensor_info.grid_pad;
    case 'full_plane_anterior'
        % handled after geometry resolution (no expansion)
    case 'full_plane_lateral'
        % handled after geometry resolution (no expansion)
    otherwise
        error('study_gamma_index_convergence:BadSensorMethod', ...
            'Unsupported sensor_placement_method "%s" (use determine_sensor_mask or full_plane_*).', ...
            CONFIG.sensor_placement_method);
end

% Resolve the final grid geometry (expansion offset + FFT-optimal padded dims).
[off, expdims, paddims] = resolve_grid_geometry([Nx, Ny, Nz], gp, ...
    CONFIG.use_grid_padding, CONFIG.pml_size);
if gp.expanded
    fprintf('  [SENSOR] Grid expanded to [%d %d %d] (offset [%d %d %d]).\n', ...
        expdims, off);
end
fprintf('  [SENSOR] Final k-Wave grid: [%d %d %d].\n', paddims);

% Build the sensor mask on the final padded grid.
sensor = struct();
sensor.mask = zeros(paddims);
switch CONFIG.sensor_placement_method
    case 'determine_sensor_mask'
        m = min(paddims, size(sensor_mask_orig));
        sensor.mask(1:m(1), 1:m(2), 1:m(3)) = double(sensor_mask_orig(1:m(1), 1:m(2), 1:m(3)));
    case 'full_plane_anterior'
        sensor.mask(off(1) + CONFIG.sensor_x_index, :, :) = 1;
    case 'full_plane_lateral'
        sensor.mask(:, off(2) + CONFIG.sensor_y_index, :) = 1;
end
numSensorPts = sum(sensor.mask(:));
fprintf('  [SENSOR] %d active points.\n', numSensorPts);
if numSensorPts == 0
    error('study_gamma_index_convergence:EmptySensor', 'Sensor mask is empty.');
end

%% ========================= EMBED ONTO FINAL GRID =======================
% Original-grid data (media, p0) sits at offset `off` in the padded grid; water
% (medium) / zero (p0, dose) fills the expansion + FFT padding. This reproduces
% the two-stage expand-then-pad flow of run_standalone_simulation exactly.

% Padded k-Wave media (density/sound/alpha in water fill). Gruneisen is NOT a
% k-Wave wave-propagation property, so it is carried as a SEPARATE padded array
% (needed only for p0 and the pressure->dose conversion), never inside kmedium.
[kmed1, grun1_pad] = embed_medium(medium1, paddims, off);
[kmed3, ~]         = embed_medium(medium3, paddims, off);   % CT_3 gruneisen unused (recon is CT_1)

% Padded initial-pressure sources (zero fill)
p0_A = embed_fill(p0_A_orig, paddims, off, 0);
p0_B = embed_fill(p0_B_orig, paddims, off, 0);

% Reconstruction (model) medium is CT_1 for BOTH studies. Crop to the expanded
% grid for the per-iteration pressure->dose conversion.
grun_recon_exp = grun1_pad(1:expdims(1),     1:expdims(2), 1:expdims(3));
dens_recon_exp = kmed1.density(1:expdims(1), 1:expdims(2), 1:expdims(3));

% Reference truth (D1) and the per-study low-dose masks on the expanded grid.
ref_dose_exp  = embed_fill(D1, expdims, off, 0);   % common reference for both studies
doseMaskA_exp = embed_mask(doseMaskA, expdims, off);
doseMaskB_exp = embed_mask(doseMaskB, expdims, off);

% Body mask (CT_1, the recon anatomy) for masking the reconstructed dose.
if isfield(sct1, 'bodyMask')
    body_exp = embed_mask(sct1.bodyMask, expdims, off);
else
    body_exp = true(expdims);
end

%% ========================= k-WAVE GRID (shared) ========================
% One kgrid for both studies. dt satisfies CFL for the faster of the two media;
% Nt covers the slower of the two, so the shared grid is stable for either
% forward/recon medium combination.

fprintf('\n[KGRID] Setting up shared k-Wave grid...\n');
kgrid = kWaveGrid(paddims(1), dx, paddims(2), dy, paddims(3), dz);

maxC = max(max(kmed1.sound_speed(:)), max(kmed3.sound_speed(:)));
minC = min(min(kmed1.sound_speed(kmed1.sound_speed > 0)), ...
           min(kmed3.sound_speed(kmed3.sound_speed > 0)));
dt   = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;

gridDiag = sqrt((paddims(1)*dx)^2 + (paddims(2)*dy)^2 + (paddims(3)*dz)^2);
simTime  = 2.5 * gridDiag / minC;
Nt       = ceil(simTime / dt);

% Air (~343 m/s) sets minC and inflates Nt; shorten the recording when it does.
if isfield(CONFIG, 'Nt_scaling') && CONFIG.Nt_scaling ~= 0
    air_c = 343;
    if isfield(CONFIG.tissue_tables, 'threshold_2') && ...
            isfield(CONFIG.tissue_tables.threshold_2, 'air')
        air_c = CONFIG.tissue_tables.threshold_2.air.sound_speed;
    end
    if minC <= air_c
        Nt = max(1, ceil(Nt / CONFIG.Nt_scaling));
    end
end

kgrid.dt = dt;
kgrid.Nt = Nt;
fprintf('  dt = %.2e s, Nt = %d, maxC = %.0f m/s, minC = %.0f m/s\n', dt, Nt, maxC, minC);

if CONFIG.use_gpu
    try
        gpuDevice;
        dataCast = 'gpuArray-single';
        fprintf('  Compute: GPU\n');
    catch
        dataCast = 'single';
        fprintf('  Compute: CPU (GPU unavailable)\n');
    end
else
    dataCast = 'single';
    fprintf('  Compute: CPU\n');
end
inputArgs = {'Smooth', false, 'PMLInside', false, 'PMLSize', CONFIG.pml_size, ...
             'DataCast', dataCast, 'PlotSim', false};

%% ========================= STUDY A: matched CT_1 =======================

fprintf('\n=========================================================\n');
fprintf('  STUDY A - matched CT_1 (fwd CT_1, recon CT_1, gamma vs D1)\n');
fprintf('=========================================================\n');
[pass_A, snap_A, noise_amp_A] = run_tr_convergence('A/CT1', kgrid, inputArgs, ...
    kmed1, kmed1, p0_A, sensor, ...
    grun_recon_exp, dens_recon_exp, ref_dose_exp, doseMaskA_exp, body_exp, ...
    num_pulses, expdims, spacing_mm, CONFIG.gamma_n, CONFIG.num_time_reversal_iter, CONFIG, []);

%% ========================= STUDY B: cross-CT ===========================

fprintf('\n=========================================================\n');
fprintf('  STUDY B - cross-CT (fwd CT_3, recon CT_1, gamma vs D1)\n');
fprintf('=========================================================\n');
[pass_B, snap_B] = run_tr_convergence('B/CT3->CT1', kgrid, inputArgs, ...
    kmed3, kmed1, p0_B, sensor, ...
    grun_recon_exp, dens_recon_exp, ref_dose_exp, doseMaskB_exp, body_exp, ...
    num_pulses, expdims, spacing_mm, CONFIG.gamma_n, CONFIG.num_time_reversal_iter, CONFIG, []);

%% ========================= STUDY N: noise-floor null ===================
% Null hypothesis: reconstruct from PURE ELECTRONIC NOISE (no coherent signal),
% injected at the same absolute amplitude Study A saw, through the identical
% deconvolution + TR-convergence + gamma pipeline (gamma vs the same truth D1
% over the same 10% dose region). The A/B signal curves should sit far above
% this floor. Only meaningful when the pulse model actually injects noise; when
% CONFIG.null_realizations > 1 the pass rates are averaged over independent
% noise draws (the snapshot panel uses the first realization).

run_null = CONFIG.convolution_kernel > 0 && CONFIG.conv_noise_level > 0 && ...
           isfinite(noise_amp_A) && noise_amp_A > 0;
pass_N = [];
snap_N = [];
if run_null
    n_real = max(1, round(CONFIG.null_realizations));
    fprintf('\n=========================================================\n');
    fprintf('  STUDY N - noise floor (noise-only, recon CT_1, gamma vs D1)\n');
    fprintf('  noise amp = %.3e Pa (matched to Study A), %d realization(s)\n', ...
        noise_amp_A, n_real);
    fprintf('=========================================================\n');
    opts_null  = struct('noise_only', true, 'forced_noise_amp', noise_amp_A);
    passN_acc  = zeros(CONFIG.num_time_reversal_iter, 1);
    for r = 1:n_real
        [pn, sn] = run_tr_convergence(sprintf('N/noise r%d', r), kgrid, inputArgs, ...
            kmed1, kmed1, p0_A, sensor, ...
            grun_recon_exp, dens_recon_exp, ref_dose_exp, doseMaskA_exp, body_exp, ...
            num_pulses, expdims, spacing_mm, CONFIG.gamma_n, ...
            CONFIG.num_time_reversal_iter, CONFIG, opts_null);
        passN_acc = passN_acc + pn(:);
        if r == 1, snap_N = sn; end
    end
    pass_N = passN_acc / n_real;
else
    fprintf(['\n[STUDY N] Skipped noise-floor null: requires a pulse model with ' ...
             'noise (CONFIG.convolution_kernel > 0 and conv_noise_level > 0).\n']);
end

%% ========================= CONSOLE SUMMARY =============================

N_iter   = CONFIG.num_time_reversal_iter;
have_null = ~isempty(pass_N);
fprintf('\n==================== GAMMA %g%%/%g mm CONVERGENCE ====================\n', ...
    CONFIG.gamma_n, CONFIG.gamma_n);
if have_null
    fprintf('  %-6s  %-14s  %-14s  %-16s\n', 'iter', 'Study A (%)', 'Study B (%)', 'Noise floor (%)');
    for it = 1:N_iter
        fprintf('  %-6d  %-14.2f  %-14.2f  %-16.2f\n', it, pass_A(it), pass_B(it), pass_N(it));
    end
else
    fprintf('  %-6s  %-14s  %-14s\n', 'iter', 'Study A (%)', 'Study B (%)');
    for it = 1:N_iter
        fprintf('  %-6d  %-14.2f  %-14.2f\n', it, pass_A(it), pass_B(it));
    end
end
fprintf('====================================================================\n');

%% ========================= PLOTS =======================================

green = [0.15, 0.60, 0.20];
red   = [0.80, 0.15, 0.15];
gray  = [0.45, 0.45, 0.45];   % noise floor (null)
null_lbl = 'Noise floor (null)';

if CONFIG.plot_results
    ptxt = strrep(CONFIG.patient_id, '_', '\_');
    stxt = strrep(CONFIG.session,    '_', '\_');
    ftxt = strrep(file_tag,          '_', '\_');

    % --- Truth-vs-recon dose + gamma-index comparison panels (iters 1 & 10) ---
    % Study A (matched CT_1), Study B (cross-CT), and the noise-floor null.
    fhA = plot_dose_gamma_panel(snap_A, ...
        sprintf('Study A (matched CT\\_1)  |  %s', ftxt), ptxt, stxt);
    save_and_maybe_close(fhA, CONFIG, [file_tag '_panel_studyA']);

    fhB = plot_dose_gamma_panel(snap_B, ...
        sprintf('Study B (CT\\_3\\rightarrowCT\\_1)  |  %s', ftxt), ptxt, stxt);
    save_and_maybe_close(fhB, CONFIG, [file_tag '_panel_studyB']);

    if have_null
        fhN = plot_dose_gamma_panel(snap_N, ...
            sprintf('Noise floor / null (noise-only)  |  %s', ftxt), ptxt, stxt);
        save_and_maybe_close(fhN, CONFIG, [file_tag '_panel_noise']);
    end

    % --- Gamma pass rate vs TR iteration ---
    % The noise floor (gray) is overlaid on every convergence graph so the A/B
    % signal curves can be seen sitting above it.

    % Figure - Study A alone (+ noise floor).
    sA = {pass_A}; lA = {'Recon CT\_1 vs Truth D1'}; cA = {green};
    if have_null, sA{end+1} = pass_N; lA{end+1} = null_lbl; cA{end+1} = gray; end
    fh1 = plot_gamma_convergence(sA, lA, cA, CONFIG.gamma_n, ...
        sprintf('Study A: Matched CT\\_1 Reconstruction  |  %s', ftxt), ...
        CONFIG.patient_id, CONFIG.session);
    save_and_maybe_close(fh1, CONFIG, [file_tag '_convergence_studyA']);

    % Figure - Study B alone (+ noise floor).
    sB = {pass_B}; lB = {'Recon CT\_3\rightarrowCT\_1 vs Truth D1'}; cB = {red};
    if have_null, sB{end+1} = pass_N; lB{end+1} = null_lbl; cB{end+1} = gray; end
    fh2 = plot_gamma_convergence(sB, lB, cB, CONFIG.gamma_n, ...
        sprintf('Study B: Cross-CT (CT\\_3 delivery, CT\\_1 model)  |  %s', ftxt), ...
        CONFIG.patient_id, CONFIG.session);
    save_and_maybe_close(fh2, CONFIG, [file_tag '_convergence_studyB']);

    % Figure - both studies overlaid (+ noise floor).
    s3 = {pass_A, pass_B};
    l3 = {'A: Recon CT\_1 vs D1', 'B: Recon CT\_3\rightarrowCT\_1 vs D1'};
    c3 = {green, red};
    if have_null, s3{end+1} = pass_N; l3{end+1} = null_lbl; c3{end+1} = gray; end
    fh3 = plot_gamma_convergence(s3, l3, c3, CONFIG.gamma_n, ...
        sprintf('Gamma Convergence: Signal vs Noise Floor  |  %s', ftxt), ...
        CONFIG.patient_id, CONFIG.session);
    save_and_maybe_close(fh3, CONFIG, [file_tag '_convergence_overlay']);

    % Figure - noise floor alone (its own convergence, mirroring A/B).
    if have_null
        fhN2 = plot_gamma_convergence({pass_N}, {null_lbl}, {gray}, CONFIG.gamma_n, ...
            sprintf('Noise Floor (null): Noise-only Reconstruction  |  %s', ftxt), ...
            CONFIG.patient_id, CONFIG.session);
        save_and_maybe_close(fhN2, CONFIG, [file_tag '_convergence_noise']);
    end
end

%% ========================= SAVE RESULTS ================================

total_runtime = toc(study_timer);

% Always assemble the RESULTS output (also written to disk when save_results).
RESULTS = struct();
RESULTS.config          = CONFIG;
RESULTS.file_tag        = file_tag;
RESULTS.gamma_crit      = CONFIG.gamma_n;
RESULTS.num_iterations  = N_iter;
RESULTS.iterations      = (1:N_iter)';
RESULTS.pass_rate_A     = pass_A;   % matched CT_1, recon vs D1
RESULTS.pass_rate_B     = pass_B;   % cross-CT (CT_3 fwd, CT_1 recon), recon vs D1
RESULTS.pass_rate_N     = pass_N;   % noise-floor null (noise-only), recon vs D1; [] if skipped
RESULTS.snapshot_A      = snap_A;   % iter-1/10 recon dose + gamma map (Study A)
RESULTS.snapshot_B      = snap_B;   % iter-1/10 recon dose + gamma map (Study B)
RESULTS.snapshot_N      = snap_N;   % iter-1/10 recon dose + gamma map (noise floor); [] if skipped
RESULTS.noise_amp_A     = noise_amp_A;   % absolute noise amplitude (Pa) used for the null
RESULTS.dose_filename_ct1 = CONFIG.dose_filename_ct1;
RESULTS.dose_filename_ct3 = CONFIG.dose_filename_ct3;
RESULTS.spacing_mm      = spacing_mm;
RESULTS.grid_size       = expdims;
RESULTS.total_runtime_s = total_runtime;

if CONFIG.save_results
    if ~isempty(CONFIG.results_dir) && ~isfolder(CONFIG.results_dir)
        mkdir(CONFIG.results_dir);
    end
    save(CONFIG.output_file, '-struct', 'RESULTS', '-v7.3');
    fprintf('\nResults saved to: %s\n', CONFIG.output_file);
end

fprintf('\nTotal runtime: %.1f s (%.2f min).\n', total_runtime, total_runtime / 60);
fprintf('Gamma index convergence study complete.\n');

end   % ===== end of main function study_gamma_index_convergence =====


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function [pass_rates, snap, noise_amp] = run_tr_convergence(tag, kgrid, inputArgs, ...
        kmed_fwd, kmed_recon, p0, sensor, ...
        grun_recon_exp, dens_recon_exp, ref_dose_exp, doseMask_exp, body_exp, ...
        num_pulses, expdims, spacing_mm, gamma_crit, N_iter, CONFIG, opts)
%RUN_TR_CONVERGENCE Forward sim (through kmed_fwd) + iterative time reversal
%  (through kmed_recon), recording the gamma pass rate of the reconstructed dose
%  against ref_dose_exp AFTER EVERY iteration.
%
%  The measurement is simulated once through the TRUE medium (kmed_fwd); every
%  back-propagation and residual re-simulation uses the RECON (model) medium
%  (kmed_recon). For Study A both are CT_1; for Study B kmed_fwd is CT_3 and
%  kmed_recon is CT_1 (model/anatomy mismatch).
%
%  Returns:
%    pass_rates - [N_iter x 1], gamma %/mm pass rate per TR iteration.
%    snap       - snapshot struct for the dose/gamma comparison panel, with the
%                 reconstructed dose and gamma-index map captured at the iterations
%                 in CONFIG.snapshot_iters (clamped to [1, N_iter]). Fields:
%                   .iters      [1 x M] captured iteration numbers
%                   .recon      {1 x M} recon dose volumes (expdims)
%                   .gmap       {1 x M} gamma-index maps (expdims; [] on failure)
%                   .pass       [1 x M] pass rate at each captured iteration
%                   .ref_dose   truth D1 on the expanded grid
%                   .slice      axial (Z) index of the truth max, for display
%                   .spacing_mm, .gamma_crit
%    noise_amp  - absolute electronic-noise amplitude (Pa) injected by the
%                 measurement chain (0 when no pulse model). Feeds the null study.
%
%  NULL / NOISE-FLOOR MODE (opts.noise_only == true): no coherent forward source
%  is simulated. Instead the sensor traces are pure electronic noise at
%  opts.forced_noise_amp (matched to a real study), pushed through the identical
%  deconvolution + TR-convergence + per-iteration gamma pipeline. The resulting
%  pass_rates are the noise floor: real studies should sit far above them.

    if nargin < 19 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'noise_only'),       opts.noise_only = false;      end
    if ~isfield(opts, 'forced_noise_amp'), opts.forced_noise_amp = [];   end

    ex = expdims(1); ey = expdims(2); ez = expdims(3);

    % Snapshot iterations for the comparison panel (default [1 10]).
    snap_iters = CONFIG.snapshot_iters(:)';
    snap_iters = unique(min(max(round(snap_iters), 1), N_iter));

    snap = struct();
    snap.iters      = [];
    snap.recon      = {};
    snap.gmap       = {};
    snap.pass       = [];
    snap.ref_dose   = ref_dose_exp;
    snap.spacing_mm = spacing_mm;
    snap.gamma_crit = gamma_crit;
    [~, truth_max_idx] = max(ref_dose_exp(:));
    [~, ~, snap.slice] = ind2sub([ex, ey, ez], truth_max_idx);

    %% ---- Forward simulation (measurement through the true medium) ----
    % GPU note: with DataCast='gpuArray-single', kspaceFirstOrder3D returns a
    % gpuArray and clears the GPU device on exit. Retaining that gpuArray and
    % re-feeding it into the next call invalidates it ("The data no longer exists
    % on the device"). So EVERY kspaceFirstOrder3D output is gathered to CPU here
    % and all persistent loop state (sensorData, reconPressure) is kept on CPU;
    % k-Wave re-casts the CPU inputs to GPU internally each call.
    FS = 1 / kgrid.dt;
    if opts.noise_only
        % ---- NULL: no coherent source; sensor sees only electronic noise ----
        % Build zero-signal traces of the SAME shape a real forward sim would
        % produce, then let the measurement chain inject noise at the matched
        % absolute amplitude and Wiener-deconvolve it exactly as for real data.
        nSens = nnz(sensor.mask);
        Nt    = kgrid.Nt;
        fprintf('[%s] Noise-only null: %d sensors x %d samples, noise amp = %.3e Pa.\n', ...
            tag, nSens, Nt, opts.forced_noise_amp);
        sensorData = zeros(nSens, Nt);
        [sensorData, sensorData_measured, noise_amp] = apply_measurement_chain( ...
            sensorData, kgrid.dt, FS, CONFIG, true, opts.forced_noise_amp);
    else
        fprintf('[%s] Forward simulation...\n', tag);
        source_fwd = struct('p0', p0);
        fwd_tic    = tic;
        sensorData = gather(kspaceFirstOrder3D(kgrid, kmed_fwd, source_fwd, sensor, inputArgs{:}));
        fprintf('[%s] Forward complete (%.1f s). Sensor data: [%d x %d]\n', ...
            tag, toc(fwd_tic), size(sensorData, 1), size(sensorData, 2));

        sensorData = gather(smooth(sensorData));

        %% ---- Physical measurement chain (pulse / freq response / noise / deconv) ----
        [sensorData, sensorData_measured, noise_amp] = apply_measurement_chain( ...
            sensorData, kgrid.dt, FS, CONFIG);
    end

    %% ---- Pressure->dose conversion factor for the recon medium ----
    conversionFactor = grun_recon_exp .* dens_recon_exp;
    conversionFactor(conversionFactor == 0) = 1;
    air_recon_exp = dens_recon_exp < 100;   % air voxels carry no real PA dose signal

    mask_to_region = ~(isfield(CONFIG, 'mask_recon_to_dose_region') && ...
                       ~CONFIG.mask_recon_to_dose_region);

    pass_rates = nan(N_iter, 1);

    %% ---- Time-reversal loop with per-iteration gamma ----
    tr_tic = tic;
    for it = 1:N_iter
        % Back-propagate the (residual-corrected) sensor data through the model medium.
        source_tr        = struct();
        source_tr.p_mask = sensor.mask;
        source_tr.p      = fliplr(sensorData);
        source_tr.p_mode = 'dirichlet';

        sensor_tr        = struct();
        sensor_tr.mask   = ones(size(p0));
        sensor_tr.record = {'p_final'};

        p0_recon = kspaceFirstOrder3D(kgrid, kmed_recon, source_tr, sensor_tr, inputArgs{:});
        if isstruct(p0_recon) && isfield(p0_recon, 'p_final')
            reconPressure = reshape(gather(p0_recon.p_final), size(p0));
        else
            reconPressure = reshape(gather(p0_recon), size(p0));
        end
        reconPressure = max(reconPressure, 0);   % CPU

        % ---- Convert THIS iteration's estimate to dose, then gamma vs truth ----
        rp = reconPressure(1:ex, 1:ey, 1:ez) * CONFIG.correction_factor;
        recon_dose = (rp ./ conversionFactor) * num_pulses;
        if mask_to_region
            recon_dose = recon_dose .* double(doseMask_exp);
        end
        recon_dose = recon_dose .* double(body_exp);
        recon_dose(air_recon_exp) = 0;

        % LSQ normalization of the recon to the truth (from study_pass_rates_allsegments).
        if isfield(CONFIG, 'normalize') && CONFIG.normalize
            recon_dose = recon_dose * least_squares_gain(ref_dose_exp, recon_dose);
        end

        want_map = ismember(it, snap_iters);
        [pass_rates(it), gmap_it] = gamma_pass_rate(ref_dose_exp, recon_dose, ...
            gamma_crit, spacing_mm, want_map);
        fprintf('[%s] iter %2d/%d:  gamma %g%%/%g mm pass = %.2f%%  (max recon p = %.3e)\n', ...
            tag, it, N_iter, gamma_crit, gamma_crit, pass_rates(it), max(rp(:)));

        % Capture the recon dose + gamma map at the requested snapshot iterations.
        if want_map
            snap.iters(end + 1) = it;
            snap.recon{end + 1} = recon_dose;
            snap.gmap{end + 1}  = gmap_it;
            snap.pass(end + 1)  = pass_rates(it);
        end

        % ---- Residual correction for the next iteration (through the model medium) ----
        % reconPressure is CPU; k-Wave re-casts it to GPU internally. Gather the
        % result so sensorData stays a CPU array across iterations.
        if it < N_iter
            source_resid    = struct('p0', reconPressure);
            sensorDataRecon = gather(kspaceFirstOrder3D(kgrid, kmed_recon, source_resid, sensor, inputArgs{:}));
            sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
        end
    end
    fprintf('[%s] TR convergence complete (%.1f s, %d iterations).\n', tag, toc(tr_tic), N_iter);
end


function [sensorData, sensorData_measured, noise_amp] = apply_measurement_chain( ...
        sensorData, dt, FS, CONFIG, noise_only, forced_noise_amp)
%APPLY_MEASUREMENT_CHAIN Model the physical acquisition chain, in order:
%  (1) convolve each sensor trace with a Gaussian pulse kernel, (2) apply the
%  sensor frequency response, (3) add electronic noise, (4) Wiener-deconvolve the
%  pulse kernel. With CONFIG.convolution_kernel == 0 only the frequency response
%  is applied. Ported from run_standalone_simulation.m.
%
%  noise_only (default false): the NULL path. The coherent (post-frequency-
%  response) signal is zeroed, so the only content is freshly drawn electronic
%  noise; everything downstream (noise draw + Wiener deconvolution) is identical
%  to the real path, so the noise "appears in the same way".
%  forced_noise_amp (default []): use this absolute noise amplitude (Pa) instead
%  of the signal-scaled default -- lets the null match a real study's noise level.
%
%  Third output noise_amp is the absolute electronic-noise amplitude used (Pa),
%  0 when there is no pulse model (no noise-injection point).
    if nargin < 5 || isempty(noise_only), noise_only = false; end
    if nargin < 6,                        forced_noise_amp = []; end
    noise_amp = 0;

    if CONFIG.convolution_kernel > 0
        conv_kernel_sigma  = CONFIG.convolution_kernel;
        conv_noise_level   = CONFIG.conv_noise_level;
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

        sensorData_conv   = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));
        sensorData_resp   = gaussianFilter(sensorData_conv, FS, 0.35e6, 100, true);

        if noise_only
            % Null: strip the coherent signal, keep only the noise draw.
            sensorData_resp = zeros(size(sensorData_resp));
        end

        if ~isempty(forced_noise_amp)
            noise_amp = forced_noise_amp;
        else
            noise_amp = conv_noise_level * max(abs(sensorData_resp(:)));
        end
        sensorData_noisy  = sensorData_resp + noise_amp * randn(size(sensorData_resp));
        sensorData_deconv = real(ifft( ...
            fft(sensorData_noisy, [], 2) .* H_conj ./ (H_power + conv_deconv_lambda), [], 2));

        sensorData          = single(sensorData_deconv);
        sensorData_measured = single(sensorData_deconv);
    else
        % No pulse model => no electronic-noise injection point. In noise-only
        % mode there is nothing to reconstruct; the caller guards against this.
        if noise_only
            sensorData = zeros(size(sensorData), 'single');
        else
            % gather so the returned sensor data is CPU-resident (gaussianFilter
            % preserves the input's device; the caller keeps sensorData on the CPU).
            sensorData = gather(gaussianFilter(sensorData, FS, 0.35e6, 100, true));
        end
        sensorData_measured = sensorData;
    end
end


function [pr, gmap_out] = gamma_pass_rate(ref, tgt, crit, spacing, want_map)
%GAMMA_PASS_RATE Global gamma pass rate (%) of tgt against ref at crit%/crit mm.
%  Reference = ref (truth); the 10% reference cutoff sets the eval mask. CalcGamma
%  console output is suppressed via evalc. Returns NaN on failure / empty ref.
%
%  Second output gmap_out is the full 3D gamma-index map (same size as ref); it is
%  only returned when want_map is true (default false), otherwise []. Voxels below
%  the 10% eval cutoff are NaN'd in gmap_out so they don't skew the display scale.
    if nargin < 5, want_map = false; end
    gmap_out = [];
    ref = double(ref);
    tgt = double(tgt);
    if max(ref(:)) <= 0
        pr = NaN;
        return;
    end
    mask = ref >= 0.10 * max(ref(:));
    ref_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', ref);   %#ok<NASGU>
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing, 'data', tgt);   %#ok<NASGU>
    gmap = [];
    try
        evalc(['gmap = CalcGamma(ref_struct, tgt_struct, crit, crit, ', ...
               '''local'', 0, ''limit'', crit*2, ''restrict'', 1);']);
        pr = 100 * mean(gmap(mask) <= 1);
        if want_map
            gmap_out = gmap;
            gmap_out(~mask) = NaN;   % show gamma only in the evaluated region
        end
    catch
        pr = NaN;
    end
end


function g = least_squares_gain(rs_truth, recon)
%LEAST_SQUARES_GAIN Scalar gain aligning a recon to its truth (relative norm).
%  g = sum(rs.*recon)/sum(recon.^2) over the truth's 10% region -- the g
%  minimizing ||rs - g*recon||^2 there. Falls back to 1 for empty/zero inputs.
%  (Ported from study_pass_rates_allsegments.m.)
    rs_truth = double(rs_truth);
    recon    = double(recon);
    if max(rs_truth(:)) > 0
        mask = rs_truth >= 0.10 * max(rs_truth(:));
    else
        mask = true(size(rs_truth));
    end
    r = recon(mask);
    denom = sum(r .^ 2);
    if denom > 0
        g = sum(rs_truth(mask) .* r) / denom;
    else
        g = 1;
    end
end


function medium = build_medium(sct, CONFIG)
%BUILD_MEDIUM Acoustic medium from a CBCT struct + the standalone uniform/bath
%  overrides. Mirrors the medium setup in run_standalone_simulation.m. Guarantees
%  a full 3D alpha_coeff so the padding/embedding helpers work uniformly.
    medium = create_acoustic_medium(sct, CONFIG);
    gs = medium.grid_size;

    if CONFIG.force_uniform_density
        medium.density = ones(gs) * CONFIG.uniform_density;
    end
    if CONFIG.force_uniform_sound_speed
        medium.sound_speed = ones(gs) * CONFIG.uniform_sound_speed;
    end
    if CONFIG.force_uniform_attenuation
        medium.alpha_coeff = ones(gs) * CONFIG.uniform_alpha_coeff;
        medium.alpha_power = 1.1;
    end
    if CONFIG.force_uniform_gruneisen
        medium.gruneisen = ones(gs) * CONFIG.uniform_gruneisen;
    end

    % Force the coupling bath (outside body / on couch) to the uniform medium.
    if isfield(sct, 'bodyMask')
        outside_body = ~logical(sct.bodyMask);
        if isfield(sct, 'couchMask')
            outside_body = outside_body | logical(sct.couchMask);
        end
        medium.density(outside_body)     = CONFIG.uniform_density;
        medium.sound_speed(outside_body) = CONFIG.uniform_sound_speed;
        medium.alpha_coeff(outside_body) = CONFIG.uniform_alpha_coeff;
        medium.gruneisen(outside_body)   = CONFIG.uniform_gruneisen;
    end

    if numel(medium.alpha_coeff) == 1
        medium.alpha_coeff = medium.alpha_coeff * ones(gs);
    end
end


function p0 = compute_p0(doseGrid, medium, num_pulses)
%COMPUTE_P0 Initial pressure: p0 = (dose/num_pulses) .* gruneisen .* density,
%  k-Wave-smoothed. Mirrors run_standalone_simulation.m.
    dose_per_pulse = doseGrid / num_pulses;
    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    p0 = smooth(p0);
end


function [off, expdims, paddims] = resolve_grid_geometry(orig_dims, gp, use_pad, pml)
%RESOLVE_GRID_GEOMETRY Final grid geometry from a determine_sensor_mask grid_pad.
%  Because expansion pre-pads with water and FFT padding post-pads with water,
%  original-grid data lives at offset `off` within the FFT-optimal `paddims`; the
%  intermediate `expdims` (expansion only) is where the recon dose is compared.
%    off     = [dim1_pre, dim2_pre, dim3_pre] (0 when not expanded)
%    expdims = orig_dims + pre + post
%    paddims = FFT-optimal( expdims )   (== expdims when use_pad is false)
%  Mapping (per run_standalone_simulation): dim1<-y, dim2<-x, dim3<-z.
    if isstruct(gp) && isfield(gp, 'expanded') && gp.expanded
        off  = [gp.y_pre, gp.x_pre, gp.z_pre];
        post = [gp.y_post, gp.x_post, gp.z_post];
    else
        off  = [0, 0, 0];
        post = [0, 0, 0];
    end
    expdims = orig_dims + off + post;
    if use_pad
        paddims = [find_optimal_kwave_size(expdims(1), pml), ...
                   find_optimal_kwave_size(expdims(2), pml), ...
                   find_optimal_kwave_size(expdims(3), pml)];
    else
        paddims = expdims;
    end
end


function [kmed, gruneisen_pad] = embed_medium(medium, dims, off)
%EMBED_MEDIUM Place a medium's arrays into the padded grid with water fill.
%  Returns:
%    kmed          - the k-Wave medium struct (density 1000, sound 1540, alpha 0
%                    in the fill; alpha_power = 1.1). Contains ONLY the fields
%                    kspaceFirstOrder3D uses -- NO gruneisen.
%    gruneisen_pad - the padded Gruneisen array (0 in the fill), carried
%                    separately because Gruneisen is a photoacoustic source /
%                    dose-conversion property, not a wave-propagation property.
    kmed = struct();
    kmed.density     = embed_fill(medium.density,     dims, off, 1000);
    kmed.sound_speed = embed_fill(medium.sound_speed, dims, off, 1540);
    kmed.alpha_coeff = embed_fill(medium.alpha_coeff, dims, off, 0);
    kmed.alpha_power = 1.1;

    gruneisen_pad = embed_fill(medium.gruneisen, dims, off, 0);
end


function out = embed_fill(src, dims, off, fillval)
%EMBED_FILL Place src into a dims-sized array of fillval at 1-based offset `off`.
    out = fillval * ones(dims);
    out(off(1) + (1:size(src, 1)), ...
        off(2) + (1:size(src, 2)), ...
        off(3) + (1:size(src, 3))) = src;
end


function out = embed_mask(src, dims, off)
%EMBED_MASK Place a logical src into a dims-sized false array at offset `off`.
    out = false(dims);
    out(off(1) + (1:size(src, 1)), ...
        off(2) + (1:size(src, 2)), ...
        off(3) + (1:size(src, 3))) = logical(src);
end


function [doseGrid, spacing_mm, meterset, gantry_angle] = load_field_dose(filepath, CONFIG)
%LOAD_FIELD_DOSE Load a step15 field/total dose .mat into a 3D array (+metadata).
%  Auto-detects the step15 formats (sparse/dense field_dose, sparse/dense total).
%  Ported from run_standalone_simulation.m. spacing/meterset/gantry are [] when
%  the file does not carry them.
    spacing_mm   = [];
    meterset     = [];
    gantry_angle = [];

    if ~isfile(filepath)
        error('load_field_dose:NotFound', 'Dose file not found: %s', filepath);
    end
    dose_data   = load(filepath);
    dose_fields = fieldnames(dose_data);

    if isfield(dose_data, 'field_dose')
        fd = dose_data.field_dose;
        if ~isfield(fd, 'dose_Gy')
            error('load_field_dose:NoDose', 'field_dose struct missing dose_Gy field.');
        end
        if (isfield(fd, 'is_sparse') && fd.is_sparse) || issparse(fd.dose_Gy)
            if ~isfield(fd, 'dose_dims')
                error('load_field_dose:NoDims', 'field_dose.dose_dims missing (sparse dose).');
            end
            doseGrid = reshape(full(fd.dose_Gy), fd.dose_dims);
        else
            doseGrid = double(fd.dose_Gy);
        end
        if isfield(fd, 'spacing') && ~isempty(fd.spacing),   spacing_mm   = fd.spacing(:)';   end
        if isfield(fd, 'meterset') && ~isempty(fd.meterset), meterset     = fd.meterset;      end
        if isfield(fd, 'gantry_angle'),                      gantry_angle = fd.gantry_angle;  end
    elseif isfield(dose_data, 'total_rs_dose_sparse')
        if ~isfield(dose_data, 'total_rs_dose_dims')
            error('load_field_dose:NoDims', 'total_rs_dose_dims missing (sparse total dose).');
        end
        doseGrid = reshape(full(dose_data.total_rs_dose_sparse), dose_data.total_rs_dose_dims);
    elseif isfield(dose_data, 'total_rs_dose')
        doseGrid = dose_data.total_rs_dose;
    elseif isfield(dose_data, 'dose_Gy')
        doseGrid = dose_data.dose_Gy;
    elseif numel(dose_fields) == 1
        doseGrid = dose_data.(dose_fields{1});
    else
        error('load_field_dose:AutoDetect', ...
            'Cannot auto-detect dose variable. Fields: %s', strjoin(dose_fields, ', '));
    end

    doseGrid = double(doseGrid);
    if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
        error('load_field_dose:Not3D', 'Dose data must be a 3D numeric array: %s', filepath);
    end
end


function sct = load_cbct(filepath)
%LOAD_CBCT Load a CBCT{1,3}_resampled struct (cubeHU + spacing/bodyMask/...).
    if ~isfile(filepath)
        error('load_cbct:NotFound', 'CBCT file not found: %s', filepath);
    end
    cbct_data = load(filepath);
    if isfield(cbct_data, 'CBCT1_resampled')
        sct = cbct_data.CBCT1_resampled;
    elseif isfield(cbct_data, 'CBCT3_resampled')
        sct = cbct_data.CBCT3_resampled;
    else
        error('load_cbct:NoVar', 'CBCT1_resampled / CBCT3_resampled not found in %s', filepath);
    end
    if ~isfield(sct, 'cubeHU')
        error('load_cbct:NoCube', 'CBCT resampled struct missing required field: cubeHU (%s)', filepath);
    end
end


function assert_same_grid(name, ref_dims, got_dims)
%ASSERT_SAME_GRID Error unless got_dims matches ref_dims (shared dose grid check).
    if ~isequal(ref_dims, got_dims)
        error('study_gamma_index_convergence:GridMismatch', ...
            '%s grid [%s] does not match the reference dose grid [%s]. All four volumes must share one grid.', ...
            name, num2str(got_dims), num2str(ref_dims));
    end
end


function fh = plot_gamma_convergence(pass_cell, labels, colors, gamma_crit, ttl, patient_id, session)
%PLOT_GAMMA_CONVERGENCE Gamma pass rate (%) vs TR iteration for one or more series.
%  pass_cell - cell of [N x 1] pass-rate vectors; labels/colors index-matched.
%  Returns the figure handle fh (for saving by the caller).
    N = numel(pass_cell{1});
    fh = figure('Name', ttl, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [120, 120, 780, 520]);
    hold on;
    markers = {'-o', '-s', '-^', '-d'};
    h = gobjects(1, numel(pass_cell));
    for i = 1:numel(pass_cell)
        m = markers{mod(i - 1, numel(markers)) + 1};
        h(i) = plot(1:N, pass_cell{i}(:)', m, 'Color', colors{i}, ...
            'MarkerFaceColor', colors{i}, 'MarkerSize', 6, 'LineWidth', 1.8, ...
            'DisplayName', labels{i});
    end
    yline(90, 'k--', '90%', 'LineWidth', 1.0, 'FontSize', 8, ...
        'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    hold off;
    grid on; box on;
    xlim([0.5, N + 0.5]);
    ylim([0, 101]);
    xticks(1:N);
    xlabel('Time-reversal iteration');
    ylabel(sprintf('Gamma pass rate (%%)  @ %g%%/%g mm', gamma_crit, gamma_crit));
    title(sprintf('%s   |   %s / %s', ttl, ...
        strrep(patient_id, '_', '\_'), strrep(session, '_', '\_')), ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'tex');
    legend(h, 'Location', 'southeast', 'FontSize', 9, 'Interpreter', 'tex');
    drawnow;
end


function fh = plot_dose_gamma_panel(snap, study_label, patient_txt, session_txt)
%PLOT_DOSE_GAMMA_PANEL Truth-vs-recon dose + gamma-index comparison panel.
%  Renders the axial slice through the truth-dose maximum. Layout (2 x 3):
%     [ Truth D1      | Recon iter s1  | Recon iter s2 ]
%     [ (blank)       | Gamma iter s1  | Gamma iter s2 ]
%  Top row: dose in 'hot' with a shared color scale (0..max truth). Bottom row:
%  gamma index (green<=1 pass, red>1 fail) under the matching recon column, with
%  the pass rate in the title. Returns the figure handle fh.
    sl    = snap.slice;
    truth = snap.ref_dose;
    tmax  = max(truth(:));
    if ~(tmax > 0), tmax = 1; end

    nsnap = min(2, numel(snap.iters));   % panel has room for two snapshot columns
    gcmap = gamma_colormap();

    fh = figure('Name', sprintf('Dose/Gamma panel  %s', study_label), ...
        'Color', 'w', 'NumberTitle', 'off', 'Position', [80, 80, 1220, 760]);
    tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % (1,1) truth dose
    ax = nexttile(tl, 1);
    imagesc(ax, squeeze(truth(:, :, sl))');
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, 'hot');
    caxis(ax, [0, tmax]); colorbar(ax);
    xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
    title(ax, sprintf('Truth D1   (Z=%d)', sl), 'FontWeight', 'bold');

    % (1,2..3) recon dose at each snapshot iteration (shared dose scale)
    for k = 1:nsnap
        ax = nexttile(tl, 1 + k);
        imagesc(ax, squeeze(snap.recon{k}(:, :, sl))');
        axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, 'hot');
        caxis(ax, [0, tmax]); colorbar(ax);
        xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
        title(ax, sprintf('Recon dose   (iter %d)', snap.iters(k)), 'FontWeight', 'bold');
    end

    % (2,1) intentionally blank so gamma maps sit under their recon columns
    ax = nexttile(tl, 4); axis(ax, 'off');

    % (2,2..3) gamma-index map at each snapshot iteration
    for k = 1:nsnap
        ax = nexttile(tl, 4 + k);
        gmap = snap.gmap{k};
        if isempty(gmap)
            axis(ax, 'off');
            title(ax, sprintf('Gamma iter %d (unavailable)', snap.iters(k)));
            continue;
        end
        gslice = squeeze(gmap(:, :, sl))';
        him = imagesc(ax, gslice);
        set(him, 'AlphaData', ~isnan(gslice));   % transparent outside eval region
        axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gcmap);
        caxis(ax, [0, 2]); cb = colorbar(ax); cb.Label.String = '\gamma index';
        xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
        title(ax, sprintf('Gamma %g%%/%g mm  iter %d   (pass %.1f%%)', ...
            snap.gamma_crit, snap.gamma_crit, snap.iters(k), snap.pass(k)), ...
            'FontWeight', 'bold');
    end

    title(tl, sprintf('%s      |      %s / %s', study_label, patient_txt, session_txt), ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'tex');
    drawnow;
end


function cmap = gamma_colormap()
%GAMMA_COLORMAP 256-level map over gamma in [0, 2]: green for pass (<1), red for
%  fail (>1), yellow near the gamma=1 threshold. Used with caxis([0 2]).
    n = 256;
    x = linspace(0, 2, n)';
    r = min(1, max(0, (x - 0.7) / 0.6));      % ramps up to 1 approaching gamma=1
    g = min(1, max(0, (1.3 - x) / 0.6));      % ramps down past gamma=1
    b = 0.10 * ones(n, 1);
    cmap = [r, g, b];
end


function out_path = save_and_maybe_close(fh, CONFIG, name)
%SAVE_AND_MAYBE_CLOSE Save fh to CONFIG.figure_dir/<name>.png (when save_figures)
%  and close it (when close_figures). Returns the written path ('' if not saved).
    out_path = '';
    if isempty(fh) || ~ishandle(fh)
        return;
    end
    if isfield(CONFIG, 'save_figures') && CONFIG.save_figures
        if ~isempty(CONFIG.figure_dir) && ~isfolder(CONFIG.figure_dir)
            mkdir(CONFIG.figure_dir);
        end
        out_path = fullfile(CONFIG.figure_dir, [name '.png']);
        try
            exportgraphics(fh, out_path, 'Resolution', 150);
        catch
            print(fh, out_path, '-dpng', '-r150');   % fallback for older MATLAB
        end
        fprintf('  [FIG] Saved %s\n', out_path);
    end
    if isfield(CONFIG, 'close_figures') && CONFIG.close_figures
        close(fh);
    end
end


function tag = dose_file_tag(dose_filename)
%DOSE_FILE_TAG Short beam/segment tag parsed from a field-dose filename.
%  'dose_..._B15_112.mat' -> 'B15_112'. Falls back to the base filename (sans
%  extension) when the _B<beam>_<seg> pattern is absent.
    [~, base] = fileparts(char(dose_filename));
    tok = regexp(base, '_B(\d+)_(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        tag = sprintf('B%s_%s', tok{1}, tok{2});
    else
        tag = base;
    end
end


function C = merge_config(base, override)
%MERGE_CONFIG Shallow overlay: every field of override replaces base's field.
    C = base;
    if isempty(override) || ~isstruct(override)
        return;
    end
    fn = fieldnames(override);
    for i = 1:numel(fn)
        C.(fn{i}) = override.(fn{i});
    end
end


function CONFIG = default_config()
%DEFAULT_CONFIG Default parameters for the gamma-index TR-convergence study.
%  Any field may be overridden by the CONFIG struct passed to the main function.
    % Shared simulation defaults from get_default_config (single source of truth).
    CONFIG = get_default_config();

    % CT_1 field dose. The CT_3 counterpart filename is derived from this by
    % swapping the '_CT_1_' token for '_CT_3_' (same beam/segment, other CT image).
    CONFIG.dose_filename_ct1 = 'dose_1194203_Session_1_reference_CT_1_B15_112.mat';
    CONFIG.dose_filename_ct3 = '';   % '' => derive from the CT_1 name (recommended)

    CONFIG.cbct_filename_ct1 = 'CBCT1_resampled.mat';
    CONFIG.cbct_filename_ct3 = 'CBCT3_resampled.mat';

    % Sensor geometry, tissue model, forward/reconstruction params, the
    % convolution/noise chain, use_grid_padding and normalize are all inherited
    % from get_default_config (they matched this study's previous values). Only
    % the study-specific overrides and extra output fields are set below.

    % --- Convergence study: 20 time-reversal iterations, gamma after each ---
    CONFIG.num_time_reversal_iter = 20;  % more TR iterations than the default (1)
    CONFIG.gamma_n                = 3;    % gamma criterion, evaluated as n%/n mm

    % TR iterations snapshotted for the truth-vs-recon dose/gamma panels.
    CONFIG.snapshot_iters = [1, 10];

    % Noise-floor null hypothesis: reconstruct pure electronic noise (matched to
    % Study A's noise amplitude) through the same TR-convergence pipeline, so the
    % A/B signal curves can be shown sitting far above it. null_realizations > 1
    % averages the noise-floor pass rates over independent noise draws (the panel
    % uses the first realization). Requires convolution_kernel > 0 to be active.
    CONFIG.null_realizations = 1;

    % --- Output ---
    CONFIG.output_file   = '';   % '' => AnalysisResults/.../gamma_convergence per dose
    CONFIG.results_dir   = '';   % '' => same as figure_dir
    CONFIG.save_figures  = false; % batch driver sets true
    CONFIG.close_figures = false; % batch driver sets true (headless)
    CONFIG.figure_dir    = '';   % '' => AnalysisResults/<pt>/<sess>/gamma_convergence
end


function tables = define_tissue_tables()
%DEFINE_TISSUE_TABLES Tissue property lookup tables for HU thresholding.
%  (Ported from run_standalone_simulation.m so this study is self-contained.)

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
    tables.threshold_2.air_hu_threshold = -300;
    tables.threshold_2.air = struct( ...
        'density',     1.2, ...
        'sound_speed', 343, ...
        'alpha_coeff', 0, ...
        'alpha_power', 1.0, ...
        'gruneisen',   0);
end
