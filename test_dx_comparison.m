%% =========================================================================
%  TEST_DX_COMPARISON.m
%  Voxel-scale (dx,dy,dz) sweep for a SINGLE beam/segment pair, with a full
%  run_standalone_comparison-style VISUAL verification at every dx value.
%
%  WHAT IT DOES
%  ------------
%  This is the voxel_scale (dx) leg of study_optimization_sweeps.m, run on its
%  own and WITHOUT parallelization, so every dx value can be visually verified:
%
%   1. Runs the full two-dose (reference-CT + counterpart-CT) forward +
%      time-reversal pipeline at each dx value in DX.downscale_values.
%      dx is swept via CONFIG.downscale_factor, which scales the sim voxel size
%      dx,dy,dz UP by that factor (coarser grid). dt and Nt follow automatically
%      from the CFL condition, exactly as in study_optimization_sweeps.
%   2. The downscale_factor == 1 run is the REFERENCE. Every other dx is scored
%      against it for STABILITY (does the coarse recon reproduce the trusted
%      full-resolution recon?) and RUNTIME, using the same success criteria as
%      study_optimization_sweeps.
%   3. For EACH dx it also reproduces run_standalone_comparison.m's figures so
%      the pass rate can be inspected by eye:
%        - 2x3 dose-comparison panels (two reconstructed doses)
%        - 1x3 axial gamma / |error| / normalized-recon panel per gamma pair
%        - (optional) sensor-vs-dose view and p0 convergence history
%      The four standalone gamma pairs (each 3%/3mm) are computed at that dx's
%      OWN (downscaled) grid, so the number printed on a gamma panel is exactly
%      the pass rate scored for that dx.
%   4. Everything is printed to the command window and drawn to figure windows.
%      No parallel pool, no report/.mat files.
%
%  RELATION TO study_optimization_sweeps.m
%  ---------------------------------------
%   * STABILITY score: the swept recon is resampled BACK to the native dose
%     dims and gamma-compared to the reference (dx=1) recon, on the native grid
%     (identical to the sweep's stability metric).
%   * FIGURE / accuracy score: gamma is computed at the dx's own grid (identical
%     to how run_standalone_comparison would render it with that downscale
%     factor), so the figures and the accuracy column agree.
%
%  NOTE: run_standalone_comparison.m contains a hard-coded "dx = dx*2" line that
%  doubles the spacing on every run; that is NOT part of the optimization sweep
%  and is deliberately NOT reproduced here. dx here is controlled purely by
%  CONFIG.downscale_factor, matching study_optimization_sweeps.
%
%  Remote/HIPAA execution only. Do not run locally.
%% =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path (run from root).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% ======================= dx VALUES TO SWEEP =============================
%  EDIT THIS LINE to change which voxel scales are tested. Each value multiplies
%  dx,dy,dz. 1 = native resolution (always used as the reference baseline; it is
%  auto-added if you remove it).
DX.downscale_values = 2.^[1/4, 1/3, 1/2, 1, 2, 3];

%% ========================= BASE CONFIGURATION ===========================
%  Identical defaults to study_optimization_sweeps.m / run_standalone_comparison.m.
%  Every dx run uses these; only CONFIG.downscale_factor changes per run.

% Shared defaults from get_default_config (single source of truth); identical to
% study_optimization_sweeps.m. Only the fields below differ or are study-specific.
% downscale_factor stays at the default (1) and is overwritten per dx in the loop.
CONFIG = get_default_config();

% Single beam/segment pair to study (CT_3 dose; the default is the CT_1 dose).
CONFIG.dose_filename = 'dose_1194203_Session_1_reference_CT_3_B15_112.mat';
CONFIG.ct_pair      = [1, 3];
CONFIG.reference_ct = 1;

CONFIG.num_time_reversal_iter = 10;      % more TR iterations than the default (1)

% Deterministic electronic-noise draw so pass-rate differences between dx values
% reflect the grid change, not run-to-run noise variance (same as the sweep).
CONFIG.rng_seed              = 42;

% Success criteria (same thresholds as study_optimization_sweeps).
CONFIG.success_stability_pct = 98;       % criterion-1 stability threshold (%)
CONFIG.runtime_tol           = 0.02;     % allow 2% timing jitter for "<= default"

% Gaussian sigma (voxels) applied to dose slices for DISPLAY ONLY in the dose
% panels (fills recon speckle). 0 = disable.
CONFIG.viz_smooth_sigma = 1.0;

%% ========================= WHAT TO DISPLAY ==============================
%  Everything is shown in figure windows. Toggle individual figure groups here.
DISPLAY.dose_panels   = true;    % 2x3 two-recon dose comparison, per dx
DISPLAY.gamma_panels  = true;    % axial gamma/error/recon, per gamma pair, per dx
DISPLAY.sensor_view   = false;   % sensor-vs-dose-mask 1x3 view, per dx
DISPLAY.convergence   = false;   % p0 convergence history (reference-CT recon), per dx

% Figure saving. Files go to AnalysisResults/[PatientID]/[Session]/dx_comparison.
DISPLAY.save_figures  = true;
DISPLAY.figure_format = 'png';   % 'png' | 'fig' | 'both'
DISPLAY.figure_dir    = fullfile(CONFIG.working_dir, 'AnalysisResults', ...
                                 CONFIG.patient_id, CONFIG.session, 'dx_comparison');

%% ========================= PRIMARY GAMMA CRITERION ======================
GAMMA_CRITERIA = {3, 3, '3%/3mm'};       % single criterion, as in both scripts
PRIMARY = 1;

%% ===================== RESOLVE DOSE PAIR & CBCT PATHS ====================
%  (Identical resolution logic to run_standalone_comparison / the sweep.)

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

ct_tok = regexp(dose_basename_A, '_CT_(\d+)_', 'tokens', 'once');
if isempty(ct_tok)
    error('test_dx_comparison:NoCTtoken', ...
        'Dose filename "%s" has no _CT_k token.', dose_basename_A);
end
ct_a = str2double(ct_tok{1});
if ~ismember(ct_a, CONFIG.ct_pair)
    error('test_dx_comparison:CTnotInPair', ...
        'Dose CT index %d is not in CONFIG.ct_pair = [%s].', ct_a, num2str(CONFIG.ct_pair));
end
ct_b_all = CONFIG.ct_pair(CONFIG.ct_pair ~= ct_a);
ct_b     = ct_b_all(1);
ct_list  = [ct_a, ct_b];

if ~ismember(CONFIG.reference_ct, CONFIG.ct_pair)
    error('test_dx_comparison:BadReferenceCT', ...
        'CONFIG.reference_ct = %d must be one of CONFIG.ct_pair = [%s].', ...
        CONFIG.reference_ct, num2str(CONFIG.ct_pair));
end
recon_cbct_filepath = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', CONFIG.reference_ct));

dose_basename_B = regexprep(dose_basename_A, '_CT_\d+_', sprintf('_CT_%d_', ct_b));
dose_filepath_B = fullfile(processed_dir, dose_basename_B);

if ~isempty(CONFIG.cbct_file_override)
    cbct_filepath_A = CONFIG.cbct_file_override;
else
    cbct_filepath_A = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_a));
end
cbct_filepath_B = fullfile(processed_dir, sprintf('CBCT%d_resampled.mat', ct_b));

PATHS = struct();
PATHS.dose         = {dose_filepath_A, dose_filepath_B};
PATHS.cbct         = {cbct_filepath_A, cbct_filepath_B};
PATHS.ct_list      = ct_list;
PATHS.reference_ct = CONFIG.reference_ct;
PATHS.recon_cbct   = recon_cbct_filepath;
PATHS.label        = {sprintf('CT_%d (listed)', ct_a), ...
                      sprintf('CT_%d (counterpart)', ct_b)};

%% ========================= LOAD PLAN BEAM METADATA =======================

if ~isfield(CONFIG, 'beam_metadata') || isempty(CONFIG.beam_metadata)
    metadata_filepath = fullfile(processed_dir, 'metadata.mat');
    if isfile(metadata_filepath)
        try
            md = load(metadata_filepath, 'metadata');
            if isfield(md, 'metadata') && isfield(md.metadata, 'beam_metadata') ...
                    && ~isempty(md.metadata.beam_metadata)
                CONFIG.beam_metadata = md.metadata.beam_metadata;
                fprintf('  Loaded beam_metadata for %d beams.\n', numel(CONFIG.beam_metadata));
            else
                fprintf('  [WARN] metadata.mat present but no beam_metadata field.\n');
                CONFIG.beam_metadata = [];
            end
        catch ME
            fprintf('  [WARN] Failed to load metadata.mat: %s\n', ME.message);
            CONFIG.beam_metadata = [];
        end
    else
        fprintf('  [WARN] No metadata.mat in processed_dir; exclusion zone will be empty.\n');
        CONFIG.beam_metadata = [];
    end
end

% Shared tissue tables (same builder the pipeline uses).
CONFIG.tissue_tables = define_tissue_tables();
CONFIG.tissue_tables.uniform = struct( ...
    'density',     CONFIG.uniform_density, ...
    'sound_speed', CONFIG.uniform_sound_speed, ...
    'alpha_coeff', CONFIG.uniform_alpha_coeff, ...
    'alpha_power', 1.1, ...
    'gruneisen',   CONFIG.uniform_gruneisen);

%% ===================== BUILD THE dx RUN LIST ============================
%  Ensure the native (downscale_factor == 1) baseline is present, and put it
%  FIRST so it is available as the reference before any coarse run is scored.

dx_values = DX.downscale_values(:)';
if ~any(abs(dx_values - 1) < 1e-9)
    dx_values = [1, dx_values];
    fprintf('  [INFO] downscale_factor = 1 baseline added to the sweep list.\n');
end
% Sort ascending and move the 1.0 baseline to the front.
dx_values = unique(dx_values, 'stable');
[~, base_pos] = min(abs(dx_values - 1));
dx_values = [dx_values(base_pos), dx_values(setdiff(1:numel(dx_values), base_pos))];
nDx = numel(dx_values);

fprintf('=========================================================\n');
fprintf('  dx (voxel-scale) comparison + visual verification\n');
fprintf('  Patient %s / %s   dose: %s\n', CONFIG.patient_id, CONFIG.session, CONFIG.dose_filename);
fprintf('  Reference CT: CT_%d   (CT pair [%s])\n', CONFIG.reference_ct, num2str(CONFIG.ct_pair));
fprintf('  %d dx values: %s\n', nDx, mat2str(round(dx_values, 4)));
fprintf('  (downscale_factor = 1 is the reference baseline; no parallel pool)\n');
if DISPLAY.save_figures
    if ~isfolder(DISPLAY.figure_dir), mkdir(DISPLAY.figure_dir); end
    fprintf('  Saving figures (%s) to: %s\n', DISPLAY.figure_format, DISPLAY.figure_dir);
end
fprintf('=========================================================\n\n');

%% ========================= MAIN dx LOOP =================================
%  Serial. For each dx: run both doses, compute the four standalone gamma pairs
%  at the dx grid (for figures + accuracy), score stability/runtime vs the
%  baseline (native-resampled), draw the requested figures, and record a row.

ROWS = struct('dx', {}, 'label', {}, 'dims', {}, 'runtime', {}, 'iters', {}, ...
    'stability', {}, 'recon1_vs_dose1', {}, 'recon2_vs_dose2', {}, ...
    'crit1', {}, 'crit2', {}, 'success', {}, 'ok', {}, 'error', {});

baseline = [];   % filled by the dx == 1 run

for k = 1:nDx
    dx_val = dx_values(k);
    is_baseline = abs(dx_val - 1) < 1e-9;
    dx_lbl = fmt_dx(dx_val);

    fprintf('\n#################### dx run %d/%d : downscale_factor = %s%s ####################\n', ...
        k, nDx, dx_lbl, tern(is_baseline, '  [REFERENCE]', ''));

    row = struct('dx', dx_val, 'label', dx_lbl, 'dims', [NaN NaN NaN], ...
        'runtime', NaN, 'iters', [NaN NaN], 'stability', NaN, ...
        'recon1_vs_dose1', NaN, 'recon2_vs_dose2', NaN, ...
        'crit1', false, 'crit2', false, 'success', false, 'ok', false, 'error', '');

    try
        cfg = CONFIG;
        cfg.downscale_factor = dx_val;

        % --- run both doses at this dx ---
        o1 = test_field_once(PATHS.dose{1}, PATHS.cbct{1}, PATHS.ct_list(1), ...
                             PATHS.reference_ct, PATHS.recon_cbct, cfg);
        o2 = test_field_once(PATHS.dose{2}, PATHS.cbct{2}, PATHS.ct_list(2), ...
                             PATHS.reference_ct, PATHS.recon_cbct, cfg);
        o1.label = PATHS.label{1};
        o2.label = PATHS.label{2};

        if ~isequal(o1.gridSize, o2.gridSize)
            error('test_dx_comparison:GridMismatch', ...
                'The two doses ended on different grids ([%s] vs [%s]) at dx=%s.', ...
                num2str(o1.gridSize), num2str(o2.gridSize), dx_lbl);
        end

        % --- map to reference-CT (idx1) vs counterpart-CT (idx2) ---
        idx1 = find(PATHS.ct_list == PATHS.reference_ct, 1);
        idx2 = find(PATHS.ct_list ~= PATHS.reference_ct, 1);
        if isempty(idx2), tmp = setdiff([1, 2], idx1); idx2 = tmp(1); end
        outs = {o1, o2};
        oref = outs{idx1};
        ooth = outs{idx2};

        spacing_mm = oref.spacing_mm;

        dose1  = double(oref.doseGrid);      recon1 = double(oref.recon_dose);
        dose2  = double(ooth.doseGrid);      recon2 = double(ooth.recon_dose);
        label1 = oref.label;                 label2 = ooth.label;
        sensor1 = oref.sensor_mask;          sensor2 = ooth.sensor_mask;

        % recon_A / recon_B feed the dose-panel figure (listed = di 1).
        recon_A = double(o1.recon_dose);     recon_B = double(o2.recon_dose);

        % --- least-squares normalization (recon -> own-CT truth) ---
        if CONFIG.normalize
            recon1  = recon1 * least_squares_gain(dose1, recon1);
            recon2  = recon2 * least_squares_gain(dose2, recon2);
            recon_A = recon_A * least_squares_gain(double(o1.doseGrid), recon_A);
            recon_B = recon_B * least_squares_gain(double(o2.doseGrid), recon_B);
        end

        % --- the four standalone gamma pairs (this dx's own grid) ---
        gamma_pairs = {
            'dose1_vs_dose2',  dose1,  dose2,  sprintf('%s truth  vs  %s truth', label1, label2), sensor1;
            'recon1_vs_dose1', dose1,  recon1, sprintf('%s recon  vs  %s truth', label1, label1), sensor1;
            'recon2_vs_dose2', dose2,  recon2, sprintf('%s recon  vs  %s truth', label2, label2), sensor2;
            'recon2_vs_dose1', dose1,  recon2, sprintf('%s recon  vs  %s truth', label2, label1), sensor2;
        };
        gamma_results = struct();
        fprintf('  [Gamma @ dx=%s]  (10%% low-dose cutoff on reference)\n', dx_lbl);
        for pp = 1:size(gamma_pairs, 1)
            gr = run_gamma_pair(gamma_pairs{pp, 2}, gamma_pairs{pp, 3}, spacing_mm, GAMMA_CRITERIA);
            gr.name = gamma_pairs{pp, 1};  gr.title = gamma_pairs{pp, 4};
            gr.ref  = gamma_pairs{pp, 2};  gr.tgt   = gamma_pairs{pp, 3};
            gr.sensor_mask = gamma_pairs{pp, 5};
            gamma_results.(gamma_pairs{pp, 1}) = gr;
            fprintf('     %-34s %s = %s\n', gr.title, GAMMA_CRITERIA{PRIMARY, 3}, ...
                pass_str(gr.pass_rates(PRIMARY)));
        end

        runtime = o1.fwd_time + o1.tr_time + o2.fwd_time + o2.tr_time;
        row.dims            = oref.native_dims;
        row.runtime         = runtime;
        row.iters           = [o1.num_iters_done, o2.num_iters_done];
        row.recon1_vs_dose1 = gamma_results.recon1_vs_dose1.pass_rates(PRIMARY);
        row.recon2_vs_dose2 = gamma_results.recon2_vs_dose2.pass_rates(PRIMARY);
        row.ok              = true;

        % --- reference bookkeeping / stability scoring ---
        if is_baseline
            baseline = struct();
            baseline.recon1_native = single(oref.recon_native);
            baseline.spacing_native = oref.spacing_native;
            baseline.runtime = runtime;
            baseline.recon1_vs_dose1 = row.recon1_vs_dose1;
            row.stability = 100;        % vs itself
            row.crit1 = true;           % baseline trivially reproduces itself
            row.crit2 = false;
            row.success = true;
        elseif ~isempty(baseline)
            r1 = double(oref.recon_native);
            if ~isequal(size(r1), size(baseline.recon1_native))
                r1 = max(imresize3(r1, size(baseline.recon1_native)), 0);
            end
            stab = run_gamma_pair(double(baseline.recon1_native), r1, ...
                                  baseline.spacing_native, GAMMA_CRITERIA);
            row.stability = stab.pass_rates(PRIMARY);

            rt_ok = runtime <= baseline.runtime * (1 + CONFIG.runtime_tol);
            row.crit1 = (row.stability >= CONFIG.success_stability_pct) && rt_ok;
            row.crit2 = (row.recon1_vs_dose1 > baseline.recon1_vs_dose1) && rt_ok;
            row.success = row.crit1 || row.crit2;
        else
            fprintf('  [WARN] Baseline (dx=1) not available yet; stability not scored.\n');
        end

        % --- FIGURES: reproduce run_standalone_comparison's results ---
        fig_tag = sprintf('dx=%s', dx_lbl);
        dtag    = safe_tag(dx_lbl);   % filesystem-safe dx label

        if DISPLAY.sensor_view
            dose_mask_vis = double(o1.doseGrid) >= 0.10 * max(double(o1.doseGrid(:)));
            plot_sensor_dose_planes(dose_mask_vis, logical(o1.sensor_mask), ...
                spacing_mm, o1.density, CONFIG, fig_tag);
            save_fig(DISPLAY, sprintf('dx_%s_sensor_view', dtag));
        end

        if DISPLAY.dose_panels
            plot_dose_panels(recon_A, recon_B, o1.sensor_mask, o1.density, spacing_mm, ...
                sprintf('Dose Comparison: Two Reconstructed Doses  [%s]', fig_tag), ...
                {o1.label, o2.label}, CONFIG.viz_smooth_sigma);
            save_fig(DISPLAY, sprintf('dx_%s_dose_panels', dtag));
        end

        if DISPLAY.convergence
            plot_convergence_history(oref.conv_max_pressure, oref.conv_rel_change, ...
                oref.num_iters_done, CONFIG.convergence_tol, max(oref.p0(:)), fig_tag);
            save_fig(DISPLAY, sprintf('dx_%s_convergence', dtag));
        end

        if DISPLAY.gamma_panels
            pair_fn = fieldnames(gamma_results);
            for pf = 1:numel(pair_fn)
                gr = gamma_results.(pair_fn{pf});
                [~, mdi] = max(gr.ref(:));
                [~, ~, cz_g] = ind2sub(size(gr.ref), mdi);
                plot_gamma_and_error_axial(gr, gr.ref, gr.tgt, gr.sensor_mask, cz_g, ...
                    sprintf('%s   [%s]', gr.title, fig_tag));
                save_fig(DISPLAY, sprintf('dx_%s_gamma_%s', dtag, gr.name));
            end
        end
        drawnow;

        fprintf('  runtime = %.1f s   dims [%s]   iters %s\n', ...
            runtime, num2str(row.dims), mat2str(row.iters));

    catch ME
        row.ok = false;
        row.error = sprintf('%s | %s', ME.identifier, ME.message);
        fprintf('  [FAIL] dx=%s : %s\n', dx_lbl, ME.message);
    end

    ROWS(end+1) = row; %#ok<SAGROW>
end

%% ========================= SUMMARY TABLE ================================

plabel = GAMMA_CRITERIA{PRIMARY, 3};
fprintf('\n\n================================================================================\n');
fprintf('  dx COMPARISON SUMMARY   (criterion: %s ; success thresholds from study_optimization_sweeps)\n', plabel);
fprintf('================================================================================\n');
fprintf('  SUCCESS if EITHER (vs the dx=1 reference):\n');
fprintf('    (1) stability gamma(default recon vs this recon) >= %g%%  AND  runtime <= default\n', ...
    CONFIG.success_stability_pct);
fprintf('    (2) recon1-vs-truth gamma > default''s               AND  runtime <= default\n');
fprintf('    (runtime tolerance for "<= default": %.0f%%)\n\n', 100 * CONFIG.runtime_tol);

fprintf('  %-9s %-13s %8s %8s %9s %10s %10s %5s %5s %8s\n', ...
    'dx', 'dims', 'runtime', 'dRun%', 'iters', ['stab-' plabel], ['R1vT-' plabel], 'c1', 'c2', 'RESULT');
fprintf('  %s\n', repmat('-', 1, 96));

base_runtime = NaN;
for i = 1:numel(ROWS)
    if ROWS(i).ok && abs(ROWS(i).dx - 1) < 1e-9, base_runtime = ROWS(i).runtime; end
end
for i = 1:numel(ROWS)
    r = ROWS(i);
    if ~r.ok
        fprintf('  %-9s %-13s %8s   ERROR: %s\n', r.label, '', 'FAIL', r.error);
        continue;
    end
    if ~isnan(base_runtime) && base_runtime > 0
        dRun = 100 * (r.runtime - base_runtime) / base_runtime;
    else
        dRun = NaN;
    end
    res_str = 'no';
    if r.success, res_str = 'SUCCESS'; end
    fprintf('  %-9s %-13s %7.1fs %+7.1f %9s %9s %9s %5s %5s %8s\n', ...
        r.label, ['[' num2str(r.dims) ']'], r.runtime, dRun, mat2str(r.iters), ...
        pass_str(r.stability), pass_str(r.recon1_vs_dose1), ...
        yn(r.crit1), yn(r.crit2), res_str);
end
fprintf('  %s\n', repmat('-', 1, 96));
fprintf('  stab-%s : gamma(reference dx=1 recon  vs  this dx recon), native grid  -> STABILITY\n', plabel);
fprintf('  R1vT-%s : gamma(reference-CT recon  vs  its RayStation truth), this dx grid -> ACCURACY (matches the gamma figures)\n', plabel);
fprintf('\ndx comparison complete. %d/%d dx runs succeeded.\n', sum([ROWS.ok]), numel(ROWS));


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function s = fmt_dx(v)
%FMT_DX  Compact label for a downscale factor.
    if abs(v - round(v)) < 1e-9
        s = sprintf('%d', round(v));
    else
        s = sprintf('%.4g', v);
    end
end

function s = pass_str(v)
%PASS_STR  Format a gamma pass rate (NaN -> FAIL).
    if isnan(v), s = 'FAIL'; else, s = sprintf('%.2f%%', v); end
end

function s = yn(tf)
    if tf, s = 'Y'; else, s = 'n'; end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

function s = safe_tag(lbl)
%SAFE_TAG  Turn a dx label (e.g. '1.189') into a filesystem-safe token ('1p189').
    s = regexprep(lbl, '\.', 'p');
    s = regexprep(s, '[^A-Za-z0-9_]', '_');
end

function save_fig(DISPLAY, basename)
%SAVE_FIG  Save the current figure to DISPLAY.figure_dir (png / fig / both).
    if ~isfield(DISPLAY, 'save_figures') || ~DISPLAY.save_figures, return; end
    if ~isfolder(DISPLAY.figure_dir), mkdir(DISPLAY.figure_dir); end
    fig = gcf;
    fmt = 'png';
    if isfield(DISPLAY, 'figure_format') && ~isempty(DISPLAY.figure_format)
        fmt = lower(DISPLAY.figure_format);
    end
    try
        if any(strcmpi(fmt, {'png', 'both'}))
            exportgraphics(fig, fullfile(DISPLAY.figure_dir, [basename '.png']), ...
                'Resolution', 150);
        end
        if any(strcmpi(fmt, {'fig', 'both'}))
            savefig(fig, fullfile(DISPLAY.figure_dir, [basename '.fig']));
        end
    catch ME
        fprintf('  [WARN] Could not save figure "%s": %s\n', basename, ME.message);
    end
end


function out = test_field_once(dose_filepath, cbct_filepath, ct_this, ...
        reference_ct, recon_cbct_filepath, CONFIG)
%TEST_FIELD_ONCE  Single-dose forward + iterative TR reconstruction.
%  Same pipeline as run_standalone_comparison.m / study_optimization_sweeps.m
%  simulate_field_once, headless. Returns BOTH the working-grid quantities (for
%  the run_standalone_comparison-style figures at this dx) AND the recon
%  resampled back to native dims (for the cross-dx stability score).
%
%  OUTPUT struct fields:
%    Working grid (this dx, = gridSize_orig, water-expanded if the sensor grew):
%      .recon_dose  .doseGrid  .doseMask  .sensor_mask  .density  .p0
%      .spacing_mm  .gridSize  .num_pulses
%      .conv_max_pressure  .conv_rel_change  .num_iters_done
%      .fwd_time  .tr_time
%    Native (pre-downscale):
%      .recon_native  .dose_native  .spacing_native  .native_dims

    out = struct();
    blind_recon = (ct_this ~= reference_ct);

    %% ---- LOAD DOSE ----
    if ~isfile(dose_filepath)
        error('test_field_once:NoDose', 'Dose file not found: %s', dose_filepath);
    end
    dose_data   = load(dose_filepath);
    dose_fields = fieldnames(dose_data);

    fd_gantry_angle   = [];
    step15_spacing_mm = [];

    if isfield(dose_data, 'field_dose')
        fd = dose_data.field_dose;
        if ~isfield(fd, 'dose_Gy'), error('field_dose struct missing dose_Gy field.'); end
        if (isfield(fd, 'is_sparse') && fd.is_sparse) || issparse(fd.dose_Gy)
            if ~isfield(fd, 'dose_dims'), error('field_dose.dose_dims missing.'); end
            doseGrid = reshape(full(fd.dose_Gy), fd.dose_dims);
        else
            doseGrid = double(fd.dose_Gy);
        end
        if isfield(fd, 'spacing') && ~isempty(fd.spacing)
            step15_spacing_mm = fd.spacing(:)';
        end
        if isfield(fd, 'meterset') && ~isempty(fd.meterset) && fd.meterset > 0
            CONFIG.meterset = fd.meterset;
        end
        if isfield(fd, 'gantry_angle')
            fd_gantry_angle = fd.gantry_angle;
        end
    elseif isfield(dose_data, 'total_rs_dose_sparse')
        if ~isfield(dose_data, 'total_rs_dose_dims'), error('total_rs_dose_dims missing.'); end
        doseGrid = reshape(full(dose_data.total_rs_dose_sparse), dose_data.total_rs_dose_dims);
    elseif isfield(dose_data, 'total_rs_dose')
        doseGrid = dose_data.total_rs_dose;
    elseif isfield(dose_data, 'dose_Gy')
        doseGrid = dose_data.dose_Gy;
    elseif numel(dose_fields) == 1
        doseGrid = dose_data.(dose_fields{1});
    else
        error('Cannot auto-detect dose variable. Fields: %s', strjoin(dose_fields, ', '));
    end
    doseGrid = double(doseGrid);
    if ~isnumeric(doseGrid) || ndims(doseGrid) ~= 3
        error('Dose data must be a 3D numeric array.');
    end

    gridSize = size(doseGrid);
    Nx = gridSize(1); Ny = gridSize(2); Nz = gridSize(3);

    %% ---- LOAD CBCT ----
    if ~isfile(cbct_filepath)
        error('test_field_once:NoCBCT', 'CBCT file not found: %s', cbct_filepath);
    end
    cbct_data = load(cbct_filepath);
    if isfield(cbct_data, 'CBCT1_resampled')
        sct = cbct_data.CBCT1_resampled;
    elseif isfield(cbct_data, 'CBCT3_resampled')
        sct = cbct_data.CBCT3_resampled;
    else
        error('CBCT1_resampled / CBCT3_resampled not found in %s', cbct_filepath);
    end
    if ~isfield(sct, 'cubeHU'), error('CBCT struct missing cubeHU.'); end

    if isfield(sct, 'spacing') && ~isempty(sct.spacing)
        spacing_mm = sct.spacing(:)';
    elseif ~isempty(step15_spacing_mm)
        spacing_mm = step15_spacing_mm;
    else
        error('No spacing available (CBCT or field_dose).');
    end

    if ~isequal(gridSize, size(sct.cubeHU))
        error('Dose grid [%s] != CBCT grid [%s].', ...
            num2str(gridSize), num2str(size(sct.cubeHU)));
    end

    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end

    % Native reference (BEFORE any downscale) used for the stability comparison.
    dose_native    = doseGrid;
    spacing_native = spacing_mm;
    nativeDims     = [Nx, Ny, Nz];

    %% ---- GRID DOWNSCALING (scales dx,dy,dz up by downscale_factor) ----
    if CONFIG.downscale_factor ~= 1
        df     = CONFIG.downscale_factor;
        new_Nx = max(2, round(Nx / df));
        new_Ny = max(2, round(Ny / df));
        new_Nz = max(2, round(Nz / df));

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
        sct.spacing = spacing_mm;
        Nx = new_Nx; Ny = new_Ny; Nz = new_Nz;
        gridSize = [Nx, Ny, Nz];
    end
    dx = spacing_mm(1) / 1000;
    dy = spacing_mm(2) / 1000;
    dz = spacing_mm(3) / 1000;

    % Downscaled-native dims: the anatomy block size before FFT padding / sensor
    % grid expansion (equals nativeDims when downscale_factor == 1).
    dsDims = [Nx, Ny, Nz];

    %% ---- ACOUSTIC MEDIUM (forward = TRUE geometry) ----
    medium = create_acoustic_medium(sct, CONFIG);
    medium = apply_standalone_medium_overrides(medium, sct, CONFIG);

    if blind_recon
        if ~isfile(recon_cbct_filepath)
            error('Reference CBCT for recon not found: %s', recon_cbct_filepath);
        end
        rec_cbct = load(recon_cbct_filepath);
        if isfield(rec_cbct, 'CBCT1_resampled')
            sct_recon = rec_cbct.CBCT1_resampled;
        elseif isfield(rec_cbct, 'CBCT3_resampled')
            sct_recon = rec_cbct.CBCT3_resampled;
        else
            error('CBCT1_resampled / CBCT3_resampled not found in %s', recon_cbct_filepath);
        end
        if ~isfield(sct_recon, 'cubeHU'), error('Reference CBCT missing cubeHU.'); end
        if ~isequal(size(sct_recon.cubeHU), nativeDims) && CONFIG.downscale_factor == 1
            error('Reference CBCT grid [%s] != dose grid [%s].', ...
                num2str(size(sct_recon.cubeHU)), num2str(nativeDims));
        end
        if CONFIG.downscale_factor ~= 1
            sct_recon.cubeHU = imresize3(sct_recon.cubeHU, size(sct.cubeHU));
            if isfield(sct_recon, 'bodyMask')
                sct_recon.bodyMask = imresize3(single(sct_recon.bodyMask), size(sct.cubeHU), 'nearest') > 0.5;
            end
            if isfield(sct_recon, 'couchMask')
                sct_recon.couchMask = imresize3(single(sct_recon.couchMask), size(sct.cubeHU), 'nearest') > 0.5;
            end
        end
        sct_recon.spacing = spacing_mm;
        medium_recon = create_acoustic_medium(sct_recon, CONFIG);
        medium_recon = apply_standalone_medium_overrides(medium_recon, sct_recon, CONFIG);
    else
        medium_recon = medium;
    end

    %% ---- INITIAL PRESSURE p0 ----
    if isfield(sct, 'bodyMask')
        doseGrid = doseGrid .* double(sct.bodyMask);
    end
    meterset       = CONFIG.meterset;
    num_pulses     = ceil(meterset / CONFIG.dose_per_pulse_cGy);
    dose_per_pulse = doseGrid / num_pulses;

    p0 = dose_per_pulse .* medium.gruneisen .* medium.density;
    p0 = smooth(p0);

    doseThreshold = 0.01 * max(doseGrid(:));
    doseMask      = doseGrid > doseThreshold;
    if ~any(doseMask(:)) || max(p0(:)) == 0
        error('No significant dose or zero initial pressure.');
    end

    %% ---- OPTIMAL GRID PADDING ----
    Nx_orig = Nx; Ny_orig = Ny; Nz_orig = Nz;
    gridSize_orig = gridSize;
    medium_orig   = medium;

    % Core index ranges that locate the downscaled-native anatomy block within
    % the final cropped recon grid. Overwritten if the sensor expands the grid.
    core_xr = 1:dsDims(1);
    core_yr = 1:dsDims(2);
    core_zr = 1:dsDims(3);

    if CONFIG.use_grid_padding
        Nx_pad = find_optimal_kwave_size(Nx, CONFIG.pml_size);
        Ny_pad = find_optimal_kwave_size(Ny, CONFIG.pml_size);
        Nz_pad = find_optimal_kwave_size(Nz, CONFIG.pml_size);
    else
        Nx_pad = Nx; Ny_pad = Ny; Nz_pad = Nz;
    end

    did_pad = ~isequal([Nx_pad, Ny_pad, Nz_pad], [Nx, Ny, Nz]);
    if did_pad
        [medium, p0] = pad_forward(medium, p0, Nx, Ny, Nz, Nx_pad, Ny_pad, Nz_pad);
        if blind_recon
            medium_recon = pad_medium_to(medium_recon, [Nx_pad, Ny_pad, Nz_pad]);
        else
            medium_recon = medium;
        end
        Nx = Nx_pad; Ny = Ny_pad; Nz = Nz_pad;
        gridSize = [Nx, Ny, Nz];
    end

    %% ---- SENSOR PLACEMENT ----
    sensor      = struct();
    sensor.mask = zeros(Nx, Ny, Nz);

    switch CONFIG.sensor_placement_method
        case 'full_plane_anterior'
            sensor.mask(CONFIG.sensor_x_index, :, :) = 1;
        case 'full_plane_lateral'
            sensor.mask(:, CONFIG.sensor_y_index, :) = 1;
        case 'determine_sensor_mask'
            sct_for_sensor = sct;
            if ~isfield(sct_for_sensor, 'couchMask')
                sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
            end
            if ~isfield(sct_for_sensor, 'origin')
                sct_for_sensor.origin = [0, 0, 0];
            end
            sct_for_sensor.spacing = spacing_mm;

            field_dose_for_sensor = struct();
            field_dose_for_sensor.dose_Gy = doseGrid;
            if ~isempty(fd_gantry_angle)
                field_dose_for_sensor.gantry_angle = fd_gantry_angle;
            else
                field_dose_for_sensor.gantry_angle = 0;
            end
            field_dose_for_sensor.origin     = sct_for_sensor.origin;
            field_dose_for_sensor.spacing    = spacing_mm;
            field_dose_for_sensor.dimensions = [Nx_orig, Ny_orig, Nz_orig];

            beam_meta = [];
            if isfield(CONFIG, 'beam_metadata') && ~isempty(CONFIG.beam_metadata)
                beam_meta = CONFIG.beam_metadata;
            end

            [sensor_mask_orig, sensor_info_orig] = determine_sensor_mask( ...
                sct_for_sensor, field_dose_for_sensor, beam_meta, CONFIG);

            % --- GRID EXPANSION HANDLING (water-pad + re-FFT-pad) ---
            gp = sensor_info_orig.grid_pad;
            if gp.expanded
                density_unp    = medium.density(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                soundSpeed_unp = medium.sound_speed(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                if numel(medium.alpha_coeff) > 1
                    alphaCoeff_unp = medium.alpha_coeff(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                else
                    alphaCoeff_unp = medium.alpha_coeff;
                end
                gruneisen_unp  = medium.gruneisen(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
                p0_unp         = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

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

                medium_orig = struct('density', density_exp, 'sound_speed', soundSpeed_exp, ...
                    'alpha_coeff', alphaCoeff_exp, 'gruneisen', gruneisen_exp);

                if blind_recon
                    medium_recon = crop_medium(medium_recon, [Nx_orig, Ny_orig, Nz_orig]);
                    medium_recon = expand_medium(medium_recon, [Nx_exp, Ny_exp, Nz_exp], xr, yr, zr);
                end

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

                % The anatomy block now lives at (xr,yr,zr) in the expanded grid,
                % which survives the FFT-pad -> crop step unchanged.
                core_xr = xr; core_yr = yr; core_zr = zr;

                Nx_orig = Nx_exp; Ny_orig = Ny_exp; Nz_orig = Nz_exp;
                gridSize_orig = [Nx_orig, Ny_orig, Nz_orig];

                if CONFIG.use_grid_padding
                    Nx_pad2 = find_optimal_kwave_size(Nx_orig, CONFIG.pml_size);
                    Ny_pad2 = find_optimal_kwave_size(Ny_orig, CONFIG.pml_size);
                    Nz_pad2 = find_optimal_kwave_size(Nz_orig, CONFIG.pml_size);
                else
                    Nx_pad2 = Nx_orig; Ny_pad2 = Ny_orig; Nz_pad2 = Nz_orig;
                end

                if ~isequal([Nx_pad2, Ny_pad2, Nz_pad2], [Nx_orig, Ny_orig, Nz_orig])
                    [medium, p0] = pad_forward(medium, p0, Nx_orig, Ny_orig, Nz_orig, ...
                                               Nx_pad2, Ny_pad2, Nz_pad2);
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

            m1 = min(Nx, size(sensor_mask_orig, 1));
            m2 = min(Ny, size(sensor_mask_orig, 2));
            m3 = min(Nz, size(sensor_mask_orig, 3));
            sensor.mask(1:m1, 1:m2, 1:m3) = double(sensor_mask_orig(1:m1, 1:m2, 1:m3));
        otherwise
            error('test_field_once:UnknownSensor', ...
                'This study supports determine_sensor_mask / full_plane_* only (got "%s").', ...
                CONFIG.sensor_placement_method);
    end

    if sum(sensor.mask(:)) == 0
        error('Sensor mask is empty.');
    end

    %% ---- k-WAVE GRID & MEDIUM SETUP ----
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    maxC = max(medium.sound_speed(:));
    minC = min(medium.sound_speed(medium.sound_speed > 0));
    dt   = CONFIG.cfl_number * min([dx, dy, dz]) / maxC;

    gridDiag = sqrt((Nx*dx)^2 + (Ny*dy)^2 + (Nz*dz)^2);
    simTime  = 2.5 * gridDiag / minC;
    Nt       = ceil(simTime / dt);

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

    kmedium_fwd             = struct();
    kmedium_fwd.density     = medium.density;
    kmedium_fwd.sound_speed = medium.sound_speed;
    kmedium_fwd.alpha_coeff = medium.alpha_coeff;
    kmedium_fwd.alpha_power = 1.1;

    kmedium_recon             = struct();
    kmedium_recon.density     = medium_recon.density;
    kmedium_recon.sound_speed = medium_recon.sound_speed;
    kmedium_recon.alpha_coeff = medium_recon.alpha_coeff;
    kmedium_recon.alpha_power = 1.1;

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

    %% ---- FORWARD SIMULATION ----
    source_fwd    = struct();
    source_fwd.p0 = p0;

    fwd_tic    = tic;
    sensorData = kspaceFirstOrder3D(kgrid, kmedium_fwd, source_fwd, sensor, inputArgs{:});
    fwd_time   = toc(fwd_tic);

    sensorData          = smooth(sensorData);
    FS                  = 1 / kgrid.dt;
    sensorData_measured = sensorData;

    %% ---- PULSE / FREQUENCY RESPONSE / NOISE / DECONV ----
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

        sensorData_conv = real(ifft(fft(sensorData_cpu, [], 2) .* H, [], 2));
        sensorData_resp = gaussianFilter(sensorData_conv, FS, 0.35e6, 100, true);

        % Deterministic electronic-noise draw so the stability comparison isolates
        % the swept parameter (not run-to-run RNG variance).
        if isfield(CONFIG, 'rng_seed') && ~isempty(CONFIG.rng_seed)
            rng(CONFIG.rng_seed);
        end
        noise_amp        = conv_noise_level * max(abs(sensorData_resp(:)));
        sensorData_noisy = sensorData_resp + noise_amp * randn(size(sensorData_resp));

        sensorData_deconv = real(ifft( ...
            fft(sensorData_noisy, [], 2) .* H_conj ./ (H_power + conv_deconv_lambda), [], 2));

        sensorData          = single(sensorData_deconv);
        sensorData_measured = single(sensorData_deconv);
    else
        sensorData          = gaussianFilter(sensorData, FS, 0.35e6, 100, true);
        sensorData_measured = sensorData;
    end

    %% ---- ITERATIVE TIME REVERSAL ----
    reconPressure      = zeros(gridSize);
    reconPressure_prev = zeros(gridSize);
    num_iters_done     = 0;

    % Convergence tracking (for the p0-convergence figure).
    conv_max_pressure = zeros(CONFIG.num_time_reversal_iter, 1);
    conv_rel_change   = nan(CONFIG.num_time_reversal_iter, 1);

    tr_total_tic = tic;
    for tr_iter = 1:CONFIG.num_time_reversal_iter
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
        reconPressure  = max(reconPressure, 0);
        num_iters_done = tr_iter;

        conv_max_pressure(tr_iter) = max(reconPressure(:));

        converged = false;
        if tr_iter > 1
            norm_prev = norm(reconPressure_prev(:));
            if norm_prev > 0
                rel_change = norm(reconPressure(:) - reconPressure_prev(:)) / norm_prev;
            else
                rel_change = Inf;
            end
            conv_rel_change(tr_iter) = rel_change;
            if rel_change < CONFIG.convergence_tol
                converged = true;
            end
        end
        reconPressure_prev = reconPressure;

        if converged, break; end

        if tr_iter < CONFIG.num_time_reversal_iter
            source_resid    = struct();
            source_resid.p0 = reconPressure;
            sensorDataRecon = kspaceFirstOrder3D(kgrid, kmedium_recon, source_resid, sensor, inputArgs{:});
            sensorData      = sensorData + (sensorData_measured - sensorDataRecon);
        end
    end
    reconPressure = gather(reconPressure) * CONFIG.correction_factor;
    tr_time = toc(tr_total_tic);

    %% ---- CROP TO ORIGINAL (padded -> gridSize_orig) ----
    if did_pad || ~isequal(size(reconPressure), gridSize_orig)
        reconPressure = reconPressure(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
        Nx = Nx_orig; Ny = Ny_orig; Nz = Nz_orig;
        gridSize = gridSize_orig;
        medium   = medium_orig;
        sensor.mask = sensor.mask(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);
    end

    if blind_recon
        medium = crop_medium(medium_recon, gridSize);
    end

    p0 = p0(1:Nx_orig, 1:Ny_orig, 1:Nz_orig);

    %% ---- PRESSURE SCALE CORRECTION ----
    if CONFIG.use_pressure_scale_correction
        p0_max_orig = max(p0, [], 'all');
        recon_max   = max(reconPressure(:));
        if recon_max > 0
            reconPressure = reconPressure * (p0_max_orig / recon_max);
        end
    end

    %% ---- PRESSURE -> DOSE ----
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
    recon_dose(medium.density < 100) = 0;

    %% ---- NATIVE-RESAMPLED RECON (for the cross-dx stability score) ----
    % Extract the anatomy core (dsDims) from the working grid, then resample to
    % native dims so the stability gamma runs at original resolution.
    recon_core = recon_dose(core_xr, core_yr, core_zr);
    if ~isequal(size(recon_core), nativeDims)
        recon_native = max(imresize3(recon_core, nativeDims), 0);
    else
        recon_native = recon_core;
    end

    %% ---- PACK OUTPUT ----
    % Working grid (this dx) for the run_standalone_comparison-style figures:
    out.recon_dose        = recon_dose;
    out.doseGrid          = doseGrid;
    out.doseMask          = doseMask;
    out.sensor_mask       = sensor.mask;
    out.density           = medium.density;
    out.p0                = p0;
    out.spacing_mm        = spacing_mm;
    out.gridSize          = gridSize_orig;
    out.num_pulses        = num_pulses;
    out.conv_max_pressure = conv_max_pressure;
    out.conv_rel_change   = conv_rel_change;
    out.num_iters_done    = num_iters_done;
    out.fwd_time          = fwd_time;
    out.tr_time           = tr_time;
    % Native (for stability scoring):
    out.recon_native      = recon_native;
    out.dose_native       = dose_native;
    out.spacing_native    = spacing_native;
    out.native_dims       = nativeDims;
end


function [medium, p0] = pad_forward(medium, p0, Nx, Ny, Nz, Nxp, Nyp, Nzp)
%PAD_FORWARD  Water-pad the forward medium + p0 into an [Nxp Nyp Nzp] grid.
    density_pad    = ones(Nxp, Nyp, Nzp) * 1000;
    soundSpeed_pad = ones(Nxp, Nyp, Nzp) * 1540;
    alphaCoeff_pad = zeros(Nxp, Nyp, Nzp);
    gruneisen_pad  = zeros(Nxp, Nyp, Nzp);
    p0_pad         = zeros(Nxp, Nyp, Nzp);

    density_pad(1:Nx, 1:Ny, 1:Nz)    = medium.density;
    soundSpeed_pad(1:Nx, 1:Ny, 1:Nz) = medium.sound_speed;
    if numel(medium.alpha_coeff) > 1
        alphaCoeff_pad(1:Nx, 1:Ny, 1:Nz) = medium.alpha_coeff;
    else
        alphaCoeff_pad(:) = medium.alpha_coeff;
    end
    gruneisen_pad(1:Nx, 1:Ny, 1:Nz) = medium.gruneisen;
    p0_pad(1:Nx, 1:Ny, 1:Nz)        = p0;

    medium.density     = density_pad;
    medium.sound_speed = soundSpeed_pad;
    medium.alpha_coeff = alphaCoeff_pad;
    medium.gruneisen   = gruneisen_pad;
    p0 = p0_pad;
end


function medium = apply_standalone_medium_overrides(medium, sct, config)
%APPLY_STANDALONE_MEDIUM_OVERRIDES  Force-uniform toggles + coupling bath.
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
%DEFINE_TISSUE_TABLES  Tissue property lookup tables for HU thresholding.
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


function m = pad_medium_to(m, sz)
%PAD_MEDIUM_TO  Water-pad the four medium fields into an sz grid.
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
%EXPAND_MEDIUM  Place cropped medium fields into a water-filled sz grid.
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
%CROP_MEDIUM  Crop the four medium fields to the leading sz sub-block.
    m.density     = m.density(1:sz(1),     1:sz(2),     1:sz(3));
    m.sound_speed = m.sound_speed(1:sz(1), 1:sz(2),     1:sz(3));
    if numel(m.alpha_coeff) > 1
        m.alpha_coeff = m.alpha_coeff(1:sz(1), 1:sz(2), 1:sz(3));
    end
    m.gruneisen   = m.gruneisen(1:sz(1),   1:sz(2),     1:sz(3));
end


function g = least_squares_gain(rs_truth, recon)
%LEAST_SQUARES_GAIN  Scalar gain aligning a recon to its RS truth (relative norm).
%  g = sum(rs.*recon)/sum(recon.^2) over the truth's 10% low-dose region; the g
%  minimizing ||rs - g*recon||^2 there. Falls back to 1 for empty/zero inputs.
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


function gr = run_gamma_pair(ref_vol, tgt_vol, spacing_mm, criteria)
%RUN_GAMMA_PAIR  Gamma index of tgt_vol against ref_vol at each criterion row.
%  ref_vol (ground truth) sets the 10% low-dose evaluation mask. Returns .maps
%  (per-criterion gamma maps, retained for plotting), .pass_rates, .criteria,
%  .cutoff_Gy, .eval_mask. Mirrors run_standalone_comparison's run_gamma_pair.
    nCrit = size(criteria, 1);
    gr = struct('maps', {cell(nCrit, 1)}, 'pass_rates', nan(nCrit, 1), ...
                'criteria', {criteria}, 'cutoff_Gy', 0, 'eval_mask', []);

    if isempty(ref_vol) || isempty(tgt_vol) || max(ref_vol(:)) <= 0
        gr.eval_mask = false(size(ref_vol));
        return;
    end

    ref_struct.start = [0, 0, 0]; ref_struct.width = spacing_mm; ref_struct.data = double(ref_vol);
    tgt_struct.start = [0, 0, 0]; tgt_struct.width = spacing_mm; tgt_struct.data = double(tgt_vol);

    cutoff    = 0.10 * max(ref_vol(:));
    eval_mask = ref_vol >= cutoff;
    gr.cutoff_Gy = cutoff;
    gr.eval_mask = eval_mask;

    for gc = 1:nCrit
        pct = criteria{gc, 1};
        dta = criteria{gc, 2};
        try
            gmap = CalcGamma(ref_struct, tgt_struct, pct, dta, ...
                'local', 0, 'limit', dta * 2, 'restrict', 1);
            gr.maps{gc}       = gmap;
            gr.pass_rates(gc) = 100 * mean(gmap(eval_mask) <= 1);
        catch ME
            warning('Gamma failed: %s', ME.message);
            gr.maps{gc}       = [];
            gr.pass_rates(gc) = NaN;
        end
    end
end


function plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, density, config, tag)
%PLOT_SENSOR_DOSE_PLANES  1x3 anatomical view of sensor geometry vs dose mask.
    if nargin < 6, tag = ''; end
    if ~isequal(size(sensor_mask), size(dose_mask))
        sz = size(dose_mask);
        sensor_mask = sensor_mask(1:sz(1), 1:sz(2), 1:sz(3));
    end
    [Nx3, Ny3, Nz3] = size(dose_mask);

    x_ax = (1:Nx3) * spacing_mm(1);
    y_ax = (1:Ny3) * spacing_mm(2);
    z_ax = (1:Nz3) * spacing_mm(3);

    dose_cor  = squeeze(any(dose_mask,   2));
    dose_sag  = squeeze(any(dose_mask,   1));
    dose_axi  = squeeze(any(dose_mask,   3));
    sens_cor  = squeeze(any(sensor_mask, 2));
    sens_sag  = squeeze(any(sensor_mask, 1));
    sens_axi  = squeeze(any(sensor_mask, 3));

    have_density = ~isempty(density) && isequal(size(density), size(dose_mask));
    wl_center = 1050; wl_width = 350;
    wl_min    = wl_center - wl_width / 2;
    if have_density
        ct_projs = {
            squeeze(mean(double(density), 2))', ...
            squeeze(mean(double(density), 1))', ...
            squeeze(mean(double(density), 3))'
        };
    end

    figure('Name', ['Sensor Placement vs Dose Mask  ' tag], 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [80, 80, 1300, 420]);
    sgtitle(sprintf('Sensor Placement vs Dose Mask  (\\geq10%% max)   |   %s', tag), ...
        'FontWeight', 'bold', 'FontSize', 11);

    view_data = {
        dose_cor', sens_cor', x_ax, z_ax, 'X (mm)', 'Z (mm)', 'Coronal  (mean-proj along Y)';
        dose_sag', sens_sag', y_ax, z_ax, 'Y (mm)', 'Z (mm)', 'Sagittal  (mean-proj along X)';
        dose_axi', sens_axi', x_ax, y_ax, 'X (mm)', 'Y (mm)', 'Axial  (mean-proj along Z)';
    };

    dose_color   = [0.20, 0.50, 0.90];
    sensor_color = [0.90, 0.10, 0.10];

    for col = 1:3
        ax  = subplot(1, 3, col);
        d2d = double(view_data{col, 1});
        s2d = double(view_data{col, 2});
        xv  = view_data{col, 3};
        yv  = view_data{col, 4};

        if have_density
            dn     = (ct_projs{col} - wl_min) / wl_width;
            dn     = max(0, min(1, dn));
            bg_rgb = repmat(dn, [1, 1, 3]);
        else
            bg_rgb = ones([size(d2d), 3]);
        end
        image(ax, xv, yv, bg_rgb);
        hold(ax, 'on');

        if any(d2d(:))
            dose_rgb = repmat(reshape(dose_color, [1,1,3]), [size(d2d), 1]);
            h_dose   = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = d2d * 0.45;
        end
        if any(s2d(:))
            sens_rgb = repmat(reshape(sensor_color, [1,1,3]), [size(s2d), 1]);
            h_sens   = image(ax, xv, yv, sens_rgb);
            h_sens.AlphaData = s2d * 0.85;
        end

        hold(ax, 'off');
        axis(ax, 'xy'); axis(ax, 'image');
        xlabel(ax, view_data{col, 5}); ylabel(ax, view_data{col, 6});
        title(ax, view_data{col, 7}, 'FontSize', 10);

        hold(ax, 'on');
        p1 = patch(ax, NaN, NaN, dose_color,   'FaceAlpha', 0.45, 'EdgeColor', 'none');
        p2 = patch(ax, NaN, NaN, sensor_color, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
        hold(ax, 'off');
        if col == 3
            legend(ax, [p1, p2], {'Dose mask (>=10%)', 'Sensor'}, ...
                'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
        end
    end
    drawnow;
end


function plot_dose_panels(original, recon, sensor_mask, density, spacing_mm, titleStr, rowLabels, smooth_sigma)
%PLOT_DOSE_PANELS  2x3 dose comparison: coronal, sagittal, axial.
%  Row 1 = listed-dose recon, Row 2 = counterpart-dose recon. Dose (jet,
%  semi-transparent) over density (grayscale). Sensor contour in red. Both rows
%  share the clim of `original`. smooth_sigma (voxels) is DISPLAY-only smoothing.
    if nargin < 8 || isempty(smooth_sigma), smooth_sigma = 0; end

    gridSize = size(original);
    if ~isequal(size(sensor_mask), gridSize)
        sensor_mask = sensor_mask(1:gridSize(1), 1:gridSize(2), 1:gridSize(3));
    end
    have_density = ~isempty(density) && isequal(size(density), gridSize);

    [~, max_idx] = max(original(:));
    [cx, cy, cz] = ind2sub(gridSize, max_idx);

    max_dose = max(original(:));
    if max_dose == 0, max_dose = 1; end

    x_ax = (1:gridSize(1)) * spacing_mm(1);
    y_ax = (1:gridSize(2)) * spacing_mm(2);
    z_ax = (1:gridSize(3)) * spacing_mm(3);

    cmap_jet = jet(256);
    wl_center = 1050; wl_width = 350;
    wl_min    = wl_center - wl_width / 2;

    figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
        'Position', [50, 50, 1380, 700]);
    sgtitle(sprintf('%s\nIsocenter (max dose): X=%d  Y=%d  Z=%d voxel  |  Dose clim [0, %.4f] Gy', ...
        titleStr, cx, cy, cz, max_dose), 'FontWeight', 'bold', 'FontSize', 11);

    if nargin < 7 || isempty(rowLabels)
        row_labels = {'Original', 'Reconstructed'};
    else
        row_labels = rowLabels;
    end
    doses = {original, recon};

    for row = 1:2
        d   = doses{row};
        lbl = row_labels{row};

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
            xv = xvecs{col};
            yv = yvecs{col};

            if have_density
                dn     = (density_slices{col} - wl_min) / wl_width;
                dn     = max(0, min(1, dn));
                bg_rgb = repmat(dn, [1, 1, 3]);
            else
                bg_rgb = zeros([size(dose_slice), 3]);
            end
            image(ax, xv, yv, bg_rgb);
            hold(ax, 'on');

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

            s = sensor_slices{col};
            if any(s(:))
                contour(ax, xv, yv, s, [0.5, 0.5], 'r-', 'LineWidth', 2);
            end
            hold(ax, 'off');

            axis(ax, 'xy'); axis(ax, 'image');
            colormap(ax, cmap_jet);
            caxis(ax, [0, max_dose]);
            cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';

            xlabel(ax, xlbls{col}); ylabel(ax, ylbls{col});
            title(ax, sprintf('%s  %s', lbl, tsuffs{col}));
        end
    end
    drawnow;
end


function plot_gamma_and_error_axial(gamma_results, original, recon, sensor_mask, cz, titleStr)
%PLOT_GAMMA_AND_ERROR_AXIAL  1x(nCrit+2) axial figure:
%  gamma map(s) | absolute dose error | normalized reconstruction.
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
    nCol       = nCrit + 2;

    figure('Name', ['Gamma / Error / Normalized Recon - ' titleStr], 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [50, 300, 360 * nCol, 370]);
    sgtitle(sprintf('%s\nAxial Plane (Z = %d voxel)  Gamma Index, |Target - Reference|, Normalized Recon', titleStr, cz), ...
        'FontWeight', 'bold', 'FontSize', 11);

    gamma_clim   = [0, 2];
    sensor_slice = squeeze(sensor_mask(:, :, cz))';

    for g = 1:nCrit
        ax = subplot(1, nCol, g);
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

    orig_slice  = squeeze(original(:, :, cz))';
    recon_slice = squeeze(recon(:, :, cz))';

    ax = subplot(1, nCol, nCrit + 1);
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

    ax = subplot(1, nCol, nCrit + 2);
    recon_max = max(recon_slice(:));
    if recon_max == 0, recon_max = 1; end
    imagesc(ax, recon_slice, [0, recon_max]);
    axis(ax, 'xy'); axis(ax, 'image');
    colormap(ax, 'parula');
    cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
    hold(ax, 'on');
    if any(sensor_slice(:))
        contour(ax, sensor_slice, [0.5, 0.5], 'r-', 'LineWidth', 2);
    end
    hold(ax, 'off');
    xlabel(ax, 'X (voxel)'); ylabel(ax, 'Y (voxel)');
    title(ax, sprintf('Normalized Recon\nMax: %.4f Gy', recon_max));

    drawnow;
end


function cmap = gamma_colormap_gyr()
%GAMMA_COLORMAP_GYR  Green-yellow-red colormap for gamma index.
    n    = 256;
    half = round(n / 2);
    rest = n - half;
    r = [linspace(0, 1, half)'; ones(rest, 1)];
    g = [linspace(0.85, 1, half)'; linspace(1, 0, rest)'];
    b = zeros(n, 1);
    cmap = [r, g, b];
end


function plot_convergence_history(conv_max_pressure, conv_rel_change, num_iters, tol, p0_max_orig, tag)
%PLOT_CONVERGENCE_HISTORY  p0 convergence over TR iterations.
    if nargin < 6, tag = ''; end
    iters   = 1:num_iters;
    p_vals  = conv_max_pressure(iters);
    rc_vals = conv_rel_change(iters);

    figure('Name', ['p0 Convergence  ' tag], 'Color', 'w', ...
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
    ylabel('Relative Change ||p_n - p_{n-1}|| / ||p_{n-1}||', 'Color', [0.8, 0.1, 0.1]);

    xlabel('TR Iteration');
    title(sprintf('p_0 Convergence  (%d/%d iterations)   %s', num_iters, numel(conv_max_pressure), tag), ...
        'FontWeight', 'bold');
    xlim([0.5, num_iters + 0.5]);
    grid on;
    drawnow;
end
