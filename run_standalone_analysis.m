%% =========================================================================
%  RUN_STANDALONE_ANALYSIS.m
%  Analysis-only standalone pipeline  (NO k-Wave).
%
%  Reproduces the dose / gamma analysis and ALL plots of
%  run_standalone_simulation.m, run_standalone_comparison.m, and
%  study_gamma_pass_rates_chart.m WITHOUT running any forward simulation or
%  time-reversal reconstruction. Already-reconstructed doses are loaded from
%  disk via load_recon_dose_data.m (written by pipeline_simulate /
%  run_single_field_simulation), and only the cheap, deterministic pieces are
%  recomputed: tissue density (HU->tissue lookup), sensor placement
%  (determine_sensor_mask, pure geometry), gamma analysis, and the gamma
%  pass-rate sweep.
%
%  Omitted (require k-Wave and so cannot be reproduced here):
%    - p0 convergence-history plot (TR-iteration metric)
%    - noise-ensemble error bars on the pass-rate chart
%  =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path when this script
% is run from the repo root (pipeline flows already genpath PipelineScripts).
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% ========================= CONFIGURATION ================================

CONFIG.working_dir    = '/mnt/weka/home/80030361/ETHOS_Simulations';
CONFIG.patient_id     = '1194203';
CONFIG.session        = 'Session_1';
CONFIG.treatment_site = 'Pancreas';
CONFIG.gruneisen_method = 'threshold_2';

% --- Analysis mode ---
%   'single'     : single recon vs its RayStation truth (run_standalone_simulation)
%   'comparison' : recon_A vs recon_B, CT_1 vs CT_3 of one field (run_standalone_comparison)
%   'sweep'      : 'comparison' PLUS the gamma pass-rate-vs-criterion chart
%                  (study_gamma_pass_rates_chart, no noise ensemble)
CONFIG.analysis_mode = 'single';

% --- Field selection ---
% Leave beam/segment empty to auto-pick the first matched field ("any dose").
CONFIG.beam    = [];      % scalar beam number, or [] for first match
CONFIG.segment = [];      % scalar segment number, or [] for no segment filter
CONFIG.ct_pair = [1, 3];  % CT labels paired for comparison / sweep modes

% Gamma pass-rate sweep grid (n%/n mm). Used only in 'sweep' mode.
CONFIG.criterion_n = (1:0.5:5)';

% Explicit recon config-hash override ('' => auto-discover on disk).
CONFIG.config_hash = '';

% --- Tissue model (needed by create_medium for the density background) ---
CONFIG.force_uniform_density     = false;
CONFIG.force_uniform_sound_speed = false;
CONFIG.force_uniform_attenuation = false;
CONFIG.force_uniform_gruneisen   = false;

CONFIG.uniform_density      = 1000;
CONFIG.uniform_sound_speed  = 1540;
CONFIG.uniform_alpha_coeff  = 0;
CONFIG.uniform_alpha_power  = 1.1;
CONFIG.uniform_gruneisen    = 1.0;

% --- Sensor geometry (needed by determine_sensor_mask) ---
CONFIG.sensor_placement_method = 'determine_sensor_mask';
CONFIG.elements_per_side = 32;
CONFIG.element_pitch_mm  = 4.35;   % Kerf derived inside determine_sensor_mask as (pitch - size)
CONFIG.element_size_mm   = 3.65;
CONFIG.sensor_standoff_mm = 5;
CONFIG.jaw_margin_mm      = 10;
CONFIG.sensor_placement   = 'anterior';
CONFIG.pml_size           = 10;
CONFIG.aim_at_iso         = true;
CONFIG.force_turn_angle   = 290;

% --- Normalization ---
% Least-squares relative normalization (method from study_pass_rates_allsegments.m):
% each reconstruction is rescaled by the scalar gain that best fits it to its
% OWN RayStation truth over that truth's 10% low-dose region, before the gamma
% analysis. RS truth volumes are left in absolute Gy. Enabled by default.
CONFIG.normalize = true;

% --- Visualization ---
% Gaussian sigma (voxels) applied to dose slices for DISPLAY ONLY in the dose
% comparison panels. Fills speckle gaps in a spotty recon. 0 disables.
CONFIG.viz_smooth_sigma = 1.0;

%% ========================= PRINT CONFIGURATION ===========================

fprintf('=========================================================\n');
fprintf('  Standalone ANALYSIS  (no k-Wave)\n');
fprintf('=========================================================\n');
fprintf('  Patient:        %s / %s\n', CONFIG.patient_id, CONFIG.session);
fprintf('  Mode:           %s\n', CONFIG.analysis_mode);
fprintf('  Tissue model:   %s\n', CONFIG.gruneisen_method);
fprintf('  Sensor:         %s\n', CONFIG.sensor_placement_method);
fprintf('=========================================================\n');

%% ========================= LOAD RECON DOSES =============================
% Reference (ref) and target (tgt) doses depend on the mode:
%   single      ref = RayStation truth, tgt = reconstructed dose
%   comparison  ref = CT_1 recon,       tgt = CT_3 recon
%   sweep       same as comparison, plus the pass-rate chart

switch lower(CONFIG.analysis_mode)
    case 'single'
        % Pick the first field (or the CONFIG-specified beam/segment) and load
        % it via Mode='set' (no single-match requirement). A given beam/segment
        % may exist for both CT_1 and CT_3; for single mode we just take the
        % first returned field ("any dose").
        if isempty(CONFIG.beam)
            [beam, seg] = pick_first_field(CONFIG);
        else
            beam = CONFIG.beam;
            seg  = CONFIG.segment;
        end

        out = load_field_set(CONFIG, beam, seg);
        fld = out.fields(1);

        refDose      = double(fld.rs_dose);
        tgtDose      = double(fld.recon_dose);
        cbct         = fld.cbct;
        metadata     = out.metadata;
        gantry_angle = getfield_def(fld.rtplan, 'gantry_angle', 0);
        labelRef     = 'RayStation (truth)';
        labelTgt     = sprintf('Reconstructed (%s)', field_ct_label(fld));
        panelTitle   = sprintf('Dose Comparison: %s vs %s', labelRef, labelTgt);

        % Least-squares normalization: scale the recon (target) to its RS truth
        % (the reference); the truth is left unscaled. (See CONFIG.normalize.)
        if isfield(CONFIG, 'normalize') && CONFIG.normalize
            gT = least_squares_gain(refDose, tgtDose);
            tgtDose = tgtDose * gT;
            fprintf('  [Normalize] LS gain recon->truth: x%.4g\n', gT);
        end

    case {'comparison', 'sweep'}
        ct1 = sprintf('CT_%d', CONFIG.ct_pair(1));
        ct2 = sprintf('CT_%d', CONFIG.ct_pair(2));

        % Find a beam/segment that has BOTH CT recons, then load that set and
        % select each CT field by its own ct_label (robust to filename quirks).
        if isempty(CONFIG.beam)
            [beam, seg] = pick_first_paired_field(CONFIG, ct1, ct2);
        else
            beam = CONFIG.beam;
            seg  = CONFIG.segment;
        end

        out = load_field_set(CONFIG, beam, seg);
        iA  = find_field_by_ct(out.fields, ct1);
        iB  = find_field_by_ct(out.fields, ct2);
        if isempty(iA) || isempty(iB)
            found = strjoin(arrayfun(@(f) field_ct_label(f), out.fields, ...
                'UniformOutput', false), ', ');
            error('run_standalone_analysis:NoPair', ...
                ['Beam %s / segment %s does not have both %s and %s recon ', ...
                 'doses (found CT labels: %s). Set CONFIG.beam/segment to a ', ...
                 'field reconstructed on both CTs.'], ...
                mat2str(beam), mat2str(seg), ct1, ct2, found);
        end
        fldA = out.fields(iA);
        fldB = out.fields(iB);

        refDose      = double(fldA.recon_dose);
        tgtDose      = double(fldB.recon_dose);
        cbct         = fldA.cbct;
        metadata     = out.metadata;
        gantry_angle = getfield_def(fldA.rtplan, 'gantry_angle', 0);
        labelRef     = sprintf('%s recon', ct1);
        labelTgt     = sprintf('%s recon', ct2);
        panelTitle   = sprintf('Recon vs Recon: %s vs %s', labelRef, labelTgt);

        % Least-squares normalization: scale EACH recon to its OWN-CT RS truth
        % over that truth's 10% region. (See CONFIG.normalize.)
        if isfield(CONFIG, 'normalize') && CONFIG.normalize
            gA = least_squares_gain(double(fldA.rs_dose), refDose);
            gB = least_squares_gain(double(fldB.rs_dose), tgtDose);
            refDose = refDose * gA;
            tgtDose = tgtDose * gB;
            fprintf('  [Normalize] LS gains: %s x%.4g, %s x%.4g (each -> own truth)\n', ...
                ct1, gA, ct2, gB);
        end

    otherwise
        error('run_standalone_analysis:BadMode', ...
            'CONFIG.analysis_mode must be single | comparison | sweep.');
end

spacing_mm = metadata.spacing(:)';
gridSize   = size(refDose);
fprintf('  Grid: [%d x %d x %d]  spacing [%.2f %.2f %.2f] mm\n', ...
    gridSize(1), gridSize(2), gridSize(3), spacing_mm(1), spacing_mm(2), spacing_mm(3));

beam_meta = [];
if isfield(metadata, 'beam_metadata') && ~isempty(metadata.beam_metadata)
    beam_meta = metadata.beam_metadata;
end

%% ========================= TISSUE DENSITY ===============================
% Density background for the dose / sensor panels. Pure HU->tissue lookup
% (no k-Wave). Falls back to [] (black/white background) if CBCT is missing.

density = [];
if ~isempty(cbct) && isfield(cbct, 'cubeHU')
    sct = cbct;
    sct.spacing = spacing_mm;
    medium  = create_medium(sct, CONFIG);
    density = medium.density;
    if ~isequal(size(density), gridSize)
        fprintf('  [WARN] Density grid %s != dose grid %s; dropping density bg.\n', ...
            mat2str(size(density)), mat2str(gridSize));
        density = [];
    end
else
    fprintf('  [WARN] No CBCT cubeHU available; plotting without density background.\n');
end

%% ============= SUMMED PLAN DOSE FOR SENSOR PLACEMENT ====================
% Sensor placement uses the step15 SUMMED PLAN DOSE (total_rs_dose.mat) so the
% exclusion zone and array center are DETERMINISTIC and identical across every
% field/segment, instead of shifting with whichever single recon is loaded.
% determine_sensor_mask reads this path via config.total_dose_file.
processed_dir   = fullfile(CONFIG.working_dir, 'RayStationFiles', ...
    CONFIG.patient_id, CONFIG.session, 'processed');
total_dose_file = fullfile(processed_dir, 'total_rs_dose.mat');
if isfile(total_dose_file)
    CONFIG.total_dose_file = total_dose_file;
    fprintf('  Sensor placement dose: %s (summed plan dose)\n', total_dose_file);
else
    CONFIG.total_dose_file = '';
    fprintf(['  [WARN] total_rs_dose.mat not found in %s;\n' ...
             '         sensor placement falls back to the per-field dose.\n'], ...
             processed_dir);
end

%% ========================= SENSOR PLACEMENT =============================
% Geometric placement only (no k-Wave). Cropped to the dose grid for overlay.

sensor_mask = compute_sensor_mask(cbct, refDose, spacing_mm, gantry_angle, ...
    beam_meta, CONFIG);
fprintf('  Sensor voxels (on dose grid): %d\n', sum(sensor_mask(:)));

%% ========================= GAMMA ANALYSIS ===============================

gamma_results = compute_gamma_results(refDose, tgtDose, spacing_mm);

%% ========================= PLOTS =======================================

% Axial slice index at the reference max-dose voxel.
[~, max_idx] = max(refDose(:));
[~, ~, cz_gamma] = ind2sub(gridSize, max_idx);

% 2x3 dose panels (original/reference vs recon/target).
plot_dose_panels(refDose, tgtDose, sensor_mask, density, spacing_mm, ...
    panelTitle, {labelRef, labelTgt}, CONFIG.viz_smooth_sigma);

% 1x3 sensor-vs-dose anatomical view.
dose_mask = refDose >= 0.10 * max(refDose(:));
plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, density, CONFIG);

% 1x4 gamma + absolute-error axial figure.
if ~isempty(gamma_results)
    plot_gamma_and_error_axial(gamma_results, refDose, tgtDose, sensor_mask, cz_gamma);
end

% Gamma pass-rate vs criterion chart (sweep mode only; no ensemble bars).
if strcmpi(CONFIG.analysis_mode, 'sweep')
    fprintf('\n[Sweep] Computing gamma pass-rate sweep...\n');
    sweep_pr = gamma_sweep_pass_rates(refDose, tgtDose, CONFIG.criterion_n, spacing_mm);
    plot_gamma_pass_rate_curve(CONFIG.criterion_n, sweep_pr, labelRef, labelTgt, ...
        CONFIG.patient_id, [], []);
end

fprintf('\n[run_standalone_analysis] Done.\n');


%% =========================================================================
%  LOADING HELPERS
%% =========================================================================

function out = load_field_set(CONFIG, beam, seg)
%LOAD_FIELD_SET Load matched fields via load_recon_dose_data (Mode='set').
%  'set' does NOT require a unique match, so a beam/segment present for several
%  CTs / plan types loads them all; the caller selects which to use.
    args = {'Mode', 'set', 'IncludeEthos', false};
    if ~isempty(beam), args = [args, {'Beam', beam}]; end
    if ~isempty(seg),  args = [args, {'Segment', seg}]; end
    if isfield(CONFIG, 'config_hash') && ~isempty(CONFIG.config_hash)
        args = [args, {'Hash', CONFIG.config_hash}];
    end
    out = load_recon_dose_data(CONFIG.patient_id, CONFIG.session, CONFIG, args{:});
end

function [beam, seg] = pick_first_field(CONFIG)
%PICK_FIRST_FIELD Beam/segment of the first indexed processed field dose.
    ldcfg.working_dir = CONFIG.working_dir;
    fidx = list_processed_field_doses(CONFIG.patient_id, CONFIG.session, ldcfg);
    beam = fidx(1).beam_index;
    seg  = fidx(1).segment;
end

function [beam, seg] = pick_first_paired_field(CONFIG, ct1, ct2)
%PICK_FIRST_PAIRED_FIELD First beam/segment present on BOTH ct1 and ct2.
%  Uses the lightweight file index (no array loads) and a lenient CT parse, so
%  it does not depend on the strict dose-filename regex.
    ldcfg.working_dir = CONFIG.working_dir;
    fidx = list_processed_field_doses(CONFIG.patient_id, CONFIG.session, ldcfg);

    n = numel(fidx);
    keys = strings(1, n);          % "beam|seg" identity per file
    cts  = strings(1, n);          % CT label per file
    for i = 1:n
        keys(i) = sprintf('%d|%d', fidx(i).beam_index, fidx(i).segment);
        cts(i)  = string(ct_label_from_name(fidx(i).source_mat_filename));
    end

    seen = unique(keys, 'stable');
    for k = 1:numel(seen)
        m = keys == seen(k);
        if any(cts(m) == ct1) && any(cts(m) == ct2)
            idx  = find(m, 1);
            beam = fidx(idx).beam_index;
            seg  = fidx(idx).segment;
            return;
        end
    end

    error('run_standalone_analysis:NoPairedField', ...
        ['No field is reconstructed on both %s and %s. Set CONFIG.beam / ', ...
         'CONFIG.segment, or adjust CONFIG.ct_pair.'], ct1, ct2);
end

function idx = find_field_by_ct(fields, want_ct)
%FIND_FIELD_BY_CT Index of the first loaded field whose CT label is want_ct.
    idx = [];
    for i = 1:numel(fields)
        if strcmpi(field_ct_label(fields(i)), want_ct)
            idx = i;
            return;
        end
    end
end

function lbl = field_ct_label(fld)
%FIELD_CT_LABEL CT label of a loaded field: rtplan.ct_label, else filename.
    lbl = '';
    if isfield(fld, 'rtplan') && isfield(fld.rtplan, 'ct_label') ...
            && ~isempty(fld.rtplan.ct_label)
        lbl = strrep(char(fld.rtplan.ct_label), '-', '_');
    end
    if isempty(lbl) && isfield(fld, 'source_mat_filename')
        lbl = ct_label_from_name(fld.source_mat_filename);
    end
    if isempty(lbl), lbl = 'unknown'; end
end

function ct_label = ct_label_from_name(name)
%CT_LABEL_FROM_NAME Lenient CT-label parse: matches CT_1, CT-1, CT1, ...
    ct_label = '';
    tok = regexp(char(name), 'CT[_-]?(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        ct_label = sprintf('CT_%s', tok{1});
    end
end

function v = getfield_def(s, f, default)
%GETFIELD_DEF Struct field value, or default when missing/empty/NaN.
    v = default;
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        candidate = s.(f);
        if isnumeric(candidate) && all(isnan(candidate(:)))
            return;
        end
        v = candidate;
    end
end


%% =========================================================================
%  SENSOR PLACEMENT (geometric; no k-Wave)
%% =========================================================================

function sensor_mask = compute_sensor_mask(cbct, refDose, spacing_mm, ...
        gantry_angle, beam_meta, CONFIG)
%COMPUTE_SENSOR_MASK Place the 2D array via determine_sensor_mask and crop to
%  the dose grid. The k-Wave grid-expansion/padding logic is intentionally
%  skipped (only relevant to the forward sim); for an overlay we just crop any
%  expansion back to the dose-grid dimensions. Falls back to an empty mask on
%  any failure so plots still render.

    gridSize    = size(refDose);
    sensor_mask = false(gridSize);

    if ~strcmpi(CONFIG.sensor_placement_method, 'determine_sensor_mask')
        return;
    end
    if isempty(cbct) || ~isfield(cbct, 'bodyMask')
        fprintf('  [WARN] No CBCT bodyMask; skipping sensor placement.\n');
        return;
    end

    try
        sct_for_sensor = cbct;
        if ~isfield(sct_for_sensor, 'couchMask') || isempty(sct_for_sensor.couchMask)
            sct_for_sensor.couchMask = false(size(sct_for_sensor.bodyMask));
        end
        if ~isfield(sct_for_sensor, 'origin') || isempty(sct_for_sensor.origin)
            sct_for_sensor.origin = [0, 0, 0];
        end
        sct_for_sensor.spacing = spacing_mm;

        field_dose_for_sensor             = struct();
        field_dose_for_sensor.dose_Gy     = refDose;
        field_dose_for_sensor.gantry_angle = gantry_angle;
        field_dose_for_sensor.origin      = sct_for_sensor.origin;
        field_dose_for_sensor.spacing     = spacing_mm;
        field_dose_for_sensor.dimensions  = [gridSize(1), gridSize(2), gridSize(3)];

        sm = determine_sensor_mask(sct_for_sensor, field_dose_for_sensor, ...
            beam_meta, CONFIG);

        m1 = min(size(sm, 1), gridSize(1));
        m2 = min(size(sm, 2), gridSize(2));
        m3 = min(size(sm, 3), gridSize(3));
        sensor_mask(1:m1, 1:m2, 1:m3) = logical(sm(1:m1, 1:m2, 1:m3));
    catch ME
        fprintf('  [WARN] Sensor placement failed (%s); no sensor overlay.\n', ...
            ME.message);
        sensor_mask = false(gridSize);
    end
end


%% =========================================================================
%  GAMMA ANALYSIS  (copied from run_standalone_simulation.m:1270-1330)
%% =========================================================================

function g = least_squares_gain(rs_truth, recon)
%LEAST_SQUARES_GAIN Scalar gain aligning a recon to its RS truth (relative norm).
%  g = sum(rs.*recon)/sum(recon.^2) over the truth's 10% low-dose region; the g
%  minimizing ||rs - g*recon||^2 there. Falls back to 1 for empty/zero inputs.
%  (Method mirrored from study_pass_rates_allsegments.m.)
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


function gamma_results = compute_gamma_results(refDose, tgtDose, spacing_mm)
    gamma_results = struct();

    if exist('CalcGamma', 'file') ~= 2
        warning('CalcGamma not found. Skipping gamma analysis.');
        gamma_results = [];
        return;
    end

    fprintf('\n[Gamma] Running gamma analysis...\n');

    ref_struct.start = [0, 0, 0];
    ref_struct.width = spacing_mm;
    ref_struct.data  = double(refDose);

    tgt_struct.start = [0, 0, 0];
    tgt_struct.width = spacing_mm;
    tgt_struct.data  = double(tgtDose);

    low_dose_cutoff = 0.10 * max(refDose(:));
    gamma_eval_mask = refDose >= low_dose_cutoff;

    gamma_criteria = {10, 10, '10%/10mm'; 5, 5, '5%/5mm'; 3, 3, '3%/3mm'};
    gamma_maps     = cell(size(gamma_criteria, 1), 1);
    pass_rates     = zeros(size(gamma_criteria, 1), 1);

    for gc = 1:size(gamma_criteria, 1)
        pct_crit = gamma_criteria{gc, 1};
        dta_crit = gamma_criteria{gc, 2};
        lbl      = gamma_criteria{gc, 3};

        fprintf('  [%s] Computing...', lbl);
        try
            gmap = CalcGamma(ref_struct, tgt_struct, pct_crit, dta_crit, ...
                'local', 0, 'limit', dta_crit * 2, 'restrict', 1);
            gamma_maps{gc} = gmap;
            eval_vals       = gmap(gamma_eval_mask);
            pass_rate       = 100 * mean(eval_vals <= 1);
            pass_rates(gc)  = pass_rate;
            fprintf('  Pass rate: %.2f%%\n', pass_rate);
        catch ME
            warning('Gamma [%s] failed: %s', lbl, ME.message);
            gamma_maps{gc} = [];
            pass_rates(gc) = NaN;
            fprintf('  FAILED\n');
        end
    end

    gamma_results.maps        = gamma_maps;
    gamma_results.pass_rates  = pass_rates;
    gamma_results.criteria    = gamma_criteria;
    gamma_results.cutoff_Gy   = low_dose_cutoff;
    gamma_results.eval_mask   = gamma_eval_mask;

    fprintf('\n  ------ Gamma Pass Rates (10%% low-dose cutoff) ------\n');
    for gc = 1:size(gamma_criteria, 1)
        if isnan(pass_rates(gc))
            fprintf('  %-12s  FAILED\n', gamma_criteria{gc, 3});
        else
            fprintf('  %-12s  %.2f%%\n', gamma_criteria{gc, 3}, pass_rates(gc));
        end
    end
end


%% =========================================================================
%  PLOTTING  (copied from run_standalone_comparison.m / study_*.m)
%% =========================================================================

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


function plot_sensor_dose_planes(dose_mask, sensor_mask, spacing_mm, density, config)
%PLOT_SENSOR_DOSE_PLANES  1x3 anatomical view of sensor geometry vs dose mask.
%  Shows three orthogonal projections (coronal, sagittal, axial).
%  CT density is rendered as a grayscale background (mean-projection).
%  Dose mask (dose >= 10% max) drawn as a filled semi-transparent blue region.
%  Sensor drawn as a solid red line/region  computed via max-projection so it
%  always appears regardless of which slice the dose centroid falls on.

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

    figure('Name', 'Gamma & Absolute Error  Axial', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [50, 300, 1400, 370]);
    sgtitle(sprintf('Axial Plane (Z = %d voxel)  Gamma Index & Absolute Error', cz), ...
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


function plot_gamma_pass_rate_curve(crit_n, pass_rates, label_ref, label_tgt, patient_id, ens_mean, ens_std)
%PLOT_GAMMA_PASS_RATE_CURVE  Pass rate vs gamma criterion (n%/n mm).
%  crit_n      column vector of criterion values
%  pass_rates  nominal (single-draw) pass rates in percent (NaN entries skipped)
%  ens_mean,ens_std  optional noise-ensemble mean and std per criterion. When
%                    provided (non-empty), error bars are overlaid as the
%                    primary series and the nominal curve is shown faintly.
%  A red horizontal reference line is drawn at the 90% pass rate.

    if nargin < 6, ens_mean = []; end
    if nargin < 7, ens_std  = []; end

    crit_n     = crit_n(:);
    pass_rates = pass_rates(:);
    have_ens   = ~isempty(ens_mean) && ~isempty(ens_std);

    figure('Name', 'Gamma Pass Rate vs Criterion', 'Color', 'w', ...
        'NumberTitle', 'off', 'Position', [120, 120, 760, 470]);
    hold on;

    legend_handles = gobjects(0);
    legend_labels  = {};

    valid = ~isnan(pass_rates);
    if have_ens
        % Nominal single-draw curve shown faintly for reference.
        h_nom = plot(crit_n(valid), pass_rates(valid), '--o', ...
            'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.0, 'MarkerSize', 4);
        legend_handles(end+1) = h_nom; %#ok<AGROW>
        legend_labels{end+1}  = 'Nominal (single draw)';

        % Ensemble mean +/- std error bars (primary series).
        em = ens_mean(:); es = ens_std(:);
        vv = ~isnan(em);
        h_ens = errorbar(crit_n(vv), em(vv), es(vv), 'b-o', 'LineWidth', 1.8, ...
            'MarkerSize', 6, 'MarkerFaceColor', [0.2, 0.4, 1.0], 'CapSize', 8);
        legend_handles(end+1) = h_ens; %#ok<AGROW>
        legend_labels{end+1}  = 'Ensemble mean \pm std';

        % Annotate ensemble points with mean +/- std
        for k = 1:numel(crit_n)
            if vv(k)
                text(crit_n(k), em(k) + es(k) + 1.5, ...
                    sprintf('%.1f\\pm%.1f', em(k), es(k)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    else
        h_nom = plot(crit_n(valid), pass_rates(valid), 'b-o', 'LineWidth', 1.8, ...
            'MarkerSize', 6, 'MarkerFaceColor', [0.2, 0.4, 1.0]);
        legend_handles(end+1) = h_nom; %#ok<AGROW>
        legend_labels{end+1}  = 'Pass rate';

        for k = 1:numel(crit_n)
            if valid(k)
                text(crit_n(k), pass_rates(k) + 1.5, sprintf('%.1f%%', pass_rates(k)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    end

    % 90% pass-rate reference line (red)
    yline(90, 'r-', '90% pass rate', 'LineWidth', 1.8, ...
        'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'bottom', 'FontSize', 9);

    hold off;
    grid on;
    xlim([min(crit_n) - 0.5, max(crit_n) + 0.5]);
    ylim([0, 105]);
    xticks(crit_n);
    xlabel('Gamma criterion  (n%/n mm)');
    ylabel('Gamma pass rate (%)');
    title(sprintf('Gamma Pass Rate vs Criterion   |   Patient %s\n%s vs %s', ...
        patient_id, label_ref, label_tgt), 'FontWeight', 'bold', 'FontSize', 11);
    legend(legend_handles, legend_labels, 'Location', 'southeast', 'FontSize', 8);
    drawnow;
end


function pr = gamma_sweep_pass_rates(reconA, reconB, crit_vec, spacing_mm)
%GAMMA_SWEEP_PASS_RATES  Pass rate for each n%/n mm criterion between two doses.
%  Serial over criteria. Reference = reconA, target = reconB, with a 10%
%  low-dose cutoff on the reference.
    ref_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(reconA));
    tgt_struct = struct('start', [0, 0, 0], 'width', spacing_mm, 'data', double(reconB));

    cutoff    = 0.10 * max(reconA(:));
    eval_mask = reconA >= cutoff;

    pr = nan(numel(crit_vec), 1);
    for k = 1:numel(crit_vec)
        c = crit_vec(k);
        try
            gmap  = CalcGamma(ref_struct, tgt_struct, c, c, ...
                'local', 0, 'limit', c * 2, 'restrict', 1);
            pr(k) = 100 * mean(gmap(eval_mask) <= 1);
        catch
            pr(k) = NaN;
        end
    end
end


%% =========================================================================
%  TISSUE MEDIUM  (copied from run_standalone_simulation.m:1904-1998)
%% =========================================================================

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
end
