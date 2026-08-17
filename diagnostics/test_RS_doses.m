%% =========================================================================
%  TEST_RS_DOSES.m
%  ETHOS Photoacoustic Pipeline — RayStation Dose Verification
%  =========================================================================
%
%  PURPOSE:
%  Verifies that RayStation-computed field doses agree with the ETHOS system
%  truth dose (RTDOSE DICOM export).  Run this after pipeline_compress.m has
%  successfully processed a session.
%
%  The ETHOS RTDOSE is treated as the reference.  Field doses from
%  step15_process_doses are split by CT label (CT_1 vs CT_3), summed
%  separately, then each sum is compared to the truth via global 3%/3mm
%  gamma analysis.  Some disagreement is expected because RayStation
%  recomputes dose on different CT tables; large failures indicate a
%  miscalibration or export error.
%
%  INPUTS (via CONFIG struct below):
%    patient_id         - Patient identifier string (e.g. '1194203')
%    session            - Session string (e.g. 'Session_1')
%    treatment_site     - Treatment site name matching EthosExports tree
%    working_dir        - Repository root (default: directory of this file)
%    gamma_dose_pct     - Dose difference criterion (%, default 3)
%    gamma_dist_mm      - Distance-to-agreement criterion (mm, default 3)
%    gamma_dose_cutoff_pct - Low-dose exclusion threshold (% of max, default 10)
%    save_figures       - Save PNG figures to figure_dir (default true)
%    figure_dir         - Output directory; '' = auto-resolve under AnalysisResults
%
%  OUTPUTS:
%    Console summary table with pass rate / mean gamma / max gamma
%    PNG figures (when save_figures = true):
%      dose_comparison.png   — orthogonal views of ETHOS, CT_1 sum, CT_3 sum
%      dose_difference.png   — CT_1−ETHOS and CT_3−ETHOS difference maps
%      gamma_maps.png        — gamma index maps with pass rates in title
%
%  ALGORITHM:
%    0. Initialise, addpath
%    1. Load ETHOS truth dose from RTDOSE*.dcm in sct directory
%    2. Load processed field doses via load_processed_data
%    3. Sum field doses by ct_label (CT_1 / CT_3)
%    4. Resample ETHOS to RS dose grid if dimensions differ
%    5. Run gamma analysis for each available CT sum vs ETHOS
%    6. Print summary table
%    7. Generate and optionally save figures
%    8. Report elapsed time
%
%  EXAMPLE:
%    % Edit the CONFIG block below, then run:
%    test_RS_doses
%
%  DEPENDENCIES:
%    load_processed_data.m
%    CalcGamma.m (Mark Geurts, University of Wisconsin)
%    Image Processing Toolbox (imresize3)
%
%  PREREQUISITES:
%    pipeline_compress.m must have run successfully for the same patient/session.
%
%  AUTHOR: ETHOS Pipeline Team
%  DATE: May 2026
%  =========================================================================

clear; clc; close all;

% Ensure the moved helper functions in utils/ are on the path when run from root.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% ========================= CONFIGURATION =================================

CONFIG.patient_id            = '1194203';
CONFIG.session               = 'Session_1';
CONFIG.treatment_site        = 'Pancreas';
CONFIG.working_dir           = fileparts(mfilename('fullpath'));
CONFIG.gamma_dose_pct        = 3.0;
CONFIG.gamma_dist_mm         = 3.0;
CONFIG.gamma_dose_cutoff_pct = 10.0;
CONFIG.save_figures          = true;
CONFIG.figure_dir            = '';   % '' = AnalysisResults/<id>/<session>/rs_verification/
CONFIG.recompute_sums        = false; % true = re-accumulate from field doses even if cache exists
CONFIG.n_fractions           = 10;   % number of fractions; accumulated dose is divided by this

%% ========================= STEP 0: INITIALISE ============================

script_timer = tic;
patient_id   = char(CONFIG.patient_id);
session      = char(CONFIG.session);

addpath(CONFIG.working_dir);

fprintf('=========================================================\n');
fprintf('  TEST_RS_DOSES — RayStation Dose Verification\n');
fprintf('=========================================================\n');
fprintf('  Patient:  %s\n', patient_id);
fprintf('  Session:  %s\n', session);
fprintf('  Started:  %s\n', datetime('now'));
fprintf('=========================================================\n\n');

%% ========================= STEP 1: ETHOS TRUTH DOSE ======================

fprintf('[STEP 1] Loading ETHOS truth dose...\n');

sct_dir  = fullfile(CONFIG.working_dir, 'EthosExports', patient_id, ...
                    CONFIG.treatment_site, session, 'sct');
rd_path = fullfile(sct_dir, 'RD_reference.dcm');

if ~isfile(rd_path)
    error('test_RS_doses:FileNotFound', ...
        'RD_reference.dcm not found in:\n  %s\nCheck CONFIG.treatment_site and that DICOM exports are present.', ...
        sct_dir);
end
dose_info = dicominfo(rd_path);
ethos_dose = double(squeeze(dicomread(rd_path)));

if isfield(dose_info, 'DoseGridScaling')
    ethos_dose = ethos_dose * dose_info.DoseGridScaling;
end

fprintf('  File:    %s\n', 'RD_reference.dcm');
fprintf('  Size:    [%d x %d x %d]\n', size(ethos_dose));
fprintf('  Max dose: %.4f Gy\n\n', max(ethos_dose(:)));

%% ========================= STEP 2 & 3: LOAD / ACCUMULATE DOSES ===========

processed_dir = fullfile(CONFIG.working_dir, 'RayStationFiles', patient_id, session, 'processed');
sums_cache    = fullfile(processed_dir, 'dose_sums_cache.mat');

cache_loaded = false;
if ~CONFIG.recompute_sums && isfile(sums_cache)
    fprintf('[STEP 2] Skipping load_processed_data — using cached dose sums.\n');
    fprintf('[STEP 3] Loading cached dose sums from:\n  %s\n', sums_cache);
    loaded_cache  = load(sums_cache, 'ct1_dose', 'ct3_dose', 'n_ct1', 'n_ct3', 'spacing');
    ct1_dose      = loaded_cache.ct1_dose;
    ct3_dose      = loaded_cache.ct3_dose;
    n_ct1         = loaded_cache.n_ct1;
    n_ct3         = loaded_cache.n_ct3;
    spacing       = loaded_cache.spacing;
    rs_size       = size(ct1_dose);
    sct_resampled = struct();
    cache_loaded  = true;
    fprintf('  CT_1 fields: %d  (max sum dose: %.4f Gy)\n', n_ct1, max(ct1_dose(:)));
    fprintf('  CT_3 fields: %d  (max sum dose: %.4f Gy)\n\n', n_ct3, max(ct3_dose(:)));
end

if ~cache_loaded
    fprintf('[STEP 2] Loading processed RayStation field doses...\n');

    [field_doses, sct_resampled, ~, dose_metadata] = ...
        load_processed_data(patient_id, session, CONFIG);

    rs_size = dose_metadata.dimensions;   % [ny, nx, nz]
    spacing = dose_metadata.spacing;      % [dx, dy, dz] mm

    fprintf('  RS dose grid: [%d x %d x %d], spacing [%.2f x %.2f x %.2f] mm\n\n', ...
        rs_size, spacing);

    fprintf('[STEP 3] Accumulating field doses by CT label...\n');

    ct1_dose = zeros(rs_size, 'double');
    ct3_dose = zeros(rs_size, 'double');
    n_ct1 = 0;  n_ct3 = 0;  n_other = 0;
    warned_other = false;

    for i = 1:numel(field_doses)
        fd = field_doses{i};
        if isempty(fd), continue; end

        label = upper(strtrim(fd.ct_label));

        switch label
            case 'CT_1'
                ct1_dose = ct1_dose + fd.dose_Gy;
                n_ct1 = n_ct1 + 1;
            case 'CT_3'
                ct3_dose = ct3_dose + fd.dose_Gy;
                n_ct3 = n_ct3 + 1;
            otherwise
                n_other = n_other + 1;
                if ~warned_other
                    warning('test_RS_doses:UnknownCTLabel', ...
                        'Field %d has ct_label ''%s'' (not CT_1 or CT_3); excluded from both sums.', ...
                        i, fd.ct_label);
                    warned_other = true;
                end
        end
    end

    fprintf('  CT_1 fields: %d  (max sum dose: %.4f Gy)\n', n_ct1, max(ct1_dose(:)));
    fprintf('  CT_3 fields: %d  (max sum dose: %.4f Gy)\n', n_ct3, max(ct3_dose(:)));
    if n_other > 0
        fprintf('  Other/unlabelled fields excluded: %d\n', n_other);
    end
    fprintf('\n');

    fprintf('  Saving dose sums cache to:\n  %s\n\n', sums_cache);
    save(sums_cache, 'ct1_dose', 'ct3_dose', 'n_ct1', 'n_ct3', 'spacing', '-v7.3');
end

% Scale to per-fraction dose
fprintf('  Scaling by 1/%d fractions...\n\n', CONFIG.n_fractions);
ct1_dose = ct1_dose / CONFIG.n_fractions;
ct3_dose = ct3_dose / CONFIG.n_fractions;

do_ct1 = (n_ct1 > 0);
do_ct3 = (n_ct3 > 0);

if ~do_ct1
    warning('test_RS_doses:NoCT1', 'No CT_1 field doses found — skipping CT_1 analysis.');
end
if ~do_ct3
    warning('test_RS_doses:NoCT3', 'No CT_3 field doses found — skipping CT_3 analysis.');
end

%% ========================= STEP 4: RESAMPLE ETHOS IF NEEDED ==============

fprintf('[STEP 4] Checking ETHOS/RS grid alignment...\n');

if ~isequal(size(ethos_dose), rs_size)
    fprintf('  [WARNING] ETHOS size [%d x %d x %d] differs from RS size [%d x %d x %d]\n', ...
        size(ethos_dose), rs_size);
    fprintf('  Resampling ETHOS to RS grid...\n');
    ethos_dose = rs_resample(ethos_dose, rs_size);
    fprintf('  Resampled. New max: %.4f Gy\n\n', max(ethos_dose(:)));
else
    fprintf('  Grids match. No resampling needed.\n\n');
end

%% ========================= STEP 5: GAMMA ANALYSIS ========================

fprintf('[STEP 5] Running gamma analysis (%.0f%%/%.0fmm)...\n', ...
    CONFIG.gamma_dose_pct, CONFIG.gamma_dist_mm);

gamma_ct1    = [];
gamma_ct3    = [];
gamma_ct1_ct3 = [];

if do_ct1
    fprintf('\n  --- ETHOS vs CT_1 ---\n');
    gamma_ct1 = rs_gamma(ethos_dose, ct1_dose, spacing, CONFIG);
end

if do_ct3
    fprintf('\n  --- ETHOS vs CT_3 ---\n');
    gamma_ct3 = rs_gamma(ethos_dose, ct3_dose, spacing, CONFIG);
end

if do_ct1 && do_ct3
    fprintf('\n  --- CT_1 vs CT_3 ---\n');
    gamma_ct1_ct3 = rs_gamma(ct1_dose, ct3_dose, spacing, CONFIG);
end

%% ========================= STEP 6: SUMMARY TABLE =========================

fprintf('\n=========================================================\n');
fprintf('  Dose Verification Summary (%.0f%%/%.0fmm global gamma)\n', ...
    CONFIG.gamma_dose_pct, CONFIG.gamma_dist_mm);
fprintf('=========================================================\n');
fprintf('  %-22s  %9s  %8s  %7s  %s\n', ...
    'Comparison', 'Pass Rate', 'Mean γ', 'Max γ', 'N evaluated');
fprintf('  %s\n', repmat('-', 1, 65));

if do_ct1
    fprintf('  %-22s  %8.1f%%  %8.3f  %7.3f  %d\n', ...
        'ETHOS vs CT_1', gamma_ct1.pass_rate, gamma_ct1.mean_gamma, ...
        gamma_ct1.max_gamma, gamma_ct1.num_evaluated);
else
    fprintf('  %-22s  %s\n', 'ETHOS vs CT_1', 'SKIPPED (no CT_1 fields found)');
end

if do_ct3
    fprintf('  %-22s  %8.1f%%  %8.3f  %7.3f  %d\n', ...
        'ETHOS vs CT_3', gamma_ct3.pass_rate, gamma_ct3.mean_gamma, ...
        gamma_ct3.max_gamma, gamma_ct3.num_evaluated);
else
    fprintf('  %-22s  %s\n', 'ETHOS vs CT_3', 'SKIPPED (no CT_3 fields found)');
end

if do_ct1 && do_ct3
    fprintf('  %-22s  %8.1f%%  %8.3f  %7.3f  %d\n', ...
        'CT_1 vs CT_3', gamma_ct1_ct3.pass_rate, gamma_ct1_ct3.mean_gamma, ...
        gamma_ct1_ct3.max_gamma, gamma_ct1_ct3.num_evaluated);
elseif ~do_ct1 || ~do_ct3
    fprintf('  %-22s  %s\n', 'CT_1 vs CT_3', 'SKIPPED (need both CT_1 and CT_3)');
end

fprintf('=========================================================\n\n');

%% ========================= STEP 7: FIGURES ===============================

if CONFIG.save_figures || nargout == 0

    fprintf('[STEP 7] Generating figures...\n');

    % Resolve figure output directory
    if isempty(CONFIG.figure_dir)
        fig_dir = fullfile(CONFIG.working_dir, 'AnalysisResults', ...
                           patient_id, session, 'rs_verification');
    else
        fig_dir = CONFIG.figure_dir;
    end

    do_save = CONFIG.save_figures;
    if do_save
        try
            if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
        catch ME_mkdir
            warning('test_RS_doses:MkdirFailed', ...
                'Could not create figure directory %s: %s', fig_dir, ME_mkdir.message);
            do_save = false;
        end
    end

    % Orthogonal slice indices at max ETHOS dose voxel
    [r0, c0, s0] = ind2sub(size(ethos_dose), find(ethos_dose == max(ethos_dose(:)), 1));

    body_mask = [];
    if isfield(sct_resampled, 'bodyMask')
        body_mask = logical(sct_resampled.bodyMask);
    end

    % Figure 1 — Dose comparison
    rs_fig1_dose_comparison(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir);

    % Figure 2 — Dose difference maps
    rs_fig2_dose_difference(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir);

    % Figure 3 — Gamma maps (ETHOS vs CT_1, ETHOS vs CT_3)
    rs_fig3_gamma_maps(gamma_ct1, gamma_ct3, do_ct1, do_ct3, ...
        r0, c0, s0, CONFIG, do_save, fig_dir);

    % Figure 4 — CT_1 vs CT_3 gamma map
    if do_ct1 && do_ct3
        rs_fig4_ct1_vs_ct3_gamma(gamma_ct1_ct3, ct1_dose, r0, c0, s0, CONFIG, do_save, fig_dir);
    end

    if do_save
        fprintf('  Figures saved to: %s\n\n', fig_dir);
    end
end

%% ========================= STEP 8: ELAPSED ===============================

fprintf('[DONE] test_RS_doses complete in %.1f sec\n', toc(script_timer));


%% =========================================================================
%  HELPER FUNCTIONS
%% =========================================================================

function resampled = rs_resample(dose, target_size)
% Trilinear resample of 3D dose to target_size.  Uses imresize3 when
% available (Image Processing Toolbox R2017a+), otherwise falls back to
% interp3.  Clamps negatives introduced by interpolation.

    if isequal(size(dose), target_size)
        resampled = dose;
        return;
    end

    if exist('imresize3', 'file')
        resampled = imresize3(double(dose), target_size, 'linear');
    else
        src = size(dose);
        [Xi, Yi, Zi] = ndgrid(...
            linspace(1, src(1), target_size(1)), ...
            linspace(1, src(2), target_size(2)), ...
            linspace(1, src(3), target_size(3)));
        resampled = interp3(double(dose), Yi, Xi, Zi, 'linear', 0);
    end

    resampled(resampled < 0) = 0;
end


function gamma = rs_gamma(reference, evaluated, spacing, config)
% Compute global 3D gamma index between two dose arrays on the same grid.
% Wraps CalcGamma (Mark Geurts).  Applies low-dose cutoff mask.
%
% Returns struct with: pass_rate, mean_gamma, max_gamma, gamma_map,
%                      num_evaluated, num_passed, computation_time_sec

    if ~isequal(size(reference), size(evaluated))
        error('test_RS_doses:SizeMismatch', ...
            'Reference [%s] and evaluated [%s] sizes differ.', ...
            mat2str(size(reference)), mat2str(size(evaluated)));
    end

    dose_pct   = config.gamma_dose_pct;
    dist_mm    = config.gamma_dist_mm;
    cutoff_pct = config.gamma_dose_cutoff_pct;

    max_ref       = max(reference(:));
    threshold     = max_ref * cutoff_pct / 100;
    mask          = (reference >= threshold) | (evaluated >= threshold);

    fprintf('    Dose threshold: %.4f Gy (%.0f%% of %.4f Gy)\n', ...
        threshold, cutoff_pct, max_ref);
    fprintf('    Voxels above threshold: %d / %d (%.1f%%)\n', ...
        sum(mask(:)), numel(mask), 100 * sum(mask(:)) / numel(mask));

    ref_s.start = [0, 0, 0];
    ref_s.width = spacing(:)';
    ref_s.data  = reference;

    tar_s.start = [0, 0, 0];
    tar_s.width = spacing(:)';
    tar_s.data  = evaluated;

    fprintf('    Running CalcGamma (%.0f%%/%.0fmm)...\n', dose_pct, dist_mm);
    t0 = tic;

    gamma_raw = CalcGamma(ref_s, tar_s, dose_pct, dist_mm, ...
        'local', 0, 'restrict', 1, 'res', 20, 'limit', 2);

    elapsed = toc(t0);
    fprintf('    CalcGamma done in %.1f sec\n', elapsed);

    gamma_map = gamma_raw;
    gamma_map(~mask) = NaN;

    valid = gamma_map(mask);

    gamma.pass_rate            = 100 * sum(valid <= 1) / numel(valid);
    gamma.mean_gamma           = mean(valid);
    gamma.max_gamma            = max(valid);
    gamma.gamma_map            = gamma_map;
    gamma.num_evaluated        = numel(valid);
    gamma.num_passed           = sum(valid <= 1);
    gamma.computation_time_sec = elapsed;
end


function rs_fig1_dose_comparison(ethos, ct1, ct3, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir)
% Figure 1: orthogonal dose views for ETHOS, CT_1 sum, CT_3 sum (3 rows x 3 cols).

    all_vals = ethos(:);
    if do_ct1, all_vals = [all_vals; ct1(:)]; end %#ok<AGROW>
    if do_ct3, all_vals = [all_vals; ct3(:)]; end %#ok<AGROW>
    dose_clim = [0, max(all_vals)];

    row_labels = {'ETHOS (truth)', 'CT\_1 sum', 'CT\_3 sum'};
    volumes    = {ethos, ct1, ct3};
    do_row     = [true, do_ct1, do_ct3];

    col_titles = {'Axial', 'Coronal', 'Sagittal'};

    fig = figure('Position', [50 50 1100 800], 'Visible', 'off');

    for row = 1:3
        if ~do_row(row), continue; end
        vol = volumes{row};

        slices = { vol(:, :, s0), ...           % axial
                   squeeze(vol(:, c0, :)), ...   % coronal
                   squeeze(vol(r0, :, :))' };    % sagittal (transpose for correct orientation)

        bm_slices = {};
        if ~isempty(body_mask)
            bm_slices = { body_mask(:, :, s0), ...
                          squeeze(body_mask(:, c0, :)), ...
                          squeeze(body_mask(r0, :, :))' };
        end

        for col = 1:3
            subplot(3, 3, (row-1)*3 + col);
            imagesc(slices{col}, dose_clim);
            colormap(gca, 'jet'); colorbar;
            axis image; axis off;
            if ~isempty(bm_slices)
                hold on;
                contour(double(bm_slices{col}), [0.5 0.5], 'g', 'LineWidth', 1);
                hold off;
            end
            if row == 1
                title(col_titles{col}, 'FontWeight', 'bold');
            end
            if col == 1
                ylabel(row_labels{row});
            end
        end
    end

    sgtitle(sprintf('Dose Comparison — Patient %s / %s', ...
        strrep(strtrim(inputname(1)), '_', '\_'), ''), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'dose_comparison.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
    set(fig, 'Visible', 'on');
end


function rs_fig2_dose_difference(ethos, ct1, ct3, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir)
% Figure 2: dose difference maps CT_1−ETHOS and CT_3−ETHOS (up to 2 rows x 3 cols).

    diffs    = {};
    titles   = {};
    if do_ct1, diffs{end+1} = ct1 - ethos; titles{end+1} = 'CT\_1 \minus ETHOS'; end
    if do_ct3, diffs{end+1} = ct3 - ethos; titles{end+1} = 'CT\_3 \minus ETHOS'; end

    if isempty(diffs), return; end

    all_diff = cellfun(@(d) max(abs(d(:))), diffs);
    diff_lim = max(all_diff);
    if diff_lim < 1e-10, diff_lim = 1; end

    col_titles = {'Axial', 'Coronal', 'Sagittal'};
    n_rows = numel(diffs);

    fig = figure('Position', [50 50 1100, 350*n_rows + 100], 'Visible', 'off');

    for row = 1:n_rows
        d = diffs{row};
        sl = { d(:, :, s0), squeeze(d(:, c0, :)), squeeze(d(r0, :, :))' };

        bm_slices = {};
        if ~isempty(body_mask)
            bm_slices = { body_mask(:, :, s0), ...
                          squeeze(body_mask(:, c0, :)), ...
                          squeeze(body_mask(r0, :, :))' };
        end

        for col = 1:3
            subplot(n_rows, 3, (row-1)*3 + col);
            imagesc(sl{col}, [-diff_lim, diff_lim]);
            colormap(gca, rdbu_colormap()); colorbar;
            axis image; axis off;
            if ~isempty(bm_slices)
                hold on;
                contour(double(bm_slices{col}), [0.5 0.5], 'k', 'LineWidth', 0.8);
                hold off;
            end
            if row == 1
                title(col_titles{col}, 'FontWeight', 'bold');
            end
            if col == 1
                ylabel(titles{row});
            end
        end
    end

    sgtitle('Dose Difference Maps (Gy)', 'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'dose_difference.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
    set(fig, 'Visible', 'on');
end


function rs_fig3_gamma_maps(gamma_ct1, gamma_ct3, do_ct1, do_ct3, ...
        r0, c0, s0, config, do_save, fig_dir)
% Figure 3: gamma maps for CT_1 and CT_3 vs ETHOS (up to 2 rows x 3 cols).

    gammas = {};
    labels = {};
    if do_ct1, gammas{end+1} = gamma_ct1; labels{end+1} = sprintf('CT\\_1 (%.1f%% pass)', gamma_ct1.pass_rate); end
    if do_ct3, gammas{end+1} = gamma_ct3; labels{end+1} = sprintf('CT\\_3 (%.1f%% pass)', gamma_ct3.pass_rate); end

    if isempty(gammas), return; end

    gamma_clim  = [0, 2];
    col_titles  = {'Axial', 'Coronal', 'Sagittal'};
    n_rows      = numel(gammas);
    bg_gray     = [0.85 0.85 0.85];

    fig = figure('Position', [50 50 1100, 350*n_rows + 100], 'Visible', 'off');

    for row = 1:n_rows
        gmap = gammas{row}.gamma_map;
        sl = { gmap(:, :, s0), squeeze(gmap(:, c0, :)), squeeze(gmap(r0, :, :))' };

        for col = 1:3
            ax = subplot(n_rows, 3, (row-1)*3 + col);
            set(ax, 'Color', bg_gray);
            imagesc(sl{col}, gamma_clim);
            colormap(ax, gamma_colormap()); colorbar;
            axis image; axis off;
            if row == 1
                title(col_titles{col}, 'FontWeight', 'bold');
            end
            if col == 1
                ylabel(labels{row});
            end
        end
    end

    pass_str = '';
    if do_ct1
        pass_str = sprintf('CT\\_1: %.1f%%', gamma_ct1.pass_rate);
    end
    if do_ct3
        if ~isempty(pass_str), pass_str = [pass_str, '  |  ']; end %#ok<AGROW>
        pass_str = [pass_str, sprintf('CT\\_3: %.1f%%', gamma_ct3.pass_rate)]; %#ok<AGROW>
    end

    sgtitle(sprintf('Gamma Maps (%.0f%%/%.0fmm) — %s', ...
        config.gamma_dose_pct, config.gamma_dist_mm, pass_str), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'gamma_maps.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
    set(fig, 'Visible', 'on');
end


function rs_fig4_ct1_vs_ct3_gamma(gamma_result, ct1_dose, r0, c0, s0, config, do_save, fig_dir)
% Figure 4: gamma map between CT_1 and CT_3 dose sums (1 row x 3 cols).

    gmap = gamma_result.gamma_map;
    sl   = { gmap(:, :, s0), squeeze(gmap(:, c0, :)), squeeze(gmap(r0, :, :))' };

    gamma_clim = [0, 2];
    col_titles = {'Axial', 'Coronal', 'Sagittal'};
    bg_gray    = [0.85 0.85 0.85];

    fig = figure('Position', [50 50 1100 400], 'Visible', 'off');

    for col = 1:3
        ax = subplot(1, 3, col);
        set(ax, 'Color', bg_gray);
        imagesc(sl{col}, gamma_clim);
        colormap(ax, gamma_colormap()); colorbar;
        axis image; axis off;
        title(col_titles{col}, 'FontWeight', 'bold');
        if col == 1
            ylabel(sprintf('CT\\_1 vs CT\\_3 (%.1f%% pass)', gamma_result.pass_rate));
        end
    end

    sgtitle(sprintf('CT\\_1 vs CT\\_3 Gamma (%.0f%%/%.0fmm) — pass rate %.1f%%', ...
        config.gamma_dose_pct, config.gamma_dist_mm, gamma_result.pass_rate), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'gamma_ct1_vs_ct3.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
    set(fig, 'Visible', 'on');
end


function cmap = rdbu_colormap()
% Diverging red-white-blue colormap for dose differences.
    n = 256; half = n / 2;
    r = [linspace(0.1, 1, half)'; linspace(1, 0.8, half)'];
    g = [linspace(0.2, 1, half)'; linspace(1, 0.1, half)'];
    b = [linspace(0.7, 1, half)'; linspace(1, 0.1, half)'];
    cmap = [r, g, b];
end


function cmap = gamma_colormap()
% Green -> yellow -> red colormap for gamma index display.
    n = 256; half = round(n / 2); rem_ = n - half;
    r = [linspace(0, 1, half)';    ones(rem_, 1)];
    g = [ones(half, 1) * 0.8;     linspace(0.8, 0, rem_)'];
    b = [linspace(0.2, 0, half)';  zeros(rem_, 1)];
    cmap = [r, g, b];
end
