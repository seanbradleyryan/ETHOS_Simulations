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
% Original (legacy) scheme: search for any RTDOSE*.dcm in sct directory
% rd_files = dir(fullfile(sct_dir, 'RTDOSE*.dcm'));
%
% if isempty(rd_files)
%     error('test_RS_doses:FileNotFound', ...
%         'No RTDOSE*.dcm file found in:\n  %s\nCheck CONFIG.treatment_site and that DICOM exports are present.', ...
%         sct_dir);
% end
%
% rd_path   = fullfile(sct_dir, rd_files(1).name);

% Current scheme: use canonical RD_reference.dcm as the ETHOS truth sample
rd_files = dir(fullfile(sct_dir, 'RD_reference.dcm'));

if isempty(rd_files)
    error('test_RS_doses:FileNotFound', ...
        'RD_reference.dcm not found in:\n  %s\nCheck CONFIG.treatment_site and that the reference dose export is present.', ...
        sct_dir);
end

rd_path   = fullfile(sct_dir, rd_files(1).name);
dose_info = dicominfo(rd_path);
ethos_dose = double(squeeze(dicomread(rd_path)));

if isfield(dose_info, 'DoseGridScaling')
    ethos_dose = ethos_dose * dose_info.DoseGridScaling;
end

fprintf('  File:           %s\n', rd_files(1).name);
fprintf('  Size:           [%d x %d x %d]\n', size(ethos_dose));
fprintf('  Class:          %s\n', class(ethos_dose));
if isfield(dose_info, 'DoseGridScaling')
    fprintf('  Dose scaling:   %.6g (applied)\n', dose_info.DoseGridScaling);
end
if isfield(dose_info, 'PixelSpacing') && isfield(dose_info, 'SliceThickness')
    fprintf('  DICOM spacing:  [%.2f x %.2f x %.2f] mm\n', ...
        dose_info.PixelSpacing(1), dose_info.PixelSpacing(2), dose_info.SliceThickness);
end
fprintf('  Empty?          %d   NaN count: %d   Inf count: %d\n', ...
    isempty(ethos_dose), sum(isnan(ethos_dose(:))), sum(isinf(ethos_dose(:))));
fprintf('  Stats (Gy):     max=%.4f   mean=%.4f   median=%.4f   std=%.4f\n', ...
    max(ethos_dose(:)), mean(ethos_dose(:)), median(ethos_dose(:)), std(ethos_dose(:)));
fprintf('  Nonzero voxels: %d / %d (%.2f%%)   sum=%.2f Gy-vox\n\n', ...
    nnz(ethos_dose), numel(ethos_dose), ...
    100 * nnz(ethos_dose) / numel(ethos_dose), sum(ethos_dose(:)));

%% ========================= STEP 2: LOAD PROCESSED FIELD DOSES ============

fprintf('[STEP 2] Loading processed RayStation field doses...\n');

[field_doses, sct_resampled, ~, dose_metadata] = ...
    load_processed_data(patient_id, session, CONFIG);

% --- Explicit geometry extraction (do not assume any upstream state) -----
if ~isstruct(dose_metadata)
    error('test_RS_doses:BadMetadata', 'dose_metadata is not a struct.');
end
if ~isfield(dose_metadata, 'dimensions') || numel(dose_metadata.dimensions) ~= 3
    error('test_RS_doses:BadMetadata', ...
        'dose_metadata.dimensions must be a 3-element vector.');
end
if ~isfield(dose_metadata, 'spacing') || numel(dose_metadata.spacing) ~= 3
    error('test_RS_doses:BadMetadata', ...
        'dose_metadata.spacing must be a 3-element vector.');
end

rs_size = double(dose_metadata.dimensions(:)');   % row vector [ny, nx, nz]
spacing = double(dose_metadata.spacing(:)');      % row vector [dx, dy, dz] mm

if any(rs_size <= 0) || any(spacing <= 0)
    error('test_RS_doses:BadMetadata', ...
        'Grid dimensions or spacing contain non-positive values.');
end

fprintf('  Field doses loaded:  %d entries\n', numel(field_doses));
n_empty = 0; n_with_dose = 0;
for ii = 1:numel(field_doses)
    if isempty(field_doses{ii})
        n_empty = n_empty + 1;
    elseif isfield(field_doses{ii}, 'dose_Gy') && ~isempty(field_doses{ii}.dose_Gy)
        n_with_dose = n_with_dose + 1;
    end
end
fprintf('  Non-empty entries:   %d   (empty: %d, with dose_Gy: %d)\n', ...
    numel(field_doses) - n_empty, n_empty, n_with_dose);
fprintf('  RS dose grid:        [%d x %d x %d]\n', rs_size);
fprintf('  RS dose spacing:     [%.3f x %.3f x %.3f] mm\n', spacing);
fprintf('  RS voxel volume:     %.4f mm^3\n', prod(spacing));
fprintf('  RS FOV:              [%.1f x %.1f x %.1f] mm\n\n', ...
    rs_size .* spacing);

%% ========================= STEP 3: ACCUMULATE CT_1 AND CT_3 ==============

fprintf('[STEP 3] Accumulating field doses by CT label...\n');

ct1_dose = zeros(rs_size, 'double');
ct3_dose = zeros(rs_size, 'double');
n_ct1 = 0;  n_ct3 = 0;  n_other = 0;
warned_other = false;

fprintf('  %-5s  %-8s  %-12s  %-10s  %-10s  %-10s\n', ...
    'idx', 'label', 'source_file', 'max(Gy)', 'mean(Gy)', 'nnz');
fprintf('  %s\n', repmat('-', 1, 70));

for i = 1:numel(field_doses)
    fd = field_doses{i};
    if isempty(fd), continue; end

    % --- Resolve ct_label explicitly within this script -------------------
    % Do NOT rely on upstream step15 having stamped fd.ct_label; this script
    % must be runnable standalone against the processed/ directory.
    % Resolution order:
    %   1. fd.ct_label field (if present and non-empty)
    %   2. Parse from fd.source_file filename: ..._<plan_type>_<CT_n>_B<beam>_<seg>.mat
    %   3. Empty -> treated as 'other'
    raw_label = '';
    if isfield(fd, 'ct_label') && ~isempty(fd.ct_label)
        raw_label = fd.ct_label;
    elseif isfield(fd, 'source_file') && ~isempty(fd.source_file)
        ct_tokens = regexp(fd.source_file, '_(CT[_]?\d+)_B\d+', ...
                           'tokens', 'once', 'ignorecase');
        if ~isempty(ct_tokens)
            raw_label = ct_tokens{1};
        end
    end
    label = upper(strtrim(char(raw_label)));

    % --- Validate dose volume size before accumulation --------------------
    if ~isfield(fd, 'dose_Gy') || isempty(fd.dose_Gy)
        warning('test_RS_doses:MissingDose', ...
            'Field %d has no dose_Gy; skipping.', i);
        continue;
    end
    if ~isequal(size(fd.dose_Gy), rs_size)
        warning('test_RS_doses:DoseSizeMismatch', ...
            'Field %d dose size [%s] differs from RS grid [%s]; skipping.', ...
            i, mat2str(size(fd.dose_Gy)), mat2str(rs_size));
        continue;
    end

    src_short = '';
    if isfield(fd, 'source_file') && ~isempty(fd.source_file)
        [~, src_short, ~] = fileparts(fd.source_file);
        if numel(src_short) > 12, src_short = src_short(1:12); end
    end

    fprintf('  %-5d  %-8s  %-12s  %-10.4f  %-10.4g  %-10d\n', ...
        i, label, src_short, max(fd.dose_Gy(:)), mean(fd.dose_Gy(:)), nnz(fd.dose_Gy));

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
                    i, label);
                warned_other = true;
            end
    end
end

do_ct1 = (n_ct1 > 0);
do_ct3 = (n_ct3 > 0);

if ~do_ct1
    warning('test_RS_doses:NoCT1', 'No CT_1 field doses found — skipping CT_1 analysis.');
end
if ~do_ct3
    warning('test_RS_doses:NoCT3', 'No CT_3 field doses found — skipping CT_3 analysis.');
end

fprintf('\n  --- Accumulated sums ---\n');
fprintf('  CT_1: %d fields  max=%.4f  mean=%.4f  nnz=%d (%.2f%%)  sum=%.2f\n', ...
    n_ct1, max(ct1_dose(:)), mean(ct1_dose(:)), nnz(ct1_dose), ...
    100*nnz(ct1_dose)/numel(ct1_dose), sum(ct1_dose(:)));
fprintf('  CT_3: %d fields  max=%.4f  mean=%.4f  nnz=%d (%.2f%%)  sum=%.2f\n', ...
    n_ct3, max(ct3_dose(:)), mean(ct3_dose(:)), nnz(ct3_dose), ...
    100*nnz(ct3_dose)/numel(ct3_dose), sum(ct3_dose(:)));
if n_other > 0
    fprintf('  Other/unlabelled fields excluded: %d\n', n_other);
end

% Sanity: NaN/Inf in accumulated sums
if any(isnan(ct1_dose(:))) || any(isinf(ct1_dose(:)))
    warning('test_RS_doses:CT1BadValues', 'CT_1 sum contains NaN/Inf values.');
end
if any(isnan(ct3_dose(:))) || any(isinf(ct3_dose(:)))
    warning('test_RS_doses:CT3BadValues', 'CT_3 sum contains NaN/Inf values.');
end
fprintf('\n');

%% ========================= STEP 4: RESAMPLE ETHOS IF NEEDED ==============

fprintf('[STEP 4] Checking ETHOS/RS grid alignment...\n');

if ~isequal(size(ethos_dose), rs_size)
    fprintf('  [WARNING] ETHOS size [%d x %d x %d] differs from RS size [%d x %d x %d]\n', ...
        size(ethos_dose), rs_size);
    pre_max  = max(ethos_dose(:));
    pre_mean = mean(ethos_dose(:));
    pre_sum  = sum(ethos_dose(:));
    fprintf('  Pre-resample:  max=%.4f  mean=%.4f  sum=%.2f\n', pre_max, pre_mean, pre_sum);
    fprintf('  Resampling ETHOS to RS grid (trilinear)...\n');
    ethos_dose = rs_resample(ethos_dose, rs_size);
    fprintf('  Post-resample: max=%.4f  mean=%.4f  sum=%.2f  size=[%d x %d x %d]\n', ...
        max(ethos_dose(:)), mean(ethos_dose(:)), sum(ethos_dose(:)), size(ethos_dose));
    fprintf('  Max changed by %.2f%%, mean by %.2f%%\n\n', ...
        100*(max(ethos_dose(:))-pre_max)/pre_max, ...
        100*(mean(ethos_dose(:))-pre_mean)/pre_mean);
else
    fprintf('  Grids match. No resampling needed.\n\n');
end

%% ========================= STEP 5: GAMMA ANALYSIS ========================

fprintf('[STEP 5] Running gamma analysis (%.0f%%/%.0fmm)...\n', ...
    CONFIG.gamma_dose_pct, CONFIG.gamma_dist_mm);

gamma_ct1 = [];
gamma_ct3 = [];

if do_ct1
    fprintf('\n  --- ETHOS vs CT_1 ---\n');
    gamma_ct1 = rs_gamma(ethos_dose, ct1_dose, spacing, CONFIG);
end

if do_ct3
    fprintf('\n  --- ETHOS vs CT_3 ---\n');
    gamma_ct3 = rs_gamma(ethos_dose, ct3_dose, spacing, CONFIG);
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

    fprintf('  Slice indices at ETHOS max-dose voxel: row=%d, col=%d, slice=%d\n', r0, c0, s0);
    fprintf('  ETHOS max dose at this voxel: %.4f Gy\n', ethos_dose(r0, c0, s0));

    % Figure 1 — Dose comparison (absolute, shared color limits)
    rs_fig1_dose_comparison(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir, patient_id, session, spacing);

    % Figure 1b — Globally normalized (each shown as fraction of GLOBAL max across all 3 volumes)
    rs_fig1b_global_norm(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir, patient_id, session);

    % Figure 1c — Locally normalized (each shown as fraction of ITS OWN max)
    rs_fig1c_local_norm(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir, patient_id, session);

    % Figure 2 — Dose difference maps
    rs_fig2_dose_difference(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir);

    % Figure 3 — Gamma maps
    rs_fig3_gamma_maps(gamma_ct1, gamma_ct3, do_ct1, do_ct3, ...
        r0, c0, s0, CONFIG, do_save, fig_dir);

    % Figure 4 — Dose-volume histograms (cumulative DVH)
    rs_fig4_dvh(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        do_save, fig_dir);

    % Figure 5 — Line profiles through max-dose voxel
    rs_fig5_profiles(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        r0, c0, s0, spacing, do_save, fig_dir);

    % Figure 6 — Voxelwise scatter (ETHOS vs CT_1/CT_3) above low-dose cutoff
    rs_fig6_scatter(ethos_dose, ct1_dose, ct3_dose, do_ct1, do_ct3, ...
        CONFIG, do_save, fig_dir);

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
        r0, c0, s0, body_mask, do_save, fig_dir, patient_id, session, spacing)
% Figure 1: orthogonal dose views for ETHOS, CT_1 sum, CT_3 sum (3 rows x 3 cols).
% Absolute dose, single shared color scale (global max across the 3 volumes).

    all_vals = ethos(:);
    if do_ct1, all_vals = [all_vals; ct1(:)]; end %#ok<AGROW>
    if do_ct3, all_vals = [all_vals; ct3(:)]; end %#ok<AGROW>
    dose_clim = [0, max(all_vals)];

    row_labels = { ...
        sprintf('ETHOS (truth)\nmax=%.3f Gy', max(ethos(:))), ...
        sprintf('CT\\_1 sum\nmax=%.3f Gy', max(ct1(:))), ...
        sprintf('CT\\_3 sum\nmax=%.3f Gy', max(ct3(:))) };
    volumes = {ethos, ct1, ct3};
    do_row  = [true, do_ct1, do_ct3];

    col_titles = { ...
        sprintf('Axial (slice %d, z = %.1f mm)', s0, (s0-1)*spacing(3)), ...
        sprintf('Coronal (col %d, x = %.1f mm)', c0, (c0-1)*spacing(1)), ...
        sprintf('Sagittal (row %d, y = %.1f mm)', r0, (r0-1)*spacing(2)) };

    fig = figure('Name', 'Fig 1: Dose Comparison (absolute)', ...
                 'Position', [50 50 1200 850], 'Visible', 'on');

    for row = 1:3
        if ~do_row(row), continue; end
        vol = volumes{row};

        slices = { vol(:, :, s0), ...
                   squeeze(vol(:, c0, :)), ...
                   squeeze(vol(r0, :, :))' };

        bm_slices = {};
        if ~isempty(body_mask)
            bm_slices = { body_mask(:, :, s0), ...
                          squeeze(body_mask(:, c0, :)), ...
                          squeeze(body_mask(r0, :, :))' };
        end

        for col = 1:3
            subplot(3, 3, (row-1)*3 + col);
            imagesc(slices{col}, dose_clim);
            colormap(gca, 'jet');
            cb = colorbar; cb.Label.String = 'Dose (Gy)';
            axis image; axis off;
            if ~isempty(bm_slices)
                hold on;
                contour(double(bm_slices{col}), [0.5 0.5], 'g', 'LineWidth', 1);
                hold off;
            end
            if row == 1
                title(col_titles{col}, 'FontWeight', 'bold', 'FontSize', 10);
            end
            if col == 1
                ylabel(row_labels{row}, 'FontWeight', 'bold', 'Visible', 'on');
            end
        end
    end

    sgtitle(sprintf('Fig 1 — Absolute Dose (shared colour scale, global max = %.3f Gy)\nPatient %s / %s', ...
        dose_clim(2), strrep(patient_id,'_','\_'), strrep(session,'_','\_')), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'dose_comparison.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
end


function rs_fig1b_global_norm(ethos, ct1, ct3, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir, patient_id, session)
% Figure 1b: orthogonal dose views, each volume normalized to the GLOBAL max
% (max across ETHOS, CT_1 sum, CT_3 sum). Lets you compare absolute amplitude.

    global_max = max(ethos(:));
    if do_ct1, global_max = max(global_max, max(ct1(:))); end
    if do_ct3, global_max = max(global_max, max(ct3(:))); end
    if global_max <= 0, global_max = 1; end

    eN  = ethos / global_max;
    c1N = ct1   / global_max;
    c3N = ct3   / global_max;

    row_labels = { ...
        sprintf('ETHOS (truth)\nlocal max / global = %.2f', max(ethos(:))/global_max), ...
        sprintf('CT\\_1 sum\nlocal max / global = %.2f', max(ct1(:))/global_max), ...
        sprintf('CT\\_3 sum\nlocal max / global = %.2f', max(ct3(:))/global_max) };
    vols    = {eN, c1N, c3N};
    do_row  = [true, do_ct1, do_ct3];
    col_titles = {'Axial', 'Coronal', 'Sagittal'};

    fig = figure('Name', 'Fig 1b: Globally normalized', ...
                 'Position', [80 80 1200 850], 'Visible', 'on');

    for row = 1:3
        if ~do_row(row), continue; end
        v = vols{row};
        slices = { v(:, :, s0), squeeze(v(:, c0, :)), squeeze(v(r0, :, :))' };

        bm_slices = {};
        if ~isempty(body_mask)
            bm_slices = { body_mask(:, :, s0), ...
                          squeeze(body_mask(:, c0, :)), ...
                          squeeze(body_mask(r0, :, :))' };
        end

        for col = 1:3
            subplot(3, 3, (row-1)*3 + col);
            imagesc(slices{col}, [0, 1]);
            colormap(gca, 'jet');
            cb = colorbar; cb.Label.String = 'Dose / global max';
            axis image; axis off;
            if ~isempty(bm_slices)
                hold on;
                contour(double(bm_slices{col}), [0.5 0.5], 'g', 'LineWidth', 1);
                hold off;
            end
            if row == 1, title(col_titles{col}, 'FontWeight', 'bold'); end
            if col == 1, ylabel(row_labels{row}, 'FontWeight', 'bold', 'Visible', 'on'); end
        end
    end

    sgtitle(sprintf('Fig 1b — Globally Normalized (global max = %.3f Gy)\nPatient %s / %s', ...
        global_max, strrep(patient_id,'_','\_'), strrep(session,'_','\_')), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'dose_comparison_globalnorm.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
end


function rs_fig1c_local_norm(ethos, ct1, ct3, do_ct1, do_ct3, ...
        r0, c0, s0, body_mask, do_save, fig_dir, patient_id, session)
% Figure 1c: orthogonal dose views, each volume normalized to its OWN max.
% Lets you compare spatial distribution / shape regardless of absolute level.

    local_e  = max(ethos(:)); if local_e  <= 0, local_e  = 1; end
    local_c1 = max(ct1(:));   if local_c1 <= 0, local_c1 = 1; end
    local_c3 = max(ct3(:));   if local_c3 <= 0, local_c3 = 1; end

    eN  = ethos / local_e;
    c1N = ct1   / local_c1;
    c3N = ct3   / local_c3;

    row_labels = { ...
        sprintf('ETHOS\nself-max = %.3f Gy', local_e), ...
        sprintf('CT\\_1 sum\nself-max = %.3f Gy', local_c1), ...
        sprintf('CT\\_3 sum\nself-max = %.3f Gy', local_c3) };
    vols    = {eN, c1N, c3N};
    do_row  = [true, do_ct1, do_ct3];
    col_titles = {'Axial', 'Coronal', 'Sagittal'};

    fig = figure('Name', 'Fig 1c: Locally normalized', ...
                 'Position', [110 110 1200 850], 'Visible', 'on');

    for row = 1:3
        if ~do_row(row), continue; end
        v = vols{row};
        slices = { v(:, :, s0), squeeze(v(:, c0, :)), squeeze(v(r0, :, :))' };

        bm_slices = {};
        if ~isempty(body_mask)
            bm_slices = { body_mask(:, :, s0), ...
                          squeeze(body_mask(:, c0, :)), ...
                          squeeze(body_mask(r0, :, :))' };
        end

        for col = 1:3
            subplot(3, 3, (row-1)*3 + col);
            imagesc(slices{col}, [0, 1]);
            colormap(gca, 'jet');
            cb = colorbar; cb.Label.String = 'Dose / self max';
            axis image; axis off;
            if ~isempty(bm_slices)
                hold on;
                contour(double(bm_slices{col}), [0.5 0.5], 'g', 'LineWidth', 1);
                hold off;
            end
            if row == 1, title(col_titles{col}, 'FontWeight', 'bold'); end
            if col == 1, ylabel(row_labels{row}, 'FontWeight', 'bold', 'Visible', 'on'); end
        end
    end

    sgtitle(sprintf('Fig 1c — Locally Normalized (shape comparison)\nPatient %s / %s', ...
        strrep(patient_id,'_','\_'), strrep(session,'_','\_')), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'dose_comparison_localnorm.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
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

    % Annotate titles with per-row min/max/mean for context
    for kk = 1:numel(diffs)
        d = diffs{kk};
        titles{kk} = sprintf('%s\nmin=%.3f max=%.3f mean=%.3f Gy', ...
            titles{kk}, min(d(:)), max(d(:)), mean(d(:))); %#ok<AGROW>
    end

    col_titles = {'Axial', 'Coronal', 'Sagittal'};
    n_rows = numel(diffs);

    fig = figure('Name', 'Fig 2: Dose Difference', ...
                 'Position', [140 50 1100, 350*n_rows + 120], 'Visible', 'on');

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
            colormap(gca, rdbu_colormap());
            cb = colorbar; cb.Label.String = 'Dose diff (Gy)';
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
                ylabel(titles{row}, 'FontWeight', 'bold', 'Visible', 'on');
            end
        end
    end

    sgtitle(sprintf('Fig 2 — Dose Difference Maps (Gy), symmetric scale ±%.3f Gy', diff_lim), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'dose_difference.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
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

    fig = figure('Name', 'Fig 3: Gamma Maps', ...
                 'Position', [170 50 1100, 350*n_rows + 120], 'Visible', 'on');

    for row = 1:n_rows
        gmap = gammas{row}.gamma_map;
        sl = { gmap(:, :, s0), squeeze(gmap(:, c0, :)), squeeze(gmap(r0, :, :))' };

        for col = 1:3
            ax = subplot(n_rows, 3, (row-1)*3 + col);
            set(ax, 'Color', bg_gray);
            imagesc(sl{col}, gamma_clim);
            colormap(ax, gamma_colormap());
            cb = colorbar; cb.Label.String = '\gamma index';
            axis image; axis off;
            if row == 1
                title(col_titles{col}, 'FontWeight', 'bold');
            end
            if col == 1
                ylabel(labels{row}, 'FontWeight', 'bold', 'Visible', 'on');
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

    sgtitle(sprintf('Fig 3 — Gamma Maps (%.0f%%/%.0fmm, global) — %s', ...
        config.gamma_dose_pct, config.gamma_dist_mm, pass_str), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'gamma_maps.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
end


function rs_fig4_dvh(ethos, ct1, ct3, do_ct1, do_ct3, do_save, fig_dir)
% Figure 4: Cumulative dose-volume histograms (whole-grid).

    fig = figure('Name', 'Fig 4: DVH', ...
                 'Position', [200 100 900 600], 'Visible', 'on');
    hold on;

    legends = {};
    [d, v] = local_cdvh(ethos);
    plot(d, v, 'k-', 'LineWidth', 2);
    legends{end+1} = sprintf('ETHOS (max=%.3f Gy)', max(ethos(:)));

    if do_ct1
        [d, v] = local_cdvh(ct1);
        plot(d, v, 'b-', 'LineWidth', 1.5);
        legends{end+1} = sprintf('CT\\_1 sum (max=%.3f Gy)', max(ct1(:)));
    end
    if do_ct3
        [d, v] = local_cdvh(ct3);
        plot(d, v, 'r-', 'LineWidth', 1.5);
        legends{end+1} = sprintf('CT\\_3 sum (max=%.3f Gy)', max(ct3(:)));
    end

    grid on; box on;
    xlabel('Dose (Gy)', 'FontWeight', 'bold');
    ylabel('Fractional volume \geq Dose', 'FontWeight', 'bold');
    title('Fig 4 — Cumulative Dose-Volume Histogram (whole grid)', 'FontSize', 13, 'FontWeight', 'bold');
    legend(legends, 'Location', 'northeast', 'Interpreter', 'tex');
    ylim([0, 1]);
    hold off;

    if do_save
        out = fullfile(fig_dir, 'dvh.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
end


function [dose_axis, frac_vol] = local_cdvh(vol)
% Cumulative DVH helper. Returns sorted dose axis and fractional volume >= each dose.
    v = sort(vol(:), 'descend');
    n = numel(v);
    frac_vol = (1:n)' / n;
    dose_axis = v;
    % Downsample for plotting if huge
    if n > 5000
        idx = round(linspace(1, n, 5000));
        dose_axis = dose_axis(idx);
        frac_vol  = frac_vol(idx);
    end
end


function rs_fig5_profiles(ethos, ct1, ct3, do_ct1, do_ct3, ...
        r0, c0, s0, spacing, do_save, fig_dir)
% Figure 5: 1D line profiles through the ETHOS max-dose voxel along each axis.

    ny = size(ethos, 1); nx = size(ethos, 2); nz = size(ethos, 3);

    yy = ((1:ny) - 1) * spacing(2);
    xx = ((1:nx) - 1) * spacing(1);
    zz = ((1:nz) - 1) * spacing(3);

    e_y = squeeze(ethos(:, c0, s0));   % vary row (Y)
    e_x = squeeze(ethos(r0, :, s0));   % vary col (X)
    e_z = squeeze(ethos(r0, c0, :));   % vary slice (Z)

    fig = figure('Name', 'Fig 5: Line Profiles', ...
                 'Position', [230 130 1300 420], 'Visible', 'on');

    axes_data = { ...
        struct('ax', xx, 'e', e_x, 'lbl', 'X (mm)', 'ttl', sprintf('X profile (row %d, slice %d)', r0, s0), 'getC1', @() squeeze(ct1(r0, :, s0)), 'getC3', @() squeeze(ct3(r0, :, s0))), ...
        struct('ax', yy, 'e', e_y, 'lbl', 'Y (mm)', 'ttl', sprintf('Y profile (col %d, slice %d)', c0, s0), 'getC1', @() squeeze(ct1(:, c0, s0)), 'getC3', @() squeeze(ct3(:, c0, s0))), ...
        struct('ax', zz, 'e', e_z, 'lbl', 'Z (mm)', 'ttl', sprintf('Z profile (row %d, col %d)', r0, c0), 'getC1', @() squeeze(ct1(r0, c0, :)), 'getC3', @() squeeze(ct3(r0, c0, :))) };

    for k = 1:3
        subplot(1, 3, k); hold on;
        plot(axes_data{k}.ax, axes_data{k}.e, 'k-', 'LineWidth', 2);
        legs = {'ETHOS'};
        if do_ct1
            plot(axes_data{k}.ax, axes_data{k}.getC1(), 'b-', 'LineWidth', 1.5);
            legs{end+1} = 'CT\_1 sum'; %#ok<AGROW>
        end
        if do_ct3
            plot(axes_data{k}.ax, axes_data{k}.getC3(), 'r-', 'LineWidth', 1.5);
            legs{end+1} = 'CT\_3 sum'; %#ok<AGROW>
        end
        grid on; box on;
        xlabel(axes_data{k}.lbl, 'FontWeight', 'bold');
        ylabel('Dose (Gy)', 'FontWeight', 'bold');
        title(axes_data{k}.ttl);
        legend(legs, 'Location', 'best', 'Interpreter', 'tex');
        hold off;
    end

    sgtitle('Fig 5 — Line Profiles Through ETHOS Max-Dose Voxel', ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'profiles.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
end


function rs_fig6_scatter(ethos, ct1, ct3, do_ct1, do_ct3, config, do_save, fig_dir)
% Figure 6: Voxelwise scatter ETHOS vs each CT sum, restricted to voxels
% above the gamma low-dose cutoff. Annotated with Pearson r and slope.

    threshold = max(ethos(:)) * config.gamma_dose_cutoff_pct / 100;
    n_panels = double(do_ct1) + double(do_ct3);
    if n_panels == 0, return; end

    fig = figure('Name', 'Fig 6: Voxelwise Scatter', ...
                 'Position', [260 160 600*n_panels 550], 'Visible', 'on');

    panel = 0;
    if do_ct1
        panel = panel + 1;
        subplot(1, n_panels, panel);
        local_scatter_panel(ethos, ct1, threshold, 'CT\_1 sum', 'b');
    end
    if do_ct3
        panel = panel + 1;
        subplot(1, n_panels, panel);
        local_scatter_panel(ethos, ct3, threshold, 'CT\_3 sum', 'r');
    end

    sgtitle(sprintf('Fig 6 — Voxelwise Dose Correlation (voxels > %.0f%% of ETHOS max = %.3f Gy)', ...
        config.gamma_dose_cutoff_pct, threshold), ...
        'FontSize', 13, 'FontWeight', 'bold');

    if do_save
        out = fullfile(fig_dir, 'scatter.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('  Saved: %s\n', out);
    end
end


function local_scatter_panel(ref, eval, threshold, label, color)
    mask = (ref >= threshold) | (eval >= threshold);
    x = ref(mask); y = eval(mask);

    % Subsample if very large for plotting performance
    n_max_plot = 50000;
    if numel(x) > n_max_plot
        idx = randperm(numel(x), n_max_plot);
        xp = x(idx); yp = y(idx);
    else
        xp = x; yp = y;
    end

    scatter(xp, yp, 4, color, 'filled', 'MarkerFaceAlpha', 0.2);
    hold on;
    lim = max([max(x), max(y), 1e-3]);
    plot([0, lim], [0, lim], 'k--', 'LineWidth', 1.2);   % identity line

    % Linear fit (no intercept fudge — slope only via least squares)
    if numel(x) > 1
        p = polyfit(x, y, 1);
        plot([0, lim], polyval(p, [0, lim]), 'g-', 'LineWidth', 1.5);
        r = corr(x, y);
        text(0.05*lim, 0.92*lim, sprintf('n=%d\nr=%.4f\nslope=%.3f\nintercept=%.3f Gy', ...
            numel(x), r, p(1), p(2)), 'FontSize', 10, 'BackgroundColor', 'w', ...
            'EdgeColor', 'k', 'VerticalAlignment', 'top');
    end
    hold off;
    grid on; box on; axis equal; xlim([0 lim]); ylim([0 lim]);
    xlabel('ETHOS dose (Gy)', 'FontWeight', 'bold');
    ylabel(sprintf('%s dose (Gy)', label), 'FontWeight', 'bold', 'Interpreter', 'tex');
    title(sprintf('ETHOS vs %s', label), 'Interpreter', 'tex');
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
