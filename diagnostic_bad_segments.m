%% =========================================================================
%  DIAGNOSTIC_BAD_SEGMENTS.m
%  Segments whose two truth doses (RS CT_1 vs RS CT_3) disagree the most
%
%  Stored data only (nothing recomputed): pass rates from the Step 2.5
%  summary, gamma maps + LS recon gains from the segment_metrics folded into
%  each CT_1 recon file.
%    1. Histogram of the truth1_vs_truth3 pass rate over all segments (2.5% bins).
%    2. The N_WORST lowest-pass segments, printed.
%    3. One figure per worst segment, three tabs, orthogonal views at the
%       truth CT_1 max-dose voxel, CT_1 body contour in green:
%         Doses          - truth CT_1, truth CT_3, recon CT_1, recon CT_3
%         Gamma          - stored gamma map of each Step 2.5 comparison
%         Anatomy change - truth CT_3 - truth CT_1 and density CT_3 - CT_1,
%                          with failing truth1_vs_truth3 voxels (gamma > 1)
%                          outlined in black. This shows whether the dose
%                          disagreement follows an anatomy change between the
%                          two CBCTs (gas, body contour, organ motion).
%
%  Summary rows have no segment numbers, so they are taken from the CT_1
%  recon filenames. Step 2.5 orders a beam's rows the same way: sorted
%  segments with both CT_1 and CT_3 recons. Each picked segment's stored pass
%  rate is checked against its summary row.
%  Figures (.fig, tabs kept) go to
%  AnalysisResults/<pid>/<session>/<label>/diagnostic_bad_segments/.
%  =========================================================================

clear; clc; close all;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% ========================= CONFIGURATION ================================

CONFIG = get_default_config();     % patient_id, session, working_dir, sim settings
CONFIG.config_hash = '';           % '' = hash of this CONFIG; or e.g. 'a9a3e1e6'

N_WORST          = 5;              % segments to plot
BIN_WIDTH_PCT    = 2.5;            % histogram bin width (%)
DENSITY_CLIM     = 300;            % +/- color limit of the density difference (kg/m^3)

%% ========================= LOCATE STORED RESULTS ========================

pid   = CONFIG.patient_id;
ses   = CONFIG.session;
hash8 = CONFIG.config_hash;
if isempty(hash8)
    hash8 = compute_sim_config_hash(CONFIG);
end
sim_dir = get_simulation_dir(CONFIG.working_dir, pid, ses, hash8);

run_cfg = hashtoconfig(hash8, CONFIG);
if isfield(run_cfg, 'reconstruction_mode') && strcmp(run_cfg.reconstruction_mode, 'beam')
    error('diagnostic_bad_segments:BeamMode', ...
        'Hash %s is a beam-mode run; this script needs per-segment results.', hash8);
end
out_dir = fullfile(CONFIG.working_dir, 'AnalysisResults', pid, ses, ...
    config_label(run_cfg), 'diagnostic_bad_segments');
if ~isfolder(out_dir), mkdir(out_dir); end

summary_file = fullfile(sim_dir, sprintf('segment_metrics_summary_%s.mat', hash8));
if ~isfile(summary_file)
    error('diagnostic_bad_segments:NoSummary', ...
        'No Step 2.5 summary:\n  %s\nRun step25_segment_metrics (or step3_analysis) first.', summary_file);
end
SM        = load(summary_file);
plan_type = char(SM.plan_type);
c_tt      = find(strcmp(SM.comparisons, 'truth1_vs_truth3'));
fprintf('[DIAG] %s %s | hash %s | plan %s | gamma %g%%/%gmm\n', pid, ses, hash8, ...
    plan_type, SM.gamma_dose_pct, SM.gamma_dist_mm);

%% ========================= ONE ROW PER SEGMENT ==========================
% Segment numbers from the CT_1 recon filenames (..._CT_1_B<beam>_<seg>_recon_<hash>.mat)
% that also have a CT_3 recon.

files = dir(fullfile(sim_dir, sprintf('*_CT_1_B*_recon_%s.mat', hash8)));
names = {files.name};
if ~strcmpi(plan_type, 'any')
    names = names(contains(names, ['_', plan_type, '_']));
end
tok   = regexp(names, '_B(\d+)_(\d+)_recon_', 'tokens', 'once');
keep  = ~cellfun(@isempty, tok);
names = names(keep);
tok   = tok(keep);
file_beam = cellfun(@(t) str2double(t{1}), tok);
file_seg  = cellfun(@(t) str2double(t{2}), tok);
has_ct3   = cellfun(@(n) isfile(fullfile(sim_dir, strrep(n, '_CT_1_', '_CT_3_'))), names);

beam = []; seg = []; pass = [];
for n = 1:numel(SM.beams)
    b = SM.beams(n);
    p = SM.gamma.seg_pass{n}(:, c_tt);
    s = unique(file_seg(file_beam == b & has_ct3));
    if numel(s) ~= numel(p)
        warning('diagnostic_bad_segments:RowCount', ...
            'Beam %d: %d recon pairs on disk vs %d summary rows; its segments are left out of the worst list.', ...
            b, numel(s), numel(p));
        s = nan(size(p));
    end
    beam = [beam; repmat(b, numel(p), 1)]; %#ok<AGROW>
    seg  = [seg; s(:)];                    %#ok<AGROW>
    pass = [pass; p(:)];                   %#ok<AGROW>
end

valid    = find(~isnan(pass) & ~isnan(seg));
[~, ord] = sort(pass(valid));
worst    = valid(ord(1:min(N_WORST, numel(valid))));
if isempty(worst)
    error('diagnostic_bad_segments:NoSegments', 'No segment has both a truth pass rate and a segment number.');
end

fprintf('[DIAG] %d segments with a truth pass rate. Worst %d:\n', nnz(~isnan(pass)), numel(worst));
for r = worst'
    fprintf('    beam %3d  seg %4d  truth1 vs truth3 %6.2f %%\n', beam(r), seg(r), pass(r));
end

%% ========================= HISTOGRAM ====================================

p  = pass(~isnan(pass));
lo = min(floor(min(p) / BIN_WIDTH_PCT) * BIN_WIDTH_PCT, 100 - BIN_WIDTH_PCT);
fig = figure('Name', 'Truth CT_1 vs CT_3 pass rates', 'Color', 'w', 'NumberTitle', 'off');
histogram(p, 'BinEdges', lo:BIN_WIDTH_PCT:100);
xline(pass(worst(end)), 'r--', sprintf('%d worst', numel(worst)), 'LineWidth', 1.5);
grid on;
xlabel(sprintf('Truth CT_1 vs truth CT_3 gamma pass rate (%%, %g%%/%gmm)', ...
    SM.gamma_dose_pct, SM.gamma_dist_mm), 'Interpreter', 'none');
ylabel('Segments');
title(sprintf('%s %s | %d segments | median %.1f%%', pid, ses, numel(p), median(p)), ...
    'Interpreter', 'none');
savefig(fig, fullfile(out_dir, 'truth_pass_rate_histogram.fig'));

%% ========================= ONE FIGURE PER WORST SEGMENT =================

for r = worst'
    b = beam(r);
    s = seg(r);
    fprintf('[DIAG] Plotting beam %d seg %d...\n', b, s);
    S = load_recon_dose_data(pid, ses, CONFIG, 'Mode', 'set', 'Beam', b, 'Segment', s, ...
        'PlanType', plan_type, 'Hash', hash8, 'IncludeEthos', false);
    cts = arrayfun(@(f) strrep(char(f.rtplan.ct_label), '-', '_'), S.fields, 'UniformOutput', false);
    f1  = S.fields(find(strcmp(cts, 'CT_1'), 1));
    f3  = S.fields(find(strcmp(cts, 'CT_3'), 1));

    L  = load(f1.recon_file, 'segment_metrics');
    sm = L.segment_metrics;
    if abs(sm.comparison(c_tt).gamma_pass_rate - pass(r)) > 1e-6
        error('diagnostic_bad_segments:RowMismatch', ...
            ['Beam %d seg %d: stored pass rate %.2f%% but summary row %.2f%%. The summary ' ...
             'does not match the recon files; re-run step25_segment_metrics.'], ...
            b, s, sm.comparison(c_tt).gamma_pass_rate, pass(r));
    end

    t1 = double(f1.rs_dose);
    t3 = double(f3.rs_dose);
    r1 = double(f1.recon_dose) * sm.recon_ct1_gain;   % same scaling Step 2.5 scored
    r3 = double(f3.recon_dose) * sm.recon_ct3_gain;
    body = f1.cbct.bodyMask;   % Step 1.5 CBCTs are on the dose grid
    drho = double(f3.cbct.cubeDensity) - double(f1.cbct.cubeDensity);
    ddose_pct = 100 * (t3 - t1) / max(t1(:));

    nComp = numel(sm.comparison);
    gmaps = cell(1, nComp);
    for c = 1:nComp
        gmaps{c} = nan(sm.vol_size);
        if ~isempty(sm.comparison(c).gamma_vals)
            gmaps{c}(sm.comparison(c).mask_idx) = double(sm.comparison(c).gamma_vals);
        end
    end
    gamma_titles = arrayfun(@(c) sprintf('%s %.1f%%', c.name, c.gamma_pass_rate), ...
        sm.comparison, 'UniformOutput', false);
    fail_mask = gmaps{c_tt} > 1;   % NaN (outside the eval mask) counts as not failing

    [~, imax]    = max(t1(:));
    [iy, ix, iz] = ind2sub(size(t1), imax);
    ijk  = [iy, ix, iz];
    dmax = max([t1(:); t3(:)]);
    dlim = max(abs(ddose_pct(:)));
    if dlim == 0, dlim = 1; end

    fig = figure('Name', sprintf('Beam %d seg %d | truth1 vs truth3 %.1f%%', b, s, pass(r)), ...
        'Color', 'w', 'NumberTitle', 'off', 'Position', [40, 40, 1500, 950]);
    tg = uitabgroup(fig);
    draw_tab(uitab(tg, 'Title', 'Doses'), {t1, t3, r1, r3}, ...
        {'Truth CT_1 (Gy)', 'Truth CT_3 (Gy)', 'Recon CT_1 (Gy)', 'Recon CT_3 (Gy)'}, ...
        {jet(256), jet(256), jet(256), jet(256)}, repmat({[0, dmax]}, 1, 4), ...
        ijk, sm.spacing, body, []);
    draw_tab(uitab(tg, 'Title', 'Gamma'), gmaps, gamma_titles, ...
        repmat({blue_white_red()}, 1, nComp), repmat({[0, 2]}, 1, nComp), ...
        ijk, sm.spacing, body, []);
    draw_tab(uitab(tg, 'Title', 'Anatomy change'), {ddose_pct, drho}, ...
        {'Truth CT_3 - CT_1 (% of CT_1 max)', 'Density CT_3 - CT_1 (kg/m^3)'}, ...
        {blue_white_red(), blue_white_red()}, {[-dlim, dlim], [-DENSITY_CLIM, DENSITY_CLIM]}, ...
        ijk, sm.spacing, body, fail_mask);
    savefig(fig, fullfile(out_dir, sprintf('bad_segment_B%d_S%d.fig', b, s)));
    clear S f1 f3 L sm;
end
fprintf('[DIAG] Figures saved to %s\n', out_dir);


%% ========================= LOCAL HELPERS ================================

function draw_tab(tab, vols, labels, cmaps, clims, ijk, spacing, body, outline)
%DRAW_TAB 3 rows (transverse / coronal / sagittal) x one column per volume.
%  NaN voxels are blank. Green = body contour; black = outline mask (skipped if []).
    views = {'Transverse', 'Coronal', 'Sagittal'};
    tl = tiledlayout(tab, 3, numel(vols), 'TileSpacing', 'compact', 'Padding', 'compact');
    for v = 1:3
        for c = 1:numel(vols)
            ax = nexttile(tl);
            [img, xv, yv] = ortho_slice(vols{c}, v, ijk, spacing);
            h = imagesc(ax, xv, yv, img);
            set(h, 'AlphaData', ~isnan(img));
            axis(ax, 'image');
            if v > 1, set(ax, 'YDir', 'normal'); end   % superior (higher iz) at the top
            colormap(ax, cmaps{c}); clim(ax, clims{c}); colorbar(ax);
            hold(ax, 'on');
            m = ortho_slice(body, v, ijk, spacing);
            if any(m(:)), contour(ax, xv, yv, double(m), [0.5, 0.5], 'g', 'LineWidth', 1); end
            if ~isempty(outline)
                m = ortho_slice(outline, v, ijk, spacing);
                if any(m(:)), contour(ax, xv, yv, double(m), [0.5, 0.5], 'k', 'LineWidth', 1); end
            end
            hold(ax, 'off');
            xlabel(ax, 'mm'); ylabel(ax, 'mm');
            title(ax, sprintf('%s | %s', labels{c}, views{v}), 'Interpreter', 'none');
        end
    end
end

function [img, xv, yv] = ortho_slice(vol, v, ijk, spacing)
%ORTHO_SLICE Slice through voxel ijk = [iy ix iz]; v = 1 transverse, 2 coronal,
%  3 sagittal. Axes in mm from the grid corner; spacing = [dx dy dz].
    [ny, nx, nz] = size(vol);
    x = (0:nx - 1) * spacing(1);
    y = (0:ny - 1) * spacing(2);
    z = (0:nz - 1) * spacing(3);
    switch v
        case 1, img = vol(:, :, ijk(3));           xv = x; yv = y;
        case 2, img = squeeze(vol(ijk(1), :, :))'; xv = x; yv = z;
        case 3, img = squeeze(vol(:, ijk(2), :))'; xv = y; yv = z;
    end
end

function cmap = blue_white_red()
%BLUE_WHITE_RED Diverging map: blue (low) -> white (middle) -> red (high).
%  Gamma on [0 2]: white = 1, red = fail. Differences on [-L L]: white = 0.
    t = linspace(0, 1, 128)';
    cmap = [t, t, ones(128, 1); ones(128, 1), flipud(t), flipud(t)];
end
