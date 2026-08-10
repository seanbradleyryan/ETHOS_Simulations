classdef SimPlotter
%SIMPLOTTER  All ETHOS result figures, one static method per figure type.
%
%   PURPOSE:
%   The single home for the plotting primitives that were copy-pasted (as file-
%   local functions, hence un-callable) across run_*/study_* scripts. Every
%   method takes an ethos.SimResult (and an optional gamma struct from
%   ethos.Analysis), so the same figures render whether the dose was simulated
%   or loaded.
%
%   METHODS:
%       all(result, G)                 - reconVsTruth + sensorPlacement (+convergence).
%       reconVsTruth(result, G, ...)   - truth | recon | signed diff [| gamma] axial.
%       changePanels(rA, rB, G, ...)   - reconA | reconB | diff [| gamma] axial.
%       doseDifference(result, G)      - alias to reconVsTruth (diff panel included).
%       gammaMap(result, G)            - standalone gamma-index map at truth max slice.
%       sensorPlacement(result)        - 3 orthogonal max-projections + array overlay.
%       convergence(result)            - TR convergence history (sim-only).
%
%   See also: ethos.SimResult, ethos.Analysis, ethos.StudyRunner

    methods (Static)
        function all(result, G)
        %ALL  Full standalone-style panel set for one result.
            if nargin < 2, G = []; end
            ethos.SimPlotter.reconVsTruth(result, G);
            ethos.SimPlotter.sensorPlacement(result);
            if result.isSimulated()
                ethos.SimPlotter.convergence(result);
            end
        end

        function reconVsTruth(result, G, labels, smooth_sigma)
        %RECONVSTRUTH  Axial truth | recon | (truth - recon) [| gamma].
            if nargin < 2, G = []; end
            if nargin < 3 || isempty(labels), labels = {'RayStation truth', 'Recon'}; end
            if nargin < 4 || isempty(smooth_sigma), smooth_sigma = 0; end

            volA = result.rs_truth;
            volB = result.recon_dose;
            if isempty(volA), volA = volB; labels{1} = labels{2}; end
            sliceRef = result.rs_truth;
            if isempty(sliceRef), sliceRef = volB; end

            [gmap, gcrit, gpass] = ethos.SimPlotter.unpackGamma(G);
            ttl = ethos.SimPlotter.titleFor(result, 'Recon vs Truth');
            ethos.SimPlotter.diffAxial(volA, volB, sliceRef, result.density, ...
                result.sensor_mask, result.spacing(), labels, ttl, smooth_sigma, ...
                gmap, gcrit, gpass);
        end

        function changePanels(resultA, resultB, G, smooth_sigma)
        %CHANGEPANELS  Axial reconA | reconB | (reconA - reconB) [| gamma].
            if nargin < 3, G = []; end
            if nargin < 4 || isempty(smooth_sigma), smooth_sigma = 0; end
            sliceRef = resultA.rs_truth;
            if isempty(sliceRef), sliceRef = resultA.recon_dose; end
            [gmap, gcrit, gpass] = ethos.SimPlotter.unpackGamma(G);
            ttl = ethos.SimPlotter.titleFor(resultA, 'Change detection (CT_1 vs CT_3)');
            ethos.SimPlotter.diffAxial(resultA.recon_dose, resultB.recon_dose, ...
                sliceRef, resultA.density, resultA.sensor_mask, resultA.spacing(), ...
                {'Recon CT_1', 'Recon CT_3'}, ttl, smooth_sigma, gmap, gcrit, gpass);
        end

        function doseDifference(result, G)
        %DOSEDIFFERENCE  Convenience: the truth|recon|diff panel (diff included).
            if nargin < 2, G = []; end
            ethos.SimPlotter.reconVsTruth(result, G);
        end

        function gammaMap(result, G)
        %GAMMAMAP  Standalone gamma-index map at the truth max-dose slice.
            if nargin < 2 || isempty(G) || ~isfield(G, 'gamma_map')
                error('ethos:SimPlotter:NoGamma', ...
                    'gammaMap requires a gamma struct from ethos.Analysis.gamma.');
            end
            gmap  = G.gamma_map;
            ref   = result.rs_truth;
            if isempty(ref), ref = result.recon_dose; end
            gridSize = size(gmap);
            [~, midx]  = max(ref(:));
            [~, ~, cz] = ind2sub(size(ref), midx);
            cz = min(max(cz, 1), gridSize(3));

            sp   = result.spacing();
            x_ax = (1:gridSize(1)) * sp(1);
            y_ax = (1:gridSize(2)) * sp(2);
            g_sl = squeeze(gmap(:, :, cz))';

            figure('Name', 'Gamma Index Map', 'Color', 'w', 'NumberTitle', 'off');
            ax = axes;
            ethos.SimPlotter.drawBackground(ax, result.density, cz, x_ax, y_ax);
            hold(ax, 'on');
            h = imagesc(ax, x_ax, y_ax, g_sl, [0, 2]);
            h.AlphaData = ~isnan(g_sl);
            colormap(ax, ethos.SimPlotter.gammaGYR(256));
            ethos.SimPlotter.overlaySensor(ax, result.sensor_mask, cz, x_ax, y_ax, 'k');
            hold(ax, 'off');
            axis(ax, 'xy'); axis(ax, 'image'); caxis(ax, [0, 2]);
            cb = colorbar(ax); cb.Label.String = '\gamma index';
            xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)');
            if isfield(G, 'criteria')
                title(ax, sprintf('\\gamma  %g%%/%gmm  (pass %.1f%%)  Z=%d', ...
                    G.criteria(1), G.criteria(2), G.pass_rate, cz));
            else
                title(ax, sprintf('\\gamma index  Z=%d', cz));
            end
            drawnow;
        end

        function sensorPlacement(result)
        %SENSORPLACEMENT  1x3 orthogonal max-projection view of array vs dose mask.
            dose = result.rs_truth;
            if isempty(dose), dose = result.recon_dose; end
            if isempty(dose)
                warning('ethos:SimPlotter:NoDose', 'No dose volume to plot.');
                return;
            end
            dose_mask = dose >= 0.10 * max(dose(:));
            sensor_mask = result.sensor_mask;
            if isempty(sensor_mask) || ~isequal(size(sensor_mask), size(dose_mask))
                warning('ethos:SimPlotter:NoSensor', ...
                    ['No co-gridded sensor mask on this result; use ' ...
                     'StudyRunner.ensureSensor first. Showing dose mask only.']);
                sensor_mask = false(size(dose_mask));
            end
            method = 'unknown';
            if isfield(result.sensor_info, 'method'), method = result.sensor_info.method; end
            ethos.SimPlotter.sensorPlanes(dose_mask, sensor_mask, result.spacing(), ...
                result.density, method);
        end

        function convergence(result)
        %CONVERGENCE  Time-reversal convergence history (sim-only).
            prov = result.provenance;
            has_p = isfield(prov, 'conv_max_pressure') && ~isempty(prov.conv_max_pressure);
            has_r = isfield(prov, 'conv_rel_change')   && ~isempty(prov.conv_rel_change);
            if ~has_p && ~has_r
                return;   % nothing recorded (loaded result or DAS)
            end
            figure('Name', 'TR Convergence', 'Color', 'w', 'NumberTitle', 'off');
            if has_p
                subplot(1, 2, 1);
                plot(prov.conv_max_pressure, '-o', 'LineWidth', 1.5);
                grid on; xlabel('Iteration'); ylabel('Max pressure (Pa)');
                title('Max reconstructed pressure');
            end
            if has_r
                subplot(1, 2, 2);
                semilogy(prov.conv_rel_change, '-o', 'LineWidth', 1.5);
                grid on; xlabel('Iteration'); ylabel('Relative change');
                title('TR relative change');
            end
            drawnow;
        end
    end

    % =====================================================================
    %  PRIVATE RENDERING PRIMITIVES (ported from the canonical study helpers)
    % =====================================================================
    methods (Static, Access = private)

        function diffAxial(volA, volB, slice_ref_vol, density, sensor_mask, ...
                spacing_mm, labels, titleStr, smooth_sigma, gamma_map, gamma_crit, gamma_pass)
        %DIFFAXIAL  volA | volB | (volA - volB) [| gamma] at the truth max slice.
            if nargin < 9  || isempty(smooth_sigma), smooth_sigma = 0; end
            if nargin < 10, gamma_map = []; end
            if nargin < 11, gamma_crit = []; end
            if nargin < 12, gamma_pass = NaN; end

            if ~isequal(size(volA), size(volB))
                warning('ethos:SimPlotter:SizeMismatch', ...
                    'volA %s and volB %s differ; skipping panel.', ...
                    mat2str(size(volA)), mat2str(size(volB)));
                return;
            end
            gridSize   = size(volA);
            have_gamma = ~isempty(gamma_map) && isequal(size(gamma_map), gridSize);
            n_panels   = 3 + double(have_gamma);
            [~, midx]  = max(slice_ref_vol(:));
            [~, ~, cz] = ind2sub(size(slice_ref_vol), midx);
            cz = min(max(cz, 1), gridSize(3));

            x_ax = (1:gridSize(1)) * spacing_mm(1);
            y_ax = (1:gridSize(2)) * spacing_mm(2);

            cmap_jet  = jet(256);
            wl_center = 1050; wl_width = 350;
            wl_min    = wl_center - wl_width / 2;

            have_density = ~isempty(density) && isequal(size(density), gridSize);
            if have_density, dens_sl = squeeze(density(:, :, cz))'; else, dens_sl = []; end
            if ~isempty(sensor_mask) && isequal(size(sensor_mask), gridSize)
                sensor_sl = squeeze(sensor_mask(:, :, cz))';
            else
                sensor_sl = [];
            end

            a_sl = squeeze(volA(:, :, cz))';
            b_sl = squeeze(volB(:, :, cz))';

            figure('Name', titleStr, 'Color', 'w', 'NumberTitle', 'off', ...
                'Position', [60, 80, 460 * n_panels, 460]);
            sgtitle(sprintf('%s\nAxial slice at truth max-dose voxel (Z=%d)', ...
                titleStr, cz), 'FontWeight', 'bold', 'FontSize', 11, ...
                'Interpreter', 'none');

            rmA = ethos.SimPlotter.robustDoseMax(volA, 100);
            rmB = ethos.SimPlotter.robustDoseMax(volB, 100);

            ax1 = subplot(1, n_panels, 1);
            ethos.SimPlotter.renderDosePanel(ax1, a_sl, dens_sl, sensor_sl, x_ax, y_ax, ...
                rmA, cmap_jet, wl_min, wl_width, smooth_sigma, ...
                sprintf('%s  Axial (Z=%d)', labels{1}, cz));

            ax2 = subplot(1, n_panels, 2);
            ethos.SimPlotter.renderDosePanel(ax2, b_sl, dens_sl, sensor_sl, x_ax, y_ax, ...
                rmB, cmap_jet, wl_min, wl_width, smooth_sigma, ...
                sprintf('%s  Axial (Z=%d)', labels{2}, cz));

            % Panel 3: signed difference.
            ax3 = subplot(1, n_panels, 3);
            diff_sl = a_sl - b_sl;
            if smooth_sigma > 0, diff_sl = imgaussfilt(diff_sl, smooth_sigma); end
            dmax = ethos.SimPlotter.robustDoseMax(abs(diff_sl), 100);
            if ~(dmax > 0), dmax = max(abs(diff_sl(:))); end
            if ~(dmax > 0), dmax = 1; end

            bg_rgb = ethos.SimPlotter.densBG(dens_sl, size(diff_sl), wl_min, wl_width);
            image(ax3, x_ax, y_ax, bg_rgb); hold(ax3, 'on');
            cmap_div = ethos.SimPlotter.blueWhiteRed(256);
            norm_d   = max(-1, min(1, diff_sl / dmax));
            idx      = max(1, min(256, round((norm_d + 1) / 2 * 255) + 1));
            diff_rgb = reshape(cmap_div(idx(:), :), [size(diff_sl), 3]);
            h_diff   = image(ax3, x_ax, y_ax, diff_rgb);
            h_diff.AlphaData = min(1, abs(diff_sl) / (0.30 * dmax));
            if ~isempty(sensor_sl) && any(sensor_sl(:))
                contour(ax3, x_ax, y_ax, double(sensor_sl), [0.5, 0.5], 'k-', 'LineWidth', 1.5);
            end
            hold(ax3, 'off'); axis(ax3, 'xy'); axis(ax3, 'image');
            colormap(ax3, cmap_div); caxis(ax3, [-dmax, dmax]);
            cb = colorbar(ax3); cb.Label.String = 'Difference (Gy)';
            xlabel(ax3, 'X (mm)'); ylabel(ax3, 'Y (mm)');
            title(ax3, sprintf('%s - %s  Axial (Z=%d)', labels{1}, labels{2}, cz));

            % Panel 4: gamma.
            if have_gamma
                ax4  = subplot(1, n_panels, 4);
                g_sl = squeeze(gamma_map(:, :, cz))';
                bg   = ethos.SimPlotter.densBG(dens_sl, size(g_sl), wl_min, wl_width);
                image(ax4, x_ax, y_ax, bg); hold(ax4, 'on');
                h_g  = imagesc(ax4, x_ax, y_ax, g_sl, [0, 2]);
                h_g.AlphaData = ~isnan(g_sl);
                colormap(ax4, ethos.SimPlotter.gammaGYR(256));
                if ~isempty(sensor_sl) && any(sensor_sl(:))
                    contour(ax4, x_ax, y_ax, double(sensor_sl), [0.5, 0.5], 'k-', 'LineWidth', 1.5);
                end
                hold(ax4, 'off'); axis(ax4, 'xy'); axis(ax4, 'image'); caxis(ax4, [0, 2]);
                cb4 = colorbar(ax4); cb4.Label.String = '\gamma index';
                xlabel(ax4, 'X (mm)'); ylabel(ax4, 'Y (mm)');
                if ~isempty(gamma_crit)
                    title(ax4, sprintf('\\gamma  %g%%/%gmm  (pass %.1f%%)  Z=%d', ...
                        gamma_crit(1), gamma_crit(end), gamma_pass, cz));
                else
                    title(ax4, sprintf('\\gamma index  (pass %.1f%%)  Z=%d', gamma_pass, cz));
                end
            end
            drawnow;
        end

        function renderDosePanel(ax, dose_slice, density_slice, sensor_slice, xv, yv, ...
                row_max, cmap_jet, wl_min, wl_width, smooth_sigma, ttl)
            dose_slice = double(dose_slice);
            if smooth_sigma > 0, dose_slice = imgaussfilt(dose_slice, smooth_sigma); end
            bg_rgb = ethos.SimPlotter.densBG(density_slice, size(dose_slice), wl_min, wl_width);
            image(ax, xv, yv, bg_rgb); hold(ax, 'on');

            norm_d   = max(0, min(1, dose_slice / row_max));
            idx      = max(1, min(256, round(norm_d * 255) + 1));
            dose_rgb = reshape(cmap_jet(idx(:), :), [size(dose_slice), 3]);
            above    = dose_slice >= 0.10 * row_max;
            ramp     = 0.45 + 0.35 * min(1, (dose_slice - 0.10 * row_max) / max(0.40 * row_max, 1e-12));
            h_dose   = image(ax, xv, yv, dose_rgb);
            h_dose.AlphaData = above .* ramp;

            if ~isempty(sensor_slice) && any(sensor_slice(:))
                contour(ax, xv, yv, double(sensor_slice), [0.5, 0.5], 'r-', 'LineWidth', 2);
            end
            hold(ax, 'off'); axis(ax, 'xy'); axis(ax, 'image');
            colormap(ax, cmap_jet); caxis(ax, [0, row_max]);
            cb = colorbar(ax); cb.Label.String = 'Dose (Gy)';
            xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); title(ax, ttl);
        end

        function sensorPlanes(dose_mask, sensor_mask, spacing_mm, density, method)
            [Nx3, Ny3, Nz3] = size(dose_mask);
            x_ax = (1:Nx3) * spacing_mm(1);
            y_ax = (1:Ny3) * spacing_mm(2);
            z_ax = (1:Nz3) * spacing_mm(3);

            dose_cor = squeeze(any(dose_mask, 2)); dose_sag = squeeze(any(dose_mask, 1)); dose_axi = squeeze(any(dose_mask, 3));
            sens_cor = squeeze(any(sensor_mask, 2)); sens_sag = squeeze(any(sensor_mask, 1)); sens_axi = squeeze(any(sensor_mask, 3));

            have_density = ~isempty(density) && isequal(size(density), size(dose_mask));
            wl_center = 1050; wl_width = 350; wl_min = wl_center - wl_width / 2;
            if have_density
                ct_projs = {squeeze(mean(double(density), 2))', ...
                            squeeze(mean(double(density), 1))', ...
                            squeeze(mean(double(density), 3))'};
            end

            figure('Name', 'Sensor Placement vs Dose Mask', 'Color', 'w', ...
                'NumberTitle', 'off', 'Position', [80, 80, 1300, 420]);
            sgtitle(sprintf('Sensor Placement vs Dose Mask  (\\geq10%% max)   |   Sensor: %s', ...
                method), 'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');

            view_data = {
                dose_cor', sens_cor', x_ax, z_ax, 'X (mm)', 'Z (mm)', 'Coronal  (mean-proj along Y)';
                dose_sag', sens_sag', y_ax, z_ax, 'Y (mm)', 'Z (mm)', 'Sagittal  (mean-proj along X)';
                dose_axi', sens_axi', x_ax, y_ax, 'X (mm)', 'Y (mm)', 'Axial  (mean-proj along Z)'};
            dose_color = [0.20, 0.50, 0.90]; sensor_color = [0.90, 0.10, 0.10];

            for col = 1:3
                ax  = subplot(1, 3, col);
                d2d = double(view_data{col, 1}); s2d = double(view_data{col, 2});
                xv  = view_data{col, 3}; yv = view_data{col, 4};
                if have_density
                    dn = max(0, min(1, (ct_projs{col} - wl_min) / wl_width));
                    bg_rgb = repmat(dn, [1, 1, 3]);
                else
                    bg_rgb = ones([size(d2d), 3]);
                end
                image(ax, xv, yv, bg_rgb); hold(ax, 'on');
                if any(d2d(:))
                    h_dose = image(ax, xv, yv, repmat(reshape(dose_color, [1,1,3]), [size(d2d), 1]));
                    h_dose.AlphaData = d2d * 0.45;
                end
                if any(s2d(:))
                    h_sens = image(ax, xv, yv, repmat(reshape(sensor_color, [1,1,3]), [size(s2d), 1]));
                    h_sens.AlphaData = s2d * 0.85;
                end
                hold(ax, 'off'); axis(ax, 'xy'); axis(ax, 'image');
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

        function bg_rgb = densBG(dens_sl, target_size, wl_min, wl_width)
            if ~isempty(dens_sl) && isequal(size(dens_sl), target_size)
                dn = max(0, min(1, (dens_sl - wl_min) / wl_width));
                bg_rgb = repmat(dn, [1, 1, 3]);
            else
                bg_rgb = zeros([target_size, 3]);
            end
        end

        function drawBackground(ax, density, cz, x_ax, y_ax)
            if ~isempty(density) && cz <= size(density, 3)
                dens_sl = squeeze(density(:, :, cz))';
                bg = ethos.SimPlotter.densBG(dens_sl, size(dens_sl), 875, 350);
            else
                bg = zeros([numel(y_ax), numel(x_ax), 3]);
            end
            image(ax, x_ax, y_ax, bg);
        end

        function overlaySensor(ax, sensor_mask, cz, x_ax, y_ax, color)
            if ~isempty(sensor_mask) && cz <= size(sensor_mask, 3)
                s_sl = squeeze(sensor_mask(:, :, cz))';
                if any(s_sl(:))
                    contour(ax, x_ax, y_ax, double(s_sl), [0.5, 0.5], [color '-'], 'LineWidth', 1.5);
                end
            end
        end

        function [gmap, gcrit, gpass] = unpackGamma(G)
            gmap = []; gcrit = []; gpass = NaN;
            if ~isempty(G) && isstruct(G) && isfield(G, 'gamma_map')
                gmap = G.gamma_map;
                if isfield(G, 'criteria'),  gcrit = G.criteria;  end
                if isfield(G, 'pass_rate'), gpass = G.pass_rate; end
            end
        end

        function s = titleFor(result, prefix)
            id = strtrim(sprintf('%s / %s', result.patient_id, result.session));
            lbl = '';
            if isfield(result.rtplan, 'beam_num') && ~isempty(result.rtplan.beam_num) ...
                    && ~isnan(result.rtplan.beam_num)
                lbl = sprintf('  B%g', result.rtplan.beam_num);
                if isfield(result.rtplan, 'seg_num') && ~isnan(result.rtplan.seg_num)
                    lbl = sprintf('%s_%g', lbl, result.rtplan.seg_num);
                end
            end
            s = sprintf('%s   %s%s', prefix, id, lbl);
        end

        function cmap = gammaGYR(n)
            if nargin < 1 || isempty(n), n = 256; end
            h = floor(n / 2);
            green = [0.00, 0.60, 0.20]; yellow = [1.00, 0.95, 0.20]; red = [0.80, 0.10, 0.10];
            half1 = [linspace(green(1), yellow(1), h)', linspace(green(2), yellow(2), h)', linspace(green(3), yellow(3), h)'];
            half2 = [linspace(yellow(1), red(1), n - h)', linspace(yellow(2), red(2), n - h)', linspace(yellow(3), red(3), n - h)'];
            cmap = [half1; half2];
        end

        function cmap = blueWhiteRed(n)
            if nargin < 1 || isempty(n), n = 256; end
            bot = [0.230, 0.299, 0.754]; mid = [1, 1, 1]; top = [0.706, 0.016, 0.150];
            h = floor(n / 2);
            half1 = [linspace(bot(1), mid(1), h)', linspace(bot(2), mid(2), h)', linspace(bot(3), mid(3), h)'];
            half2 = [linspace(mid(1), top(1), n - h)', linspace(mid(2), top(2), n - h)', linspace(mid(3), top(3), n - h)'];
            cmap = [half1; half2];
        end

        function m = robustDoseMax(d, clip_pct)
            if nargin < 2 || isempty(clip_pct), clip_pct = 100; end
            vals = sort(double(d(d > 0)));
            if isempty(vals), m = 1; return; end
            if clip_pct >= 100
                m = vals(end);
            else
                k = max(1, min(numel(vals), ceil(clip_pct / 100 * numel(vals))));
                m = vals(k);
            end
            if ~(m > 0), m = vals(end); end
            if ~(m > 0), m = 1; end
        end
    end
end
