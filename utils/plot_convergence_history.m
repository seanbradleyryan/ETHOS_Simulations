function plot_convergence_history(conv_max_pressure, conv_rel_change, num_iters, tol, p0_max_orig, parent)
%PLOT_CONVERGENCE_HISTORY p0 convergence over TR iterations.
%  Left axis:  max reconstructed p0 per iteration (blue).
%              Dashed black line = max of original p0 distribution (reference).
%  Right axis: relative change between iterations (red, log-scale).
%  Optional parent (figure / uitab / uipanel) hosts the plot; when omitted a
%  new figure is created. Passing a uitab lets the caller collect several
%  panels as tabs of one window.
%  (Shared copy; formerly a local in the standalone drivers.)

    if nargin < 6, parent = []; end

    iters   = 1:num_iters;
    p_vals  = conv_max_pressure(iters);
    rc_vals = conv_rel_change(iters);

    if isempty(parent)
        figure('Name', 'p0 Convergence', 'Color', 'w', ...
            'NumberTitle', 'off', 'Position', [150, 520, 720, 390]);
    else
        axes('Parent', parent);   % draw into the supplied tab / panel
    end

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
