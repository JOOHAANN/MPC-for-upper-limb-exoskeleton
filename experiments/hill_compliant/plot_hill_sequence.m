function plot_hill_sequence(r, out)
    % PLOT_HILL_SEQUENCE All motion and idle samples, without clipping/rescaling.
    t = r.sequence;
    p = r.predictions{1};
    g = r.predictions{2};
    colors = [.80 .15 .10; .05 .35 .75; .05 .55 .35];
    labels = {'Experimental + stance activation', 'MAP + stance activation'};
    if isfield(p, 'motionIntent')
        labels = {'Experimental + neural feedback', 'MAP + neural feedback'};
    end
    plot_torque_theta(r, out, colors, labels);
    plot_hill_abdominal_torque(r, out);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1700 1200]);
    tiledlayout(4, 2, 'TileSpacing', 'compact');
    ax = nexttile([1 2]);
    hold(ax, 'on');
    plot(t.t, t.torque, '--', 'Color', colors(1, :), 'LineWidth', 1.2);
    plot(t.t, p.torque, 'Color', colors(2, :), 'LineWidth', 1.4);
    plot(t.t, g.torque, 'Color', colors(3, :), 'LineWidth', 1.4);
    ylabel('Net Hill torque (Nm)');
    legend('Original reference (idle 20 Nm is filled)', labels{1}, labels{2}, ...
           'Location', 'northoutside', 'Orientation', 'horizontal');
    grid on;
    show_protocol(ax, t, true);
    ax = nexttile([1 2]);
    hold(ax, 'on');
    plot(t.t, t.q, '--', 'Color', colors(1, :), 'LineWidth', 1.2);
    plot(t.t, p.q, 'Color', colors(2, :), 'LineWidth', 1.4);
    plot(t.t, g.q, 'Color', colors(3, :), 'LineWidth', 1.4);
    ylabel('Trunk angle (rad)');
    ylim([min(0, min([p.q; g.q])), max(1.9, max([p.q; g.q]))]);
    legend('Original Scale reference', labels{1}, labels{2}, ...
           'Location', 'northoutside', 'Orientation', 'horizontal');
    grid on;
    show_protocol(ax, t, false);
    for source = 1:2
        ax = nexttile;
        p = r.predictions{source};
        hold(ax, 'on');
        plot(t.t, r.inputs{source}(:, 1), ':', 'Color', [.35 .35 .35], 'LineWidth', 1);
        plot(t.t, p.effectiveActivation(:, 1), 'Color', [.05 .35 .75], 'LineWidth', 1.1);
        plot(t.t, p.effectiveActivation(:, 2), 'Color', [.85 .35 .10], 'LineWidth', 1.1);
        plot(t.t, p.balanceWeight, '--', 'Color', [.50 .25 .65], 'LineWidth', .8);
        title(sprintf('Source %d: input and actual Hill activation', source));
        ylabel('Activation / stance weight');
        legend('Back input', 'Back Hill activation', 'Front Hill activation', 'Stance blend', ...
               'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 9);
        grid on;
        show_protocol(ax, t, false);
    end
    for source = 1:2
        ax = nexttile;
        p = r.predictions{source};
        hold(ax, 'on');
        plot(t.t, p.gravityTorque, 'Color', [.35 .35 .35], 'LineWidth', 1);
        plot(t.t, p.dampingTorque, 'Color', [.05 .55 .35], 'LineWidth', 1);
        plot(t.t, p.boundaryTorque, 'Color', [.50 .25 .65], 'LineWidth', 1);
        plot(t.t, p.externalBalanceTorque, 'k--', 'LineWidth', 1.2);
        ylabel('Other joint torques (Nm)');
        xlabel('Time (s)');
        legend('Gravity', 'Joint damping', 'Soft limit', 'External balance = 0', ...
               'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 9);
        grid on;
        show_protocol(ax, t, false);
    end
    sgtitle(sprintf('Continuous 5 / 10 / 15 kg sequence %d | Original waits | Same Hill model | No state resets', r.group));
    exportgraphics(fig, fullfile(out, sprintf('sequence_group%d.png', r.group)), 'Resolution', 170);
    exportgraphics(fig, fullfile(out, sprintf('sequence_group%d.pdf', r.group)), 'ContentType', 'vector');
    savefig(fig, fullfile(out, sprintf('sequence_group%d.fig', r.group)));
    close(fig);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [30 30 1550 950]);
    tiledlayout(2, 2, 'TileSpacing', 'compact');
    for source = 1:2
        ax = nexttile;
        p = r.predictions{source};
        b = r.unbalanced{source};
        hold(ax, 'on');
        plot(t.t, p.q, 'LineWidth', 1.4, 'DisplayName', 'Muscle feedback');
        if ~b.failed
            plot(t.t, b.q, '--', 'Color', [.55 .55 .55], 'LineWidth', 1.1, 'DisplayName', 'Feedback disabled');
        else
            text(ax, .02, .06, 'Disabled-feedback rollout failed; curve omitted', 'Units', 'normalized', 'Color', [.7 .1 .1]);
        end
        yline(t.standingAngle, ':', 'Original standing target', 'HandleVisibility', 'off');
        title(sprintf('Source %d: posture controller enabled / disabled', source));
        ylabel('Angle (rad)');
        legend('show', 'Location', 'northoutside');
        grid on;
        show_protocol(ax, t, false);
    end
    for source = 1:2
        ax = nexttile;
        p = r.predictions{source};
        hold(ax, 'on');
        plot(t.t, p.neuralCommands, 'LineWidth', 1.1);
        plot(t.t, r.inputs{source}, ':', 'LineWidth', .8);
        title(sprintf('Source %d: added postural commands are visible', source));
        ylabel('Neural excitation');
        xlabel('Time (s)');
        legend('Back command', 'Front command', 'Back base input', 'Front mapped base input', ...
               'Location', 'northoutside', 'NumColumns', 2);
        grid on;
        show_protocol(ax, t, false);
    end
    sgtitle('Standing balance acts through Hill activation; external PD balance torque is zero');
    exportgraphics(fig, fullfile(out, sprintf('standing_balance_group%d.png', r.group)), 'Resolution', 170);
    exportgraphics(fig, fullfile(out, sprintf('standing_balance_group%d.pdf', r.group)), 'ContentType', 'vector');
    close(fig);
end

function plot_torque_theta(r, out, colors, labels)
    % A compact two-panel figure follows the supplied paper figure layout.
    t = r.sequence;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1550 800]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    ax = nexttile;
    hold(ax, 'on');
    plot(t.t, t.torque, '--', 'Color', colors(1, :), 'LineWidth', 1.2);
    plot(t.t, r.predictions{1}.torque, 'Color', colors(2, :), 'LineWidth', 1.5);
    plot(t.t, r.predictions{2}.torque, 'Color', colors(3, :), 'LineWidth', 1.5);
    ylabel('Net Hill torque (N m)');
    legend('Original reference (idle 20 N m is filled)', labels{1}, labels{2}, ...
           'Location', 'northoutside', 'Orientation', 'horizontal');
    grid on;
    show_protocol(ax, t, true);
    ax = nexttile;
    hold(ax, 'on');
    plot(t.t, t.q, '--', 'Color', colors(1, :), 'LineWidth', 1.2);
    plot(t.t, r.predictions{1}.q, 'Color', colors(2, :), 'LineWidth', 1.5);
    plot(t.t, r.predictions{2}.q, 'Color', colors(3, :), 'LineWidth', 1.5);
    ylabel('Trunk angle (rad)');
    xlabel('Time (s)');
    ylim([min(0, min([r.predictions{1}.q; r.predictions{2}.q])), max(1.9, max([r.predictions{1}.q; r.predictions{2}.q]))]);
    legend('Original Scale reference', labels{1}, labels{2}, ...
           'Location', 'northoutside', 'Orientation', 'horizontal');
    grid on;
    show_protocol(ax, t, false);
    sgtitle(sprintf('Continuous 5 / 10 / 15 kg | Group %d | Shaded: original standing waits | No state resets', r.group));
    save_outputs(fig, out, sprintf('torque_theta_group%d', r.group));
end

function save_outputs(fig, out, name)
    exportgraphics(fig, fullfile(out, [name '.png']), 'Resolution', 170);
    exportgraphics(fig, fullfile(out, [name '.pdf']), 'ContentType', 'vector');
    savefig(fig, fullfile(out, [name '.fig']));
    close(fig);
end

function show_protocol(ax, t, labels)
    hill_sequence_protocol(ax, t, labels);
end
