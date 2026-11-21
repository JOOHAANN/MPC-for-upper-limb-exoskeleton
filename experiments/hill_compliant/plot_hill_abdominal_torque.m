function plot_hill_abdominal_torque(r, out)
    % PLOT_HILL_ABDOMINAL_TORQUE Flexor torque and its Hill components in Nm.
    % The displayed total is positive; its signed joint contribution is negative.
    t = r.sequence;
    colors = [.05 .35 .75; .05 .55 .35];
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1700 950]);
    tiledlayout(2, 2, 'TileSpacing', 'compact');
    ax = nexttile([1 2]);
    hold(ax, 'on');
    for source = 1:2
        plot(ax, t.t, r.predictions{source}.frontTorque, 'Color', colors(source, :), 'LineWidth', 1.5);
    end
    ylabel(ax, 'Abdominal flexion torque (N m)');
    legend(ax, 'Experimental input + muscle feedback', 'MAP input + muscle feedback', ...
           'Location', 'northoutside', 'Orientation', 'horizontal');
    grid(ax, 'on');
    hill_sequence_protocol(ax, t, true);
    labels = {'Experimental activation', 'MAP activation'};
    for source = 1:2
        ax = nexttile;
        p = r.predictions{source};
        hold(ax, 'on');
        plot(ax, t.t, p.frontTorque, 'k', 'LineWidth', 1.5);
        plot(ax, t.t, p.frontActive, 'Color', [.05 .35 .75], 'LineWidth', 1.1);
        plot(ax, t.t, p.frontPassive, 'Color', [.85 .35 .10], 'LineWidth', 1.1);
        plot(ax, t.t, p.frontFiberDamping, 'Color', [.50 .25 .65], 'LineWidth', 1.1);
        title(ax, [labels{source} ': abdominal torque components']);
        ylabel(ax, 'Torque (N m)');
        xlabel(ax, 'Time (s)');
        legend(ax, 'Total abdominal torque', 'Active CE', 'Passive PE', 'Fiber damping (signed)', ...
               'Location', 'northoutside', 'NumColumns', 2);
        grid(ax, 'on');
        hill_sequence_protocol(ax, t, false);
    end
    sgtitle(sprintf('Group %d | Abdominal Hill torque | Positive flexion magnitude | Net joint uses the opposite sign', r.group));
    name = sprintf('abdominal_torque_group%d', r.group);
    exportgraphics(fig, fullfile(out, [name '.png']), 'Resolution', 170);
    exportgraphics(fig, fullfile(out, [name '.pdf']), 'ContentType', 'vector');
    savefig(fig, fullfile(out, [name '.fig']));
    close(fig);
end
