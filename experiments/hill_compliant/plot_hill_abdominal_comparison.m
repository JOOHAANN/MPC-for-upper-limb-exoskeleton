function plot_hill_abdominal_comparison(r, models, out)
    % PLOT_HILL_ABDOMINAL_COMPARISON Overlay both sources on common axes.
    % Torques are cached Hill predictions; no fitting, torque reference or new
    % simulation is introduced. Base inputs precede excitation feedback.
    % frontTorque is positive flexion magnitude; the joint uses -frontTorque.
    t = r.sequence;
    p = r.predictions{1};
    g = r.predictions{2};
    experimentalBase = hill_pair_inputs(t, 1, models{1}.generator);
    mapBase = hill_pair_inputs(t, 2, models{2}.generator);
    colors = [.05 .35 .75; .05 .55 .35];
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1700 1450]);
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile([1 2]);
    comparison(ax, t, p.frontTorque, g.frontTorque, colors, ...
               'Abdominal flexion torque (positive magnitude; opposite sign in net joint torque)', 'Torque (N m)', true);
    ax = nexttile([1 2]);
    comparison(ax, t, p.effectiveActivation(:, 2), g.effectiveActivation(:, 2), colors, ...
               'Abdominal Hill activation after muscle excitation control (not verified %MVC)', 'Activation', false);
    ax = nexttile;
    comparison(ax, t, p.frontActive, g.frontActive, colors, ...
               'Active CE torque', 'Torque (N m)', false);
    ax = nexttile;
    comparison(ax, t, p.frontPassive, g.frontPassive, colors, ...
               'Passive PE torque', 'Torque (N m)', false);
    ax = nexttile;
    comparison(ax, t, p.frontFiberDamping, g.frontFiberDamping, colors, ...
               'Signed fiber-damping torque', 'Torque (N m)', false);
    xlabel(ax, 'Time (s)');
    ax = nexttile;
    hold(ax, 'on');
    plot(ax, t.t, experimentalBase(:, 2), 'Color', colors(1, :), 'LineWidth', 1.2);
    plot(ax, t.t, mapBase(:, 2), 'Color', colors(2, :), 'LineWidth', 1.2);
    title(ax, 'Abdominal base inputs before recruitment and muscle feedback');
    ylabel(ax, 'Input activation');
    xlabel(ax, 'Time (s)');
    legend(ax, 'Acquisition-derived rect_abd_r input', 'Auxiliary prediction from MAP back input', ...
           'Location', 'northoutside', 'FontSize', 9, 'Interpreter', 'none');
    grid(ax, 'on');
    hill_sequence_protocol(ax, t, false);
    sgtitle(sprintf('Abdominal torque | Group %d | Experimental vs MAP | Positive flexion magnitude | Model predictions', r.group));
    name = sprintf('abdominal_torque_comparison_group%d', r.group);
    exportgraphics(fig, fullfile(out, [name '.png']), 'Resolution', 170);
    exportgraphics(fig, fullfile(out, [name '.pdf']), 'ContentType', 'vector');
    savefig(fig, fullfile(out, [name '.fig']));
    close(fig);
    % Preserve the exact plotted values, including small passive/damping terms.
    traces = table(t.t, t.loadKg, t.waiting, p.frontTorque, g.frontTorque, ...
                   p.effectiveActivation(:, 2), g.effectiveActivation(:, 2), ...
                   p.frontActive, g.frontActive, p.frontPassive, g.frontPassive, ...
                   p.frontFiberDamping, g.frontFiberDamping, experimentalBase(:, 2), mapBase(:, 2), ...
                   -p.frontTorque, -g.frontTorque, ...
                   'VariableNames', {'time', 'nominalLoadKg', 'waiting', ...
                          'experimentalAbdominalTorqueNm', 'mapAbdominalTorqueNm', ...
                          'experimentalHillActivation', 'mapHillActivation', ...
                          'experimentalActiveCETorqueNm', 'mapActiveCETorqueNm', ...
                          'experimentalPassivePETorqueNm', 'mapPassivePETorqueNm', ...
                          'experimentalFiberDampingTorqueNm', 'mapFiberDampingTorqueNm', ...
                          'experimentalBaseFrontInput', 'mapBaseFrontInput', ...
                          'experimentalAbdominalJointContributionNm', 'mapAbdominalJointContributionNm'});
    writetable(traces, fullfile(out, [name '.csv']));
end

function comparison(ax, t, experimental, mapped, colors, label, units, labels)
    hold(ax, 'on');
    plot(ax, t.t, experimental, 'Color', colors(1, :), 'LineWidth', 1.2);
    plot(ax, t.t, mapped, 'Color', colors(2, :), 'LineWidth', 1.2);
    title(ax, label);
    ylabel(ax, units);
    legend(ax, 'Experimental input + muscle feedback', 'MAP input + muscle feedback', ...
           'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 10);
    grid(ax, 'on');
    if strcmp(units, 'Activation')
        if max([experimental(:); mapped(:)]) <= .15
            ylim(ax, [0 .15]);
        else
            ylim(ax, [0 1]);
        end
    end
    hill_sequence_protocol(ax, t, labels);
end
