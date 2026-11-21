function summary = plot_hill_soft_boundary(root)
    % PLOT_HILL_SOFT_BOUNDARY Separate the elastic limit law from actual use.
    % All torques come from hill_soft_limit, the function used by the model.
    % Neither angle clipping nor a parameter adjustment is performed here.
    cache = load(fullfile(root, 'data', 'processed', 'hill_compliant', 'neural_models.mat'));
    out = fullfile(root, 'results', 'hill_neural_sequence');
    replayFile = fullfile(out, 'hill_sequence.mat');
    if ~isfile(replayFile)
        options = struct('modelFile', 'neural_models.mat', ...
                         'outputDirectory', 'hill_neural_sequence', 'renderFigures', false);
        result = run_hill_sequence(root, options);
    else
        replay = load(replayFile, 'result');
        result = replay.result;
    end
    for source = 1:2
        assert(isequal(cache.models{source}.parameters, result.models{source}.parameters) && ...
               isequal(cache.models{source}.neuralParameters, result.models{source}.neuralParameters), ...
               'Cached replay differs from the selected model. Run demo_hill_neural first.');
    end
    par = cache.models{1}.parameters;
    assert(isequal(par(29:32), cache.models{2}.parameters(29:32)), ...
           'The two activation sources must share the soft-limit settings.');
    stiffness = par(29);
    cap = par(30);
    lower = par(31);
    upper = par(32);
    sourceNames = {'Experimental activation', 'MAP activation'};

    % Recompute the elastic torque at every stored angle and check the replay.
    records = struct([]);
    for group = 1:3
        for source = 1:2
            p = result.groups{group}.predictions{source};
            expected = hill_soft_limit(p.q, cache.models{source}.parameters);
            residual = max(abs(expected(:) - p.boundaryTorque(:)));
            assert(residual < 1e-10, 'Stored soft-limit torque does not match the model law.');
            rec = struct('group', group, 'source', sourceNames{source}, ...
                         'minimumAngleRad', min(p.q), 'maximumAngleRad', max(p.q), ...
                         'peakSoftTorqueNm', max(abs(p.boundaryTorque)), ...
                         'activeSamples', nnz(abs(p.boundaryTorque) > 1e-10), ...
                         'lawResidualNm', residual);
            records = [records; rec]; %#ok<AGROW>
        end
    end
    summary = struct2table(records);
    angleMin = min(summary.minimumAngleRad);
    angleMax = max(summary.maximumAngleRad);

    % A penetration is measured beyond an onset angle, not from upright.
    penetrationRad = [.01; .02; .05; .10; .20; .30; .50; .80];
    lowerAngleRad = lower - penetrationRad;
    upperAngleRad = upper + penetrationRad;
    lowerTorqueNm = hill_soft_limit(lowerAngleRad, par);
    upperTorqueNm = hill_soft_limit(upperAngleRad, par);
    assert(max(abs(lowerTorqueNm + upperTorqueNm)) < 1e-10);
    samples = table(penetrationRad, lowerAngleRad, lowerTorqueNm, upperAngleRad, upperTorqueNm);
    writetable(samples, fullfile(out, 'soft_boundary_samples.csv'));
    writetable(summary, fullfile(out, 'soft_boundary_actual_summary.csv'));

    blue = [.10 .35 .75];
    orange = [.90 .38 .10];
    green = [.18 .55 .32];
    characteristic = figure('Name', 'Soft boundary: torque versus angle', ...
                            'Color', 'w', 'Position', [80 80 1100 790]);
    layout = tiledlayout(characteristic, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile(layout);
    hold(ax, 'on');
    deadZone = patch(ax, [lower upper upper lower], [-90 -90 90 90], [.94 .94 .94], ...
                     'EdgeColor', 'none');
    actualRange = patch(ax, [angleMin angleMax angleMax angleMin], [-90 -90 90 90], ...
                        [.84 .94 .86], 'EdgeColor', 'none');
    theta = linspace(lower - .80, upper + .80, 1801);
    law = plot(ax, theta, hill_soft_limit(theta, par), 'Color', blue, 'LineWidth', 2.4);
    yline(ax, cap, ':', 'Color', [.5 .5 .5], 'HandleVisibility', 'off');
    yline(ax, -cap, ':', 'Color', [.5 .5 .5], 'HandleVisibility', 'off');
    xline(ax, lower, '--', sprintf('Lower onset: %.2f rad', lower), ...
          'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
    xline(ax, upper, '--', sprintf('Upper onset: %.2f rad', upper), ...
          'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
    text(ax, lower - .75, 83, sprintf('Saturation: +%.0f N m', cap), 'FontSize', 10);
    text(ax, upper + .22, -72, sprintf('Saturation: -%.0f N m', cap), 'FontSize', 10);
    xlabel(ax, 'Trunk angle, \theta (rad)');
    ylabel(ax, 'Elastic soft-limit torque (N m)');
    title(ax, sprintf('Current elastic boundary: k = %.0f N m/rad; cap = %.0f N m', stiffness, cap));
    xlim(ax, [theta(1) theta(end)]);
    ylim(ax, [-90 90]);
    legend(ax, [law deadZone actualRange], {'Model soft-limit law', ...
           'Zero-torque interval', sprintf('Actual angle range: %.3f-%.3f rad', angleMin, angleMax)}, ...
           'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
    style_axes(ax);

    ax = nexttile(layout);
    hold(ax, 'on');
    d = linspace(0, .35, 601);
    plot(ax, d, hill_soft_limit(lower - d, par), 'Color', blue, 'LineWidth', 2.4);
    selected = [.05 .10 .20 .30];
    values = hill_soft_limit(lower - selected, par);
    scatter(ax, selected, values, 45, orange, 'filled');
    for k = 1:numel(selected)
        text(ax, selected(k) + .004, values(k) + 3, ...
             sprintf('%.2f rad: %.2f N m', selected(k), values(k)), 'FontSize', 10);
    end
    text(ax, .008, 62, 'Below lower onset: positive torque; above upper onset: negative torque.', ...
         'FontSize', 10);
    xlabel(ax, 'Angle beyond the onset, d (rad)');
    ylabel(ax, 'Restoring torque magnitude (N m)');
    title(ax, 'Smooth onset: the torque increases with penetration on either side');
    xlim(ax, [0 .35]);
    ylim(ax, [0 75]);
    style_axes(ax);
    save_figure(characteristic, out, 'soft_boundary_characteristic');

    % Each sequence contains all three loads and their standing gaps.
    actual = figure('Name', 'Actual soft-boundary torque: both activation sources', ...
                    'Color', 'w', 'Position', [100 60 1180 900]);
    layout = tiledlayout(actual, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    colors = {blue, orange};
    lineStyles = {'-', '--'};
    for group = 1:3
        g = result.groups{group};
        ax = nexttile(layout);
        hold(ax, 'on');
        for source = 1:2
            plot(ax, g.sequence.t, g.predictions{source}.q, 'Color', colors{source}, ...
                 'LineStyle', lineStyles{source}, 'LineWidth', 1.6);
        end
        yline(ax, lower, ':', sprintf('%.2f rad', lower), 'HandleVisibility', 'off');
        yline(ax, upper, ':', sprintf('%.2f rad', upper), 'HandleVisibility', 'off');
        title(ax, sprintf('Sequence %d: 5 / 10 / 15 kg and standing gaps', group));
        ylabel(ax, '\theta (rad)');
        ylim(ax, [min(lower -.08, angleMin -.05), max(upper + .08, angleMax + .05)]);
        xlim(ax, [g.sequence.t(1) g.sequence.t(end)]);
        if group == 1
            legend(ax, sourceNames, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
        end
        if group == 3
            xlabel(ax, 'Time (s)');
        end
        style_axes(ax);
        ax = nexttile(layout);
        hold(ax, 'on');
        for source = 1:2
            plot(ax, g.sequence.t, g.predictions{source}.boundaryTorque, ...
                 'Color', colors{source}, 'LineStyle', lineStyles{source}, 'LineWidth', 1.8);
        end
        peak = max(summary.peakSoftTorqueNm(summary.group == group));
        span = max(1, 1.15 * peak);
        title(ax, sprintf('Actual soft-limit torque: peak |\tau| = %.3g N m', peak));
        if peak == 0
            text(ax, 1, .52 * span, 'Both sources overlap at zero.', 'FontSize', 11, 'Color', green);
        end
        text(ax, 1, .25 * span, 'Elastic boundary only; joint damping and muscle torque are separate.', 'FontSize', 10);
        ylabel(ax, '\tau_{soft} (N m)');
        ylim(ax, [-span span]);
        xlim(ax, [g.sequence.t(1) g.sequence.t(end)]);
        if group == 3
            xlabel(ax, 'Time (s)');
        end
        style_axes(ax);
    end
    save_figure(actual, out, 'soft_boundary_actual_torque');
    fprintf('Current actual soft-limit peak: %.6g N m.\n', max(summary.peakSoftTorqueNm));
    fprintf('Figures and numeric samples saved to: %s\n', out);
end

function style_axes(ax)
    set(ax, 'FontName', 'Arial', 'FontSize', 11, 'Layer', 'top');
    grid(ax, 'on');
    box(ax, 'on');
end

function save_figure(fig, out, name)
    exportgraphics(fig, fullfile(out, [name '.png']), 'Resolution', 180);
    exportgraphics(fig, fullfile(out, [name '.pdf']), 'ContentType', 'vector');
end
