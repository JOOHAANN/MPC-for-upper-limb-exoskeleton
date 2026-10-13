function plot_hill_compliant(r, out, mode)
    % PLOT_HILL_COMPLIANT Display all conditions and expose every force component.
    for speed = [35 45 55]
        indices = find([r.trials.speed] == speed);
        [~, order] = sort([r.trials(indices).loadKg]);
        indices = indices(order);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [25 25 1800 1100]);
        tiledlayout(3, 4, 'TileSpacing', 'compact');
        for k = indices
            t = r.trials(k);
            p = r.predictions{1, k};
            g = r.predictions{2, k};
            generatedModel = select_model(r.models, r.heldSpeeds, 2, t.speed);
            rawGenerated = hill_pair_inputs(t, 2, generatedModel.generator);
            nexttile;
            plot(t.t, t.a, 'LineWidth', 1.2);
            hold on;
            plot(t.t, [t.front, rawGenerated(:, 2)], '--');
            grid on;
            title(sprintf('%g kg / speed %g', t.loadKg, speed));
            ylabel('Activation');
            legend('Experimental back', 'MAP back', 'Experimental front (filtered)', 'Auxiliary generated front', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
            nexttile;
            colororder(gca, [0 .35 .7; .85 .25 .1; .15 .6 .3]);
            plot(t.t, [t.q, p.q, g.q] * 180 / pi, 'LineWidth', 1.3);
            hold on;
            names = {'Original Scale reference', 'Hill: experimental', 'Hill: MAP'};
            if isfield(r, 'archived')
                plot(t.t, r.archived{k}(:, 2) * 180 / pi, '--', 'Color', [.55 .55 .55]);
                names{4} = 'Selected archived BSE';
            end
            ylabel('Torso angle (deg)');
            grid on;
            legend(names, 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
            nexttile;
            colororder(gca, [0 .35 .7; .85 .25 .1; .15 .6 .3]);
            plot(t.t, [t.torque, p.torque, g.torque], 'LineWidth', 1.2);
            hold on;
            torqueNames = {'Reference', 'Hill: experimental', 'Hill: MAP'};
            if isfield(r, 'archived')
                plot(t.t, r.archived{k}(:, 1), '--', 'Color', [.55 .55 .55]);
                torqueNames{4} = 'Selected archived BSE';
            end
            ylabel('Net Hill torque (Nm)');
            grid on;
            legend(torqueNames, 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
            nexttile;
            colororder(gca, [0 .35 .7; .85 .25 .1; .15 .6 .3; .55 .2 .65]);
            plot(t.t, [p.backPassive, -p.frontPassive, g.backPassive, -g.frontPassive], 'LineWidth', 1.1);
            hold on;
            plot(t.t, p.boundaryTorque, 'k--', 'LineWidth', 1.6);
            plot(t.t, g.boundaryTorque, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.3);
            grid on;
            ylabel('Passive / boundary torque (Nm)');
            legend('Exp back PE', 'Exp front PE', 'MAP back PE', 'MAP front PE', 'Exp soft limit', 'MAP soft limit', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
        end
        sgtitle([mode ' | Same compliant Hill pair | Force-driven coupled integration']);
        exportgraphics(fig, fullfile(out, sprintf('comparison_speed%d.png', speed)), 'Resolution', 160);
        exportgraphics(fig, fullfile(out, sprintf('comparison_speed%d.pdf', speed)), 'ContentType', 'vector');
        savefig(fig, fullfile(out, sprintf('comparison_speed%d.fig', speed)));
        close(fig);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [30 30 1400 1000]);
        tiledlayout(3, 2, 'TileSpacing', 'compact');
        for k = indices
            for source = 1:2
                t = r.trials(k);
                p = r.predictions{source, k};
                nexttile;
                plot(t.t, [p.torque, p.backTorque, -p.frontTorque], 'LineWidth', 1.2);
                hold on;
                plot(t.t, t.torque, 'k--');
                grid on;
                title(sprintf('%g kg / source %d', t.loadKg, source));
                ylabel('Torque (Nm)');
                legend('Net Hill', 'Extensor Hill', 'Flexor Hill (signed)', 'Reference', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
            end
        end
        sgtitle([mode ' | Opposing Hill torques; soft-limit force excluded from net Hill output']);
        exportgraphics(fig, fullfile(out, sprintf('torque_components_speed%d.png', speed)), 'Resolution', 160);
        close(fig);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 1300 950]);
        tiledlayout(3, 2, 'TileSpacing', 'compact');
        for k = indices
            for source = 1:2
                t = r.trials(k);
                p = r.predictions{source, k};
                model = select_model(r.models, r.heldSpeeds, source, t.speed);
                raw = hill_pair_inputs(t, source, model.generator);
                nexttile;
                plot(t.t, [raw, p.effectiveActivation], 'LineWidth', 1.1);
                grid on;
                title(sprintf('%g kg / source %d', t.loadKg, source));
                ylabel('Activation');
                legend('Back input', 'Front input', 'Back entering Hill', 'Estimated front entering Hill', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
            end
        end
        sgtitle([mode ' | Inputs versus calibrated Hill activation; recruitment is estimated']);
        exportgraphics(fig, fullfile(out, sprintf('activation_calibration_speed%d.png', speed)), 'Resolution', 160);
        close(fig);
        % Show exactly when the soft spring acts and its actual force budget.
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 1400 1000]);
        tiledlayout(3, 2, 'TileSpacing', 'compact');
        for k = indices
            t = r.trials(k);
            p = r.predictions{1, k};
            g = r.predictions{2, k};
            model = select_model(r.models, r.heldSpeeds, 1, t.speed);
            par = model.parameters;
            nexttile;
            plot(t.t, [t.q, p.q, g.q], 'LineWidth', 1.3);
            hold on;
            yline(par(31), 'k--', 'Lower soft limit', 'HandleVisibility', 'off');
            yline(par(32), 'k--', 'Upper soft limit', 'HandleVisibility', 'off');
            grid on;
            ylabel('Angle (rad)');
            title(sprintf('%g kg | limits %.2f / %.2f rad', t.loadKg, par(31), par(32)));
            legend('Original Scale reference', 'Experimental input', 'MAP input', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
            nexttile;
            plot(t.t, [p.boundaryTorque, g.boundaryTorque], 'LineWidth', 1.4);
            hold on;
            yline(0, 'k:', 'HandleVisibility', 'off');
            grid on;
            ylabel('Soft-limit torque (Nm)');
            title(sprintf('Peak exp %.1f / MAP %.1f Nm | k %.0f, cap %.0f', max(abs(p.boundaryTorque)), max(abs(g.boundaryTorque)), par(29), par(30)));
            legend('Experimental input', 'MAP input', 'Location', 'northoutside', 'Orientation', 'horizontal');
        end
        sgtitle([mode ' | Elastic soft limits: no reference tracking, no boundary damping']);
        exportgraphics(fig, fullfile(out, sprintf('soft_limits_speed%d.png', speed)), 'Resolution', 160);
        exportgraphics(fig, fullfile(out, sprintf('soft_limits_speed%d.pdf', speed)), 'ContentType', 'vector');
        close(fig);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [45 45 1400 1000]);
        tiledlayout(3, 2, 'TileSpacing', 'compact');
        for k = indices
            for source = 1:2
                t = r.trials(k);
                p = r.predictions{source, k};
                nexttile;
                plot(t.t, [p.torque, p.gravityTorque, p.dampingTorque, p.boundaryTorque], 'LineWidth', 1.2);
                hold on;
                yline(0, 'k:', 'HandleVisibility', 'off');
                grid on;
                title(sprintf('%g kg / source %d', t.loadKg, source));
                ylabel('Signed torque (Nm)');
                legend('Net Hill', 'Gravity', 'Joint damping', 'Soft limit', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 8);
            end
        end
        sgtitle([mode ' | Equation budget: I * acceleration = Hill + gravity + damping + soft limit']);
        exportgraphics(fig, fullfile(out, sprintf('equation_budget_speed%d.png', speed)), 'Resolution', 160);
        close(fig);
        % Show the two main valleys at a useful scale; the full archived
        % angle remains visible in the separate four-column comparison.
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [45 45 1400 950]);
        tiledlayout(3, 2, 'TileSpacing', 'compact');
        for k = indices
            for source = 1:2
                t = r.trials(k);
                p = r.predictions{source, k};
                nexttile;
                plot(t.t, [t.q, p.q] * 180 / pi, 'LineWidth', 1.4);
                hold on;
                [~, ref] = findpeaks(-t.q, 'MinPeakProminence', .12, 'MinPeakDistance', 50);
                [~, pred] = findpeaks(-p.q, 'MinPeakProminence', .12, 'MinPeakDistance', 50);
                plot(t.t(ref), t.q(ref) * 180 / pi, 'o', 'Color', [0 .35 .7], 'HandleVisibility', 'off');
                plot(t.t(pred), p.q(pred) * 180 / pi, 's', 'Color', [.85 .25 .1], 'HandleVisibility', 'off');
                grid on;
                ylim([min(0, min(p.q)), max(1.8, max(p.q))] * 180 / pi);
                xlabel('Time (s)');
                ylabel('Angle (deg)');
                title(sprintf('%g kg | source %d | RMSE %.1f deg | valleys %d', t.loadKg, source, sqrt(mean((p.q - t.q).^2)) * 180 / pi, numel(pred)));
                legend('Original Scale reference', 'Hill prediction', 'Location', 'northoutside', 'Orientation', 'horizontal');
            end
        end
        sgtitle([mode ' | Angle trends and valley locations | Source 1: experimental, source 2: MAP']);
        exportgraphics(fig, fullfile(out, sprintf('theta_focus_speed%d.png', speed)), 'Resolution', 160);
        exportgraphics(fig, fullfile(out, sprintf('theta_focus_speed%d.pdf', speed)), 'ContentType', 'vector');
        close(fig);
    end
end

function model = select_model(models, heldSpeeds, source, speed)
    if isempty(heldSpeeds)
        model = models{source};
    else
        model = models{source, find(heldSpeeds == speed, 1)};
    end
end
