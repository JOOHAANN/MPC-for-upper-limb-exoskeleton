function plot_bse_results(results, outputDir)
%PLOT_BSE_RESULTS Plot all outer folds and fixed, preselected example trials.
    colors = [0.55 0.58 0.63; 0.12 0.42 0.72; 0.89 0.35 0.15];
    protocols = {'load', 'speed'};
    metrics = results.metrics;
    fig = figure('Color', 'w', 'Position', [80 80 1250 780]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:2
        nexttile(p);
        means = zeros(1, 3);
        for m = 1:3
            mask = strcmp(metrics.protocol, protocols{p}) & metrics.model == m;
            means(m) = mean(metrics.rmse(mask));
        end
        b = bar(means, 'FaceColor', 'flat'); b.CData = colors;
        xticklabels({'Lumped', 'Physics', 'Hybrid'});
        ylabel('Mean trial RMSE (Nm)'); grid on;
        title(sprintf('Leave-%s-out: all 25 trials', protocols{p}));
        for m = 1:3
            text(m, means(m), sprintf(' %.2f', means(m)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
        ylim([0 max(means)*1.2]);
        nexttile(p+2); hold on;
        fm = results.foldMetrics;
        for m = 1:3
            mask = strcmp(fm.protocol, protocols{p}) & fm.model == m;
            plot(fm.heldOut(mask), fm.macroRMSE(mask), '-o', ...
                'Color', colors(m, :), 'LineWidth', 1.5);
        end
        ylabel('Mean held-out trial RMSE (Nm)');
        if p == 1, xlabel('Held-out load (kg)'); else, xlabel('Held-out speed label'); end
        grid on; legend(results.modelNames, 'Location', 'best');
        title('Each fold retrains and tunes on training groups only');
    end
    sgtitle('Offline torque estimation | nested grouped validation | no random sample split');
    save_plot(fig, outputDir, 'generalization_summary');

    fig = figure('Color', 'w', 'Position', [90 90 1250 780]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    examples = [12.5 45 1; 15 45 1; 10 65 2; 5 25 2];
    counts = arrayfun(@(t) numel(t.torque), results.trials);
    offsets = [0, cumsum(counts)];
    for k = 1:size(examples, 1)
        trial = find([results.trials.loadKg] == examples(k, 1) & ...
            [results.trials.speedLabel] == examples(k, 2));
        t = results.trials(trial);
        indices = offsets(trial)+1:offsets(trial+1);
        pred = results.predictions{examples(k, 3)}(indices, :);
        nexttile; hold on;
        plot(t.time, t.torque, 'k', 'LineWidth', 1.6);
        for m = 1:3
            plot(t.time, pred(:, m), 'Color', colors(m, :), 'LineWidth', 1.1);
        end
        xlabel('Time in recording (s)'); ylabel('Reference torque (Nm)'); grid on;
        title(sprintf('%g kg / speed %g | leave-%s-out', t.loadKg, t.speedLabel, protocols{examples(k, 3)}));
        if k == 1, legend([{'Torque.sto'}, results.modelNames], 'Location', 'best'); end
    end
    sgtitle('Preselected held-out examples | measured kinematics and EMG are inputs');
    save_plot(fig, outputDir, 'held_out_predictions');

    fig = figure('Color', 'w', 'Position', [100 100 1200 520]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    for p = 1:2
        delta = nan(5);
        loads = unique([results.trials.loadKg]); speeds = unique([results.trials.speedLabel]);
        for i = 1:5
            for j = 1:5
                mask = strcmp(metrics.protocol, protocols{p}) & ...
                    metrics.loadKg == loads(i) & metrics.speedLabel == speeds(j);
                base = metrics.rmse(mask & metrics.model == 2);
                hybrid = metrics.rmse(mask & metrics.model == 3);
                delta(i, j) = 100*(base-hybrid)/base;
            end
        end
        nexttile; imagesc(speeds, loads, delta); set(gca, 'YDir', 'normal');
        xlabel('Speed label'); ylabel('Load (kg)'); colorbar; clim([-50 50]);
        title(sprintf('Leave-%s-out: hybrid RMSE reduction (%%)', protocols{p}));
        for i = 1:5
            for j = 1:5
                text(speeds(j), loads(i), sprintf('%+.0f', delta(i, j)), ...
                    'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
            end
        end
    end
    sgtitle('All conditions shown | positive = better than physics basis; negative = worse');
    save_plot(fig, outputDir, 'condition_improvement');
end

function save_plot(fig, folder, name)
    exportgraphics(fig, fullfile(folder, [name '.png']), 'Resolution', 160);
    exportgraphics(fig, fullfile(folder, [name '.pdf']), 'ContentType', 'vector');
end
