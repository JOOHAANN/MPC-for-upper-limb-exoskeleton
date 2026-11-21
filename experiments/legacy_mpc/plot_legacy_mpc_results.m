function summary = plot_legacy_mpc_results(root, reports)
    % PLOT_LEGACY_MPC_RESULTS Show actual full replays without hiding overshoot.
    % The archived unassisted curves are stored baselines, not new simulations.
    out = fullfile(root, 'results', 'legacy_mpc');
    if nargin < 2
        sources = {'control_experimental', 'control_map'};
        reports = struct([]);
        for i = 1:2
            reports = [reports; jsondecode(fileread(fullfile(out, [sources{i} '_checked_status.json'])))]; %#ok<AGROW>
        end
    end
    assert(numel(reports) == 2 && all([reports.completed]), ...
           'Both complete replays are required; a partial tail must not be plotted.');
    labels = {'Experimental activation', 'MAP activation'};
    records = struct([]);
    main = figure('Name', 'Archived MPC: torque, theta and assistance', 'Color', 'w', 'Position', [80 40 1180 950]);
    primary = tiledlayout(main, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    components = figure('Name', 'Archived MPC: external and abdominal terms', 'Color', 'w', 'Position', [100 50 1180 900]);
    secondary = tiledlayout(components, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    blue = [.09 .4 .68]; red = [.83 .12 .16]; grey = [.52 .52 .52];
    for source = 1:2
        d = load(fullfile(reports(source).workDirectory, 'legacy_execution.mat'));
        count = d.NumberOfData;
        assert(reports(source).lastStep == count - 1);
        t = (0:count - 1) * d.Ts;
        back = d.MTF(1:count) * d.coeffs(1);
        assist = d.Tac(1:count);
        used = isfinite(d.legacyExitFlag);
        rec = struct('source', labels{source}, 'minimumAngleRad', min(d.theta), ...
                     'maximumAngleRad', max(d.theta), 'assistancePeakNm', max(assist), ...
                     'softBoundaryAbsPeakNm', max(abs(d.F_k_up)), ...
                     'externalBalanceAbsPeakNm', max(abs(d.T_balance)), ...
                     'solverCalls', nnz(used), 'nonpositiveSolverExits', nnz(d.legacyExitFlag(used) <= 0));
        records = [records; rec]; %#ok<AGROW>
        ax = nexttile(primary, source); hold(ax, 'on');
        plot(ax, t, d.Torque_Ceinms, '--', 'Color', red, 'LineWidth', 1.1);
        plot(ax, t, d.Torque_real, 'Color', grey, 'LineWidth', 1.0);
        plot(ax, t, back, 'Color', blue, 'LineWidth', 1.4);
        plot(ax, t, d.Threshold, ':k', 'LineWidth', 1.0);
        ylabel(ax, 'Legacy back torque (N m)'); title(ax, labels{source});
        legend(ax, {'Reference', 'Archived unassisted', 'Legacy MPC assisted', 'Threshold'}, ...
               'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8, 'Box', 'off');
        style_axes(ax, t);
        ax = nexttile(primary, source + 2); hold(ax, 'on');
        plot(ax, t, d.theta_Ceinms, '--', 'Color', red, 'LineWidth', 1.1);
        plot(ax, t, d.theta_real, 'Color', grey);
        plot(ax, t, d.theta, 'Color', blue, 'LineWidth', 1.4);
        ylabel(ax, '\theta (rad)');
        title(ax, sprintf('Actual range: %.2f to %.2f rad', min(d.theta), max(d.theta)));
        style_axes(ax, t);
        ax = nexttile(primary, source + 4); hold(ax, 'on');
        plot(ax, t, assist, 'Color', [.15 .52 .3], 'LineWidth', 1.4);
        plot(ax, t, d.u_max, '--', 'Color', red);
        ylabel(ax, 'Applied assistance (N m)'); xlabel(ax, 'Time (s)');
        title(ax, sprintf('Peak: %.2f N m; intended cap: %.0f N m', max(assist), max(d.u_max)));
        style_axes(ax, t);
        ax = nexttile(secondary, source); hold(ax, 'on');
        plot(ax, t, d.F_k_up, 'Color', [.55 .28 .65]);
        plot(ax, t, d.T_balance, 'Color', [.83 .52 .08]);
        title(ax, labels{source}); ylabel(ax, 'External torque (N m)');
        legend(ax, {'Old soft boundary', 'External low-activation balance'}, 'Location', 'southoutside', 'Box', 'off');
        style_axes(ax, t);
        ax = nexttile(secondary, source + 2);
        plot(ax, t, d.T_front, 'Color', red); yline(ax, 0, ':k');
        ylabel(ax, 'Old abdominal term (N m)');
        title(ax, 'Signed angle polynomial in the old model');
        style_axes(ax, t);
        ax = nexttile(secondary, source + 4);
        scatter(ax, t(used), d.legacyExitFlag(used), 4, blue, 'filled');
        yline(ax, 0, '--', 'Color', red); yticks(ax, [-2 0 1 2]);
        ylabel(ax, 'fmincon exit flag'); xlabel(ax, 'Time (s)');
        style_axes(ax, t);
    end
    summary = struct2table(records);
    writetable(summary, fullfile(out, 'live_summary.csv'));
    title(primary, 'Archived BSE and MPC settings; separate from the selected neural/Hill model');
    title(secondary, 'Legacy torque decomposition and solver diagnostics');
    exportgraphics(main, fullfile(out, 'legacy_mpc_torque_theta.png'), 'Resolution', 170);
    exportgraphics(components, fullfile(out, 'legacy_mpc_components.png'), 'Resolution', 170);
    disp(summary);
end

function style_axes(ax, t)
    set(ax, 'FontName', 'Arial', 'FontSize', 10);
    xlim(ax, [t(1) t(end)]); grid(ax, 'on'); box(ax, 'on');
end
