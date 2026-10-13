function audit = audit_bse_targets(root, outputDir)
%AUDIT_BSE_TARGETS Separate synthetic calibration targets from recorded torque.
% Runs the archived general_model_op_scale BSE exactly as its check script.
% Its saved coefficients were calibrated on these cases: this is an apparent
% fit diagnostic, never an independent generalization result.
    config = jsondecode(fileread(fullfile(root, 'config', 'workflows', ...
        'training_map__muscle_tendon__general_model_op_scale.json')));
    filenames = {'Combination_5_35_10_55_15_45.mat', ...
        'Combination_5_45_10_35_15_55.mat', 'Combination_5_55_10_45_15_35.mat'};
    arrays = cell(3, 1);
    for k = 1:3
        data = load(asset(root, config, filenames{k}), 'Combination');
        arrays{k} = data.Combination;
    end
    combined = vertcat(arrays{:});
    d = load(asset(root, config, 'coeffs.mat'), 'coeffs');
    c = num2cell(d.coeffs);
    assert(numel(c) == 12);
    source = asset(root, config, 'Torque_function_constant.m');
    previousPath = path;
    cleanup = onCleanup(@() path(previousPath));
    addpath(fileparts(source), '-begin');
    clear Torque_function_constant;
    assert(strcmp(which('Torque_function_constant'), source), 'Wrong BSE variant on path.');
    x = Torque_function_constant(zeros(size(combined, 1), 1), combined(:, 3), ...
        combined(:, 3), combined(:, 2), combined(:, 4), 0.5, 0.24086, 0.001, c{:});
    windows = [5321 11740; 22037 26059; 36732 41109; ...
        52761 58381; 66452 72417; 84099 87884; ...
        98004 102357; 113394 118437; 127817 133717];
    conditions = [5 35; 10 55; 15 45; 5 45; 10 35; 15 55; 5 55; 10 45; 15 35];
    audit = struct([]);
    fig = figure('Color', 'w', 'Position', [60 60 1400 850]);
    tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:9
        idx = windows(k, 1)+1:windows(k, 2)-1;
        synthetic = combined(idx, 1);
        recorded = combined(idx, 4);
        predicted = x(idx, 1);
        audit(k).loadKg = conditions(k, 1);
        audit(k).speedLabel = conditions(k, 2);
        audit(k).syntheticTargetIdentityMaxError = max(abs(synthetic-310*combined(idx, 2)));
        audit(k).syntheticVsRecordedRMSE = sqrt(mean((synthetic-recorded).^2));
        audit(k).bseVsSyntheticRMSE = sqrt(mean((predicted-synthetic).^2));
        audit(k).bseVsRecordedRMSE = sqrt(mean((predicted-recorded).^2));
        audit(k).angleRMSE = sqrt(mean((x(idx, 2)-combined(idx, 3)).^2));
        nexttile; hold on;
        t = (0:numel(idx)-1)*0.001;
        plot(t, recorded, 'k', 'LineWidth', 1.4);
        plot(t, synthetic, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
        plot(t, predicted, 'Color', [0.12 0.42 0.72], 'LineWidth', 1.2);
        title(sprintf('%g kg / %g | ref RMSE %.1f Nm', conditions(k, :), audit(k).bseVsRecordedRMSE));
        xlabel('Time in segment (s)'); ylabel('Torque (Nm)'); grid on;
        if k == 1, legend('Recorded reference (col 4)', '310 x activation (col 1)', 'Archived BSE', 'Location', 'best'); end
    end
    sgtitle('Current BSE audit: calibration uses synthetic torque | archived fit, NOT held-out validation');
    writetable(struct2table(audit), fullfile(outputDir, 'original_bse_audit.csv'));
    exportgraphics(fig, fullfile(outputDir, 'original_bse_target_audit.png'), 'Resolution', 160);
    exportgraphics(fig, fullfile(outputDir, 'original_bse_target_audit.pdf'), 'ContentType', 'vector');
    save(fullfile(outputDir, 'original_bse_audit.mat'), 'audit', 'windows', 'conditions');
end

function filename = asset(root, config, name)
    index = find(strcmp({config.files.target}, name));
    assert(isscalar(index), 'Missing or ambiguous manifest entry: %s', name);
    filename = fullfile(root, config.files(index).source);
end
