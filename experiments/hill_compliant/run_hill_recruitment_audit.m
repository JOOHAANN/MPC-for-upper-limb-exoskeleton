function audit = run_hill_recruitment_audit(root)
    % RUN_HILL_RECRUITMENT_AUDIT Separate input evidence from inferred activity.
    % Raw project activation is not independently verified MVC-normalized EMG.
    % Counterfactuals change three input-mapping gains, holding mechanics fixed;
    % their results cannot establish a physiological need for cocontraction.
    cache = fullfile(root, 'results', 'hill_sequence', 'hill_sequence.mat');
    if ~isfile(cache)
        run_hill_sequence(root, struct('renderFigures', false));
    end
    data = load(cache, 'result');
    result = data.result;
    out = fullfile(root, 'results', 'hill_recruitment_audit');
    if ~isfolder(out)
        mkdir(out);
    end
    records = struct([]);
    counterfactual = struct([]);
    parameters = struct([]);
    sources = {'Experimental activation', 'MAP activation'};
    audit = struct('groups', {cell(1, 3)});
    for group = 1:3
        r = result.groups{group};
        t = r.sequence;
        moving = ~logical(t.waiting);
        details = cell(1, 2);
        for source = 1:2
            model = result.models{source};
            par = model.parameters;
            p = r.predictions{source};
            base = hill_pair_inputs(t, source, model.generator);
            d = memory(base(:, 1), .1, t.dt) - memory(base(:, 1), .35, t.dt);
            observed = par(15) * base(:, 2);
            inferred = par(16) * base(:, 1) + par(17) * max(0, d) + par(18) * max(0, -d);
            recruitment = observed + inferred;
            assert(max(abs(-expm1(-max(0, recruitment)) - r.inputs{source}(:, 2))) < 1e-10);
            posterior = p.backTorque(moving);
            abdominal = p.frontTorque(moving);
            rec = struct('group', group, 'source', sources{source}, ...
                         'backInputMean', mean(base(moving, 1)), 'frontBaseInputMean', mean(base(moving, 2)), ...
                         'frontBaseInputPeak', max(base(moving, 2)), ...
                         'backHillActivationMean', mean(p.effectiveActivation(moving, 1)), ...
                         'frontHillActivationMean', mean(p.effectiveActivation(moving, 2)), ...
                         'backHillActivationPeak', max(p.effectiveActivation(moving, 1)), ...
                         'frontHillActivationPeak', max(p.effectiveActivation(moving, 2)), ...
                         'backTorqueMeanNm', mean(posterior), 'frontTorqueMeanNm', mean(abdominal), ...
                         'backTorquePeakNm', max(posterior), 'frontTorquePeakNm', max(abdominal), ...
                         'frontToBackMeanTorqueRatio', sum(abdominal) / sum(posterior), ...
                         'frontGreaterThanBackFramePercent', 100 * mean(abdominal > posterior), ...
                         'frontRecruitmentFromBackPercent', 100 * sum(inferred(moving)) / sum(recruitment(moving)), ...
                         'backPositiveWorkJ', t.dt * sum(max(0, posterior .* p.v(moving))), ...
                         'frontPositiveWorkJ', t.dt * sum(max(0, -abdominal .* p.v(moving))), ...
                         'frontResistiveWorkDuringRiseJ', t.dt * sum(abdominal .* max(0, p.v(moving))));
            records = append(records, rec);
            qVariants = cell(1, 4);
            fractions = [1 .75 .5 0];
            failed = false(1, 4);
            for k = 1:numel(fractions)
                changed = par;
                changed(16:18) = fractions(k) * changed(16:18);
                if k == 1
                    prediction = p;
                else
                    inputs = hill_pair_inputs(t, source, model.generator, changed);
                    stance = struct('waiting', t.waiting, 'parameters', [t.standingAngle, result.options.controller]);
                    prediction = hill_pair_rollout(changed, inputs, t.m, [t.q(1), t.v0], t.dt, 10, true, stance);
                end
                failed(k) = prediction.failed;
                rec = struct('group', group, 'source', sources{source}, ...
                             'backDrivenFrontGainFraction', fractions(k), 'failed', double(prediction.failed), ...
                             'thetaRMSEDegrees', NaN, 'torqueRMSENm', NaN, ...
                             'minimumAngleRad', NaN, 'maximumAngleRad', NaN, ...
                             'softLimitPeakNm', NaN, 'frontTorqueMeanNm', NaN);
                if ~prediction.failed
                    qVariants{k} = prediction.q;
                    rec.thetaRMSEDegrees = sqrt(mean((prediction.q(moving) - t.q(moving)).^2)) * 180 / pi;
                    rec.torqueRMSENm = sqrt(mean((prediction.torque(moving) - t.torque(moving)).^2));
                    rec.minimumAngleRad = min(prediction.q);
                    rec.maximumAngleRad = max(prediction.q);
                    rec.softLimitPeakNm = max(abs(prediction.boundaryTorque));
                    rec.frontTorqueMeanNm = mean(prediction.frontTorque(moving));
                end
                counterfactual = append(counterfactual, rec);
            end
            details{source} = struct('base', base, 'observed', observed, 'inferred', inferred, ...
                                     'qVariants', {qVariants}, 'fractions', fractions, 'failed', failed);
            traces = table(t.t, t.waiting, base(:, 1), base(:, 2), ...
                           p.effectiveActivation(:, 1), p.effectiveActivation(:, 2), ...
                           p.backTorque, p.frontTorque, p.torque, observed, inferred, ...
                           'VariableNames', {'time', 'waiting', 'backBaseInput', 'frontBaseInput', ...
                                  'backHillActivation', 'frontHillActivation', ...
                                  'backTorqueNm', 'abdominalTorqueNm', 'netHillTorqueNm', ...
                                  'frontRecruitmentFromFrontInput', 'frontRecruitmentFromBackAndHistory'});
            writetable(traces, fullfile(out, sprintf('group%d_source%d.csv', group, source)));
            if group == 1
                rec = struct('source', sources{source}, 'backCapacityNm', par(1), 'frontCapacityNm', par(9), ...
                             'frontInputGain', par(15), 'backCoactivationGain', par(16), ...
                             'risingGain', par(17), 'fallingGain', par(18), ...
                             'backExponent', par(6), 'frontExponent', par(10), ...
                             'backMomentArmM', par(21), 'frontMomentArmM', par(22), ...
                             'derivedBackF0N', par(1) / par(21), 'derivedFrontF0N', par(9) / par(22));
                parameters = append(parameters, rec);
            end
        end
        audit.groups{group} = struct('baseline', r, 'details', {details});
        plot_audit(r, details, sources, out);
        plot_counterfactual(r, details, sources, out);
    end
    audit.metrics = struct2table(records);
    audit.parameters = struct2table(parameters);
    audit.counterfactuals = struct2table(counterfactual);
    writetable(audit.metrics, fullfile(out, 'recruitment_metrics.csv'));
    writetable(audit.parameters, fullfile(out, 'parameter_audit.csv'));
    writetable(audit.counterfactuals, fullfile(out, 'counterfactual_metrics.csv'));
    save(fullfile(out, 'recruitment_audit.mat'), 'audit');
    disp(audit.metrics(:, {'group', 'source', 'frontBaseInputMean', 'frontHillActivationMean', ...
                          'backTorqueMeanNm', 'frontTorqueMeanNm', 'frontRecruitmentFromBackPercent'}));
    fprintf('Audit completed: source evidence, model activation, opposing torques and fixed-mechanics counterfactuals.\n');
end

function plot_audit(r, details, sources, out)
    t = r.sequence;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1700 1400]);
    tiledlayout(4, 2, 'TileSpacing', 'compact');
    for row = 1:4
        for source = 1:2
            ax = nexttile;
            p = r.predictions{source};
            d = details{source};
            hold(ax, 'on');
            if row == 1
                plot(ax, t.t, d.base, 'LineWidth', 1);
                title(ax, [sources{source} ': base inputs']);
                ylabel(ax, 'Input activation');
                ylim(ax, [0 1]);
                legend(ax, 'Back input', 'Front base input', 'Location', 'northoutside');
            elseif row == 2
                plot(ax, t.t, p.effectiveActivation, 'LineWidth', 1);
                title(ax, 'Actual effective Hill activation (includes stance)');
                ylabel(ax, 'Model activation');
                ylim(ax, [0 1]);
                legend(ax, 'Back Hill activation', 'Abdominal Hill activation', 'Location', 'northoutside');
            elseif row == 3
                plot(ax, t.t, p.backTorque, 'LineWidth', 1);
                plot(ax, t.t, p.frontTorque, 'LineWidth', 1);
                plot(ax, t.t, p.torque, 'k--', 'LineWidth', 1);
                title(ax, 'Same axes: positive back / abdominal magnitudes and net torque');
                ylabel(ax, 'Torque (N m)');
                legend(ax, 'Back extensor', 'Abdominal flexor', 'Net = back - abdomen', 'Location', 'northoutside', 'FontSize', 9);
            else
                plot(ax, t.t, d.observed, 'LineWidth', 1);
                plot(ax, t.t, d.inferred, 'LineWidth', 1);
                title(ax, 'Contributions to the abdominal recruitment mapping');
                ylabel(ax, 'Recruitment signal');
                xlabel(ax, 'Time (s)');
                legend(ax, 'From front base input', 'From back input and history', 'Location', 'northoutside', 'FontSize', 9);
            end
            grid(ax, 'on');
            hill_sequence_protocol(ax, t, false);
        end
    end
    sgtitle(sprintf('Recruitment audit | Group %d | Project activation is not verified %%MVC | Model output is not a measurement', r.group));
    save_plot(fig, out, sprintf('recruitment_audit_group%d', r.group));
end

function plot_counterfactual(r, details, sources, out)
    t = r.sequence;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1550 750]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    for source = 1:2
        ax = nexttile;
        d = details{source};
        hold(ax, 'on');
        plot(ax, t.t, t.q, 'k--', 'LineWidth', 1, 'DisplayName', 'Original reference');
        for k = 1:numel(d.fractions)
            if ~d.failed(k)
                plot(ax, t.t, d.qVariants{k}, 'LineWidth', 1, ...
                     'DisplayName', sprintf('Back-driven front gain x %g', d.fractions(k)));
            end
        end
        title(ax, sources{source});
        ylabel(ax, 'Angle (rad)');
        xlabel(ax, 'Time (s)');
        legend(ax, 'Location', 'northoutside', 'FontSize', 9);
        grid(ax, 'on');
        hill_sequence_protocol(ax, t, false);
        if any(d.failed)
            text(ax, .01, .04, sprintf('Failed variants (not plotted): %s', num2str(d.fractions(d.failed))), ...
                 'Units', 'normalized', 'Color', [.7 .1 .1], 'FontSize', 10);
        end
    end
    sgtitle('Input-mapping ablation with fixed mechanics | Not a refit or evidence of physiological necessity');
    save_plot(fig, out, sprintf('mapping_ablation_group%d', r.group));
end

function save_plot(fig, out, name)
    exportgraphics(fig, fullfile(out, [name '.png']), 'Resolution', 170);
    exportgraphics(fig, fullfile(out, [name '.pdf']), 'ContentType', 'vector');
    savefig(fig, fullfile(out, [name '.fig']));
    close(fig);
end

function y = memory(x, tau, dt)
    r = exp(-dt / tau);
    y = filter(1 - r, [1 -r], x, r * x(1));
end

function records = append(records, rec)
    if isempty(records)
        records = rec;
    else
        records(end + 1) = rec;
    end
end
