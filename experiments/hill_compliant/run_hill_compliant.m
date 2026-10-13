function result = run_hill_compliant(root, cache, mode, renderFigures)
    % RUN_HILL_COMPLIANT Replay transparent, coupled muscle-driven calibration.
    if nargin < 2
        cache = fullfile(root, 'data', 'processed', 'hill_compliant', 'calibrated_models.mat');
    end
    assert(isfile(cache), 'Missing calibrated Hill package.');
    if nargin < 3
        mode = 'Calibration';
    end
    if nargin < 4
        renderFigures = true;
    end
    d = load(cache);
    models = d.models;
    trials = d.trials;
    heldSpeeds = [];
    if iscell(trials)
        trials = [trials{:}];
    end
    out = fullfile(root, 'results', 'hill_compliant');
    if strcmp(mode, 'Whole-speed holdout')
        models = d.validationModels;
        heldSpeeds = d.heldSpeeds;
        out = fullfile(out, 'validation');
    end
    if ~isfolder(out)
        mkdir(out);
    end
    labels = {'Experimental activation', 'MAP activation'};
    predictions = cell(2, numel(trials));
    inputs = predictions;
    records = struct([]);
    for source = 1:2
        for k = 1:numel(trials)
            t = trials(k);
            model = select_model(models, heldSpeeds, source, t.speed);
            if ~isempty(heldSpeeds)
                assert(~ismember(k, model.trainIds) && ismember(k, model.testIds), 'Training/test overlap.');
            end
            a = hill_pair_inputs(t, source, model.generator, model.parameters);
            p = hill_pair_rollout(model.parameters, a, t.m, [t.q(1) t.v0], .01, 10);
            assert(~p.failed, 'Hill rollout failed; no angle clipping or fallback is applied.');
            predictions{source, k} = p;
            inputs{source, k} = a;
            rec = score(t, p);
            rec.source = labels{source};
            rec.loadKg = t.loadKg;
            rec.speed = t.speed;
            if isempty(records)
                records = rec;
            else
                records(end + 1) = rec;
            end %#ok<AGROW>
            sample = table(t.t, t.q, p.q, t.torque, p.torque, p.backTorque, p.frontTorque, p.backPassive, p.frontPassive, ...
                           p.gravityTorque, p.dampingTorque, p.boundaryTorque, p.equilibriumResidual, p.momentumResidual, ...
                           p.backActive, p.frontActive, p.backFiberDamping, p.frontFiberDamping, p.effectiveActivation(:, 1), p.effectiveActivation(:, 2), ...
                           'VariableNames', {'time', 'referenceAngle', 'predictedAngle', 'referenceTorque', 'netHillTorque', ...
                                           'backHillTorque', 'frontHillTorque', 'backPassiveTorque', 'frontPassiveTorque', 'gravityTorque', ...
                                           'jointDampingTorque', 'boundaryTorque', 'forceEquilibriumResidual', 'momentumResidual', ...
                                           'backActiveTorque', 'frontActiveTorque', 'backFiberDampingTorque', 'frontFiberDampingTorque', 'backHillActivation', 'frontHillActivation'});
            writetable(sample, fullfile(out, sprintf('source%d_load%d_speed%d.csv', source, t.loadKg, t.speed)));
        end
    end
    result = struct('models', {models}, 'trials', trials, 'heldSpeeds', heldSpeeds, 'predictions', {predictions}, 'inputs', {inputs}, 'metrics', struct2table(records));
    old = fullfile(root, 'data', 'processed', 'hill_compliant', 'archived_bse.mat');
    if isfile(old)
        d = load(old, 'archived');
        result.archived = d.archived;
    end
    save(fullfile(out, 'hill_compliant.mat'), 'result');
    writetable(result.metrics, fullfile(out, 'metrics.csv'));
    if renderFigures
        plot_hill_compliant(result, out, mode);
    end
    fprintf('\n%s\n', mode);
    disp(groupsummary(result.metrics, 'source', 'mean', {'thetaRMSE', 'torqueRMSE', 'correlation', 'Wmatched', 'outsidePercent'}));
    if strcmp(mode, 'Calibration')
        result.validation = run_hill_compliant(root, cache, 'Whole-speed holdout', renderFigures);
        save(fullfile(out, 'hill_compliant.mat'), 'result');
    end
end

function model = select_model(models, heldSpeeds, source, speed)
    if isempty(heldSpeeds)
        model = models{source};
    else
        model = models{source, find(heldSpeeds == speed, 1)};
    end
end

function r = score(t, p)
    [~, ref] = findpeaks(-t.q, 'MinPeakProminence', .12, 'MinPeakDistance', 50);
    [~, pred] = findpeaks(-p.q, 'MinPeakProminence', .12, 'MinPeakDistance', 50);
    c = corrcoef(t.q, p.q);
    match = numel(ref) == 2 && numel(pred) == 2;
    timing = NaN;
    valley = NaN;
    if match
        timing = mean(abs(ref - pred)) * .01;
        valley = mean(abs(t.q(ref) - p.q(pred)));
    end
    r = struct('thetaRMSE', sqrt(mean((p.q - t.q).^2)), 'torqueRMSE', sqrt(mean((p.torque - t.torque).^2)), ...
               'correlation', c(1, 2), 'Wmatched', double(match), 'predictedValleys', numel(pred), ...
               'minimumError', abs(min(p.q) - min(t.q)), 'maximumError', abs(max(p.q) - max(t.q)), ...
               'valleyTimeMAE', timing, 'valleyAngleMAE', valley, 'outsidePercent', 100 * mean(p.q < 0 | p.q > 1.8), ...
               'minAngle', min(p.q), 'maxAngle', max(p.q), 'maxAbsTorque', max(abs(p.torque)), ...
               'backPassiveMax', max(p.backPassive), 'frontPassiveMax', max(p.frontPassive), ...
               'boundaryMax', max(abs(p.boundaryTorque)), 'boundaryRMS', sqrt(mean(p.boundaryTorque.^2)), ...
               'gravityMax', max(abs(p.gravityTorque)), 'dampingMax', max(abs(p.dampingTorque)), ...
               'torqueCorrelation', corr(p.torque, t.torque), 'equilibriumResidualMax', max(p.equilibriumResidual), ...
               'momentumResidualMax', max(abs(p.momentumResidual)));
end
