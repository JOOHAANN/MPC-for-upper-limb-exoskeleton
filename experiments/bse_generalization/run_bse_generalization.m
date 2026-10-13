function results = run_bse_generalization(root, rebuildData)
%RUN_BSE_GENERALIZATION Nested leave-load/leave-speed-out torque benchmark.
% Hyperparameters and all feature scaling are selected on outer-training data.
% No saved full-data BSE coefficients enter the cross-validation benchmark.
    if nargin < 2
        rebuildData = false;
    end
    setup_project();
    outputDir = fullfile(root, 'results', 'bse_generalization');
    if ~isfolder(outputDir), mkdir(outputDir); end
    cacheDir = fullfile(root, 'data', 'processed', 'bse_generalization');
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    cache = fullfile(cacheDir, 'raw_sensor_trials.mat');
    if rebuildData || ~isfile(cache)
        trials = build_bse_trials(root);
        preprocessingVersion = 1;
        save(cache, 'trials', 'preprocessingVersion');
    else
        data = load(cache, 'trials', 'preprocessingVersion');
        assert(data.preprocessingVersion == 1, 'Rebuild outdated data cache.');
        trials = data.trials;
    end
    assert(numel(trials) == 25);
    [lumped, physics, residual, y, group] = stack_trials(trials);
    lambdas = [0.001, 0.01, 0.1, 1, 10];
    protocols = {'load', 'speed'};
    rows = struct([]);
    foldRows = struct([]);
    predictions = cell(2, 1);
    selected = struct([]);
    for protocol = 1:2
        if protocol == 1
            levels = [trials.loadKg];
        else
            levels = [trials.speedLabel];
        end
        uniqueLevels = unique(levels);
        predicted = nan(numel(y), 3);
        for fold = 1:numel(uniqueLevels)
            testTrials = find(levels == uniqueLevels(fold));
            trainTrials = find(levels ~= uniqueLevels(fold));
            train = ismember(group, trainTrials);
            test = ismember(group, testTrials);
            assert(~any(train & test));
            innerLevels = unique(levels(trainTrials));
            score = zeros(numel(lambdas), 2);
            % Inner folds also leave out entire loads or speeds.
            for inner = 1:numel(innerLevels)
                validation = train & ismember(group, find(levels == innerLevels(inner)));
                fitting = train & ~validation;
                for l = 1:numel(lambdas)
                    for model = 1:2
                        if model == 1, x = lumped; else, x = physics; end
                        fitted = bse_ridge('fit', x(fitting, :), y(fitting), group(fitting), lambdas(l));
                        estimate = bse_ridge('predict', fitted, x(validation, :));
                        score(l, model) = score(l, model) + macro_mse(y(validation), estimate, group(validation));
                    end
                end
            end
            [~, best] = min(score, [], 1);
            lambdaLumped = lambdas(best(1));
            lambdaPhysics = lambdas(best(2));
            residualScore = zeros(numel(lambdas)+1, 1);
            for inner = 1:numel(innerLevels)
                validation = train & ismember(group, find(levels == innerLevels(inner)));
                fitting = train & ~validation;
                base = bse_ridge('fit', physics(fitting, :), y(fitting), group(fitting), lambdaPhysics);
                errorTrain = y(fitting) - bse_ridge('predict', base, physics(fitting, :));
                baseValidation = bse_ridge('predict', base, physics(validation, :));
                residualScore(end) = residualScore(end) + macro_mse(y(validation), baseValidation, group(validation));
                for l = 1:numel(lambdas)
                    correction = bse_ridge('fit', residual(fitting, :), errorTrain, group(fitting), lambdas(l));
                    estimate = baseValidation + bse_ridge('predict', correction, residual(validation, :));
                    residualScore(l) = residualScore(l) + macro_mse(y(validation), estimate, group(validation));
                end
            end
            [~, bestResidual] = min(residualScore);
            fitLumped = bse_ridge('fit', lumped(train, :), y(train), group(train), lambdaLumped);
            fitPhysics = bse_ridge('fit', physics(train, :), y(train), group(train), lambdaPhysics);
            predicted(test, 1) = bse_ridge('predict', fitLumped, lumped(test, :));
            predicted(test, 2) = bse_ridge('predict', fitPhysics, physics(test, :));
            predicted(test, 3) = predicted(test, 2);
            lambdaResidual = Inf;
            if bestResidual <= numel(lambdas)
                lambdaResidual = lambdas(bestResidual);
                errorTrain = y(train) - bse_ridge('predict', fitPhysics, physics(train, :));
                correction = bse_ridge('fit', residual(train, :), errorTrain, group(train), lambdaResidual);
                predicted(test, 3) = predicted(test, 3) + bse_ridge('predict', correction, residual(test, :));
            end
            k = numel(selected) + 1;
            selected(k).protocol = protocols{protocol};
            selected(k).heldOut = uniqueLevels(fold);
            selected(k).trainTrials = trainTrials;
            selected(k).testTrials = testTrials;
            selected(k).lambdaLumped = lambdaLumped;
            selected(k).lambdaPhysics = lambdaPhysics;
            selected(k).lambdaResidual = lambdaResidual;
            for model = 1:3
                k = numel(foldRows)+1;
                foldRows(k).protocol = protocols{protocol};
                foldRows(k).heldOut = uniqueLevels(fold);
                foldRows(k).model = model;
                foldRows(k).macroRMSE = mean_trial_rmse(y(test), predicted(test, model), group(test));
            end
            fprintf('%s held out %g: selected residual lambda %g\n', protocols{protocol}, uniqueLevels(fold), lambdaResidual);
        end
        assert(all(isfinite(predicted), 'all'));
        predictions{protocol} = predicted;
        for trial = 1:numel(trials)
            mask = group == trial;
            for model = 1:3
                k = numel(rows)+1;
                rows(k).protocol = protocols{protocol};
                rows(k).trial = trials(trial).id;
                rows(k).loadKg = trials(trial).loadKg;
                rows(k).speedLabel = trials(trial).speedLabel;
                rows(k).model = model;
                rows(k).rmse = sqrt(mean((predicted(mask, model)-y(mask)).^2));
                rows(k).mae = mean(abs(predicted(mask, model)-y(mask)));
                rows(k).rSquared = 1-sum((predicted(mask, model)-y(mask)).^2)/sum((y(mask)-mean(y(mask))).^2);
            end
        end
    end
    metrics = struct2table(rows);
    foldMetrics = struct2table(foldRows);
    modelNames = {'Lumped activation', 'Physics basis', 'Physics + residual'};
    writetable(metrics, fullfile(outputDir, 'trial_metrics.csv'));
    writetable(foldMetrics, fullfile(outputDir, 'fold_metrics.csv'));
    results = struct('trials', trials, 'metrics', metrics, 'foldMetrics', foldMetrics, ...
        'predictions', {predictions}, 'selected', selected, 'modelNames', {modelNames});
    save(fullfile(outputDir, 'benchmark_results.mat'), 'results');
    plot_bse_results(results, outputDir);
    audit_bse_targets(root, outputDir);
    fprintf('\nMean per-trial RMSE (Nm):\n');
    for p = 1:2
        for m = 1:3
            mask = strcmp(metrics.protocol, protocols{p}) & metrics.model == m;
            fprintf('%s / %s: %.3f\n', protocols{p}, modelNames{m}, mean(metrics.rmse(mask)));
        end
    end
    fprintf('\nOutputs: %s\n', outputDir);
end

function [a, b, r, y, group] = stack_trials(trials)
    aa = cell(numel(trials), 1); bb = aa; rr = aa; yy = aa; gg = aa;
    for k = 1:numel(trials)
        [aa{k}, bb{k}, rr{k}] = bse_features(trials(k));
        yy{k} = trials(k).torque;
        gg{k} = k*ones(size(yy{k}));
    end
    a = vertcat(aa{:}); b = vertcat(bb{:}); r = vertcat(rr{:});
    y = vertcat(yy{:}); group = vertcat(gg{:});
end

function value = macro_mse(y, estimate, group)
    ids = unique(group);
    errors = zeros(numel(ids), 1);
    for k = 1:numel(ids)
        mask = group == ids(k);
        errors(k) = mean((estimate(mask)-y(mask)).^2);
    end
    value = mean(errors);
end

function value = mean_trial_rmse(y, estimate, group)
    ids = unique(group);
    errors = zeros(numel(ids), 1);
    for k = 1:numel(ids)
        mask = group == ids(k);
        errors(k) = sqrt(mean((estimate(mask)-y(mask)).^2));
    end
    value = mean(errors);
end
