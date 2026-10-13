function test_bse_generalization()
%TEST_BSE_GENERALIZATION Check target isolation and grouped evaluation records.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'experiments', 'bse_generalization'));
    d = load(fullfile(root, 'data', 'processed', 'bse_generalization', ...
        'raw_sensor_trials.mat'), 'trials');
    trials = d.trials;
    assert(numel(trials) == 25 && numel(unique({trials.id})) == 25);
    % Poisoning the target must not change any predictor.
    trial = trials(1);
    [a0, p0, r0] = bse_features(trial);
    trial.torque(:) = NaN;
    [a1, p1, r1] = bse_features(trial);
    assert(isequal(a0, a1) && isequal(p0, p1) && isequal(r0, r1));
    % Repeating a trial does not change the trial-balanced fitted model.
    x = [0; 1; 2; 3]; y = [1; 3; 4; 8]; g = [1; 1; 2; 2];
    model0 = bse_ridge('fit', x, y, g, 0.1);
    model1 = bse_ridge('fit', [x; x(1:2)], [y; y(1:2)], [g; g(1:2)], 0.1);
    assert(max(abs(model0.beta-model1.beta)) < 1e-10);
    saved = fullfile(root, 'results', 'bse_generalization', 'benchmark_results.mat');
    if isfile(saved)
        d = load(saved, 'results'); results = d.results;
        for protocol = {'load', 'speed'}
            coverage = zeros(1, 25);
            for fold = results.selected
                if ~strcmp(fold.protocol, protocol{1}), continue; end
                assert(isempty(intersect(fold.trainTrials, fold.testTrials)));
                assert(numel(union(fold.trainTrials, fold.testTrials)) == 25);
                coverage(fold.testTrials) = coverage(fold.testTrials) + 1;
                if strcmp(protocol{1}, 'load')
                    trainLevels = [trials(fold.trainTrials).loadKg];
                    testLevels = [trials(fold.testTrials).loadKg];
                else
                    trainLevels = [trials(fold.trainTrials).speedLabel];
                    testLevels = [trials(fold.testTrials).speedLabel];
                end
                assert(isempty(intersect(trainLevels, testLevels)));
            end
            assert(all(coverage == 1));
        end
    end
    fprintf('PASS: target-free features, equal trial weighting, disjoint grouped folds.\n');
end
