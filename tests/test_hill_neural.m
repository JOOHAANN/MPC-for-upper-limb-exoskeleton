function test_hill_neural(modelFile)
    % TEST_HILL_NEURAL Physics parity and causality, not physiological proof.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'experiments', 'hill_compliant'));
    if nargin < 1
        modelFile = fullfile(root, 'data', 'processed', 'hill_compliant', 'neural_models.mat');
    end
    d = load(modelFile);
    controller = [200 60 0 120 2 .15];
    if isfield(d, 'controller')
        controller = reshape(d.controller, 1, 6);
    end
    data = load(fullfile(root, 'data', 'processed', 'hill_compliant', 'paired_sequences.mat'));
    assert(isequal(d.models{1}.parameters([13 14 19 20 29:32]), d.models{2}.parameters([13 14 19 20 29:32])));
    for source = 1:2
        model = d.models{source};
        t = data.sequences(1);
        a = hill_pair_inputs(t, source, model.generator, model.parameters);
        % Include the loaded phase so tonic load compensation is exercised.
        ids = 601:1050;
        assert(any(t.m(ids) > 1));
        initial = [t.q(ids(1)), 0];
        stance = struct('waiting', t.waiting(ids), 'parameters', [t.standingAngle controller], ...
                        'neuralParameters', model.neuralParameters);
        p = hill_pair_rollout(model.parameters, a(ids, :), t.m(ids), initial, .01, 2, true, stance);
        readable = hill_pair_rollout(model.parameters, a(ids, :), t.m(ids), initial, .01, 2, false, stance);
        assert(~p.failed && ~readable.failed);
        assert(max(abs(p.q - readable.q)) < 1e-5 && max(abs(p.torque - readable.torque)) < .01);
        assert(max(abs(p.motionIntent - readable.motionIntent)) < 1e-10);
        assert(max(abs(p.torque - p.backTorque + p.frontTorque)) < 1e-9);
        assert(max(abs(p.frontTorque - p.frontActive - p.frontPassive - p.frontFiberDamping)) < 1e-6);
        assert(all(p.externalBalanceTorque == 0));
        if numel(model.neuralParameters) == 10 && model.neuralParameters(10) > 1e-4
            uncompensated = stance;
            uncompensated.neuralParameters(10) = 0;
            noTonic = hill_pair_rollout(model.parameters, a(ids, :), t.m(ids), initial, .01, 2, true, uncompensated);
            assert(max(p.motionIntent - noTonic.motionIntent) > 1e-6, 'Load compensation had no effect on intent.');
        end
        cut = 150;
        changed = a(ids, :);
        changed(cut + 1:end, :) = .7;
        future = hill_pair_rollout(model.parameters, changed, t.m(ids), initial, .01, 2, true, stance);
        assert(isequal(p.q(1:cut), future.q(1:cut)), 'Future EMG changed past states.');
        t.q(:) = NaN;
        t.torque(:) = NaN;
        withoutReferences = hill_pair_inputs(t, source, model.generator, model.parameters);
        assert(isequal(a, withoutReferences), 'Reference entered excitation mapping.');
        fprintf('Source %d: C/MATLAB parity, torque balance, no external PD, causal inputs.\n', source);
    end
    fprintf('PASS: neural/Hill numerical checks; these are not independent muscle-force validation.\n');
end
