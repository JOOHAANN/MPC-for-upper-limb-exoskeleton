function test_hill_sequence()
    % TEST_HILL_SEQUENCE Continuous mechanics, idle balance and reference audits.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'experiments', 'hill_compliant'));
    data = load(fullfile(root, 'results', 'hill_sequence', 'hill_sequence.mat'));
    package = load(fullfile(root, 'data', 'processed', 'hill_compliant', 'calibrated_models.mat'), 'models');
    originals = {'5_35_10_55_15_45', '5_45_10_35_15_55', '5_55_10_45_15_35'};
    convergence = struct([]);
    stress = struct([]);
    for group = 1:3
        r = data.result.groups{group};
        t = r.sequence;
        raw = load(fullfile(root, 'data', 'processed', ['training_map__muscle_tendon__' originals{group} '_done__build_data'], 'Combination_raw.mat'));
        assert(isequal(t.q, raw.Combination(1:10:end, 3)), 'Original full angle reference changed.');
        assert(isequal(t.torque, raw.Combination(1:10:end, 1)), 'Original torque reference changed.');
        for source = 1:2
            model = package.models{source};
            par = model.parameters;
            p = r.predictions{source};
            a = r.inputs{source};
            stance = struct('waiting', t.waiting, 'parameters', [t.standingAngle, data.result.options.controller]);
            assert(~p.failed && all(isfinite(p.q)) && all(p.externalBalanceTorque == 0));
            assert(max(abs(p.torque - p.backTorque + p.frontTorque)) < 1e-9);
            assert(max(abs(p.backTorque - p.backActive - p.backPassive - p.backFiberDamping)) < 1e-6);
            assert(max(abs(p.frontTorque - p.frontActive - p.frontPassive - p.frontFiberDamping)) < 1e-6);
            % Independently recover abdominal tendon force from geometry and
            % tendon strain; this checks the newton-valued plotted/exported data.
            lo = par(22) / par(11);
            width = sin(par(28));
            slack = par(24) - lo * sqrt(par(12)^2 - width^2);
            lmt = par(24) + par(22) * (p.q - 1.4);
            tendon = lmt - lo * sqrt(max(0, p.frontLengthNormalized.^2 - width^2));
            independentForce = (par(9) / par(22)) * (max(0, (tendon - slack) / (par(26) * slack))).^2;
            assert(max(abs(independentForce - p.muscleForces.frontTendonN)) < 1e-5);
            assert(min(p.muscleForces.frontTendonN) >= 0);
            forceResidual = p.muscleForces.frontTendonN - p.muscleForces.frontActiveN - ...
                            p.muscleForces.frontPassiveN - p.muscleForces.frontFiberDampingN;
            assert(max(abs(forceResidual)) < 1e-5);
            assert(max(p.equilibriumResidual) < 1e-7 && max(abs(p.momentumResidual)) < 1e-4);
            inertia = t.m * par(13)^2 + 4 * par(19) * par(20)^2 / 3;
            assert(max(abs(inertia .* p.acceleration - p.torque - p.gravityTorque - p.dampingTorque - p.boundaryTorque)) < 1e-8);
            assert(all(p.neuralCommands(:) >= 0 & p.neuralCommands(:) <= 1));
            assert(all(p.effectiveActivation(:) >= 0 & p.effectiveActivation(:) <= 1));
            % These bounds follow activation dynamics and forbid activation resets.
            bounds = 1 - exp(-t.dt ./ par([8 27]));
            assert(all(max(abs(diff(p.effectiveActivation)), [], 1) <= bounds + 1e-9));
            for stage = 1:3
                edges = [t.segments(stage).startIndex, t.segments(stage).endIndex + 1];
                for k = edges
                    assert(abs(p.q(k) - p.q(k - 1)) < t.dt * max(abs(p.v)) * 1.25 + 1e-8);
                    % Check sampled momentum near switches independently of the
                    % solver's residual; an instantaneous state reset would fail.
                    rhs = p.torque + p.gravityTorque + p.dampingTorque + p.boundaryTorque;
                    residual = inertia(k) * (p.v(k) - p.v(k - 1)) / t.dt - (rhs(k) + rhs(k - 1)) / 2;
                    assert(abs(residual) < 5, 'Mechanical state jumped at a protocol boundary.');
                end
            end
            fine = hill_pair_rollout(par, a, t.m, [t.q(1), t.v0], t.dt, 20, true, stance);
            assert(~fine.failed);
            rec = struct('group', group, 'source', source, 'angleDifferenceRad', max(abs(p.q - fine.q)), ...
                         'torqueDifferenceNm', max(abs(p.torque - fine.torque)));
            convergence = append(convergence, rec);
            % Future input edits and reference poisoning cannot change a prefix.
            cut = floor(numel(t.q) / 2);
            changedTrial = t;
            changedTrial.q(2:end) = NaN;
            changedTrial.torque(:) = NaN;
            unchanged = hill_pair_inputs(changedTrial, source, model.generator, par);
            assert(isequal(a, unchanged));
            changedTrial.a(cut + 1:end, :) = .8;
            changedTrial.front(cut + 1:end) = .2;
            changed = hill_pair_inputs(changedTrial, source, model.generator, par);
            again = hill_pair_rollout(par, changed, t.m, [t.q(1), t.v0], t.dt, 10, true, stance);
            assert(isequal(p.q(1:cut), again.q(1:cut)), 'Future input leaked into a prefix.');
            % Same initial state and parameters: controller disabled must reduce
            % exactly to the original compliant Hill forward solver.
            n = 250;
            disabled = stance;
            disabled.waiting = zeros(n, 1);
            base = hill_pair_rollout(par, a(1:n, :), t.m(1:n), [t.q(1), t.v0], t.dt, 10);
            passThrough = hill_pair_rollout(par, a(1:n, :), t.m(1:n), [t.q(1), t.v0], t.dt, 10, true, disabled);
            assert(max(abs(base.q - passThrough.q)) < 1e-12 && max(abs(base.torque - passThrough.torque)) < 1e-9);
            % Append 30 seconds of unloaded standing without reinitialization.
            count = 3000;
            held = t;
            held.a = [t.a; zeros(count, 2)];
            held.front = [t.front; repmat(t.front(end), count, 1)];
            held.m = [t.m; zeros(count, 1)];
            held.loadKg = [t.loadKg; repmat(t.loadKg(end), count, 1)];
            held.waiting = [t.waiting; ones(count, 1)];
            extendedInputs = hill_pair_inputs(held, source, model.generator, par);
            extendedStance = stance;
            extendedStance.waiting = held.waiting;
            extended = hill_pair_rollout(par, extendedInputs, held.m, [t.q(1), t.v0], t.dt, 10, true, extendedStance);
            assert(~extended.failed && isequal(p.q, extended.q(1:numel(t.q))));
            tail = numel(extended.q) - 999:numel(extended.q);
            assert(max(abs(extended.q(tail) - t.standingAngle)) < .001);
            assert(max(abs(extended.v(tail))) < .001);
            rec = struct('group', group, 'source', source, 'durationSeconds', (numel(extended.q) - 1) * t.dt, ...
                         'finalAngleErrorRad', max(abs(extended.q(tail) - t.standingAngle)), ...
                         'finalSpeedRadPerSecond', max(abs(extended.v(tail))));
            stress = append(stress, rec);
        end
    end
    % Compare readable MATLAB and C with both movement and wait transitions.
    for source = 1:2
        par = package.models{source}.parameters;
        n = 80;
        x = linspace(0, 1, n)';
        inputs = [.03 + .15 * sin(pi * x).^2, .06 + .2 * sin(pi * x).^2];
        stance = struct('waiting', [zeros(30, 1); ones(50, 1)], ...
                        'parameters', [1.5, data.result.options.controller]);
        p1 = hill_pair_rollout(par, inputs, zeros(n, 1), [1.25, -.3], .01, 2, false, stance);
        p2 = hill_pair_rollout(par, inputs, zeros(n, 1), [1.25, -.3], .01, 2, true, stance);
        assert(~p1.failed && ~p2.failed);
        assert(max(abs(p1.q - p2.q)) < 1e-6 && max(abs(p1.torque - p2.torque)) < .001);
        assert(max(abs(p1.neuralCommands(:) - p2.neuralCommands(:))) < 1e-6);
    end
    checks = struct2table(convergence);
    assert(max(checks.angleDifferenceRad) < .01 && max(checks.torqueDifferenceNm) < 3);
    assert(max(data.result.waitingMetrics.maxSettledAngleErrorRad) < .01);
    assert(all(data.result.movementMetrics.Wmatched == 1));
    writetable(checks, fullfile(root, 'results', 'hill_sequence', 'convergence.csv'));
    writetable(struct2table(stress), fullfile(root, 'results', 'hill_sequence', 'stability.csv'));
    fprintf(['PASS: six continuous sequences, original references, abdominal tendon forces, Hill-only posture balance, ' ...
             'no resets, causality, step convergence, long waits and MATLAB/C parity.\n']);
end

function records = append(records, record)
    if isempty(records)
        records = record;
    else
        records(end + 1) = record;
    end
end
