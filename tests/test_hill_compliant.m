function test_hill_compliant()
    % TEST_HILL_COMPLIANT Independent force, geometry, momentum and causality checks.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'experiments', 'hill_compliant'));
    d = load(fullfile(root, 'results', 'hill_compliant', 'hill_compliant.mat'));
    inputAudit = load(fullfile(root, 'data', 'processed', 'hill_compliant', 'calibrated_models.mat'), 'calibrationInputAudit', 'validationInputAudit');
    legacy = load(fullfile(root, 'data', 'processed', 'hill_forward', 'calibrated_models.mat'), 'trials');
    for k = 1:9
        assert(isequal(d.result.trials(k).q, legacy.trials(k).q), 'Reference was modified.');
        assert(d.result.trials(k).v0 == legacy.trials(k).v0, 'Initial velocity no longer matches the original reference.');
    end
    checks = struct([]);
    allResults = {d.result, d.result.validation};
    for mode = 1:2
        r = allResults{mode};
        for source = 1:2
            for k = 1:9
                model = model_for_trial(r, source, k);
                par = model.parameters;
                dynamics = par;
                t = r.trials(k);
                p = r.predictions{source, k};
                a = r.inputs{source, k};
                expectedInputs = inputAudit.calibrationInputAudit;
                if mode == 2
                    expectedInputs = inputAudit.validationInputAudit;
                end
                assert(max(abs(a(:) - reshape(expectedInputs{source, k}, [], 1))) < 1e-12, 'Python/MATLAB input mapping differs.');
                assert(~p.failed && all(isfinite(p.q)));
                inside = p.q >= par(31) & p.q <= par(32);
                assert(all(p.boundaryTorque(inside) == 0));
                assert(all(p.boundaryTorque(p.q < par(31)) >= 0) && all(p.boundaryTorque(p.q > par(32)) <= 0));
                assert(max(abs(p.boundaryTorque)) <= par(30) + 1e-9);
                assert(max(abs(p.torque - p.backTorque + p.frontTorque)) < 1e-9);
                assert(max(abs(p.backTorque - p.backActive - p.backPassive - p.backFiberDamping)) < 1e-6);
                assert(max(abs(p.frontTorque - p.frontActive - p.frontPassive - p.frontFiberDamping)) < 1e-6);
                % Tendon stretch must close the length geometry of BOTH units.
                units = [par([1 2 3 7]); par([9 11 12 26])];
                arms = par(21:22);
                lengthRef = par(23:24);
                for u = 1:2
                    cap = units(u, 1);
                    ratio = units(u, 2);
                    lref = units(u, 3);
                    strainAtF0 = units(u, 4);
                    lo = arms(u) / ratio;
                    ln = p.backLengthNormalized;
                    torque = p.backTorque;
                    if u == 2
                        ln = p.frontLengthNormalized;
                        torque = p.frontTorque;
                    end
                    assert(all(ln >= sin(par(28)) - 1e-12), 'Fiber length violates constant-height geometry.');
                    projected = lo * sqrt(max(0, ln.^2 - sin(par(28))^2));
                    slack = lengthRef(u) - lo * sqrt(lref^2 - sin(par(28))^2);
                    Lmt = lengthRef(u) - (3 - 2 * u) * arms(u) * (p.q - 1.4);
                    tendon = Lmt - projected;
                    assert(all(tendon > 0), 'Tendon length must remain positive.');
                    % A slack tendon has zero force and can be shorter than its
                    % slack length; use the forward law rather than invert F=0.
                    expectedForce = max(0, (tendon - slack) / (strainAtF0 * slack)).^2;
                    assert(max(abs(expectedForce - torque / cap)) < 1e-8, 'Tendon force and length geometry disagree.');
                end
                assert(max(p.equilibriumResidual) < 1e-7 && max(abs(p.momentumResidual)) < 1e-4);
                I = t.m * par(13)^2 + 4 * par(19) * par(20)^2 / 3;
                assert(max(abs(I .* p.acceleration - p.torque - p.gravityTorque - p.dampingTorque - p.boundaryTorque)) < 1e-8);
                q0 = [t.q(1), t.v0];
                t.q(:) = NaN;
                t.torque(:) = NaN;
                t.qArchived(:) = NaN;
                a2 = hill_pair_inputs(t, source, model.generator, par);
                assert(isequal(a, a2));
                again = hill_pair_rollout(dynamics, a2, t.m, q0, .01, 10);
                assert(isequaln(p, again));
                fine = hill_pair_rollout(dynamics, a2, t.m, q0, .01, 20);
                rec = struct('mode', mode, 'source', source, 'trial', k, 'angleDifference', max(abs(p.q - fine.q)), ...
                             'torqueDifference', max(abs(p.torque - fine.torque)));
                if isempty(checks)
                    checks = rec;
                else
                    checks(end + 1) = rec;
                end %#ok<AGROW>
            end
            % Future EMG changes must leave the already-computed prefix unchanged.
            model = model_for_trial(r, source, 1);
            par = model.parameters;
            dynamics = par;
            t = r.trials(1);
            cut = floor(numel(t.q) / 2);
            a = hill_pair_inputs(t, source, model.generator, par);
            q0 = [t.q(1), t.v0];
            p = hill_pair_rollout(dynamics, a, t.m, q0, .01, 10);
            t.a(cut + 1:end, :) = .8;
            t.front(cut + 1:end) = .2;
            a2 = hill_pair_inputs(t, source, model.generator, par);
            changed = hill_pair_rollout(dynamics, a2, t.m, q0, .01, 10);
            assert(isequal(p.q(1:cut), changed.q(1:cut)), 'Future input leaked into the predicted prefix.');
            % Readable MATLAB fallback and compiled solver must agree.
            n = 80;
            p1 = hill_pair_rollout(dynamics, a(1:n, :), t.m(1:n), q0, .01, 2, false);
            p2 = hill_pair_rollout(dynamics, a(1:n, :), t.m(1:n), q0, .01, 2, true);
            assert(max(abs(p1.q - p2.q)) < 1e-6 && max(abs(p1.torque - p2.torque)) < .001);
        end
    end
    checks = struct2table(checks);
    writetable(checks, fullfile(root, 'results', 'hill_compliant', 'convergence.csv'));
    assert(max(checks.angleDifference) < .01 && max(checks.torqueDifference) < 3);
    stressRecords = struct([]);
    for source = 1:2
        model = model_for_trial(d.result, source, 1);
        representativeTrials = [1 8 6];
        for profile = 1:3
            trialId = representativeTrials(profile);
            t = d.result.trials(trialId);
            repeat = ceil(6001 / numel(t.q));
            t.a = repmat(t.a, repeat, 1);
            t.front = repmat(t.front, repeat, 1);
            t.m = repmat(t.m, repeat, 1);
            inputs = hill_pair_inputs(t, source, model.generator, model.parameters);
            p = hill_pair_rollout(model.parameters, inputs, t.m, [t.q(1), t.v0], .01, 10);
            assert(~p.failed && all(isfinite(p.q)), 'Long numerical stress failed.');
            rec = struct('source', source, 'trial', trialId, 'durationSeconds', (numel(t.m) - 1) * .01, ...
                         'minimumDegrees', min(p.q) * 180 / pi, 'maximumDegrees', max(p.q) * 180 / pi, ...
                         'outside18RadPercent', 100 * mean(p.q < 0 | p.q > 1.8));
            if isempty(stressRecords)
                stressRecords = rec;
            else
                stressRecords(end + 1) = rec;
            end %#ok<AGROW>
        end
    end
    writetable(struct2table(stressRecords), fullfile(root, 'results', 'hill_compliant', 'stability.csv'));
    % Body/load geometry and the same boundary mechanism belong to one subject.
    common = [13 14 19 20 29:32];
    assert(isequal(d.result.models{1}.parameters(common), d.result.models{2}.parameters(common)));
    for fold = 1:3
        assert(isequal(d.result.validation.models{1, fold}.parameters(common), d.result.validation.models{2, fold}.parameters(common)));
    end
    % Independently exercise both branches of the spring, even if a specific
    % fitted trial remains entirely inside the free-motion interval.
    par = d.result.models{1}.parameters;
    q = [par(31) - .2, par(31), mean(par(31:32)), par(32), par(32) + .2];
    torque = hill_soft_limit(q, par);
    assert(torque(1) > 0 && all(torque(2:4) == 0) && torque(5) < 0);
    assert(max(abs(torque)) < par(30));
    eps = 1e-6;
    assert(abs(hill_soft_limit(par(31) - eps, par) / eps) < .1);
    assert(abs(hill_soft_limit(par(32) + eps, par) / eps) < .1);
    fprintf('PASS: both Hill forces, tendon geometry, momentum convergence, elastic soft limits, shared geometry, reference isolation, causal inputs and C/MATLAB parity.\n');
end

function model = model_for_trial(r, source, k)
    if isempty(r.heldSpeeds)
        model = r.models{source};
    else
        model = r.models{source, find(r.heldSpeeds == r.trials(k).speed, 1)};
    end
end
