function result = run_hill_sequence(root, options)
    % RUN_HILL_SEQUENCE Continuous 5/10/15-kg movements and standing waits.
    % All mechanical, activation and input-memory states persist throughout.
    if nargin < 2
        options = struct();
    end
    defaults = struct('groups', 1:3, 'steps', 10, 'useMex', true, 'renderFigures', true, ...
                      'controller', [], ...
                      'modelFile', 'calibrated_models.mat', 'outputDirectory', 'hill_sequence');
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(options, names{k})
            options.(names{k}) = defaults.(names{k});
        end
    end
    models = load(fullfile(root, 'data', 'processed', 'hill_compliant', options.modelFile));
    if isempty(options.controller)
        if isfield(models, 'controller')
            options.controller = reshape(models.controller, 1, 6);
        else
            options.controller = [200 60 0 120 2 .15];
        end
    end
    data = load(fullfile(root, 'data', 'processed', 'hill_compliant', 'paired_sequences.mat'));
    assert(all(ismember(options.groups, 1:3)));
    out = fullfile(root, 'results', options.outputDirectory);
    if ~isfolder(out)
        mkdir(out);
    end
    labels = {'Experimental activation', 'MAP activation'};
    records = struct([]);
    waitingRecords = struct([]);
    forceRecords = struct([]);
    torqueRecords = struct([]);
    result = struct('groups', {cell(1, 3)}, 'options', options, 'models', {models.models});
    for group = options.groups
        t = data.sequences(group);
        predictions = cell(1, 2);
        unbalanced = predictions;
        inputs = predictions;
        for source = 1:2
            model = models.models{source};
            a = hill_pair_inputs(t, source, model.generator, model.parameters);
            stance = struct('waiting', t.waiting, 'parameters', [t.standingAngle, options.controller]);
            if isfield(model, 'neuralParameters')
                stance.neuralParameters = model.neuralParameters;
                if isfield(model, 'intentDecoder'), stance.intentDecoder = model.intentDecoder; end
            end
            p = hill_pair_rollout(model.parameters, a, t.m, [t.q(1), t.v0], t.dt, options.steps, options.useMex, stance);
            assert(~p.failed, 'Continuous Hill rollout failed. No curve replacement is applied.');
            p.muscleForces = hill_muscle_forces(model.parameters, p);
            predictions{source} = p;
            inputs{source} = a;
            disabled = stance;
            disabled.waiting(:) = 0;
            if isfield(disabled, 'neuralParameters')
                disabled.neuralParameters(7) = 0;
            end
            unbalanced{source} = hill_pair_rollout(model.parameters, a, t.m, [t.q(1), t.v0], t.dt, options.steps, options.useMex, disabled);
            % A disabled-control counterfactual may fail; never show its
            % sentinel tail as a physically valid predicted trajectory.
            settled = logical(t.waiting);
            for stage = 1:3
                segment = t.segments(stage);
                selected = segment.startIndex:segment.endIndex;
                [~, valleys] = findpeaks(-p.q(selected), 'MinPeakProminence', .12, 'MinPeakDistance', 50);
                referenceValleys = findpeaks(-t.q(selected), 'MinPeakProminence', .12, 'MinPeakDistance', 50);
                rec = struct('group', group, 'source', labels{source}, 'loadKg', segment.loadKg, 'speedCode', segment.speedCode, ...
                             'thetaRMSEDegrees', sqrt(mean((p.q(selected) - t.q(selected)).^2)) * 180 / pi, ...
                             'torqueRMSENm', sqrt(mean((p.torque(selected) - t.torque(selected)).^2)), ...
                             'predictedValleys', numel(valleys), 'Wmatched', double(numel(valleys) == 2 && numel(referenceValleys) == 2));
                records = append_record(records, rec);
                firstWait = segment.endIndex + 1;
                settled(firstWait:min(numel(settled), segment.endIndex + round(2 / t.dt))) = false;
            end
            rec = struct('group', group, 'source', labels{source}, 'standingAngleRad', t.standingAngle, ...
                         'settlingAllowanceSeconds', 2, 'maxSettledAngleErrorRad', max(abs(p.q(settled) - t.standingAngle)), ...
                         'meanSettledAngleErrorRad', mean(abs(p.q(settled) - t.standingAngle)), ...
                         'maxSettledSpeedRadPerSecond', max(abs(p.v(settled))), ...
                         'meanSettledHillTorqueNm', mean(p.torque(settled)), ...
                         'expectedUnloadedGravityNm', 9.81 * model.parameters(19) * model.parameters(20) * cos(t.standingAngle), ...
                         'minAngleRad', min(p.q), 'maxAngleRad', max(p.q), ...
                         'softLimitPeakNm', max(abs(p.boundaryTorque)), 'externalBalancePeakNm', max(abs(p.externalBalanceTorque)), ...
                         'forceResidualMax', max(p.equilibriumResidual), 'momentumResidualMax', max(abs(p.momentumResidual)));
            waitingRecords = append_record(waitingRecords, rec);
            f = p.muscleForces;
            active = ~logical(t.waiting);
            rec = struct('group', group, 'source', labels{source}, ...
                         'frontMomentArmM', model.parameters(22), 'frontF0N', model.parameters(9) / model.parameters(22), ...
                         'frontTendonPeakN', max(f.frontTendonN), 'frontMovementMeanN', mean(f.frontTendonN(active)), ...
                         'frontSettledStandingMeanN', mean(f.frontTendonN(settled)), ...
                         'frontActivePeakN', max(f.frontActiveN), 'frontPassivePeakN', max(f.frontPassiveN), ...
                         'frontFiberDampingAbsPeakN', max(abs(f.frontFiberDampingN)), ...
                         'backTendonPeakN', max(f.backTendonN), ...
                         'frontComponentResidualN', max(abs(f.frontTendonN - f.frontActiveN - f.frontPassiveN - f.frontFiberDampingN)));
            forceRecords = append_record(forceRecords, rec);
            rec = struct('group', group, 'source', labels{source}, ...
                         'frontMomentArmM', model.parameters(22), ...
                         'abdominalFlexionTorquePeakNm', max(p.frontTorque), ...
                         'abdominalMovementMeanTorqueNm', mean(p.frontTorque(active)), ...
                         'abdominalSettledStandingMeanTorqueNm', mean(p.frontTorque(settled)), ...
                         'abdominalActivePeakNm', max(p.frontActive), 'abdominalPassivePeakNm', max(p.frontPassive), ...
                         'abdominalFiberDampingAbsPeakNm', max(abs(p.frontFiberDamping)), ...
                         'abdominalJointContributionMinNm', min(-p.frontTorque), ...
                         'abdominalJointContributionMaxNm', max(-p.frontTorque), ...
                         'backHillTorquePeakNm', max(p.backTorque), ...
                         'torqueComponentResidualNm', max(abs(p.frontTorque - p.frontActive - p.frontPassive - p.frontFiberDamping)));
            torqueRecords = append_record(torqueRecords, rec);
            traces = table(t.t, t.loadKg, t.m, t.waiting, t.referenceIsFilled, t.q, p.q, p.v, t.torque, p.torque, ...
                           p.backTorque, p.frontTorque, p.gravityTorque, p.dampingTorque, p.boundaryTorque, p.externalBalanceTorque, ...
                           p.balanceWeight, a(:, 1), a(:, 2), p.neuralCommands(:, 1), p.neuralCommands(:, 2), ...
                           p.effectiveActivation(:, 1), p.effectiveActivation(:, 2), p.backLengthNormalized, p.frontLengthNormalized, ...
                           p.requestedHillTorque, p.postureCorrectionRequest, p.equilibriumResidual, p.momentumResidual, ...
                           'VariableNames', {'time', 'nominalLoadKg', 'effectiveLoadKg', 'waiting', 'referenceIsFilled', ...
                                  'referenceAngle', 'predictedAngle', 'predictedSpeed', 'referenceTorque', 'netHillTorque', ...
                                  'backHillTorque', 'frontHillTorque', 'gravityTorque', 'jointDampingTorque', 'softLimitTorque', 'externalBalanceTorque', ...
                                  'stanceBlend', 'recordedOrGeneratedBackInput', 'mappedFrontInput', 'commandBack', 'commandFront', ...
                                  'backHillActivation', 'frontHillActivation', 'backFiberLengthNormalized', 'frontFiberLengthNormalized', ...
                                  'requestedHillTorque', 'postureCorrectionRequest', 'forceEquilibriumResidual', 'momentumResidual'});
            traces.backTendonForceN = f.backTendonN;
            traces.frontTendonForceN = f.frontTendonN;
            traces.backActiveForceN = f.backActiveN;
            traces.backPassiveForceN = f.backPassiveN;
            traces.backFiberDampingForceN = f.backFiberDampingN;
            traces.frontActiveForceN = f.frontActiveN;
            traces.frontPassiveForceN = f.frontPassiveN;
            traces.frontFiberDampingForceN = f.frontFiberDampingN;
            traces.abdominalFlexionTorqueNm = p.frontTorque;
            traces.abdominalJointContributionNm = -p.frontTorque;
            traces.abdominalActiveTorqueNm = p.frontActive;
            traces.abdominalPassiveTorqueNm = p.frontPassive;
            traces.abdominalFiberDampingTorqueNm = p.frontFiberDamping;
            if isfield(p, 'motionIntent')
                traces.motionIntentAngleRad = p.motionIntent;
            end
            writetable(traces, fullfile(out, sprintf('group%d_source%d.csv', group, source)));
        end
        r = struct('group', group, 'sequence', t, 'predictions', {predictions}, 'unbalanced', {unbalanced}, 'inputs', {inputs});
        result.groups{group} = r;
        if options.renderFigures
            plot_hill_sequence(r, out);
            plot_hill_abdominal_comparison(r, models.models, out);
        end
    end
    result.movementMetrics = struct2table(records);
    result.waitingMetrics = struct2table(waitingRecords);
    result.muscleForceMetrics = struct2table(forceRecords);
    result.muscleTorqueMetrics = struct2table(torqueRecords);
    writetable(result.movementMetrics, fullfile(out, 'movement_metrics.csv'));
    writetable(result.waitingMetrics, fullfile(out, 'waiting_metrics.csv'));
    writetable(result.muscleForceMetrics, fullfile(out, 'muscle_force_metrics.csv'));
    writetable(result.muscleTorqueMetrics, fullfile(out, 'muscle_torque_metrics.csv'));
    save(fullfile(out, 'hill_sequence.mat'), 'result');
    fprintf('\nContinuous sequence: one initialization per source/sequence; no state reset.\n');
    disp(groupsummary(result.movementMetrics, 'source', 'mean', {'thetaRMSEDegrees', 'torqueRMSENm', 'Wmatched'}));
    disp(result.waitingMetrics(:, {'group', 'source', 'maxSettledAngleErrorRad', 'meanSettledHillTorqueNm', 'externalBalancePeakNm'}));
end

function records = append_record(records, rec)
    if isempty(records)
        records = rec;
    else
        records(end + 1) = rec;
    end
end
