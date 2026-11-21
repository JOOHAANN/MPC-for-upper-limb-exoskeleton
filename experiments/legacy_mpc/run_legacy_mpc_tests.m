function reports = run_legacy_mpc_tests(root, sources, mode)
    % RUN_LEGACY_MPC_TESTS Execute archived simulators in isolated workspaces.
    % Model equations, controller settings and archived inputs are preserved.
    % Any exception is recorded together with the actual completed trajectory.
    if nargin < 2
        sources = {'control_experimental', 'control_map'};
    end
    if nargin < 3
        mode = 'checked';
    end
    assert(ismember(mode, {'original', 'checked', 'reference_only'}));
    addpath(root);
    oldDirectory = pwd;
    oldVisibility = get(groot, 'DefaultFigureVisible');
    cleanup = onCleanup(@() restore_environment(oldDirectory, oldVisibility)); %#ok<NASGU>
    set(groot, 'DefaultFigureVisible', 'off');
    reports = struct([]);
    out = fullfile(root, 'results', 'legacy_mpc');
    if ~isfolder(out)
        mkdir(out);
    end
    for source = 1:numel(sources)
        workDir = prepare_workflow(sources{source});
        cd(workDir);
        clear nonlinear_system_period mpc_cost_function_increment;
        assert(strcmp(which('nonlinear_system_period'), fullfile(workDir, 'nonlinear_system_period.m')), ...
               'The selected legacy predictor must resolve inside its own workflow.');
        manifest = jsondecode(fileread('workflow_manifest.json'));
        if strcmp(sources{source}, 'control_experimental')
            entry = 'offline_simulator_NMPC_quasimodo_incremental_MPC';
        else
            entry = 'offline_simulator_NMPC_quasimodo_incremental_MPC_with_MAP';
        end
        originalEntry = entry;
        if ~strcmp(mode, 'original')
            entry = make_checked_entry(workDir, entry, mode);
        end
        diary(fullfile(workDir, 'legacy_run.log'));
        rng(0, 'twister');
        timer = tic;
        report = execute_entry(entry);
        report.elapsedSeconds = toc(timer);
        report.source = sources{source};
        report.originalFolder = manifest.original_folder;
        report.workDirectory = workDir;
        report.simulator = originalEntry;
        report.executedEntry = entry;
        report.mode = mode;
        report.predictor = fullfile(workDir, 'nonlinear_system_period.m');
        report.modelPolicy = 'Archived BSE/controller parameters; not the selected neural/Hill model';
        report.executionChanges = {};
        if ~strcmp(mode, 'original')
            report.executionChanges = {'Serial finite differences; same fmincon SQP settings', ...
                'Hold the final recorded sample at the prediction-window tail', ...
                'Use zero initial assistance instead of indexing U(0)', ...
                'Record solver diagnostics; omit legacy plotting block'};
            if strcmp(mode, 'reference_only')
                report.executionChanges{end + 1} = 'Control target uses unchanged Combination(:,3), not archived Theta_Tac0';
            end
        end
        diary off;
        write_json(fullfile(workDir, 'execution_status.json'), report);
        write_json(fullfile(out, [sources{source} '_' mode '_status.json']), report);
        reports = [reports; report]; %#ok<AGROW>
        fprintf('%s: complete=%d; last sample=%d; wall time=%.1f s\n', ...
                report.source, report.completed, report.lastStep, report.elapsedSeconds);
    end
end

function report = execute_entry(entry)
    % The legacy script clears its own workspace, so keep bookkeeping outside.
    try
        run(entry);
        report = struct('completed', true, 'errorIdentifier', '', 'errorMessage', '');
    catch legacyError
        report = struct('completed', false, 'errorIdentifier', legacyError.identifier, ...
                        'errorMessage', getReport(legacyError, 'extended', 'hyperlinks', 'off'));
        fprintf('%s\n', report.errorMessage);
    end
    report.lastStep = 0;
    report.expectedSteps = 0;
    if exist('k', 'var')
        report.lastStep = k;
    end
    if exist('NumberOfData', 'var')
        report.expectedSteps = NumberOfData - 1;
    end
    variables = {'Combination', 'coeffs', 'Ts', 'PH', 'CH', 'n', 'Threshold', 'u_max', ...
                 'activation_sew', 'theta', 'dtheta', 'MTF', 'Tac', 'U', 't', ...
                 'theta_Ceinms', 'Torque_Ceinms', 'theta_real', 'Torque_real', ...
                 'T_balance', 'F_k_up', 'T_front', 'k', 'NumberOfData', ...
                 'legacyExitFlag', 'legacySolveSeconds', 'legacyIterations', ...
                 'legacyEvaluations', 'legacyObjective', 'legacyPredictionMax'};
    present = false(size(variables));
    for index = 1:numel(variables)
        present(index) = exist(variables{index}, 'var') == 1;
    end
    save('legacy_execution.mat', variables{present}, 'report');
end

function entry = make_checked_entry(directory, originalEntry, mode)
    % Execution-only patches are applied to a generated copy, not source code.
    code = fileread(fullfile(directory, [originalEntry '.m']));
    cut = strfind(code, 'model_torque = MTF(1:end - 1)');
    assert(numel(cut) == 1, 'Unexpected legacy plotting boundary.');
    lineEnd = find(code(cut:end) == newline, 1);
    code = code(1:cut + lineEnd - 1);
    code = strrep(code, '''UseParallel'', true', '''UseParallel'', false');
    code = strrep(code, 'u0 = U(k - 1);', 'u0 = U(max(1, k - 1));');
    if strcmp(mode, 'reference_only')
        assert(contains(code, 'theta_real = Theta_Tac0;'));
        code = strrep(code, 'theta_real = Theta_Tac0;', 'theta_real = theta_Ceinms;');
    end
    names = {'dtheta_real', 'theta_real', 'activation_sew'};
    for i = 1:numel(names)
        code = strrep(code, [names{i} '(k:k + PH)'], ...
                     [names{i} '(min(k:k + PH, NumberOfData))']);
    end
    init = ['legacyExitFlag = NaN(NumberOfData, 1);' newline ...
            'legacySolveSeconds = zeros(NumberOfData, 1);' newline ...
            'legacyIterations = zeros(NumberOfData, 1);' newline ...
            'legacyEvaluations = zeros(NumberOfData, 1);' newline ...
            'legacyObjective = NaN(NumberOfData, 1);' newline ...
            'legacyPredictionMax = NaN(NumberOfData, 1);' newline ...
            'legacyRunTimer = tic; legacyProgressTimer = tic;' newline];
    code = strrep(code, 'for k = 1:(NumberOfData - 1)', ...
                 [init 'for k = 1:(NumberOfData - 1)']);
    code = strrep(code, 'if mod(k, 100) == 0', 'if mod(k, 1000) == 0');
    code = strrep(code, 'if max(Predict_torque(:, 1)) > Threshold(k)', ...
                 ['legacyPredictionMax(k) = max(Predict_torque(:, 1));' newline ...
                  '        if max(Predict_torque(:, 1)) > Threshold(k)']);
    code = strrep(code, 'du = fmincon(', ...
                 ['legacySolveTimer = tic;' newline ...
                  '            [du, legacyValue, legacyFlag, legacyOutput] = fmincon(']);
    log = ['legacyExitFlag(k) = legacyFlag;' newline ...
           '            legacySolveSeconds(k) = toc(legacySolveTimer);' newline ...
           '            legacyIterations(k) = legacyOutput.iterations;' newline ...
           '            legacyEvaluations(k) = legacyOutput.funcCount;' newline ...
           '            legacyObjective(k) = legacyValue;' newline ...
           '            if toc(legacyProgressTimer) > 15' newline ...
           '                fprintf(''Sample %d/%d; wall %.1f s; solve %.4f s; exit %d\n'', k, NumberOfData, toc(legacyRunTimer), legacySolveSeconds(k), legacyFlag);' newline ...
           '                legacyProgressTimer = tic;' newline ...
           '            end' newline ...
           '            Tac(k + 1) = u0 + du(1);'];
    code = strrep(code, 'Tac(k + 1) = u0 + du(1);', log);
    if contains(originalEntry, 'with_MAP')
        entry = ['legacy_mpc_map_' mode];
    else
        entry = ['legacy_mpc_experimental_' mode];
    end
    file = fopen(fullfile(directory, [entry '.m']), 'w', 'n', 'UTF-8');
    assert(file >= 0);
    cleanup = onCleanup(@() fclose(file)); %#ok<NASGU>
    fprintf(file, '%s', code);
end

function write_json(path, value)
    file = fopen(path, 'w', 'n', 'UTF-8');
    assert(file >= 0, 'Could not write the execution report.');
    cleanup = onCleanup(@() fclose(file)); %#ok<NASGU>
    fprintf(file, '%s\n', jsonencode(value, 'PrettyPrint', true));
end

function restore_environment(directory, visibility)
    cd(directory);
    set(groot, 'DefaultFigureVisible', visibility);
end
