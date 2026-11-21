% Run the fixed MPC on both activation sources and all nine conditions.
% Requires Python with numpy/scipy/matplotlib and a GCC C compiler.
% Angle/torque references are not used by the human plant or MPC preview.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'hill_compliant'));
addpath(fullfile(root, 'experiments', 'legacy_mpc'));
model = load(fullfile(root, 'data', 'processed', 'hill_compliant', 'neural_models.mat'));
replayFile = fullfile(root, 'results', 'hill_neural_sequence', 'hill_sequence.mat');
replayMatches = false;
if isfile(replayFile)
    replay = load(replayFile, 'result');
    replayMatches = true;
    for source = 1:2
        replayMatches = replayMatches && ...
            isequal(model.models{source}.parameters, replay.result.models{source}.parameters) && ...
            isequal(model.models{source}.neuralParameters, replay.result.models{source}.neuralParameters);
    end
end
if ~replayMatches
    run_hill_sequence(root, struct('modelFile', 'neural_models.mat', ...
        'outputDirectory', 'hill_neural_sequence', 'renderFigures', false));
end

pythonExecutable = find_hill_mpc_python();
for script = {'run_hill_mpc.py', 'evaluate_hill_mpc.py'}
    scriptPath = fullfile(root, 'experiments', 'legacy_mpc', script{1});
    [status, output] = system(sprintf('"%s" "%s"', pythonExecutable, scriptPath), '-echo');
    assert(status == 0, 'Hill MPC backend failed: %s', output);
end
hillMPCResults = load(fullfile(root, 'results', 'hill_mpc', 'hill_mpc.mat'));
demo_hill_mpc_plots(false);
fprintf('Fixed MPC results: %s\n', fullfile(root, 'results', 'hill_mpc'));
