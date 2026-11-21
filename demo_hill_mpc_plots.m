function demo_hill_mpc_plots(regenerate)
% Display all saved optimized-model/MPC comparisons in one tabbed window.
% No model fitting or simulation is performed. Optional false skips rendering.
if nargin < 1, regenerate = true; end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'legacy_mpc'));
outputDirectory = fullfile(root, 'results', 'hill_mpc');
if regenerate
    assert(isfile(fullfile(outputDirectory, 'hill_mpc.mat')), ...
        'Saved MPC trajectories are missing. Run demo_hill_mpc first.');
    pythonExecutable = find_hill_mpc_python();
    scriptPath = fullfile(root, 'experiments', 'legacy_mpc', 'plot_hill_mpc_results.py');
    [status, output] = system(sprintf('"%s" "%s"', pythonExecutable, scriptPath), '-echo');
    assert(status == 0, 'Hill/MPC plotting failed: %s', output);
end
names = {'mpc_torque_theta_group1', 'mpc_torque_theta_group2', 'mpc_torque_theta_group3', ...
         'theta_all9_experimental', 'torque_all9_experimental', ...
         'theta_all9_map', 'torque_all9_map', ...
         'mpc_muscle_and_boundary_torque_group1', 'mpc_muscle_and_boundary_torque_group2', ...
         'mpc_muscle_and_boundary_torque_group3', 'theta_15kg_all_speeds'};
labels = {'Sequence 1', 'Sequence 2', 'Sequence 3', ...
          'Experimental: theta', 'Experimental: torque', 'MAP: theta', 'MAP: torque', ...
          'Muscles 1', 'Muscles 2', 'Muscles 3', '15 kg detail'};
if ~regenerate && ~all(cellfun(@(name) isfile(fullfile(outputDirectory, [name '.png'])), names))
    % A fresh checkout includes published figures but no ignored result files.
    outputDirectory = fullfile(root, 'docs', 'hill_mpc_evaluation');
end
for index = 1:numel(names)
    assert(isfile(fullfile(outputDirectory, [names{index} '.png'])), ...
        'Missing figure %s. Run demo_hill_mpc_plots to render all figures.', names{index});
end
window = figure('Name', 'Optimized Hill model and MPC: all nine conditions', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'Position', [.04 .06 .92 .84]);
tabs = uitabgroup(window);
for index = 1:numel(names)
    panel = uitab(tabs, 'Title', labels{index});
    axesHandle = axes('Parent', panel, 'Position', [.005 .005 .99 .99]);
    image(axesHandle, imread(fullfile(outputDirectory, [names{index} '.png'])));
    axis(axesHandle, 'image');
    axis(axesHandle, 'off');
end
fprintf('All figures: %s\n', outputDirectory);
fprintf('Combined PDF: %s\n', fullfile(outputDirectory, 'optimized_model_and_mpc_recovery.pdf'));
end
