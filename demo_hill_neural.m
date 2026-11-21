% Replay fixed-parameter neural/Hill models across all nine load/speed cases.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'hill_compliant'));
options = struct('modelFile', 'neural_models.mat', 'outputDirectory', 'hill_neural_sequence');
neuralSequenceResults = run_hill_sequence(root, options);
figure('Name', 'Nine conditions: same neural/Hill model', 'Color', 'w');
image(imread(fullfile(root, 'results', 'hill_neural_sequence', 'torque_theta_group1.png')));
axis image off;
figure('Name', 'Abdominal torque from the same Hill system', 'Color', 'w');
image(imread(fullfile(root, 'results', 'hill_neural_sequence', 'abdominal_torque_comparison_group1.png')));
axis image off;
