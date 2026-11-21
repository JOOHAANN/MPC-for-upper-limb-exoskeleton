% Continuous original 5/10/15-kg actions, with muscle-driven standing waits.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'hill_compliant'));
options = struct();
outputDirectory = 'hill_sequence';
if isfile(fullfile(root, 'data', 'processed', 'hill_compliant', 'neural_models.mat'))
    % Prefer the screened fixed-parameter neural/Hill model; keep the old
    % calibration's outputs available as an explicit historical baseline.
    options.modelFile = 'neural_models.mat';
    options.outputDirectory = 'hill_neural_sequence';
    outputDirectory = options.outputDirectory;
end
sequenceResults = run_hill_sequence(root, options);
figure('Name', 'Continuous 5 / 10 / 15 kg Hill sequence', 'Color', 'w');
image(imread(fullfile(root, 'results', outputDirectory, 'torque_theta_group1.png')));
axis image off;
figure('Name', 'Abdominal torque: experimental versus MAP activation', 'Color', 'w');
image(imread(fullfile(root, 'results', outputDirectory, 'abdominal_torque_comparison_group1.png')));
axis image off;
