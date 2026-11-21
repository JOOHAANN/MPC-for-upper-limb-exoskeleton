% Plot the current soft-limit law and its actual torque in all nine cases.
% This uses the selected model without fitting or changing any parameters.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'hill_compliant'));
softBoundarySummary = plot_hill_soft_boundary(root);
disp(softBoundarySummary);
