% DEMO_BSE_TORQUE_THETA Compare archived BSE and a new two-output surrogate.
% Unlike the earlier torque estimator, this experiment rolls theta forward.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'experiments', 'bse_generalization'));
jointResults = compare_bse_torque_theta(root);
