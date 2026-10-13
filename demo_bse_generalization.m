% DEMO_BSE_GENERALIZATION Audit BSE targets and test torque estimators.
% This is an offline, grouped-validation experiment, not an NMPC controller.
% The original BSE source, parameters, and input assets are never overwritten.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot, fullfile(projectRoot, 'experiments', 'bse_generalization'));
bseResults = run_bse_generalization(projectRoot);
