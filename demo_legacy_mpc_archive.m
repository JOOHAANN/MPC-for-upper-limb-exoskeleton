% Diagnostic reproduction of the archived BSE/MPC settings; known unstable.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'legacy_mpc'));
legacyMPCReports = run_legacy_mpc_tests(root);
legacyMPCTable = struct2table(legacyMPCReports);
disp(legacyMPCTable(:, {'source', 'completed', 'lastStep', 'elapsedSeconds'}));
legacyMPCSummary = plot_legacy_mpc_results(root, legacyMPCReports);
