% Two compliant Hill units drive both torque and angle from one force balance.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'hill_compliant'));
compliantResults = run_hill_compliant(root);
% Expose the parameters used in calibration and each held-speed experiment.
p = compliantResults.models{1}.parameters;
fprintf('\nAll-trial shared geometry: M=%.4f kg, L1=%.5f m, L2=%.5f m, phi=%.5f rad\n', p(19), p(20), p(13), p(14));
fprintf('All-trial soft limits: k=%.1f Nm/rad, cap=%.1f Nm, onset=[%.3f %.3f] rad\n', p(29), p(30), p(31), p(32));
for fold = 1:3
    p = compliantResults.validation.models{1, fold}.parameters;
    fprintf('Held speed %g: M=%.3f, L1=%.4f, L2=%.4f, phi=%.4f; soft k=%.1f, cap=%.1f, onset=[%.3f %.3f]\n', ...
            compliantResults.validation.heldSpeeds(fold), p(19), p(20), p(13), p(14), p(29), p(30), p(31), p(32));
end
figure('Name', 'Compliant Hill muscle comparison', 'Color', 'w');
image(imread(fullfile(root, 'results', 'hill_compliant', 'validation', 'comparison_speed45.png')));
axis image off;
