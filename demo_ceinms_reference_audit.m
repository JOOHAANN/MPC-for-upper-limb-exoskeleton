function demo_ceinms_reference_audit()
% Display torque recalculation and its limits; no independent calibration.
root = fileparts(mfilename('fullpath'));
directory = fullfile(root, 'docs', 'ceinms_reference_audit');
file = fullfile(directory, 'reference_torque_diagnostics.png');
assert(isfile(file), 'Run experiments/multi_subject/audit_ceinms_reference.py first.');
figure('Name', 'CEINMS reference audit: uncalibrated candidates', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [60 40 1400 1000]);
ax = axes('Position', [0.005 0.005 0.99 0.99]);
image(ax, imread(file));
axis(ax, 'image');
axis(ax, 'off');
fprintf('No independent torque/force target is available for these records.\n');
fprintf('Candidates are sensitivity estimates, not calibrated ground truth.\n');
fprintf('Audit: %s\n', fullfile(directory, 'audit.json'));
fprintf('Report: %s\n', fullfile(root, 'docs', 'CEINMS_REFERENCE_AUDIT.zh-CN.md'));
end
