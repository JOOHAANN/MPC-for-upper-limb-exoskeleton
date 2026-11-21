% Audit excessive inferred abdominal recruitment and compare it with back activity.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'experiments', 'hill_compliant'));
recruitmentAudit = run_hill_recruitment_audit(root);
figure('Name', 'Back versus abdominal recruitment audit', 'Color', 'w');
image(imread(fullfile(root, 'results', 'hill_recruitment_audit', 'recruitment_audit_group1.png')));
axis image off;
