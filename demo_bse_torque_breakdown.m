% Diagnose torque components without changing or refitting either model.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'experiments','bse_angle_priority'));
torqueDiagnostics=diagnose_angle_torque(root);
figure('Name','Boundary torque and torque components','Color','w');
image(imread(fullfile(root,'results','bse_torque_breakdown','components_speed45.png')));
axis image off;
