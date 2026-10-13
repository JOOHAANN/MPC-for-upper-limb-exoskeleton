% Compare paired activation sources with the SAME bounded angle model.
warning('BSE:HistoricalSurrogate','Historical angle surrogate, NOT Hill forward dynamics. Use demo_hill_forward for the muscle-driven model.');
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'experiments','bse_angle_priority'));
angleResults = run_angle_priority(root);
% Open one held-out comparison; all speeds are exported as PNG and PDF.
figure('Name','Angle-priority: held-out speed 45','Color','w');
image(imread(fullfile(root,'results','bse_angle_priority','comparison_mode2_speed45.png')));
axis image off;
