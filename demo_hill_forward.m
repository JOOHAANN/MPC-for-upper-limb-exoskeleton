% Genuine muscle-driven forward dynamics: no reference tracking or angle target.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'experiments','hill_forward'));
hillResults=run_hill_forward(root);
fprintf('Development baseline: W-shape and range requirements are NOT yet satisfied. See docs/HILL_FORWARD.zh-CN.md.\n');
figure('Name','Hill muscle forward dynamics','Color','w');
image(imread(fullfile(root,'results','hill_forward','heldout_speed45.png'))); axis image off;
