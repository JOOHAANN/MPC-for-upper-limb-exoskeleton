function demo_exp_activation_A01()
% View both Exp groups after recomputing activation with shapeFactor = +0.1.
% This viewer does not retrain MAP/BSE or recompute torque and MPC.
root = fileparts(mfilename('fullpath'));
directory = fullfile(root, 'docs', 'exp_activation_A01');
files = {'activation_before_after.png', '240725Exp_A01_grid.png', ...
    '240730Exp_A01_grid.png'};
titles = {'Back and abdomen: before / after', '240725Exp: A = 0.1', ...
    '240730Exp: A = 0.1'};
assert(all(cellfun(@(f) isfile(fullfile(directory, f)), files)), ...
    'Run experiments/multi_subject/recompute_exp_activation.py to create the figures.');
window = figure('Name', 'Exp activation: A = +0.1', 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', [70 50 1350 950]);
tabs = uitabgroup(window);
for index = 1:numel(files)
    tab = uitab(tabs, 'Title', titles{index});
    ax = axes('Parent', tab, 'Position', [0.005 0.005 0.99 0.99]);
    image(ax, imread(fullfile(directory, files{index})));
    axis(ax, 'image');
    axis(ax, 'off');
end
fprintf('A = 0.1 input profile: %s\n', ...
    fullfile(root, 'data', 'processed', 'exp_activation_A01', 'profile.json'));
fprintf('Report: %s\n', fullfile(root, 'docs', 'EXP_ACTIVATION_A01.zh-CN.md'));
end
