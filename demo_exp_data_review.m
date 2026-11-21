function demo_exp_data_review()
% Display the Exp data audit; this does not fit a model or replay MPC.
% One complete M-shaped pair is retained per condition where available.
root = fileparts(mfilename('fullpath'));
directory = fullfile(root, 'docs', 'exp_one_cycle_review');
files = {'240730Exp_EMG_activation_levels.png', '240725Exp_keep_first_M.png', ...
    '240725Exp_one_M_grid.png', '240730Exp_one_M_grid.png'};
titles = {'EMG vs activation levels', 'Keep one complete M', ...
    '240725Exp conditions', '240730Exp conditions'};
assert(all(cellfun(@(f) isfile(fullfile(directory, f)), files)), ...
    ['Review figures are missing. Run experiments/multi_subject/review_exp_one_cycle.py ' ...
    'and experiments/multi_subject/compare_exp_signal_levels.py.']);
window = figure('Name', 'Exp data review: existing IK retained', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [70, 50, 1350, 950]);
tabs = uitabgroup(window);
for index = 1:numel(files)
    tab = uitab(tabs, 'Title', titles{index});
    ax = axes('Parent', tab, 'Position', [0.005, 0.005, 0.99, 0.99]);
    image(ax, imread(fullfile(directory, files{index})));
    axis(ax, 'image');
    axis(ax, 'off');
end
fprintf('Data decisions: %s\n', fullfile(directory, 'conditions.csv'));
fprintf('Runtime explanation: %s\n', fullfile(root, 'docs', 'PROJECT_RUNTIME.zh-CN.md'));
end
