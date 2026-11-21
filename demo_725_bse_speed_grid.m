function demo_725_bse_speed_grid(mode)
% Show 725 BSE cases with CEINMS support; five unsupported cases are excluded.
% This viewer does not reuse older MPC plots as results of the new BSE model.
% demo_725_bse_speed_grid        - show saved results without training.
% demo_725_bse_speed_grid('run') - finish a missing grid fit, skip it if saved.
if nargin < 1, mode = 'plots'; end
mode = validatestring(mode, {'plots','run'});
root = fileparts(mfilename('fullpath'));
if strcmp(mode, 'run')
    addpath(fullfile(root, 'experiments', 'legacy_mpc'));
    executable = find_hill_mpc_python();
    backend = fullfile(root, 'experiments', 'multi_subject', 'run_exp_speed_grid.py');
    [status, output] = system(sprintf('"%s" "%s" --fit', executable, backend), '-echo');
    assert(status == 0, '725 grid workflow failed: %s', output);
end
directory = fullfile(root, 'docs', 'exp_speed_grid');
assert(isfile(fullfile(directory,'evaluation.json')), ...
    'Saved grid results missing. Run demo_725_bse_speed_grid(''run'').');
window = figure('Name','725 BSE | 10 retained conditions | no MPC assistance', ...
    'NumberTitle','off','Color','w','Position',[30 40 1450 900]);
tabs = uitabgroup(window);
files = {'trend_summary.png','reference_grid.png','activation_timing_grid.png', ...
    'experimental_theta.png','map_theta.png','experimental_torque.png', ...
    'map_torque.png','experimental_muscles.png','map_muscles.png','soft_boundary.png'};
labels = {'Trend overview','Reference audit','Activation timing','Experimental theta','MAP theta', ...
    'Experimental torque','MAP torque','Experimental back / abdomen', ...
    'MAP back / abdomen','Soft boundary torque'};
for index = 1:numel(files)
    path = fullfile(directory, files{index});
    if ~isfile(path), continue; end
    tab = uitab(tabs,'Title',labels{index});
    ax = axes('Parent',tab,'Position',[0.005 0.005 0.99 0.99]);
    image(ax,imread(path));
    axis(ax,'image'); axis(ax,'off');
end
fprintf('10 retained conditions are calibration; 5 without torque are excluded.\n');
fprintf('Independent generalization is not tested in this grid.\n');
fprintf('Grid metrics: %s\n',fullfile(directory,'evaluation.json'));
fprintf('Method: %s\n',fullfile(root,'docs','EXP_SPEED_GRID.zh-CN.md'));
end
