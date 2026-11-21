function demo_exp_subjects(mode)
% Show newly fitted Exp results, or complete only missing training stages.
% demo_exp_subjects          - show saved plots; does not fit any model.
% demo_exp_subjects('run')   - skip completed subjects and run missing stages.
% demo_exp_subjects('status')- print completion status without changing models.
if nargin < 1, mode = 'plots'; end
mode = validatestring(mode, {'plots','run','status'});
root = fileparts(mfilename('fullpath'));
if any(strcmp(mode, {'run','status'}))
    addpath(fullfile(root, 'experiments', 'legacy_mpc'));
    executable = find_hill_mpc_python();
    backend = fullfile(root, 'experiments', 'multi_subject', 'run_exp_subjects.py');
    command = sprintf('"%s" "%s"', executable, backend);
    if strcmp(mode,'status'), command = [command ' --status']; end
    [status, output] = system(command, '-echo');
    assert(status == 0, 'Exp subject pipeline failed: %s', output);
    if strcmp(mode,'status'), return; end
end
directory = fullfile(root, 'docs', 'exp_subjects_A01');
threeLoadDirectory = fullfile(root, 'docs', 'exp_three_loads');
threeLoadRun = isfile(fullfile(threeLoadDirectory,'completion_status.json'));
if threeLoadRun, directory = threeLoadDirectory; end
assert(isfile(fullfile(directory, 'completion_status.json')), ...
    'Results missing. Run demo_exp_subjects(''run'').');
window = figure('Name','Exp MAP, BSE and NMPC results','NumberTitle','off', ...
    'Color','w','Position',[50 30 1420 950]);
tabs = uitabgroup(window);
files = {'completion_overview.png', ...
    '240725Exp_A01_map.png','240730Exp_A01_map.png', ...
    '240725Exp_A01_exclusion_evidence.png', ...
    '240730Exp_A01_continuous.png', ...
    '240730Exp_A01_cases_experimental.png','240730Exp_A01_cases_map.png'};
labels = {'Completion / quality','725 MAP','730 MAP','725 data limitations', ...
    '730 BSE + NMPC','730 experimental cases','730 MAP cases'};
if threeLoadRun
    files = {'completion_overview.png','240725Exp_A01_map.png', ...
        '240730Exp_A01_map.png','three_anchor_inputs.png', ...
        '240725Exp_A01_calibration.png', ...
        '240725Exp_A01_continuous.png','240725Exp_A01_cases_experimental.png', ...
        '240725Exp_A01_cases_map.png','240730Exp_A01_continuous.png', ...
        '240730Exp_A01_cases_experimental.png','240730Exp_A01_cases_map.png'};
    labels = {'Three-load status','725 MAP','730 MAP','725 5/10/15 inputs', ...
        '725 three fitted actions', ...
        '725 BSE + NMPC','725 experimental cases','725 MAP cases', ...
        '730 BSE + NMPC','730 experimental cases','730 MAP cases'};
end
for index = 1:numel(files)
    path = fullfile(directory, files{index});
    if ~isfile(path), continue; end
    tab = uitab(tabs, 'Title', labels{index});
    ax = axes('Parent',tab,'Position',[0.005 0.005 0.99 0.99]);
    image(ax,imread(path));
    axis(ax,'image'); axis(ax,'off');
end
fprintf('Subject status: %s\n',fullfile(directory,'completion_status.json'));
if threeLoadRun, method = 'EXP_THREE_LOADS.zh-CN.md'; else, method = 'EXP_SUBJECTS_A01.zh-CN.md'; end
fprintf('Method and limitations: %s\n',fullfile(root,'docs',method));
end
