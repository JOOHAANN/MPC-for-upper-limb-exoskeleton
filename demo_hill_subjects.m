function demo_hill_subjects(mode)
% View or reproduce personalized MAP, BSE/Hill and MPC subject results.
% demo_hill_subjects          - display the latest quality-gated figures.
% demo_hill_subjects('archive')- view the historical, less strictly gated run.
% demo_hill_subjects('replay') - replay frozen models on the local recordings.
% demo_hill_subjects('train')  - recalibrate MAP/BSE, tune MPC and validate.
% demo_hill_subjects('tune')   - retune only MPC for minimal assistance.
% demo_hill_subjects('threshold-replay') - replay the saved near-threshold MPC.
if nargin < 1, mode = 'plots'; end
mode = validatestring(mode, {'plots', 'archive', 'replay', 'train', 'tune', 'threshold-replay'});
root = fileparts(mfilename('fullpath'));
if any(strcmp(mode, {'replay', 'train', 'tune', 'threshold-replay'}))
    addpath(fullfile(root, 'experiments', 'legacy_mpc'));
    pythonExecutable = find_hill_mpc_python();
    backend = fullfile(root, 'experiments', 'multi_subject', 'run_subjects.py');
    if any(strcmp(mode, {'tune', 'threshold-replay'}))
        backend = fullfile(root, 'experiments', 'multi_subject', 'run_threshold.py');
    end
    command = sprintf('"%s" "%s"', pythonExecutable, backend);
    if strcmp(mode, 'train'), command = [command ' --retrain']; end
    if strcmp(mode, 'tune'), command = [command ' --tune']; end
    [status, output] = system(command, '-echo');
    assert(status == 0, 'Subject pipeline failed: %s', output);
end
directory = fullfile(root, 'docs', 'multi_subject_evaluation');
qualityDirectory = fullfile(root, 'docs', 'subject_quality_v2');
qualityRun = ~any(strcmp(mode, {'archive','replay'})) && ...
    isfile(fullfile(qualityDirectory, 'multi_subject_overview.png'));
if qualityRun, directory = qualityDirectory; end
trackingDirectory = fullfile(root, 'docs', 'mpc_threshold_tracking');
trackingFile = dir(fullfile(trackingDirectory, 'multi_subject_overview.png'));
qualityFile = dir(fullfile(qualityDirectory, 'multi_subject_overview.png'));
trackingRun = ~any(strcmp(mode, {'archive','replay'})) && ~isempty(trackingFile) && ...
    (isempty(qualityFile) || trackingFile.datenum >= qualityFile.datenum);
if trackingRun
    directory = trackingDirectory;
    qualityRun = true;
end
assert(isfile(fullfile(directory, 'multi_subject_overview.png')), ...
    'Published figures are missing. Run demo_hill_subjects(''replay'').');
window = figure('Name', 'Personalized MAP, BSE/Hill and MPC', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [60, 40, 1400, 950]);
tabs = uitabgroup(window);
overview = uitab(tabs, 'Title', 'Overview');
showImage(overview, fullfile(directory, 'multi_subject_overview.png'));
if qualityRun && isfile(fullfile(directory, 'quality_overview.png'))
    qualityTab = uitab(tabs, 'Title', 'Data quality first');
    showImage(qualityTab, fullfile(directory, 'quality_overview.png'));
end
if qualityRun && ~trackingRun && isfile(fullfile(directory, 'boundary_audit.png'))
    boundaryTab = uitab(tabs, 'Title', 'Soft boundary audit');
    showImage(boundaryTab, fullfile(directory, 'boundary_audit.png'));
end
subjects = {'240617Pilot', '240627Pilot', '240719Pilot', '240723Pilot', '240725Pilot'};
slugs = {'before_after', 'continuous', 'activation', 'cases_experimental', 'cases_map'};
titles = {'Controller before / after', 'Continuous actions', 'Activation inputs', 'Experimental cases', 'MAP cases'};
for index = 1:numel(subjects)
    subjectTab = uitab(tabs, 'Title', subjects{index});
    pages = uitabgroup(subjectTab);
    if qualityRun && ~isfile(fullfile(directory, [subjects{index} '_continuous.png']))
        tab = uitab(pages, 'Title', 'Not calibrated');
        uicontrol('Parent', tab, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.08 0.35 0.84 0.25], 'BackgroundColor', 'white', ...
            'FontSize', 16, 'String', sprintf('%s\nQuality / coverage gate not met.\nSee data_quality_audit_v2.pdf and QUALITY_FIRST.zh-CN.md.', subjects{index}));
        continue;
    end
    for page = 1:numel(slugs)
        file = fullfile(directory, [subjects{index} '_' slugs{page} '.png']);
        if isfile(file)
            tab = uitab(pages, 'Title', titles{page});
            showImage(tab, file);
        end
    end
end
if trackingRun
    fprintf('Full report: %s\n', fullfile(directory, 'threshold_tracking_results.pdf'));
    fprintf('Controller method: %s\n', fullfile(root, 'docs', 'MPC_THRESHOLD_TRACKING.zh-CN.md'));
elseif qualityRun
    fprintf('Full report: %s\n', fullfile(directory, 'quality_first_results.pdf'));
    fprintf('Data quality and method: %s\n', fullfile(root, 'docs', 'QUALITY_FIRST.zh-CN.md'));
else
    fprintf('Full report: %s\n', fullfile(directory, 'multi_subject_results.pdf'));
    fprintf('Data quality and method: %s\n', fullfile(root, 'docs', 'MULTI_SUBJECT.zh-CN.md'));
end
end

function showImage(parent, file)
axesHandle = axes('Parent', parent, 'Position', [0.005, 0.005, 0.99, 0.99]);
image(axesHandle, imread(file));
axis(axesHandle, 'image');
axis(axesHandle, 'off');
end
