function window = demo_nine_case_75pct(subject)
% View saved BSE/NMPC results; this does not retrain or change a model.
% Examples: demo_nine_case_75pct; demo_nine_case_75pct('725');
if nargin < 1, subject = 'all'; end
root = fileparts(mfilename('fullpath'));
folder = fullfile(root, 'docs', 'nine_case_75pct');
manifestPath = fullfile(folder, 'report_manifest.json');
assert(isfile(manifestPath), 'The nine-case report has not been generated.');
manifest = jsondecode(fileread(manifestPath));
files = cellstr(string(manifest.figures));
subject = char(string(subject));
keep = contains(files, 'peak_tracking_overview');
if strcmpi(subject, 'all')
    keep(:) = true;
else
    assert(any(strcmp(subject, {'528','725','730','240528Pilot', ...
        '240725Exp_A01','240730Exp_A01'})), 'Choose 528, 725, 730, or all.');
    keep = keep | contains(files, subject);
end
files = files(keep);
window = figure('Name', 'Nine conditions: fixed BSE and 75% peak NMPC', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', ...
    'Position', [.025 .04 .95 .89]);
tabs = uitabgroup(window);
for index = 1:numel(files)
    [~,name,extension] = fileparts(files{index});
    filePath = fullfile(folder, [name extension]);
    assert(isfile(filePath), 'Missing published figure: %s', filePath);
    label = strrep(name, '240725Exp_A01', '725');
    label = strrep(label, '240730Exp_A01', '730');
    label = strrep(label, '240528Pilot', '528');
    panel = uitab(tabs, 'Title', strrep(label, '_', ' '));
    ax = axes('Parent', panel, 'Position', [.005 .005 .99 .99]);
    image(ax, imread(filePath));
    axis(ax, 'image'); axis(ax, 'off');
end
fprintf('Black dashed line: 75%% of the unassisted back-muscle peak.\n');
fprintf('Green band: 0-5 N m below that threshold.\n');
fprintf('Full report: %s\n', fullfile(root, 'output', 'pdf', 'nine_case_bse_nmpc_75pct.pdf'));
end
