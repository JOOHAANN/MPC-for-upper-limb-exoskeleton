function subjects = bse_subjects(names, modelVersion)
%BSE_SUBJECTS Load one fixed model per subject and activation source.
% The demonstration contains only the current four-subject parameter set.
if nargin < 2, modelVersion = 'current'; end
if nargin < 1 || isempty(names)
    names = ["528", "617", "725", "730"];
end
names = string(names);
catalog = ["528", "617", "725", "730"];
fullNames = ["240528Pilot", "240617Pilot", "240725Exp_A01", "240730Exp_A01"];
folder = fullfile(fileparts(mfilename('fullpath')), 'data');
assert(any(strcmp(modelVersion,{'current','pearson'})), ...
    'BSE:Version','The Demo contains only the current calibrated parameters.');
subjects = cell(1, numel(names));
for k = 1:numel(names)
    index = find(names(k) == catalog | names(k) == fullNames, 1);
    assert(~isempty(index), 'BSE:Subject', 'Choose 528, 617, 725 or 730.');
    cache = load(fullfile(folder, "subject_"+catalog(index)+".mat"), 'subject');
    subjects{k} = cache.subject;
    subjects{k}.modelVersion = "current";
end
end
