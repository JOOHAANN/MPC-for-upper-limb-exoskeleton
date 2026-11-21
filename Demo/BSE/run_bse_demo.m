function results = run_bse_demo(subjectNames, showFigures, modelVersion)
%RUN_BSE_DEMO Run the readable BSE script for one or all four subjects.
%   run_bse_demo
%   results = run_bse_demo("617",false);
if nargin < 1 || isempty(subjectNames)
    subjectNames = {'528','617','725','730'};
end
if nargin < 2
    showFigures = true;
end
if nargin >= 3
    assert(any(strcmp(modelVersion,{'current','pearson'})), ...
        'BSE:Version','The Demo contains only the current calibrated parameters.');
end
Subject_list = subjectNames;
Show_figures = showFigures;
scriptFile = fullfile(fileparts(mfilename('fullpath')), ...
    'Lumped_model_check_with_coeffs_thesis_combination.m');
run(scriptFile);
end
