"""Expose the existing demo calculations in two readable procedural scripts.

Only code organization changes here. The calibrated Hill equations, controller
algorithm, inputs and plots are reused without numerical substitutions.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEMO = ROOT / "Demo"


def main():
    bse = DEMO / "BSE/run_bse_demo.m"
    text = bse.read_text()
    if "Lumped_model_check_with_coeffs_thesis_combination" in text:
        print("The procedural entry scripts already exist.")
        return
    body = text[text.index("outputFolder = fullfile(folder, 'output');") :]
    body = body.replace(
        "if strcmp(modelVersion,'pearson'), outputFolder = fullfile(outputFolder,'pearson'); end\n",
        "",
    )
    body = body.replace(
        "subjects = bse_subjects(subjectNames,modelVersion);",
        "subjects = bse_subjects(Subject_list);",
    )
    body = body.replace(
        "sourceNames = ",
        "%% Read the subjects and allocate result arrays\nsourceNames = ",
        1,
    )
    body = body.replace(
        "for k = 1:numel(subjects)",
        "%% Run one fixed model per subject and activation source\nfor k = 1:numel(subjects)",
        1,
    )
    body = body.replace(
        "    cases = {};", "    cases = cell(1,9);\n    Case_count = 0;"
    )
    body = body.replace(
        "            plant = bse_model(model, subject.balanceController, sequence.standingAngle, useNative);\n            traces{source} = plant.simulate(inputs, sequence.mass, sequence.waiting, ...\n                sequence.initialState, sequence.dt, 10);",
        """            % The coefficients stay fixed for every condition of this source.
            coeffs = model.parameters;
            neural_coeffs = model.neuralParameters;
            M = coeffs(19);                 % Trunk mass (kg)
            L1 = coeffs(20);                % Trunk centre-of-mass lever (m)
            L2 = coeffs(13);                % External load lever (m)
            R_back = coeffs(21);            % Back muscle moment arm (m)
            R_front = coeffs(22);           % Abdominal moment arm (m)
            Ts = sequence.dt;
            Model_output = Torque_function(model, subject.balanceController, ...
                sequence, inputs, useNative);
            traces{source} = Model_output;""",
    )
    body = body.replace(
        "            cases{end+1} = item; %#ok<AGROW>",
        """            if isfield(segment,'cadenceBPM'), item.cadenceBPM = segment.cadenceBPM; end
            if isfield(segment,'recordId'), item.recordId = segment.recordId; end
            Case_count = Case_count+1;
            cases{Case_count} = item;""",
    )
    body = body.replace(
        "    conditions = cellfun(@(c) [c.loadKg,c.speedCode], cases, 'UniformOutput',false);\n    assert(numel(cases) == 9 && size(unique(vertcat(conditions{:}),'rows'),1) == 9, ...\n        'BSE:Coverage', 'Nine distinct recorded conditions are required.');",
        "    assert(Case_count == 9,'BSE:Coverage','Nine recorded action bouts are required.');",
    )
    body = body.replace(
        "results.metrics = struct2table(rows);",
        "%% Save the traces, correlations and errors\nresults.metrics = struct2table(rows);",
        1,
    )
    split = body.index("\nend\n\nfunction plot_cases")
    body = body[:split] + body[split + 4 :]
    body = body.replace(
        "function plot_cases",
        "%% Plot torque, theta and activation for all nine recorded bouts\nfunction plot_cases",
        1,
    )
    body = body.replace(
        "    heading = sprintf('%g kg | speed code %g',item.loadKg,item.speedCode);",
        """    heading = sprintf('%g kg | speed code %g',item.loadKg,item.speedCode);
    if item.speedCode < 0 && isfield(item,'cadenceBPM')
        heading = sprintf('%g kg | %.2f BPM | bout %g',item.loadKg,item.cadenceBPM,k);
    end""",
    )
    header = """%% Check the lumped BSE model: experimental and trained MAP activation
% Run this script directly, or call run_bse_demo for a selected subject.
% The coupled Hill model produces both torque and theta. Reference signals
% are used only for the plots and error calculations.

%% Choose the subjects and figure display
if ~exist('Subject_list','var')
    Subject_list = {'528','617','725','730'};
end
if ~exist('Show_figures','var')
    Show_figures = true;
end
showFigures = Show_figures;
folder = fileparts(mfilename('fullpath'));
modelVersion = 'current';

%% Load the current calibrated coefficients and prepare the Hill solver
"""
    (DEMO / "BSE/Lumped_model_check_with_coeffs_thesis_combination.m").write_text(
        header + body
    )
    bse.write_text(
        """function results = run_bse_demo(subjectNames, showFigures, modelVersion)
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
Lumped_model_check_with_coeffs_thesis_combination;
end
"""
    )
    nmpc = DEMO / "NMPC/run_nmpc_demo.m"
    text = nmpc.read_text()
    body = text[text.index("folder = fileparts") :]
    body = body.replace(
        "subjects = bse_subjects(subjectNames,'frozen');",
        "subjects = bse_subjects(Subject_list);",
    )
    body = body.replace(
        'results.modelVersion = "frozen";', 'results.modelVersion = "current";'
    )
    body = body.replace(
        "for k = 1:numel(subjects)",
        "%% Simulate each subject with its own controller settings\nfor k = 1:numel(subjects)",
        1,
    )
    body = body.replace(
        "        options = subject.nmpc;",
        """        % Keep one controller setting per subject across all three actions.
        options = subject.nmpc;
        Ts = sequence.dt;
        PH = options.points;             % Prediction points
        CH = options.blocks;             % Control blocks
        u_max = options.capNm;           % Actuator torque limit (N m)""",
        1,
    )
    body = body.replace(
        "        [threshold,limits] = torque_limits(sequence,baseline);",
        """        Torque_Tac0 = baseline(:,3);
        Theta_Tac0 = baseline(:,1);
        [threshold,limits] = torque_limits(sequence,baseline);
        Threshold = threshold;           % 75% of the unassisted back peak""",
    )
    body = body.replace(
        "        actual = runs{source}.actual;",
        """        actual = runs{source}.actual;
        Torque_Tac = actual(:,3);
        Theta_Tac = actual(:,1);
        Tac = actual(:,25);              % Applied assistance torque""",
        1,
    )
    body = body.replace(
        "results.metrics = struct2table(rows);",
        "%% Save the continuous traces and threshold tracking metrics\nresults.metrics = struct2table(rows);",
        1,
    )
    split = body.index("\nend\n\nfunction [threshold, limits]")
    body = body[:split] + body[split + 4 :]
    header = """%% Offline NMPC simulation with trained MAP activation
% Run this script directly, or call run_nmpc_demo for selected subjects.
% All plant states remain continuous through 5, 10 and 15 kg actions and waits.

%% Choose the subjects and figure display
if ~exist('Subject_list','var')
    Subject_list = {'528','617','725','730'};
end
if ~exist('Show_figures','var')
    Show_figures = true;
end
showFigures = Show_figures;

%% Load the current BSE model and subject-specific NMPC parameters
"""
    (
        DEMO / "NMPC/offline_simulator_NMPC_quasimodo_incremental_MPC_with_MAP.m"
    ).write_text(header + body)
    nmpc.write_text("""function results = run_nmpc_demo(subjectNames, showFigures)
%RUN_NMPC_DEMO Run the readable continuous NMPC script for four subjects.
%   run_nmpc_demo
%   results = run_nmpc_demo("617",false);
if nargin < 1 || isempty(subjectNames)
    subjectNames = {'528','617','725','730'};
end
if nargin < 2
    showFigures = true;
end
Subject_list = subjectNames;
Show_figures = showFigures;
offline_simulator_NMPC_quasimodo_incremental_MPC_with_MAP;
end
""")
    print(
        "Installed the two procedural scripts; run_bse_demo and run_nmpc_demo remain available."
    )


if __name__ == "__main__":
    main()
