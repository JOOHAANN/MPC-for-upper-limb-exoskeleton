%% Check the lumped BSE model: experimental and trained MAP activation
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
outputFolder = fullfile(folder, 'output');
if ~isfolder(outputFolder), mkdir(outputFolder); end
useNative = build_bse_solver;
subjects = bse_subjects(Subject_list);
%% Read the subjects and allocate result arrays
sourceNames = ["Experimental activation", "MAP activation"];
rows = struct([]);
results.subjects = cell(size(subjects));
results.modelVersion = string(modelVersion);

%% Run one fixed model per subject and activation source
for k = 1:numel(subjects)
    subject = subjects{k};
    fprintf('\nBSE | %s | nine conditions, two fixed parameter vectors\n', subject.id);
    if isfield(subject,'calibrationTargetsPassed') && ~all(subject.calibrationTargetsPassed)
        warning('BSE:UnmetTargets', ...
            '%s includes a nonpassing calibration candidate. See data/fit_summary.csv.',subject.id);
    end
    cases = cell(1,9);
    Case_count = 0;
    for group = 1:numel(subject.sequences)
        sequence = subject.sequences{group};
        traces = cell(1, 2);
        for source = 1:2
            model = subject.models{source};
            % The coefficients stay fixed for every condition of this source.
            coeffs = model.parameters;
            neural_coeffs = model.neuralParameters;
            M = coeffs(19);                 % Trunk mass (kg)
            L1 = coeffs(20);                % Trunk centre-of-mass lever (m)
            L2 = coeffs(13);                % External load lever (m)
            R_back = coeffs(21);            % Back muscle moment arm (m)
            R_front = coeffs(22);           % Abdominal moment arm (m)
            coeffs([19 20 13 21 22]) = [M L1 L2 R_back R_front];
            model.parameters = coeffs;
            model.neuralParameters = neural_coeffs;
            inputs = bse_inputs(sequence, model, source);
            Ts = sequence.dt;
            Model_output = Torque_function(model, subject.balanceController, ...
                sequence, inputs, useNative);
            traces{source} = Model_output;
        end
        for stage = 1:numel(sequence.segments)
            segment = sequence.segments{stage};
            index = segment.startIndex:segment.endIndex;
            item = struct('loadKg',segment.loadKg, 'speedCode',segment.speedCode, ...
                'time',sequence.time(index)-sequence.time(index(1)), ...
                'referenceTheta',sequence.referenceTheta(index), ...
                'referenceTorque',sequence.referenceTorque(index), ...
                'backActivation',sequence.backActivation(index,:), ...
                'outputs',cat(3,traces{1}(index,:),traces{2}(index,:)));
            item.referenceTorque(~sequence.torqueValid(index)) = NaN;
            if isfield(segment,'cadenceBPM'), item.cadenceBPM = segment.cadenceBPM; end
            if isfield(segment,'recordId'), item.recordId = segment.recordId; end
            Case_count = Case_count+1;
            cases{Case_count} = item;
            for source = 1:2
                out = item.outputs(:,:,source);
                valid = isfinite(item.referenceTorque);
                correlation = corrcoef(out(:,1), item.referenceTheta);
                torqueCorrelation = corrcoef(out(valid,3), item.referenceTorque(valid));
                row = struct('Subject',string(subject.id), 'Source',sourceNames(source), ...
                    'LoadKg',segment.loadKg, 'SpeedCode',segment.speedCode, 'CaseNumber',Case_count, ...
                    'ThetaCorrelation',correlation(1,2), ...
                    'ThetaRMSEDeg',rad2deg(sqrt(mean((out(:,1)-item.referenceTheta).^2))), ...
                    'TorqueRMSENm',sqrt(mean((out(valid,3)-item.referenceTorque(valid)).^2)), ...
                    'TorqueCorrelation',torqueCorrelation(1,2), ...
                    'TorqueReferenceFraction',mean(valid), ...
                    'BoundaryPeakNm',max(abs(out(:,26))), ...
                    'FiberResidual',max(abs(out(:,16))), ...
                    'MomentumResidualNm',max(abs(out(:,17))));
                rows = [rows; row]; %#ok<AGROW>
            end
        end
    end
    assert(Case_count == 9,'BSE:Coverage','Nine recorded action bouts are required.');
    plot_cases(subject, cases, 'torque', outputFolder, showFigures);
    plot_cases(subject, cases, 'theta', outputFolder, showFigures);
    plot_cases(subject, cases, 'activation', outputFolder, showFigures);
    results.subjects{k} = struct('id',subject.id, 'models',{subject.models}, 'cases',{cases});
    fprintf('Completed %s: 9/9 conditions.\n', subject.id);
end
%% Save the traces, correlations and errors
results.metrics = struct2table(rows);
writetable(results.metrics, fullfile(outputFolder, 'bse_metrics.csv'));
save(fullfile(outputFolder, 'bse_results.mat'), 'results', '-v7');
disp(results.metrics(:,{'Subject','Source','LoadKg','SpeedCode','TorqueCorrelation','ThetaCorrelation','TorqueRMSENm'}));
fprintf('\nBSE figures and numerical results: %s\n', outputFolder);

%% Plot torque, theta and activation for all nine recorded bouts
function plot_cases(subject, cases, signal, folder, showFigures)
visible = 'off';
if showFigures, visible = 'on'; end
fig = figure('Name',sprintf('%s | BSE %s',subject.id,signal), ...
    'Color','w', 'Position',[80 80 1440 960], 'Visible',visible);
layout = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
colors = [0 .447 .741; .494 .184 .556];
for k = 1:9
    item = cases{k};
    ax = nexttile(layout);
    hold(ax,'on');
    switch signal
        case 'torque'
            plot(ax,item.time,item.referenceTorque,'-','Color',[.85 .2 .15],'LineWidth',1.5);
            for source = 1:2
                plot(ax,item.time,item.outputs(:,3,source),'Color',colors(source,:),'LineWidth',1.4);
            end
            ylabel(ax,'Net muscle torque (N m)');
            shade_missing_reference(ax,item.time,isfinite(item.referenceTorque));
        case 'theta'
            plot(ax,item.time,item.referenceTheta,'--','Color',[.85 .2 .15],'LineWidth',1.5);
            for source = 1:2
                plot(ax,item.time,item.outputs(:,1,source),'Color',colors(source,:),'LineWidth',1.4);
            end
            ylabel(ax,'Theta (rad)');
        case 'activation'
            for source = 1:2
                plot(ax,item.time,item.backActivation(:,source),'Color',colors(source,:),'LineWidth',1.4);
            end
            ylabel(ax,'Back activation');
    end
    heading = sprintf('%g kg | speed code %g',item.loadKg,item.speedCode);
    if item.speedCode < 0 && isfield(item,'cadenceBPM')
        heading = sprintf('%g kg | %.2f BPM | bout %g',item.loadKg,item.cadenceBPM,k);
    end
    if strcmp(signal,'torque')
        heading = sprintf('%s | reference support %.0f%%',heading, ...
            100*mean(isfinite(item.referenceTorque)));
    end
    title(ax,heading);
    xlabel(ax,'Time within recorded action (s)');
    grid(ax,'on'); box(ax,'on');
    ax.FontSize = 10;
    xlim(ax,[item.time(1),item.time(end)]);
end
if strcmp(signal,'activation')
    labels = {'Experimental activation','Trained MAP activation'};
else
    if strcmp(signal,'torque')
        referenceLabel = 'CEINMS estimate';
        if endsWith(subject.id,'A01'), referenceLabel = 'CEINMS estimate, A=0.1'; end
        labels = {referenceLabel,'BSE torque: experimental input','BSE torque: MAP input'};
    else
        labels = {'Original IK reference','BSE theta: experimental input','BSE theta: MAP input'};
    end
end
key = legend(ax,labels,'Orientation','horizontal');
key.Layout.Tile = 'south';
title(layout,sprintf('%s | %s BSE %s | fixed parameters across all nine conditions',subject.id,subject.modelVersion,signal), ...
    'Interpreter','none','FontSize',15);
if strcmp(signal,'torque')
    subtitle(layout,'Net = back minus abdomen. Gray regions have no supported CEINMS estimate; solid red is available reference.');
elseif strcmp(signal,'theta')
    subtitle(layout,'Original angle scale; smaller angle means deeper flexion.');
else
    subtitle(layout,'Prepared experimental activation and frozen trained MAP output; no fitting in this demo.');
end
if isfield(subject,'calibrationTargetsPassed') && ~all(subject.calibrationTargetsPassed)
    subtitle(layout,'Experimental candidate: mean torque r = 0.693, theta r = 0.653; requested targets remain unmet.');
end
name = fullfile(folder,sprintf('%s_bse_%s',subject.shortName,signal));
exportgraphics(fig,[name '.png'],'Resolution',160);
savefig(fig,[name '.fig']);
if ~showFigures, close(fig); end
end

function shade_missing_reference(ax,time,valid)
% A line style must not be confused with a missing acquisition interval.
edges = diff([false; ~valid(:); false]);
first = find(edges == 1);
last = find(edges == -1)-1;
limits = ylim(ax);
for k = 1:numel(first)
    left = time(first(k));
    right = time(min(last(k)+1,numel(time)));
    area = patch(ax,[left right right left],limits([1 1 2 2]),[.6 .6 .6], ...
        'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');
    uistack(area,'bottom');
end
ylim(ax,limits);
end
