clc
clearvars -except Subject_list Show_figures Use_native Output_label
close all

%% Direct-activation BSE: two MTUs, standing PD only, no angle limits
% One fixed parameter vector per subject and activation source, across nine
% recorded conditions. Experimental and frozen trained-MAP inputs are used.
% No external assistance, movement PD, decoded angle or boundary torque.
if ~exist('Subject_list','var'),Subject_list={'528'};end
if ~exist('Show_figures','var'),Show_figures=true;end
if ~exist('Use_native','var'),Use_native=true;end
if ~exist('Output_label','var'),Output_label='';end
folder=fileparts(mfilename('fullpath'));addpath(folder,fileparts(fileparts(folder)));
outputFolder=fullfile(folder,'output',Output_label);
if ~isfolder(outputFolder),mkdir(outputFolder);end
sourceNames={'Experimental activation','MAP activation'};
results=struct;results.subjects=cell(1,numel(Subject_list));rows=struct([]);
parameterRows=struct([]);

%% Read each subject and keep coefficients fixed throughout its sequences
for subjectIndex=1:numel(Subject_list)
    data=load(fullfile(folder,'data',['subject_' char(Subject_list{subjectIndex}) '.mat']));
    subject=data.subject;cases=cell(1,9);caseCount=0;groups=cell(size(subject.sequences));
    firstCaseRow=numel(rows)+1;
    assert(isfield(subject.calibrationProtocol,'status') && strcmp(subject.calibrationProtocol.status,'calibrated'), ...
        'The input revision must be physically recalibrated before running this demo.');
    fprintf('\nDirect dual-MTU BSE | %s\n',subject.id);
    for source=1:2
        parameterRow=struct('Subject',string(subject.shortName), ...
            'Source',string(sourceNames{source}),'Parameters',{subject.models{source}.parameters}, ...
            'AbdomenScale',subject.models{source}.abdomenScale);
        parameterRows=[parameterRows;parameterRow]; %#ok<AGROW>
    end
    for group=1:numel(subject.sequences)
        q=subject.sequences{group};runs=cell(1,2);
        for source=1:2
            model=subject.models{source};
            inputs=prepare_activation(q,model,source);
            plant=dual_MTU_model(model.parameters,subject.balanceController,q.standingAngle,Use_native);
            runs{source}=simulate_sequence(plant,q,inputs);
            verify_feedback(runs{source},q);
            fprintf('  Sequence %d | %s | %d/%d samples\n',group,sourceNames{source}, ...
                runs{source}.completedSamples,numel(q.time));
        end
        groups{group}=struct('sequence',q,'runs',{runs});
        for stage=1:numel(q.segments)
            seg=q.segments{stage};caseCount=caseCount+1;
            ix=seg.startIndex:seg.endIndex;
            assert(all(q.waiting(ix)==0),'A recorded movement contains standing feedback.');
            % Include real neighboring samples to show standing-PD switching.
            context=max(1,seg.startIndex-round(.5/q.dt)):min(numel(q.time),seg.endIndex+round(1/q.dt));
            item=struct('loadKg',seg.loadKg,'speedCode',seg.speedCode, ...
                'group',group,'stage',stage, ...
                'time',q.time(context)-q.time(seg.startIndex),'actionMask',ismember(context,ix)', ...
                'referenceTheta',q.referenceTheta(context),'referenceTorque',q.referenceTorque(context), ...
                'torqueValid',q.torqueValid(context),'waiting',q.waiting(context), ...
                'outputs',cat(3,runs{1}.output(context,:),runs{2}.output(context,:)), ...
                'inputs',cat(3,runs{1}.inputs(context,:),runs{2}.inputs(context,:)), ...
                'recordedAbdomen',q.abdomenActivation(context), ...
                'predictedAbdomen',q.abdomenMAPActivation(context), ...
                'rawMAPAbdomen',q.rawMAPActivation(context,2));
            item.rawMAPBack=q.rawMAPActivation(context,1);
            if isfield(seg,'cadenceBPM'),item.cadenceBPM=seg.cadenceBPM;end
            if isfield(seg,'recordId'),item.recordId=seg.recordId;end
            item.referenceTorque(~item.torqueValid)=NaN;cases{caseCount}=item;
            for source=1:2
                out=runs{source}.output(ix,:);valid=isfinite(out(:,1));
                supported=valid & q.torqueValid(ix) & isfinite(q.referenceTorque(ix));
                complete=all(valid);
                row=struct('Subject',string(subject.shortName),'Source',string(sourceNames{source}), ...
                    'CaseNumber',caseCount,'LoadKg',seg.loadKg,'SpeedCode',seg.speedCode, ...
                    'Complete',complete,'CompletedFraction',mean(valid), ...
                    'ThetaCorrelation',pearson(out(valid,1),q.referenceTheta(ix(valid))), ...
                    'ThetaRMSEDeg',rad2deg(rmse(out(valid,1),q.referenceTheta(ix(valid)))), ...
                    'TorqueCorrelation',pearson(out(supported,3),q.referenceTorque(ix(supported))), ...
                    'TorqueRMSENm',rmse(out(supported,3),q.referenceTorque(ix(supported))), ...
                    'MinimumThetaRad',finite_min(out(:,1)),'MaximumThetaRad',finite_max(out(:,1)), ...
                    'BackPeakNm',finite_max(out(:,4)),'AbdomenPeakNm',finite_max(out(:,5)), ...
                    'BackPeakN',finite_max(out(:,4))/subject.models{source}.parameters(17), ...
                    'AbdomenPeakN',finite_max(out(:,5))/subject.models{source}.parameters(18), ...
                    'AbdomenInputScale',subject.models{source}.abdomenScale, ...
                    'BackMeanNm',mean(out(valid,4)),'AbdomenMeanNm',mean(out(valid,5)), ...
                    'AbdomenBackMeanRatio',mean(out(valid,5))/max(mean(out(valid,4)),1), ...
                    'MovementPDPeakNm',finite_max(abs(out(q.waiting(ix)<.5,22))), ...
                    'BoundaryPeakNm',finite_max(abs(out(:,26))), ...
                    'AssistancePeakNm',finite_max(abs(out(:,25))), ...
                    'MomentumResidualPeakNm',finite_max(abs(out(:,17))), ...
                    'FiberResidualPeak',finite_max(abs(out(:,16))));
                rows=[rows;row]; %#ok<AGROW>
            end
        end
    end
    assert(caseCount==9,'Expected nine recorded conditions per subject.');
    plot_direct_BSE(subject,cases,outputFolder,Show_figures);
    plot_528_scaling(subject,outputFolder,Show_figures);
    if strcmp(char(subject.shortName),'528'),plot_raw_MAP_validation(Show_figures);end
    caseMetrics=struct2table(rows(firstCaseRow:end));
    continuous=continuous_example(subject,groups,cases,caseMetrics,Use_native);
    correlationSummary=plot_requested_figures(subject,continuous,caseMetrics,outputFolder,Show_figures);
    results.subjects{subjectIndex}=struct('id',subject.id,'models',{subject.models}, ...
        'parameterPolicy',subject.parameterPolicy, ...
        'balance',subject.balanceController,'groups',{groups},'cases',{cases}, ...
        'continuousExample',continuous,'correlationSummary',correlationSummary);
end

%% Save raw outputs, metrics and the exact coefficients actually used
results.metrics=struct2table(rows);
results.protocol=struct('muscleUnits',2,'assistanceNm',0,'movementPD',false, ...
    'standingPD',true,'angleBoundaries',false,'referenceTrackingDuringMovement',false, ...
    'parameterPolicy','See each subject package and training report for calibration provenance.');
writetable(results.metrics,fullfile(outputFolder,'case_metrics.csv'));
save(fullfile(outputFolder,'BSE_results.mat'),'results','parameterRows','-v7');
disp(results.metrics(:,{'Subject','Source','LoadKg','SpeedCode','Complete','TorqueCorrelation','ThetaCorrelation','ThetaRMSEDeg'}));
fprintf('\nResults: %s\n',outputFolder);

function run=simulate_sequence(plant,q,inputs)
% There is one initialization per recorded sequence, never between its loads.
n=numel(q.time);out=nan(n,27);last=1;failure='';
[state,out(1,:)]=plant.initialize(inputs(1,:),q.mass(1),q.waiting(1),q.initialState);
for k=2:n
    try
        [state,out(k,:)]=plant.step(state,inputs(k,:),q.mass(k),q.waiting(k),q.dt,10);
        last=k;
    catch exception
        if ~strcmp(exception.identifier,'DirectMTU:Equilibrium'),rethrow(exception);end
        failure=exception.message;
        warning('DirectMTU:Stopped','Stopped at %.3f s: %s',q.time(k),failure);
        break
    end
end
run=struct('output',out,'inputs',inputs,'parameters',plant.parameters, ...
    'balance',plant.balance,'completedSamples',last,'failure',failure);
end

function verify_feedback(run,q)
out=run.output;valid=isfinite(out(:,1));moving=valid & q.waiting<.5;
assert(all(out(valid,24:26)==0,'all'),'Unexpected intent, assistance or angle boundary.');
assert(all(out(moving,20:22)==0,'all'),'Standing feedback entered an action.');
if any(moving)
    assert(max(abs(out(moving,18:19)-run.inputs(moving,:)),[],'all')<1e-12, ...
        'Movement excitation differs from the prepared input.');
end
assert(max(abs(out(valid,3)-out(valid,4)+out(valid,5)))<1e-6);
end
function value=pearson(a,b)
value=NaN;if numel(a)<3 || std(a)<1e-12 || std(b)<1e-12,return,end
c=corrcoef(a,b);value=c(1,2);
end
function value=rmse(a,b)
value=NaN;if ~isempty(a),value=sqrt(mean((a-b).^2));end
end
function value=finite_min(a)
a=a(isfinite(a));value=NaN;if ~isempty(a),value=min(a);end
end
function value=finite_max(a)
a=a(isfinite(a));value=NaN;if ~isempty(a),value=max(a);end
end

function example=continuous_example(subject,groups,cases,metrics,useNative)
% Select one recorded case per load using both torque and angle correlations.
% Each source contributes equally. No model parameter or reference is changed.
score=-inf(1,numel(cases));
for k=1:numel(cases)
    m=metrics(metrics.CaseNumber==k,:);
    if height(m)==2 && all(m.Complete)
        values=[m.TorqueCorrelation;m.ThetaCorrelation];
        if all(isfinite(values)),score(k)=mean(values);end
    end
end
best=-inf;chosen=[];
for a=find(cellfun(@(x)x.loadKg==5,cases))
    for b=find(cellfun(@(x)x.loadKg==10,cases))
        for c=find(cellfun(@(x)x.loadKg==15,cases))
            trial=[a b c];speeds=cellfun(@(x)x.speedCode,cases(trial));
            if numel(unique(speeds))==3 && all(isfinite(score(trial)))
                value=mean(score(trial));
                if value>best,best=value;chosen=trial;end
            end
        end
    end
end
assert(~isempty(chosen), ...
    'Three complete 5/10/15-kg cases with distinct speeds and finite correlations are required.');

% Join only in existing standing intervals. Run the Hill model once across
% the entire new timeline, so no muscle or joint state is reset between loads.
first=groups{cases{chosen(1)}.group}.sequence;
q=struct('time',[],'mass',[],'waiting',[],'referenceTheta',[], ...
    'referenceTorque',[],'torqueValid',[],'referenceFilled',[], ...
    'dt',first.dt,'standingAngle',first.standingAngle,'initialState',first.initialState, ...
    'segments',{cell(1,3)});
input=cell(1,2);input{1}=zeros(0,2);input{2}=zeros(0,2);
parts=cell(1,4);
for k=1:3
    item=cases{chosen(k)};original=groups{item.group}.sequence;
    seg=original.segments{item.stage};
    assert(abs(original.dt-q.dt)<1e-12 && abs(original.standingAngle-q.standingAngle)<1e-12, ...
        'Selected cases must use the same sample interval and standing angle.');
    previous=1;
    if item.stage>1,previous=original.segments{item.stage-1}.endIndex;end
    following=numel(original.time);
    if item.stage<numel(original.segments),following=original.segments{item.stage+1}.startIndex;end
    firstIndex=floor((previous+seg.startIndex)/2);
    lastIndex=floor((seg.endIndex+following)/2);
    if k==1 && item.stage==1,firstIndex=1;end
    if k==3 && item.stage==numel(original.segments),lastIndex=numel(original.time);end
    assert(original.waiting(firstIndex)>.5 && original.waiting(lastIndex)>.5, ...
        'Continuous-example joins must lie inside standing intervals.');
    sourceIndex=(firstIndex:lastIndex)';
    offset=numel(q.mass);
    fields={'mass','waiting','referenceTheta','referenceTorque','torqueValid','referenceFilled'};
    for j=1:numel(fields)
        name=fields{j};values=original.(name)(sourceIndex);
        q.(name)=[q.(name);values(:)]; %#ok<AGROW>
    end
    for source=1:2
        input{source}=[input{source};groups{item.group}.runs{source}.inputs(sourceIndex,:)]; %#ok<AGROW>
    end
    selected=seg;
    selected.startIndex=offset+seg.startIndex-firstIndex+1;
    selected.endIndex=offset+seg.endIndex-firstIndex+1;
    selected.originalGroup=item.group;
    selected.originalStage=item.stage;
    selected.originalStartIndex=seg.startIndex;
    selected.originalEndIndex=seg.endIndex;
    selected.selectionMeanR=score(chosen(k));
    q.segments{k}=selected;
    parts{k}=struct('originalGroup',item.group,'originalStartIndex',firstIndex, ...
        'originalEndIndex',lastIndex,'startIndex',offset+1,'endIndex',numel(q.mass));
end
% Use two more seconds of an existing standing interval to show recovery
% after the last lift. These are defined waiting inputs, not measured data.
waitSource=groups{cases{chosen(1)}.group}.sequence;
waitIndex=(1:round(2/q.dt))';
assert(all(waitSource.waiting(waitIndex)>.5) && all(waitSource.referenceFilled(waitIndex)), ...
    'The final extension must contain only the defined standing interval.');
offset=numel(q.mass);
for j=1:numel(fields)
    name=fields{j};values=waitSource.(name)(waitIndex);
    q.(name)=[q.(name);values(:)]; %#ok<AGROW>
end
for source=1:2
    originalInput=groups{cases{chosen(1)}.group}.runs{source}.inputs;
    input{source}=[input{source};originalInput(waitIndex,:)]; %#ok<AGROW>
end
parts{4}=struct('originalGroup',cases{chosen(1)}.group,'originalStartIndex',waitIndex(1), ...
    'originalEndIndex',waitIndex(end),'startIndex',offset+1,'endIndex',numel(q.mass));
q.time=(0:numel(q.mass)-1)'*q.dt;
runs=cell(1,2);
for source=1:2
    model=subject.models{source};
    plant=dual_MTU_model(model.parameters,subject.balanceController,q.standingAngle,useNative);
    runs{source}=simulate_sequence(plant,q,input{source});
    verify_feedback(runs{source},q);
    assert(runs{source}.completedSamples==numel(q.time), ...
        'The selected continuous simulation stopped before the end.');
end

% Report correlations on the new continuous replay, not on spliced outputs.
rows=struct([]);
for k=1:3
    seg=q.segments{k};ix=(seg.startIndex:seg.endIndex)';
    row=struct('LoadKg',seg.loadKg,'SpeedCode',seg.speedCode, ...
        'OriginalGroup',seg.originalGroup,'SelectionMeanR',seg.selectionMeanR, ...
        'TorqueRExperimental',NaN,'TorqueRMAP',NaN,'ThetaRExperimental',NaN,'ThetaRMAP',NaN);
    for source=1:2
        out=runs{source}.output(ix,:);
        observed=~logical(q.referenceFilled(ix));
        validTorque=observed & logical(q.torqueValid(ix)) & isfinite(q.referenceTorque(ix)) & isfinite(out(:,3));
        validTheta=observed & isfinite(q.referenceTheta(ix)) & isfinite(out(:,1));
        torqueR=pearson(out(validTorque,3),q.referenceTorque(ix(validTorque)));
        thetaR=pearson(out(validTheta,1),q.referenceTheta(ix(validTheta)));
        if source==1
            row.TorqueRExperimental=torqueR;row.ThetaRExperimental=thetaR;
        else
            row.TorqueRMAP=torqueR;row.ThetaRMAP=thetaR;
        end
    end
    rows=[rows;row]; %#ok<AGROW>
end
example=struct('sequence',q,'runs',{runs},'parts',{parts}, ...
    'caseMetrics',struct2table(rows),'selectionMeanR',best);
fprintf('\nSelected continuous example | %s\n',subject.shortName);
disp(example.caseMetrics);
end

function summary=plot_requested_figures(subject,example,metrics,folder,showFigures)
% The reference in a waiting interval is the original defined standing baseline.
visible='off';if showFigures,visible='on';end
colors=[0 .40 .72;.53 .19 .65];
q=example.sequence;
fig=figure('Name','Continuous BSE example','Color','w','Visible',visible, ...
    'Position',[40 40 1500 850]);
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
for panel=1:2
    ax=nexttile(layout);hold(ax,'on');
    if panel==1
        reference=q.referenceTorque;reference(~logical(q.torqueValid))=NaN;
        column=3;ylabel(ax,'Net muscle torque (N m)');
    else
        reference=q.referenceTheta;column=1;ylabel(ax,'Trunk angle (rad)');
    end
    recorded=reference;recorded(logical(q.referenceFilled))=NaN;
    standing=reference;standing(~logical(q.referenceFilled) | q.waiting<.5)=NaN;
    h(1)=plot(ax,q.time,recorded,'r--','LineWidth',1.35);
    h(2)=plot(ax,q.time,example.runs{1}.output(:,column),'Color',colors(1,:),'LineWidth',1.4);
    h(3)=plot(ax,q.time,example.runs{2}.output(:,column),'Color',colors(2,:),'LineWidth',1.4);
    h(4)=plot(ax,q.time,standing,':','Color',[.55 .55 .55],'LineWidth',1.1);
    grid(ax,'on');box(ax,'on');xlim(ax,[q.time(1) q.time(end)]);
    set(ax,'FontName','Times New Roman','FontSize',11);
    if panel==1
        limits=ylim(ax);span=max(diff(limits),1);
        ylim(ax,[min(0,limits(1)) limits(2)+.15*span]);limits=ylim(ax);
        for k=1:3
            seg=q.segments{k};at=mean(q.time([seg.startIndex seg.endIndex]));
            text(ax,at,limits(2)-.06*diff(limits), ...
                sprintf('%g kg | %g BPM',seg.loadKg,seg.speedCode/5), ...
                'HorizontalAlignment','center','Color',[.75 .15 .15],'FontSize',11);
        end
    else
        xlabel(ax,'Time (s)');
    end
end
key=legend(h,{'Reference: CEINMS torque / IK angle','BSE: experimental activation', ...
    'BSE: MAP activation','Defined standing baseline'}, ...
    'Orientation','horizontal','NumColumns',2,'FontSize',10);
key.Layout.Tile='south';
title(layout,{'Trunk Musculoskeletal Model Evaluation with', ...
    'CEINMS Kinematics and Dynamics Reference'}, ...
    'FontName','Times New Roman','FontSize',16,'Interpreter','none');
stem=fullfile(folder,sprintf('%s_continuous_torque_theta',subject.shortName));
format_demo_figure(fig);
% Enlarge this title after the shared formatter has set its standard size.
layout.Title.FontSize=1.5*layout.Title.FontSize;
drawnow;
exportgraphics(fig,[stem '.png'],'Resolution',180);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig(fig,[stem '.fig']);
export_paper_figure(fig,[stem '_paper'],2.8);
if ~showFigures,close(fig);end
writetable(example.caseMetrics,fullfile(folder,sprintf('%s_continuous_example_cases.csv',subject.shortName)));

% The bar height is the equal-weight mean of the nine per-case Pearson R values.
names={'Experimental activation','MAP activation'};
means=nan(2,2);counts=zeros(2,2);
for source=1:2
    m=metrics(metrics.Source==string(names{source}),:);
    values={m.TorqueCorrelation,m.ThetaCorrelation};
    for quantity=1:2
        x=values{quantity};x=x(isfinite(x));counts(quantity,source)=numel(x);
        if ~isempty(x),means(quantity,source)=mean(x);end
    end
end
summary=table(string(names(:)),means(1,:)',counts(1,:)',means(2,:)',counts(2,:)', ...
    'VariableNames',{'Source','TorqueMeanR','TorqueCases','ThetaMeanR','ThetaCases'});
writetable(summary,fullfile(folder,sprintf('%s_correlation_summary.csv',subject.shortName)));
fig=figure('Name','Mean BSE correlation','Color','w','Visible',visible, ...
    'Position',[80 80 900 620]);
ax=axes(fig);hold(ax,'on');
bars=bar(ax,1:2,means,.65,'grouped');
for source=1:2
    bars(source).FaceColor=colors(source,:);
    for quantity=1:2
        if isfinite(means(quantity,source))
            text(ax,bars(source).XEndPoints(quantity),means(quantity,source)+.025, ...
                sprintf('%.3f\nn = %d',means(quantity,source),counts(quantity,source)), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
        end
    end
end
xlim(ax,[.35 2.65]);finite=means(isfinite(means));
if isempty(finite),lower=0;else,lower=min(0,min(finite)-.1);end
ylim(ax,[lower 1.12]);xticks(ax,[1 2]);xticklabels(ax,{'Torque','Theta'});
ylabel(ax,'Mean Pearson correlation, R');grid(ax,'on');box(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',12);
legend(ax,bars,{'Experimental activation','MAP activation'},'Location','northoutside', ...
    'Orientation','horizontal','Box','off');
title(ax,sprintf('%s | Mean correlation across nine conditions',subject.shortName), ...
    'Interpreter','none','FontWeight','normal');
subtitle(ax,'Equal-weight mean of available per-case correlations; n is shown for each bar.');
stem=fullfile(folder,sprintf('%s_mean_correlation',subject.shortName));
format_demo_figure(fig);
exportgraphics(fig,[stem '.png'],'Resolution',180);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig(fig,[stem '.fig']);
if ~showFigures,close(fig);end
end
