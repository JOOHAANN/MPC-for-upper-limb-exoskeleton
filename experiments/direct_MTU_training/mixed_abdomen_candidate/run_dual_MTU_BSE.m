clc
clearvars -except Subject_list Show_figures Use_native Output_label
close all

%% Direct-activation BSE: two MTUs, standing PD only, no angle limits
% One fixed parameter vector per subject and activation source, across nine
% recorded conditions. Experimental and frozen trained-MAP inputs are used.
% No external assistance, movement PD, decoded angle or boundary torque.
if ~exist('Subject_list','var'),Subject_list={'528','617','725','730'};end
if ~exist('Show_figures','var'),Show_figures=true;end
if ~exist('Use_native','var'),Use_native=true;end
if ~exist('Output_label','var'),Output_label='';end
folder=fileparts(mfilename('fullpath'));addpath(folder);
outputFolder=fullfile(folder,'output',Output_label);
if ~isfolder(outputFolder),mkdir(outputFolder);end
sourceNames={'Experimental activation','MAP activation'};
results=struct;results.subjects=cell(1,numel(Subject_list));rows=struct([]);
parameterRows=struct([]);

%% Read each subject and keep coefficients fixed throughout its sequences
for subjectIndex=1:numel(Subject_list)
    data=load(fullfile(folder,'data',['subject_' char(Subject_list{subjectIndex}) '.mat']));
    subject=data.subject;cases=cell(1,9);caseCount=0;groups=cell(size(subject.sequences));
    fprintf('\nDirect dual-MTU BSE | %s\n',subject.id);
    for source=1:2
        parameterRow=struct('Subject',string(subject.shortName), ...
            'Source',string(sourceNames{source}),'Parameters',{subject.models{source}.parameters});
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
            % Include real neighboring samples to show standing-PD switching.
            context=max(1,seg.startIndex-round(.5/q.dt)):min(numel(q.time),seg.endIndex+round(1/q.dt));
            item=struct('loadKg',seg.loadKg,'speedCode',seg.speedCode, ...
                'time',q.time(context)-q.time(seg.startIndex),'actionMask',ismember(context,ix)', ...
                'referenceTheta',q.referenceTheta(context),'referenceTorque',q.referenceTorque(context), ...
                'torqueValid',q.torqueValid(context),'waiting',q.waiting(context), ...
                'outputs',cat(3,runs{1}.output(context,:),runs{2}.output(context,:)));
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
    results.subjects{subjectIndex}=struct('id',subject.id,'models',{subject.models}, ...
        'parameterPolicy',subject.parameterPolicy, ...
        'balance',subject.balanceController,'groups',{groups},'cases',{cases});
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
