% Isolate movement Kp/Kd without changing inputs, allocation blend or physics.
% This is an unassisted sensitivity experiment, not a retrained replacement.
folder=fileparts(mfilename('fullpath'));root=fileparts(folder);
addpath(fullfile(root,'Demo','BSE'));
subjects=bse_subjects({'528','617','725','730'});rows=struct([]);
for k=1:numel(subjects)
    subject=subjects{k};
    for source=1:2
        model=subject.models{source};removed=model;
        removed.neuralParameters(4:5)=0;
        activeErrors=[];removedErrors=[];differenceErrors=[];
        qmin=inf;qmax=-inf;complete=true;completed=0;
        for group=1:numel(subject.sequences)
            q=subject.sequences{group};a=bse_inputs(q,model,source);
            before=Torque_function(model,subject.balanceController,q,a,true);
            plant=bse_model(removed,subject.balanceController,q.standingAngle,true);
            [state,first]=plant.initialize(a(1,:),q.mass(1),q.waiting(1),q.initialState);
            after=nan(size(before));after(1,:)=first;last=1;
            for j=2:numel(q.time)
                try
                    [state,after(j,:)]=plant.step(state,a(j,:),q.mass(j),q.waiting(j),0,q.dt,10);
                    last=j;
                catch exception
                    if ~strcmp(exception.identifier,'BSE:Equilibrium'),rethrow(exception);end
                    complete=false;break
                end
            end
            qmin=min(qmin,min(after(1:last,1)));qmax=max(qmax,max(after(1:last,1)));
            for stage=1:numel(q.segments)
                seg=q.segments{stage};ix=seg.startIndex:seg.endIndex;
                if seg.endIndex>last,continue,end
                activeErrors(end+1)=rad2deg(sqrt(mean((before(ix,1)-q.referenceTheta(ix)).^2))); %#ok<SAGROW>
                removedErrors(end+1)=rad2deg(sqrt(mean((after(ix,1)-q.referenceTheta(ix)).^2))); %#ok<SAGROW>
                differenceErrors(end+1)=rad2deg(sqrt(mean((after(ix,1)-before(ix,1)).^2))); %#ok<SAGROW>
                completed=completed+1;
            end
        end
        names={'Experimental','MAP'};
        row=struct('Subject',string(subject.shortName),'Source',string(names{source}), ...
            'CompletedCases',completed,'CompleteSimulation',complete, ...
            'OriginalMeanThetaRMSEDeg',mean(activeErrors), ...
            'MovementKpKdZeroMeanThetaRMSEDeg',mean(removedErrors), ...
            'MeanThetaChangeRMSEDeg',mean(differenceErrors), ...
            'MovementKpKdZeroMinimumThetaRad',qmin,'MovementKpKdZeroMaximumThetaRad',qmax, ...
            'OriginalFeedbackBlend',model.neuralParameters(7));
        rows=[rows;row]; %#ok<AGROW>
    end
end
audit=struct2table(rows);disp(audit);
writetable(audit,fullfile(folder,'movement_PD_ablation_audit.csv'));
fprintf('Only movement Kp/Kd were set to zero in memory. All files, inputs, allocation blends, standing gains and physical parameters were unchanged.\n');
