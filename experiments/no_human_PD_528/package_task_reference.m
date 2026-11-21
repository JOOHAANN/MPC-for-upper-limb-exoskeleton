% Freeze the previously calibrated unassisted MAP movement as a declared task.
% This is an offline reference generator, not a feedback loop in the ablated
% human plant. Neither controller receives future experimental EMG or IK.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'));
subjects=bse_subjects('528'); subject=subjects{1};
tasks=cell(1,3);
for group=1:3
    q=subject.sequences{group}; model=subject.models{2};
    planner=bse_model(model,subject.balanceController,q.standingAngle,true);
    activation=bse_inputs(q,model,2);
    nominal=planner.simulate(activation,q.mass,q.waiting,q.initialState,q.dt,10);
    thresholds=zeros(numel(q.time),2);
    for source=1:2
        sourceModel=subject.models{source};
        sourcePlant=bse_model(sourceModel,subject.balanceController,q.standingAngle,true);
        sourceInput=bse_inputs(q,sourceModel,source);
        baseline=sourcePlant.simulate(sourceInput,q.mass,q.waiting,q.initialState,q.dt,10);
        limits=zeros(1,3);
        for stage=1:3
            seg=q.segments{stage}; ix=seg.startIndex:seg.endIndex;
            limits(stage)=.75*max(baseline(ix,4)); thresholds(ix,source)=limits(stage);
        end
        thresholds(logical(q.waiting),source)=min(limits);
    end
    tasks{group}=struct('theta',nominal(:,1),'velocity',nominal(:,2), ...
        'torque',nominal(:,4),'thresholds',thresholds);
end
modelSHA256=subject.modelSHA256;
provenance='Offline unassisted rollout of the frozen calibrated MAP model; original posture gains exist only in this reference-generation step, never in the new actual or prediction plants.';
save(fullfile(root,'Demo','NMPC','data','no_human_PD_task_528.mat'), ...
    'tasks','modelSHA256','provenance','-v7');
