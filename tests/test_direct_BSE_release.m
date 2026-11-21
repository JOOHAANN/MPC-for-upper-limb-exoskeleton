function test_direct_BSE_release
%TEST_DIRECT_BSE_RELEASE Verify recorded inputs, fixed vectors and no hidden PD.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');addpath(folder);
current=load(fullfile(folder,'data','subject_528.mat'));
if isfield(current.subject,'MAPScaling')
    test_load_scaled_528;return
end
data=load(fullfile(folder,'output','BSE_results.mat'));results=data.results;
assert(height(results.metrics)==18*numel(results.subjects));
assert(results.protocol.muscleUnits==2 && ~results.protocol.movementPD ...
    && results.protocol.standingPD && ~results.protocol.angleBoundaries);
assert(~results.protocol.referenceTrackingDuringMovement && results.protocol.assistanceNm==0);
assert(all(results.metrics.Complete));
for k=1:numel(results.subjects)
    subject=results.subjects{k};
    if contains(subject.id,'528')
        short='528';
    elseif contains(subject.id,'617')
        short='617';
    elseif contains(subject.id,'725')
        short='725';
    else
        short='730';
    end
    before=load(fullfile(root,'experiments','direct_MTU_training','before_training','data',['subject_' short '.mat']));
    assert(isequaln(subject.balance,before.subject.balanceController));
    package=load(fullfile(folder,'data',['subject_' short '.mat']));map=package.subject.abdominalMAP;
    assert(~map.usesIK && ~map.usesTorque && ~map.usesBSEOutput && ~map.usesBackToAbdomenRecruitment);
    for source=1:2
        model=subject.models{source};assert(numel(model.parameters)==24);
        assert(~isfield(model,'neuralParameters'));
        assert(~isfield(model,'generator') && model.abdomenScale>0);
        for group=1:numel(subject.groups)
            q=subject.groups{group}.sequence;run=subject.groups{group}.runs{source};a=run.output;
            recorded=rmfield(q,'abdomenMAPActivation');
            assert(isequaln(recorded,before.subject.sequences{group}), ...
                'Recorded signals, references or frozen back MAP changed.');
            for stage=1:numel(q.segments)
                seg=q.segments{stage};ix=seg.startIndex:seg.endIndex;
                assert(all(q.waiting(ix)==0),'Standing feedback is enabled inside an action.');
                duration=q.time(ix(end))-q.time(ix(1));
                predicted=predict_abdominal_MAP(map,seg.loadKg,60/duration, ...
                    (q.time(ix)-q.time(ix(1)))/duration);
                assert(max(abs(predicted-q.abdomenMAPActivation(ix)))<1e-12);
            end
            assert(isequaln(run.parameters(:),model.parameters(:)),'Parameters changed between records.');
            assert(all(isfinite(a),'all') && run.completedSamples==numel(q.time));
            assert(all(a(:,24:26)==0,'all'),'Unexpected intent, assistance or angle boundaries.');
            moving=q.waiting<.5;
            assert(all(a(moving,20:22)==0,'all'),'Human feedback entered a movement interval.');
            input=prepare_activation(q,model,source);
            assert(isequaln(input,run.inputs));
            channel=q.abdomenActivation;
            if source==2,channel=q.abdomenMAPActivation;end
            assert(max(abs(input(:,2)-max(0,min(1,model.abdomenScale*channel))))<1e-12);
            changed=q;changed.backActivation=1-q.backActivation;
            changed.referenceTheta=zeros(size(q.time));changed.referenceTorque=1000+zeros(size(q.time));
            independent=prepare_activation(changed,model,source);
            assert(isequaln(independent(:,2),input(:,2)), ...
                'Abdominal command depends on back activation or angle/torque references.');
            assert(max(abs(a(moving,18:19)-input(moving,:)),[],'all')<1e-12);
            assert(max(abs(a(:,3)-a(:,4)+a(:,5)))<1e-6);
            assert(max(abs(a(:,17)))<1e-5 && max(abs(a(:,16)))<1e-4);
        end
    end
end
% Report accuracy separately from physics/input integrity. A failed fitting
% target must not be hidden by weakening a physics test or modifying outputs.
t=results.metrics;
for source={'Experimental activation','MAP activation'}
    rows=t.Subject=="528" & t.Source==source{1};assert(sum(rows)==9);
    torqueTarget=.8;thetaTarget=.7;
    if strcmp(source{1},'MAP activation'),torqueTarget=.75;thetaTarget=.6;end
    fprintf('%s: mean torque r %.4f (target %.2f), theta r %.4f (target %.2f), abdomen/back %.3f.\n', ...
        source{1},mean(t.TorqueCorrelation(rows)),torqueTarget, ...
        mean(t.ThetaCorrelation(rows)),thetaTarget,mean(t.AbdomenBackMeanRatio(rows)));
end
fprintf('PASS: completed replays, unchanged recordings/references/back MAP, strict abdominal sources, one fixed vector/gain per source, zero movement feedback/assistance/boundary torque, and consistent dual-MTU mechanics.\n');
end
