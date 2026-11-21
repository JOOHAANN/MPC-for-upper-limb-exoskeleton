function test_load_scaled_528(checkReplay)
%TEST_LOAD_SCALED_528 Legacy MAP parity, shared scaling and reference integrity.
if nargin<1,checkReplay=true;end
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');addpath(folder);
data=load(fullfile(folder,'data','subject_528.mat'));subject=data.subject;
old=load(fullfile(root,'experiments','load_scaled_528','before_load_rescale','data','subject_528.mat'));
mapData=load(fullfile(folder,'data','legacy_MAP_528.mat'));s=subject.MAPScaling;
assert(numel(mapData.records)==25 && isequal(s.speedCodes,[25 35 45 55 65]));
assert(isequal(s.loadKg,[5;10;15]));
assert(~subject.MAP.abdomen.usesBackToAbdomenRecruitment);
assert(~subject.MAP.abdomen.usesIK && ~subject.MAP.abdomen.usesTorque);
for loadIndex=1:3
    rows=cellfun(@(r)r.loadKg==s.loadKg(loadIndex),mapData.records);assert(sum(rows)==5);
    values=cellfun(@(r)r.experimental,mapData.records(rows),'UniformOutput',false);
    average=mean(cat(3,values{:}),3);
    assert(max(abs(average-s.meanExperimental(:,:,loadIndex)),[],'all')<1e-12);
    affine=s.meanRawMAP(:,:,loadIndex).*s.gain(loadIndex,:)+s.offset(loadIndex,:);
    assert(max(abs(min(affine)-min(average)))<1e-12);
    assert(max(abs(max(affine)-max(average)))<1e-12);
end

%% Compare every legacy back MAP profile with the actual original function
previous=pwd;cleanup=onCleanup(@()cd(previous));cd(subject.MAP.backSource);
for weight=[5 10 15]
    for speed=[25 35 45 55 65]
        [~,~,MEP]=MAP_model(speed,weight,'SM','baseline_data.mat');
        expected=(MEP(:,1)+MEP(:,2))/12+(MEP(:,3)+MEP(:,4))/4 ...
            +(MEP(:,1)+MEP(:,2))/6;
        actual=predict_528_MAP(subject.MAP,weight,speed,subject.MAP.phase);
        assert(max(abs(actual(:,1)-expected))<1e-12,'Original MAP fusion changed.');
    end
end
cd(previous);clear cleanup

protected={'time','dt','referenceTheta','referenceTorque','torqueValid','referenceFilled', ...
    'mass','nominalLoad','waiting','initialState','standingAngle','segments'};
for group=1:numel(subject.sequences)
    q=subject.sequences{group};before=old.subject.sequences{group};
    for k=1:numel(protected)
        name=protected{k};assert(isequaln(q.(name),before.(name)),'A reference or timeline changed.');
    end
    for source=1:2
        assert(subject.models{source}.abdomenScale==1);
        input=prepare_activation(q,subject.models{source},source);
        changed=q;changed.referenceTheta(:)=100;changed.referenceTorque(:)=-100;
        assert(isequaln(input,prepare_activation(changed,subject.models{source},source)));
    end
    for stage=1:3
        seg=q.segments{stage};ix=seg.startIndex:seg.endIndex;
        phase=(q.time(ix)-q.time(ix(1)))/(q.time(ix(end))-q.time(ix(1)));
        raw=predict_528_MAP(subject.MAP,seg.loadKg,seg.speedCode,phase);
        expected=rescale_MAP_by_load(raw,seg.loadKg,s);
        assert(max(abs(raw-q.rawMAPActivation(ix,:)),[],'all')<1e-12);
        assert(max(abs(expected-[q.backActivation(ix,2) q.abdomenMAPActivation(ix)]),[],'all')<1e-12);
        assert(all(q.waiting(ix)==0));
    end
end
assert(isequaln(subject.balanceController,old.subject.balanceController));
fprintf('PASS: original MAP output/fusion, five-speed experimental mean extrema, load-only scaling, unchanged references/timelines/standing balance.\n');

if ~checkReplay,return,end
data=load(fullfile(folder,'output','BSE_results.mat'));results=data.results;
assert(height(results.metrics)==18 && all(results.metrics.Complete));
for group=1:3
    for source=1:2
        run=results.subjects{1}.groups{group}.runs{source};out=run.output;
        q=subject.sequences{group};moving=q.waiting<.5;
        assert(isequaln(run.parameters(:),subject.models{source}.parameters(:)));
        assert(isequaln(run.inputs,prepare_activation(q,subject.models{source},source)));
        assert(all(isfinite(out),'all') && all(out(:,24:26)==0,'all'));
        assert(all(out(moving,20:22)==0,'all'));
        assert(max(abs(out(moving,18:19)-run.inputs(moving,:)),[],'all')<1e-12);
        assert(max(abs(out(:,3)-out(:,4)+out(:,5)))<1e-6);
        assert(max(abs(out(:,16)))<1e-4 && max(abs(out(:,17)))<1e-5);
    end
end
fprintf('PASS: 18 completed replays, one fixed physical vector per source, no fitted input gains, zero movement PD/assistance/boundary torque.\n');
end
