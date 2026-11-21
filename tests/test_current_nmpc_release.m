function test_current_nmpc_release
%TEST_CURRENT_NMPC_RELEASE Audit the actual controller runs, inputs and model.
root=fileparts(fileparts(mfilename('fullpath')));folder=fullfile(root,'Demo','NMPC');
bse=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');addpath(folder,bse);
data=load(fullfile(folder,'output','current_BSE_528','current_BSE_NMPC_528.mat'));
r=data.results;settings=data.settings;
model=load(fullfile(bse,'data','subject_528.mat'));subject=model.subject;
before=load(fullfile(bse,'output','BSE_results.mat'));
assert(strcmp(r.modelSHA256,nmpc_file_hash(fullfile(bse,'data','subject_528.mat'))));
assert(height(r.metrics)==54 && isequaln(r.models,subject.models));
assert(isequaln(r.balanceController,subject.balanceController) && isequaln(r.MAPScaling,subject.MAPScaling));
assert(isequal(settings.actuator,[80 200]) && isequal(settings.humanPD(1:4),[5 3 .005 .3]));
assert(settings.protocol.thresholdFraction==.75 && ~settings.protocol.matchedEffort);
assert(iscell(settings.pid) && numel(settings.pid)==2, ...
    'Exactly one PID tuning per activation source is required.');
for group=1:3
    q=r.groups{group}.sequence;
    assert(isequaln(q,subject.sequences{group}));
    for source=1:2
        pair=r.groups{group}.runs{source};inputs=prepare_activation(q,subject.models{source},source);
        assert(isequaln(pair.pid.pid,settings.pid{source}), ...
            'Each source must use the same PID tuning for all nine cases and both plots.');
        expected=before.results.subjects{1}.groups{group}.runs{source}.output;
        assert(max(abs(pair.baseline-expected),[],'all')<1e-7);
        assert(isequaln(pair.options,settings.controllers{source}));
        for stage=1:3
            seg=q.segments{stage};ix=seg.startIndex:seg.endIndex;
            assert(max(abs(pair.threshold(ix)-.75*max(pair.baseline(ix,4))))<1e-10);
        end
        for name={'pid','nmpc'}
            run=pair.(name{1});a=run.actual;moving=q.waiting<.5;off=moving & abs(a(:,25))<=1e-9;
            on=moving & abs(a(:,25))>1e-9;
            if strcmp(name{1},'nmpc'),assert(isequaln(run.inputs,inputs));end
            assert(all(isfinite(a),'all') && min(a(:,25))>=-1e-9 && max(a(:,25))<=80+1e-8);
            assert(max(diff(a(:,25)))/q.dt<=200+1e-6);
            assert(all(a(:,26)==0) && all(a(moving,20:22)==0,'all'));
            assert(max(abs(a(off,18:19)-inputs(off,:)),[],'all')<1e-12 && all(a(off,24)==0));
            assert(max(abs(a(moving,19)-inputs(moving,2)))<1e-12,'Abdominal motion input was altered by PD.');
            assert(max(abs(a(on,24)-pair.baseline(on,1)))<1e-7,'Movement PD is not targeting its causal unassisted shadow.');
            assert(max(abs(a(:,3)-a(:,4)+a(:,5)))<1e-6);
            assert(max(abs(a(:,16)))<1e-4 && max(abs(a(:,17)))<1e-5);
        end
        assert(max(abs(pair.nmpc.target-min(pair.baseline(:,4),pair.threshold)))<1e-8);
    end
end
fprintf('PASS: all 18 PID/NMPC pairs use the current BSE, fixed source parameters, unchanged data/scaling, matched actuator, causal gated back PD, unchanged abdominal motion input and zero boundaries.\n');
end
