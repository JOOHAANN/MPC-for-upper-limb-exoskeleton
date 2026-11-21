function test_current_bse_control
%TEST_CURRENT_BSE_CONTROL Unassisted parity and exact applied-assistance gating.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
addpath(folder,fullfile(root,'Demo','NMPC'));
data=load(fullfile(folder,'data','subject_528.mat'));subject=data.subject;
for source=1:2
    q=subject.sequences{1};model=subject.models{source};a=prepare_activation(q,model,source);
    original=dual_MTU_model(model.parameters,subject.balanceController,q.standingAngle,true);
    current=current_bse_plant(model,subject.balanceController,q.standingAngle,[5 3 .005 .3 .012],[80 200]);
    [s,x]=original.initialize(a(1,:),q.mass(1),q.waiting(1),q.initialState);
    [t,y]=current.initialize(a(1,:),q.mass(1),q.waiting(1),q.initialState);
    assert(max(abs(x-y))<1e-9);
    for k=2:numel(q.time)
        [s,x]=original.step(s,a(k,:),q.mass(k),q.waiting(k),q.dt,10);
        [t,y]=current.step(t,a(k,:),q.mass(k),q.waiting(k),0,q.dt,10);
        assert(max(abs(x-y))<1e-7,'Zero-assistance model no longer matches calibrated BSE.');
        assert(t(69)==0 && y(25)==0 && y(26)==0);
        if ~q.waiting(k),assert(all(y(20:22)==0) && max(abs(y(18:19)-a(k,:)))<1e-12);end
        assert(max(abs(t(1:2)-t(35:36)))<1e-8);
    end
    % Exercise nonzero assistance during a recorded action, then withdraw it.
    ix=q.segments{1}.startIndex;
    [t,~]=current.initialize(a(1,:),q.mass(1),q.waiting(1),q.initialState);
    for k=2:ix+150
        command=0;if k>=ix+50 && k<ix+110,command=30;end
        [t,y]=current.step(t,a(k,:),q.mass(k),q.waiting(k),command,.01,10);
        if command>0
            assert(y(25)>0 && y(25)<=30 && y(24)==t(35));
        else
            assert(y(25)==0 && t(69)==0 && y(24)==0);
            if ~q.waiting(k),assert(max(abs(y(18:19)-a(k,:)))<1e-12);end
        end
        assert(y(26)==0 && abs(y(17))<1e-5);
    end
end
fprintf('PASS: latest BSE parity at zero assistance, causal shadow, movement PD gated by applied torque, instant withdrawal and zero boundaries.\n');
end
