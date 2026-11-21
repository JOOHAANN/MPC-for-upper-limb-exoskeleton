function test_direct_dual_MTU
%TEST_DIRECT_DUAL_MTU Physics parity, exact movement gating and no boundaries.
root=fileparts(fileparts(mfilename('fullpath')));
parent=fullfile(root,'Demo','BSE');folder=fullfile(parent,'Direct_Activation_No_Limits');
addpath(parent,folder);
data=load(fullfile(folder,'data','subject_528.mat'));subject=data.subject;
p=subject.models{1}.parameters;balance=subject.balanceController;
native=dual_MTU_model(p,balance,1.5,true);
readable=dual_MTU_model(p,balance,1.5,false);
for waiting=[0 1]
    a=[.14 .025];initial=[1.2 .03];mass=5;
    [s,n]=native.initialize(a,mass,waiting,initial);
    [t,m]=readable.initialize(a,mass,waiting,initial);
    assert(max(abs(n-m))<1e-7);
    for k=1:40
        a=[.14+.02*sin(k/10),.025];flag=waiting;
        if k>20,flag=1-waiting;end
        [s,n]=native.step(s,a,mass,flag,.01,10);
        [t,m]=readable.step(t,a,mass,flag,.01,10);
        assert(max(abs(n-m))<1e-6,'Readable and native physical models differ.');
        assert(all(n(24:26)==0));
        if flag==0
            assert(all(n(20:22)==0) && max(abs(n(18:19)-a))<1e-12);
        end
    end
end
% The same existing dual-MTU physics, with the old neural path/boundary off.
originalData=load(fullfile(parent,'data','subject_528.mat'));
model=originalData.subject.models{1};model.parameters([1:14 19:28])=p;
model.parameters(15:18)=[1 0 0 0];model.parameters(29)=0;
model.neuralParameters=zeros(1,16);model.neuralParameters([2 3 6])=[1 .1 1];
original=bse_model(model,balance,1.5,true);
for waiting=[0 1]
    a=[.14 .025];[s,n]=native.initialize(a,5,waiting,[1.2 .03]);
    [t,m]=original.initialize(a,5,waiting,[1.2 .03]);
    columns=[1:19 25:27];assert(max(abs(n(columns)-m(columns)))<1e-7);
    for k=1:40
        [s,n]=native.step(s,a,5,waiting,.01,10);
        [t,m]=original.step(t,a,5,waiting,0,.01,10);
        assert(max(abs(n(columns)-m(columns)))<1e-6,'Original dual-MTU mechanics changed.');
    end
end
% Angles outside the previous soft limits remain raw states, without a spring.
for angle=[-.2 2.2]
    [~,out]=native.initialize([.02 .01],0,0,[angle 0]);
    assert(out(1)==angle && out(26)==0);
end
% Whole-sequence acceleration must preserve the step-by-step calculation.
inputs=repmat([.14 .025],60,1);mass=5*ones(60,1);waiting=[ones(20,1);zeros(20,1);ones(20,1)];
[fast,completed]=native.simulate(inputs,mass,waiting,[1.2 .03],.01,10);
slow=zeros(60,27);[state,slow(1,:)]=native.initialize(inputs(1,:),5,1,[1.2 .03]);
for k=2:60
    [state,slow(k,:)]=native.step(state,inputs(k,:),5,waiting(k),.01,10);
end
assert(completed==60 && max(abs(fast-slow),[],'all')<1e-9);
fprintf('PASS: original dual-MTU mechanics, readable/native parity, standing-only feedback, direct movement input and zero angle-boundary torque.\n');
end
