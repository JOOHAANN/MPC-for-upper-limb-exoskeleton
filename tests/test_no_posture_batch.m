function test_no_posture_batch
%TEST_NO_POSTURE_BATCH Scalar, batched and nonuniform prediction must agree.
root=fileparts(fileparts(mfilename('fullpath')));demo=fullfile(root,'Demo','NMPC');
addpath(demo,fullfile(root,'Demo','BSE'));
subjects=bse_subjects('528');subject=subjects{1};q=subject.sequences{1};
programs=load(fullfile(demo,'data','no_human_PD_motor_programs_528.mat'));
a=programs.programs{1,2}(1,:);
plant=bse_without_posture(subject.models{2},subject.balanceController,q.standingAngle,.8);
plant=no_posture_blocked_plant(plant);
[state,~]=plant.initialize(a,0,1,q.initialState);
o=struct('predictionDt',.01,'predictionSubsteps',10,'rateNmPerSec',200, ...
    'blockPoints',[5 5 5 5 10 10 15 15 15 15]);
forecast=struct('activation',repmat(a,100,1),'mass',zeros(100,1),'waiting',ones(100,1));
controls=[zeros(10,1),linspace(-4,4,10)',linspace(6,-3,10)',2*ones(10,1)];
batch=plant.previewBatch(state,controls,forecast,o);
for j=1:4
    scalar=plant.preview(state,controls(:,j),forecast,o);
    assert(~isempty(scalar) && max(abs(scalar-batch(:,:,j)),[],'all')<1e-7);
end
% Replay the same nonuniform plan through actual 10-ms calls.
request=repelem(controls(:,2),o.blockPoints(:)); actual=zeros(100,27);
for k=1:100
    previous=state(29);if previous*request(k)<0,previous=0;end
    u=sign(request(k))*min(abs(request(k)),abs(previous)+2);
    [state,actual(k,:)]=plant.step(state,a,0,1,u,.01,10);
end
assert(max(abs(actual-batch(:,:,2)),[],'all')<1e-7);
fprintf('PASS: parallel batch equals scalar prediction and nonuniform actual replay.\n');
end
