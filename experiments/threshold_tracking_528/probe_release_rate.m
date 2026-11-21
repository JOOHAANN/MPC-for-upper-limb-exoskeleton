% Smooth release within the unchanged 200-N-m/s physical actuator bound.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
source=1; group=3; sequence=subject.sequences{group}; model=subject.models{source};
plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
inputs=bse_inputs(sequence,model,source); planned=bse_inputs(sequence,model,2);
cache=load(fullfile(root,'Demo','NMPC','output','pid_nmpc_528','nmpc_group3_source1.mat'));
count=1600; sequence.time=sequence.time(1:count);
sequence.mass=sequence.mass(1:count); sequence.waiting=sequence.waiting(1:count);
inputs=inputs(1:count,:); planned=planned(1:count,:);
for rate=[20 40 60 100]
    options=cache.options; options.blocks=8;
    options.activationUncertainty=zeros(16,2); options.activationEdgeTimesSeconds=[];
    options.predictionErrorGuardScale=.3; options.predictionErrorGuardFloorNm=1;
    options.backTorqueMarginNm=.5; options.releaseRateNmPerSec=rate;
    options.backTorqueUnderWeight=200; options.trackingHeadroomNm=.5;
    options.recoveryBackTargetNm=0; options.earlyRecovery=false;
    fprintf('\nRELEASE TRIAL | rate %.3g\n',rate);
    run=simulate_threshold_nmpc(plant,sequence,inputs,planned,cache.threshold(1:count),options);
    seg=sequence.segments{1}; ix=seg.startIndex:seg.endIndex;
    demand=ix(cache.baseline(ix,4)>cache.threshold(ix));
    gap=cache.threshold(demand)-run.actual(demand,4);
    rows=[seg.loadKg seg.speedCode max(run.actual(ix,4)-cache.threshold(ix)), ...
        mean(max(gap,0)) max(run.actual(ix,25))];
    fprintf('RELEASE RESULT rate %g [load speed maximum-excess mean-demand-gap assist-peak]\n',rate);
    disp(rows);
    save(fullfile(folder,sprintf('release_%03d.mat',rate)),'run','options','rows','source','group','-v7');
end
