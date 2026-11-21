% Causal curvature correction for a fast abdominal-activation rise.
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
trials=[1 .2 200;1 .3 200;.5 .2 200;1 .2 100;1 .2 60;1 .2 40];
for trial=1:size(trials,1)
    options=cache.options; options.blocks=8;
    pred=options.activationResidualPredictor; local=pred.coefficients;
    lag=4; row1=find(pred.lags==lag); row2=find(pred.lags==2*lag); h=lag*sequence.dt;
    for point=1:options.points
        lead=min(point*options.predictionDt,trials(trial,2)); ratio=lead/h;
        col=2*point; local(:,col)=0;
        rows=[2,2*row1,2*row2]; coefficients=[1+1.5*ratio+.5*ratio^2,-2*ratio-ratio^2,.5*ratio+.5*ratio^2];
        scale=pred.scale(:); local(rows,col)=scale(rows).*coefficients(:);
    end
    pred.coefficients=(1-trials(trial,1))*pred.coefficients+trials(trial,1)*local;
    options.activationResidualPredictor=pred;
    options.predictionErrorGuardScale=.3; options.predictionErrorGuardFloorNm=1;
    options.backTorqueMarginNm=.5; options.releaseRateNmPerSec=trials(trial,3);
    options.backTorqueUnderWeight=200; options.trackingHeadroomNm=.5;
    options.recoveryBackTargetNm=0; options.earlyRecovery=false;
    fprintf('\nABDOMEN TRIAL %d | blend %g | lead cap %g | release %g\n',trial,trials(trial,:));
    run=simulate_threshold_nmpc(plant,sequence,inputs,planned,cache.threshold(1:count),options);
    seg=sequence.segments{1}; ix=seg.startIndex:seg.endIndex;
    demand=ix(cache.baseline(ix,4)>cache.threshold(ix));
    gap=cache.threshold(demand)-run.actual(demand,4);
    rows=[seg.loadKg seg.speedCode max(run.actual(ix,4)-cache.threshold(ix)), ...
        mean(max(gap,0)) max(run.actual(ix,25))];
    fprintf('ABDOMEN RESULT trial %d [load speed maximum-excess mean-demand-gap assist-peak]\n',trial);
    disp(rows);
    save(fullfile(folder,sprintf('abdomen_%02d.mat',trial)),'run','options','rows','source','group','-v7');
end
