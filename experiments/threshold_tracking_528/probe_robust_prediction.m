% Screen activation-uncertainty previews on a continuous action prefix.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
source=1; group=3; sequence=subject.sequences{group}; model=subject.models{source};
plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
inputs=bse_inputs(sequence,model,source); planned=bse_inputs(sequence,model,2);
cache=load(fullfile(root,'Demo','NMPC','output','pid_nmpc_528','nmpc_group3_source1.mat'));
uncertainty=load(fullfile(folder,'activation_uncertainty.mat'));
count=1600; sequence.time=sequence.time(1:count);
sequence.mass=sequence.mass(1:count); sequence.waiting=sequence.waiting(1:count);
inputs=inputs(1:count,:); planned=planned(1:count,:);
trials=[8 .25;8 .5;8 1;8 1.5;4 1];
for trial=1:size(trials,1)
    options=cache.options; options.blocks=trials(trial,1);
    options.activationUncertainty=trials(trial,2)*uncertainty.uncertainty;
    options.predictionErrorGuardScale=0; options.predictionErrorGuardFloorNm=1;
    options.backTorqueMarginNm=.5;
    options.backTorqueUnderWeight=200; options.trackingHeadroomNm=.5;
    options.recoveryBackTargetNm=0; options.earlyRecovery=false;
    fprintf('\nROBUST TRIAL %d | blocks %d | uncertainty %.3g\n',trial,trials(trial,:));
    run=simulate_threshold_nmpc(plant,sequence,inputs,planned,cache.threshold(1:count),options);
    seg=sequence.segments{1}; ix=seg.startIndex:seg.endIndex;
    demand=ix(cache.baseline(ix,4)>cache.threshold(ix));
    gap=cache.threshold(demand)-run.actual(demand,4);
    rows=[seg.loadKg seg.speedCode max(run.actual(ix,4)-cache.threshold(ix)), ...
        mean(max(gap,0)) max(run.actual(ix,25))];
    fprintf('ROBUST RESULT trial %d [load speed maximum-excess mean-demand-gap assist-peak]\n',trial);
    disp(rows);
    save(fullfile(folder,sprintf('robust_%02d.mat',trial)),'run','options','rows','source','group','-v7');
end
