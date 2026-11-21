% Screen move blocking on the difficult continuous experimental sequence.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
source=1; group=3; sequence=subject.sequences{group}; model=subject.models{source};
plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
inputs=bse_inputs(sequence,model,source); planned=bse_inputs(sequence,model,2);
cache=load(fullfile(root,'Demo','NMPC','output','pid_nmpc_528','nmpc_group3_source1.mat'));
trials=[1 .3 1;2 .3 1;8 .3 1;16 .3 1;1 0 2;2 0 2];
for trial=1:size(trials,1)
    options=cache.options; options.blocks=trials(trial,1);
    options.predictionErrorGuardScale=trials(trial,2);
    options.predictionErrorGuardFloorNm=trials(trial,3);
    options.backTorqueUnderWeight=200; options.trackingHeadroomNm=.5;
    options.recoveryBackTargetNm=0; options.earlyRecovery=false;
    fprintf('\nBLOCK TRIAL %d | blocks %d | guard %.3g | floor %.3g\n',trial,trials(trial,:));
    run=simulate_nmpc(plant,sequence,inputs,planned,cache.threshold,options);
    rows=[];
    for stage=1:3
        seg=sequence.segments{stage}; ix=seg.startIndex:seg.endIndex;
        demand=ix(cache.baseline(ix,4)>cache.threshold(ix));
        gap=cache.threshold(demand)-run.actual(demand,4);
        rows=[rows;seg.loadKg seg.speedCode max(run.actual(ix,4)-cache.threshold(ix)), ...
            mean(max(gap,0)) max(run.actual(ix,25))]; %#ok<AGROW>
    end
    fprintf('BLOCK RESULT trial %d [load speed maximum-excess mean-demand-gap assist-peak]\n',trial);
    disp(rows);
    save(fullfile(folder,sprintf('blocks_%02d.mat',trial)),'run','options','rows','source','group','-v7');
end
