% Validate one unchanged MAP controller on all three continuous sequences.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
source=2;
for group=1:3
    sequence=subject.sequences{group}; model=subject.models{source};
    plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
    inputs=bse_inputs(sequence,model,source); planned=inputs;
    cache=load(fullfile(root,'Demo','NMPC','output','pid_nmpc_528', ...
        sprintf('nmpc_group%d_source%d.mat',group,source)));
    options=cache.options; options.blocks=8;
    options.activationUncertainty=zeros(16,2); options.activationEdgeTimesSeconds=[];
    options.predictionErrorGuardScale=0; options.predictionErrorGuardFloorNm=2.5;
    options.backTorqueMarginNm=.5;
    options.backTorqueUnderWeight=200; options.trackingHeadroomNm=.5;
    options.recoveryBackTargetNm=0; options.earlyRecovery=false;
    fprintf('\nMAP TRACKING group %d\n',group);
    run=simulate_threshold_nmpc(plant,sequence,inputs,planned,cache.threshold,options);
    rows=[];
    for stage=1:3
        seg=sequence.segments{stage}; ix=seg.startIndex:seg.endIndex;
        demand=ix(cache.baseline(ix,4)>cache.threshold(ix));
        gap=cache.threshold(demand)-run.actual(demand,4);
        rows=[rows;seg.loadKg seg.speedCode max(run.actual(ix,4)-cache.threshold(ix)), ...
            mean(max(gap,0)) max(run.actual(ix,25))]; %#ok<AGROW>
    end
    fprintf('MAP RESULT group %d [load speed maximum-excess mean-demand-gap assist-peak]\n',group);
    disp(rows);
    save(fullfile(folder,sprintf('map_tracking_group%d.mat',group)),'run','options','rows','source','group','-v7');
end
