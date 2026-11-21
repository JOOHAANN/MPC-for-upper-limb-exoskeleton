% Screen source-specific NMPC options without changing the original demo.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
settings=load(fullfile(root,'Demo','NMPC','data','controller_528.mat'));
out=fullfile(root,'Demo','NMPC','output','pid_nmpc_528');
% Source, group, guard scale, guard floor, guard cap, under-weight, headroom.
trials=[2 1 0 1 24 200 .5; 2 2 0 1 24 200 .5; 2 3 0 1 24 200 .5; ...
        1 3 .85 1 24 200 .5; 1 3 1.5 1 40 200 .5; ...
        1 3 2 1 50 200 .5; 1 3 .85 3 24 500 1];
for trial=1:size(trials,1)
    v=trials(trial,:); source=v(1); group=v(2); sequence=subject.sequences{group};
    model=subject.models{source};
    plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
    inputs=bse_inputs(sequence,model,source);
    if source==1, planned=bse_inputs(sequence,model,2); else, planned=inputs; end
    cache=load(fullfile(out,sprintf('nmpc_group%d_source%d.mat',group,source)));
    options=cache.options;
    options.predictionErrorGuardScale=v(3); options.predictionErrorGuardFloorNm=v(4);
    options.predictionErrorGuardCapNm=v(5); options.backTorqueUnderWeight=v(6);
    options.trackingHeadroomNm=v(7); options.recoveryBackTargetNm=0;
    fprintf('\nTRIAL %d | source %d | group %d | guard %.3g floor %.3g cap %.3g under %.3g\n', ...
        trial,source,group,v(3:6));
    run=simulate_nmpc(plant,sequence,inputs,planned,cache.threshold,options);
    rows=[];
    for stage=1:3
        seg=sequence.segments{stage}; ix=seg.startIndex:seg.endIndex;
        demand=ix(cache.baseline(ix,4)>cache.threshold(ix));
        gap=cache.threshold(demand)-run.actual(demand,4);
        rows=[rows;seg.loadKg seg.speedCode max(run.actual(ix,4)-cache.threshold(ix)), ...
            mean(max(gap,0)) max(run.actual(ix,25))]; %#ok<AGROW>
    end
    fprintf('RESULT trial %d [load speed maximum-excess mean-demand-gap assist-peak]\n',trial);
    disp(rows);
    save(fullfile(folder,sprintf('probe_%02d.mat',trial)),'run','options','rows','source','group','-v7');
end
