% Preserve MAP preview while improving short-term residual continuity.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
out=fullfile(root,'Demo','NMPC','output','pid_nmpc_528');
% Source, group, local blend, blend decay, guard scale, guard floor.
trials=[1 3 1 .15 .3 1; 1 3 1 .15 .85 1; 1 3 .5 .15 .3 1; ...
        1 3 1 .30 .3 1; 1 3 1 .30 0 2; ...
        2 1 0 .15 0 2.5; 2 2 0 .15 0 2.5; 2 3 0 .15 0 2.5];
for trial=1:size(trials,1)
    v=trials(trial,:); source=v(1); group=v(2); sequence=subject.sequences{group};
    model=subject.models{source};
    plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
    inputs=bse_inputs(sequence,model,source);
    if source==1, planned=bse_inputs(sequence,model,2); else, planned=inputs; end
    cache=load(fullfile(out,sprintf('nmpc_group%d_source%d.mat',group,source)));
    options=cache.options; pred=options.activationResidualPredictor;
    local=zeros(size(pred.coefficients)); lag=6; row=find(pred.lags==lag);
    times=(1:options.points)*options.predictionDt;
    for point=1:options.points
        lead=.12*(1-exp(-times(point)/.12)); gain=lead/(lag*sequence.dt);
        for channel=1:2
            col=2*(point-1)+channel;
            local(channel,col)=pred.scale(channel)*(1+gain);
            offset=2*(row-1)+channel;
            local(offset,col)=-pred.scale(offset)*gain;
        end
    end
    blend=repelem(v(3)*exp(-times/v(4)),2);
    pred.coefficients=pred.coefficients.*(1-blend)+local.*blend;
    options.activationResidualPredictor=pred;
    options.predictionErrorGuardScale=v(5); options.predictionErrorGuardFloorNm=v(6);
    options.backTorqueUnderWeight=200; options.trackingHeadroomNm=.5;
    options.recoveryBackTargetNm=0;
    fprintf('\nLOCAL TRIAL %d | source %d | group %d | blend %.3g tau %.3g guard %.3g floor %.3g\n', ...
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
    fprintf('LOCAL RESULT trial %d [load speed maximum-excess mean-demand-gap assist-peak]\n',trial);
    disp(rows);
    save(fullfile(folder,sprintf('local_%02d.mat',trial)),'run','options','rows','source','group','-v7');
end
