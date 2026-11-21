% Inspect forecasts on a replayed state; future observations are diagnostic only.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
subjects=bse_subjects('528'); subject=subjects{1}; useNative=build_bse_solver;
source=1; sequence=subject.sequences{3}; model=subject.models{source};
plant=bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
inputs=bse_inputs(sequence,model,source); planned=bse_inputs(sequence,model,2);
cache=load(fullfile(root,'Demo','NMPC','output','pid_nmpc_528','nmpc_group3_source1.mat'));
[state,~]=plant.initialize(inputs(1,:),sequence.mass(1),sequence.waiting(1),sequence.initialState);
o=cache.options;
u=load(fullfile(folder,'activation_uncertainty.mat'));
for k=2:1060
    [state,~]=plant.step(state,inputs(k,:),sequence.mass(k),sequence.waiting(k), ...
        cache.run.actual(k,25),sequence.dt,o.actualSubsteps);
    if any(k==[1020 1030 1040 1050])
        f=nmpc_forecast(sequence,inputs,planned,k,o);
        z=plant.preview(state,zeros(4,1),f,o);
        rising=plant.preview(state,[0;20;40;40],f,o);
        flat=plant.preview(state,40*ones(4,1),f,o);
        fprintf('\nPrediction at %.2f s; current assist %.2f; threshold %.2f\n', ...
            sequence.time(k),state(29),cache.threshold(k));
        disp([(1:16)'*.025 z(:,4) rising(:,4) flat(:,4)]);
        fprintf('Scenarios: uncertainty multiplier, max zero-assist forecast torque\n');
        for gain=[-2 -1 -.5 .5 1 1.5 2]
            test=f; test.activation=max(0,min(1,f.activation+gain*u.uncertainty.*[1 -1]));
            pred=plant.preview(state,zeros(4,1),test,o);
            fprintf('%g %.3f\n',gain,max(pred(:,4)));
        end
        query=k+(1:16)'*.025/sequence.dt;
        oracle=f; oracle.activation=interp1((1:size(inputs,1))',inputs,query);
        pred=plant.preview(state,zeros(4,1),oracle,o);
        fprintf('Diagnostic future-observation preview peak %.3f\n',max(pred(:,4)));
    end
end
