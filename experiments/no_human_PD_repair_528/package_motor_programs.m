% Separate the calibrated nominal motor program from online posture feedback.
% The frozen source-specific unassisted drive is replayed as feedforward;
% future experimental drive is never supplied to the predictive controller.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','NMPC');
addpath(fullfile(root,'Demo','BSE'));
data=load(fullfile(demo,'output','pid_nmpc_528','threshold_tracking','pid_nmpc_comparison_528.mat'));
subjects=bse_subjects('528'); subject=subjects{1};
programs=cell(3,2);
plannedPrograms=cell(3,2);calibrationParameters=cell(1,2);
for group=1:3
    for source=1:2
        programs{group,source}=data.results.groups{group}.runs{source}.baseline(:,18:19);
        model=subject.models{source};q=subject.sequences{group};
        calibrationParameters{source}=[model.parameters(:);model.neuralParameters(:)];
        if source==2
            plannedPrograms{group,source}=programs{group,source};
        else
            % Predicted MAP activation must pass through this route's own
            % calibration, not the differently calibrated MAP plant.
            planner=bse_model(model,subject.balanceController,q.standingAngle,true);
            mapInput=bse_inputs(q,model,2);
            nominal=planner.simulate(mapInput,q.mass,q.waiting,q.initialState,q.dt,10);
            plannedPrograms{group,source}=nominal(:,18:19);
        end
    end
end
% A causal ARX model is refit in the new excitation-command space.
% Its offline labels may be future residuals; online features are past-only.
lags=[0 1 2 4 8 16 32]; points=100; predictionDt=.01;
residuals=cell(1,3);
for group=1:3, residuals{group}=programs{group,1}-plannedPrograms{group,1}; end
allResidual=vertcat(residuals{:}); scale=max(std(allResidual,0,1),1e-4);
X=[]; Y=[];
for group=1:3
    q=subject.sequences{group}; residual=residuals{group}; n=size(residual,1);
    for current=33:2:n-100
        history=current-lags;
        feature=reshape(residual(history,:)',1,[])./repmat(scale,1,numel(lags));
        future=interp1((1:n)',residual,current+(1:points)'*predictionDt/q.dt,'linear');
        X=[X;feature]; Y=[Y;reshape(future',1,[])]; %#ok<AGROW>
    end
end
beta=(X'*X+20*eye(size(X,2)))\(X'*Y);
predictors=cell(1,2);
for source=1:2
    coefficients=beta; if source==2, coefficients(:)=0; end
    predictors{source}=struct('lags',lags,'scale',scale,'coefficients',coefficients, ...
        'points',points,'predictionDt',predictionDt);
end
modelSHA256=subject.modelSHA256;
provenance=['Source-specific calibrated unassisted excitation commands, preserved as a frozen feedforward motor program. ' ...
    'They are generated offline from the original calibrated experimental/MAP model; no assisted posture error changes this program. ' ...
    'MAP forecasts are transformed through the same route-specific calibration as the actual plant input. ' ...
    'ARX command-residual prediction is calibrated on all nine recorded cases, not an independent test set.'];
save(fullfile(demo,'data','no_human_PD_motor_programs_528.mat'), ...
    'programs','plannedPrograms','calibrationParameters','predictors','modelSHA256','provenance','-v7');
