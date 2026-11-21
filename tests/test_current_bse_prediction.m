function test_current_bse_prediction
%TEST_CURRENT_BSE_PREDICTION Match the applied plant to its nonlinear preview.
root=fileparts(fileparts(mfilename('fullpath')));
bse=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');addpath(bse,fullfile(root,'Demo','NMPC'));
data=load(fullfile(bse,'data','subject_528.mat'));subject=data.subject;
settings=current_528_settings(subject);
for source=1:2
    q=subject.sequences{1};inputs=prepare_activation(q,subject.models{source},source);
    planned=prepare_activation(q,subject.models{2},2);o=settings.controllers{source};
    plant=current_bse_plant(subject.models{source},subject.balanceController,q.standingAngle,settings.humanPD,settings.actuator);
    [state,~]=plant.initialize(inputs(1,:),q.mass(1),q.waiting(1),q.initialState);
    current=q.segments{1}.startIndex+40;
    for k=2:current
        [state,~]=plant.step(state,inputs(k,:),q.mass(k),q.waiting(k),0,q.dt,10);
    end
    forecast=nmpc_forecast(q,inputs,planned,current,o);
    changed=q;changed.referenceTheta(:)=100;changed.referenceTorque(:)=-100;
    future=inputs;future(current+1:end,:)=1-future(current+1:end,:);
    same=nmpc_forecast(changed,future,planned,current,o);
    assert(isequaln(forecast,same),'The forecast reads future experimental inputs or reference labels.');
    controls=[5;10;15;20;20;10;5;0];
    expected=plant.preview(state,controls,forecast,o);assert(~isempty(expected));
    actual=zeros(o.points,27);
    for k=1:o.points
        block=floor((k-1)*numel(controls)/o.points)+1;
        [state,actual(k,:)]=plant.step(state,forecast.activation(k,:),forecast.mass(k), ...
            forecast.waiting(k),controls(block),o.predictionDt,o.predictionSubsteps);
    end
    assert(max(abs(actual-expected),[],'all')<1e-9,'Prediction and applied plant use different dynamics.');
end
fprintf('PASS: nonlinear preview/applied-plant parity and no future experimental EMG/IK/torque reference in online prediction.\n');
end
