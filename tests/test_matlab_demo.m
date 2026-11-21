function test_matlab_demo()
%TEST_MATLAB_DEMO Check shared mechanics and causal activation prediction.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
build_bse_solver;
subjects = bse_subjects;
maximumDifference = 0;
for k = 1:numel(subjects)
    subject = subjects{k};
    sequence = subject.exampleSequence;
    for source = 1:2
        model = subject.models{source};
        inputs = bse_inputs(sequence,model,source);
        native = bse_model(model,subject.balanceController,sequence.standingAngle,true);
        matlab = bse_model(model,subject.balanceController,sequence.standingAngle,false);
        % The readable procedural loop must reproduce the native full replay.
        loopOutput = Torque_function(model,subject.balanceController,sequence,inputs,true);
        rolloutOutput = native.simulate(inputs,sequence.mass,sequence.waiting, ...
            sequence.initialState,sequence.dt,10);
        assert(max(abs(loopOutput(:)-rolloutOutput(:))) < 1e-8, ...
            'The procedural BSE loop changed the original simulation.');
        [state,~] = native.initialize(inputs(1,:),sequence.mass(1),sequence.waiting(1),sequence.initialState);
        first = double(sequence.segments{1}.startIndex);
        last = min(size(inputs,1),first+round(1.8/sequence.dt));
        for j = 2:last
            [state,~] = native.step(state,inputs(j,:),sequence.mass(j),sequence.waiting(j),0,sequence.dt,10);
        end
        nativeState = state;
        matlabState = state;
        for j = last+(1:20)
            assistance = min(20,2*(j-last));
            [nativeState,a] = native.step(nativeState,inputs(j,:),sequence.mass(j),sequence.waiting(j),assistance,sequence.dt,10);
            [matlabState,b] = matlab.step(matlabState,inputs(j,:),sequence.mass(j),sequence.waiting(j),assistance,sequence.dt,10);
            maximumDifference = max(maximumDifference,max(abs(a-b)));
            [difference,column] = max(abs(a-b));
            assert(difference < 1e-6,'MATLAB/native mismatch: %s source %g sample %g column %g difference %.9g.', ...
                subject.id,source,j,column,difference);
            assert(abs(a(3)-a(4)+a(5)) < 1e-10,'Muscle torque accounting failed.');
        end
    end
    settings = load(fullfile(root,'Demo','NMPC','data',['controller_' subject.shortName '.mat']));
    options = settings.controller;
    observed = bse_inputs(sequence,subject.models{1},1);
    planned = bse_inputs(sequence,subject.models{1},2);
    current = double(sequence.segments{1}.startIndex)+80;
    baseline = nmpc_forecast(sequence,observed,planned,current,options);
    observed(current+1:end,:) = .99;
    altered = nmpc_forecast(sequence,observed,planned,current,options);
    assert(isequal(baseline,altered),'Forecast reads future experimental activation.');
end
fprintf('MATLAB/native motion and assistance checks passed; maximum difference %.3g.\n',maximumDifference);
fprintf('Causal forecast checks passed for all four subjects.\n');
end
