function test_bse_multisignal()
%TEST_BSE_MULTISIGNAL Compare both solvers with all added intent terms active.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Demo','BSE'));
build_bse_solver;
subjects = bse_subjects;
maximumDifference = 0;
for k = 1:numel(subjects)
    subject = subjects{k}; sequence = subject.exampleSequence;
    index = double(sequence.segments{1}.startIndex)+80;
    for source = 1:2
        model = subject.models{source};
        model.neuralParameters = [model.neuralParameters(1:10) .25 -.15 .08 -.12 .20 .10];
        inputs = bse_inputs(sequence,model,source);
        native = bse_model(model,subject.balanceController,sequence.standingAngle,true);
        matlab = bse_model(model,subject.balanceController,sequence.standingAngle,false);
        [a,x] = native.initialize(inputs(index,:),sequence.mass(index),0,[1.1 0]);
        [b,y] = matlab.initialize(inputs(index,:),sequence.mass(index),0,[1.1 0]);
        assert(max(abs(a-b)) < 1e-7 && max(abs(x-y)) < 1e-7);
        for j = index+(1:20)
            [a,x] = native.step(a,inputs(j,:),sequence.mass(j),0,10,sequence.dt,10);
            [b,y] = matlab.step(b,inputs(j,:),sequence.mass(j),0,10,sequence.dt,10);
            maximumDifference = max(maximumDifference,max(abs(x-y)));
        end
    end
end
assert(maximumDifference < 1e-6,'Multisignal MATLAB/native parity failed.');
fprintf('All six intent terms active: maximum MATLAB/native output difference %.3g.\n',maximumDifference);
end
