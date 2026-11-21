function Model_output = Torque_function(model, balance_coeffs, sequence, activation, useNative)
%TORQUE_FUNCTION Calculate theta and both muscle torques with the same model.
% Columns: theta, velocity, net torque, back torque, abdominal torque, ...
% The remaining columns contain fiber states, activation and joint balance.
% This routine has no reference-angle or reference-torque tracking input.

%% Create the Hill model with the fixed coefficients
plant = bse_model(model,balance_coeffs,sequence.standingAngle,useNative);

%% Preallocate the output arrays and initialize the Hill state
Ts = sequence.dt;
NumberOfSubsteps = 10;
NumberOfData = size(activation,1);
Model_output = zeros(NumberOfData,27);
[X,Model_output(1,:)] = plant.initialize(activation(1,:),sequence.mass(1), ...
    sequence.waiting(1),sequence.initialState);

%% Advance the same state through every action and standing interval
for k = 2:NumberOfData
    Tac = 0;
    [X,Model_output(k,:)] = plant.step(X,activation(k,:),sequence.mass(k), ...
        sequence.waiting(k),Tac,Ts,NumberOfSubsteps);
end
end
