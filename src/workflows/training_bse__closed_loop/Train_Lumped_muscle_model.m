% TRAIN_LUMPED_MUSCLE_MODEL Calibrate model parameters using the selected experiment and optimizer settings.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

load('Combination.mat');

Ts = 0.001;

damping = 0.5;

alpha0 = 0.24086; % obtain from the paper

Combination(:, 4) = zeros(length(Combination(:, 1)), 1)';

for i = 1:length(Combination(:, 4))
    if Combination(i, 2) > 0.02
        Combination(i, 4) = 20;
    end
end

%% using the least square method to obtain the optimal parameters

objective = @(coeffs) Cost_function(Combination(:, 3), Combination(:, 1), Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), 0.7);

%%%%%%%%%%% R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0%%%%%%%%%%%

lb = [0 0 0 0 15 0.2 0.345 0 0 0 -20 -20 -20];
ub = [0.3 1.2 1.2 50000 45 0.35 0.85 1.1 2 100 10 20 20];

A = [];
b = [];

Aeq = [];
beq = [];
nonlcon = [];
initialConditions = [
                     0.0380196602941307 0.647096851480744 0.204039384775119 10553.7289735908 22.1256426143705 0.327435310461100 0.845758789157293 0.772835212459554 1.19616359542222 84.6780523634433];
options = optimoptions('ga', 'MaxGenerations', 7000, 'PopulationSize', 600, 'FunctionTolerance', 1e-12, 'UseParallel', true, 'PlotFcn', @gaplotbestf, 'HybridFcn', @fmincon, 'InitialPopulationMatrix', initialConditions);

coeffs = ga(objective, 10, A, b, Aeq, beq, lb, ub, nonlcon, options);

coeffs = [coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), 0.7];

Tac = zeros(length(Combination(:, 1)), 1)';
Model_torque = Torque_function_constant(Tac, Combination(:, 3), Combination(:, 3), Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));

Torque_Tac0 = Model_torque(:, 1);
Theta_Tac0 = Model_torque(:, 2);

save('Torque_Tac0.mat', 'Torque_Tac0');
save('Theta_Tac0.mat', 'Theta_Tac0');

disp(['The optimal coefficients are: ', num2str(coeffs)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the joint torque obtain by one rebuild muscle

t = linspace(0, (length(Combination(:, 1)) - 1) * 0.001, length(Combination(:, 1)));

t_interp = t';

figure;

plot(t_interp, Combination(:, 1));

hold on;

plot(t_interp, Model_torque(:, 1));

ylabel('L5S1_torque');

legend('Real-torque', 'Model-torque');
title('Real-torque, Model-torque');

figure;
plot(t_interp, Combination(:, 3));

hold on;
plot(t_interp, Model_torque(:, 2));

ylabel('L5/S1 flex-extension angle');

legend('Target-angle', 'Model-angle');
title('Target-angle, Model-angle');
