% TRAIN_LUMPED_MUSCLE_MODEL_EMG_AND_LMT_TO_MTF Calibrate model parameters using the selected experiment and optimizer settings.
% Use prepare_workflow to select matching code and input data.
%%%%%%%% --------Program Description--------%%%%%%%%
% with predefined Geometric parameters M, L1, L2, phi and LMT to fit the reference torque
%%%%%%%% --------Program Description--------%%%%%%%%
%%%%%%%%%%% R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0%%%%%%%%%%%
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

load('Combination_5_35_10_55_15_45.mat', 'Combination');
Combination_354 = Combination;

load('Combination_5_45_10_35_15_55.mat', 'Combination');
Combination_435 = Combination;

load('Combination_5_55_10_45_15_35.mat', 'Combination');
Combination_543 = Combination;

Combination = [Combination_354; Combination_435; Combination_543];

Ts = 0.001;

damping = 0.5;

alpha0 = 0.24086; % obtain from the paper

% for i = 1:length(Combination(:,4))
%     if Combination(i,2) > 0.02
%         Combination(i,4) = 20;
%     end
% end

%% using the least square method to obtain the optimal parameters

%% 22.9819        0.345     0.756835     0.950692      92.7424          0.7
%% 24.3502        0.345     0.735338     0.958016      83.5524          0.7
%% 32.2981        0.345     0.628236     0.880302      86.0247          0.7

M = 24.3502;
L1 = 0.345;
L2 = 0.735338;
phi = 0.958016;

objective = @(coeff) Cost_function(Combination(:, 3), Combination(:, 1), Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeff(1), coeff(2), coeff(3), coeff(4), M, L1, L2, phi, coeff(5), coeff(6), 0.7);

%%%%%%%%%%% R, LoM_nominal, LsT, FoM, theta0, Torque0%%%%%%%%%%%

%%%%%%%%%%% 0.0211600119837552   0.768787457631810   0.146635739203492   24650.4840294265    0.998343837940520   37.0769468313711%%%%%%%%%%%

lb = [0.02 0.52 0.1 7000 0 0];
ub = [0.04 0.72 0.24 17000 2 100];

A = [];
b = [];

Aeq = [];
beq = [];
nonlcon = [];
initialConditions = [
                     0.0380196602941307 0.647096851480744 0.204039384775119 10553.7289735908 0 0];
options = optimoptions('ga', 'MaxGenerations', 7000, 'PopulationSize', 600, 'FunctionTolerance', 1e-12, 'UseParallel', true, 'PlotFcn', @gaplotbestf, 'HybridFcn', @fmincon, 'InitialPopulationMatrix', initialConditions);

coeff = ga(objective, 6, A, b, Aeq, beq, lb, ub, nonlcon, options);

coeffs = [coeff(1), coeff(2), coeff(3), coeff(4), M, L1, L2, phi, coeff(5), coeff(6), 0.7];

save('coeffs.mat', 'coeffs');

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

plot(t_interp, Combination(:, 4));

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
