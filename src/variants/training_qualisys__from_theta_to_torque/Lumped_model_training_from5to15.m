clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

parentFolder = exo_data_root();
childFolder_1 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_45\';
childFolder_2 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_45\';
[Reference_1, Scale_1] = Dynamic_kine_reference(childFolder_1, 153, 267);
[Reference_2, Scale_2] = Dynamic_kine_reference(childFolder_2, 154, 200);

start_stamp = 8260; % 1062 9033 7549
end_stamp = 12862; % 4177、12248 14096
start_stamp_2 = 10342; % 1062 9033
end_stamp_2 = 11146; % 4177

%% \240730Exp\CalibrationData\dynamics_stoop_5kg_45_3\246

Torque_real_1 = Reference_1(:, 1);
Smooth_activation_1 = Reference_1(:, 2);
theta_real_1 = Reference_1(:, 3);
Target_LMT_1 = Reference_1(:, 4);

Torque_real_2 = Reference_2(:, 1);
Smooth_activation_2 = Reference_2(:, 2);
theta_real_2 = Reference_2(:, 3);
Target_LMT_2 = Reference_2(:, 4);

%% %%%%%%%%%%%%%%  edit the activations, theta, and torque to fit the dynamic mode-2 %%%%%%%%%%%

% found = false;  % 用于标记是否已经找到了第一次满足条件的i
% for i = 1:length(theta_real_1)
%     if i <= start_stamp || i >= end_stamp
%         Torque_real_1(i) = 0;
%         Smooth_activation_1(i) = 0;
%         theta_real_1(i) = 0.9;
%     end
% end

Torque_real_1 = smoothdata(Torque_real_1, "gaussian", 300);
Smooth_activation_1 = smoothdata(Smooth_activation_1, "gaussian", 300);
% Torque_real_1(10336:11177) = smoothdata(Torque_real_1(10336:11177), "gaussian", 1200);
% Smooth_activation_1(10336:11177) = smoothdata(Smooth_activation_1(10336:11177), "gaussian", 1200);
Torque_real_2 = smoothdata(Torque_real_2, "gaussian", 300);
Smooth_activation_2 = smoothdata(Smooth_activation_2, "gaussian", 300);
%% Exp240405
%% 240408Youhan
%% \240319Pilot\dynamic_stoop_5kg\
%% \stoop_5kg_3\
%% \240430Youhan\CalibrationData\pilot_noexo_stoop_5_2\
%% \240507Pilot\dynamic_stoop_5kg_3\
%% \240528Pilot\CalibrationData\dynamics_stoop_15kg_45\
%% \240514Pilot\CalibrationData\dynamic_stoop_5_0\
%% 240725Pilot\dynamics_stoop_5kg
%% 240730Exp\CalibrationData\dynamics_stoop_5kg_45_3
Ts = 0.001;

damping = 0.5;

alpha0 = 0.24086; % obtain from the paper

m1 = 5;

m2 = 10;

% theta_real_1(1) = 1.3;

%% using the least square method to obtain the optimal parameters

objective = @(coeffs) Cost_function_from5to15(m1, Torque_real_1, theta_real_1, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4));
% objective_value = Torque_fitness_function(Torque_real(206:758), theta, activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6));
% objective = @(coeffs) Cost_function_from5to15(10, Torque_real_1, theta_real_1, Smooth_activation_1, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8));

%%%%%%%%%%% M, L1, L2, phi%%%%%%%%%%%
% lb = [0.02 0.07 0 0 20 0.3 0.17 0 0];
% ub = [0.2 0.3 1.2 40000 70 1.3 1.3 1.2 1.2];

lb = [20 0.2 0.345 0];
ub = [70 0.35 0.85 1.1];

% lb = [0.02 0 0 0 20 0.3 0.1 0 0];
% ub = [0.5 2 2 30000 80 2 2 1.2 1.2];
% A = [-1 0 1 1 0 0 0];
% b = 0.4;
A = [];
b = [];

Aeq = [];
beq = [];
nonlcon = [];
options = optimoptions('ga', 'MaxGenerations', 7000, 'PopulationSize', 600, 'FunctionTolerance', 1e-12, 'UseParallel', true, 'PlotFcn', @gaplotbestf, 'HybridFcn', @fmincon);
% options = optimoptions('ga','MaxGenerations', 7000,'PopulationSize', 600, 'FunctionTolerance', 1e-12 ,'UseParallel', true,'PlotFcn', @gaplotbestf,'HybridFcn', @fmincon);

coeffs = ga(objective, 4, A, b, Aeq, beq, lb, ub, nonlcon, options);

coeffs = [coeffs(1), coeffs(2), coeffs(3), coeffs(4)];

Zero_activation = zeros(length(Smooth_activation_2), 1)';
Model_torque = Torque_function_constant(m1, Torque_real_1, theta_real_1, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4));

Model_torque_2 = Torque_function_constant(m1, 50, theta_real_1, Torque_real_1, Smooth_activation_1, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));

Model_torque_3 = Torque_function_constant(m1, 0, theta_real_2, Torque_real_2, Zero_activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));

LMT = LMT_function(Smooth_activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the joint torque obtain by one rebuild muscle

figure;

plot(t_interp, Torque_real);

hold on;

plot(t_interp, Model_torque);

ylabel('L5S1_torque');

legend('Real-torque', 'Model-torque');
title('Real-torque, Model-torque');

figure;
plot(t_interp, activation);
hold on;

title('Average activation, EMGfilt-IL-r');

figure;
plot(t_interp, LMT);
title('LMT');
figure;
plot(t_interp, BE);
title('MA');

%% plot the model_LMT and target_LMT
figure;

plot(t_interp, Target_LMT);

hold on;

plot(t_interp, LMT);

ylabel('Muscle tendon length');

legend('Target-LMT', 'Model-LMT');
title('Target-LMT, Model-LMT');
