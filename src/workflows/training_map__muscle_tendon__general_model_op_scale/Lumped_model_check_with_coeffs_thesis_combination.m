% LUMPED_MODEL_CHECK_WITH_COEFFS_THESIS_COMBINATION Compare model predictions with the selected experimental reference signals.
% Use prepare_workflow to select matching code and input data.
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

% load('Combination_total.mat', 'Combination');
save('Combination_total.mat', 'Combination');

start_stamp = 6861; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg) 7532
end_stamp = 11737; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg) 11203
start_stamp_2 = 17033; % 1062 9033 17033
end_stamp_2 = 19809; % 4177、12248 19809

%% \240730Exp\CalibrationData\dynamics_stoop_5kg_45_3\246

%% %%%%%%%%%%%%%%  edit the activations, theta, and torque to fit the dynamic mode-2 %%%%%%%%%%%

%% Exp240405
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
%%% 22.6006        0.345     0.760606     0.936847      99.2107          0.7
%%% 354 0.0359648470384669   0.603618539025773   0.125422132401744   8596.00550568983    24.3502000000000    0.345000000000000   0.735338000000000   0.958016000000000   0.440229814274935   50.0683306850439    0.700000000000000
%%% 435 0.0380196672     0.661745096       0.2040389      10616.9823         24.3502           0.345        0.735338        0.958016     0.415039394          14.625             0.7
%%% 543 0.0278391986031940   0.650763116340046   0.195954593548366   15065.7095159712    24.3502000000000    0.345000000000000   0.735338000000000   0.958016000000000   0.978322545814145   3.51623129905564    0.700000000000000
coeffs = [0.0328531574526041    0.603633931100612   0.184185669889918   11039.3040599535    22.6006000000000    0.345000000000000   0.760606000000000   0.936847000000000   0.360148103553726   49.8813833223633    0.0605843948558919  0.700000000000000];

Zero_activation = zeros(length(Combination(:, 1)), 1)';
One_activation = ones(length(Combination(:, 1)), 1)';

Tac = zeros(length(Combination(:, 1)), 1)';
Model_torque = Torque_function_constant(Tac, Combination(:, 3), Combination(:, 3), Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11), coeffs(12));

Torque_Tac0 = Model_torque(:, 1);
Theta_Tac0 = Model_torque(:, 2);

save('Torque_Tac0.mat', 'Torque_Tac0');
save('Theta_Tac0.mat', 'Theta_Tac0');

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

% ==== Correlation for torque ====
C_torque = corrcoef(Combination(:, 4), Model_torque(:, 1));
R_torque = C_torque(1, 2);
disp(['Correlation (Torque) R = ', num2str(R_torque)]);

figure;
plot(t_interp, Combination(:, 3));
hold on;
plot(t_interp, Model_torque(:, 2));

ylabel('L5/S1 flex-extension angle');
legend('Target-angle', 'Model-angle');
title('Target-angle, Model-angle');

% ==== Correlation for angle ====
C_angle = corrcoef(Combination(:, 3), Model_torque(:, 2));
R_angle = C_angle(1, 2);
disp(['Correlation (Angle) R = ', num2str(R_angle)]);
