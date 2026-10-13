% MAP_OUTPUT_MODEL_FUSION Combine predicted muscle activation with the reference movement data.
% Use prepare_workflow to select matching code and input data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAP_output_model_fusion.m
%%%
%%% Combine every code together and derive and plot the MAP predictor
%%%
%%% Input: pv_coefficients.mat and pw_coefficients.mat
%%% Output: Combination activation from MAP
%%%
%%% Youhan Li (2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

set(0, 'DefaultFigureWindowStyle', 'docked');

%%%%%%%%%%%%% step 1%%%%%%%%%%%%%% Import and select the segment of the data
parentFolder = exo_data_root();
childFolder_1 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_55\';
childFolder_2 = '\240528Pilot\CalibrationData\dynamics_stoop_10kg_45\';
childFolder_3 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_35\';
[Reference_1, Scale_1] = Dynamic_kine_reference(childFolder_1, 1, 300);
[Reference_2, Scale_2] = Dynamic_kine_reference(childFolder_2, 1, 270);
[Reference_3, Scale_3] = Dynamic_kine_reference(childFolder_3, 50, 300);
%%%%%%%%%%%%% step 1%%%%%%%%%%%%%% Import and select the segment of the data

load('baseline_data.mat');
load('EMG_matrix.mat');
load('Combination_raw.mat');
load('Baselines.mat');

speed = [25, 35, 45, 55, 65];
weight = [5, 7.5, 10, 12.5, 15];

Ts = 0.001;

damping = 0.5;

alpha0 = 0.24086; % obtain from the paper

coeffs = [0.0380196562634957    0.645017822853626   0.204039417115808   10557.7108614069    22.1256427169235    0.329388435617747   0.789068624062216   0.610105856389750   1.19953597779333    86.7539190477979    0.700000000000000   -1.99676909726852   1.03495249411787    1.50672917083905];

baseline_mode = 'SM';
muscle_titles = {'IL_r', 'IL_l', 'LTPT_r', 'LTPT_l', 'LTPL_r', 'LTPL_l'};

[~, ~, MEP_1] = MAP_model(speed(4), weight(1), baseline_mode, 'baseline_data.mat');
[~, ~, MEP_2] = MAP_model(speed(3), weight(3), baseline_mode, 'baseline_data.mat');
[~, ~, MEP_3] = MAP_model(speed(2), weight(5), baseline_mode, 'baseline_data.mat');

zero_activation = zeros(20000, 1);

k = 1:200;
t_interp_1 = linspace(1, 200, 4265);
t_interp_2 = linspace(1, 200, 5018);
t_interp_3 = linspace(1, 200, 5867);

Smooth_activation_1 = (MEP_1(:, 1) + MEP_1(:, 2)) * (1 / 12) + (MEP_1(:, 3) + MEP_1(:, 4)) * (1 / 4)  + (MEP_1(:, 1) + MEP_1(:, 2)) * (1 / 6);
Smooth_activation_2 = (MEP_2(:, 1) + MEP_2(:, 2)) * (1 / 12) + (MEP_2(:, 3) + MEP_2(:, 4)) * (1 / 4)  + (MEP_2(:, 1) + MEP_2(:, 2)) * (1 / 6);
Smooth_activation_3 = (MEP_3(:, 1) + MEP_3(:, 2)) * (1 / 12) + (MEP_3(:, 3) + MEP_3(:, 4)) * (1 / 4)  + (MEP_3(:, 1) + MEP_3(:, 2)) * (1 / 6);

% Smooth_activation_1 = normalize_max_min(Smooth_activation_1, max(Reference_1(:,2)), min(Reference_1(:,2)));
% Smooth_activation_2 = normalize_max_min(Smooth_activation_2, max(Reference_2(1:14000,2)), min(Reference_2(:,2)));
% Smooth_activation_3 = normalize_max_min(Smooth_activation_3, max(Reference_3(:,2)), min(Reference_3(:,2)));

Smooth_activation_1 = normalize_max_min(Smooth_activation_1, 0.40, min(Reference_1(:, 2)));
Smooth_activation_2 = normalize_max_min(Smooth_activation_2, 0.54, min(Reference_2(:, 2)));
Smooth_activation_3 = normalize_max_min(Smooth_activation_3, 0.68, min(Reference_3(:, 2)));

Smooth_activation_1(50:150) = normalize_max_min(Smooth_activation_1(50:150), max(Smooth_activation_1(50:150)), min(Reference_1(9527:10816, 2)));
Smooth_activation_2(50:150) = normalize_max_min(Smooth_activation_2(50:150), max(Smooth_activation_2(50:150)), min(Reference_2(10086:11973, 2)));
Smooth_activation_3(50:150) = normalize_max_min(Smooth_activation_3(50:150), max(Smooth_activation_3(50:150)), min(Reference_3(7453:10250, 2)));
% Smooth_activation_3(50:150) = normalize_max_min(Smooth_activation_3(50:150), max(Smooth_activation_3(50:150)), 0.25);

Smooth_activation_1 = interp1(k', Smooth_activation_1, t_interp_1');
Smooth_activation_2 = interp1(k', Smooth_activation_2, t_interp_2');
Smooth_activation_3 = interp1(k', Smooth_activation_3, t_interp_3');

Smooth_activation_sew = [zero_activation(1:8102); Smooth_activation_1; zero_activation(1:11133); Smooth_activation_2; zero_activation(1:9412); Smooth_activation_3; zero_activation(1:1049)];

% Smooth_activation_sew = smoothdata(Smooth_activation_sew, "gaussian", 10);

Tac = zeros(length(Smooth_activation_sew(:, 1)), 1)';
% Model_torque = Torque_function_constant(Tac, Combination(:, 3), Combination(:, 3), Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));
Combination(:, 2) = Smooth_activation_sew;
Combination(:, 4) = Combination(:, 1);
Combination(:, 1) = 310 * Smooth_activation_sew;
save('Combination.mat', 'Combination');
