% LUMPED_MODEL_CHECK_WITH_COEFFS_THESIS_COMBINATION Compare model predictions with the selected experimental reference signals.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

load('Combination.mat');

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

m1 = 15;

m2 = 15;

coeffs = [0.137158002711631 0.509929292299829   0.197540449462042   1829.74428550154    27.1832128387185    0.342633377194440   0.708672958988020   1.09999778331619    1.55838113142714    87.4125776644667    0.700000000000000];

Zero_activation = zeros(length(Combination(:, 1)), 1)';
One_activation = ones(length(Combination(:, 1)), 1)';

% for i = 1:length(Smooth_activation_1)
%     if i < 15000 && i > 5000
%         One_activation(i) = 0.5;
%     end
% end

Tac = zeros(length(Combination(:, 1)), 1)';
Model_torque = Torque_function_constant(Tac, Combination(:, 3), Combination(:, 3), Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));
% Model_torque = Torque_function_constant(m1, Tac, theta_real_1, Torque_real_1, Smooth_activation_1, Smooth_activation_front, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11), coeffs(12));
for i = 1:length(Tac)
    if (i > 8700 && i < 10000) || (i > 11088 && i < 12924)
        Tac(i) = 50;
        % if Combination(i,1) > 100
        %     Tac(i) = 50;
    end
end

Torque_Tac0 = Model_torque(:, 1);
Theta_Tac0 = Model_torque(:, 2);

save('Torque_Tac0.mat', 'Torque_Tac0');
save('Theta_Tac0.mat', 'Theta_Tac0');

Model_torque_2 = Torque_function_constant(m1, Tac, Combination(:, 3), Theta_Tac0, Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11), coeffs(12), coeffs(13), coeffs(14));

Model_torque_3 = Torque_function_constant(m1, Tac, Combination(:, 3), Combination(:, 1), Combination(:, 2), Combination(:, 4), damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11), coeffs(12), coeffs(13), coeffs(14));

k = linspace(0, 100, length(Torque_Tac0));

R_T = corrcoef(Combination(:, 1), Model_torque(:, 1));
r_T = R_T(1, 2);

R_th = corrcoef(Combination(:, 3), Model_torque(:, 2));
r_th = R_th(1, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

plot(k, Combination(:, 1), 'b', 'LineWidth', 4, 'DisplayName', 'Experimental L5-S1 Torque');

text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Movement Cycle (%)', 'FontSize', 24);
ylabel('Torque (NM)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
hold off;

Combination(:, 4) = normalize_max_min(Combination(:, 4), 0.85, 0.7);

% normalize_max_min(Smooth_activation_1(50:150), max(Smooth_activation_1(50:150)), min(Reference_1((9846:11812),2)))
figure;

plot(k(1:34000), Combination(1:34000, 4), 'r--', 'LineWidth', 3, 'DisplayName', 'Normalized Experimental L_{MT}');
hold on;
plot(k(1:34000), Model_torque(1:34000, 4), 'b', 'LineWidth', 4, 'DisplayName', 'Model generated L_{MT}');

text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Movement Cycle (%)', 'FontSize', 24);
ylabel('L_{MT} (mm)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
hold off;

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
