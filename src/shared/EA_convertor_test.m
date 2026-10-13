clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

Ts = 0.001;

damping = 0.5;

alpha0 = 0.24086; % obtain from the paper

coeffs = [0.0380196562634957    0.645017822853626   0.204039417115808   10557.7108614069    22.1256427169235    0.329388435617747   0.789068624062216   0.610105856389750   1.19953597779333    86.7539190477979    0.700000000000000   -1.99676909726852   1.03495249411787    1.50672917083905];

load('Combination_raw.mat');

acitvation_raw = Combination(:, 2);
Torque_raw = Combination(:, 1);

load('Combination.mat');

excitation = Combination(:, 2);

theta = Combination(:, 3);

activation_EA = ExcAct(Combination(:, 2), -0.1);

activation_exp_EA = ExcAct(acitvation_raw, -0.1);

k = linspace(0, 100, length(activation_EA));

Combination(:, 2) = activation_EA;

save('Combination_EA.mat', 'Combination');

Tac = zeros(length(activation_EA), 1)';

Model_torque = Torque_function_constant(0, Tac, theta, theta, excitation, excitation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11), coeffs(12), coeffs(13), coeffs(14));

figure;

plot(k, acitvation_raw, 'r--', 'LineWidth', 2, 'DisplayName', 'Experimental Excitation');
hold on;
plot(k, excitation, 'g--', 'LineWidth', 2, 'DisplayName', 'MEP Generated Excitation');
plot(k, activation_EA, 'b', 'LineWidth', 3, 'DisplayName', 'MAP Activation');

text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Movement Cycle (%)', 'FontSize', 24);
ylabel('Amplitude (0-1)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
hold off;

%%%%%%%%%%%%%%%%%%% Plot the Kinematic properties of the Subject A%%%%%%%%%%%%%%%%%%

figure;

plot(k, theta, 'b', 'LineWidth', 2, 'DisplayName', 'Experimental Theta');

text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Movement Cycle (%)', 'FontSize', 24);
ylabel('Theta (rad)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
hold off;

%%%%%%%%%%%%%%%%%%% Plot the Kinematic properties of the Subject A%%%%%%%%%%%%%%%%%%
k = linspace(0, 100, length(activation_exp_EA(7598:13944)));
figure;

subplot(2, 1, 1);

plot(k, Torque_raw(7598:13944), 'b', 'LineWidth', 5, 'DisplayName', 'Experimental L5-S1 torque');

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

subplot(2, 1, 2);

plot(k, theta(7598:13944), 'b', 'LineWidth', 5, 'DisplayName', 'Experimental \theta');

text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Movement Cycle (%)', 'FontSize', 24);
ylabel('\theta (rad)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
hold off;
