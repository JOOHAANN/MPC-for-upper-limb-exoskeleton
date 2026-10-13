% OFFLINE_SIMULATOR_NMPC_QUASIMODO_INCREMENTAL_MPC_WITH_MAP Run the offline NMPC simulation for the selected activation source.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

load('Combination.mat');

Torque_Ceinms = Combination(:, 1);
activation_sew = Combination(:, 2);
theta_Ceinms = Combination(:, 3);

%% %%%%%%%%%%%%%%  edit the activations, theta, and torque to fit the dynamic mode-2 %%%%%%%%%%%

%%%%%%%%% start with designing the NMPC %%%%%%%%%%%%%%
%% scale muscle model

coeffs = [0.0380196562634957    0.645017822853626   0.204039417115808   10557.7108614069    22.1256427169235    0.329388435617747   0.789068624062216   0.610105856389750   1.19953597779333    86.7539190477979    0.700000000000000   -1.99676909726852   1.03495249411787    1.50672917083905];
%%% R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0,
%%% front_a, front_b, front_c

R = coeffs(1);
LoM_nominal = coeffs(2);
LsT = coeffs(3);
FoM = coeffs(4);
M = coeffs(5);
L1 = coeffs(6);
L2 = coeffs(7);
phi = coeffs(8);
theta0 = coeffs(9);
Torque0 = coeffs(10);
LMT0 = coeffs(11);
front_a = coeffs(12);
front_b = coeffs(13);
front_c = coeffs(14);

%% real-time NMPC

%% Preallocations

NumberOfData = length(activation_sew);

theta = zeros(NumberOfData, 1)';
dtheta = zeros(NumberOfData, 1)';
ddtheta = zeros(NumberOfData, 1)';

LoM = zeros(NumberOfData, 1)';
alpha = zeros(NumberOfData, 1)';
dalpha = zeros(NumberOfData, 1)';
V_max = zeros(NumberOfData, 1)';
LMT = zeros(NumberOfData, 1)';
MA = R * ones(NumberOfData, 1)';
LM = zeros(NumberOfData, 1)';
dLM = zeros(NumberOfData, 1)';
LT = zeros(NumberOfData, 1)';
dLMT = zeros(NumberOfData, 1)';
k_t = zeros(NumberOfData, 1)';
stiffness = zeros(NumberOfData, 1)';
F_L = zeros(NumberOfData, 1)';
F_PE = zeros(NumberOfData, 1)';
t = zeros(NumberOfData, 1)';
F_k_up = zeros(NumberOfData, 1)';
T_balance = zeros(NumberOfData, 1)';
T_front = zeros(NumberOfData, 1)';
t = zeros(NumberOfData, 1)';

MTF = zeros(NumberOfData + 1, 1)';
Tac = zeros(NumberOfData + 1, 1)';
U = zeros(NumberOfData + 1, 1)';
Y = zeros(NumberOfData + 1, 1)';
X = zeros(3, NumberOfData + 1);
dFt = zeros(NumberOfData, 1)';
% Tac = zeros(NumberOfData,1)';
strain = zeros(NumberOfData, 1)';

%% initial the state

PercentageChange = 0.15;

Ts = 0.001;

g = 9.81;

a = 0.9;

kt_spline = [6.07450471822638e+42   -1.17093401235721e+42   9.78861307802160e+40    -4.68161016668128e+39   1.41111335336212e+38    -2.75023728531764e+36   3.32535388381590e+34    -1.97635568280683e+32   -4.74852712911846e+29   1.68812705982292e+28    -9.30681300424510e+25   -1.69227796153188e+23   3.64048821058250e+21    -7.76789223923884e+18   -4.91004587706162e+16   201470178531589 239829087287.116    -1720265398.11745   -15739.3275509684   6708.90359529906    3.90612013684262];

damping = 0.5;

alpha0 = 0.24086;

w = LoM_nominal * sin(alpha0);
c = 10;
Kp = 500; % 200
Kd = 200; % 80
Ki = 0;
InT = 0;
P = 1;
D = 0.6;

PH = 100; % prediction horizon

CH = 5; % control horizon

n = 10; % block Length

Threshold = 120 * ones(NumberOfData, 1);

u_max = 80 * ones(NumberOfData, 1);

InitialConditions0 = [theta0, Torque0, LMT0];

theta(1) = InitialConditions0(1); % state x1

dtheta(1) = 0;

MTF(1) = InitialConditions0(2) / R; % state x2

LMT(1) = InitialConditions0(3); % state x3

dLMT(1) = 0;

Tac(1) = 0;

load('Theta_Tac0.mat');

load('Torque_Tac0.mat');

Torque_real = Torque_Tac0;

theta_real = Theta_Tac0;

dtheta_real = [0; diff(theta_real) ./ Ts];

du_initials = 0.5 * rand(CH, 1);  % initial the du initials

for k = 1:(NumberOfData - 1)

    %% compute LM for time k
    % display the process of the program
    if mod(k, 100) == 0
        disp([num2str((k / (NumberOfData - 1)) * 100), '%']);
    end

    if k < 17000
        Threshold(k) = 90;
    elseif k < 32000 && k > 17000
        Threshold(k) = 120;
    else
        Threshold(k) = 150;
    end

    if k == 1
        k_t(k) = 10 * FoM / LsT;
    else
        strain(k) = (LT(k - 1) - LsT) / LsT;
        if strain(k) > 0.0199
            strain(k) = 0.0199;
            stiffness(k) = polyval(kt_spline, strain(k));
        else
            stiffness(k) = polyval(kt_spline, strain(k));
        end
        k_t(k) = stiffness(k) * FoM / LsT;
    end

    LoM(k) = LoM_nominal * (1 + PercentageChange * (1 - activation_sew(k)));

    V_max(k) = 10 * LoM(k);

    LT(k) = LsT + MTF(k) / k_t(k);                                                % tendon length
    alpha(k) = acos(sqrt((LMT(k) - LT(k))^2 / (w^2 + (LMT(k) - LT(k))^2)));         % compute the pennation angles for the time k
    LM(k) = (LMT(k) - LT(k)) / cos(alpha(k));                                         % muscle length

    %% compute dLM

    F_L(k) = exp(-17.33 * (abs((LM(k) / LoM(k)) - 1)^3)); % equ (42)
    F_PE(k) = FoM * 3 * 10000 * exp(6 * ((LM(k) / LoM(k)) - 3.2)); % equ (12)
    % F_V(k)=(MTF(k) / cos(alpha(k))-F_PE(k)) / (activation(k) * FoM * exp(-17.33 * abs((LM(k) / LoM(k)) - 1)^3));%solve directly

    %% 1. directive method inverse solving

    %% 2. regression method

    e = [2.2 1];                                                             %  Fitted exactly to the linear part of Fv (k2 in the "Testing_Fv" file
    Z1 = MTF(k) / (FoM * cos(alpha(k))) - F_PE(k) / FoM - activation_sew(k) * F_L(k) * e(2);
    Z2 = activation_sew(k) * F_L(k) * e(1) / (V_max(k)) + damping / V_max(k);
    % dLM(k)=Vw_max * LoM * (0.5+exp(-3*(F_V(k)-0.02))+0.1*exp(12*(F_V(k)-1.5)));

    dLM(k) = Z1 / Z2;

    %% 3. linearized method
    % if k > 1
    % dLM(k)=(LM(k)-LM(k-1))/Ts;
    % end

    dalpha(k) = -dLM(k) / LM(k) * tan(alpha(k));

    % ddtheta(k) = 0;
    % if MTF(k)*MA(k) > 50
    %     ddtheta(k) = (MTF(k)*MA(k) + 50 - M) / I;
    % else
    %     ddtheta(k) = (MTF(k)*MA(k) + 0 - M) / I;
    % end

    % if theta(k) > 1.7
    %     % activation(k+1) = 0;
    %     k_lock_up = 2000;
    % else
    %     k_lock_up = 0;
    % end
    %
    %
    % F_k_up(k) =  k_lock_up*(abs(theta(k) - 1.7)) + 0.2*k_lock_up*(abs(dtheta(k) - 0));
    % F_k(k) =  0;

    % %%%%%%%%%%%%%%%%%%%PID using to track the theta = 1.5%%%%%%%%%%%%%%%%%%%%%%%

    % if (k < 7757) || (k > 13853 && k < 20956) || (k > 27639 && k < 36908) || (k > 41734) || (theta(k) > 1.5)
    % if k < 7757 || k > 13853
    % || (theta(k) > 1.5)

    if activation_sew(k) < 0.05
        err = theta_real(k) - theta(k);
        InT = InT + err * Ts;
        T_balance(k) = Kp * (err) + Kd * (dtheta_real(k) - dtheta(k)) + Ki * InT;
    else
        if (theta(k) > 1)  && (activation_sew(k) > 0.05)
            ee = 1 - theta(k);
            InT = InT + ee * Ts;
            F_k_up(k) = 0.2 * Kp * (ee) + 0.2 * Kd * (0 - dtheta(k)) + Ki * InT;
        elseif (theta(k) > 1.5)  && (activation_sew(k) > 0.05)
            ee = 1 - theta(k);
            InT = InT + ee * Ts;
            F_k_up(k) = 1 * Kp * (ee) + 1 * Kd * (0 - dtheta(k)) + Ki * InT;
        end
    end

    % %%%%%%%%%%%%%%%%%%%PID using to track the theta = 1.5%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%PID using to track the real theta when Tac(k) ~= 0%%%%%%%%%%%%%%%%%%%%%%%
    if Tac(k) ~= 0
        dactivation_sew = -P * (theta(k) - theta_real(k)) - D * (dtheta(k) - dtheta_real(k));
        if theta(k) < 1
            if dactivation_sew > 0.005
                dactivation_sew = 0.005;
            elseif dactivation_sew < -0.3
                dactivation_sew = -0.3;
            end
            activation_sew(k + 1) = activation_sew(k + 1) + dactivation_sew;
        else
            activation_sew(k + 1) = activation_sew(k + 1) + dactivation_sew;
        end

        % if (activation(k+1) - activation(k))/Ts > 50
        %     activation(k+1) = 50*Ts + activation(k);
        % end

        % elseif (activation(k+1) - activation(k))/Ts < -0.5
        %     activation(k+1) = -0.5*Ts + activation(k);
        if activation_sew(k + 1) < 0
            activation_sew(k + 1) = 0;
        elseif activation_sew(k + 1) > 1
            activation_sew(k + 1) = 1;
        end
    end

    %% the activation filter to smooth the data
    activation_sew(k + 1) = a * activation_sew(k) + (1 - a) * activation_sew(k + 1);
    %% the activation filter to smooth the data

    % %%%%%%%%%%%%%%%%%%PID using to track the real theta when Tac(k) ~= 0%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%data combination%%%%%%%%%%%%%%%%%%%%%%%
    if k > 9050 && k < 12537
        m = 5;
        I = m * L2^2 + (4 * M * L1^2) / 3;
    elseif k > 22662 && k < 26159
        m = 10;
        I = m * L2^2 + (4 * M * L1^2) / 3;
    elseif k > 37966 && k < 40544
        m = 15;
        I = m * L2^2 + (4 * M * L1^2) / 3;
    else
        m = 0;
        I = m * L2^2 + (4 * M * L1^2) / 3;
    end
    % %%%%%%%%%%%%%%%%%%data combination%%%%%%%%%%%%%%%%%%%%%%%

    % % %%%%%%%%%%%%%%%%%%Front rect muscle torque%%%%%%%%%%%%%%%%%%%%%%%

    T_front(k) = front_a * (theta(k) - 0.5)^2 + front_b * (theta(k) - 0.5) + front_c;

    % % %%%%%%%%%%%%%%%%%%Front rect muscle torque%%%%%%%%%%%%%%%%%%%%%%%

    ddtheta(k) = (-c * dtheta(k) + F_k_up(k)  + MTF(k) * MA(k) + Tac(k) + T_balance(k) - T_front(k) - m * g * L2 * cos(theta(k) - phi) - M * g * L1 * cos(theta(k))) / I;

    dtheta(k + 1) = dtheta(k) + ddtheta(k) * Ts;
    theta(k + 1) = theta(k) + dtheta(k) * Ts;

    %% integration and compute Ft

    LMT(k + 1) = LMT(k) + dLMT(k) * Ts;

    dLMT(k + 1) = -R * dtheta(k + 1);

    dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k))); % handnotes equ (37) and (16)

    MTF(k + 1) = MTF(k) + dFt(k) * Ts;

    if MTF(k + 1) < 0
        MTF(k + 1) = 0;
    end

    if isnan(MTF(k + 1))
        MTF(k + 1) = MTF(k);
    end

    if MTF(k + 1) > 50000
        MTF(k + 1) = MTF(k);
    end

    t(k + 1) = k * Ts;

    torque = MTF(k + 1) * MA(k + 1);
    %% optimization by cost function
    %%%%%%%%%%%%%%%%%%%%%%%%% NMPC%%%%%%%%%%%%%%

    if torque > 5

        Q = 12; % state weight
        RR = 3; % control weight
        W = 3; % increment of control weight
        InitialConditions = [MTF(k + 1), theta(k + 1), dtheta(k + 1), LMT(k + 1), dLMT(k + 1), LT(k)]; % initial state

        Predict_torque = nonlinear_system_period(m, zeros(PH, 1), InitialConditions, activation_sew(k:k + PH), theta_real(k:k + PH), dtheta_real(k:k + PH), damping, alpha0, Ts, coeffs);
        if max(Predict_torque(:, 1)) > Threshold(k)
            %%%%%%%%%%%%% activate the MPC%%%%%%%%%%%%%%

            for i = 1:PH
                if Predict_torque(i, 1) > Threshold(k)
                    Predict_torque(i, 1) = Threshold(k);
                end
            end

            r = Predict_torque(:, 1); % ref

            u0 = U(k - 1);

            cost_function = @(du) mpc_cost_function_increment(m, u0, du, n, InitialConditions, PH, CH, Q, W, RR, r, damping, alpha0, Ts, activation_sew(k:k + PH), theta_real(k:k + PH), dtheta_real(k:k + PH), coeffs);
            du0 = du_initials;
            lb = 0 * ones(CH, 1); % lower boundary of the incremental u
            ub = 0.5 * ones(CH, 1); % upper boundary of the incremental u

            A_lin = ones(1, CH); % linear constraint A matrix
            b_lin = u_max(k) - u0; % linear constraint b

            % using fmincon
            if mod(k, 100) == 0
                disp(k);
                options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'OptimalityTolerance', 1e-6, 'UseParallel', true);
            else
                options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'OptimalityTolerance', 1e-6, 'UseParallel', true);
            end

            % using GA
            % options = optimoptions('ga', 'Display', 'off', 'UseParallel', true);

            % using fmincon
            du = fmincon(cost_function, du0, A_lin, b_lin, [], [], lb, ub, [], options); % @(z) state_constraints(z, InitialConditions, N, Ts, 0, activation, coeffs, k_t(k))
            % using GA
            % z = ga(cost_function, N, [], [], [], [], lb, ub, [], options);

            % extract the first control input
            % du_initials = du;
            % if du(1)>5
            %     disp(du(1));
            % end

            Tac(k + 1) = u0 + du(1);

            % % apply the control input
            % [MTF_next, theta_next, dtheta_next] = nonlinear_system(k_t(k), InitialConditions, activation(k), u_Tac, Ts, coeffs);
            %
            x = [MTF(k + 1) * MA(k + 1), theta(k + 1), dtheta(k + 1)];
            % saving the result
            X(:, k) = x;
            U(k) = Tac(k + 1);
            Y(k) = x(1); % output the state x1
        else
            Tac(k + 1) = 0;
        end
    else
        Tac(k + 1) = 0;
        u0 = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
model_torque = MTF(1:end - 1)' .* MA';

% %%%%%%%correlation_coefficient%%%%%%%%
% R1 = corrcoef(Torque_real_1, Torque_Tac0); %torque coefficient
% CC_torque = R1(1, 2);
%
% R2 = corrcoef(theta_real_1, Theta_Tac0); %theta coefficient
% CC_theta = R2(1, 2);
%
% disp(CC_torque);
% disp(CC_theta);
% %%%%%%%correlation_coefficient%%%%%%%%

%%%%% plot the assistance/No assistance torque and the reference%%%%%
figure;
valid_idx = (k > 7757 & k < 13853) | (k > 20956 & k < 27639) | (k > 36908 & k < 41734);
subplot(2, 1, 1);
plot(t, Torque_Ceinms, 'k--', 'LineWidth', 2, 'DisplayName', 'Torque^{L5-S1}_{Ceinms}');
hold on;
plot(t, Torque_real, 'r', 'LineWidth', 3, 'DisplayName', 'Torque^{L5-S1}_{No assistance}');
text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('Torque (N.M)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show', 'Location', 'northwest', 'Orientation', 'horizontal'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
hold off;

subplot(2, 1, 2);
% plot(t, Threshold * ones(length(t), 1), 'k--', 'LineWidth', 2, 'DisplayName', 'F^{T}_{Thresh}');
plot(t, Threshold, '--k', 'LineWidth', 2, 'DisplayName', 'Torque^{L5-S1}_{Thresh}');
hold on;
plot(t, Torque_real, 'r--', 'LineWidth', 2, 'DisplayName', 'Torque^{L5-S1}_{No assistance}');
plot(t, model_torque, 'b', 'LineWidth', 3, 'DisplayName', 'Torque^{L5-S1}_{Assisted}');
plot(t, U(1:end - 1), 'g', 'LineWidth', 3, 'DisplayName', 'Torque^{Actuator}');
% plot(t, u_max * ones(length(t), 1), 'g--', 'LineWidth', 2, 'DisplayName', 'F^{Act}_{Thresh}');
plot(t, u_max, '--g', 'LineWidth', 2, 'DisplayName', 'Torque^{Actuator}_{Thresh}');
text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('Torque (N.M)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show', 'Location', 'northwest', 'Orientation', 'horizontal', 'NumColumns', 3); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
text(-0.5, Threshold(2), num2str(Threshold(2)), 'FontSize', 20, 'VerticalAlignment', 'middle');
text(-0.5, u_max(2), num2str(u_max(2)), 'FontSize', 20, 'VerticalAlignment', 'middle');

hold off;

%%%%% plot the bending angle theta and the reference%%%%%
figure;
subplot(2, 1, 1);
plot(t, theta_Ceinms, 'k--', 'LineWidth', 2, 'DisplayName', '\theta_{Ceinms}');
hold on;
plot(t, theta_real, 'r-', 'LineWidth', 3, 'DisplayName', '\theta_{No assistance}');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('\theta (rad)', 'FontSize', 24);
text(10000, 2, 'Stoop 5kg speed 45');
text(23000, 2, 'Stoop 10kg speed 35');
text(39000, 2, 'Stoop 15kg speed 45');
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 24;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;

subplot(2, 1, 2);
plot(t, theta_real, 'r--', 'LineWidth', 2, 'DisplayName', '\theta_{No assistance}');
hold on;
plot(t, theta, 'b-', 'LineWidth', 3, 'DisplayName', '\theta_{Assisted}');
text(10000, 2, 'Stoop 5kg speed 45');
text(23000, 2, 'Stoop 10kg speed 35');
text(39000, 2, 'Stoop 15kg speed 45');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('\theta (rad)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 24;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;

hold off;

%%%%% plot the torque and bending angle theta compare to the reference%%%%%
figure;
subplot(2, 1, 1);
plot(t, Torque_Ceinms, 'r--', 'LineWidth', 2, 'DisplayName', 'Torque^{L5-S1}_{Ceinms}');
hold on;
plot(t, Torque_real, 'b', 'LineWidth', 3, 'DisplayName', 'Torque^{L5-S1}_{No assistance}');
text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('Torque (N.M)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
hold off;

subplot(2, 1, 2);
plot(t, theta_Ceinms, 'r--', 'LineWidth', 2, 'DisplayName', '\theta_{Ceinms}');
hold on;
plot(t, Theta_Tac0, 'b-', 'LineWidth', 3, 'DisplayName', '\theta_{No assistance}');
text(10000, 2, 'Stoop 5kg speed 45');
text(23000, 2, 'Stoop 10kg speed 35');
text(39000, 2, 'Stoop 15kg speed 45');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('\theta (rad)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;

hold off;

%%%%% plot the assistance/No assistance torque and the reference, and theta%%%%%
figure;

subplot(2, 1, 1);
% plot(t, Threshold * ones(length(t), 1), 'k--', 'LineWidth', 2, 'DisplayName', 'F^{T}_{Thresh}');
plot(t, Threshold, '--k', 'LineWidth', 2, 'DisplayName', 'Torque^{L5-S1}_{Threshold}');
hold on;
plot(t, Torque_real, 'r--', 'LineWidth', 2, 'DisplayName', 'Torque^{L5-S1}_{No assistance}');
plot(t, model_torque, 'b', 'LineWidth', 3, 'DisplayName', 'Torque^{L5-S1}_{Assisted}');
plot(t, U(1:end - 1), 'g', 'LineWidth', 3, 'DisplayName', 'Torque^{Actuator}');
% plot(t, u_max * ones(length(t), 1), 'g--', 'LineWidth', 2, 'DisplayName', 'F^{Act}_{Thresh}');
plot(t, u_max, '--g', 'LineWidth', 2, 'DisplayName', 'Torque^{Actuator}_{Threshold}');
text(10000, 300, 'Stoop 5kg speed 45');
text(23000, 300, 'Stoop 10kg speed 35');
text(39000, 300, 'Stoop 15kg speed 45');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('Torque (N.M)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show', 'Location', 'northwest', 'Orientation', 'horizontal', 'NumColumns', 3); % display the legends
leg.FontSize = 14;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;
text(-0.5, Threshold(2), num2str(Threshold(2)), 'FontSize', 20, 'VerticalAlignment', 'middle');
text(-0.5, u_max(2), num2str(u_max(2)), 'FontSize', 20, 'VerticalAlignment', 'middle');

subplot(2, 1, 2);
plot(t, theta_real, 'r--', 'LineWidth', 2, 'DisplayName', '\theta_{No assistance}');
hold on;
plot(t, theta, 'b-', 'LineWidth', 3, 'DisplayName', '\theta_{Assisted}');
text(10000, 2, 'Stoop 5kg speed 45');
text(23000, 2, 'Stoop 10kg speed 35');
text(39000, 2, 'Stoop 15kg speed 45');
xlabel('Time (sec)', 'FontSize', 24);
ylabel('\theta (rad)', 'FontSize', 24);
% title('Sine and Cosine Waves');
leg = legend('show'); % display the legends
leg.FontSize = 24;
legend boxoff;  % remove the box of legends
ax = gca;
ax.FontSize = 24;

hold off;
