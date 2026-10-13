function x = Torque_function_constant(m, Torque_real_1, theta_real_1, Ts, M, L1, L2, phi)
    % TORQUE_FUNCTION_CONSTANT Integrate the lumped muscle and exoskeleton dynamics.
    % Use prepare_workflow to select matching code and input data.

    g = 9.81;

    NumberOfData = length(theta_real_1);

    %% Preallocations

    LoM = zeros(NumberOfData, 1)';
    alpha = zeros(NumberOfData, 1)';
    dalpha = zeros(NumberOfData, 1)';
    V_max = zeros(NumberOfData, 1)';
    LMT = zeros(NumberOfData, 1)';
    MA = zeros(NumberOfData, 1)';
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
    F_k_down = zeros(NumberOfData, 1)';
    T_balance = zeros(NumberOfData, 1)';
    Torque = zeros(NumberOfData, 1)';

    MTF = zeros(NumberOfData + 1, 1)';
    dFt = zeros(NumberOfData, 1)';
    Torque_gravity = zeros(NumberOfData, 1)';
    dTorque = zeros(NumberOfData, 1)';

    % Tac = zeros(NumberOfData,1)';

    %% compute the pennation angle for the first time

    I = m * L2^2 + (4 * M * L1^2) / 3;

    c = 10;

    P = 0.001;

    D = 0.0001;

    dtheta = [0; diff(theta_real_1) ./ Ts];
    ddtheta = [0; diff(dtheta) ./ Ts];
    dTorque_real_1 = [0; diff(Torque_real_1) ./ Ts];

    for k = 1:(NumberOfData - 1)

        %% compute LM for time k

        % if theta(k) < 0.6
        %     k_lock_down = 1000;
        % else
        %     k_lock_down = 0;
        % end
        %
        % if theta(k) > 1.8
        %     % activation(k+1) = 0;
        %     k_lock_up = 1000;
        % else
        %     k_lock_up = 0;
        % end
        %
        % F_k_down(k) =  k_lock_down*(abs(theta(k)-0.6));
        % F_k_up(k) =  k_lock_up*(abs(theta(k)-1.8));
        % F_k(k) =  0;

        % %%%%%%%%%%%%%%%%%%%PD using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%
        % if k < 8260 || k > 12862
        % T_balance(k) = - P*(theta(k)-theta_real(k)) - D*(dtheta(k) - dtheta_real(k));
        % end
        % %%%%%%%%%%%%%%%%%%%PD using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%

        % %%%%%%%%%%%%%%%%%%PD using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%
        %
        % if k == 1
        %     dTorque(k) = 0;
        % else
        %     dTorque(k) = (Torque(k) - Torque(k-1))/Ts;
        % end
        % theta_real_1(k+1) = theta_real_1(k+1) - P*(Torque(k) - Torque_real_1(k)) - D*(dTorque(k) - dTorque_real_1(k));
        %
        %
        % %%%%%%%%%%%%%%%%%%PD using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%

        Torque(k) = I * ddtheta(k) + m * g * L2 * cos(theta_real_1(k) - phi) + M * g * L1 * cos(theta_real_1(k)) + c * dtheta(k);
        Torque_gravity(k) = m * g * L2 * cos(theta_real_1(k) - phi) + M * g * L1 * cos(theta_real_1(k));

        %% integration and compute Ft

        t(k + 1) = k * Ts;

    end

    x(:, 1) = Torque;
    x(:, 2) = theta_real_1;
    x(:, 3) = dtheta;
    x(:, 4) = ddtheta;
    x(:, 5) = Torque_gravity;

    % %% plot figure
    %
    % MTF = MTF';
    %
    % figure;
    % plot(t(1:750),MTF(1:750));
    % ylim([-40, 90]);
    % hold on
    % plot(t(1:750), MTF_real(1:750));
    % plot(t(1:750), activation(1:750));
    % xlabel('t');
    % ylabel('Force');
    % ylim([-200, 1000]);
    % legend('MTF', 'F_desired','activation');
    % title('Ft(t) , F desired force, and 1000*activation');
    % grid on;
