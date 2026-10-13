%% this function can compute the LMT of pulley model

function [LMT, BE] = simplify_pulley_model(L2, L3, R, theta)

    % theta = linspace(40, 90, 300);

    %% when LMT = 0.33
    % L2 = 0.33;
    % L3 = 0.003;
    % R = 0.03;

    %% the when LMT = 0.17
    % L2 = 0.17280;
    % L3 = 0.006;
    % R = 0.052;

    %% the big muscle IL
    % L2 = 0.102;
    % L3 = 0.005;
    % R = 0.0412;

    LMT = sqrt(L3^2 + L2^2 + 2 * L2 * L3 * cos(theta) - R^2);

    % figure
    % plot(theta,LMT);
    % xlabel('theta');
    % ylabel('LMT');
    % title('simplify model');

    %% the length of moment arm

    AC = L3^2 + L2^2 + 2 * L2 * L3 * cos(theta);

    C1 = (AC + L2^2 - L3^2) ./ (2 * L2 * sqrt(AC));

    Sa1a2 = sqrt(1 - C1.^2) .* LMT ./ sqrt(AC) + C1 .* R ./ sqrt(AC);

    BE = Sa1a2 * L2;

    % figure
    % plot(theta,BE);
    % xlabel('theta');
    % ylabel('Moment arm');
    % title('Moment arm against theta');

    %% dynamics
end
