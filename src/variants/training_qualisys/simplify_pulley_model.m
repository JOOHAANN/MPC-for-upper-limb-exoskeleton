%% this function can compute the LMT of pulley model

function [LMT, MA] = simplify_pulley_model(L1, R, theta)

    LMT = sqrt(L1^2 - R^2) + R * (pi / 2 - theta);

    % figure
    % plot(theta,LMT);
    % xlabel('theta');
    % ylabel('LMT');
    % title('simplify model');

    %% the length of moment arm

    MA = R;

    % figure
    % plot(theta,BE);
    % xlabel('theta');
    % ylabel('Moment arm');
    % title('Moment arm against theta');

    %% dynamics
end
