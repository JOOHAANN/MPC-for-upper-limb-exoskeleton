%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% getXP.mat
%%%
%%% Estimates the exitacion primitive (XP) for a specified speed and
%%% elevation. There are 2 baseline modes: 'GM', which is the generic mode
%%% and no file should be indicated and 'SM', which is the subject-specific
%%% mode and a file should be indicated in the next argument (baseline_file)
%%%
%%% Jose Gonzalez (2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [weights] = getWeights(speed, elevation, baseline_mode, baseline_file)

    Alpha = 1;
    Beta = 1;

    % load the regression constants
    % load(strcat('regression_const'),strcat('regression_const'));

    % model calculated only with speeds: 1 3 5 and elevations -20 0 20
    % load(strcat('withExceptionsreg_2step_Add_all'),strcat('reg_2step_Add_all'));
    % regression_const = reg_2step_Add_all;

    % model calculated with all data. Should be the same as the original one.
    load(strcat('regression_const'));
    % regression_const = reg_2step_Add_all;

    %%%% get the subject-specific baseline
    if strcmp(baseline_mode, 'SM')
        load(baseline_file);
        load('pv_coefficients.mat');
        load('pw_coefficients.mat');
        baseline_val = baseline_data;

        %%%%% get the generic mode baseline
    elseif strcmp(baseline_mode, 'GM')
        load('weights_generic', 'weights_generic');
        baseline_val = weights_generic';
        clear weight_generic;
    end

    % solves the equations for each muscle
    weights = zeros(2, 6); %%%%%%%%%%%%%%%%%%%%
    for c = 1:2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for m = 1:6

            coef = pw_coefficients(:, c, m);
            temp1 = coef(1) * 10^2;
            temp1 = temp1 + coef(2) * 10;
            temp1 = temp1 + coef(3);

            % coef = pw_coefficients(:,c,m);
            % temp1 = coef(1)*5;
            % temp1 = temp1 + coef(2);

            coef = pv_coefficients(:, c, m);
            temp2 = coef(1) * 45^2;
            temp2 = temp2 + coef(2) * 45;
            temp2 = temp2 + coef(3);

            correction = Alpha * temp1 + Beta * temp2;

            coef = pw_coefficients(:, c, m);
            temp1 = coef(1) * elevation^2;
            temp1 = temp1 + coef(2) * elevation;
            temp1 = temp1 + coef(3);

            % coef = pw_coefficients(:,c,m);
            % temp1 = coef(1)*elevation;
            % temp1 = temp1 + coef(2);

            coef = pv_coefficients(:, c, m);
            temp2 = coef(1) * speed^2;
            temp2 = temp2 + coef(2) * speed;
            temp2 = temp2 + coef(3);

            % weights(c,m) = (baseline_val(c,m) + temp1) + temp2-correction;
            weights(c, m) = Alpha * temp1 + Beta * temp2 + baseline_val(c, m) - correction;

            if weights(c, m) < 0
                weights(c, m) = 0; % deals with negative values.
            end
        end

        if c == 1
            cc = 2;
        elseif c == 2
            cc = 1;
        else
            cc = c;
        end
    end

end
