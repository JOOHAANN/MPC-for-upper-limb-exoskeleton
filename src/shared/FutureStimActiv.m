function    [Stimulations, STIM, FutureActivations, NextActivation] = ...
    FutureStimActiv(ii, activation, Period, StartingStep, DutyPercent, StepsToStabilize)
    % This function's output is the muscle excitaion of the next time steps;
    % both cases: next one and the future (up to horizon). So the number of its
    % elements are the same as the horizon (1:InitialHorizon)

    % Excitation (which I might have called it stimulation in the codes) isn't
    % important for us. The important thing is the activation which we should
    % be sure to know the timing of it. Excitation is only used to spot changes
    % and counting in some part of the code.

    global   dt InitialHorizon

    ZeroExcit = 0.05;
    ZeroAct = 0.05;

    %% Preallocation
    Stimulations = zeros(1, InitialHorizon);
    tau = zeros(1, InitialHorizon);
    FutureActivations = zeros(1, InitialHorizon);

    %% Main code
    for jj = 1:InitialHorizon

        if ii + jj <= StartingStep

            Stimulations(jj) = ZeroExcit;
            FutureActivations(jj) = ZeroAct;

        elseif ii + jj >= StartingStep + 1

            CurrCycle = floor((ii + jj - StartingStep) / Period);

            if (ii + jj - (StartingStep + 1)) - CurrCycle * Period <= DutyPercent * Period

                Stimulations(jj) = (floor(CurrCycle / StepsToStabilize) + 1) * 0.05;
                tau(jj) = 0.015;

            else

                Stimulations(jj) = ZeroExcit;
                tau(jj) = 0.091;

            end

            if jj == 1

                FutureActivations(jj) = dt * (Stimulations(jj) - activation) / tau(jj) + activation;

            else

                FutureActivations(jj) = dt * (Stimulations(jj) - FutureActivations(jj - 1)) / tau(jj) + FutureActivations(jj - 1);

            end

        end

        NextActivation = FutureActivations(1);
        STIM = Stimulations(1);

    end
