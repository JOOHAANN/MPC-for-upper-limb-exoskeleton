function [NextActivation, STIM] = ActivationDynamics(StartingStep, DutyPercent, Period, StepsToStabilize, activation, i)

    persistent EndOfPeriod1 EndOfPeriod2 EndOfPeriod3 EndOfPeriod4 EndOfPeriod5 EndOfPeriod6 EndOfPeriod7 EndOfPeriod8 EndOfPeriod9 EndOfPeriod10 EndOfPeriod11 EndOfPeriod12 EndOfPeriod13 EndOfPeriod14 EndOfPeriod15 EndOfPeriod16 EndOfPeriod17 EndOfPeriod18 EndOfPeriod19 EndOfPeriod20

    global dt

    if isempty(i)
        EndOfPeriod1 = 0;
        EndOfPeriod2 = 0;
        EndOfPeriod3 = 0;
        EndOfPeriod4 = 0;
        EndOfPeriod5 = 0;
        EndOfPeriod6 = 0;
        EndOfPeriod7 = 0;
        EndOfPeriod8 = 0;
        EndOfPeriod9 = 0;
        EndOfPeriod10 = 0;
        EndOfPeriod11 = 0;
        EndOfPeriod12 = 0;
        EndOfPeriod13 = 0;
        EndOfPeriod14 = 0;
        EndOfPeriod15 = 0;
        EndOfPeriod16 = 0;
        EndOfPeriod17 = 0;
        EndOfPeriod18 = 0;
        EndOfPeriod19 = 0;
        EndOfPeriod20 = 0;
    end

    if floor(i / (StartingStep + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.05;
        deactivate = 0.0;
        EndOfPeriod1 = i + 1;

    elseif floor(i / (EndOfPeriod1 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.1;
        deactivate = 0.0;
        EndOfPeriod2 = i + 1;

    elseif floor(i / (EndOfPeriod2 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.15;
        deactivate = 0.0;
        EndOfPeriod3 = i + 1;

    elseif floor(i / (EndOfPeriod3 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.2;
        deactivate = 0.0;
        EndOfPeriod4 = i + 1;

    elseif floor(i / (EndOfPeriod4 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.25;
        deactivate = 0.0;
        EndOfPeriod5 = i + 1;

    elseif floor(i / (EndOfPeriod5 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.3;
        deactivate = 0.0;
        EndOfPeriod6 = i + 1;

    elseif floor(i / (EndOfPeriod6 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.35;
        deactivate = 0.0;
        EndOfPeriod7 = i + 1;

    elseif floor(i / (EndOfPeriod7 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.4;
        deactivate = 0.0;
        EndOfPeriod8 = i + 1;

    elseif floor(i / (EndOfPeriod8 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.45;
        deactivate = 0.0;
        EndOfPeriod9 = i;

    elseif floor(i / (EndOfPeriod9 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.5;
        deactivate = 0.0;
        EndOfPeriod10 = i;

    elseif floor(i / (EndOfPeriod10 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.55;
        deactivate = 0.0;
        EndOfPeriod11 = i;

    elseif floor(i / (EndOfPeriod11 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.6;
        deactivate = 0.0;
        EndOfPeriod12 = i;

    elseif floor(i / (EndOfPeriod12 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.65;
        deactivate = 0.0;
        EndOfPeriod13 = i;

    elseif floor(i / (EndOfPeriod13 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.7;
        deactivate = 0.0;
        EndOfPeriod14 = i;

    elseif floor(i / (EndOfPeriod14 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.75;
        deactivate = 0.0;
        EndOfPeriod15 = i;

    elseif floor(i / (EndOfPeriod15 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.8;
        deactivate = 0.0;
        EndOfPeriod16 = i;

    elseif floor(i / (EndOfPeriod16 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.85;
        deactivate = 0.0;
        EndOfPeriod17 = i;

    elseif floor(i / (EndOfPeriod17 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.9;
        deactivate = 0.0;
        EndOfPeriod18 = i;

    elseif floor(i / (EndOfPeriod18 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.95;
        deactivate = 0.0;
        EndOfPeriod19 = i;

    elseif floor(i / (EndOfPeriod19 + StepsToStabilize * Period)) == 0
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 1;
        deactivate = 0.0;
        EndOfPeriod20 = i;

    else
        PeriodSteps = Period;
        DutySteps = DutyPercent * PeriodSteps;
        activate = 0.0;
        deactivate = 0.0;
    end

    CurrentPlaceInPeriod = i - StartingStep - floor((i - StartingStep) / PeriodSteps) * PeriodSteps;

    if i <= StartingStep || CurrentPlaceInPeriod > DutySteps
        STIM = deactivate;
        tau(i) = 0.091;
    else
        STIM = activate;
        tau(i) = 0.033;
    end

    NextActivation = dt * (STIM - activation) / tau(i) + activation;

end
