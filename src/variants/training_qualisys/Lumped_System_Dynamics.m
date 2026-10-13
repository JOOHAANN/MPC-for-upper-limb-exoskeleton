% States:  x= [FT, LMT, dLMT, Fac];

% The difference of this code wrt the previous (NMPC for the paper) is that
% in here we assume that we know the change is going to be happening in the
% activation in the future (because we know the treadmill speed).

clc;
clear;
close all;

global LoM_nominal PercentageChange FoM LsT w dt g kt_spline damping  InitialHorizon InputIncrement

% NumberOfData= 26000;

%% Soleus parameters
PercentageChange = 0.15;
damping = 0.4;
% damping= 0.8;

% kt_spline is extracted from the file "Stiffness.m"
kt_spline = [6.07450471822638e+42    -1.17093401235721e+42   9.78861307802160e+40    -4.68161016668128e+39   1.41111335336212e+38    -2.75023728531764e+36   3.32535388381590e+34    -1.97635568280683e+32   -4.74852712911846e+29   1.68812705982292e+28    -9.30681300424510e+25   -1.69227796153188e+23   3.64048821058250e+21    -7.76789223923884e+18   -4.91004587706162e+16   201470178531589 239829087287.116    -1720265398.11745   -15739.3275509684   6708.90359529906    3.90612013684262];

%% Lumped model based on CEINMS. Soleus is dominant.
FoM = 3500;

LoM_nominal = 0.103;
LsT = 0.24044;
alpha0 = 20.14 * pi / 180;

w = LoM_nominal * sin(alpha0);
g = 9.81;

%% Activation dynamics parameters
StartingStep = 50;
DutyPercent = 0.1;
Period = 800;
StepsToStabilize = 4;

NumberOfData = StartingStep + 20 * Period * StepsToStabilize;

%% Preallocations
LoM = zeros(NumberOfData, 1)';
alpha = zeros(NumberOfData, 1)';
dalpha = zeros(NumberOfData, 1)';
Vmax = zeros(NumberOfData, 1)';
LMT = zeros(NumberOfData, 1)';
LM = zeros(NumberOfData, 1)';
LT = zeros(NumberOfData, 1)';
dLMT = zeros(NumberOfData, 1)';
kT = zeros(NumberOfData, 1)';
stiffness = zeros(NumberOfData, 1)';
FL = zeros(NumberOfData, 1)';
FPE = zeros(NumberOfData, 1)';
ddLMT = zeros(NumberOfData, 1)';
t = zeros(NumberOfData, 1)';
ddtheta = zeros(NumberOfData, 1)';
dtheta = zeros(NumberOfData, 1)';
theta = zeros(NumberOfData, 1)';
MA = zeros(NumberOfData, 1)';

FT = zeros(NumberOfData + 1, 1)';
dFT = zeros(NumberOfData, 1)';
Fac = zeros(NumberOfData, 1)';
strain = zeros(NumberOfData, 1)';

STIM = zeros(NumberOfData, 1)';
Factuator = zeros(NumberOfData, 1)';
W = zeros(NumberOfData, 1)';
activation = zeros(NumberOfData + 1, 1)';
ActPredict = zeros(InitialHorizon, NumberOfData + 1)';
%% 0.352956
InitialConditions = [FT(1)    LMT(1)   dLMT(1)    LM(1) Fac(1)    activation(1)     LT(1)];          % Initial conditions for assigning initial conditions for the muscle and actuator
InitialConditions = [150    0.352956   -0.0000    0.0621  694.8291    0.0100    0.2380];        % Initial conditions for assigning initial conditions for the muscle and actuator
PeriodCounter = ones(20 * StepsToStabilize + 1, 1)';
FTmax = ones(20 * StepsToStabilize + 1, 1)';
FutureStims = zeros(NumberOfData, InitialHorizon);

% ControlHorizon= InitialHorizon*ones(NumberOfData,1);
FutureDesired_FT = zeros(NumberOfData, InitialHorizon);
FutureFT = zeros(NumberOfData, InitialHorizon);

FT(1) = InitialConditions(1);
LMT(1) = InitialConditions(2);
dLMT(1) = InitialConditions(3);
LM(1) = InitialConditions(4);
Fac(1) = InitialConditions(5);
activation(1) = InitialConditions(6);
LT(1) = InitialConditions(7);
theta(1) = 0.3;
dtheta(1) = 0;
MA(1) = 0.032;

ii = 1;
dt = 0.001; % Sampling time
%%

load('activation.mat');

NumberOfData = length(Smooth_activation);

activation = Smooth_activation;

% load('FT0.mat');
for i = 1:NumberOfData

    t(i) = (i - 1) * dt;
    W(i) = 35 * g * 3;

    %% First order activation dynamics

    % if i==1
    %     activation(i)= 0.0;
    % end
    %
    % [NextActivation, STIM(i)]= ActivationDynamics(StartingStep, DutyPercent, Period, StepsToStabilize, activation(i), i);
    %  activation(i+1)= NextActivation;
    %  if i > 0
    % activation(i+1)= 0;
    %  end
    if i > 1
        if STIM(i) > STIM(i - 1)
            ii = ii + 1;
            PeriodCounter(ii) = i;
        end
    end

    Fac(i) = 0;

    %% Muscle and hopping dynamics
    LoM(i) = LoM_nominal * (1 + PercentageChange * (1 - activation(i)));         % LoM = Length of Muscle optimal
    Vmax(i) = 10 * LoM(i);

    if i == 1
        kT(i) = 10 * FoM / LsT;
    else
        strain(i) = (LT(i - 1) - LsT) / LsT;
        stiffness(i) = polyval(kt_spline, strain(i));
        kT(i) = stiffness(i) * FoM / LsT;
    end
    LT(i) = LsT + FT(i) / kT(i);                                               % Tendon length
    alpha(i) = acos(sqrt((LMT(i) - LT(i))^2 / (w^2 + (LMT(i) - LT(i))^2)));      % Pennation angle
    LM(i) = (LMT(i) - LT(i)) / cos(alpha(i));                                  % Muscle length

    FL(i) = exp(-17.33 * (abs(LM(i) / LoM(i) - 1))^3);                           % Active Force-length relation
    FPE(i) = FoM * 3e4 * exp(6 * (LM(i) / LoM(i) - 3.2));                            % Passive Force-length relation

    k = [2.2 1];                                                             %  Fitted exactly to the linear part of Fv (k2 in the "Testing_Fv" file
    Z1 = FT(i) / (FoM * cos(alpha(i))) - FPE(i) / FoM - activation(i) * FL(i) * k(2);
    Z2 = activation(i) * FL(i) * k(1) / (Vmax(i)) + damping / Vmax(i);

    dLM(i) = Z1 / Z2;
    % if i == 1
    % dLM(i) = 0;
    % else
    % dLM(i) = LM(i) - LM(i-1);
    % end

    dalpha(i) = -dLM(i) / LM(i) * tan(alpha(i));                                 % Pennation angle derivative

    dFT(i) = kT(i) * (dLMT(i) - dLM(i) * cos(alpha(i))  + LM(i) * sin(alpha(i)) * dalpha(i));
    ddLMT(i) = -g / W(i) * (FT(i) + Fac(i) - W(i));

    %% Integration
    FT(i + 1) =     FT(i) + dt * dFT(i);
    dLMT(i + 1) = dLMT(i) + dt * ddLMT(i);
    LMT(i + 1) =   LMT(i) + dt * dLMT(i);
    %
end

%% plot the MTF
figure;
plot(t, FT(1:64050));

% FT0 = FT;
% save('lumped.mat', 'FT0');
