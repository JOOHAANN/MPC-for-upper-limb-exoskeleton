function plant = bse_model(model, balanceController, standingAngle, useNative)
%BSE_MODEL Two compliant Hill muscles coupled to one trunk joint.
%   plant = bse_model(model, balanceController, standingAngle)
%   [state, output] = plant.initialize(activation, mass, waiting, initialState)
%   [state, output] = plant.step(state, activation, mass, waiting, u, dt, substeps)
%   output = plant.simulate(activation, mass, waiting, initialState, dt, substeps)
%
% Output columns are listed in plant.columns. Both torque and theta come from
% the same fiber/tendon/joint equilibrium. Reference trajectories are absent.
% Set useNative=false to inspect or execute the equivalent MATLAB equations.
if nargin < 4
    useNative = exist('bse_hill_mex_v2', 'file') == 3;
end
%% Fixed coefficients: physical model, muscle recruitment and neural control
p = model.parameters(:)';
r = model.neuralParameters(:)';
assert(any(numel(r)==[10 15 16]),'BSE:Parameters','Expected 10, 15 or 16 neural parameters.');
r = [r zeros(1,16-numel(r))];
% Neural order: intent amplitude, exponent, time constant, posture Kp/Kd,
% negative-request limit, feedback blend, coactivation, EMG lead, tonic gain.
% Optional terms: front, front*load, back*load, back^2, derivative*load, intent bias.
cfg = [standingAngle, balanceController(:)'];
assert(numel(p) == 32 && numel(r) == 16 && numel(cfg) == 7, ...
    'BSE:Parameters', 'Expected 32 physical, up to 16 neural and 7 standing parameters.');
assert(all(isfinite([p r cfg])), 'BSE:Parameters', 'Parameters must be finite.');
assert(strcmp(model.intentDecoder, 'subtractive'), ...
    'BSE:Decoder', 'This demo uses the subtractive intent decoder.');
physics = physical_parameters(p);
plant.initialize = @initialize;
plant.step = @advance;
plant.preview = @preview;
plant.simulate = @simulate;
plant.physics = physics;
plant.useNative = useNative;
plant.columns = {'theta','velocity','netTorque','backTorque','abdomenTorque', ...
    'backActive','backPassive','backFiberDamping','abdomenActive', ...
    'abdomenPassive','abdomenFiberDamping','backFiberLength','abdomenFiberLength', ...
    'backActivation','abdomenActivation','fiberResidual','momentumResidual', ...
    'backCommand','abdomenCommand','standingBlend','requestedHillTorque', ...
    'postureRequest','postureIntegral','motionIntent','assistanceTorque', ...
    'boundaryTorque','gravityTorque'};

    function [state, output] = initialize(a, mass, waiting, initial)
        if useNative
            [state, output] = bse_hill_mex_v2('init', p, r, cfg, ...
                double(a), double(mass), double(waiting), double(initial));
        else
            [state, output] = initial_state(physics, r, cfg, double(a), ...
                double(mass), double(waiting), double(initial));
        end
    end

    function [state, output] = advance(state, a, mass, waiting, u, dt, steps)
        if useNative
            [state, output] = bse_hill_mex_v2('step', p, r, cfg, state, ...
                double(a), double(mass), double(waiting), double(u), dt, double(steps));
        else
            [state, output] = advance_state(physics, r, cfg, state, double(a), ...
                double(mass), double(waiting), double(u), dt, steps);
        end
    end

    function output = simulate(a, mass, waiting, initial, dt, steps)
        if useNative
            output = bse_hill_mex_v2('rollout', p, r, cfg, a, ...
                double(mass), double(waiting), initial, dt, double(steps));
            return
        end
        output = zeros(size(a, 1), 27);
        [state, output(1,:)] = initialize(a(1,:), mass(1), waiting(1), initial);
        for k = 2:size(a, 1)
            [state, output(k,:)] = advance(state, a(k,:), mass(k), waiting(k), 0, dt, steps);
        end
    end

    function output = preview(state, controls, forecast, options)
        if useNative
            output = bse_hill_mex_v2('preview', p, r, cfg, state, controls, ...
                forecast.activation, forecast.mass, forecast.waiting, ...
                options.predictionDt, double(options.predictionSubsteps), options.rateNmPerSec);
            return
        end
        count = size(forecast.activation, 1);
        output = zeros(count, 27);
        try
            for k = 1:count
                block = floor((k-1)*numel(controls)/count) + 1;
                limit = options.rateNmPerSec*options.predictionDt;
                u = max(state(29)-limit, min(state(29)+limit, controls(block)));
                [state, output(k,:)] = advance(state, forecast.activation(k,:), ...
                    forecast.mass(k), forecast.waiting(k), u, ...
                    options.predictionDt, options.predictionSubsteps);
            end
        catch exception
            if strcmp(exception.identifier, 'BSE:Equilibrium')
                output = [];
            else
                rethrow(exception)
            end
        end
    end
end

function h = physical_parameters(p)
% Parameters keep their original order in the MAT files. Use names below.
h.bodyMass = p(19);
h.bodyLever = p(20);
h.loadLever = p(13);
h.loadOffset = p(14);
h.jointDamping = p(4);
h.pennationSine = sin(p(28));
h.boundaryStiffness = p(29);
h.boundaryCapacity = p(30);
h.boundaryLower = p(31);
h.boundaryUpper = p(32);
h.back = muscle(p(1), p(2), p(3), p(5), p(6), p(7), p(8), p(21), p(23), 1);
h.abdomen = muscle(p(9), p(11), p(12), p(25), p(10), p(26), p(27), p(22), p(24), -1);
end

function unit = muscle(capacity, ratio, referenceLength, vmax, exponent, strain, tau, arm, lengthMT, direction)
unit = struct('capacity',capacity, 'ratio',ratio, 'referenceLength',referenceLength, ...
    'maximumVelocity',vmax, 'activationExponent',exponent, 'tendonStrain',strain, ...
    'activationTime',tau, 'momentArm',arm, 'referenceMTLength',lengthMT, ...
    'optimalLength',arm/ratio, 'direction',direction);
end

function [s, out] = initial_state(h, r, cfg, a, mass, waiting, initial)
%% Initial joint angle, activation, compliant fibers and intent memory
% State layout follows the native kernel, with all 34 states retained.
% 1:2 joint; 3:4 activation; 5:6 fiber lengths; 7:8 stance/integral;
% 9:11 intent memories; 12:15 requests/commands; 16:27 fiber results;
% 28 momentum residual; 29 assist; 30 neutral input; 31:34 prior inputs.
s = zeros(34, 1);
s(1:2) = initial(:);
s(7) = clip(waiting, 0, 1);
s(9:10) = a(1);
s(30) = excitation(h.back, h, cfg(1), max(0, -gravity(h, cfg(1), 0)));
tonic = load_tonic(h, r, cfg, mass, s(30));
signal = intent_input(a(1),a(2),mass,0,r)-r(10)*tonic;
s(11) = cfg(1) + r(16) - r(1)*clip(signal, 0, 1)^r(2);
[command, s(12), s(13)] = human_request(h, r, cfg, s(1), s(2), ...
    mass, s(11), a(2), 0, 0, s(7) > 1e-8);
units = [h.back, h.abdomen];
for u = 1:2
    unit = units(u);
    feedback = s(7) + (1-s(7))*r(7);
    s(13+u) = (1-feedback)*clip(a(u), 0, 1) + feedback*command(u);
    s(2+u) = clip(s(13+u), 0, 1)^unit.activationExponent;
    f = fiber(unit, h, s(1), s(2+u), 0, .001, true);
    s = store_fiber(s, u, f, h);
end
s(31:34) = [a(:); mass; waiting];
out = state_output(h, s);
end

function [s, out] = advance_state(h, r, cfg, s, a, massEnd, waitEnd, assistEnd, dt, steps)
%% Coupled Hill fiber/tendon equilibrium and joint momentum balance
units = [h.back, h.abdomen];
step = dt/steps;
assistStart = s(29);
for k = 1:steps
    fraction = k/steps;
    mass = s(33) + fraction*(massEnd-s(33));
    assist = assistStart + fraction*(assistEnd-assistStart);
    targetBlend = clip(s(34) + fraction*(waitEnd-s(34)), 0, 1);
    s(7) = targetBlend + (s(7)-targetBlend)*exp(-step/cfg(7));
    s(8) = clip(s(8)+step*s(7)*(cfg(1)-s(1)), -.5, .5);
    input = clip(s(31:32)' + fraction*(a(:)'-s(31:32)'), 0, 1);
    s(9) = input(1) + (s(9)-input(1))*exp(-step/.06);
    s(10) = input(1) + (s(10)-input(1))*exp(-step/.25);
    tonic = load_tonic(h, r, cfg, mass, s(30));
    signal = clip(intent_input(input(1),input(2),mass,s(9)-s(10),r) ...
        +r(9)*(s(9)-s(10))/(.25-.06)-r(10)*tonic, 0, 1);
    targetIntent = cfg(1)+r(16)-r(1)*signal^r(2);
    s(11) = targetIntent + (s(11)-targetIntent)*exp(-step/r(3));

    post = [0 0]; moving = [0 0]; dp = 0; cp = 0; dm = 0; cm = 0;
    if s(7) > 1e-8
        [post, dp, cp] = human_request(h, r, cfg, s(1), s(2), mass, ...
            s(11), input(2), s(8), assist, true);
    end
    if s(7) < 1-1e-8
        [moving, dm, cm] = human_request(h, r, cfg, s(1), s(2), mass, ...
            s(11), input(2), s(8), assist, false);
    end
    s(12) = s(7)*dp + (1-s(7))*dm;
    s(13) = s(7)*cp + (1-s(7))*cm;
    for u = 1:2
        s(13+u) = s(7)*post(u) + (1-s(7))*((1-r(7))*input(u)+r(7)*moving(u));
        target = clip(s(13+u), 0, 1)^units(u).activationExponent;
        s(2+u) = target+(s(2+u)-target)*exp(-step/units(u).activationTime);
    end

    % Backward Euler solves muscle equilibrium and momentum together.
    inertia = mass*h.loadLever^2 + 4*h.bodyMass*h.bodyLever^2/3;
    qold = s(1); vold = s(2); guess = qold+step*vold;
    for iteration = 1:9
        back = fiber(h.back, h, guess, s(3), s(5), step, false);
        abdomen = fiber(h.abdomen, h, guess, s(4), s(6), step, false);
        torque = h.back.capacity*back(2)-h.abdomen.capacity*abdomen(2);
        residual = inertia*(guess-qold-step*vold)/step^2 ...
            + h.jointDamping*(guess-qold)/step-torque ...
            - gravity(h, guess, mass)-boundary(h, guess)-assist;
        if abs(residual) < 1e-7
            break
        end
        epsilon = 1e-6;
        bp = fiber(h.back, h, guess+epsilon, s(3), s(5), step, false);
        ap = fiber(h.abdomen, h, guess+epsilon, s(4), s(6), step, false);
        torquePlus = h.back.capacity*bp(2)-h.abdomen.capacity*ap(2);
        derivative = inertia/step^2+h.jointDamping/step-(torquePlus-torque)/epsilon ...
            -(gravity(h, guess+epsilon, mass)-gravity(h, guess, mass))/epsilon ...
            -(boundary(h, guess+epsilon)-boundary(h, guess))/epsilon;
        guess = guess-clip(residual/derivative, -.1, .1);
    end
    s(1) = guess;
    s(2) = (guess-qold)/step;
    back = fiber(h.back, h, guess, s(3), s(5), step, false);
    abdomen = fiber(h.abdomen, h, guess, s(4), s(6), step, false);
    s = store_fiber(s, 1, back, h);
    s = store_fiber(s, 2, abdomen, h);
    s(28) = inertia*(s(2)-vold)/step-h.back.capacity*s(16) ...
        +h.abdomen.capacity*s(17)-gravity(h, guess, mass)-boundary(h, guess) ...
        +h.jointDamping*s(2)-assist;
    if ~all(isfinite(s)) || abs(s(1)) > 4 || abs(s(2)) > 50 ...
            || sum(abs(s(26:27))) > 1e-4 || abs(s(28)) > 1
        error('BSE:Equilibrium', 'The coupled Hill solve did not converge.');
    end
end
s(29) = assistEnd;
s(31:34) = [a(:); massEnd; waitEnd];
out = state_output(h, s);
end

function f = fiber(unit, h, q, activation, oldLength, dt, initial)
% Tensile tendon: FT/F0 = positive(strain / referenceStrain)^2.
% Fiber equilibrium: FT = (a*fL*fV + fPE + .03*vN)*cos(pennation).
w = h.pennationSine;
referenceProjection = sqrt(unit.referenceLength^2-w^2);
slack = unit.referenceMTLength-unit.optimalLength*referenceProjection;
compliance = unit.tendonStrain*slack/unit.optimalLength;
lengthMT = referenceProjection-unit.direction*unit.ratio*(q-1.4);
lower = 0;
upper = max(1e-5, lengthMT+slack/unit.optimalLength);
if initial
    x = (lower+upper)/2;
else
    x = sqrt(max(0, oldLength^2-w^2));
end
x = clip(x, lower+1e-10, upper-1e-10);
for iteration = 1:40
    length = sqrt(x*x+w*w);
    cosine = x/length;
    dcosine = w*w/length^3;
    velocity = 0; dvelocity = 0;
    if ~initial
        velocity = (length-oldLength)/(dt*unit.maximumVelocity);
        dvelocity = cosine/(dt*unit.maximumVelocity);
    end
    stretch = max(0, lengthMT-x);
    tendon = stretch^2/compliance^2;
    dtendon = -2*stretch/compliance^2;
    fl = exp(-((length-1)/.45)^2);
    dfl = -2*(length-1)/.45^2*fl*cosine;
    fv = force_velocity(velocity);
    if velocity <= -1
        dfv = 0;
    elseif velocity <= 0
        dfv = 5/(1-4*velocity)^2;
    else
        dfv = 2/(1+4*velocity)^2;
    end
    fp = passive_force(length);
    dfp = 0;
    if length > 1 && 5*(length-1)/.6 < 35
        dfp = (5/.6)*exp(5*(length-1)/.6)*cosine/expm1(5);
    end
    active = activation*fl*fv;
    dactive = activation*(dfl*fv+fl*dfv*dvelocity);
    error = tendon-(active+fp+.03*velocity)*cosine;
    derivative = dtendon-(dactive+dfp+.03*dvelocity)*cosine ...
        -(active+fp+.03*velocity)*dcosine;
    if abs(error) < 3e-11
        break
    end
    if error > 0, lower = x; else, upper = x; end
    next = x-error/derivative;
    if ~isfinite(next) || next <= lower || next >= upper
        next = (lower+upper)/2;
    end
    x = next;
end
length = sqrt(x*x+w*w); cosine = x/length;
velocity = 0;
if ~initial, velocity = (length-oldLength)/(dt*unit.maximumVelocity); end
tendon = (max(0, lengthMT-x)/compliance)^2;
active = activation*exp(-((length-1)/.45)^2)*force_velocity(velocity)*cosine;
passive = passive_force(length)*cosine;
f = [x, tendon, active, passive, velocity, tendon-active-passive-.03*velocity*cosine];
end

function [command, desired, correction] = human_request(h, r, cfg, q, v, mass, intent, frontInput, integral, assist, standing)
% Postural feedback changes excitation; it is never an extra joint torque.
units = [h.back, h.abdomen];
minimum = zeros(1, 2);
for u = 1:2
    f = fiber(units(u), h, q, 0, 0, .001, true);
    minimum(u) = units(u).capacity*f(2);
end
if standing
    correction = cfg(5)*tanh((cfg(2)*(cfg(1)-q)-cfg(3)*v+cfg(4)*integral)/cfg(5));
    desired = -gravity(h, q, mass)+correction-assist;
    difference = desired-minimum(1)+minimum(2);
    targets = minimum+cfg(6)+[max(difference,0), max(-difference,0)];
else
    f = fiber(h.abdomen, h, q, clip(frontInput,0,1)^h.abdomen.activationExponent, 0, .001, true);
    front = h.abdomen.capacity*f(2);
    correction = r(4)*(intent-q)-r(5)*v;
    desired = -gravity(h, q, mass)+correction-assist;
    if desired < 0, desired = r(6)*tanh(desired/r(6)); end
    targets = [max(minimum(1),desired+front), max(front,minimum(1)-desired)]+r(8);
end
command = [excitation(h.back,h,q,targets(1)), excitation(h.abdomen,h,q,targets(2))];
end

function a = excitation(unit, h, q, torque)
% Static inversion only requests excitation; dynamic forces remain free.
w = h.pennationSine;
slack = unit.referenceMTLength-unit.optimalLength*sqrt(unit.referenceLength^2-w*w);
lengthMT = unit.referenceMTLength-unit.direction*unit.momentArm*(q-1.4);
force = max(0, torque/unit.capacity);
tendonLength = slack*(1+unit.tendonStrain*sqrt(force));
x = (lengthMT-tendonLength)/unit.optimalLength;
if x <= 0, a = 0; return, end
length = sqrt(x*x+w*w);
activation = (force/(x/length)-passive_force(length))/max(exp(-((length-1)/.45)^2),1e-8);
a = clip(activation,0,1)^(1/unit.activationExponent);
end

function tonic = load_tonic(h, r, cfg, mass, neutral)
tonic = 0;
if r(10) > 0
    tonic = max(0, excitation(h.back,h,cfg(1),max(0,-gravity(h,cfg(1),mass)))-neutral);
end
end

function value = intent_input(back,front,mass,difference,r)
% Fixed cross-muscle/load features. No reference angle or trial identity.
load = mass/15;
value = back+r(11)*front+r(12)*front*load+r(13)*back*load ...
    +r(14)*back^2+r(15)*difference*load;
end

function s = store_fiber(s, u, f, h)
s(4+u) = sqrt(f(1)^2+h.pennationSine^2);
s(15+u) = f(2); s(17+u) = f(3); s(19+u) = f(4);
s(21+u) = f(5); s(23+u) = f(1); s(25+u) = f(6);
end

function out = state_output(h, s)
back = h.back.capacity; abdomen = h.abdomen.capacity;
out = [s(1:2)', back*s(16)-abdomen*s(17), back*s(16), abdomen*s(17), ...
    back*s(18), back*s(20), back*.03*s(22)*s(24)/s(5), ...
    abdomen*s(19), abdomen*s(21), abdomen*.03*s(23)*s(25)/s(6), ...
    s(5:6)', s(3:4)', max(abs(s(26:27))), s(28), s(14:15)', s(7), ...
    s(12:13)', s(8), s(11), s(29), boundary(h,s(1)), gravity(h,s(1),s(33))];
end

function torque = gravity(h, q, mass)
torque = -9.81*h.bodyMass*h.bodyLever*cos(q) ...
    -9.81*mass*h.loadLever*cos(q-h.loadOffset);
end

function torque = boundary(h, q)
% Elastic soft boundary only. There is no hidden boundary damping.
low = max(h.boundaryLower-q,0); high = max(q-h.boundaryUpper,0);
low = low^2/(low+.05); high = high^2/(high+.05);
torque = h.boundaryCapacity*(tanh(h.boundaryStiffness*low/h.boundaryCapacity) ...
    -tanh(h.boundaryStiffness*high/h.boundaryCapacity));
end

function value = passive_force(length)
value = max(0, expm1(min(35,5*(length-1)/.6))/expm1(5));
end

function value = force_velocity(velocity)
if velocity <= 0
    value = max(0,1+velocity)/(1-velocity/.25);
else
    value = (1+1.5*velocity/.25)/(1+velocity/.25);
end
end

function value = clip(value, lower, upper)
value = max(lower, min(upper, value));
end
