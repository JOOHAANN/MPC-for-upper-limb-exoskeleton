function plant = dual_MTU_model(parameters,balance,standingAngle,useNative)
%DUAL_MTU_MODEL Direct activation, two compliant MTUs, standing balance only.
% The complete MATLAB equations are below. The optional MEX solver performs
% the same calculation faster. There is no intent decoder, movement PD,
% angle boundary torque or angle clipping.
if nargin<4,useNative=true;end
p=parameters(:)';assert(numel(p)==24);
cfg=[standingAngle balance(:)'];h=physical_parameters(p);
if useNative,build_solver;end
plant.initialize=@initialize;plant.step=@advance;plant.simulate=@simulate;
plant.parameters=p;plant.balance=cfg;
    function [state,out]=initialize(input,mass,waiting,initial)
        if useNative
            [state,out]=direct_mtu_mex('init',p,cfg,double(input),double(mass),double(waiting),double(initial));
        else
            [state,out]=initial_state(h,cfg,input,mass,waiting,initial);
        end
    end
    function [state,out]=advance(state,input,mass,waiting,dt,steps)
        if useNative
            [state,out]=direct_mtu_mex('step',p,cfg,state,double(input),double(mass),double(waiting),dt,double(steps));
        else
            [state,out]=advance_state(h,cfg,state,input,mass,waiting,dt,steps);
        end
    end
    function [out,completed]=simulate(inputs,mass,waiting,initial,dt,steps)
        if useNative
            [out,completed]=direct_mtu_mex('rollout',p,cfg,double(inputs),double(mass), ...
                double(waiting),double(initial),dt,double(steps));
        else
            n=size(inputs,1);out=nan(n,27);completed=1;
            [state,out(1,:)]=initialize(inputs(1,:),mass(1),waiting(1),initial);
            for k=2:n
                try
                    [state,out(k,:)]=advance(state,inputs(k,:),mass(k),waiting(k),dt,steps);
                    completed=k;
                catch exception
                    if ~strcmp(exception.identifier,'DirectMTU:Equilibrium'),rethrow(exception);end
                    break
                end
            end
        end
    end
end

function build_solver
folder=fileparts(mfilename('fullpath'));
source=fullfile(folder,'private','native','direct_mtu_mex.c');
header=fullfile(folder,'private','native','mtu_mechanics.h');
target=fullfile(folder,'private',['direct_mtu_mex.' mexext]);
binary=dir(target);dependencies=[dir(source);dir(header)];
if isempty(binary) || binary.datenum<max([dependencies.datenum])
    clear direct_mtu_mex
    mex('-R2017b','-O','-outdir',fullfile(folder,'private'), ...
        '-output','direct_mtu_mex',source);
end
end

function [s,out]=initial_state(h,cfg,a,mass,waiting,initial)
s=zeros(34,1);s(1:2)=initial(:);s(7)=double(waiting>.5);
command=a(:)';
if s(7)>0
    [command,s(12),s(13)]=standing_request(h,cfg,s(1),s(2),mass,0);
end
units=[h.back h.abdomen];
for u=1:2
    s(13+u)=clip(command(u),0,1);
    s(2+u)=s(13+u)^units(u).activationExponent;
    f=fiber(units(u),h,s(1),s(2+u),0,.001,true);
    s=store_fiber(s,u,f,h);
end
s(31:34)=[a(:);mass;waiting];out=state_output(h,s);
end

function [s,out]=advance_state(h,cfg,s,a,massEnd,waitEnd,dt,steps)
units=[h.back h.abdomen];step=dt/steps;
for k=1:steps
    fraction=k/steps;mass=s(33)+fraction*(massEnd-s(33));
    input=clip(s(31:32)'+fraction*(a(:)'-s(31:32)'),0,1);
    command=input;request=0;correction=0;
    if waitEnd>.5
        s(7)=1+(s(7)-1)*exp(-step/cfg(7));
        s(8)=clip(s(8)+step*(cfg(1)-s(1)),-.5,.5);
        [command,request,correction]=standing_request(h,cfg,s(1),s(2),mass,s(8));
    else
        s(7)=0; % Standing feedback exits at the first action sample.
    end
    s(12)=s(7)*request;s(13)=s(7)*correction;
    for u=1:2
        s(13+u)=s(7)*command(u)+(1-s(7))*input(u);
        target=clip(s(13+u),0,1)^units(u).activationExponent;
        s(2+u)=target+(s(2+u)-target)*exp(-step/units(u).activationTime);
    end
    % Solve the same compliant-fiber and joint momentum equations together.
    inertia=mass*h.loadLever^2+4*h.bodyMass*h.bodyLever^2/3;
    qold=s(1);vold=s(2);guess=qold+step*vold;
    for iteration=1:9
        back=fiber(h.back,h,guess,s(3),s(5),step,false);
        abdomen=fiber(h.abdomen,h,guess,s(4),s(6),step,false);
        torque=h.back.capacity*back(2)-h.abdomen.capacity*abdomen(2);
        residual=inertia*(guess-qold-step*vold)/step^2 ...
            +h.jointDamping*(guess-qold)/step-torque-gravity(h,guess,mass);
        if abs(residual)<1e-7,break,end
        epsilon=1e-6;
        bp=fiber(h.back,h,guess+epsilon,s(3),s(5),step,false);
        ap=fiber(h.abdomen,h,guess+epsilon,s(4),s(6),step,false);
        torquePlus=h.back.capacity*bp(2)-h.abdomen.capacity*ap(2);
        derivative=inertia/step^2+h.jointDamping/step-(torquePlus-torque)/epsilon ...
            -(gravity(h,guess+epsilon,mass)-gravity(h,guess,mass))/epsilon;
        % A Newton increment bound protects iteration, not the joint angle.
        guess=guess-clip(residual/derivative,-.1,.1);
    end
    s(1)=guess;s(2)=(guess-qold)/step;
    back=fiber(h.back,h,guess,s(3),s(5),step,false);
    abdomen=fiber(h.abdomen,h,guess,s(4),s(6),step,false);
    s=store_fiber(s,1,back,h);s=store_fiber(s,2,abdomen,h);
    s(28)=inertia*(s(2)-vold)/step-h.back.capacity*s(16) ...
        +h.abdomen.capacity*s(17)-gravity(h,guess,mass)+h.jointDamping*s(2);
    if ~all(isfinite(s)) || sum(abs(s(26:27)))>1e-4 || abs(s(28))>1
        error('DirectMTU:Equilibrium','Fiber/joint equilibrium failed; no output was clipped.');
    end
end
s(31:34)=[a(:);massEnd;waitEnd];out=state_output(h,s);
end

function [command,desired,correction]=standing_request(h,cfg,q,v,mass,integral)
units=[h.back h.abdomen];minimum=zeros(1,2);
for u=1:2
    f=fiber(units(u),h,q,0,0,.001,true);
    minimum(u)=units(u).capacity*f(2);
end
correction=cfg(5)*tanh((cfg(2)*(cfg(1)-q)-cfg(3)*v+cfg(4)*integral)/cfg(5));
desired=-gravity(h,q,mass)+correction;
difference=desired-minimum(1)+minimum(2);
targets=minimum+cfg(6)+[max(difference,0) max(-difference,0)];
command=[excitation(h.back,h,q,targets(1)) excitation(h.abdomen,h,q,targets(2))];
end

function out=state_output(h,s)
back=h.back.capacity;abdomen=h.abdomen.capacity;
out=[s(1:2)',back*s(16)-abdomen*s(17),back*s(16),abdomen*s(17), ...
    back*s(18),back*s(20),back*.03*s(22)*s(24)/s(5), ...
    abdomen*s(19),abdomen*s(21),abdomen*.03*s(23)*s(25)/s(6), ...
    s(5:6)',s(3:4)',max(abs(s(26:27))),s(28),s(14:15)',s(7), ...
    s(12:13)',s(8),0,0,0,gravity(h,s(1),s(33))];
end

function h = physical_parameters(p)
% Parameters keep their original order in the MAT files. Use names below.
h.bodyMass = p(15);
h.bodyLever = p(16);
h.loadLever = p(13);
h.loadOffset = p(14);
h.jointDamping = p(4);
h.pennationSine = sin(p(24));
h.back = muscle(p(1), p(2), p(3), p(5), p(6), p(7), p(8), p(17), p(19), 1);
h.abdomen = muscle(p(9), p(11), p(12), p(21), p(10), p(22), p(23), p(18), p(20), -1);
end

function unit = muscle(capacity, ratio, referenceLength, vmax, exponent, strain, tau, arm, lengthMT, direction)
unit = struct('capacity',capacity, 'ratio',ratio, 'referenceLength',referenceLength, ...
    'maximumVelocity',vmax, 'activationExponent',exponent, 'tendonStrain',strain, ...
    'activationTime',tau, 'momentArm',arm, 'referenceMTLength',lengthMT, ...
    'optimalLength',arm/ratio, 'direction',direction);
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

function s = store_fiber(s, u, f, h)
s(4+u) = sqrt(f(1)^2+h.pennationSine^2);
s(15+u) = f(2); s(17+u) = f(3); s(19+u) = f(4);
s(21+u) = f(5); s(23+u) = f(1); s(25+u) = f(6);
end

function torque = gravity(h, q, mass)
torque = -9.81*h.bodyMass*h.bodyLever*cos(q) ...
    -9.81*mass*h.loadLever*cos(q-h.loadOffset);
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
