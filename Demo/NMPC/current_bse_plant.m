function plant=current_bse_plant(model,balance,standingAngle,humanPD,actuator)
%CURRENT_BSE_PLANT Latest calibrated dual MTUs plus explicitly gated human PD.
% Human PD corrects back activation only when applied assistance is nonzero.
% Its target is the current causal zero-assistance shadow, not experimental IK.
% The shadow has identical physical coefficients and activation input.
% State: 1:34 physical plant, 35:68 unassisted shadow, 69 PD correction memory.
p=model.parameters(:)';assert(numel(p)==24);
cfg=[standingAngle balance(:)'];pd=humanPD(:)';hardware=actuator(:)';
assert(numel(pd)==5 && numel(hardware)==2);
build_solver;
plant.initialize=@initialize;plant.step=@advance;plant.preview=@preview;plant.simulate=@simulate;
plant.parameters=p;plant.balance=cfg;plant.humanPD=pd;plant.actuator=hardware;
plant.physics=struct('back',struct('capacity',p(1),'momentArm',p(17)), ...
    'abdomen',struct('capacity',p(9),'momentArm',p(18)));
    function [state,out]=initialize(input,mass,waiting,initial)
        [state,out]=current_mtu_control_mex('init',p,cfg,pd,hardware, ...
            double(input),double(mass),double(waiting),double(initial));
    end
    function [state,out]=advance(state,input,mass,waiting,command,dt,steps)
        [state,out]=current_mtu_control_mex('step',p,cfg,pd,hardware,state, ...
            double(input),double(mass),double(waiting),double(command),double(dt),double(steps));
    end
    function out=preview(state,controls,forecast,options)
        assert(options.capNm==hardware(1) && options.rateNmPerSec==hardware(2), ...
            'NMPC:Actuator','Prediction and applied actuator constraints must agree.');
        out=current_mtu_control_mex('preview',p,cfg,pd,hardware,state,double(controls), ...
            double(forecast.activation),double(forecast.mass),double(forecast.waiting), ...
            options.predictionDt,double(options.predictionSubsteps),options.rateNmPerSec);
    end
    function out=simulate(inputs,mass,waiting,initial,dt,steps)
        n=size(inputs,1);out=zeros(n,27);
        [state,out(1,:)]=initialize(inputs(1,:),mass(1),waiting(1),initial);
        for k=2:n
            [state,out(k,:)]=advance(state,inputs(k,:),mass(k),waiting(k),0,dt,steps);
        end
    end
end

function build_solver
folder=fileparts(mfilename('fullpath'));
source=fullfile(folder,'private','native','current_mtu_control.c');
header=fullfile(folder,'..','BSE','Direct_Activation_No_Limits','private','native','mtu_mechanics.h');
target=fullfile(folder,'private',['current_mtu_control_mex.' mexext]);
binary=dir(target);dependencies=[dir(source);dir(header)];
if isempty(binary) || binary.datenum<max([dependencies.datenum])
    clear current_mtu_control_mex
    mex('-R2017b','-O','-outdir',fullfile(folder,'private'),'-output','current_mtu_control_mex',source);
end
end
