function plant = bse_without_posture(model,balance,standingAngle)
%BSE_WITHOUT_POSTURE Disable human posture feedback in actual and predicted BSE.
% Physical Hill parameters, gravity, passive damping, soft stops and the
% assistance-dependent excitation response remain unchanged. There is no
% human position/velocity/integral tracking torque, including during waits.
folder=fileparts(mfilename('fullpath'));
source=fullfile(folder,'private','native','bse_no_posture_mex.c');
target=fullfile(folder,'private',['bse_no_posture_mex.' mexext]);
binary=dir(target); code=dir(source);
if isempty(binary) || binary.datenum<code.datenum
    clear bse_no_posture_mex
    mex('-R2017b','-O','-outdir',fullfile(folder,'private'), ...
        '-output','bse_no_posture_mex',source);
end
model.neuralParameters(4:5)=0; % Movement posture Kp and Kd.
balance(1:3)=0;               % Standing posture Kp, Kd and Ki.
p=model.parameters(:)'; r=model.neuralParameters(:)';
cfg=[standingAngle balance(:)'];
plant=bse_model(model,balance,standingAngle,true);
plant.initialize=@initialize; plant.step=@advance;
plant.preview=@preview; plant.simulate=@simulate;
plant.neuralParameters=r; plant.standingParameters=cfg;

    function [state,out]=initialize(a,mass,waiting,initial)
        [state,out]=bse_no_posture_mex('init',p,r,cfg,a,mass,waiting,initial);
    end
    function [state,out]=advance(state,a,mass,waiting,u,dt,steps)
        [state,out]=bse_no_posture_mex('step',p,r,cfg,state,a,mass,waiting,u,dt,steps);
    end
    function out=simulate(a,mass,waiting,initial,dt,steps)
        out=bse_no_posture_mex('rollout',p,r,cfg,a,mass,waiting,initial,dt,steps);
    end
    function out=preview(state,controls,forecast,options)
        out=bse_no_posture_mex('preview',p,r,cfg,state,controls, ...
            forecast.activation,forecast.mass,forecast.waiting, ...
            options.predictionDt,options.predictionSubsteps,options.rateNmPerSec);
    end
end
