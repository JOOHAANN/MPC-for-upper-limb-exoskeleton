function value=no_posture_nmpc_signature(folder,mainSource)
%NO_POSTURE_NMPC_SIGNATURE Reuse NMPC when only PID calibration changes.
if nargin<2, mainSource=fullfile(folder,'offline_PID_NMPC_no_human_PD_528.m'); end
code=fileread(mainSource);
first=strfind(code,'function run=simulate_controller(');
last=strfind(code,'function [u,m]=pid_update(');
shared=code(first:last-1);
first=strfind(code,'function [u,plan,flag]=nmpc_update(');
last=strfind(code,'function check_run(');
predictive=code(first:last-1);
engine=java.security.MessageDigest.getInstance('SHA-256');
engine.update(uint8([shared predictive]));
files={'bse_without_posture.m','private/native/bse_no_posture_mex.c','nmpc_forecast.m', ...
    '../BSE/bse_inputs.m','../BSE/private/native/hill_plant.h', ...
    '../BSE/private/native/hill_mechanics.h','../BSE/private/native/hill_neural_control.h'};
for k=1:numel(files)
    f=fopen(fullfile(folder,files{k}),'rb'); bytes=fread(f,inf,'*uint8'); fclose(f); engine.update(bytes);
end
value=lower(reshape(dec2hex(typecast(engine.digest(),'uint8'),2)',1,[]));
end
