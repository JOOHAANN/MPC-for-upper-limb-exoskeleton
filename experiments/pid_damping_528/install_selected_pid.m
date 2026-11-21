% Publish the selected PID without changing either NMPC or BSE model.
folder=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','NMPC');
selection=jsondecode(fileread(fullfile(folder,'selection.json')));
g=selection.gains;
pid=struct('Kp',g(1),'Ki',g(2),'Kd',g(3),'derivativeTime',g(4));
save(fullfile(demo,'data','pid_damped_528.mat'),'pid','-v7');
name=fullfile(demo,'data','threshold_tracking_settings_528.mat');
settings=load(name);
settings.pid=pid;
settings.protocol.pidGainsUnchanged=false;
settings.protocol.samePIDAcrossCasesAndRQs=true;
settings.protocol.pidTuning='Damping-focused PID with rise-rate-aware conditional integration';
save(name,'-struct','settings','-v7');
disp(pid);
