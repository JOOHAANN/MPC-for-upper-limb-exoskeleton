% Refresh publication layout only; the verified closed-loop data are reused.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','NMPC'));
for mode={'actuator_matched','effort_matched'}
    out=fullfile(root,'Demo','NMPC','output','pid_nmpc_528',mode{1});
    data=load(fullfile(out,'pid_nmpc_comparison_528.mat'));
    plot_pid_nmpc_comparison(data.results,out,false);
end
fprintf('Publication figures refreshed from verified traces.\n');
