function redraw_RQ1_528(showFigures)
%REDRAW_RQ1_528 Refresh Fig. 7 directly from the verified MATLAB simulations.
if nargin<1,showFigures=true;end
folder=fileparts(mfilename('fullpath'));addpath(folder,fileparts(folder));
output=fullfile(folder,'output','current_BSE_528');
file=fullfile(output,'current_BSE_NMPC_528.mat');
assert(isfile(file),'Run offline_PID_NMPC_comparison_528 first.');
data=load(file,'results');results=data.results;
c=struct('nominal',[.8 .15 .15],'pid',[.88 .42 .08],'nmpc',[0 .36 .70], ...
    'pidAssist',[.08 .48 .22],'nmpcAssist',[.55 .15 .65]);
names={'Experimental activation','MAP activation'};suffixes={'experimental','MAP'};
for source=1:2
    plot_RQ1_continuous(results,source,c,output,suffixes{source},names{source},showFigures);
end
end
