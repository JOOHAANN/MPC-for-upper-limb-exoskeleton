clc
clearvars -except Group_index Source_index
close all

%% Unassisted BSE: distinguish human standing and movement feedback
% No external controller is called. Exoskeleton torque is zero throughout.
% Trace the original calibrated model at its existing 1-ms internal step.
% PD terms request muscle recruitment; they are not additional joint torques.
if ~exist('Group_index','var'), Group_index=1; end
if ~exist('Source_index','var'), Source_index=1; end
folder=fileparts(mfilename('fullpath'));
addpath(folder);
subjects=bse_subjects('528');subject=subjects{1};
q=subject.sequences{Group_index};model=subject.models{Source_index};
plant=bse_model(model,subject.balanceController,q.standingAngle,true);
inputs=bse_inputs(q,model,Source_index);
p=model.parameters;r=model.neuralParameters;
balance=subject.balanceController;
Kp=balance(1);Kd=balance(2);Ki=balance(3);limit=balance(4);tau=balance(6);
outputFolder=fullfile(folder,'output','pd_diagnostic_528');
if ~isfolder(outputFolder),mkdir(outputFolder);end

%% Replay the unassisted model and separate the two correction requests
steps=10;h=q.dt/steps;n=numel(q.time);
actual=zeros(n,27);standing=zeros(n,1);movement=zeros(n,1);
[state,actual(1,:)]=plant.initialize(inputs(1,:),q.mass(1),q.waiting(1),q.initialState);
[standing(1),movement(1)]=correction_terms(actual(1,:),state,actual(1,1:2),r,balance,q.standingAngle);
for k=2:n
    for inner=1:steps
        f=inner/steps;
        input=inputs(k-1,:)+f*(inputs(k,:)-inputs(k-1,:));
        mass=q.mass(k-1)+f*(q.mass(k)-q.mass(k-1));
        waiting=q.waiting(k-1)+f*(q.waiting(k)-q.waiting(k-1));
        previousPose=state(1:2)';
        [state,row]=plant.step(state,input,mass,waiting,0,h,1);
        [standing(k),movement(k)]=correction_terms(row,state,previousPose,r,balance,q.standingAngle);
    end
    actual(k,:)=row;
end
% Confirm that instrumentation does not change the existing BSE trajectory.
ordinary=plant.simulate(inputs,q.mass,q.waiting,q.initialState,q.dt,steps);
replayError=max(abs(actual(:,1:19)-ordinary(:,1:19)),[],'all');
requestError=max(abs(standing+movement-actual(:,22)));
assert(replayError<1e-6,'Diagnostic replay changed the original BSE.');
assert(requestError<1e-7,'Separated PD requests do not sum to the recorded request.');
assert(max(abs(actual(:,25)))==0,'This diagnostic must remain unassisted.');

%% Exact commanded mode intervals; blend transitions remain continuous
changes=[1;find(diff(q.waiting)~=0)+1;n+1];
rows=struct([]);
for j=1:numel(changes)-1
    first=changes(j);last=changes(j+1)-1;
    stop=q.time(end);if changes(j+1)<=n,stop=q.time(changes(j+1));end
    mode="Movement";if q.waiting(first)>.5,mode="Standing / waiting";end
    rows=[rows;struct('Mode',mode,'StartSeconds',q.time(first),'EndSeconds',stop)]; %#ok<AGROW>
end
intervals=struct2table(rows);disp(intervals);
fprintf('Standing target %.3f rad; Kp %.6g; Kd %.6g; Ki %.6g; correction-request limit %.3f N m; blend time constant %.3f s.\n', ...
    q.standingAngle,Kp,Kd,Ki,limit,tau);
fprintf('Movement Kp %.6g; Kd %.6g; movement recruitment blend %.6g.\n',r(4),r(5),r(7));
fprintf('Maximum replay difference %.3g; request decomposition error %.3g N m.\n',replayError,requestError);

%% Plot the original unassisted BSE, not a PID/NMPC-controlled trajectory
labels={'Experimental activation','MAP activation'};
overview=make_figure(q,actual,standing,movement,r,labels{Source_index});
save_plot(overview,outputFolder,sprintf('unassisted_PD_schedule_group%d_source%d',Group_index,Source_index));
zoom=make_figure(q,actual,standing,movement,r,labels{Source_index});
axesHandles=findall(zoom,'Type','axes');
for j=1:numel(axesHandles),xlim(axesHandles(j),[0 min(14,q.time(end))]);end
save_plot(zoom,outputFolder,sprintf('unassisted_PD_first_action_group%d_source%d',Group_index,Source_index));
trace=table(q.time,q.waiting,actual(:,20), ...
    (1-actual(:,20))*r(7),standing,movement,actual(:,22), ...
    actual(:,1),actual(:,4),actual(:,5),actual(:,25), ...
    'VariableNames',{'TimeSeconds','WaitingFlag','StandingBlend','MovementFeedbackBlend', ...
    'StandingCorrectionRequestNm','MovementCorrectionRequestNm','TotalCorrectionRequestNm', ...
    'ThetaRad','BackMuscleTorqueNm','AbdominalTorqueNm','AssistanceTorqueNm'});
writetable(trace,fullfile(outputFolder,sprintf('unassisted_PD_trace_group%d_source%d.csv',Group_index,Source_index)));
writetable(intervals,fullfile(outputFolder,sprintf('PD_mode_intervals_group%d.csv',Group_index)));
save(fullfile(outputFolder,sprintf('unassisted_PD_group%d_source%d.mat',Group_index,Source_index)), ...
    'actual','standing','movement','intervals','q','model','balance','replayError','requestError','-v7');
fprintf('Unassisted BSE diagnostic: %s\n',outputFolder);

function [standing,movement]=correction_terms(out,state,oldPose,r,balance,target)
b=out(20);q=oldPose(1);v=oldPose(2);
standing=0;movement=0;
if b>1e-8
    standing=b*balance(4)*tanh((balance(1)*(target-q)-balance(2)*v ...
        +balance(3)*state(8))/balance(4));
end
if b<1-1e-8,movement=(1-b)*(r(4)*(out(24)-q)-r(5)*v);end
end

function fig=make_figure(q,a,standing,movement,r,source)
fig=figure('Color','w','Position',[80 40 1450 1000]);
layout=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout);hold(ax,'on');shade_waits(ax,q,[-.06 1.12]);
stairs(ax,q.time,q.waiting,'k--','LineWidth',1.15);
plot(ax,q.time,a(:,20),'Color',[.82 .39 .05],'LineWidth',1.7);
plot(ax,q.time,(1-a(:,20))*r(7),'Color',[.03 .37 .69],'LineWidth',1.5);
ylim(ax,[-.06 1.12]);ylabel(ax,'Recruitment weight');
legend(ax,'Waiting flag','Standing feedback: b','Movement feedback: (1-b) x blend', ...
    'Location','eastoutside');title(ax,'(a) Human feedback scheduling');format_axis(ax,q.time);
ax=nexttile(layout);hold(ax,'on');
plot(ax,q.time,standing,'Color',[.82 .39 .05],'LineWidth',1.4);
plot(ax,q.time,movement,'Color',[.03 .37 .69],'LineWidth',1.25);
plot(ax,q.time,a(:,22),'k:','LineWidth',1.05);yline(ax,0,'Color',[.6 .6 .6]);
ylabel(ax,'Correction request (N m)');
legend(ax,'Standing PD request','Movement PD request','Total posture request', ...
    'Location','eastoutside');title(ax,'(b) PD requests to the muscle-recruitment model');format_axis(ax,q.time);
ax=nexttile(layout);hold(ax,'on');
plot(ax,q.time,a(:,1),'Color',[.1 .1 .1],'LineWidth',1.5);
yline(ax,q.standingAngle,'--','Color',[.82 .39 .05],'LineWidth',1);
ylabel(ax,'Trunk angle (rad)');legend(ax,'Unassisted BSE','Standing target: 1.5 rad', ...
    'Location','eastoutside');title(ax,'(c) Original calibrated unassisted motion');format_axis(ax,q.time);
ax=nexttile(layout);hold(ax,'on');
plot(ax,q.time,a(:,4),'Color',[.03 .37 .69],'LineWidth',1.4);
plot(ax,q.time,a(:,5),'Color',[.48 .16 .63],'LineWidth',1.1);
plot(ax,q.time,a(:,25),'k--','LineWidth',1.4);
ylabel(ax,'Actual torque (N m)');xlabel(ax,'Time (s)');
legend(ax,'Back muscle','Abdominal muscle','Exoskeleton: zero throughout', ...
    'Location','eastoutside');title(ax,'(d) Hill-muscle outputs and zero external assistance');format_axis(ax,q.time);
title(layout,['Unassisted BSE only | Subject 528 | ' source],'FontName','Times New Roman','FontSize',17);
subtitle(layout,'Shaded intervals: commanded standing/waiting mode. PD correction requests are not extra mechanical joint torques.', ...
    'FontName','Times New Roman','FontSize',10);
end

function shade_waits(ax,q,limits)
edges=[1;find(diff(q.waiting)~=0)+1;numel(q.time)+1];
for j=1:numel(edges)-1
    first=edges(j);if q.waiting(first)<.5,continue,end
    stop=q.time(end);if edges(j+1)<=numel(q.time),stop=q.time(edges(j+1));end
    patch(ax,[q.time(first) stop stop q.time(first)],limits([1 1 2 2]), ...
        [.92 .95 .91],'EdgeColor','none','HandleVisibility','off');
end
end
function format_axis(ax,time)
set(ax,'FontName','Times New Roman','FontSize',11,'Box','on', ...
    'XGrid','on','YGrid','on','GridAlpha',.12);xlim(ax,[time(1) time(end)]);
end
function save_plot(fig,folder,stem)
exportgraphics(fig,fullfile(folder,[stem '.png']),'Resolution',200);
exportgraphics(fig,fullfile(folder,[stem '.pdf']),'ContentType','vector');
savefig(fig,fullfile(folder,[stem '.fig']));
end
