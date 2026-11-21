function plot_current_PD(results,folder,showFigures)
%PLOT_CURRENT_PD Audit actual assistance, neural corrections and mode gates.
visible='off';if showFigures,visible='on';end
fig=figure('Color','w','Visible',visible,'Position',[40 30 1450 1100]);
layout=tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
q=results.groups{1}.sequence;t=q.time;names={'Experimental activation','MAP activation'};
for source=1:2
    pair=results.groups{1}.runs{source};p=pair.pid.actual;n=pair.nmpc.actual;moving=q.waiting<.5;
    ax=nexttile(layout,source);hold(ax,'on');
    plot(ax,t,p(:,25),'Color',[.85 .4 .08]);plot(ax,t,n(:,25),'Color',[0 .36 .7]);yline(ax,80,'k--');
    ylabel(ax,'Applied assistance (N m)');title(ax,names{source});
    legend(ax,'PID','NMPC','Shared actuator cap','Location','northwest');
    ax=nexttile(layout,2+source);hold(ax,'on');
    pd=p(:,18)-pair.nmpc.inputs(:,1);nd=n(:,18)-pair.nmpc.inputs(:,1);
    pd(~moving)=0;nd(~moving)=0;
    plot(ax,t,pd,'Color',[.85 .4 .08]);plot(ax,t,nd,'Color',[0 .36 .7]);
    ylabel(ax,'Back input correction');legend(ax,'PID plant human PD','NMPC plant human PD','Location','northwest');
    ax=nexttile(layout,4+source);hold(ax,'on');
    plot(ax,t,pair.baseline(:,5),'k:');plot(ax,t,p(:,5),'Color',[.85 .4 .08]);plot(ax,t,n(:,5),'Color',[0 .36 .7]);
    ylabel(ax,'Abdominal torque (N m)');legend(ax,'Unassisted','With PID','With NMPC','Location','northwest');
    ax=nexttile(layout,6+source);hold(ax,'on');
    stairs(ax,t,double(moving & abs(p(:,25))>1e-9),'Color',[.85 .4 .08]);
    stairs(ax,t,double(moving & abs(n(:,25))>1e-9),'Color',[0 .36 .7]);
    stairs(ax,t,double(~moving),'k:');ylim(ax,[-.05 1.2]);
    ylabel(ax,'Feedback enabled');xlabel(ax,'Continuous time (s)');
    legend(ax,'PID movement PD','NMPC movement PD','Standing/waiting','Location','northwest');
end
for ax=findobj(fig,'Type','axes')'
    grid(ax,'on');box(ax,'on');xlim(ax,[t(1) t(end)]);set(ax,'FontName','Times New Roman','FontSize',9);
end
title(layout,'Subject 528 | Actual-assistance-gated human PD | Latest calibrated dual-MTU model');
subtitle(layout,'Movement PD corrects back activation only. Its current target is the causal unassisted shadow. Standing balance remains active during waits.');
format_demo_figure(fig);
exportgraphics(fig,fullfile(folder,'528_assistance_PD_audit.png'),'Resolution',180);
exportgraphics(fig,fullfile(folder,'528_assistance_PD_audit.pdf'),'ContentType','vector');
savefig(fig,fullfile(folder,'528_assistance_PD_audit.fig'));
if ~showFigures,close(fig);end
end
