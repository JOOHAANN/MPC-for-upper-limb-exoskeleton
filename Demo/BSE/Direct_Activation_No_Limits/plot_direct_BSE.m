function plot_direct_BSE(subject,cases,folder,showFigures)
%PLOT_DIRECT_BSE Raw, unshifted torque/theta and feedback audit for nine cases.
colors=[0 .40 .72;.53 .19 .65];
types={'torque','theta','muscles','feedback','abdominal_input','back_input'};
for kind=1:numel(types)
    visible='off';if showFigures,visible='on';end
    fig=figure('Color','w','Visible',visible,'Position',[50 50 1500 970]);
    layout=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
    for k=1:9
        item=cases{k};ax=nexttile(layout);hold(ax,'on');h=gobjects(1,0);
        switch types{kind}
            case 'back_input'
                h(1)=plot(ax,item.time,item.inputs(:,1,1),'Color',colors(1,:),'LineWidth',1.2);
                h(2)=plot(ax,item.time,item.rawMAPBack,'--','Color',[.5 .5 .5],'LineWidth',1.1);
                h(3)=plot(ax,item.time,item.inputs(:,1,2),'Color',colors(2,:),'LineWidth',1.2);
                names={'Legacy experimental back','Original fused MAP','Load-calibrated MAP input'};
                ylabel(ax,'Back activation input');
            case 'abdominal_input'
                h(1)=plot(ax,item.time,item.inputs(:,2,1),'Color',colors(1,:),'LineWidth',1.2);
                h(2)=plot(ax,item.time,item.rawMAPAbdomen,'--','Color',[.5 .5 .5],'LineWidth',1.1);
                h(3)=plot(ax,item.time,item.inputs(:,2,2),'Color',colors(2,:),'LineWidth',1.2);
                names={'Legacy experimental abdomen','Raw abdominal MAP','Load-calibrated MAP input'};
                ylabel(ax,'Abdominal excitation input');
            case 'torque'
                h(1)=plot(ax,item.time,item.referenceTorque,'r--','LineWidth',1.2);
                for source=1:2,h(source+1)=plot(ax,item.time,item.outputs(:,3,source),'Color',colors(source,:),'LineWidth',1.2);end
                names={'CEINMS reference','Experimental input','MAP input'};ylabel(ax,'Net muscle torque (N m)');
            case 'theta'
                h(1)=plot(ax,item.time,item.referenceTheta,'r--','LineWidth',1.2);
                for source=1:2,h(source+1)=plot(ax,item.time,item.outputs(:,1,source),'Color',colors(source,:),'LineWidth',1.2);end
                names={'Original IK reference','Experimental input','MAP input'};ylabel(ax,'Trunk angle (rad)');
            case 'muscles'
                for source=1:2
                    h(2*source-1)=plot(ax,item.time,item.outputs(:,4,source),'Color',colors(source,:),'LineWidth',1.2);
                    h(2*source)=plot(ax,item.time,item.outputs(:,5,source),'--','Color',colors(source,:),'LineWidth',1.1);
                end
                names={'Experimental back','Experimental abdomen','MAP back','MAP abdomen'};ylabel(ax,'Muscle torque (N m)');
            case 'feedback'
                for source=1:2,h(source)=plot(ax,item.time,item.outputs(:,22,source),'Color',colors(source,:),'LineWidth',1.2);end
                h(3)=plot(ax,item.time,item.outputs(:,26,1),'k--');
                movementFeedback=item.outputs(:,22,1);movementFeedback(item.waiting>.5)=0;
                h(4)=plot(ax,item.time,movementFeedback,':','Color',[.2 .55 .2],'LineWidth',1.3);
                names={'Experimental standing PD request','MAP standing PD request', ...
                    'Angle boundary: measured zero','PD during movement: measured zero'};
                ylabel(ax,'Correction request (N m)');
        end
        heading=sprintf('%g kg | code %g',item.loadKg,item.speedCode);
        if item.speedCode<0 && isfield(item,'cadenceBPM')
            heading=sprintf('%g kg | %.2f BPM | bout %d',item.loadKg,item.cadenceBPM,k);
        end
        title(ax,heading,'FontWeight','normal');
        xlabel(ax,'Time from recorded action onset (s)');grid(ax,'on');box(ax,'on');
        xlim(ax,[item.time(1) item.time(end)]);set(ax,'FontName','Times New Roman','FontSize',9);
        if k==1,handles=h;legendNames=names;end
        if any(~isfinite(item.outputs(:,1,:)),'all')
            subtitle(ax,'Simulation stopped: missing output is not interpolated','FontSize',8);
        end
    end
    key=legend(handles,legendNames,'NumColumns',min(4,numel(legendNames)),'FontSize',10);key.Layout.Tile='south';
    title(layout,sprintf('%s | Direct-activation dual-MTU BSE | %s',subject.shortName,types{kind}), ...
        'FontName','Times New Roman','FontSize',16,'Interpreter','none');
    policy='Transferred coefficients, not refitted';
    if isfield(subject,'calibrationProtocol'),policy='New subject-specific calibration | Nine-record fitting replay';end
    subtitle(layout,['Standing balance only | No movement PD, assistance or angle boundaries | ' policy]);
    stem=fullfile(folder,sprintf('%s_%s',subject.shortName,types{kind}));
    format_demo_figure(fig);
    exportgraphics(fig,[stem '.png'],'Resolution',180);
    exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig(fig,[stem '.fig']);
    if ~showFigures,close(fig);end
end
if strcmp(subject.shortName,'528'),overview_528(cases,folder,showFigures);end
end

function overview_528(cases,folder,showFigures)
% Three representative loads; all-nine plots and metrics remain available.
visible='off';if showFigures,visible='on';end
fig=figure('Color','w','Visible',visible,'Position',[40 30 1500 1300]);
layout=tiledlayout(fig,4,3,'TileSpacing','compact','Padding','compact');
colors=[0 .40 .72;.53 .19 .65];
for row=1:4
    for k=1:3
        item=cases{k};ax=nexttile(layout);hold(ax,'on');
        switch row
            case 1
                plot(ax,item.time,item.referenceTorque,'r--','LineWidth',1.2);
                for s=1:2,plot(ax,item.time,item.outputs(:,3,s),'Color',colors(s,:),'LineWidth',1.3);end
                ylabel(ax,'Net torque (N m)');
                if k==1,legend(ax,'CEINMS reference','Experimental','MAP','FontSize',8, ...
                        'Location','northoutside','NumColumns',3);end
            case 2
                plot(ax,item.time,item.referenceTheta,'r--','LineWidth',1.2);
                for s=1:2,plot(ax,item.time,item.outputs(:,1,s),'Color',colors(s,:),'LineWidth',1.3);end
                ylabel(ax,'Trunk angle (rad)');
            case 3
                for s=1:2
                    plot(ax,item.time,item.outputs(:,4,s),'Color',colors(s,:),'LineWidth',1.2);
                    plot(ax,item.time,item.outputs(:,5,s),'--','Color',colors(s,:),'LineWidth',1.2);
                end
                ylabel(ax,'Muscle torque (N m)');
                if k==1,legend(ax,'Exp. back','Exp. abdomen','MAP back','MAP abdomen','FontSize',8, ...
                        'Location','northoutside','NumColumns',2);end
            case 4
                for s=1:2,plot(ax,item.time,item.outputs(:,22,s),'Color',colors(s,:),'LineWidth',1.1);end
                plot(ax,item.time,item.outputs(:,26,1),'k--','LineWidth',1.2);
                moving=item.outputs(:,22,1);moving(item.waiting>.5)=0;
                plot(ax,item.time,moving,':','Color',[.2 .55 .2],'LineWidth',1.3);
                ylabel(ax,'PD request (N m)');
                if k==1,legend(ax,'Exp. standing','MAP standing','Boundary: 0', ...
                    'Movement PD: 0','FontSize',8,'Location','northoutside','NumColumns',2);end
        end
        if row==1,title(ax,sprintf('%g kg | code %g',item.loadKg,item.speedCode));end
        if row==4,xlabel(ax,'Time from recorded action onset (s)');end
        set(ax,'FontName','Times New Roman','FontSize',10,'Box','on');grid(ax,'on');
        xlim(ax,[item.time(1) item.time(end)]);
    end
end
title(layout,'Subject 528 | Direct-activation dual-MTU BSE | Representative 5/10/15-kg bouts', ...
    'FontName','Times New Roman','FontSize',16);
subtitle(layout,'No assistance, movement PD or angle boundaries. Standing balance retained. All curves are raw model outputs.');
stem=fullfile(folder,'528_overview');
format_demo_figure(fig);
exportgraphics(fig,[stem '.png'],'Resolution',180);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig(fig,[stem '.fig']);
if ~showFigures,close(fig);end
end
