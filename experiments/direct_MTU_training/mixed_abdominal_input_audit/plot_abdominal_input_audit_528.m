clc
close all

%% Distinguish experimental provenance from model-generated abdominal drive
% This figure audits the mixed-input candidate identified in the JSON report.
% It is not a plot of directly measured internal muscle activation.
folder=fileparts(mfilename('fullpath'));
output=folder;
data=load(fullfile(output,'abdominal_audit_curves.mat'));curves=data.curves;
fig=figure('Color','w','Position',[50 40 1550 900]);
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for source=1:2
    for k=1:3
        item=curves{(source-1)*3+k};ax=nexttile(layout);hold(ax,'on');
        if source==1
            plot(ax,item.time,item.recordedInterpolated,'Color',[.65 .65 .65],'LineWidth',.9);
            plot(ax,item.time,item.experimentalStored,'k--','LineWidth',1.2);
            plot(ax,item.time,item.preparedCommand,'Color',[0 .4 .72],'LineWidth',1.4);
            plot(ax,item.time,item.actualActivation,'Color',[.1 .55 .3],'LineWidth',1.2);
            names={'CEINMS rect abd r, interpolated','Stored experimental-derived input', ...
                'Mixed MTU excitation command','MTU activation state'};
        else
            plot(ax,item.time,item.experimentalStored,'k--','LineWidth',1.1);
            plot(ax,item.time,item.sourceBase,'Color',[.53 .19 .65],'LineWidth',1.3);
            plot(ax,item.time,item.preparedCommand,'Color',[.9 .4 .08],'LineWidth',1.4);
            plot(ax,item.time,item.actualActivation,'Color',[.1 .55 .3],'LineWidth',1.2);
            names={'Experimental-derived abdominal input','Auxiliary predicted abdominal input', ...
                'Mixed MTU excitation command','MTU activation state'};
        end
        labels={'Experimental branch','MAP branch'};
        title(ax,sprintf('%s | %g kg | code %g',labels{source},item.loadKg,item.speedCode));
        ylabel(ax,'Dimensionless input / activation');xlabel(ax,'Time within recorded action (s)');
        grid(ax,'on');box(ax,'on');set(ax,'FontName','Times New Roman','FontSize',10);
        if k==1,legend(ax,names,'Location','best','FontSize',8);end
    end
end
title(layout,'Subject 528 | Abdominal signal provenance and actual model drive','FontName','Times New Roman','FontSize',16);
subtitle(layout,'Audit of the mixed-input candidate. Experimental provenance does not make the final mixed drive a measured activation.');
stem=fullfile(output,'abdominal_input_provenance_528');
exportgraphics(fig,[stem '.png'],'Resolution',180);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig(fig,[stem '.fig']);
