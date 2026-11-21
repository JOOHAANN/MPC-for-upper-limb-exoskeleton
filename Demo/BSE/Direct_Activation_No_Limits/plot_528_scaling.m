function plot_528_scaling(subject,output,showFigures)
%PLOT_528_SCALING Show the experimental mean and the shared load-only mapping.
if nargin<3,showFigures=true;end
visible='off';if showFigures,visible='on';end
s=subject.MAPScaling;fig=figure('Color','w','Visible',visible,'Position',[40 40 1450 760]);
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
labels={'Back','Abdomen'};
for muscle=1:2
    for load=1:3
        ax=nexttile(layout);hold(ax,'on');h=gobjects(1,5);
        experimental=s.meanExperimental(:,muscle,load);raw=s.meanRawMAP(:,muscle,load);
        mapped=raw*s.gain(load,muscle)+s.offset(load,muscle);
        h(1)=plot(ax,subject.MAP.phase,experimental,'k','LineWidth',1.7);
        h(2)=plot(ax,subject.MAP.phase,raw,'--','Color',[.5 .5 .5],'LineWidth',1.1);
        h(3)=plot(ax,subject.MAP.phase,mapped,'Color',[.53 .19 .65],'LineWidth',1.5);
        h(4)=yline(ax,s.experimentalMinimum(load,muscle),':','Color',[.2 .6 .3]);
        h(5)=yline(ax,s.experimentalMaximum(load,muscle),':','Color',[.2 .6 .3]);
        title(ax,sprintf('%g kg | %s',s.loadKg(load),labels{muscle}));
        subtitle(ax,sprintf('Gain %.3f | Offset %.4f',s.gain(load,muscle),s.offset(load,muscle)));
        xlabel(ax,'Normalized action phase');ylabel(ax,'Activation input');grid(ax,'on');box(ax,'on');
        set(ax,'FontName','Times New Roman','FontSize',10);
        if muscle==1 && load==1
            handles=h;
        end
    end
end
key=legend(handles,'Five-speed experimental mean','Raw MAP mean','Affine-calibrated MAP mean', ...
    'Experimental mean minimum','Experimental mean maximum','NumColumns',3,'FontSize',9);
key.Layout.Tile='south';
title(layout,'Subject 528 | MAP amplitude calibration by load | Both MTUs');
subtitle(layout,{'Average the five phase-normalized experimental curves first, then take their minimum and maximum.', ...
    'The same gain and offset apply to codes 25--65; there is no per-speed amplitude fit.'});
format_demo_figure(fig);
exportgraphics(fig,fullfile(output,'528_load_scaling.png'),'Resolution',180);
exportgraphics(fig,fullfile(output,'528_load_scaling.pdf'),'ContentType','vector');
savefig(fig,fullfile(output,'528_load_scaling.fig'));
if ~showFigures,close(fig);end
end
