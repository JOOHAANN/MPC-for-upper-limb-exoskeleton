function plot_nmpc_results(subject, sequence, runs, folder, showFigures)
%PLOT_NMPC_RESULTS Plot recalculated results without running the controller.
visible = 'off';
if showFigures, visible = 'on'; end
fig = figure('Name',[subject.id ' | NMPC'],'Color','w', ...
    'Position',[60 30 1500 1120],'Visible',visible);
layout = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
colors = [0 .447 .741; .494 .184 .556];
reference = sequence.referenceTorque;
reference(~sequence.torqueValid | sequence.referenceFilled) = NaN;
for source = 1:2
    nominal = runs{source}.nominal;
    actual = runs{source}.actual;
    ax = nexttile(layout,source);
    plot(ax,sequence.time,reference,'--','Color',[.85 .2 .15],'LineWidth',1.3); hold(ax,'on');
    plot(ax,sequence.time,nominal(:,3),':','Color',[.55 .55 .55],'LineWidth',1.3);
    plot(ax,sequence.time,actual(:,3),'Color',colors(source,:),'LineWidth',1.5);
    ylabel(ax,'Net muscle torque (N m)');
    names = {'Experimental activation','MAP activation'};
    title(ax,names{source});
    legend(ax,{'CEINMS net reference','Unassisted BSE','With NMPC'}, ...
        'Location','northoutside','Orientation','horizontal','FontSize',8);
    format_axis(ax,sequence,true);

    ax = nexttile(layout,2+source);
    plot(ax,sequence.time,nominal(:,4),':','Color',[.55 .55 .55],'LineWidth',1.3); hold(ax,'on');
    plot(ax,sequence.time,actual(:,4),'Color',colors(source,:),'LineWidth',1.5);
    limitLine = nan(size(sequence.time));
    for stage = 1:3
        s = sequence.segments{stage};
        index = s.startIndex:s.endIndex;
        limitLine(index) = runs{source}.threshold(index);
        limit = limitLine(index(1));
        patch(ax,[sequence.time(index(1)),sequence.time(index(end)), ...
            sequence.time(index(end)),sequence.time(index(1))], ...
            [limit-5,limit-5,limit,limit],[.7 .85 .7], ...
            'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
    end
    plot(ax,sequence.time,limitLine,'k--','LineWidth',1.3);
    ylabel(ax,'Back muscle torque (N m)');
    legend(ax,{'Unassisted back','Assisted back','75% peak threshold'}, ...
        'Location','northwest','FontSize',8);
    format_axis(ax,sequence,false);

    ax = nexttile(layout,4+source);
    plot(ax,sequence.time,sequence.referenceTheta,'--','Color',[.85 .2 .15],'LineWidth',1.3); hold(ax,'on');
    plot(ax,sequence.time,nominal(:,1),':','Color',[.55 .55 .55],'LineWidth',1.3);
    plot(ax,sequence.time,actual(:,1),'Color',colors(source,:),'LineWidth',1.5);
    ylabel(ax,'Theta (rad)');
    legend(ax,{'Original angle reference','Unassisted BSE','With NMPC'}, ...
        'Location','northoutside','Orientation','horizontal','FontSize',8);
    format_axis(ax,sequence,false);

    ax = nexttile(layout,6+source);
    plot(ax,sequence.time,actual(:,25),'Color',[.13 .55 .32],'LineWidth',1.4); hold(ax,'on');
    yline(ax,subject.nmpc.capNm,'r--','LineWidth',1.2);
    ylabel(ax,'Assistance torque (N m)'); xlabel(ax,'Continuous time (s)');
    legend(ax,{'Applied assistance','Actuator limit'},'Location','northwest','FontSize',8);
    format_axis(ax,sequence,false);
end
title(layout,[subject.id ' | subject-specific NMPC on the same Hill model'], ...
    'Interpreter','none','FontSize',15);
subtitle(layout,'One initialization per source. Black dashed limits apply to back muscle torque; CEINMS is a net-torque estimate.');
save_figure(fig,folder,[subject.shortName '_nmpc'],showFigures);

fig = figure('Name',[subject.id ' | muscle components'],'Color','w', ...
    'Position',[90 90 1450 680],'Visible',visible);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for source = 1:2
    actual = runs{source}.actual; nominal = runs{source}.nominal;
    ax = nexttile(layout,source);
    plot(ax,sequence.time,nominal(:,5),':','Color',[.55 .55 .55],'LineWidth',1.3); hold(ax,'on');
    plot(ax,sequence.time,actual(:,5),'Color',colors(source,:),'LineWidth',1.4);
    title(ax,names{source}); ylabel(ax,'Abdominal torque magnitude (N m)');
    legend(ax,{'Unassisted','With NMPC'},'Location','northwest');
    format_axis(ax,sequence,false);
    ax = nexttile(layout,2+source);
    plot(ax,sequence.time,nominal(:,26),':','Color',[.55 .55 .55],'LineWidth',1.3); hold(ax,'on');
    plot(ax,sequence.time,actual(:,26),'Color',[.85 .45 .1],'LineWidth',1.4);
    ylabel(ax,'Soft boundary torque (N m)'); xlabel(ax,'Continuous time (s)');
    format_axis(ax,sequence,false);
end
title(layout,[subject.id ' | abdominal and boundary torque'],'Interpreter','none','FontSize',15);
subtitle(layout,'Abdominal torque acts in flexion; the signed joint contribution is minus this magnitude.');
save_figure(fig,folder,[subject.shortName '_muscle_components'],showFigures);
end

function format_axis(ax, sequence, labelStages)
grid(ax,'on'); box(ax,'on'); ax.FontSize = 10;
xlim(ax,[sequence.time(1),sequence.time(end)]);
if labelStages
    limits = ylim(ax);
    ylim(ax,[limits(1),limits(2)+.16*diff(limits)]);
    for stage = 1:3
        s = sequence.segments{stage};
        center = mean(sequence.time([s.startIndex,s.endIndex]));
        label = sprintf('%g kg / code %g',s.loadKg,s.speedCode);
        if s.speedCode < 0 && isfield(s,'cadenceBPM')
            label = sprintf('%g kg / %.2f BPM',s.loadKg,s.cadenceBPM);
        end
        text(ax,center,limits(2)+.04*diff(limits),label, ...
            'HorizontalAlignment','center','FontSize',9,'Clipping','on');
    end
end
end

function save_figure(fig, folder, name, showFigures)
format_demo_figure(fig);
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',160);
savefig(fig,fullfile(folder,[name '.fig']));
if ~showFigures, close(fig); end
end
