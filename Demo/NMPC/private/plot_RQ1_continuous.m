function plot_RQ1_continuous(results,source,c,folder,suffix,name,showFigures)
%PLOT_RQ1_CONTINUOUS Resolve muscle loading and assistance in separate panels.
% All traces come from the saved continuous closed-loop simulation. Only
% settled standing samples are omitted from the displayed time axis.
q=results.groups{1}.sequence;pair=results.groups{1}.runs{source};
c.pidAssist=[.08 .48 .22];c.nmpcAssist=[.55 .15 .65];
[index,time,cuts,segments,omitted]=shorten_standing(q,pair);
writetable(omitted,fullfile(folder,['RQ1_display_omissions_' suffix '.csv']));
mapping=table(q.time(index),time,'VariableNames',{'OriginalTimeS','DisplayedTimeS'});
writetable(mapping,fullfile(folder,['RQ1_display_time_' suffix '.csv']));
visible='off';if showFigures,visible='on';end
fig=figure('Name',['RQ1 continuous example | ' name],'Color','w','Visible',visible, ...
    'Units','centimeters','Position',[2 2 18.5 13],'Renderer','painters', ...
    'GraphicsSmoothing','on');
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
plotTime=time;plotTime(cuts)=NaN; % Do not join lines across omitted intervals.
target=min(pair.baseline(index,4),pair.threshold(index));
threshold=nan(size(time));
for k=1:numel(segments)
    ids=segments{k}.startIndex:segments{k}.endIndex;
    threshold(ids)=pair.threshold(index(ids));
end

ax=nexttile(layout);hold(ax,'on');
h1=plot(ax,plotTime,pair.baseline(index,4),'--','Color',c.nominal);
h5=plot(ax,plotTime,target,':','Color',[.1 .45 .45]);
h2=plot(ax,plotTime,pair.pid.actual(index,4),'-','Color',c.pid);
h3=plot(ax,plotTime,pair.nmpc.actual(index,4),'-','Color',c.nmpc);
h4=plot(ax,plotTime,threshold,'k--');
style_lines([h1 h2 h3 h4 h5],[1.0 1.05 1.2 .95 .7]);
format_axis(ax,time);ylabel(ax,{'Muscle torque','(N m)'});ax.XTickLabel=[];
upper=max([pair.baseline(index,4);pair.pid.actual(index,4);pair.nmpc.actual(index,4)]);
ylim(ax,[-8 upper*1.15]);
limits=ylim(ax);
for k=1:numel(segments)
    seg=segments{k};center=mean(time([seg.startIndex seg.endIndex]));
    text(ax,center,limits(2)-.035*diff(limits), ...
        sprintf('%g kg | %g BPM',seg.loadKg,seg.speedCode/5), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8.2);
end
legend(ax,[h1 h2 h3 h4 h5],{'Unassisted muscle','PID muscle','NMPC muscle', ...
    'Muscle threshold','Upper-clipped reference'},'Location','southoutside', ...
    'NumColumns',2,'Box','off','FontSize',8.2);

ax=nexttile(layout);hold(ax,'on');
% Solid strokes resolve the fast assistance changes without smoothing them.
h1=plot(ax,plotTime,pair.pid.actual(index,25),'-','Color',c.pidAssist);
h2=plot(ax,plotTime,pair.nmpc.actual(index,25),'-','Color',c.nmpcAssist);
h3=yline(ax,pair.options.capNm,'--','Color',[.4 .4 .4],'LineWidth',.9);
style_lines([h1 h2],[1.05 1.15]);
format_axis(ax,time);ylabel(ax,{'Assistance','(N m)'});
ylim(ax,[-2 pair.options.capNm*1.06]);yticks(ax,[0 40 80]);
xlabel(ax,'Displayed time (s)');
legend(ax,[h1 h2 h3],{'PID assistance','NMPC assistance','80 N m limit'}, ...
    'Location','southoutside','Orientation','horizontal','NumColumns',3, ...
    'Box','off','FontSize',8.2);
title(layout,'Continuous lifting comparison between PID and NMPC', ...
    'FontName','Times New Roman','FontWeight','normal','FontSize',10);

format_demo_figure(fig);
stem=fullfile(folder,['RQ1_continuous_example_' suffix]);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
exportgraphics(fig,[stem '.png'],'Resolution',900);savefig(fig,[stem '.fig']);
export_paper_figure(fig,[stem '_paper'],2.85,900);
if ~showFigures,close(fig);end
fprintf('RQ1 %s: %.2f -> %.2f s displayed; movement and recovery traces retained.\n', ...
    suffix,q.time(end),time(end));
end

function style_lines(lines,widths)
for k=1:numel(lines)
    lines(k).LineWidth=widths(k);
    lines(k).UserData=struct('PaperLineWidth',widths(k));
end
end

function format_axis(ax,time)
set(ax,'FontName','Times New Roman','FontSize',8.2,'LineWidth',.6, ...
    'Box','on','TickDir','in','XGrid','on','YGrid','on','GridAlpha',.10);
xlim(ax,[time(1) time(end)]);
end

function [index,time,cuts,segments,omitted]=shorten_standing(q,pair)
% Keep every movement and unsettled recovery sample. The tolerances below
% identify stationary waits, not admissible errors for either controller.
assert(max(abs(diff(q.time)-q.dt))<1e-10);
assistance=max(abs([pair.pid.actual(:,25) pair.nmpc.actual(:,25)]),[],2);
angleError=max(abs([pair.pid.actual(:,1)-q.standingAngle ...
    pair.nmpc.actual(:,1)-q.standingAngle pair.baseline(:,1)-q.standingAngle]),[],2);
velocity=max(abs([pair.pid.actual(:,2) pair.nmpc.actual(:,2) pair.baseline(:,2)]),[],2);
torqueError=max(abs([pair.pid.actual(:,4)-pair.baseline(:,4) ...
    pair.nmpc.actual(:,4)-pair.baseline(:,4)]),[],2);
settled=q.waiting>.5 & assistance<.05 & angleError<deg2rad(.15) ...
    & velocity<deg2rad(.5) & torqueError<.4;
edges=diff([false;settled;false]);first=find(edges==1);last=find(edges==-1)-1;
keep=true(size(q.time));rows=struct([]);edgeSamples=round(.5/q.dt);
for k=1:numel(first)
    if last(k)-first(k)+1<=round(2/q.dt),continue,end
    ids=(first(k)+edgeSamples:last(k)-edgeSamples)';
    keep(ids)=false;
    row=struct('OriginalFirstS',q.time(ids(1)),'OriginalLastS',q.time(ids(end)), ...
        'RemovedDurationS',numel(ids)*q.dt,'MaximumOmittedAssistanceNm',max(assistance(ids)), ...
        'MaximumOmittedAngleErrorDeg',rad2deg(max(angleError(ids))), ...
        'MaximumOmittedMuscleDifferenceNm',max(torqueError(ids)));
    rows=[rows;row]; %#ok<AGROW>
end
assert(all(keep(~settled)) && all(keep(q.waiting<.5)));
index=find(keep);time=(0:numel(index)-1)'*q.dt;
cuts=find(diff(index)>1)+1;segments=q.segments;
for k=1:numel(segments)
    seg=segments{k};original=seg.startIndex:seg.endIndex;
    assert(all(keep(original)),'A movement sample was removed from the display.');
    seg.startIndex=find(index==seg.startIndex,1);seg.endIndex=find(index==seg.endIndex,1);
    assert(seg.endIndex-seg.startIndex+1==numel(original));
    segments{k}=seg;
end
omitted=struct2table(rows);
end
