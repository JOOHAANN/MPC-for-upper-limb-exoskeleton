function plot_pid_nmpc_comparison(results,folder,showFigures)
%PLOT_PID_NMPC_COMPARISON English publication figures from unchanged traces.
% Nine-condition panels use 0.4 s before and 0.8 s after each recorded action.
% A separate deviation figure resolves small angle changes without shifting
% or rescaling the absolute trajectories. Timing and angle are distinct.
sources = {'Experimental activation','MAP activation'};
suffixes = {'experimental','MAP'};
protocolName = protocol_label(results.mode);
colors.nominal = [.80 .15 .15];
colors.pid = [.88 .42 .08];
colors.nmpc = [0 .36 .70];
colors.pidAssist = [.15 .55 .20];
colors.nmpcAssist = [.55 .20 .65];
for source=1:2
    if isempty(results.groups{1}.runs{source}), continue, end
    summary = results.summary(results.summary.Source==string(sources{source}),:);
    if strcmp(results.mode,'threshold_tracking')
        effortNote='Same source-specific NMPC tuning and trajectories as RQ1; shared actuator limits';
    else
        effortNote = sprintf('Total squared-torque effort mismatch: %+.2f%% across nine conditions', ...
            summary.EffortMismatchPercent);
    end
    for kind=1:4
        labels = {'RQ1 torque comparison','RQ2 L5-S1 angle','RQ2 angle deviation','RQ2 demand-assistance timing'};
        showThis = showFigures && ((strcmp(results.mode,'actuator_matched') && kind==1) ...
            || (strcmp(results.mode,'effort_matched') && any(kind==[2 3])) ...
            || (strcmp(results.mode,'threshold_tracking') && any(kind==[1 2 3])));
        fig = new_figure([labels{kind} ' | ' results.mode ' | ' sources{source}],showThis,[18.5 16.5]);
        layout=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
        legendHandles=[];
        errorLimit=angle_error_limit(results,source);
        for loadIndex=1:3
            loadKg=[5 10 15]; loadKg=loadKg(loadIndex);
            for speedIndex=1:3
                speedCode=[35 45 55]; speedCode=speedCode(speedIndex);
                [sequence,pair,segment]=find_case(results,source,loadKg,speedCode);
                first=max(1,segment.startIndex-40);
                last=min(numel(sequence.time),segment.endIndex+80);
                index=first:last;
                time=sequence.time(index)-sequence.time(segment.startIndex);
                threshold=pair.threshold(segment.startIndex);
                ax=nexttile(layout); hold(ax,'on');
                switch kind
                    case 1
                        handles=draw_torques(ax,time,pair,index,threshold,colors);
                        if loadIndex==1 && speedIndex==1, legendHandles=handles; end
                        if speedIndex==1, ylabel(ax,{'Muscle / assistance','torque (N m)'}); end
                    case 2
                        handles=draw_angles(ax,time,pair,index,colors);
                        if loadIndex==1 && speedIndex==1, legendHandles=handles; end
                        if speedIndex==1, ylabel(ax,'L5/S1 angle (rad)'); end
                    case 3
                        handles=draw_deviation(ax,time,pair,index,colors);
                        ylim(ax,[-errorLimit errorLimit]);
                        if loadIndex==1 && speedIndex==1, legendHandles=handles; end
                        if speedIndex==1, ylabel(ax,'Angle change (deg)'); end
                        p=metric_row(results,source,loadKg,speedCode,"PID");
                        n=metric_row(results,source,loadKg,speedCode,"NMPC");
                        subtitle(ax,sprintf('RMSE: PID %.2f; NMPC %.2f deg', ...
                            p.ThetaDeviationRMSEDeg,n.ThetaDeviationRMSEDeg), ...
                            'FontSize',7,'FontWeight','normal','FontName','Times New Roman');
                    case 4
                        demand=max(0,pair.baseline(index,4)-threshold);
                        h1=plot(ax,time,demand/max(max(demand),eps),'k--','LineWidth',1.05);
                        u=pair.pid.actual(index,25);
                        h2=plot(ax,time,u/max(max(u),eps),'Color',colors.pid,'LineWidth',1.15);
                        u=pair.nmpc.actual(index,25);
                        h3=plot(ax,time,u/max(max(u),eps),'Color',colors.nmpc,'LineWidth',1.15);
                        handles=[h1 h2 h3];
                        if loadIndex==1 && speedIndex==1, legendHandles=handles; end
                        ylim(ax,[-.03 1.27]);
                        if speedIndex==1, ylabel(ax,'Normalized torque'); end
                        p=metric_row(results,source,loadKg,speedCode,"PID");
                        n=metric_row(results,source,loadKg,speedCode,"NMPC");
                        subtitle(ax,sprintf('Lag: PID %+.0f; NMPC %+.0f ms', ...
                            p.DemandLagMs,n.DemandLagMs),'FontSize',6.8, ...
                            'FontWeight','normal','FontName','Times New Roman');
                end
                title(ax,sprintf('%g kg | %g BPM',loadKg,speedCode/5),'FontWeight','normal');
                if loadIndex==3, xlabel(ax,'Time from action onset (s)'); end
                format_axis(ax,time);
            end
        end
        switch kind
            case 1
                names={'Unassisted muscle torque','PID muscle torque','NMPC muscle torque', ...
                    'PID assistance','NMPC assistance','Muscle threshold','Actuator limit', ...
                    'Upper-clipped reference'};
                stem='RQ1_torque_all_nine';
                heading='RQ1: muscle-torque threshold exceedance';
                note='Both actuators: ideal instantaneous torque, 0-80 N m, 20 ms updates';
                if isfield(pair.options,'instantRelease') && pair.options.instantRelease
                    note='Both actuators: 0-80 N m; increase <=200 N m/s; instantaneous withdrawal';
                elseif ~isfield(pair.options,'instantActuator') || ~pair.options.instantActuator
                    note='Both actuators: 0-80 N m, 200 N m/s, 20 ms updates';
                end
            case 2
                names={'Unassisted BSE','With PID','With NMPC'};
                stem='RQ2_theta_all_nine'; heading='RQ2: effect of assistance on L5/S1 motion';
                note=effortNote;
            case 3
                names={'PID minus unassisted','NMPC minus unassisted'};
                stem='RQ2_theta_deviation_all_nine'; heading='RQ2: angle deviation from unassisted motion';
                note=effortNote;
            case 4
                names={'Unassisted excess-torque demand','Delivered PID assistance','Delivered NMPC assistance'};
                stem='RQ2_timing_all_nine'; heading='RQ2: demand and delivered assistance timing';
                note='Each profile normalized by its own peak; no temporal shifts';
        end
        lg=legend(legendHandles,names,'Orientation','horizontal','Box','off', ...
            'NumColumns',min(3,numel(names)),'FontSize',7.4);
        lg.Layout.Tile='south';
        title(layout,{[heading ' | ' sources{source}],protocolName}, ...
            'FontSize',10,'FontWeight','normal');
        subtitle(layout,note,'FontSize',8);
        save_figure(fig,folder,[stem '_' suffixes{source}],showThis);
    end
    continuous_example(results,source,colors,folder,suffixes{source},sources{source},showFigures);
end
end

function handles=draw_torques(ax,time,pair,index,threshold,c)
h1=plot(ax,time,pair.baseline(index,4),'--','Color',c.nominal,'LineWidth',1.05);
h2=plot(ax,time,pair.pid.actual(index,4),'Color',c.pid,'LineWidth',1.1);
h3=plot(ax,time,pair.nmpc.actual(index,4),'Color',c.nmpc,'LineWidth',1.1);
h4=plot(ax,time,pair.pid.actual(index,25),'-.','Color',c.pidAssist,'LineWidth',.95);
h5=plot(ax,time,pair.nmpc.actual(index,25),'-.','Color',c.nmpcAssist,'LineWidth',.95);
h6=yline(ax,threshold,'k--','LineWidth',.95);
h7=yline(ax,pair.options.capNm,'--','Color',[.45 .45 .45],'LineWidth',.95);
target=min(pair.baseline(index,4),pair.threshold(index));
h8=plot(ax,time,target,':','Color',[.1 .45 .45],'LineWidth',1.1);
handles=[h1 h2 h3 h4 h5 h6 h7 h8];
upper=max([pair.baseline(index,4);pair.pid.actual(index,4);pair.nmpc.actual(index,4); ...
    pair.pid.actual(index,25);pair.nmpc.actual(index,25);threshold;pair.options.capNm]);
lower=min([0;pair.baseline(index,4);pair.pid.actual(index,4);pair.nmpc.actual(index,4)]);
ylim(ax,[lower-.025*(upper-lower),upper+.08*(upper-lower)]);
end

function handles=draw_angles(ax,time,pair,index,c)
h1=plot(ax,time,pair.baseline(index,1),'--','Color',c.nominal,'LineWidth',1.15);
h2=plot(ax,time,pair.pid.actual(index,1),'Color',c.pid,'LineWidth',1.1);
h3=plot(ax,time,pair.nmpc.actual(index,1),'Color',c.nmpc,'LineWidth',1.1);
handles=[h1 h2 h3];
end

function handles=draw_deviation(ax,time,pair,index,c)
yline(ax,0,'k:','LineWidth',.7,'HandleVisibility','off');
h1=plot(ax,time,rad2deg(pair.pid.actual(index,1)-pair.baseline(index,1)), ...
    'Color',c.pid,'LineWidth',1.1);
h2=plot(ax,time,rad2deg(pair.nmpc.actual(index,1)-pair.baseline(index,1)), ...
    'Color',c.nmpc,'LineWidth',1.1);
handles=[h1 h2];
end

function continuous_example(results,source,c,folder,suffix,name,showFigures)
sequence=results.groups{1}.sequence;
pair=results.groups{1}.runs{source};
time=sequence.time; index=1:numel(time);
protocolName=protocol_label(results.mode);
plot_RQ1_continuous(results,source,c,folder,suffix,name,showFigures);

fig=new_figure(['RQ2 continuous example | ' results.mode ' | ' name],showFigures,[18.5 12]);
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout); hold(ax,'on');
h=draw_angles(ax,time,pair,index,c);
format_axis(ax,time); ylabel(ax,'L5/S1 angle (rad)');
legend(ax,h,{'Unassisted BSE','With PID','With NMPC'},'Location','southoutside', ...
    'Orientation','horizontal','Box','off','FontSize',8);
label_actions(ax,sequence);
ax=nexttile(layout); hold(ax,'on');
h=draw_deviation(ax,time,pair,index,c);
format_axis(ax,time); ylabel(ax,'Angle change (deg)'); xlabel(ax,'Time (s)');
legend(ax,h,{'PID minus unassisted','NMPC minus unassisted'},'Location','southoutside', ...
    'Orientation','horizontal','Box','off','FontSize',8);
if source==2
    mainTitle=title(layout,'Motion and deviation', ...
        'FontSize',10,'FontWeight','normal');
    titleScale=2;
else
    mainTitle=title(layout,{['RQ2: motion and deviation | ' name],protocolName}, ...
        'FontSize',10,'FontWeight','normal');
    titleScale=1;
end
save_figure(fig,folder,['RQ2_continuous_example_' suffix],showFigures,mainTitle,titleScale);
end

function label=protocol_label(mode)
if strcmp(mode,'threshold_tracking')
    label='Clipped-reference NMPC | Shared actuator limits';
elseif strcmp(mode,'actuator_matched')
    label='Original NMPC | Matched actuator limits';
else
    label='Effort-matched NMPC | Soft muscle-torque threshold';
end
end

function label_actions(ax,sequence)
limits=ylim(ax); span=diff(limits);
ylim(ax,[limits(1) limits(2)+.16*span]);
for stage=1:3
    segment=sequence.segments{stage};
    center=mean(sequence.time([segment.startIndex segment.endIndex]));
    text(ax,center,limits(2)+.05*span,sprintf('%g kg | %g BPM', ...
        segment.loadKg,segment.speedCode/5),'HorizontalAlignment','center','FontSize',8);
end
end

function [sequence,pair,segment]=find_case(results,source,loadKg,speedCode)
for group=1:3
    sequence=results.groups{group}.sequence;
    for stage=1:3
        segment=sequence.segments{stage};
        if segment.loadKg==loadKg && segment.speedCode==speedCode
            pair=results.groups{group}.runs{source}; return
        end
    end
end
error('Comparison:Case','The requested load/speed combination is missing.');
end

function row=metric_row(results,source,loadKg,speedCode,controller)
t=results.metrics;
row=t(t.SourceIndex==source & t.LoadKg==loadKg & t.SpeedCode==speedCode & t.Controller==controller,:);
end

function limit=angle_error_limit(results,source)
value=0;
for group=1:3
    pair=results.groups{group}.runs{source};
    value=max(value,max(abs(pair.pid.actual(:,1)-pair.baseline(:,1))));
    value=max(value,max(abs(pair.nmpc.actual(:,1)-pair.baseline(:,1))));
end
limit=max(.25,1.15*rad2deg(value));
end

function fig=new_figure(name,showFigures,sizeCm)
visible='off'; if showFigures, visible='on'; end
fig=figure('Name',name,'Color','w','Visible',visible,'Units','centimeters', ...
    'Position',[2 2 sizeCm],'Renderer','painters');
end

function format_axis(ax,time)
set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.65,'Box','on', ...
    'TickDir','in','XGrid','on','YGrid','on','GridAlpha',.10);
xlim(ax,[time(1) time(end)]);
end

function save_figure(fig,folder,stem,showFigures,mainTitle,titleScale)
format_demo_figure(fig);
% Scale these titles after the shared formatter so its default cannot override them.
if nargin>=6
    mainTitle.FontSize=mainTitle.FontSize*titleScale;
    drawnow;
end
exportgraphics(fig,fullfile(folder,[stem '.png']),'Resolution',400);
exportgraphics(fig,fullfile(folder,[stem '.pdf']),'ContentType','vector');
savefig(fig,fullfile(folder,[stem '.fig']));
% Publication copies use actual 8-point text at one printed column wide.
if startsWith(stem,'RQ1_continuous_example_')
    export_paper_figure(fig,fullfile(folder,[stem '_paper']),2.45);
elseif startsWith(stem,'RQ2_continuous_example_')
    export_paper_figure(fig,fullfile(folder,[stem '_paper']),3.0);
end
if ~showFigures, close(fig); end
end
