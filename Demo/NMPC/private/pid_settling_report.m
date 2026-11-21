function [events,summaryTable]=pid_settling_report(results,folder,showFigures)
%PID_SETTLING_REPORT Record post-peak threshold-band entry for each action.
% This is an operational in-action metric, not an LTI step-response settling
% time. Demand windows come from the unassisted replay ONLY for evaluation;
% no future demand or unassisted torque is supplied to the reactive PID.
% A window must last >=0.30 s. After its largest positive overshoot, entry
% requires >=0.15 s continuously within +/-max(5 N m, 5% of the threshold).
% A later final settling time is reported if the trace re-enters and stays
% inside for the rest of the recorded demand interval, with >=0.15 s hold.
names={'Experimental activation','MAP activation'};
controllers={'PID','NMPC'};rows=struct([]);
for group=1:3
    q=results.groups{group}.sequence;
    for source=1:2
        pair=results.groups{group}.runs{source};
        for stage=1:3
            seg=q.segments{stage};ix=(seg.startIndex:seg.endIndex)';
            threshold=pair.threshold(ix);demand=pair.baseline(ix,4)>threshold;
            edge=diff([false;demand;false]);starts=find(edge==1);ends=find(edge==-1)-1;
            for event=1:numel(starts)
                first=starts(event);last=ends(event);
                if (last-first+1)*q.dt<.30-1e-10,continue,end
                indices=ix(first:last);limit=threshold(first);band=max(5,.05*limit);
                for method=1:2
                    run=pair.(lower(controllers{method}));error=run.actual(indices,4)-limit;
                    [peak,peakIndex]=max(error);inside=abs(error)<=band;
                    % Suppress sample-scale noise only for counting excursions;
                    % all tracking and band-entry values use the raw torque.
                    smoothError=movmean(error,max(1,round(.05/q.dt)));
                    high=smoothError>band;low=smoothError<-band;
                    overshoots=sum(diff([false;high;false])==1);
                    undershoots=sum(diff([false;low;false])==1);
                    % Include both endpoints so the observed time span is
                    % at least 0.15 s, rather than just 15 samples at 100 Hz.
                    holdSamples=max(2,ceil(.15/q.dt-1e-10)+1);
                    if numel(inside)>=holdSamples
                        windows=conv(double(inside),ones(holdSamples,1),'valid');
                        entry=find(windows>=holdSamples & (1:numel(windows))'>=peakIndex,1);
                    else
                        entry=[];
                    end
                    firstEntry=NaN;confirmed=NaN;stayed=false;
                    if ~isempty(entry)
                        firstEntry=(entry-1)*q.dt;
                        confirmed=(entry+holdSamples-2)*q.dt;
                        stayed=all(inside(entry:end));
                    end
                    lastOutside=find(~inside(peakIndex:end),1,'last');
                    if isempty(lastOutside),finalIndex=peakIndex;
                    else,finalIndex=peakIndex+lastOutside;end
                    settled=numel(inside)-finalIndex+1>=holdSamples;
                    finalTime=NaN;if settled,finalTime=(finalIndex-1)*q.dt;end
                    row=struct('Source',string(names{source}),'SourceIndex',source, ...
                        'Controller',string(controllers{method}),'Group',group,'Stage',stage, ...
                        'LoadKg',seg.loadKg,'SpeedCode',seg.speedCode,'Event',event, ...
                        'ThresholdNm',limit,'BandHalfWidthNm',band,'DemandStartIndex',indices(1), ...
                        'DemandStartTimeSeconds',q.time(indices(1)), ...
                        'DemandDurationSeconds',numel(indices)*q.dt, ...
                        'PeakOvershootNm',max(peak,0),'PeakTimeSeconds',(peakIndex-1)*q.dt, ...
                        'OvershootExcursions',overshoots,'UndershootExcursions',undershoots, ...
                        'SustainedEntry',~isempty(entry),'FirstEntrySeconds',firstEntry, ...
                        'ConfirmationTimeSeconds',confirmed,'StayedUntilDemandEnd',stayed, ...
                        'SettledByDemandEnd',settled,'FinalSettlingTimeSeconds',finalTime, ...
                        'BandFractionAfterPeak',mean(inside(peakIndex:end)));
                    rows=[rows;row]; %#ok<AGROW>
                end
            end
        end
    end
end
assert(~isempty(rows),'No eligible demand intervals were found.');
events=struct2table(rows);summary=struct([]);
for source=1:2
    for method=1:2
        x=events(events.SourceIndex==source & events.Controller==string(controllers{method}),:);
        entered=x.FirstEntrySeconds(x.SustainedEntry);
        final=x.FinalSettlingTimeSeconds(x.SettledByDemandEnd);
        row=struct('Source',string(names{source}),'Controller',string(controllers{method}), ...
            'EligibleEvents',height(x),'SustainedEntries',sum(x.SustainedEntry), ...
            'StayedUntilDemandEnd',sum(x.StayedUntilDemandEnd), ...
            'SettledByDemandEnd',sum(x.SettledByDemandEnd), ...
            'AtMostTwoOvershootsAndSettled',sum(x.OvershootExcursions<=2 & x.SettledByDemandEnd), ...
            'MeanEntrySecondsOfEntries',mean(entered,'omitnan'), ...
            'MeanFinalSettlingSecondsOfSettled',mean(final,'omitnan'), ...
            'MeanBandFractionAfterPeak',mean(x.BandFractionAfterPeak));
        summary=[summary;row]; %#ok<AGROW>
    end
end
summaryTable=struct2table(summary);
writetable(events,fullfile(folder,'PID_NMPC_band_entry_events.csv'));
writetable(summaryTable,fullfile(folder,'PID_NMPC_band_entry_summary.csv'));

%% Show the three continuous-sequence examples, including the tolerance band
visible='off';if showFigures,visible='on';end
fig=figure('Name','PID and NMPC band entry','Color','w','Visible',visible, ...
    'Position',[40 40 1400 1050]);
layout=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
q=results.groups{1}.sequence;colors={[.84 .35 .08],[0 .36 .70]};handles=[];
shortNames={'Experimental','MAP'};
for stage=1:3
    seg=q.segments{stage};ix=(max(1,seg.startIndex-40):min(numel(q.time),seg.endIndex+80))';
    t=q.time(ix)-q.time(seg.startIndex);
    for source=1:2
        pair=results.groups{1}.runs{source};ax=nexttile(layout);hold(ax,'on');
        limit=pair.threshold(seg.startIndex);band=max(5,.05*limit);
        x0=0;x1=q.time(seg.endIndex)-q.time(seg.startIndex);
        patch(ax,[x0 x1 x1 x0],[limit-band limit-band limit+band limit+band], ...
            [.7 .85 .7],'FaceAlpha',.24,'EdgeColor','none','HandleVisibility','off');
        h1=plot(ax,t,pair.baseline(ix,4),':','Color',[.5 .5 .5],'LineWidth',1.0);
        h2=plot(ax,t,pair.pid.actual(ix,4),'Color',colors{1},'LineWidth',1.15);
        h3=plot(ax,t,pair.nmpc.actual(ix,4),'Color',colors{2},'LineWidth',1.15);
        h4=plot(ax,[x0 x1],[limit limit],'k--','LineWidth',1.0);
        if isempty(handles),handles=[h1 h2 h3 h4];end
        for method=1:2
            x=events(events.Group==1 & events.Stage==stage & events.SourceIndex==source ...
                & events.Controller==string(controllers{method}) & events.SustainedEntry,:);
            for k=1:height(x)
                at=x.DemandStartIndex(k)+round(x.FirstEntrySeconds(k)/q.dt);
                run=pair.(lower(controllers{method}));
                plot(ax,q.time(at)-q.time(seg.startIndex),run.actual(at,4),'o', ...
                    'Color',colors{method},'MarkerFaceColor',colors{method}, ...
                    'MarkerSize',4,'HandleVisibility','off');
            end
        end
        a=events(events.Group==1 & events.Stage==stage & events.SourceIndex==source & events.Controller=="PID",:);
        b=events(events.Group==1 & events.Stage==stage & events.SourceIndex==source & events.Controller=="NMPC",:);
        title(ax,sprintf('%s | %g kg / %g',shortNames{source},seg.loadKg,seg.speedCode), ...
            'FontWeight','normal','FontSize',8.5);
        subtitle(ax,sprintf('PID E %d/%d S %d | NMPC E %d/%d S %d', ...
            sum(a.SustainedEntry),height(a),sum(a.SettledByDemandEnd), ...
            sum(b.SustainedEntry),height(b),sum(b.SettledByDemandEnd)));
        ylabel(ax,'Back torque (N m)');if stage==3,xlabel(ax,'Time from action onset (s)');end
        grid(ax,'on');box(ax,'on');xlim(ax,[t(1) t(end)]);set(ax,'FontName','Times New Roman','FontSize',9);
    end
end
lg=legend(handles,{'Unassisted','PID','NMPC','Threshold'},'Orientation','horizontal','Box','off');
lg.Layout.Tile='south';
title(layout,'Subject 528 | Post-peak threshold-band entry','FontWeight','normal');
subtitle(layout,{'Band: +/- max(5 N m, 5% threshold). E: post-peak entry held >=0.15 s; S: final stay until demand ends.', ...
    'Counts use demand intervals >=0.30 s. This is an in-action metric, not a fixed-step settling time.'});
format_demo_figure(fig);
exportgraphics(fig,fullfile(folder,'PID_NMPC_band_entry_528.png'),'Resolution',180);
exportgraphics(fig,fullfile(folder,'PID_NMPC_band_entry_528.pdf'),'ContentType','vector');
savefig(fig,fullfile(folder,'PID_NMPC_band_entry_528.fig'));
if ~showFigures,close(fig);end
end
