function metrics=plot_raw_MAP_validation(showFigures)
%PLOT_RAW_MAP_VALIDATION Compare the frozen raw MAP with the original EMG.
% Six back-muscle channels x 25 conditions. No load rescaling, output
% normalization, or fitted gain/offset is applied to MAP predictions here.
if nargin<1,showFigures=true;end
folder=fileparts(mfilename('fullpath'));addpath(fileparts(fileparts(folder)));
data=load(fullfile(folder,'data','legacy_MAP_528.mat'),'map','records');
reference=load(fullfile(folder,'data','raw_MAP_reference_528.mat'));
names={'IL right','IL left','LTpT right','LTpT left','LTpL right','LTpL left'};
rows=struct([]);predictions=zeros(size(reference.reference));
for condition=1:25
    loadKg=reference.loadKg(condition);speedCode=reference.speedCode(condition);
    [~,raw]=predict_528_MAP(data.map,loadKg,speedCode,linspace(0,1,200)');
    predictions(:,:,condition)=raw;
    observed=reference.reference(:,:,condition);
    for muscle=1:6
        r=corrcoef(raw(:,muscle),observed(:,muscle));
        row=struct('LoadKg',loadKg,'SpeedBPM',speedCode/5,'Muscle',string(names{muscle}), ...
            'Correlation',r(1,2),'RMSE',sqrt(mean((raw(:,muscle)-observed(:,muscle)).^2)));
        rows=[rows;row]; %#ok<AGROW>
    end
end
metrics=struct2table(rows);assert(height(metrics)==150);
output=fullfile(folder,'output','paper');if ~isfolder(output),mkdir(output);end
writetable(metrics,fullfile(output,'raw_MAP_metrics_528.csv'));
summary=table(mean(metrics.Correlation),std(metrics.Correlation),mean(metrics.RMSE),std(metrics.RMSE), ...
    'VariableNames',{'MeanR','SDR','MeanRMSE','SDRMSE'});
writetable(summary,fullfile(output,'raw_MAP_summary_528.csv'));
save(fullfile(output,'raw_MAP_validation_528.mat'),'metrics','summary','predictions','reference');

visible='off';if showFigures,visible='on';end
fig=figure('Color','w','Visible',visible,'Name','Raw MAP validation | no amplitude calibration', ...
    'Units','inches','Position',[1 1 3.4 1.7]);
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
values={metrics.Correlation,metrics.RMSE};labels={'Correlation','RMSE'};
h=gobjects(1,3);
for k=1:2
    ax=nexttile(layout);hold(ax,'on');
    if k==1,edges=.55:.05:1;else,edges=0:.04:.4;end
    h(1)=histogram(ax,values{k},edges,'FaceColor',[.35 .36 .85],'EdgeColor',[.2 .2 .3]);
    ylim(ax,[0 1.12*max(h(1).Values)]);
    mu=mean(values{k});sd=std(values{k});
    h(2)=xline(ax,mu,'--','Color',[.8 .1 .15],'LineWidth',1);
    h(3)=xline(ax,mu-sd,'--','Color',[0 .55 .2],'LineWidth',1);
    xline(ax,mu+sd,'--','Color',[0 .55 .2],'LineWidth',1,'HandleVisibility','off');
    xlabel(ax,labels{k});ylabel(ax,'Frequency');box(ax,'on');
    title(ax,sprintf('(%c) %s','a'+k-1,labels{k}),'FontWeight','normal');
    subtitle(ax,sprintf('Mean %.2f; SD %.2f',mu,sd),'FontSize',8.2);
    if k==1,xlim(ax,[.55 1]);xticks(ax,[.6 .8 1]);else,xlim(ax,[0 .4]);xticks(ax,[0 .2 .4]);end
end
key=legend(h,{'Raw MAP','Mean','Mean +/- SD'},'Orientation','horizontal', ...
    'NumColumns',3,'Box','off','FontSize',8.2);key.Layout.Tile='south';
title(layout,'Raw MAP predictions','FontSize',9,'FontWeight','normal');
export_paper_figure(fig,fullfile(output,'528_raw_MAP_histograms'),1.7);
if ~showFigures,close(fig);end
disp(summary);
end
