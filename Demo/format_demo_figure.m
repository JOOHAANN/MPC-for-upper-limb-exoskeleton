function format_demo_figure(fig)
%FORMAT_DEMO_FIGURE Apply shared typography and line widths to demo plots.
% Increase the previous figure text by 20%. Enlarge tick numbers a further
% 50% so that the numeric scales are easier to read.
% All existing plot and border line widths are scaled by 1.5, once.
bodySize=13.6*1.2;
titleSize=17*1.2;
tickSize=13.6*1.2*1.5;
fontName='Times New Roman';

% Give multi-panel layouts room for the larger labels and outside legends.
drawnow;
axesList=findall(fig,'Type','axes');
count=numel(axesList);
limits=cell(count,3);
for k=1:count
    ax=axesList(k);
    limits(k,:)={ax.XLim,ax.YLim,ax.ZLim};
    if ~isappdata(ax,'DemoOriginalLimitModes')
        setappdata(ax,'DemoOriginalLimitModes',{ax.XLimMode,ax.YLimMode,ax.ZLimMode});
    end
end
% Inside legends would cover substantially more data after enlarging text.
% Place them above their own axes; preserve tiled-layout legends below plots.
legends=findall(fig,'Type','legend');
inside={'northwest','northeast','southwest','southeast','best'};
for k=1:numel(legends)
    item=legends(k);
    if any(strcmp(item.Location,inside))
        item.Location='northoutside';
        item.Orientation='horizontal';
        labels=item.String;
        if ischar(labels),entryCount=size(labels,1);else,entryCount=numel(labels);end
        item.NumColumns=min(3,max(1,entryCount));
    end
end
layouts=findall(fig,'Type','tiledlayout');
rows=0;
for k=1:numel(layouts)
    rows=max(rows,layouts(k).GridSize(1));
end
if rows>=4 || count>=10
    minimumSize=[1600 1450];
elseif rows==3 || count>=8
    minimumSize=[1500 1150];
elseif count>=5
    minimumSize=[1350 1000];
elseif count>=2
    minimumSize=[1200 850];
else
    minimumSize=[1200 650];
end
originalUnits=fig.Units;fig.Units='pixels';
position=fig.Position;position(3:4)=max(position(3:4),minimumSize);
fig.Position=position;fig.Units=originalUnits;

% Record all widths first, because changing an axes border can update its
% rulers. Idempotence keeps a second formatting call from thickening again.
if ~isappdata(fig,'DemoLineWidthsScaled150')
    lines=findall(fig,'-property','LineWidth');
    widths=zeros(size(lines));
    for k=1:numel(lines),widths(k)=lines(k).LineWidth;end
    for k=1:numel(lines),lines(k).LineWidth=1.5*widths(k);end
    setappdata(fig,'DemoLineWidthsScaled150',true);
end

% Set size directly on each graphics object so explicit old sizes do not
% override the shared style. Preserve line data, colors and plot limits.
objects=findall(fig,'-property','FontSize');
for k=1:numel(objects)
    item=objects(k);
    if isprop(item,'FontUnits'),item.FontUnits='points';end
    if isprop(item,'FontSizeMode'),item.FontSizeMode='manual';end
    item.FontSize=bodySize;
    if isprop(item,'FontName'),item.FontName=fontName;end
end

% MATLAB stores axes and tiled-layout titles separately from regular text.
for k=1:numel(axesList)
    ax=axesList(k);
    % MATLAB restores axes subtitles to auto sizing when a FIG is opened.
    % Keep their auto size equal to the body font despite larger tick labels.
    ax.FontSize=tickSize;
    ax.XAxis.FontSize=tickSize;ax.YAxis.FontSize=tickSize;ax.ZAxis.FontSize=tickSize;
    ax.TitleFontSizeMultiplier=bodySize/tickSize;
    ax.Title.FontUnits='points';ax.Title.FontSize=titleSize;ax.Title.FontName=fontName;
    ax.Subtitle.FontUnits='points';ax.Subtitle.FontSize=bodySize;ax.Subtitle.FontName=fontName;
    labels=[ax.XLabel ax.YLabel ax.ZLabel];
    for j=1:numel(labels)
        labels(j).FontUnits='points';labels(j).FontSize=bodySize;labels(j).FontName=fontName;
    end
    % Short panels need long vertical labels on two lines to avoid adjacent
    % titles and labels, especially after the font-size increase.
    label=ax.YLabel.String;
    wrapLabel=ischar(label) && ((rows>=4 && numel(label)>15) || ...
        (rows==3 && numel(label)>24));
    if wrapLabel
        spaces=find(isspace(label));
        if ~isempty(spaces)
            [~,nearest]=min(abs(spaces-(numel(label)+1)/2));
            split=spaces(nearest);
            ax.YLabel.String={strtrim(label(1:split-1)),strtrim(label(split+1:end))};
        end
    end
end
for k=1:numel(layouts)
    layouts(k).Padding='loose';
    layouts(k).Title.FontSize=titleSize;layouts(k).Title.FontName=fontName;
    layouts(k).Subtitle.FontSize=bodySize;layouts(k).Subtitle.FontName=fontName;
end
% Resizing can otherwise change automatic limits even though no data changed.
for k=1:count
    ax=axesList(k);ax.XLim=limits{k,1};ax.YLim=limits{k,2};ax.ZLim=limits{k,3};
end
drawnow;

% Short panels can otherwise stack the enlarged y tick numbers.
% Retain evenly spaced existing ticks while leaving the plotted limits intact.
if rows>=3
    pixelsPerPoint=get(groot,'ScreenPixelsPerInch')/72;
    for k=1:count
        ax=axesList(k);oldUnits=ax.Units;ax.Units='pixels';height=ax.Position(4);ax.Units=oldUnits;
        ticks=ax.YTick;
        maxLabels=max(2,floor(height/(tickSize*pixelsPerPoint*1.1))+1);
        if numel(ticks)>maxLabels
            chosen=unique(round(linspace(1,numel(ticks),maxLabels)));
            ax.YTick=ticks(chosen);
        end
    end
end
drawnow;
end
