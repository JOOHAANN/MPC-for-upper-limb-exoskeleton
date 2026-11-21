function export_paper_figure(fig,filename,heightInches,resolution)
%EXPORT_PAPER_FIGURE Export a single-column copy with readable 8-point text.
% The original figure and its plotted data are left unchanged.
if nargin<4,resolution=450;end
paper=copyobj(fig,groot);
paper.Visible='off';paper.Units='inches';
paper.Position=[1 1 3.4 heightInches];paper.Color='w';
paper.PaperUnits='inches';paper.PaperSize=[3.4 heightInches];
paper.PaperPosition=[0 0 3.4 heightInches];paper.PaperPositionMode='manual';
objects=findall(paper,'-property','FontSize');
for k=1:numel(objects)
    item=objects(k);
    if isprop(item,'FontUnits'),item.FontUnits='points';end
    item.FontSize=8.2;
    if isprop(item,'FontName'),item.FontName='Times New Roman';end
end
axesList=findall(paper,'Type','axes');
for k=1:numel(axesList)
    ax=axesList(k);ax.FontSize=8.2;
    ax.XAxis.FontSize=8.2;ax.YAxis.FontSize=8.2;
    ax.LabelFontSizeMultiplier=1;ax.TitleFontSizeMultiplier=1;
    ax.Title.FontSize=9;ax.Subtitle.FontSize=8.2;
    ax.XLabel.FontSize=8.2;ax.YLabel.FontSize=8.2;
    if isequal(ax.YLabel.String,'Net muscle torque (N m)')
        ax.YLabel.String={'Net muscle torque','(N m)'};
    elseif isequal(ax.YLabel.String,'Trunk angle (rad)')
        ax.YLabel.String={'Trunk angle','(rad)'};
    end
    ax.LineWidth=.6;
    curves=findall(ax,'Type','line');values=[];
    for j=1:numel(curves),values=[values;curves(j).YData(:)];end %#ok<AGROW>
    values=values(isfinite(values));
    if ~isempty(values)
        span=max(max(values)-min(values),eps);limits=ax.YLim;
        ax.YLim=[min(limits(1),min(values)-.03*span) ...
            max(limits(2),max(values)+.05*span)];
    end
    % Convert archived speed labels to the physical movement speed in BPM.
    annotations=findall(ax,'Type','text');actionLabels=gobjects(0);
    for j=1:numel(annotations)
        label=annotations(j).String;
        if ischar(label) && contains(label,' kg') && (contains(label,'code') || contains(label,'BPM'))
            values=regexp(label,'^([\d.]+) kg\s*[|/]?\s*(?:speed\s*)?code\s*([\d.]+)$','tokens','once');
            if ~isempty(values)
                label=sprintf('%s kg | %g BPM',values{1},str2double(values{2})/5);
            end
            annotations(j).String=label;
            actionLabels(end+1)=annotations(j); %#ok<AGROW>
        end
    end
    if ~isempty(actionLabels)
        limits=ax.YLim;ax.YLim=[limits(1) limits(2)+.07*diff(limits)];
        for j=1:numel(actionLabels)
            actionLabels(j).VerticalAlignment='top';
            position=actionLabels(j).Position;
            position(2)=ax.YLim(2)-.025*diff(ax.YLim);
            actionLabels(j).Position=position;
        end
    end
end
layouts=findall(paper,'Type','tiledlayout');
for k=1:numel(layouts)
    layouts(k).Padding='compact';layouts(k).TileSpacing='compact';
    layouts(k).Title.FontSize=9;layouts(k).Subtitle.FontSize=8.2;
    layouts(k).Title.FontName='Times New Roman';
    layouts(k).Subtitle.FontName='Times New Roman';
end
legends=findall(paper,'Type','legend');
for k=1:numel(legends)
    item=legends(k);item.FontSize=8.2;
    item.ItemTokenSize=[13 8];
    if numel(item.String)>3,item.NumColumns=2;end
    labels=item.String;
    for j=1:numel(labels)
        if strcmp(labels{j},'Reference: CEINMS torque / IK angle')
            labels{j}=sprintf('Reference: CEINMS torque\n/ IK angle');
        elseif strcmp(labels{j},'BSE: experimental activation')
            labels{j}=sprintf('BSE: experimental\nactivation');
        elseif strcmp(labels{j},'Defined standing baseline')
            labels{j}=sprintf('Defined standing\nbaseline');
        end
    end
    item.String=labels;
end
% COPYOBJ does not retain a shared legend's tile assignment.
if numel(layouts)==1 && numel(legends)==1
    legends(1).Layout.Tile='south';
end
% Keep physical line widths suitable for the final 3.4-inch printed column.
lines=findall(paper,'Type','line');
for k=1:numel(lines)
    style=lines(k).UserData;
    if isstruct(style) && isfield(style,'PaperLineWidth')
        lines(k).LineWidth=style.PaperLineWidth;
    else
        lines(k).LineWidth=1.05;
    end
end
constants=findall(paper,'Type','constantline');
for k=1:numel(constants),constants(k).LineWidth=.9;end
drawnow;
print(paper,[filename '.pdf'],'-dpdf','-painters');
print(paper,[filename '.png'],'-dpng',['-r' num2str(resolution)]);
savefig(paper,[filename '.fig']);
close(paper);
end
