function test_demo_figure_format
% Check the shared style and the FIG files created by both demo scripts.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'Demo'));
f=figure('Visible','off');layout=tiledlayout(f,1,1);ax=nexttile(layout);
p=plot(ax,[0 1],[0 1],'LineWidth',1.2);hold(ax,'on');
limit=yline(ax,.5,'--','LineWidth',.8);ax.LineWidth=.6;
xlabel(ax,'Time');ylabel(ax,'Torque');title(ax,'Axis title');
key=legend(ax,[p limit],{'Signal','Limit'},'Location','northwest');
title(layout,'Layout title');subtitle(layout,'Layout subtitle');
drawnow;oldX=ax.XLim;oldY=ax.YLim;
format_demo_figure(f);
assert(abs(p.LineWidth-1.8)<1e-10 && abs(limit.LineWidth-1.2)<1e-10);
assert(abs(ax.LineWidth-.9)<1e-10 && isequal(ax.XLim,oldX) && isequal(ax.YLim,oldY));
assert(abs(ax.FontSize-24.48)<1e-10 && abs(ax.Title.FontSize-20.4)<1e-10);
assert(abs(ax.XAxis.FontSize-24.48)<1e-10 && abs(ax.YAxis.FontSize-24.48)<1e-10);
assert(abs(ax.XLabel.FontSize-16.32)<1e-10 && abs(ax.YLabel.FontSize-16.32)<1e-10);
assert(abs(layout.Title.FontSize-20.4)<1e-10 && abs(layout.Subtitle.FontSize-16.32)<1e-10);
assert(strcmp(key.Location,'northoutside') && abs(key.FontSize-16.32)<1e-10);
format_demo_figure(f);
assert(abs(p.LineWidth-1.8)<1e-10 && abs(limit.LineWidth-1.2)<1e-10, ...
    'Formatting an existing figure a second time must not thicken it again.');
close(f);

folders={fullfile(root,'Demo','NMPC','output','current_BSE_528'), ...
    fullfile(root,'Demo','BSE','Direct_Activation_No_Limits','output')};
expected=[16 10];
for group=1:numel(folders)
    files=dir(fullfile(folders{group},'*.fig'));
    assert(numel(files)>=expected(group),'The expected demo FIG files are missing.');
    for k=1:numel(files)
        fig=openfig(fullfile(files(k).folder,files(k).name),'invisible');
        assert(isappdata(fig,'DemoLineWidthsScaled150'),'%s was not formatted.',files(k).name);
        items=findall(fig,'-property','FontSize');
        for j=1:numel(items)
            % MATLAB may recreate empty, invisible axes subtitle objects in
            % auto mode when a FIG is opened. They do not render any text.
            if isa(items(j),'matlab.graphics.primitive.Text') && isempty(items(j).String)
                continue;
            end
            sizePt=items(j).FontSize;
            assert(abs(sizePt-16.32)<1e-6 || abs(sizePt-20.4)<1e-6 || abs(sizePt-24.48)<1e-6, ...
                '%s has an inconsistent font size (%g pt, %s).',files(k).name,sizePt,class(items(j)));
        end
        axesList=findall(fig,'Type','axes');
        for j=1:numel(axesList)
            a=axesList(j);
            assert(abs(a.FontSize-24.48)<1e-6 && abs(a.Title.FontSize-20.4)<1e-6);
            assert(abs(a.XAxis.FontSize-24.48)<1e-6 && abs(a.YAxis.FontSize-24.48)<1e-6);
            assert(abs(a.XLabel.FontSize-16.32)<1e-6 && abs(a.YLabel.FontSize-16.32)<1e-6);
        end
        layouts=findall(fig,'Type','tiledlayout');
        for j=1:numel(layouts)
            assert(abs(layouts(j).Title.FontSize-20.4)<1e-6);
            assert(abs(layouts(j).Subtitle.FontSize-16.32)<1e-6);
        end
        legends=findall(fig,'Type','legend');
        for j=1:numel(legends)
            assert(abs(legends(j).FontSize-16.32)<1e-6);
            assert(~any(strcmp(legends(j).Location,{'northwest','northeast','southwest','southeast','best'})));
        end
        close(fig);
    end
    fprintf('PASS: %d formatted figures in %s\n',numel(files),folders{group});
end
end
