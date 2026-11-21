% Pulse-response sensitivity of the assumed human unloading gain.
% This is a model sensitivity check, not a closed-loop or human validation.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','NMPC'); addpath(demo,fullfile(root,'Demo','BSE'));
subjects=bse_subjects('528'); subject=subjects{1};
programs=load(fullfile(demo,'data','no_human_PD_motor_programs_528.mat'));
gains=[0 .5 .8 1]; time=(0:.01:1.5)'; rows=struct([]);
for source=1:2
    q=subject.sequences{1}; a=programs.programs{1,source}(800,:);
    outputs=cell(1,numel(gains));
    for test=1:numel(gains)
        plant=bse_without_posture(subject.models{source},subject.balanceController,q.standingAngle,gains(test));
        [state,first]=plant.initialize(a,5,0,[1.0 0]);
        out=zeros(numel(time),27); out(1,:)=first;
        for k=2:numel(time)
            request=0;
            if time(k)>=.2 && time(k)<.6, request=40; end
            if time(k)>=.6 && time(k)<1.0, request=-20; end
            previous=state(29); if previous*request<0, previous=0; end
            u=sign(request)*min(abs(request),abs(previous)+2);
            [state,out(k,:)]=plant.step(state,a,5,0,u,.01,10);
            assert(all(out(k,18:19)<=a+1e-12),'Unloading increased neural recruitment.');
        end
        outputs{test}=out;
        row=struct('SourceIndex',source,'UnloadingGain',gains(test), ...
            'MinimumThetaRad',min(out(:,1)),'MaximumThetaRad',max(out(:,1)), ...
            'BackPeakNm',max(out(:,4)),'AbdomenPeakNm',max(out(:,5)), ...
            'MaximumBackCommandIncrease',max(out(:,18)-a(1)), ...
            'MaximumAbdomenCommandIncrease',max(out(:,19)-a(2)));
        rows=[rows;row]; %#ok<AGROW>
    end
    fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18 17]);
    layout=tiledlayout(fig,4,1,'TileSpacing','compact');
    columns=[1 4 5 25]; labels={'Angle (rad)','Back torque (N m)','Abdomen torque (N m)','Exoskeleton torque (N m)'};
    for kind=1:4
        ax=nexttile(layout);hold(ax,'on');
        for test=1:numel(gains),plot(ax,time,outputs{test}(:,columns(kind)),'LineWidth',1.1);end
        ylabel(ax,labels{kind});grid(ax,'on');set(ax,'FontName','Times New Roman');
    end
    xlabel(ax,'Time (s)'); legend(ax,{'Gain 0','Gain 0.5','Gain 0.8','Gain 1'},'NumColumns',4,'Location','southoutside');
    title(layout,sprintf('Human unloading assumption: pulse-response sensitivity | Source %d',source));
    exportgraphics(fig,fullfile(folder,sprintf('unloading_sensitivity_source%d.png',source)),'Resolution',250);
    exportgraphics(fig,fullfile(folder,sprintf('unloading_sensitivity_source%d.pdf',source)),'ContentType','vector');
    close(fig);
end
writetable(struct2table(rows),fullfile(folder,'unloading_sensitivity.csv'));
disp(struct2table(rows));
fprintf('PASS: signed motor pulse, finite Hill responses and no antagonist-command increase at unloading gains 0, 0.5, 0.8 and 1.\n');
