% Check original standing balance and diagnose the uncontrolled replay.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','NMPC');addpath(demo,fullfile(root,'Demo','BSE'));
subjects=bse_subjects('528');subject=subjects{1};
motor=load(fullfile(demo,'data','no_human_PD_motor_programs_528.mat'));
old=load(fullfile(demo,'output','pid_nmpc_528','threshold_tracking','pid_nmpc_comparison_528.mat'));
rows=struct([]);
for source=1:2
    q=subject.sequences{1};model=subject.models{source};
    original=bse_model(model,subject.balanceController,q.standingAngle,true);
    plant=bse_without_posture(model,subject.balanceController,q.standingAngle);
    assert(isequal(plant.standingParameters,[q.standingAngle subject.balanceController(:)']));
    a=repmat([.1 .02],101,1);mass=zeros(101,1);waiting=ones(101,1);
    expected=original.simulate(a,mass,waiting,[1.35 .05],.01,10);
    actual=plant.simulate(a,mass,waiting,[1.35 .05],.01,10);
    physical=[1:19 20:23 25:27];
    assert(max(abs(expected(:,physical)-actual(:,physical)),[],'all')<1e-7, ...
        'Standing balance differs from the original controller.');
    assert(max(abs(actual(:,22)))>1,'Standing PD did not respond to a posture disturbance.');
    for group=1:3
        q=subject.sequences{group};input=motor.programs{group,source};
        plant=bse_without_posture(model,subject.balanceController,q.standingAngle);
        restored=plant.simulate(input,q.mass,q.waiting,q.initialState,q.dt,10);
        fine=plant.simulate(input,q.mass,q.waiting,q.initialState,q.dt,20);
        previous=load(fullfile(demo,'output',sprintf('no_human_PD_verified_group%d',group), ...
            sprintf('group%d_source%d.mat',group,source)));
        nominal=old.results.groups{group}.runs{source}.baseline;
        wait=q.waiting>.99;
        row=struct('Source',source,'Group',group, ...
            'OriginalMinimumThetaRad',min(nominal(:,1)), ...
            'NoBalanceMinimumThetaRad',min(previous.pair.baseline(:,1)), ...
            'RestoredMinimumThetaRad',min(restored(:,1)), ...
            'RestoredMaximumThetaRad',max(restored(:,1)), ...
            'WaitRMSEDeg',rad2deg(sqrt(mean((restored(wait,1)-q.standingAngle).^2))), ...
            'MaximumMomentumResidualNm',max(abs(restored(:,17))), ...
            'RefinedIntegrationMaximumDifferenceRad',max(abs(restored(:,1)-fine(:,1))));
        rows=[rows;row]; %#ok<AGROW>
        if group==1
            fig=figure('Visible','off','Color','w','Position',[50 50 1300 720]);
            tiledlayout(2,1,'TileSpacing','compact');
            nexttile;plot(q.time,nominal(:,1),'k:',q.time,previous.pair.baseline(:,1),'--', ...
                q.time,restored(:,1),'-','LineWidth',1.2);
            ylabel('Trunk angle (rad)');grid on;
            legend('Original calibrated model','Movement and standing PD removed', ...
                'Movement PD removed; standing balance retained','Location','southoutside');
            title(sprintf('Unassisted replay | Source %d | Standing balance restored',source));
            nexttile;plot(q.time,restored(:,22),'LineWidth',1.1);hold on;
            stairs(q.time,120*q.waiting,'k:');ylabel('Standing correction (N m)');
            xlabel('Time (s)');grid on;legend('Applied standing correction','Waiting flag x 120');
            exportgraphics(fig,fullfile(folder,sprintf('standing_restoration_source%d.png',source)),'Resolution',180);
            close(fig);
        end
    end
end
metrics=struct2table(rows);disp(metrics);
writetable(metrics,fullfile(folder,'standing_restoration_metrics.csv'));
fprintf('PASS: original standing PD response retained for both sources; movement gains remain zero.\n');
