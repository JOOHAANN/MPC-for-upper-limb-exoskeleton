function diagnostics=diagnose_angle_torque(root)
%DIAGNOSE_ANGLE_TORQUE Replay archived wall torque and decompose saved predictions.
% New-model components use saved analytic acceleration, never a second
% numerical derivative of the plotted/downsampled angle.
    file=fullfile(root,'results','bse_angle_priority','angle_priority.mat');
    assert(isfile(file),'Run demo_bse_angle_priority first.');
    d=load(file,'result'); result=d.result;
    out=fullfile(root,'results','bse_torque_breakdown');
    if ~isfolder(out), mkdir(out); end
    cfg=jsondecode(fileread(fullfile(root,'config','workflows', ...
        'training_map__muscle_tendon__general_model_op_scale.json')));
    names={'Combination_5_35_10_55_15_45.mat','Combination_5_45_10_35_15_55.mat', ...
        'Combination_5_55_10_45_15_35.mat'};
    arrays=cell(3,1);
    for k=1:3, d=load(asset(root,cfg,names{k})); arrays{k}=d.Combination; end
    data=vertcat(arrays{:});
    d=load(asset(root,cfg,'coeffs.mat')); c=num2cell(d.coeffs);
    oldPath=path; cleanup=onCleanup(@()path(oldPath));
    addpath(fileparts(asset(root,cfg,'Torque_function_constant.m')),'-begin');
    clear Torque_function_constant;
    assert(strcmp(which('Torque_function_constant'),asset(root,cfg,'Torque_function_constant.m')));
    x=Torque_function_constant(zeros(size(data,1),1),data(:,3),data(:,3), ...
        data(:,2),data(:,4),.5,.24086,.001,c{:});
    fprintf('Replayed archived BSE and recorded output column 6 (boundary torque).\n');
    windows=[5321 11740;22037 26059;36732 41109;52761 58381;66452 72417; ...
        84099 87884;98004 102357;113394 118437;127817 133717];
    archived=cell(1,9); components=cell(2,9); records=struct([]);
    labels={'Experimental activation','MAP activation'}; speeds=[35 45 55];
    for k=1:9
        t=result.trials(k); ix=(windows(k,1)+1:10:windows(k,2)-1)';
        archived{k}=x(ix,[1 2 6 10]);
        assert(isequal(data(ix,3),t.q));
        assert(max(abs(wall_law(x(ix,2))-x(ix,6)))<1e-8);
        % The recorded wall uses the accepted state at k, as in the original
        % output. Inner fixed-point iterates evaluate the same wall separately.
        for source=1:2
            p=result.predictions{source,2,k};
            model=result.models{source,find(speeds==t.speed)+1};
            M=22.6006; L1=.345; L2=.760606; phi=.936847;
            inertia=t.m*L2^2+4*M*L1^2/3;
            parts=[inertia.*p.acceleration,model.damping*p.v, ...
                9.81*M*L1*cos(p.q),9.81*t.m*L2.*cos(p.q-phi),zeros(size(p.q))];
            assert(max(abs(sum(parts,2)-p.torque))<1e-9,'Torque decomposition mismatch.');
            % Expose the virtual equilibrium drive as well: zero wall torque
            % does NOT imply the reduced model has no restoring action.
            drive=inertia.*(p.equilibrium-p.q)/model.tau^2;
            dynamicDamping=-2*inertia.*p.v/model.tau;
            assert(max(abs(drive+dynamicDamping-parts(:,1)))<1e-8);
            parts=[parts,drive,dynamicDamping];
            components{source,k}=parts;
            rec=struct('source',labels{source},'loadKg',t.loadKg,'speed',t.speed, ...
                'inertiaMaxAbs',max(abs(parts(:,1))),'dampingMaxAbs',max(abs(parts(:,2))), ...
                'bodyGravityMaxAbs',max(abs(parts(:,3))),'loadGravityMaxAbs',max(abs(parts(:,4))), ...
                'newBoundaryMaxAbs',0,'newTorqueMaxAbs',max(abs(p.torque)), ...
                'archivedBoundaryMaxAbs',max(abs(x(ix,6))), ...
                'archivedBoundaryActivePercent',100*mean(abs(x(ix,6))>1e-9), ...
                'archivedBalanceMaxAbs',max(abs(x(ix,10))), ...
                'equilibriumDriveMaxAbs',max(abs(drive)), ...
                'dynamicDampingMaxAbs',max(abs(dynamicDamping)), ...
                'equivalentStiffnessMax',max(inertia/model.tau^2));
            if isempty(records), records=rec; else, records(end+1)=rec; end %#ok<AGROW>
            samples=table(t.t,t.q,p.q,t.torque,p.torque,parts(:,1),parts(:,2),parts(:,3),parts(:,4), ...
                parts(:,5),x(ix,6),x(ix,10),'VariableNames',{'time','referenceAngle','predictedAngle', ...
                'referenceTorque','predictedTorque','inertiaTorque','dampingTorque','bodyGravity', ...
                'loadGravity','newBoundaryTorque','archivedBoundaryTorque','archivedBalanceTorque'});
            samples.virtualEquilibriumDrive=drive;
            samples.dynamicDampingTorque=dynamicDamping;
            writetable(samples,fullfile(out,sprintf('components_source%d_load%d_speed%d.csv',source,t.loadKg,t.speed)));
        end
    end
    for speed=speeds
        fig=figure('Visible','off','Color','w','Position',[20 20 1800 1000]);
        tiledlayout(3,4,'TileSpacing','compact');
        ids=find([result.trials.speed]==speed);
        for k=ids
            t=result.trials(k); p=result.predictions{1,2,k}; g=result.predictions{2,2,k}; a=archived{k};
            nexttile; plot(t.t,[t.q p.q g.q a(:,2)]*180/pi,'LineWidth',1.2); grid on;
            ylabel('Angle (deg)'); title(sprintf('%g kg / speed %g',t.loadKg,speed));
            if k==ids(1), legend('Reference','New: experimental','New: MAP','Archived BSE','Location','best'); end
            nexttile; plot(t.t,a(:,3),'LineWidth',1.4); hold on;
            plot(t.t,zeros(size(t.t)),'--','LineWidth',1.4); grid on; ylabel('Boundary torque (Nm)');
            title('Actual wall contribution');
            if k==ids(1), legend('Archived BSE (MAP)','New model: both inputs = 0','Location','best'); end
            for source=1:2
                nexttile; parts=components{source,k}; p=result.predictions{source,2,k};
                plot(t.t,[p.torque,parts(:,1:4)],'LineWidth',1.1); hold on;
                plot(t.t,t.torque,'k--','LineWidth',1.1); grid on; ylabel('Torque (Nm)'); title(labels{source});
                if k==ids(1), legend('Total','Inertia','Damping','Body gravity','Load gravity','Reference','Location','best'); end
            end
        end
        sgtitle('Held-out speed | Archived wall vs new-model torque components');
        exportgraphics(fig,fullfile(out,sprintf('components_speed%d.png',speed)),'Resolution',160);
        exportgraphics(fig,fullfile(out,sprintf('components_speed%d.pdf',speed)),'ContentType','vector');
        savefig(fig,fullfile(out,sprintf('components_speed%d.fig',speed))); close(fig);
        fig=figure('Visible','off','Color','w','Position',[30 30 1300 950]);
        tiledlayout(3,2,'TileSpacing','compact');
        for k=ids
            t=result.trials(k);
            for source=1:2
                parts=components{source,k}; nexttile;
                plot(t.t,parts(:,[6 7 1]),'LineWidth',1.2); grid on;
                ylabel('Torque (Nm)'); xlabel('Time (s)');
                title(sprintf('%g kg / %s',t.loadKg,labels{source}));
                if k==ids(1), legend('Virtual equilibrium drive','Dynamic damping','Net inertia torque','Location','best'); end
            end
        end
        sgtitle(sprintf('Speed %g | Internal impedance terms (not boundary contact forces)',speed));
        exportgraphics(fig,fullfile(out,sprintf('virtual_impedance_speed%d.png',speed)),'Resolution',160);
        exportgraphics(fig,fullfile(out,sprintf('virtual_impedance_speed%d.pdf',speed)),'ContentType','vector'); close(fig);
    end
    angle=linspace(-60,155,1200)';
    fig=figure('Visible','off','Color','w','Position',[50 50 1000 500]);
    plot(angle,wall_law(angle*pi/180),'LineWidth',1.7); hold on;
    plot(angle,zeros(size(angle)),'--','LineWidth',1.2); xline(0,':'); xline(95,':'); grid on;
    xlabel('Angle (deg)'); ylabel('Boundary torque (Nm)');
    title('Archived wall: k = 1200 Nm/rad, saturation = 900 Nm, damping = 0');
    legend('Archived BSE wall law','New model: no wall torque','Location','best');
    exportgraphics(fig,fullfile(out,'boundary_law.png'),'Resolution',160);
    exportgraphics(fig,fullfile(out,'boundary_law.pdf'),'ContentType','vector'); close(fig);
    diagnostics=struct('summary',struct2table(records),'components',{components},'archived',{archived});
    writetable(diagnostics.summary,fullfile(out,'summary.csv'));
    save(fullfile(out,'torque_diagnostics.mat'),'diagnostics');
    disp(diagnostics.summary);
    fprintf('PASS: all component sums and recorded archived wall values verified.\nOutputs: %s\n',out);
end
function tau=wall_law(q)
% Exact archived static wall law. Its c_wall is zero, so velocity is absent.
    upper=q-95*pi/180; lower=-q;
    tau=-900*tanh((1200/900)*max(upper,0)).*(.5*(1+tanh(30*upper))) ...
        +900*tanh((1200/900)*max(lower,0)).*(.5*(1+tanh(30*lower)));
end
function file=asset(root,cfg,name)
    item=cfg.files(strcmp({cfg.files.target},name)); assert(isscalar(item));
    file=fullfile(root,item.source);
end
