function result=run_hill_forward(root,refit,refine)
%RUN_HILL_FORWARD Same Hill structure, separately calibrated activation sources.
    if nargin<2, refit=false; end
    if nargin<3, refine=false; end
    out=fullfile(root,'results','hill_forward'); if ~isfolder(out), mkdir(out); end
    trials=load_trials(root);
    cacheFile=fullfile(root,'data','processed','hill_forward','calibrated_models.mat');
    cached=[];
    if ~refit && isfile(cacheFile)
        cached=load(cacheFile,'models','trials');
        assert(isequaln(trials,cached.trials),'Input data changed; call run_hill_forward(root,true) to refit.');
        fprintf('Replaying calibrated Hill models. Use run_hill_forward(root,true) to refit.\n');
    end
    models=cell(2,4); predictions=cell(2,2,9); records=struct([]);
    labels={'Experimental activation','MAP activation'}; speeds=[35 45 55];
    for source=1:2
        for fold=0:3
            if fold==0, train=1:9; test=1:9; else
                test=find([trials.speed]==speeds(fold)); train=setdiff(1:9,test);
            end
            if isempty(cached), model=fit_hill_parameters(trials(train),source); else, model=cached.models{source,fold+1}; end
            model.trainIds=train; model.testIds=test;
            models{source,fold+1}=model;
            fprintf('%s fold %d: loss %.4f parameters %s\n',labels{source},fold,model.loss,mat2str(model.parameters,5));
            for k=test
                t=trials(k); p=hill_rollout(model.parameters,t.a(:,source),t.m,[t.q(1) t.v0],.01,2);
                assert(~p.failed,'A selected Hill rollout failed; do not interpret its optimizer penalty states as predictions.');
                predictions{source,1+(fold>0),k}=p;
                r=score(t,p); r.source=labels{source}; r.mode=1+(fold>0); r.loadKg=t.loadKg; r.speed=t.speed;
                if isempty(records), records=r; else, records(end+1)=r; end %#ok<AGROW>
            end
            save(fullfile(out,'checkpoint.mat'),'models','predictions');
        end
    end
    if refine
        % Transfer ONLY an initialization from the paired source fitted on
        % the SAME training trials; no held-out reference enters this step.
        for source=1:2
            for fold=0:3
                old=models{source,fold+1}; other=models{3-source,fold+1};
                assert(isequal(old.trainIds,other.trainIds));
                fitted=fit_hill_parameters(trials(old.trainIds),source,other.parameters);
                if fitted.loss<old.loss
                    fitted.trainIds=old.trainIds; fitted.testIds=old.testIds;
                    models{source,fold+1}=fitted;
                    for k=old.testIds
                        t=trials(k);
                        predictions{source,1+(fold>0),k}=hill_rollout(fitted.parameters,t.a(:,source),t.m,[t.q(1) t.v0],.01,2);
                    end
                end
                fprintf('Paired initialization %s fold %d: %.4f -> %.4f\n',labels{source},fold,old.loss,models{source,fold+1}.loss);
            end
        end
    end
    records=struct([]);
    for source=1:2
        for mode=1:2
            for k=1:9
                t=trials(k); p=predictions{source,mode,k};
                assert(~p.failed,'The candidate fails a forward rollout. Saved baseline parameters have NOT been replaced.');
                rec=score(t,p); rec.source=labels{source}; rec.mode=mode; rec.loadKg=t.loadKg; rec.speed=t.speed;
                if isempty(records), records=rec; else, records(end+1)=rec; end %#ok<AGROW>
            end
        end
    end
    result=struct('trials',trials,'models',{models},'predictions',{predictions},'metrics',struct2table(records));
    oldFile=fullfile(root,'results','bse_torque_theta','joint_comparison.mat');
    if isfile(oldFile)
        d=load(oldFile,'results');
        for k=1:9, assert(isequal(trials(k).q,d.results.trials(k).q)); end
        result.archived=d.results.archived;
    end
    save(fullfile(out,'hill_forward.mat'),'result'); writetable(result.metrics,fullfile(out,'metrics.csv'));
    parameterRecords=struct([]);
    for source=1:2
        for fold=0:3
            model=models{source,fold+1}; par=model.parameters;
            rec=struct('source',labels{source},'fold',fold,'F0',par(1),'lo',par(2),'R',par(3), ...
                'lref',par(4),'B',par(5),'vmax',par(6),'gamma',par(7),'loss',model.loss,'exitflag',model.exitflag);
            if isempty(parameterRecords), parameterRecords=rec; else, parameterRecords(end+1)=rec; end %#ok<AGROW>
        end
    end
    writetable(struct2table(parameterRecords),fullfile(out,'parameters.csv'));
    if ~isfolder(fileparts(cacheFile)), mkdir(fileparts(cacheFile)); end
    if isempty(cached) || refine, save(cacheFile,'models','trials'); end
    plot_results(result,out);
    disp(groupsummary(result.metrics,{'source','mode'},'mean',{'thetaRMSE','torqueRMSE','Wmatched','outsidePercent'}));
end
function trials=load_trials(root)
    cfg=jsondecode(fileread(fullfile(root,'config','workflows','training_map__muscle_tendon__general_model_op_scale.json')));
    names={'5_35_10_55_15_45','5_45_10_35_15_55','5_55_10_45_15_35'};
    raw=cell(3,1); generated=raw;
    for j=1:3
        item=cfg.files(strcmp({cfg.files.target},['Combination_' names{j} '.mat']));
        d=load(fullfile(root,item.source)); generated{j}=d.Combination;
        d=load(fullfile(root,'data','processed',['training_map__muscle_tendon__' names{j} '_done__build_data'],'Combination_raw.mat'));
        raw{j}=d.Combination;
        assert(isequal(raw{j}(:,3),generated{j}(:,3)) && isequal(raw{j}(:,1),generated{j}(:,4)));
    end
    raw=vertcat(raw{:}); generated=vertcat(generated{:});
    windows=[5321 11740;22037 26059;36732 41109;52761 58381;66452 72417;84099 87884;98004 102357;113394 118437;127817 133717];
    contact=[6742 10348;22885 25303;37516 40124;53975 57348;67641 71085;85027 87002;99146 101117;114510 117291;128682 132730];
    conditions=[5 35;10 55;15 45;5 45;10 35;15 55;5 55;10 45;15 35];
    for k=1:9
        ix=(windows(k,1)+1:10:windows(k,2)-1)'; a=[raw(ix,2),generated(ix,2)];
        % Fixed causal load transfer, 80 ms, identical for both sources.
        % This replaces an instantaneous attached-load switch; it is an
        % assumed effective-load schedule, not a measured contact force.
        m=conditions(k,1)*double(ix>contact(k,1)&ix<contact(k,2)); r=exp(-.01/.08); m=filter(1-r,[1 -r],m);
        % Initial velocity uses only the frame BEFORE the prediction starts.
        v0=(raw(ix(1),3)-raw(ix(1)-10,3))/.01;
        trials(k)=struct('a',a,'q',raw(ix,3),'v0',v0,'torque',raw(ix,1),'m',m, ...
            'loadKg',conditions(k,1),'speed',conditions(k,2),'t',(0:numel(ix)-1)'*.01); %#ok<AGROW>
    end
end

function r=score(t,p)
    [~,lo]=findpeaks(-p.q,'MinPeakProminence',.12,'MinPeakDistance',50);
    c=corrcoef(t.q,p.q);
    r=struct('thetaRMSE',sqrt(mean((p.q-t.q).^2)),'torqueRMSE',sqrt(mean((p.torque-t.torque).^2)), ...
        'Wmatched',double(numel(lo)==2),'valleys',numel(lo),'correlation',c(1,2), ...
        'minimumError',abs(min(p.q)-min(t.q)),'maximumError',abs(max(p.q)-max(t.q)), ...
        'outsidePercent',100*mean(p.q<0|p.q>95*pi/180),'minAngle',min(p.q),'maxAngle',max(p.q), ...
        'maxTorque',max(p.torque),'passiveTorqueMax',max(p.passiveTorque), ...
        'passiveTorqueImpulseFraction',trapz(p.passiveTorque)/max(trapz(p.torque),eps),'failed',p.failed);
end
function plot_results(r,out)
    for speed=[35 45 55]
        fig=figure('Visible','off','Color','w','Position',[30 30 1900 950]); tiledlayout(3,4,'TileSpacing','compact');
        for k=find([r.trials.speed]==speed)
            t=r.trials(k); p=r.predictions{1,2,k}; g=r.predictions{2,2,k};
            nexttile; plot(t.t,t.a,'LineWidth',1.2); grid on; ylabel('Activation');
            legend('Experimental','MAP','Location','best');
            nexttile; plot(t.t,[t.q p.q g.q]*180/pi,'LineWidth',1.2); grid on; ylabel('Angle (deg)');
            names={'Reference','Hill: experimental','Hill: MAP'};
            if isfield(r,'archived'), hold on; plot(t.t,r.archived{k}(:,2)*180/pi,'--','Color',[.5 .5 .5]); names{4}='Archived BSE (MAP)'; end
            title(sprintf('%g kg / speed %g',t.loadKg,speed)); legend(names,'Location','best');
            nexttile; plot(t.t,[t.torque p.torque g.torque],'LineWidth',1.2); grid on; ylabel('Muscle torque (Nm)');
            if isfield(r,'archived'), hold on; plot(t.t,r.archived{k}(:,1),'--','Color',[.5 .5 .5]); end
            nexttile; plot(t.t,[p.activeTorque p.passiveTorque g.activeTorque g.passiveTorque],'LineWidth',1.1); grid on;
            ylabel('Hill torque components (Nm)'); legend('Exp active','Exp passive','MAP active','MAP passive','Location','best');
        end
        sgtitle('Held-out speed | Same Hill forward model | No angle target, reference feedback or wall torque');
        exportgraphics(fig,fullfile(out,sprintf('heldout_speed%d.png',speed)),'Resolution',160);
        exportgraphics(fig,fullfile(out,sprintf('heldout_speed%d.pdf',speed)),'ContentType','vector'); savefig(fig,fullfile(out,sprintf('heldout_speed%d.fig',speed))); close(fig);
    end
end
