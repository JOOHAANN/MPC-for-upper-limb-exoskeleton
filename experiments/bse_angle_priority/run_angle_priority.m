function result = run_angle_priority(root)
%RUN_ANGLE_PRIORITY Paired-input calibration and whole-speed validation.
    out=fullfile(root,'results','bse_angle_priority'); if ~isfolder(out), mkdir(out); end
    cfg=jsondecode(fileread(fullfile(root,'config','workflows', ...
        'training_map__muscle_tendon__general_model_op_scale.json')));
    names={'5_35_10_55_15_45','5_45_10_35_15_55','5_55_10_45_15_35'};
    raw=cell(3,1); generated=raw;
    for j=1:3
        item=cfg.files(strcmp({cfg.files.target},['Combination_' names{j} '.mat']));
        d=load(fullfile(root,item.source)); generated{j}=d.Combination;
        d=load(fullfile(root,'data','processed', ...
            ['training_map__muscle_tendon__' names{j} '_done__build_data'],'Combination_raw.mat'));
        raw{j}=d.Combination;
        assert(isequal(raw{j}(:,3),generated{j}(:,3)) && isequal(raw{j}(:,1),generated{j}(:,4)), ...
            'The two activation sources do not share identical references.');
    end
    raw=vertcat(raw{:}); generated=vertcat(generated{:});
    windows=[5321 11740;22037 26059;36732 41109;52761 58381;66452 72417; ...
        84099 87884;98004 102357;113394 118437;127817 133717];
    contact=[6742 10348;22885 25303;37516 40124;53975 57348;67641 71085; ...
        85027 87002;99146 101117;114510 117291;128682 132730];
    conditions=[5 35;10 55;15 45;5 45;10 35;15 55;5 55;10 45;15 35];
    for k=1:9
        ix=(windows(k,1)+1:10:windows(k,2)-1)';
        trials(k)=struct('a',[raw(ix,2),generated(ix,2)],'q',raw(ix,3), ...
            'torque',raw(ix,1),'m',conditions(k,1)*double(ix>contact(k,1)&ix<contact(k,2)), ...
            'loadKg',conditions(k,1),'speed',conditions(k,2),'t',(0:numel(ix)-1)'*.01); %#ok<AGROW>
    end
    labels={'Experimental activation','MAP activation'};
    predictions=cell(2,2,9); models=cell(2,4); records=struct([]);
    for source=1:2
        for fold=0:3
            if fold==0, train=1:9; test=1:9; else
                speeds=[35 45 55]; test=find([trials.speed]==speeds(fold)); train=setdiff(1:9,test);
            end
            model=fit_model(trials(train),source); models{source,fold+1}=model;
            fprintf('%s, fold %d: objective %.5f, exit %d\n',labels{source},fold,model.loss,model.exitflag);
            for k=test
                t=trials(k); p=angle_rollout(model,t.a(:,source),t.m,t.loadKg,t.q(1),.01);
                assert(all(isfinite(p.q)) && min(p.q)>-1e-9 && max(p.q)<model.limit+1e-9);
                mode=1+(fold>0); predictions{source,mode,k}=p;
                rec=metrics(t,p); rec.source=labels{source}; rec.mode=mode; rec.loadKg=t.loadKg; rec.speed=t.speed;
                if isempty(records), records=rec; else, records(end+1)=rec; end %#ok<AGROW>
            end
        end
    end
    result=struct('trials',trials,'models',{models},'predictions',{predictions},'metrics',struct2table(records));
    % Optional archived BSE overlay uses its original MAP-input replay.
    baselineFile=fullfile(root,'results','bse_torque_theta','joint_comparison.mat');
    if isfile(baselineFile)
        d=load(baselineFile,'results');
        for k=1:9
            assert(isequal(d.results.trials(k).q,trials(k).q),'Baseline reference mismatch.');
        end
        result.archived=d.results.archived;
    end
    writetable(result.metrics,fullfile(out,'metrics.csv'));
    save(fullfile(out,'angle_priority.mat'),'result');
    draw_figures(result,out);
    disp(groupsummary(result.metrics,{'source','mode'},'mean', ...
        {'thetaRMSE','torqueRMSE','shapeCorrelation','minimumError','maximumError','Wmatched'}));
    fprintf('Results: %s\n',out);
end

function model=fit_model(trials,source)
    X=[]; y=[];
    for t=trials
        f=angle_features(t.a(:,source),t.m,t.loadKg,.01);
        X=[X;f]; y=[y;t.q]; %#ok<AGROW>
    end
    model=struct('limit',95*pi/180,'tau',.15,'damping',.5,'mu',mean(X), ...
        'scale',max(std(X),.01),'beta',[]);
    Z=[ones(size(y)),(X-model.mu)./model.scale];
    target=min(.995,max(.005,y/model.limit)); target=log(target./(1-target));
    beta=(Z'*Z+10*eye(size(Z,2)))\(Z'*target);
    opts=optimoptions('lsqnonlin','Display','off','MaxIterations',100, ...
        'MaxFunctionEvaluations',4000,'FunctionTolerance',1e-6);
    [model.beta,model.loss,~,model.exitflag]=lsqnonlin(@loss,beta,-20*ones(size(beta)),20*ones(size(beta)),opts);
    function r=loss(b)
        candidate=model; candidate.beta=b; r=[];
        for t=trials
            p=angle_rollout(candidate,t.a(:,source),t.m,t.loadKg,t.q(1),.01);
            % Trial-balanced angle and slope losses dominate the torque term.
            n=sqrt(numel(t.q));
            e=[(p.q-t.q)/.15; .5*(p.v-gradient(t.q,.01))/.5; ...
                .15*(p.torque-t.torque)/100]/n;
            landmarks=reference_landmarks(t.q);
            e=[e;.8*(p.q(landmarks)-t.q(landmarks))/.15/sqrt(numel(landmarks))];
            r=[r;e/sqrt(numel(trials))]; %#ok<AGROW>
        end
        r=[r;.025*b(2:end)];
    end
end

function ix=reference_landmarks(q)
    [~,lo]=findpeaks(-q,'MinPeakProminence',.12,'MinPeakDistance',50);
    [~,hi]=findpeaks(q,'MinPeakProminence',.12,'MinPeakDistance',50);
    ix=unique([1;lo;hi;numel(q)]);
end
function r=metrics(t,p)
    c=corrcoef(t.q,p.q);
    [~,ref]=findpeaks(-t.q,'MinPeakProminence',.12,'MinPeakDistance',50);
    [~,pred]=findpeaks(-p.q,'MinPeakProminence',.12,'MinPeakDistance',50);
    timing=NaN; extrema=NaN; peakError=NaN; peakTime=NaN;
    matched=numel(ref)==2 && numel(pred)==2;
    if matched
        timing=mean(abs(ref-pred))*.01; extrema=mean(abs(t.q(ref)-p.q(pred)));
        [rpeak,ri]=max(t.q(ref(1):ref(2))); [ppeak,pi]=max(p.q(pred(1):pred(2)));
        peakError=abs(rpeak-ppeak); peakTime=abs(ref(1)+ri-pred(1)-pi)*.01;
    end
    r=struct('thetaRMSE',sqrt(mean((p.q-t.q).^2)), ...
        'torqueRMSE',sqrt(mean((p.torque-t.torque).^2)),'shapeCorrelation',c(1,2), ...
        'minimumError',abs(min(p.q)-min(t.q)),'maximumError',abs(max(p.q)-max(t.q)), ...
        'minAngle',min(p.q),'maxAngle',max(p.q),'maxSpeed',max(abs(p.v)), ...
        'referenceValleys',numel(ref),'predictedValleys',numel(pred),'Wmatched',double(matched), ...
        'valleyTimeMAE',timing,'valleyAngleMAE',extrema, ...
        'middlePeakAngleError',peakError,'middlePeakTimeError',peakTime);
end
function draw_figures(result,out)
    modes={'Calibration (all nine trials)','Held-out speed validation'};
    for mode=1:2
        for speed=[35 45 55]
            fig=figure('Color','w','Position',[50 50 1500 950],'Visible','off');
            tiledlayout(3,3,'TileSpacing','compact');
            ids=find([result.trials.speed]==speed);
            for k=ids
                t=result.trials(k); p=result.predictions{1,mode,k}; g=result.predictions{2,mode,k};
                nexttile; plot(t.t,t.a,'LineWidth',1.1); ylabel('Activation'); title(sprintf('%g kg, speed %g',t.loadKg,speed)); grid on;
                if k==ids(1), legend('Experimental','MAP','Location','best'); end
                nexttile; plot(t.t,[t.q p.q g.q]*180/pi,'LineWidth',1.4); ylabel('Angle (deg)'); grid on;
                names={'Reference','Experimental input','MAP input'};
                if isfield(result,'archived')
                    hold on; plot(t.t,result.archived{k}(:,2)*180/pi,'--','Color',[.5 .5 .5]);
                    names{4}='Archived BSE (MAP)';
                end
                if k==ids(1), legend(names,'Location','best'); end
                nexttile; plot(t.t,[t.torque p.torque g.torque],'LineWidth',1.1); ylabel('Torque (Nm)'); grid on;
                if isfield(result,'archived'), hold on; plot(t.t,result.archived{k}(:,1),'--','Color',[.5 .5 .5]); end
                if k==ids(1), legend(names,'Location','best'); end
                xlabel('Time (s)');
            end
            sgtitle([modes{mode} ' | Same bounded model, separately fitted parameters']);
            exportgraphics(fig,fullfile(out,sprintf('comparison_mode%d_speed%d.png',mode,speed)),'Resolution',160);
            exportgraphics(fig,fullfile(out,sprintf('comparison_mode%d_speed%d.pdf',mode,speed)),'ContentType','vector'); close(fig);
        end
    end
end
