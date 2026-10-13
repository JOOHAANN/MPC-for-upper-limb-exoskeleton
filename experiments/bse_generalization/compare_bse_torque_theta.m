function results = compare_bse_torque_theta(root, referenceVariant)
%COMPARE_BSE_TORQUE_THETA Joint torque/angle benchmark with explicit limitations.
% Archived BSE parameters and upstream MAP preprocessing have seen these
% conditions. This is a conditional comparison, not new-subject validation.
    if nargin < 2, referenceVariant = 'auto'; end
    setup_project();
    out = fullfile(root, 'results', 'bse_torque_theta');
    if ~isfolder(out), mkdir(out); end
    ids = {'general_model', 'general_model_op_damp', 'general_model_op_scale'};
    assert(ismember(referenceVariant, [{'auto'}, ids]), 'Unknown reference variant.');
    arrays = cell(3, 1);
    cfg = get_config(root, ids{3});
    names = {'Combination_5_35_10_55_15_45.mat', ...
        'Combination_5_45_10_35_15_55.mat', 'Combination_5_55_10_45_15_35.mat'};
    for k = 1:3
        d = load(asset(root, cfg, names{k}), 'Combination'); arrays{k} = d.Combination;
    end
    data = vertcat(arrays{:});
    windows = [5321 11740;22037 26059;36732 41109;52761 58381; ...
        66452 72417;84099 87884;98004 102357;113394 118437;127817 133717];
    conditions = [5 35;10 55;15 45;5 45;10 35;15 55;5 55;10 45;15 35];
    contact = [6742 10348;22885 25303;37516 40124;53975 57348; ...
        67641 71085;85027 87002;99146 101117;114510 117291;128682 132730];
    trials = struct([]);
    for k = 1:9
        idx = (windows(k,1)+1:10:windows(k,2)-1)';
        trials(k).indices = idx;
        trials(k).t = (0:numel(idx)-1)'*0.01;
        trials(k).a = data(idx,2);
        trials(k).q = data(idx,3);
        trials(k).torque = data(idx,4);
        trials(k).loadKg = conditions(k,1);
        trials(k).speed = conditions(k,2);
        trials(k).m = conditions(k,1)*double(idx>contact(k,1) & idx<contact(k,2));
    end
    archived = cell(3,9);
    candidateMetrics = struct([]);
    oldPath = path; cleanup = onCleanup(@() path(oldPath));
    for model = 1:3
        cfg = get_config(root, ids{model});
        d = load(asset(root,cfg,'coeffs.mat'), 'coeffs'); c = num2cell(d.coeffs);
        source = asset(root,cfg,'Torque_function_constant.m');
        path(oldPath); addpath(fileparts(source), '-begin'); clear Torque_function_constant;
        assert(strcmp(which('Torque_function_constant'),source));
        x = Torque_function_constant(zeros(size(data,1),1), data(:,3),data(:,3), ...
            data(:,2),data(:,4),0.5,0.24086,0.001,c{:});
        for k = 1:9
            archived{model,k} = x(trials(k).indices,1:2);
            record = score_record(ids{model},trials(k),archived{model,k});
            if isempty(candidateMetrics), candidateMetrics = record; else, candidateMetrics(end+1) = record; end
        end
        fprintf('Replayed archived variant: %s\n',ids{model});
    end
    path(oldPath);
    candidates = struct2table(candidateMetrics);
    writetable(candidates,fullfile(out,'archived_candidates.csv'));
    % Whole speed groups are withheld. Tune dynamics using the other speed
    % group inside each outer-training set; never use test reference angles.
    speeds = [35 45 55];
    omega = [2 4 8 12]; damping = [0.7 1 1.5];
    newPrediction = cell(1,9); selectedOld = cell(1,9); folds = struct([]);
    for f = 1:3
        test = find([trials.speed]==speeds(f));
        train = find([trials.speed]~=speeds(f));
        oldScore = zeros(1,3);
        for m = 1:3
            for k = train
                metric = score_record('',trials(k),archived{m,k});
                oldScore(m) = oldScore(m)+metric.jointScore/numel(train);
            end
        end
        [~,bestOld] = min(oldScore);
        if ~strcmp(referenceVariant,'auto'), bestOld = find(strcmp(ids,referenceVariant)); end
        innerSpeeds = unique([trials(train).speed]);
        bestScore = Inf;
        for w = omega
            for zeta = damping
                value = 0;
                for inner = innerSpeeds
                    fitIds = train([trials(train).speed]~=inner);
                    valIds = train([trials(train).speed]==inner);
                    fitted = fit_surrogate(trials(fitIds),w,zeta);
                    for k = valIds
                        prediction = rollout_surrogate(fitted,trials(k));
                        metric = score_record('',trials(k),prediction);
                        value = value + metric.jointScore/numel(train);
                    end
                end
                if value<bestScore
                    bestScore = value; bestW = w; bestZeta = zeta;
                end
            end
        end
        fitted = fit_surrogate(trials(train),bestW,bestZeta);
        for k = test
            newPrediction{k} = rollout_surrogate(fitted,trials(k));
            poisoned = trials(k);
            poisoned.q(2:end) = NaN; poisoned.torque(:) = NaN;
            assert(isequaln(newPrediction{k},rollout_surrogate(fitted,poisoned)), ...
                'Future reference leaked into the new rollout.');
            selectedOld{k} = archived{bestOld,k};
        end
        folds(f).heldOutSpeed = speeds(f);
        folds(f).trainTrials = train; folds(f).testTrials = test;
        folds(f).archivedVariant = ids{bestOld};
        folds(f).omega = bestW; folds(f).zeta = bestZeta;
        folds(f).surrogate = fitted;
        fprintf('Held-out speed %g: archived=%s, omega=%g, zeta=%g\n',speeds(f),ids{bestOld},bestW,bestZeta);
    end
    records = struct([]);
    for k = 1:9
        record = score_record('Selected archived BSE',trials(k),selectedOld{k});
        if isempty(records), records = record; else, records(end+1) = record; end
        records(end+1) = score_record('New reduced state model',trials(k),newPrediction{k});
    end
    metrics = struct2table(records);
    writetable(metrics,fullfile(out,'torque_theta_metrics.csv'));
    results = struct('trials',trials,'archivedIds',{ids},'candidates',candidates, ...
        'folds',folds,'archived',{selectedOld},'new',{newPrediction},'metrics',metrics);
    save(fullfile(out,'joint_comparison.mat'),'results');
    plot_joint(results,out);
    disp(groupsummary(metrics,'model','mean',{'torqueRMSE','thetaRMSE','jointScore'}));
    fprintf('Outputs: %s\n',out);
end

function model = fit_surrogate(trials,w,zeta)
    eqX = []; eqY = []; group = [];
    for k = 1:numel(trials)
        t = trials(k); dt = 0.01;
        % Derivatives are permitted only during fitting. No test-angle
        % derivative or later test angle is passed to the rollout function.
        velocity = gradient(t.q,dt); acceleration = gradient(velocity,dt);
        eqX = [eqX; equilibrium_features(t)]; %#ok<AGROW>
        eqY = [eqY; t.q+2*zeta/w*velocity+acceleration/w^2]; %#ok<AGROW>
        group = [group; k*ones(size(t.q))]; %#ok<AGROW>
    end
    model.w = w; model.zeta = zeta;
    model.equilibrium = bse_ridge('fit',eqX,eqY,group,0.01);
    torqueX = []; torqueY = [];
    for k = 1:numel(trials)
        state = rollout_state(model,trials(k));
        torqueX = [torqueX; torque_features(trials(k),state)]; %#ok<AGROW>
        torqueY = [torqueY; trials(k).torque]; %#ok<AGROW>
    end
    model.torque = bse_ridge('fit',torqueX,torqueY,group,0.01);
end

function x = equilibrium_features(t)
    a = t.a;
    % No phase, time, condition ID, or reference trajectory is a predictor.
    x = [a,a.^2,a.^3,t.m.*a];
end

function state = rollout_state(model,t)
    desired = bse_ridge('predict',model.equilibrium,equilibrium_features(t));
    desired = min(max(desired,0),1.7); % Fixed equilibrium range, not fitted to test data.
    w = model.w; z = model.zeta; dt = 0.01;
    A = [0 1;-w*w -2*z*w]; B = [0;w*w];
    discrete = expm([A B;zeros(1,3)]*dt);
    state = zeros(numel(t.a),3);
    state(1,1) = t.q(1); % Only this initial angle is read; initial velocity is zero.
    for i = 1:numel(t.a)-1
        state(i,3) = -w*w*state(i,1)-2*z*w*state(i,2)+w*w*desired(i);
        state(i+1,1:2) = (discrete(1:2,1:2)*state(i,1:2)' + discrete(1:2,3)*desired(i))';
    end
    state(end,3) = -w*w*state(end,1)-2*z*w*state(end,2)+w*w*desired(end);
end

function x = torque_features(t,state)
    q = state(:,1); v = state(:,2); acceleration = state(:,3); a = t.a;
    % Learned observation equation; not an exact muscle-force balance law.
    x = [a,a.^2,a.^3,cos(q),t.m.*cos(q-0.936847),v,acceleration, ...
        a.*q,a.*v,t.m.*a];
end

function prediction = rollout_surrogate(model,t)
    state = rollout_state(model,t);
    torque = bse_ridge('predict',model.torque,torque_features(t,state));
    prediction = [torque,state(:,1)];
end

function metric = score_record(name,t,prediction)
    metric.model = name; metric.loadKg = t.loadKg; metric.speed = t.speed;
    metric.torqueRMSE = sqrt(mean((prediction(:,1)-t.torque).^2));
    metric.thetaRMSE = sqrt(mean((prediction(:,2)-t.q).^2));
    metric.thetaRMSEdeg = rad2deg(metric.thetaRMSE);
    % Fixed scales set before candidate selection, not test-range normalization.
    metric.jointScore = metric.torqueRMSE/100+metric.thetaRMSE/0.5;
end

function plot_joint(r,out)
    for group = 0:3
        if group == 0
            ids = [1 8 6]; % Preselected 5/35, 10/45, 15/55; all cases in other plots.
        else
            speedLevels = [35 45 55]; ids = find([r.trials.speed]==speedLevels(group));
        end
        fig = figure('Color','w','Position',[70 60 1350 850]);
        tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
        for row = 1:3
            k = ids(row); t = r.trials(k);
            for col = 1:2
                nexttile; hold on;
                if col==1, reference=t.torque; unit='Torque (Nm)'; else, reference=t.q; unit='Theta (rad)'; end
                plot(t.t,reference,'k','LineWidth',1.6);
                plot(t.t,r.archived{k}(:,col),'Color',[0.12 0.42 0.72],'LineWidth',1.3);
                plot(t.t,r.new{k}(:,col),'Color',[0.89 0.35 0.15],'LineWidth',1.3);
                xlabel('Time in segment (s)'); ylabel(unit); grid on;
                oldError=sqrt(mean((r.archived{k}(:,col)-reference).^2));
                newError=sqrt(mean((r.new{k}(:,col)-reference).^2));
                title(sprintf('%g kg / speed %g | RMSE old %.3f, new %.3f',t.loadKg,t.speed,oldError,newError));
                if row==1 && col==1, legend('Reference','Selected archived BSE','New reduced state model','Location','best'); end
            end
        end
        if group == 0
            sgtitle('Torque AND theta | old archived BSE vs new reduced state model | conditional comparison');
            name='torque_theta_comparison';
        else
            sgtitle(sprintf('Torque AND theta | withheld speed %g | conditional MAP-input comparison',r.trials(ids(1)).speed));
            name=sprintf('torque_theta_speed_%g',r.trials(ids(1)).speed);
        end
        exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',150);
        exportgraphics(fig,fullfile(out,[name '.pdf']),'ContentType','vector');
    end
end

function cfg = get_config(root,id)
    cfg=jsondecode(fileread(fullfile(root,'config','workflows', ...
        ['training_map__muscle_tendon__' id '.json'])));
end

function filename = asset(root,cfg,name)
    k=find(strcmp({cfg.files.target},name)); assert(isscalar(k));
    filename=fullfile(root,cfg.files(k).source);
end

