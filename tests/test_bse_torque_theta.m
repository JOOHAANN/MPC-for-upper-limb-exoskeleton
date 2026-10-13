function test_bse_torque_theta()
%TEST_BSE_TORQUE_THETA Verify saved state-rollout outputs and held-out groups.
    root = fileparts(fileparts(mfilename('fullpath')));
    saved = load(fullfile(root,'results','bse_torque_theta','joint_comparison.mat'),'results');
    r = saved.results;
    coverage = zeros(1,9);
    for fold = r.folds
        assert(isempty(intersect(fold.trainTrials,fold.testTrials)));
        assert(isempty(intersect([r.trials(fold.trainTrials).speed], ...
            [r.trials(fold.testTrials).speed])));
        coverage(fold.testTrials) = coverage(fold.testTrials)+1;
        poles = eig([0 1;-fold.omega^2 -2*fold.zeta*fold.omega]);
        assert(all(real(poles)<0));
    end
    assert(all(coverage==1));
    for k = 1:9
        prediction = r.new{k};
        assert(isequal(size(prediction),[numel(r.trials(k).q),2]));
        assert(all(isfinite(prediction),'all'));
        assert(prediction(1,2)==r.trials(k).q(1));
        mask = strcmp(r.metrics.model,'New reduced state model') & ...
            r.metrics.loadKg==r.trials(k).loadKg & r.metrics.speed==r.trials(k).speed;
        expected = sqrt(mean((prediction(:,2)-r.trials(k).q).^2));
        assert(abs(expected-r.metrics.thetaRMSE(mask))<1e-12);
    end
    fprintf('PASS: held-out groups, stable state poles, finite two-output rollouts, metrics.\n');
end
