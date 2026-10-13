function test_bse_angle_priority()
%TEST_BSE_ANGLE_PRIORITY Check exact integration, boundedness and split outputs.
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'experiments','bse_angle_priority'));
    d=load(fullfile(root,'results','bse_angle_priority','angle_priority.mat'));
    r=d.result; rng(104);
    for source=1:2
        model=r.models{source,1}; dt=.01;
        a=repelem(rand(600,1),100); mass=15*double(a>.5);
        p=angle_rollout(model,a,mass,15,.8,dt);
        assert(all(isfinite(p.q)) && min(p.q)>=0 && max(p.q)<=model.limit);
        % Compare transfer-function implementation with direct state iteration.
        z=exp(-dt/model.tau); A=[z 0;z*dt/model.tau z]; B=1-sum(A,2);
        x=[.8;.8]; q=zeros(size(a));
        for k=1:numel(a), q(k)=x(2); x=A*x+B*p.equilibrium(k); end
        assert(max(abs(q-p.q))<1e-10,'Exact integration mismatch.');
        assert(max(abs(eig(A)))<1 && all(B>=0));
        for mode=1:2
            for k=1:9
                t=r.trials(k); saved=r.predictions{source,mode,k};
                if mode==1, m=model; else, speeds=[35 45 55]; m=r.models{source,find(speeds==t.speed)+1}; end
                q0=t.q(1); t.q(:)=NaN; t.torque(:)=NaN;
                actual=angle_rollout(m,t.a(:,source),t.m,t.loadKg,q0,.01);
                assert(isequal(actual,saved),'Prediction changed without reference labels.');
            end
        end
    end
    assert(height(r.metrics)==36 && all(isfinite(r.metrics.thetaRMSE)));
    fprintf('PASS: paired-source rollouts, 600 s bounded stress, exact integration, and reference isolation.\n');
end
