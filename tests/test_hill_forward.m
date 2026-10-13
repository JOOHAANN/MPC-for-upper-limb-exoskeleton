function test_hill_forward()
%TEST_HILL_FORWARD Verify one force drives both outputs and integration converges.
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'experiments','hill_forward'));
    d=load(fullfile(root,'results','hill_forward','hill_forward.mat')); r=d.result;
    convergence=struct([]); speeds=[35 45 55];
    for source=1:2
        for k=1:9
            t=r.trials(k); f=find(speeds==t.speed)+1; model=r.models{source,f}; par=model.parameters;
            assert(~ismember(k,model.trainIds));
            assert(all([r.trials(model.trainIds).speed]~=t.speed));
            p=r.predictions{source,2,k};
            assert(~p.failed && all(isfinite(p.q)) && all(p.torque>=0));
            I=t.m*.760606^2+4*22.6006*.345^2/3;
            assert(max(abs(p.torque-p.activeTorque-p.passiveTorque))<1e-9);
            assert(max(abs(I.*p.acceleration-p.torque-p.gravityTorque-p.dampingTorque))<1e-8);
            assert(all(p.boundaryTorque==0));
            % The prediction interface only accepts inputs and the initial angle.
            q0=[t.q(1) t.v0]; t.q(:)=NaN; t.torque(:)=NaN;
            again=hill_rollout(par,t.a(:,source),t.m,q0,.01,2);
            assert(isequaln(again,p));
            fine=hill_rollout(par,t.a(:,source),t.m,q0,.01,4);
            rec=struct('source',source,'trial',k,'maxAngleDifference',max(abs(p.q-fine.q)), ...
                'maxTorqueDifference',max(abs(p.torque-fine.torque)));
            assert(rec.maxAngleDifference<.002 && rec.maxTorqueDifference<1);
            if isempty(convergence), convergence=rec; else, convergence(end+1)=rec; end %#ok<AGROW>
        end
        % Independently compare the optional compiled integrator to readable MATLAB.
        t=r.trials(1); par=r.models{source,2}.parameters;
        matlab=hill_rollout(par,t.a(:,source),t.m,[t.q(1) t.v0],.01,2,false);
        compiled=hill_rollout(par,t.a(:,source),t.m,[t.q(1) t.v0],.01,2,true);
        assert(max(abs(matlab.q-compiled.q))<1e-9);
        assert(max(abs(matlab.torque-compiled.torque))<1e-7);
    end
    writetable(struct2table(convergence),fullfile(root,'results','hill_forward','integration_convergence.csv'));
    fprintf('PASS: Hill torque identity, Newton-Euler balance, zero wall, reference isolation, speed splits, RK4 refinement, MEX parity.\n');
end
