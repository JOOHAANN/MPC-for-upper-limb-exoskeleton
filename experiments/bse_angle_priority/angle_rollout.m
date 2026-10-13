function p = angle_rollout(model, a, mass, nominalMass, q0, dt)
%ANGLE_ROLLOUT Causal, bounded, critically damped equilibrium dynamics.
% No reference trajectory is accepted: only the initial angle is required.
% Both activation sources use this exact function and feature definition.
    assert(isfinite(q0) && q0>=0 && q0<=model.limit);
    X = angle_features(a,mass,nominalMass,dt);
    z = [ones(size(a)), (X-model.mu)./model.scale]*model.beta;
    u = model.limit./(1+exp(-max(-40,min(40,z))));
    h = model.tau; r = exp(-dt/h);
    A = [r 0; r*dt/h r]; B = [1-r;1-r-r*dt/h];
    % Exact zero-order-hold integration, using a positive state transition.
    % B>=0 and each row of [A B] sums to one, so the angle interval is invariant.
    s = q0 + filter([0 B(1)],[1 -r],u-q0);
    q = q0 + filter([0 B(2) A(2,1)*B(1)-r*B(2)], ...
        [1 -2*r r*r],u-q0);
    v = (s-q)/h; acceleration = ((u-s)/h-v)/h;
    % Newton-Euler torque for the SAME predicted angle, not a second regressor.
    M=22.6006; L1=.345; L2=.760606; phi=.936847;
    inertia = mass*L2^2+4*M*L1^2/3;
    torque = inertia.*acceleration+model.damping*v ...
        +9.81*M*L1*cos(q)+9.81*mass*L2.*cos(q-phi);
    p=struct('q',q,'v',v,'acceleration',acceleration,'torque',torque,'equilibrium',u);
end
