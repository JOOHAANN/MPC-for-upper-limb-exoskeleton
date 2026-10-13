function X = angle_features(a,mass,nominalMass,dt)
%ANGLE_FEATURES Causal activation memory and known external-load context.
% No elapsed phase, speed label, or reference angle is a model feature.
    a=memory(a(:),.04,dt); fast=memory(a,.10,dt); slow=memory(a,.40,dt);
    % Smooth contact context causally to avoid an artificial command jump.
    d=a-fast; lag=fast-slow; m=memory(mass(:)/15,.20,dt); n=nominalMass/15;
    X=[a,a.^2,a.^3,fast,slow,d,lag,a.*d,a.*lag, ...
       n*a,n*a.^2,n*fast,m,m.*a,m.*fast];
end
function y=memory(a,tau,dt)
    r=exp(-dt/tau);
    y=filter(1-r,[1 -r],a,r*a(1));
end
