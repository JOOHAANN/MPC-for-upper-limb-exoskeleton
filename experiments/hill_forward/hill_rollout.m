function p=hill_rollout(parameters,a,m,q0,dt,substeps,useMex)
%HILL_ROLLOUT Rigid-tendon Hill muscle drives the joint through its own torque.
% Parameters: F0, optimal fiber length, moment arm, normalized fiber length
% at q=1.4, joint viscosity, maximum normalized shortening speed, input power.
% No reference torque/trajectory, equilibrium-angle predictor, or wall force.
    if nargin<6, substeps=2; end
    if nargin<7, useMex=true; end
    assert(numel(parameters)==7 && all(isfinite(parameters)) && all(parameters([1 2 3 6 7])>0));
    assert(iscolumn(a) && iscolumn(m) && numel(a)==numel(m) && numel(a)>=2);
    assert(all(isfinite(a)) && all(isfinite(m)) && all(m>=0));
    assert(all(isfinite(q0)) && ismember(numel(q0),[1 2]) && dt>0 && isfinite(dt));
    assert(substeps>=1 && substeps==floor(substeps));
    n=numel(a); h=dt/substeps; q=q0(1); v=0;
    if numel(q0)==2, v=q0(2); end
    state=zeros(n,2); state(1,:)=[q v]; failed=false;
    if useMex && exist('hill_integrate_mex','file')==3
        [state,failed]=hill_integrate_mex(parameters,a,m,q0,dt,substeps);
    else
    for k=1:n-1
        for step=1:substeps
            s=(step-1)/substeps;
            aa=a(k)+(a(k+1)-a(k))*s; mm=m(k)+(m(k+1)-m(k))*s;
            ab=a(k)+(a(k+1)-a(k))*(s+.5/substeps);
            mb=m(k)+(m(k+1)-m(k))*(s+.5/substeps);
            ac=a(k)+(a(k+1)-a(k))*(s+1/substeps);
            mc=m(k)+(m(k+1)-m(k))*(s+1/substeps);
            [~,acc1]=forces(parameters,q,v,aa,mm);
            v2=v+h*acc1/2; [~,acc2]=forces(parameters,q+h*v/2,v2,ab,mb);
            v3=v+h*acc2/2; [~,acc3]=forces(parameters,q+h*v2/2,v3,ab,mb);
            v4=v+h*acc3; [~,acc4]=forces(parameters,q+h*v3,v4,ac,mc);
            q=q+h*(v+2*v2+2*v3+v4)/6;
            v=v+h*(acc1+2*acc2+2*acc3+acc4)/6;
            % Reject numerical/geometric failure; NEVER clip a physical angle.
            if ~isfinite(q+v) || abs(q)>4 || abs(v)>50
                failed=true; break;
            end
        end
        if failed, state(k+1:end,:)=repmat([4 0],n-k,1); break; end
        state(k+1,:)=[q v];
    end
    end
    q=state(:,1); v=state(:,2);
    [torque,acceleration,active,passive,fiberLength,fiberVelocity,gravity]=forces(parameters,q,v,a,m);
    p=struct('q',q,'v',v,'torque',torque,'acceleration',acceleration, ...
        'activeTorque',active,'passiveTorque',passive,'fiberLength',fiberLength, ...
        'fiberVelocity',fiberVelocity,'gravityTorque',gravity, ...
        'dampingTorque',-parameters(5)*v,'boundaryTorque',zeros(size(q)),'failed',failed);
end
function [torque,acc,active,passive,lm,vm,gravity]=forces(p,q,v,a,m)
    F0=p(1); lo=p(2); R=p(3); lref=p(4); B=p(5); vmax=p(6);
    alpha0=.24086; width=lo*sin(alpha0);
    % Constant-height pennation and energy-consistent constant moment arm:
    % dLmt/dq=-R, so the muscle's generalized torque is +R*Ft.
    projected=lo*sqrt(lref^2-sin(alpha0)^2)-R*(q-1.4);
    lm=sqrt(projected.^2+width^2); cosine=max(projected,0)./lm;
    vm=-R*v.*cosine; ln=lm/lo; vn=vm/(vmax*lo);
    fl=exp(-((ln-1)/.45).^2);
    fv=zeros(size(vn)); shorten=vn<=0;
    fv(shorten)=max(0,1+vn(shorten))./(1-vn(shorten)/.25);
    fv(~shorten)=(1+1.5*vn(~shorten)/.25)./(1+vn(~shorten)/.25);
    % Tensile parallel elastic element. Exponent cap is numerical overflow
    % protection only; geometric-domain violations are reported separately.
    fp=max(0,(exp(min(40,5*(ln-1)/.6))-1)/(exp(5)-1));
    active=R*F0*max(0,min(1,a)).^p(7).*fl.*fv.*cosine;
    passive=R*F0*fp.*cosine; torque=active+passive;
    gravity=-9.81*22.6006*.345*cos(q)-9.81*m*.760606.*cos(q-.936847);
    inertia=m*.760606^2+4*22.6006*.345^2/3;
    acc=(torque+gravity-B*v)./inertia;
end
