function test_no_posture_feasibility
%TEST_NO_POSTURE_FEASIBILITY Compare eliminated slacks with the explicit QP.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'Demo','NMPC'));
rng(528);n=4;points=8;
fast=optimoptions('quadprog','Algorithm','active-set','Display','off', ...
    'OptimalityTolerance',1e-9,'ConstraintTolerance',1e-9,'MaxIterations',150);
reference=optimoptions('quadprog','Algorithm','interior-point-convex','Display','off', ...
    'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10,'MaxIterations',200);
errors=zeros(12,1);
for trial=1:12
    A=randn(3*points,n);b=.5*randn(3*points,1);lo=-2*ones(n,1);hi=-lo;
    [u,flag]=no_posture_feasibility(A,b,zeros(n,1),lo,hi,points,fast);
    assert(flag>0);
    limited=fast;limited.MaxIterations=1;
    [recovered,recoveredFlag]=no_posture_feasibility(A,b,zeros(n,1),lo,hi,points,limited);
    assert(all(isfinite(recovered)) && all(recovered>=lo & recovered<=hi));
    weights=[ones(points,1);100*ones(2*points,1)];
    [x,~,status]=quadprog(diag([2e-8*ones(n,1);2*weights]), ...
        zeros(n+3*points,1),[A,-eye(3*points)],b,[],[],[lo;zeros(3*points,1)], ...
        [hi;inf(3*points,1)],[],reference);
    assert(status>0);
    r=A*u-b;
    actual=sum(weights.*max(r,0).^2)+1e-8*(u'*u);
    target=sum(weights.*x(n+1:end).^2)+1e-8*sum(x(1:n).^2);
    errors(trial)=abs(actual-target)/max(1,target);
    recoveredMerit=sum(weights.*max(A*recovered-b,0).^2)+1e-8*(recovered'*recovered);
    assert(abs(recoveredMerit-target)/max(1,target)<1e-4 && recoveredFlag>0, ...
        'The limited-iteration active-set solve did not recover through fallback.');
end
disp(errors');
assert(max(errors)<1e-4,'Eliminated-slack solver differs from the explicit convex problem.');
fprintf('PASS: pointwise slacks match 12 explicit QP references, including active-set iteration-limit recovery.\n');
end
