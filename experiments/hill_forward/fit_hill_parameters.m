function model=fit_hill_parameters(trials,source,seedParameters)
    % Only F0*R and R/lo are identifiable in this one-joint rigid-tendon
    % geometry. Fix R to the archived value and optimize those combinations.
    lower=[40 .06 .75 0 1 .5]; upper=[550 1.1 1.3 8 10 2];
    if nargin<3
    rng(114+source); candidates=rand(180,6);
    candidates(1,:) = ([180 .7 .95 4 3 1]-lower)./(upper-lower);
    screening=zeros(180,1);
    for j=1:180, e=loss(candidates(j,:)); screening(j)=sum(e.^2); end
    [~,order]=sort(screening); starts=candidates(order(1:6),:);
    else
        starts=zeros(size(seedParameters,1),6);
        for i=1:size(seedParameters,1)
            p=seedParameters(i,:); x=[p(1)*p(3),p(3)/p(2),p(4:7)];
            starts(i,:)=min(1,max(0,(x-lower)./(upper-lower)));
        end
    end
    options=optimoptions('lsqnonlin','Display','off','MaxIterations',100, ...
        'MaxFunctionEvaluations',1600,'FiniteDifferenceStepSize',1e-5,'FunctionTolerance',2e-5);
    best=Inf;
    for j=1:size(starts,1)
        initial=starts(j,:);
        [z,value,~,flag]=lsqnonlin(@loss,initial,zeros(1,6),ones(1,6),options);
        fprintf('  start %d loss %.4f exit %d\n',j,value,flag);
        if value<best, best=value; model=struct('parameters',decode(z),'loss',value,'exitflag',flag); end
    end
    function par=decode(z)
        x=lower+z.*(upper-lower); R=.0328531574526041;
        par=[x(1)/R,R/x(2),R,x(3:6)];
    end
    function residual=loss(z)
        par=decode(z); residual=[];
        for t=trials
            p=hill_rollout(par,t.a(:,source),t.m,[t.q(1) t.v0],.01,2); n=sqrt(numel(t.q));
            % Fit reference only here. Forward prediction never reads labels.
            e=[(p.q-t.q)/.2; .35*(p.torque-t.torque)/70; ...
                .2*(p.v-gradient(t.q,.01))/.6; ...
                5*(max(0,-p.q)+max(0,p.q-95*pi/180)); ...
                .08*diff(p.torque)/(.01*300);0; ...
                .3*p.passiveTorque/100]/n;
            [~,lo]=findpeaks(-t.q,'MinPeakProminence',.12,'MinPeakDistance',50);
            [~,hi]=findpeaks(t.q,'MinPeakProminence',.12,'MinPeakDistance',50);
            ix=unique([1;lo;hi;numel(t.q)]);
            e=[e;.5*(p.q(ix)-t.q(ix))/.2/sqrt(numel(ix))];
            e(~isfinite(e))=1e6;
            residual=[residual;e/sqrt(numel(trials))]; %#ok<AGROW>
        end
    end
end
