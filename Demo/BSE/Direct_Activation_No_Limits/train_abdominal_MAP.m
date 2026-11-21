function model = train_abdominal_MAP(records)
%TRAIN_ABDOMINAL_MAP Repeat the abdominal MAP fit from 25 recorded profiles.
% This function does not use IK, torque, BSE outputs or back-to-abdomen terms.
% It writes a candidate model. Repackage inputs and refit BSE before publishing.
folder=fileparts(mfilename('fullpath'));
if nargin<1
    data=load(fullfile(folder,'data','legacy_MAP_528.mat'));records=data.records;
end
phase=linspace(0,1,200)';rng(528);
n=numel(records);Y=zeros(200,n);loads=zeros(n,1);speeds=loads;
for k=1:n
    r=records{k};Y(:,k)=r.profile(:);
    loads(k)=r.loadKg;speeds(k)=r.speedCode;
end
amplitude=max(prctile(Y(:),99),.01);Y=Y/amplitude;

%% Four nonnegative temporal patterns, with one shared amplitude scale
W=.1+.7*rand(200,4);H=.01+.29*rand(4,n);
for k=1:1200
    H=H.*((W'*Y)./max(W'*W*H+.0001*H,1e-12));
    W=W.*((Y*H')./max(W*H*H'+.0001*W,1e-12));
    scale=max(max(W,[],1),1e-10);W=W./scale;H=H.*scale';
end

%% Three Gaussian lobes describe each temporal pattern
options=optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',5000);
primitives=zeros(200,4);gaussianParameters=zeros(4,9);
for k=1:4
    bestCost=inf;
    for centers=[.2 .5 .8;.3 .45 .7;.15 .65 .85]'
        guess=reshape([.6*ones(1,3);centers';.15*ones(1,3)],1,9);
        [z,cost]=lsqnonlin(@(z)gaussian_curve(z,phase)-W(:,k),guess, ...
            repmat([-2 0 .025],1,3),repmat([3 1 .6],1,3),options);
        if cost<bestCost,best=z;bestCost=cost;end
    end
    gaussianParameters(k,:)=best;
    primitives(:,k)=max(0,gaussian_curve(best,phase));
end
weights=zeros(n,4);
for k=1:n,weights(k,:)=lsqnonneg(primitives,Y(:,k))';end

%% Original additive load/speed increments, referenced to 10 kg / code 45
speedIncrement=zeros(n,4);loadIncrement=zeros(n,4);
for k=1:n
    speedReference=find(loads==loads(k) & speeds==45,1);
    loadReference=find(loads==10 & speeds==speeds(k),1);
    speedIncrement(k,:)=weights(k,:)-weights(speedReference,:);
    loadIncrement(k,:)=weights(k,:)-weights(loadReference,:);
end
selected=ismember(loads,[5 10 15]);
loadCoefficients=zeros(3,4);speedCoefficients=zeros(3,4);
for k=1:4
    loadCoefficients(:,k)=polyfit(loads(selected),loadIncrement(selected,k),2)';
    speedCoefficients(:,k)=polyfit(speeds,speedIncrement(:,k),2)';
end
baseline=weights(loads==10 & speeds==45,:)';
model=struct('primitives',primitives,'gaussianParameters',gaussianParameters, ...
    'baseline',baseline,'loadCoefficients',loadCoefficients, ...
    'speedCoefficients',speedCoefficients,'phase',phase,'amplitudeScale',amplitude, ...
    'targetMuscle','rect_abd_r','architecture', ...
    '4 abdominal NMF temporal primitives, gauss3 fitting, original additive quadratic load/speed increments with 10kg/code45 baseline', ...
    'usesIK',false,'usesTorque',false,'usesBSEOutput',false,'usesBackToAbdomenRecruitment',false);
output=fullfile(folder,'output','abdominal_MAP_retraining');
if ~isfolder(output),mkdir(output);end
save(fullfile(output,'candidate_model.mat'),'model','-v7');
fprintf('Saved abdominal MAP candidate to %s\n',output);
end

function y=gaussian_curve(z,phase)
z=reshape(z,3,3)';y=zeros(size(phase));
for k=1:3,y=y+z(k,1)*exp(-((phase-z(k,2))/z(k,3)).^2);end
end
