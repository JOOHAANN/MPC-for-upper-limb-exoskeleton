function build_528_inputs
%BUILD_528_INPUTS Reproduce legacy selection and calibrate MAP by load only.
% Raw recordings and the previous torque/theta references are read-only.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(fileparts(folder)));
project=fileparts(root);legacy=fullfile(project,'YOUHAN LI code and project description','My_MAP_model - 5 10 15');
raw=fullfile(project,'YOUHAN LI code and project description','Data','240528Pilot');
support=fullfile(root,'experiments','load_scaled_528');
data=load(fullfile(support,'before_load_rescale','data','subject_528.mat'));
subject=data.subject;phase=linspace(0,1,200)';
loads=[5 7.5 10 12.5 15];tokens=[5 75 10 125 15];speeds=[25 35 45 55 65];
% Exact row selections from the original NNMF_component.m, inclusive.
S=[11797 17280;10364 14324;2808 5021;4620 6669;3486 5154; ...
    4952 9416;5727 8919;3453 6377;4075 6350;3523 5447; ...
    5420 9558;3864 7572;5330 8588;2368 3831;3723 5388; ...
    3515 7448;2416 5825;3825 7001;3945 6189;2224 3532; ...
    5774 9452;4358 7197;4552 7392;1878 3470;4468 6102];
records=cell(1,25);trials=cell(5,5);k=0;
fprintf('Reading legacy experimental selections for all 25 conditions.\n');
for i=1:5
    for j=1:5
        k=k+1;path=fullfile(raw,sprintf('dynamics_stoop_%dkg_%d',tokens(i),speeds(j)));
        trial=experimental_trial(path);trials{i,j}=trial;
        emg=importdata(fullfile(path,'emgFilt.sto'));
        lo=emg.data(S(k,1),1);hi=emg.data(S(k,2),1);
        % Average curves on normalized action phase, not on absolute time.
        query=lo+phase*(hi-lo);
        time=(trial.time(1):.001:trial.time(end))';
        signals=interp1(trial.time,trial.activation,time,'linear');
        signals(:,1)=smoothdata(signals(:,1),'gaussian',5);
        signals(:,2)=smoothdata(signals(:,2),'gaussian',5);
        signals=legacy_smoothing(signals);
        profile=interp1(time,signals,query,'linear');
        originalWindow=[lo hi];selection='Original NNMF row window';
        if any(~isfinite(profile),'all') && loads(i)==5
            if speeds(j)==25
                % The archived second pair extends beyond CEINMS output.
                % Use the complete earlier pair visible in window_audit.png.
                lo=emg.data(1,1)+7.4;hi=emg.data(1,1)+14.8;
            elseif speeds(j)==35
                original=subject.sequences{1};seg=original.segments{1};
                lo=trial.time(1)+original.time(seg.startIndex);
                hi=trial.time(1)+original.time(seg.endIndex);
            end
            query=lo+phase*(hi-lo);profile=interp1(time,signals,query,'linear');
            selection='Complete earlier M pair; archived row window exceeds CEINMS support';
            fprintf('  5 kg / code %d: selected %.3f--%.3f s from EMG onset.\n', ...
                speeds(j),lo-emg.data(1,1),hi-emg.data(1,1));
        end
        assert(all(isfinite(profile),'all'),'Activation support is incomplete: %s',path);
        records{k}=struct('loadKg',loads(i),'speedCode',speeds(j), ...
            'phase',phase,'experimental',profile,'profile',profile(:,2), ...
            'sourcePath',path,'emgRows',S(k,:),'startEpoch',lo,'endEpoch',hi, ...
            'originalWindow',originalWindow,'selection',selection);
    end
end

%% Retain the original trained back predictor and original fusion formula
x=load(fullfile(legacy,'Primitives_fit_y.mat'));v=load(fullfile(legacy,'pv_coefficients.mat'));
w=load(fullfile(legacy,'pw_coefficients.mat'));b=load(fullfile(legacy,'baseline_data.mat'));
map.back=struct('primitives',x.Primitives_fit_y,'baseline',b.baseline_data, ...
    'loadCoefficients',w.pw_coefficients,'speedCoefficients',v.pv_coefficients);
map.abdomen=train_abdominal_MAP(records);
map.phase=phase;map.backFusion='(MEP1+MEP2)/12 + (MEP3+MEP4)/4 + (MEP1+MEP2)/6';
map.backSource=legacy;map.usesBSEOutput=false;
map.abdominalTrainingUsesIK=false;map.abdominalTrainingUsesTorque=false;

%% One amplitude calibration per load, shared by every speed
scaling.loadKg=[5;10;15];scaling.speedCodes=speeds;
scaling.gain=zeros(3,2);scaling.offset=zeros(3,2);
scaling.experimentalMinimum=zeros(3,2);scaling.experimentalMaximum=zeros(3,2);
scaling.rawMAPMinimum=zeros(3,2);scaling.rawMAPMaximum=zeros(3,2);
scaling.meanExperimental=zeros(200,2,3);scaling.meanRawMAP=zeros(200,2,3);
for i=1:3
    selected=find(cellfun(@(r)r.loadKg==scaling.loadKg(i),records));
    experimental=zeros(200,2,5);prediction=experimental;
    for j=1:5
        r=records{selected(j)};experimental(:,:,j)=r.experimental;
        prediction(:,:,j)=predict_528_MAP(map,r.loadKg,r.speedCode,phase);
    end
    e=mean(experimental,3);m=mean(prediction,3);
    emin=min(e);emax=max(e);mmin=min(m);mmax=max(m);
    assert(all(mmax-mmin>1e-8),'MAP has no usable dynamic range.');
    scaling.gain(i,:)=(emax-emin)./(mmax-mmin);
    scaling.offset(i,:)=emin-scaling.gain(i,:).*mmin;
    scaling.experimentalMinimum(i,:)=emin;scaling.experimentalMaximum(i,:)=emax;
    scaling.rawMAPMinimum(i,:)=mmin;scaling.rawMAPMaximum(i,:)=mmax;
    scaling.meanExperimental(:,:,i)=e;scaling.meanRawMAP(:,:,i)=m;
end

%% Rebuild activation inputs on the unchanged three recorded timelines
for group=1:numel(subject.sequences)
    q=subject.sequences{group};fineTime=(0:.001:q.time(end))';
    n=numel(fineTime);experimental=zeros(n,2);offset=0;
    for stage=1:3
        seg=q.segments{stage};trial=trials{find(loads==seg.loadKg),find(speeds==seg.speedCode)};
        first=1;last=300;if stage==2,last=350;end;if stage==3,first=50;end
        count=floor((trial.time(last)-trial.time(first))/.001+.0003)+1;
        index=offset+(1:count);index=index(index<=n);
        epoch=trial.time(first)+(index'-offset-1)*.001;
        value=interp1(trial.time,trial.activation,epoch,'linear');
        value=smoothdata(value,'gaussian',5);
        experimental(index,:)=value;offset=offset+count;
    end
    waiting=interp1(q.time,double(q.waiting),fineTime,'nearest','extrap')>.5;
    experimental(waiting,:)=0;
    experimental=legacy_smoothing(experimental);
    input=interp1(fineTime,experimental,q.time,'linear');
    q.backActivation(:,1)=input(:,1);q.abdomenActivation=input(:,2);
    q.backActivation(:,2)=0;q.abdomenMAPActivation=zeros(size(q.time));
    q.rawMAPActivation=zeros(numel(q.time),2);
    for stage=1:3
        seg=q.segments{stage};ix=seg.startIndex:seg.endIndex;
        localPhase=(q.time(ix)-q.time(ix(1)))/(q.time(ix(end))-q.time(ix(1)));
        rawProfile=predict_528_MAP(map,seg.loadKg,seg.speedCode,localPhase);
        calibrated=rescale_MAP_by_load(rawProfile,seg.loadKg,scaling);
        q.rawMAPActivation(ix,:)=rawProfile;
        q.backActivation(ix,2)=calibrated(:,1);q.abdomenMAPActivation(ix)=calibrated(:,2);
    end
    subject.sequences{group}=q;
end
for source=1:2
    subject.models{source}.abdomenScale=1;
    subject.models{source}.inputPolicy='Legacy experimental processing; MAP amplitude calibrated by load across all five speeds.';
end
subject=rmfield(subject,'abdominalMAP');subject.MAP=map;subject.MAPScaling=scaling;
subject.modelVersion='legacy_selection_load_scaled_MAP_dual_MTU';
subject.parameterPolicy='One fixed 24-parameter physical vector per source across nine bouts; data-derived MAP scaling depends on load only.';
subject.calibrationProtocol.numberOfInputGains=0;
subject.calibrationProtocol.inputPolicy='Legacy selection and preprocessing; shared experimental-mean peak/valley calibration by load.';
subject.calibrationProtocol.status='input_rebuilt_pending_physical_fit';
save(fullfile(folder,'data','subject_528.mat'),'subject','-v7');
save(fullfile(folder,'data','legacy_MAP_528.mat'),'map','scaling','records','-v7');
save(fullfile(support,'calibration_input.mat'),'subject','-v7');
save(fullfile(support,'input_audit.mat'),'records','scaling','map','-v7');
fprintf('Built 528 inputs; preserved original MAP fusion, nine reference curves and standing schedule.\n');
disp(table(scaling.loadKg,scaling.gain(:,1),scaling.offset(:,1), ...
    scaling.gain(:,2),scaling.offset(:,2),'VariableNames', ...
    {'LoadKg','BackGain','BackOffset','AbdomenGain','AbdomenOffset'}));
end

function trial=experimental_trial(path)
% Dynamic_kine_reference: representative activations weighted by group moments.
lmt=importdata(fullfile(path,'lmt.sto'));a=importdata(fullfile(path,'Activations.sto'));
f=importdata(fullfile(path,'MuscleForces.sto'));r=importdata(fullfile(path,'ma_L5_S1_Flex_Ext.sto'));
a=unique_clock(a);f=unique_clock(f);r=unique_clock(r);
time=lmt.data(:,1);representatives={'IL_L1_r','IL_L1_l','LTpT_T7_r','LTpT_T7_l','LTpL_L1_r','LTpL_L1_l'};
moments=zeros(numel(time),6);activation=moments;
first=find(a.data(:,1)<=time(5),1,'last');last=find(a.data(:,1)<=time(end),1,'last');
for group=1:6
    side='r';if mod(group,2)==0,side='l';end
    if group<=2
        names=muscle_names('IL_L',1:4,side);
    elseif group<=4
        names=[muscle_names('LTpT_T',7:12,side),muscle_names('LTpT_R',7:12,side)];
    else
        names=[muscle_names('IL_R',5:12,side),muscle_names('LTpL_L',5:-1:1,side)];
    end
    for k=1:numel(names)
        fi=find(strcmp(f.colheaders,names{k}),1);ri=find(strcmp(r.colheaders,names{k}),1);
        assert(~isempty(fi) && ~isempty(ri),'A legacy muscle group is missing.');
        moments(:,group)=moments(:,group)+interp1(f.data(:,1),f.data(:,fi),time) ...
            .*interp1(r.data(:,1),r.data(:,ri),time);
    end
    ai=find(strcmp(a.colheaders,representatives{group}),1);
    activation(:,group)=interp1(a.data(first:last,1),a.data(first:last,ai),time);
end
back=sum(activation.*(moments./sum(moments,2)),2);back(~isfinite(back))=.015;
ai=find(strcmp(a.colheaders,'rect_abd_r'),1);
abdomen=interp1(a.data(first:last,1),a.data(first:last,ai),time);abdomen(~isfinite(abdomen))=.008;
trial=struct('time',time,'activation',[back abdomen]);
end

function signals=legacy_smoothing(signals)
signals(:,1)=smoothdata(signals(:,1),'gaussian',300);
[b,a]=butter(2,(1/1.2)/(1000/2),'low');signals(:,2)=filter(b,a,signals(:,2));
end

function names=muscle_names(prefix,numbers,side)
names=arrayfun(@(n)sprintf('%s%d_%s',prefix,n,side),numbers,'UniformOutput',false);
end

function table=unique_clock(table)
% Some exports repeat their last row exactly; retain one copy, never average.
[~,index]=unique(table.data(:,1),'stable');
if numel(index)<size(table.data,1)
    duplicate=setdiff(1:size(table.data,1),index);
    for k=duplicate
        first=find(table.data(:,1)==table.data(k,1),1);
        assert(max(abs(table.data(first,:)-table.data(k,:)))<1e-10, ...
            'Conflicting rows share a timestamp; inspect the export before fitting.');
    end
    table.data=table.data(index,:);
end
end
