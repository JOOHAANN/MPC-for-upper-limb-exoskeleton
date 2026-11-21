% Diagnostic only: assess the actuator ceiling around the remaining peak.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
folder=fullfile(root,'Demo','NMPC');bse=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
addpath(folder,bse);d=load(fullfile(bse,'data','subject_528.mat'));s=d.subject;
c=load(fullfile(folder,'output','current_528_Exp_dense_g3','group3_source1.mat'));
q=s.sequences{3};settings=current_528_settings(s);model=s.models{1};
plant=current_bse_plant(model,s.balanceController,q.standingAngle,settings.humanPD,settings.actuator);
inputs=prepare_activation(q,model,1);n=numel(q.time);rows=[];
for onset=[37 38 39 40 41]
    [state,first]=plant.initialize(inputs(1,:),q.mass(1),q.waiting(1),q.initialState);
    out=zeros(n,27);out(1,:)=first;
    for k=2:n
        command=c.pair.nmpc.actual(k,25);
        if q.time(k)>=onset && q.time(k)<=43.5,command=80;end
        [state,out(k,:)]=plant.step(state,inputs(k,:),q.mass(k),q.waiting(k),command,q.dt,10);
    end
    ix=q.time>=37 & q.time<=43.5;
    rows=[rows;onset max(max(out(ix,4)-c.pair.threshold(ix),0)) min(out(ix,1)) max(out(ix,1))]; %#ok<AGROW>
end
disp(array2table(rows,'VariableNames',{'SaturationOnsetSeconds','PeakExcessNm','MinimumThetaRad','MaximumThetaRad'}));
save(fullfile(fileparts(mfilename('fullpath')),'peak_feasibility.mat'),'rows');
