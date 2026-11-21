% DIAGNOSTIC ONLY: future experimental activation gives an unattainable oracle.
% This run is not part of the deployable controller or reported comparison.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
folder=fullfile(root,'Demo','NMPC');bse=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');addpath(folder,bse);
d=load(fullfile(bse,'data','subject_528.mat'));subject=d.subject;
q=subject.sequences{3};model=subject.models{1};s=current_528_settings(subject);o=s.controllers{1};
o.activationResidualPredictor.coefficients(:)=0;
plant=current_bse_plant(model,subject.balanceController,q.standingAngle,s.humanPD,s.actuator);
inputs=prepare_activation(q,model,1);baseline=plant.simulate(inputs,q.mass,q.waiting,q.initialState,q.dt,10);
threshold=zeros(size(q.time));limits=zeros(1,3);
for k=1:3
    z=q.segments{k};ix=z.startIndex:z.endIndex;limits(k)=.75*max(baseline(ix,4));threshold(ix)=limits(k);
end
threshold(q.waiting>.5)=min(limits);
run=simulate_threshold_nmpc(plant,q,inputs,inputs,threshold,o);
fprintf('ORACLE DIAGNOSTIC: maximum excess %.6f N m.\n',max(max(run.actual(:,4)-threshold,0)));
save(fullfile(fileparts(mfilename('fullpath')),'oracle_group3.mat'),'run','baseline','threshold','-v7');
