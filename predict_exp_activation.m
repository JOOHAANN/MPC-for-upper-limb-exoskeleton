function prediction = predict_exp_activation(subject, loadKg, cadenceBPM, timeSeconds)
% Predict activation from a TRAINED subject MAP, without recorded EMG or IK.
% p = predict_exp_activation('240730Exp_A01', 10, 9);
% plot(p.timeSeconds, p.backActivation);
% cadenceBPM counts complete action bouts per minute, not file speed codes.
arguments
    subject (1,:) char
    loadKg (1,1) double {mustBeFinite,mustBeNonnegative}
    cadenceBPM (1,1) double {mustBeFinite,mustBePositive}
    timeSeconds double = []
end
root = fileparts(mfilename('fullpath'));
file = fullfile(root,'data','processed','exp_subjects_A01',[subject '_map.mat']);
assert(isfile(file),'Trained MAP missing for %s. Run demo_exp_subjects(''run'').',subject);
cache = load(file,'mapModel');
model = cache.mapModel;
if isempty(timeSeconds), timeSeconds = (0:0.01:60/cadenceBPM)'; end
timeSeconds = timeSeconds(:);
assert(all(isfinite(timeSeconds)) && all(timeSeconds >= 0) && ...
    all(diff(timeSeconds) > 0),'Time must be finite, increasing and nonnegative.');
features = [loadKg/15, (loadKg/15)^2, cadenceBPM/10, (cadenceBPM/10)^2];
bounded = min(max(features,model.covariateMinimum(:)'),model.covariateMaximum(:)');
design = [1,(bounded-model.mu(:)')./model.scale(:)'];
% Python stores the four-by-six coefficient weights in row-major order.
weights = reshape(max(0,design*model.beta),6,4)';
profile = model.primitives*weights;
phase = min(max(timeSeconds*cadenceBPM/60,0),1);
six = interp1(model.phase(:),profile,phase,'linear');
prediction = struct('timeSeconds',timeSeconds,'phase',phase, ...
    'sixBackActivation',six, ...
    'backActivation',min(max(six*model.lumpWeights(:),0),1), ...
    'channelNames',{{'IL_L1_r','IL_L1_l','LTpT_T7_r','LTpT_T7_l','LTpL_L1_r','LTpL_L1_l'}}, ...
    'subject',subject,'loadKg',loadKg,'cadenceBPM',cadenceBPM, ...
    'covariatesOutsideTrainingDomain',any(abs(features-bounded)>1e-10), ...
    'timeOutsideOneAction',timeSeconds>60/cadenceBPM, ...
    'sixChannelRangePassed',all(six(:)>=0 & six(:)<=1), ...
    'usesRecordedEMGOrReference',false,'modelFile',file);
end
