function forecast = nmpc_forecast(sequence, observed, planned, current, options)
%NMPC_FORECAST Trained MAP preview corrected by past activation residuals.
% Future experimental EMG, reference theta and reference torque are not read.
count = size(planned,1);
query = current + (1:options.points)'*options.predictionDt/sequence.dt;
query = max(1,min(count,query));
forecast.activation = interp1((1:count)',planned,query,'linear');
predictor = options.activationResidualPredictor;
history = max(1,current-predictor.lags(:));
pastResidual = observed(history,:)-planned(history,:);
features = reshape(pastResidual',1,[])./predictor.scale(:)';
correction = features*predictor.coefficients;
forecast.activation = forecast.activation + reshape(correction,2,[])';
forecast.waiting = interp1((1:count)',double(sequence.waiting),query,'linear');
forecast.mass = interp1((1:count)',double(sequence.mass),query,'linear');
forecast.activation = max(0,min(1,forecast.activation.*(1-forecast.waiting)));
end
