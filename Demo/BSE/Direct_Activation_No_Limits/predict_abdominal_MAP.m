function [activation,outOfDomain] = predict_abdominal_MAP(model,loadKg,speedCode,phase)
%PREDICT_ABDOMINAL_MAP Original additive load/speed MAP weight formulation.
outOfDomain=loadKg<5 || loadKg>15 || speedCode<25 || speedCode>65;
weights=model.baseline;
for primitive=1:4
    loadPoly=model.loadCoefficients(:,primitive);
    speedPoly=model.speedCoefficients(:,primitive);
    weights(primitive)=weights(primitive)+polyval(loadPoly,loadKg)-polyval(loadPoly,10) ...
        +polyval(speedPoly,speedCode)-polyval(speedPoly,45);
end
profile=model.amplitudeScale*model.primitives*max(0,weights);
activation=interp1(model.phase,profile,phase(:),'linear','extrap');
end
