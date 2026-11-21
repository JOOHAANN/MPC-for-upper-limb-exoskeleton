function [raw,muscleProfiles] = predict_528_MAP(model,loadKg,speedCode,phase)
%PREDICT_528_MAP Evaluate the original back MAP and the abdominal MAP.
% The repeated IL term is intentional: preserve the author's original rule.
b=model.back;
weights=condition_weights(b,loadKg,speedCode);
MEP=b.primitives*weights;
% Optional raw six-channel output, before fusion and amplitude calibration.
muscleProfiles=interp1(linspace(0,1,200)',MEP,phase(:),'linear','extrap');
back=(MEP(:,1)+MEP(:,2))/12+(MEP(:,3)+MEP(:,4))/4 ...
    +(MEP(:,1)+MEP(:,2))/6;
back=interp1(linspace(0,1,200)',back,phase(:),'linear','extrap');
abdomen=predict_abdominal_MAP(model.abdomen,loadKg,speedCode,phase);
raw=[back abdomen];
end

function weights=condition_weights(model,loadKg,speedCode)
weights=model.baseline;
for primitive=1:size(weights,1)
    for muscle=1:size(weights,2)
        weightPoly=model.loadCoefficients(:,primitive,muscle);
        speedPoly=model.speedCoefficients(:,primitive,muscle);
        weights(primitive,muscle)=weights(primitive,muscle) ...
            +polyval(weightPoly,loadKg)-polyval(weightPoly,10) ...
            +polyval(speedPoly,speedCode)-polyval(speedPoly,45);
    end
end
weights=max(0,weights);
end
