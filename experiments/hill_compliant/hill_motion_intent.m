function signal = hill_motion_intent(value, tonic, gain, mode)
% Normalize an EMG-derived intent signal; no reference angle is an input.
compensation = gain * tonic;
switch mode
    case 'soft_tonic'
        positive = max(0, value);
        if compensation > 0
            signal = positive.^2 ./ (positive + compensation);
        else
            signal = positive;
        end
    case 'divisive'
        signal = value ./ (1 + compensation);
    case 'subtractive'
        signal = value - compensation;
    otherwise
        error('hill:decoder', 'Unknown motion intent decoder: %s', mode);
end
signal = max(0, min(1, signal));
end
