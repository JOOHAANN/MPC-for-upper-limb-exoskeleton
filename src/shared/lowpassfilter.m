function y = lowpassfilter(x, fs, fc)
%LOWPASSFILTER Apply a causal second-order Butterworth low-pass filter.
%   fs and fc are the sampling and cutoff frequencies in Hz. Preserve the
%   causal filter used by the original model; filtfilt changes its phase.
    normalizedCutoff = fc / (fs / 2);
    [numerator, denominator] = butter(2, normalizedCutoff, 'low');
    y = filter(numerator, denominator, x);
end
