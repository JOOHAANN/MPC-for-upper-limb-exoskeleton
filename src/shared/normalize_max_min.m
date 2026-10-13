function normalized_data = normalize_max_min(data, Max, Min)
%NORMALIZE_MAX_MIN Map a nonconstant vector to the requested range.
%   Preserve the original scaling equation and Max, Min argument order.
    data_max = max(data);
    data_min = min(data);
    normalized_data = Min + (data - data_min) / (data_max - data_min) * (Max - Min);
end
