function activation = rescale_MAP_by_load(raw,loadKg,scaling)
%RESCALE_MAP_BY_LOAD Apply one fixed affine map per load and MTU.
% Neither speed nor the current experimental trial is an argument.
row=find(scaling.loadKg==loadKg,1);
assert(~isempty(row),'No experimental amplitude calibration for this load.');
activation=max(0,min(1,raw.*scaling.gain(row,:)+scaling.offset(row,:)));
end
