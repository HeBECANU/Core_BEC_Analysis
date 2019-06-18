function data_out = centre_BEC(xyz_in,threshold,lim,bin_size)

% A centreing algorithm for a BEC in cartesian coordinates.
% This method builds a histogram and finds the centre of the region where
% the atom flux exceeds a given value. This is slower than using the COM
% but it is less biased by detector saturation.
% INPUTS
% For each dim: Min, Max, bin_size

    ndims = size(xyz_in,2);
    bec_cent = zeros(ndims,1);
    for axis = 1:ndims
        X = sort(xyz_in(:,axis));
        edges = lim(axis,1):bin_size(axis):lim(axis,2);
        centres = 0.5*(edges(2:end) + edges(1:end-1));
        counts = hist_adaptive_method(X,edges',1);
        bec_mask = counts(2:end-1) > threshold;
        bec_locn = centres(bec_mask);
        bec_cent(axis) = median(bec_locn);
    end
%         data_out = [];
      data_out = bec_cent';
%     data_out = xyz_in - bec_cent'; %thank you broadcasting
end

