function count_rate= count2flux(counts,nedges)

if ~iscell(nedges)
    count_rate=counts./diff(nedges);
else
    dims = size(counts);
    vol = ones(dims);%ndim volume size of each bin
    for ii = 1:size(nedges,2)
        edges = nedges{ii};
        bin_width = diff(edges);
        temp_dim = dims;
        temp_dim(ii) = 1;
        or = ones(1,size(nedges,2));
        or(ii) = length(diff(edges));
        bin_width = reshape(bin_width,or);
        vol = vol .* repmat(bin_width,temp_dim);        
    end
    count_rate=counts./vol;
end
end