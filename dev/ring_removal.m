function data = ring_removal(data,t_lim)
%small function to remove ringing from reconstructed data
out_counts_txy = cell(size(data.counts_txy));
out_num_counts = zeros(size(data.num_counts));
for ii = 1:length(data.counts_txy)
    num_counts = data.num_counts(ii);
    if num_counts == 0
        out_counts_txy{ii} = data.counts_txy{ii};
        continue
    end
    counts=data.counts_txy{ii};
    t_dif = -counts(1:end-1,1)+counts(2:end,1);
    t_dif_mask = logical([1;abs(t_dif)>t_lim]);
    out_counts_txy{ii} = counts(t_dif_mask,:);
    masked_counts = sum(~t_dif_mask);
    out_num_counts(ii) = num_counts-masked_counts;
end
data.counts_txy = out_counts_txy;
data.num_counts = out_num_counts;
end