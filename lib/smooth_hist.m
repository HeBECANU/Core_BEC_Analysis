function out_struct=smooth_hist(in_struct)

if ~isfield(in_struct,'xsorted')
    in_struct.xsorted=false;
end

if isfield(in_struct,'bin_factor') && isfield(in_struct,'bins')
    warning('%s:bins ignored using bin_factor to calculate',mfilename)
    in_struct=rmfield(in_struct,'bins');
end

if ~isfield(in_struct,'bin_factor') && ~isfield(in_struct,'bins')
    in_struct.bin_factor=10;
end


if ~isfield(in_struct,'bins')
    in_struct.bins=abs(diff([in_struct.min,in_struct.max]))/in_struct.sigma;
    in_struct.bins=round(in_struct.bin_factor*in_struct.bins);
end

in_struct.xdat=col_vec(in_struct.xdat);

%function that resturns a gaussian smoothed histogram
edges=linspace(in_struct.min,in_struct.max,in_struct.bins+1)';
centers=(edges(2:end)+edges(1:end-1))./2;
hist_counts_raw=hist_adaptive_method(in_struct.xdat,edges,in_struct.xsorted);

out_struct.counts.below=hist_counts_raw(1);
out_struct.counts.above=hist_counts_raw(end);
out_struct.counts.raw=hist_counts_raw(2:end-1);

if in_struct.sigma~=0 || ~isnan(in_struct.sigma)
    out_struct.counts.smooth=gaussfilt(centers,out_struct.counts.raw,in_struct.sigma);
else
    out_struct.counts.smooth=out_struct.counts.raw;
end

out_struct.count_rate.smooth=out_struct.counts.smooth./diff(edges);
out_struct.count_rate.raw=out_struct.counts.raw./diff(edges);

out_struct.bin.edge=edges;
out_struct.bin.centers=centers;

end

