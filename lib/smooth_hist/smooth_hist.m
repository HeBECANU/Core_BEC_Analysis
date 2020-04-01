function out_struct=smooth_hist(xdata,varargin)
%smooth_hist - a histogram algorithm based on binary search of counts

% Syntax:         bin_counts=count_search_hist(data,edges)
% Equivelent to:  bin_counts=histcounts(data,[-inf;edges;inf])
% Designed to replicate histcounts(X,edges) "The value X(i)
%is in the kth bin if edges(k) ? X(i) < edges(k+1)" 
% Inputs:
%    in_struct          - structure of all the input arguments
%    xdata     - column vector of bin edges, MUST BE ORDERED!

%
% Outputs:
%    details_struct - column vector, with length numel(edges)+1,  the first(last) element 
%                are the number of counts below(above) the first(last) edge
%    details_struct.count_rate
%    details_struct.count_rate.smooth
%
% Example: 
%     
% Other m-files required: col_vec,hist_adaptive_method
% Also See: 
% Subfunctions: 
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%  - handle no input case
% 
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2019-07-07

%------------- BEGIN CODE --------------


p = inputParser;
is_lims=@(x) isequal(size(x),[1,2]) && isnumeric(x);
is_c_logical=@(in) isequal(in,true) || isequal(in,false); %can x be cast as a logical
addOptional(p,'lims',[],is_lims);
addOptional(p,'sorted',false,is_c_logical);
addOptional(p,'sigma',nan,@isnumeric);
addOptional(p,'bin_width',nan,@isnumeric);
addOptional(p,'bin_num',nan,@isnumeric);
addOptional(p,'bin_factor',nan,@isnumeric);
addOptional(p,'max_bins',2e6,@isnumeric);
addOptional(p,'doplot',false,is_c_logical);
parse(p,varargin{:});
parsed_input= p.Results;

if nargout==0
    parsed_input.doplot=true;
end


if sum(~isnan([parsed_input.bin_factor,parsed_input.bin_width,parsed_input.bin_num]))>1
    error('must pass only one of ''bin_factor'' or ''bin_width'' or ''bin_num''  ')
end

if ~isnan(parsed_input.sigma) && sum(~isnan([parsed_input.bin_width,parsed_input.bin_num,parsed_input.bin_factor]))==0
    %using default bin factor
    parsed_input.bin_factor=5;
end

if (isnan(parsed_input.sigma) || parsed_input.sigma==0) && ~isnan(parsed_input.bin_factor)
    error('cant specify bin width using bin_factor when sigma is zero')
end

if  isempty(parsed_input.lims) || sum(isnan(parsed_input.lims)) ~= 0
    bin_limits=[nanmin(xdata),nanmax(xdata)];
else
    bin_limits=parsed_input.lims;
end

if isnan(parsed_input.sigma) && sum(~isnan([parsed_input.bin_width,parsed_input.bin_num,parsed_input.bin_factor]))==0
    %using default bin factor
    parsed_input.sigma=range(bin_limits)*5e-3;
    parsed_input.bin_factor=5;
end

if ~isnan(parsed_input.bin_factor)
    bin_width=parsed_input.sigma/parsed_input.bin_factor;
    x_bin_num=floor(range(bin_limits)/bin_width);
    bin_limits=[min(bin_limits),min(bin_limits)+x_bin_num*bin_width];
elseif ~isnan(parsed_input.bin_width)
    x_bin_num=floor(range(bin_limits)/bin_width);
    bin_limits=[min(bin_limits),min(bin_limits)+x_bin_num*bin_width];
elseif ~isnan(parsed_input.bin_num)
    x_bin_num=parsed_input.bin_num;
end

if x_bin_num==0
   error("auto binning has returned no bins, try seting the limits with the 'lim'")
end

if x_bin_num>parsed_input.max_bins
    warning('%s: number of bins exceeds %g will reduce bin number to this',mfilename,parsed_input.max_bins)
    x_bin_num=parsed_input.max_bins;
end

sigma=parsed_input.sigma;

bin_width=range(bin_limits)/x_bin_num;
if bin_width>sigma*5
     warning('%s: bin width was much larger than sigma will not smooth',mfilename)
     sigma=0;
end

if ndims(xdata)<3
    xdata=col_vec(xdata);

    %function that resturns a gaussian smoothed histogram
    edges=linspace(min(bin_limits),max(bin_limits),x_bin_num+1)';
    centers=(edges(2:end)+edges(1:end-1))./2;

    hist_counts_raw=hist_adaptive_method(xdata,edges,parsed_input.sorted,1);
    out_struct.counts.below=hist_counts_raw(1);
    out_struct.counts.above=hist_counts_raw(end);
    out_struct.counts.raw=hist_counts_raw(2:end-1);
else
    for ii = 1:ndims(xdata)
        nedges=linspace(min(bin_limits),max(bin_limits),x_bin_num+1)';
        edges{ii}=nedges;
        centers{ii}=(nedges(2:end)+nedges(1:end-1))./2;
    end

    hist_counts_raw=histcn(xdata,edges);
    out_struct.counts.raw=hist_counts_raw;
end


if sigma~=0 || ~isnan(sigma)
    out_struct.counts.smooth=gaussfiltn(centers,out_struct.counts.raw,sigma);
else
    out_struct.counts.smooth=out_struct.counts.raw;
end

out_struct.count_rate.smooth=out_struct.counts.smooth./diff(edges);
out_struct.count_rate.raw=out_struct.counts.raw./diff(edges);

out_struct.bin.edge=edges;
out_struct.bin.centers=centers;

if parsed_input.doplot
    stfig('smooth histogram')
    plot(out_struct.bin.centers,out_struct.count_rate.smooth,'k')
end


end

