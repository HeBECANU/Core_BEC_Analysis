function masked_num = mask_square_num(in,lim,varargin)
%performs a mask square and counts the total number of counts in that
%masked region

if length(varargin)<1
    invert = 0;
else
    invert = varargin{1};
end

masked_num = zeros(1,length(in));

for ii = 1:length(in)
    this_counts = in{ii};
    masked_counts = mask_square(this_counts,lim,invert);
    masked_num(ii) = size(masked_counts,1);
end