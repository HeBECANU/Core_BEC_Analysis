function rand_out=rand_interval(limits,out_size)
% rand_interval - create a vector of random numbers that is between limits

% Syntax:   rand_out=rand_interval(limits,out_size)
%
% Inputs:
%	limits      - matrix, if limits is a [1x2] (or [2x1]) vector then the same limits are applied to each element
%                 if limits has one more dimension than out_size, and that dimension is of size 2 then
%                 the limit is set on each output, this can make creation of the limits matrix somewhat
%                 more compicated (than it could be) for the case when rand_out is a vector but it allows generalization
%   out_size	- size of the output matrix, only needed if limits is a [1x2] (or [2x1]) vector (setting the 
%                 same domain on each element) as otherwise it is given by the size of the limits matrix
% Outputs:
%    rand_out	- output matrix with each element having the domain set by the limits matrix

%
% Example: 
% out_interval=[1,2];
% rand_interval(out_interval,[1,10])
%
% out_size=[2,2];
% out_interval=reshape([[1,2];[2,3];[3,4];[4,5]],cat(2,out_size,2));
% test_out=rand_interval(out_interval)

% Other m-files required: col_vec
% Also See: test_rand_interval
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%   - none
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2020-01-22

%------------- BEGIN CODE --------------

if nargin<2
    if ndims(limits)<2
        error('out size needed if limits is only of dimension 2')
    end
    out_size=size(limits);
    if out_size(end)~=2
        error('highest dimension of limits must have size 2')
    end
    out_size=out_size(1:end-1);
end


% the same domain for all outputs
out_size=col_vec(out_size);
out_size=transpose(out_size); %row vector for rands

if isequal(size(limits),[1,2]) || isequal(size(limits),[2,1])
    limits=sort(limits);
    rand_out = limits(1) + (limits(2)-limits(1)).*rand(out_size);

else
    % a different domain for each output
    size_limits=size(limits);
    if ~(isequal(size_limits(1:end-1),out_size) || isequal(size_limits(end),2) )
        error('limits size should have one extra dimension than out_size and that dimension should be of size 2')
    end
    % we need to asign the max and min limit n tensor without knowing its dimensionality
    dimensions=numel(out_size);
    idx_sub_scripts=repmat({':'},[1,dimensions]);
    idx_sub_scripts=cat(2,idx_sub_scripts,1);
    limit_min=limits(idx_sub_scripts{:});
    
    idx_sub_scripts=repmat({':'},[1,dimensions]);
    idx_sub_scripts=cat(2,idx_sub_scripts,2);
    limit_max=limits(idx_sub_scripts{:});
    if any(limit_min(:)>limit_max(:))
        error('min is greater than max limit')
    end
    if any(abs(limit_min(:)-limit_max(:))<eps*10)
        warning('limits are within eps of each other')
    end
    rand_out = limit_min + (limit_max-limit_min).*rand(out_size);
end


end