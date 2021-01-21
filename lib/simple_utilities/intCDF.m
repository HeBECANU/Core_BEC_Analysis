function I = intCDF(domain,data,varargin)
% Returns the cumulative distribution function of data x over
% a domain d.  Y = intCDF(x,d) returns a row vector of the same size as x, 
% whose entries Y(i) counts the number of elements of d that are less than x(i).
% Inputs:
%   x - data    - 1D real vector 
%   d - domain  - 1D real vector
% Example
% ```
% domain = 0:5;
% data = [1,3,3.1]
% isequal(intCDF(domain,data),[0,0,1,1,3,3])
% ```
% Optionally normalizes the output
% ```
%  domain = 0:5;
% data = [1,2,3.1]
% isequal(intCDF(domain,data,'normalize',true),[0,0,1/3,2/3,1,1])
% ```
    p = inputParser;
    addParameter(p,'normalize',false);
    addParameter(p,'reverse',false);
    parse(p,varargin{:});
    nmlz = p.Results.normalize;
    reverse = p.Results.reverse;
    domain = col_vec(domain);
    data = col_vec(data);
    I = zeros(size(domain));
    if ~issorted(domain), warning('Domain for intCDF is not sorted');end
    for xi = 1:length(domain)
        I(xi) = sum(data<domain(xi));
    end
    if nmlz
        % returns normalized CDF i.e. final value is 1
        I = I/length(data);
    end
    if reverse
        I = I(end) - I;
    end
    I = I';
end
