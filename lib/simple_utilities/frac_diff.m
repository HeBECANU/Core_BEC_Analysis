function out=frac_diff(x,y,method)
%calculate the fractional difference

% TODO
% document
% test
% more methods eg take min or max as denominator

% input
% x,y    - numerical array, any dimensionality
% method - how the denominator should be determined can be 'mean','max','min'
% output - fractional difference with the same dimensionality as the input

if nargin<3
   method='mean';
end

if size(x)~=size(y)
    fprintf('error the sizes are different')
end

if strcmp(method,'mean')
    out=(x-y)./((x+y)./2);
elseif strcmp(method,'max')
     out=(x-y)./(max(x,y));
elseif strcmp(method,'min')
    out=(x-y)./(min(x,y));
elseif strcmp(method,'x')
    out=(x-y)./x;    
elseif strcmp(method,'y')
    out=(x-y)./y;    
elseif strcmp(method,'x abs')
    out=(x-y)./abs(x);    
elseif strcmp(method,'y abs')
    out=(x-y)./abs(y);      
else
    error('invalid method input')
end 

end