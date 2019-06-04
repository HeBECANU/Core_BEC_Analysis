function out=taylor_series(x,derivs,offset)
x=x(:);
derivs=derivs(:);
if nargin==2
    offset=0;
end
%i want to calulate the result of a taylor series evaluation at the points x
% [f(offset),f'(offset),f'''(offset)...]
factor=factorial(0:(numel(derivs)-1))';
derivs=derivs./factor;
x=x-offset;
out=polyval(flip(derivs),x);
end