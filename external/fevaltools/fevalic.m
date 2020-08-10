function c = fevalic(inds, fn, varargin)
%FEVALIC
%   c = FEVALIC(inds, fn, x1, ..., xn)
%
%   Capture multiple output arguments from a function into a single cell 
%   array in any order.
%
%   This is especially useful when multiple output-argument functions need 
%   to be manipulated inside an anonymous function.
%
%   Inputs:
%       inds - numeric or logical indices denoting function output
%           argument order or masking
%       fn - function handle to be called
%       x1, ..., xn - the arguments to fn
%   Outputs:
%       c - the requested outputs from fn captured in a single cell array
%
%   Examples:
%       c = fevalic(2,@max, [3 1 2]);
%       disp(c)
%          [1]
%
%       c = fevalic([2 1] ,@max, [3 1 2]);
%       disp(c)
%          [1]    [3]
%
%       c = fevalic([false true], @max, [3 1 2]);
%       disp(c)
%          [1]
%
%       %get the max along with preceding and following number
%       prepostmax = @(x) feval(...
%                   @(x,c) x(c{1}-1:c{1}+1), ...
%                   x, fevalic(2,@max,x) ...
%               );
%       disp(prepostmax([1 3 5 4 0]))
%            3     5     4
%
%       %get element that is one index in dimension 2
%       %beyond the first negative element
%       slice2p1 = @(x) feval(...
%                   @(x,c) x(c{1},c{2}+1,c{3}), ...
%                   x, fevalic(1:3, @ind2sub, size(x), find(x<0,1)) ...
%               );
%       
%
%   See also FEVAL, FEVALI, FEVALNC.
%
%   Copyright, Benjamin P Davis, 2020

if islogical(inds)
    nout = length(inds);
else
    nout = max(inds);
end
[c{1:nout}] = feval(fn, varargin{:});
c = c(inds);