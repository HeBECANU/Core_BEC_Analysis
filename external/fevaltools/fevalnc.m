function c = fevalnc(nout, fn, varargin)
%FEVALNC
%   c = FEVALNC(nout, fn, x1, ..., xn)
%
%   Call a function with a given number of output arguments and capture the
%   outputs to a single cell array.
%
%   This is especially useful when multiple output-argument functions need 
%   to be manipulated inside an anonymous function.
%
%   Inputs:
%       nout - number of output arguments requested
%       fn - function handle to be called
%       x1, ..., xn - the arguments to fn
%   Outputs:
%       c - the requested outputs from fn captured in a single cell array
%
%   Examples:
%       c = fevalnc(2 ,@max, [3 1 2]);
%       disp(c)
%          [3]    [1]
%
%       %get the difference from the max of pre and post element
%       prepostmaxdiff = @(x) feval(...
%                   @(x,c) x([c{2}-1 c{2}+1]) - c{1}, ...
%                   x, fevalnc(2,@max,x) ...
%               );
%       disp(prepostmaxdiff([1 3 5 4 0]))
%            -2    -1
%
%       %get element that is one index in dimension 2
%       %beyond the first negative element
%       slice2p1 = @(x) feval(...
%                   @(x,c) x(c{1},c{2}+1,c{3}), ...
%                   x, fevalnc(3, @ind2sub, size(x), find(x<0,1)) ...
%               );
%       
%
%   See also FEVAL, FEVALI, FEVALIC.
%
%   Copyright, Benjamin P Davis, 2020

c = fevalic(1:nout, fn, varargin{:});