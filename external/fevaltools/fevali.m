function varargout = fevali(inds, fn, varargin)
%FEVALI
%   [y1,...,yn] = FEVALI(inds, fn, x1, ..., xn)
%
%   Return multiple arguments from a function in any order.
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
%       y1,...,yn - the outputs from fn corresponding to inds reordered 
%           as requested
%
%
%   Examples:
%       disp(fevali(2,@max, [3 1 2]));
%           1
%
%       [i, m] = fevali([2 1] ,@max, [3 1 2])
%       %i = 1, m = 3;
%
%       disp(fevali([false true], @max, [3 1 2]))
%           1
%
%       %get the max along with preceding and following number
%       prepostmax = @(x) feval(...
%                   @(x,ind) x(ind-1:ind+1), ...
%                   x, fevali(2,@max,x) ...
%               );
%       disp(prepostmax([1 3 5 4 0]))
%            3     5     4
%
%   See also FEVAL, FEVALIC, FEVALNC.
%
%   Copyright, Benjamin P Davis, 2020

varargout = fevalic(inds, fn, varargin{:});