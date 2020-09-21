function y = bound(x,lower,upper)
% Returns input x with values bounded by [lower,upper]
% Inputs:
%     x     Real array (any dimension)
%    min    Real lower bound
%    max    Real upper bound
% Example:
% ```
% min = 2
% max = 4
% isequal(bound(1:5,min,max),[2 2 3 4 4])
% ```
  y=min(max(x,lower),upper);
end