function X = midpoints(x,varargin)
% Returns an (n-1) array whose values are the midpoints of each successive
% pair within an (n)-array. Eg
% ```
% isequal(midpoints(0:3), [.5 1.5 2.5])
% ```
% 

X = 0.5*(x(2:end)+x(1:end-1));
end