function out = vert_or_horz_cat(varargin)
% a function that adaptively concatonates the right dimension
try
    out = vertcat(varargin{:});
catch
    out = horzcat(varargin{:});
end