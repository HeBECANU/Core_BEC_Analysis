function logical_out = in_intervals(A,lims,varargin)
% in_interval(A,L) returns a binary array of size(A) with entries that are
% TRUE if the corresponding element of A lies within any of the intervals
% specified by the N x 2 matrix L of interval bounds. This generalizes the 
% expression lim(1) <= A < lim(2). Note that L does not require L(:,2) >
% L(:,1).
% Optional kwargs:
%   mode       'open'   Strict inequality m < A < M
%              'closed' Inclusive limits  m <= A <= M
%              'upper'  Reverse inequality m < A <= M
% Example:
% ```
% A = 1:6;
% L = [0,3
%      5,6];
% isequal(in_intervals(A,L),[1,1,0,0,1,0])
% ```

p = inputParser;
addParameter(p,'mode','default');
parse(p,varargin{:});
mode = p.Results.mode;

if numel(lims) < 2
    error('Limits must be pairs of values');
end

if size(lims,2) ~= 2 && size(lims(1)) == 2
    lims = lims';
end

logical_out = false(size(A));
if strcmp(mode,'default')
    for l_ax = 1:size(lims,1)
        logical_out = logical_out | (A>=min(lims(l_ax,:)))&(A<max(lims(l_ax,:)));
    end
elseif strcmp(mode,'open')
    for l_ax = 1:size(lims,1)
        logical_out = logical_out | (A>min(lims(l_ax,:)))&(A<max(lims(l_ax,:)));
    end
elseif strcmp(mode,'closed')
    for l_ax = 1:size(lims,1)
        logical_out = logical_out | (A>=min(lims(l_ax,:)))&(A<=max(lims(l_ax,:)));
    end
elseif strcmp(mode,'upper')
    for l_ax = 1:size(lims,1)
        logical_out = logical_out | (A>min(lims(l_ax,:)))&(A<=max(lims(l_ax,:)));
    end
end



end