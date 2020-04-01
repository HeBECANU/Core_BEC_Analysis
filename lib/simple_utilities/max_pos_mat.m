function [pos,max_num]=max_pos_mat(mat_in)
%max_pos_mat - find the position of the maximum element in a N dimensional input
%
% Syntax:  [pos,max_num]=max_pos_mat(mat_in)
%
% Inputs:
%    mat_in            - matrix, any dimensionality
%
% Outputs:
%    pos        - vector with number of elments equal to dimension of matrix
%    max_num    - max value of the matrix
%
% Example: 
%	a=magic(5);
%	[pos,max_num]=max_pos_mat(a)
%	pos=num2cell(pos)
%   isequal(max_num,a(pos{:}))

% Other m-files required: none
% Also See: none
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%   - none
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2020-01-19

%------------- BEGIN CODE --------------

[max_num,max_idx] = max(mat_in(:));
pos=cell(size(size(mat_in)));
[pos{:}]=ind2sub(size(mat_in),max_idx);
pos=cell2mat(pos);
end