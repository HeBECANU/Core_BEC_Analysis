function mat_out=linspace_mat(x1,x2,num_steps,use_logspace,fix_vec_in)
%linspace_mat - a simple extenstion to linspace that steps between matricies of any size
%
% Syntax:  mat_out=linspace_mat(x1,x2,num_steps,use_logspace)
%
% Inputs:
%   x1           - matrix, any dimensionality, starting point
%   x2           - matrix, any dimensionality, end point
%   num_steps    - number of steps to use
%   use_logspace - step in logspace?
%   fix_vec_in   - defalut true, fix the behaviour of vector inputs for x1,x2 so that the linspace is 
%                   added to the singleton dimension instead of the 3rd dim (bc matlab does not have real vectors)
%                   cant see the use case here but it may exist
%
% Outputs:
%    pos        - vector with number of elments equal to dimension of matrix
%    max_num    - max value of the matrix
%
% Example: 
%   mat_in=magic(5);
%   a=linspace_mat(mat_in,mat_in+1,10);

% Other m-files required: none
% Also See: none
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2020-01-19

%------------- BEGIN CODE --------------


if size(x1)~=size(x2)
    error('start and end vectors must be the same size')
end

if nargin<3
    num_steps=100;
end

if nargin<4 || isempty(use_logspace) || isnan(use_logspace)
    use_logspace=false;
end

if nargin<5 || isempty(fix_vec_in) || isnan(fix_vec_in)
    fix_vec_in=true;
end



size_tmp=size(x1);
size_tmp=cat(2,size_tmp,num_steps);
mat_out=zeros(size_tmp);

% to loop over all the positions in the output matrix i will use ind2sub to convert from the single index to 
% the subscript index, this allows arb dimensions to be handled
x1_size_tmp=size(x1);
iimax=prod(x1_size_tmp);
sub_tmp=cell(size(x1_size_tmp));
if use_logspace
    for ii=1:iimax
        [sub_tmp{:}]=ind2sub(x1_size_tmp,ii);
        mat_out(sub_tmp{:},:)=logspace(x1(ii),x2(ii),num_steps);
    end   

else
    for ii=1:iimax
        [sub_tmp{:}]=ind2sub(x1_size_tmp,ii);
        mat_out(sub_tmp{:},:)=linspace(x1(ii),x2(ii),num_steps);
    end   
end


% deal with the vector case for the start end pt
if ~isscalar(x1) && iscolumn(x1) && fix_vec_in
    mat_out=permute(mat_out,[1,3,2]);
elseif ~isscalar(x1) && isrow(x1) && fix_vec_in
    mat_out=permute(mat_out,[3,2,1]);
end


end