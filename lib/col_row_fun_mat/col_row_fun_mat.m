function out=col_row_fun_mat(funhandle,mat_in,dirn)
% col_row_fun_mat - apply a function to the rows or column of a matrix
% This code alows the user to apply a function to the rows of columns of a matrix.
% It is indended that this may substantialy simplify some common codeblocks and remove some messy for loops
% Syntax:   col_row_fun_mat(funhandle,matrix_in,direction)
% Inputs:
%   funhandle  - function handle, takes a column vector as input and returns a vector (row or col) of len N
%              - if matrix is returned then is flattened with out(:)
%              - must have consistent number of outputs
%   mat_in     - [M x Q] matrix, the matrix that funhandle will operate on one dimension of
%   dirn       - direction, direction to operate on, 1 operates on columns & returns the funhandle output on columns, 2 operates
%   on rows and retuns a column vector

% Outputs:
%    out       - matrix, for dirn=1 [N*Q] for dirn=2 [M*N]
%
% Example1: 
%         test_mat=rand(5,10);
%         testfun=@(t) 100*sin(2*pi*t*100+pi)+1*sin(2*pi*t*133+pi);
%         col_row_fun_result=col_row_fun_mat(@sum,test_mat,2);
%         inbuilt_result=sum(test_mat,2);
%         isequal(inbuilt_result,col_row_fun_result)
% Example2: 
%         test_mat=rand(5,10);
%         core_function=@(x) norm((sin(x*3)+1).^2);
%         % the old way with explicit looping
%         iimax=size(test_mat,1);
%         reslut_old=zeros(iimax,1);
%         for ii=1:iimax
%             reslut_old(ii)=core_function(test_mat(ii,:));
%         end
%         % the new way is a single line solution that is nearly impossible to break
%         % and is far easier to maintain
%         col_row_fun_result=col_row_fun_mat(core_function,test_mat,2);
%         isequal(col_row_fun_result,reslut_old)

% Other m-files required: none
% Also See: test_col_row_fun_mat.m
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%    - speed improvements may be possible with less of the x(:) syntax
%    - could try to catch when the function output nonuniform len
%    - build for higher dimension operation
%
% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][col_row_fun_mat]' in the subject line OR I WILL NOT REPLY
% Last revision:2019-05-11
%------------- BEGIN CODE --------------

if dirn==1
    col_slice=mat_in(:,1);
    first_out=funhandle(col_slice(:));
    first_out=first_out(:);
    first_out_size=size(first_out);
    %fun_return_row_or_col=first_out_size(2)==1;
    %fun_return_scalar=isscalar(first_out_size);
    iimax=size(mat_in,2); 
    out=zeros(max(first_out_size),iimax); %TODO use fun_return_row_or_col to chose first_out_size
    out(:,1)=first_out;
    for ii=2:iimax
        col_slice=mat_in(:,ii);
        out_slice=funhandle(col_slice(:));
        out(:,ii)=out_slice(:);  %TODO use optional rotation here
    end
elseif dirn==2
    row_slice=mat_in(1,:); %TODO use rotation
    first_out=funhandle(row_slice(:));
    first_out=first_out(:);
    first_out_size=size(first_out);
    %fun_return_row_or_col=first_out_size(2)==1;
    %fun_return_scalar=isscalar(first_out_size);
    iimax=size(mat_in,1); 
    out=zeros(iimax,max(first_out_size)); %TODO use fun_return_row_or_col to chose first_out_size
    out(1,:)=first_out;
    for ii=2:iimax
        row_slice=mat_in(ii,:);
        out_slice=funhandle(row_slice(:));
        out(ii,:)=out_slice(:);
    end
else 
    error('dirn not 1,2')
end
end