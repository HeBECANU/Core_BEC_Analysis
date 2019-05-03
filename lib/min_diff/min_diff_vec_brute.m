function [min_val,idxy]=min_diff_vec_brute(in)
%min_diff - fing the minimum absolute difference between every element and all other elements in a vector (exluding itself)

% Syntax:  [min_val,idxy]=min_diff_vec(data_vector,tol[optional])
% Example: 


% Other m-files required: none
% Also See:
% Subfunctions: 
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%   -add option for when difference must be greater than some tolerance
%   -try a sorted vector version
% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][min_diff_vec]' in the subject line OR I WILL NOT REPLY
% Last revision:2019-05-03
%------------- BEGIN CODE --------------


%bugs/improvements
%could get a factor of 2 speedup by only looking for the min in the upper or lower diag;
% calculate the difference matrix 
% equiv to minus(in,in') for R2016b and later

diff_mat=abs(bsxfun(@minus, in ,in'));
sd=size(diff_mat);
n=sd(1);
diff_mat=diff_mat+diag(NaN(n,1));
[min_val,idx]=nanmin(diff_mat(:));
[idx1,idx2] = ind2sub(sd,idx);
idxy=sort([idx1,idx2] );
end

% sorted version ideas
% the basic implementation has complexity O(n^2) to get the differences
% then O(n^2) to find the mimimum

% a sorted version 
% O(nlog(n)) for the sort
% O(n-1) differences of neibours
% O(n) to find the minimum
% O(2) to look up the postion in the sort table



%% some other slower implementations
% i dont like that it generates the symetric difference matrix which is inherently wastefull in memory
% also the min operation is done with twice the number of inputs and in the implmentation above must use
% nanmin

%variant 1
% this masks the difference matrix with is then fed into a normal min
% ends up being slower becasue of the mask and is also very hard to recover the indicies
    % diff_mat=abs(bsxfun(@minus, in ,in'));
    % sd=size(diff_mat);
    % mask=~triu(true(sd));
    % [min_val,idx]=nanmin(diff_mat(mask));