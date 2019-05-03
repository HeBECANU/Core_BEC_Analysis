function [min_val,idxy]=min_diff_vec(in,diff_tol)
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

% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][min_diff_vec]' in the subject line OR I WILL NOT REPLY
% Last revision:2019-05-03
%------------- BEGIN CODE --------------


%bugs/improvements
% could save a tiny bit of time if idxy is not to be returned by not storing order or doing the lookup at the
% end line

[in_sorted,order]=sort(in);
diff_vec=abs(diff(in_sorted));
if nargin>1 && ~isempty(diff_tol)
    diff_vec(diff_vec<diff_tol)=inf;
end

[min_val,idx]=min(diff_vec);
idxy=sort(order([idx,idx+1]));

end