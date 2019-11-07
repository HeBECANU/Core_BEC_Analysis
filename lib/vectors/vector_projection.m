function vector_out=vector_projection(vec_a,vec_b)
% calculate the vector projection of vector a onto vector b
% Inputs
%     vec_a   - [n x m] matrix where m is the dimensionality of the vector and n is the number of independent vectors
% Outputs:
%     vector_out - projection of a onto b calculated row-wise

% Other m-files required: none
% Also See: angle_between_vec,compute_polar_angle_between_vec
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%    - 
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2019-10-30


if size(vec_a)~=size(vec_b)
    error('vectors must be the same size')
end

%algo 1
vector_out=vec_b.*dot(vec_a,vec_b,2)./(vecnorm(vec_b,2,2).^2);

% algo 2
% this seems to have very slightly better numerical perfromance
% vector_out2=(dot(vec_a,vec_b,2)./vecnorm(vec_b,2,2)).*(vec_b./vecnorm(vec_b,2,2));



end