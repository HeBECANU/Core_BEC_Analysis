function angle_out=compute_polar_angle_between_vec(in_vec,point_ref_vec,angle_ref_vec) 

% given some vector v, a pointing reference u and a angle reference w
% find the polar angle of the vector v about u such that the vector w would retun zero angle
% two ways to tackle the probjem
% both will require finding the component of w that is perpendicular to u
% w_perp(u)=w-(w·u)u

% the geometric way
% we calculate the component of v that is perpendicular to u
% v_perp(u)
% then find the angle between  v_perp(u) & w_perp(u)
% this gives the magnitude of the angle with no regards for the sign
% to get the sign we compute sign(u · ( v_perp(u) × w_perp(u)))


% cord transform
% we find a transformation such that the identity vectors i is mapped to the 
% w_perp(u) vector, k goes to u and then j is maped to the othogonal to these two
% transformation*I=vecnorm([w_perp(u);u×w_perp(u);u])
% and simply
% transformation=vecnorm([w_perp(u);u×w_perp(u);u])
% then this transformation is applied to the input 
% the resulting vectors are then converted to spherical cordinates

%handle matching input sizes
if size(point_ref_vec)~=size(angle_ref_vec)
    error('reference vectors must be the same size')
end
if size(in_vec,2)~=size(angle_ref_vec,2)
    error('inputs must have the same size in the 2nd dimension')
end

if size(in_vec,1)~=size(angle_ref_vec,1)
    if ~xor(size(in_vec,1)==1,size(angle_ref_vec,1)==1)
        error('first diemension of inputs must be matched or input_vector/reference must have length of 1')
    end
    if size(in_vec,1)==1
        in_vec=repmat(in_vec,[size(angle_ref_vec,1),1]);
    elseif size(angle_ref_vec,1)==1
        angle_ref_vec=repmat(angle_ref_vec,[size(in_vec,1),1]);
        point_ref_vec=repmat(point_ref_vec,[size(in_vec,1),1]);
    end
end

if ~isequal(size(point_ref_vec),size(angle_ref_vec),size(in_vec))
    error('inputs not matched after code should have matched them')
end

comp_angle_ref_perp_point_ref=angle_ref_vec-dot(angle_ref_vec,point_ref_vec,2).*point_ref_vec;
comp_in_vec_perp_point_ref=in_vec-dot(in_vec,point_ref_vec,2).*point_ref_vec;
angle_mag=angle_between_vec(comp_angle_ref_perp_point_ref,comp_in_vec_perp_point_ref);
angle_sign=sign(dot(point_ref_vec,cross(comp_in_vec_perp_point_ref,comp_angle_ref_perp_point_ref),2));
angle_sign(angle_sign==0)=1;
angle_out=angle_mag.*angle_sign;
angle_out=wrapTo2Pi(angle_out);

end


