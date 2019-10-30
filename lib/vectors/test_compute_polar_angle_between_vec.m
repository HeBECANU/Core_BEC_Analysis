%test_compute_polar_angle_between_vec

% life lesson
% dont test code with vectors with approx unit norm

rad2deg(compute_polar_angle_between_vec([0,1,0]*100,[1,0,0]*100,[0,0,1]*100) )
rad2deg(compute_polar_angle_between_vec([0,1,1]*100,[1,0,0]*100,[0,0,1]*100) )
rad2deg(compute_polar_angle_between_vec([0,0,-1]*100,[1,0,0]*100,[0,0,1]*100) )
rad2deg(compute_polar_angle_between_vec([0,-1,0]*100,[1,0,0]*100,[0,0,1]*100) )


%%

rad2deg(compute_polar_angle_between_vec([[0,1,0];[0,1,1];[0,0,-1];[0,-1,0]],[1,0,0],[0,0,1]) )


%%
rad2deg(compute_polar_angle_between_vec([1,1,0], [[1,0,0];[1,1,0]], [[0,0,1];[0,0,1]] ) )

%%
rad2deg(compute_polar_angle_between_vec([1,1,1], [[1,0,0];[1,1,0]], [[0,0,1];[0,0,1]] ) )

%% hangling when the angle reference and pointing reference point in the same dir
% single input
rad2deg(compute_polar_angle_between_vec([1,1,1],[1,0,0],[2,0,0]) )

% multi input
rad2deg(compute_polar_angle_between_vec([1,1,1],[[1,0,0];[1,0,0]],[[2,0,0];[0,0,1]]) )

%% zero vector

rad2deg(compute_polar_angle_between_vec([0,0,0],[0,1,0],[0,0,1]) )
rad2deg(compute_polar_angle_between_vec([1,0,0],[0,0,0],[0,0,1]) )
rad2deg(compute_polar_angle_between_vec([1,0,0],[0,1,0],[0,0,0]) )


%% vector scaling
% should be invariant under scaling
scale_factor=1e-320;
%scale_factor=1e120;
invec=[0.5+pi/1e3,sqrt(2),-1]*100;
pointrefvec=[1+48/7812,pi/1e4,pi/1e4]*100;
anglerefvec=[sqrt(7)/100,sqrt(2)/100,pi]*100;
baseline=rad2deg(compute_polar_angle_between_vec(invec,pointrefvec,anglerefvec) )
baseline-rad2deg(compute_polar_angle_between_vec(scale_factor.*invec,pointrefvec,anglerefvec) )
baseline-rad2deg(compute_polar_angle_between_vec(invec,scale_factor.*(pointrefvec),anglerefvec) )
baseline-rad2deg(compute_polar_angle_between_vec(invec,pointrefvec,scale_factor.*anglerefvec) )

% this can be slightly improved by turning on intermediate normalization but its already very good