%test_compute_polar_angle_between_vec

rad2deg(compute_polar_angle_between_vec([0,0,0],[1,0,0],[0,0,1]) )
rad2deg(compute_polar_angle_between_vec([0,1,0],[1,0,0],[0,0,1]) )
rad2deg(compute_polar_angle_between_vec([0,0,-1],[1,0,0],[0,0,1]) )
rad2deg(compute_polar_angle_between_vec([0,-1,0],[1,0,0],[0,0,1]) )


%%

rad2deg(compute_polar_angle_between_vec([[0,0,0];[0,1,0];[0,0,-1];[0,-1,0]],[1,0,0],[0,0,1]) )


%%
rad2deg(compute_polar_angle_between_vec([1,1,0], [[1,0,0];[1,1,0]], [[0,0,1];[0,0,1]] ) )

%%
rad2deg(compute_polar_angle_between_vec([1,1,1], [[1,0,0];[1,0,1]], [[0,0,1];[0,0,1]] ) )