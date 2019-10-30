%test_angle_between_vec

u=[[1,0,0];[0,1,0];[0,0,1];[1,1,0]];
v=[[1,0,0];[1,0,0];[1,0,0];[1,0,0]];

angles=angle_between_vec(u,v)/pi
isequal(angles,col_vec([0,0.5,0.5,0.25]))