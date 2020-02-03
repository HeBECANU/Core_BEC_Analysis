%test_logspace_mat

mat_in=magic(5);
a=linspace_mat(mat_in,mat_in+1,10);
isequal(col_vec(squeeze(a(1,1,:))),col_vec(linspace(mat_in(1,1),mat_in(1,1)+1,10)))


%% how it deals with vectors
% does the linspace in the singelton dimension

mat_in=[1,2,3]'
a=linspace_mat(mat_in,mat_in+1,10)


mat_in=[1,2,3]
a=linspace_mat(mat_in,mat_in+1,10)


%% if you turn off the vector fix then it does the linspace in the 3rd dim
% (which i dont see being usefull)
mat_in=[1,2,3]'
a=linspace_mat(mat_in,mat_in+1,10,[],0)
