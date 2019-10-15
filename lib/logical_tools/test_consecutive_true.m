%test_consecutive_true

in_vec=[0,0,0,1,1,1,0,0,1,0,1,1,1];
out_vec=consecutive_true(in_vec,3);

fprintf('%u ',in_vec)
fprintf('\n')
fprintf('%u ',out_vec)
fprintf('\n')

%%
in_vec=[1 1 0 0 1 1];
out_vec=consecutive_true(in_vec,3);
fprintf('%u ',in_vec)
fprintf('\n')
fprintf('%u ',out_vec)
fprintf('\n')

%%

in_vec=[1 1 0 0 1 1 1 1];
out_vec=consecutive_true(in_vec,3);
fprintf('%u ',in_vec)
fprintf('\n')
fprintf('%u ',out_vec)
fprintf('\n')


%%
in_vec=[1 1 0 0 1 1];
out_vec=consecutive_true(in_vec,3,1,1);
fprintf('%u ',in_vec)
fprintf('\n')
fprintf('%u ',out_vec)
fprintf('\n')



