%test_num_grad_vecs



test_fun=@(x) sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);


num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3)

%%
num_grad_vecs(test_fun,[1,2,3],1e-3)

num_grad_vecs(test_fun,[4,5,6],1e-3)

num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3)

isequal(num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3,1),num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3,0))
