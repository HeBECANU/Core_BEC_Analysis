%test_num_grad_vecs



test_fun=@(x) sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);


num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3)

%%
num_grad_vecs(test_fun,[1,2,3],1e-3)

num_grad_vecs(test_fun,[4,5,6],1e-3)

num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3)

isequal(num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3,1),num_grad_vecs(test_fun,[[1,2,3];[4,5,6]],1e-3,0))


%% basic test fun
grad=10;
f_simple=@(x) sum(x*grad);
num_grad=num_grad_vecs(f_simple,[1,2,3],1.324e-6,0);
grad_err=num_grad-grad;
if max(grad_err)>1e-7
    error('not right rms error %.2e',rms(grad_err))
end
