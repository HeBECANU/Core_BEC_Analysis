%test_num_grad_tensors



test_fun=@(x) sum(x(:).^2);

x_in=rand(3,2,4);
num_grad_tensors(test_fun,x_in,0.1)

%%
x_in=rand(1,1);
df_out=num_grad_tensors(test_fun,x_in,0.1)
isalmost(df_out,2*x_in,eps(10))

%%
x_in=rand(2,3);
[df_out,mean_out]=num_grad_tensors(test_fun,x_in,0.1)
diff(mean_out,1,2)
isalmost(df_out,2*x_in,eps(100))



%%

test_fun=@(x) [[1,2];[3,4]]*sum(x(:).^2);

x_in=rand(3,2,4);
num_grad_tensors(test_fun,x_in,0.1)