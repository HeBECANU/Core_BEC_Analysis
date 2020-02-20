%test_rand_interval

out_interval=[1,2];
test_out=rand_interval(out_interval,[1,10])
if sum(test_out>out_interval(2))>0 || sum(test_out<out_interval(1))>0
    error('outside limits')
end



%% constructing the limits matrix is the hardest part of using this function
out_size=[1,2];
out_interval=cat(2,cat(3,2,3), cat(3,1,2))
out_interval=reshape([[2,3];[1,2]],cat(2,out_size,2))
test_out=rand_interval(out_interval,out_size)


%% in higher dimensions it can get a bit painfull but carefully used reshape makes it easier

out_interval=cat(1,cat(2,cat(3,1,2), cat(3,2,3)),cat(2,cat(3,3,4), cat(3,4,5)));
out_size=[2,2];
out_interval=reshape([[1,2];[2,3];[3,4];[4,5]],cat(2,out_size,2));
test_out=rand_interval(out_interval,out_size)

%% you dont even need to pass the output size if you are specifying the limits for each element
out_size=[2,2];
out_interval=reshape([[1,2];[2,3];[3,4];[4,5]],cat(2,out_size,2));
test_out=rand_interval(out_interval)



%% the function will check that the limits are sensible and that max and min are not equal (will produce a warning)

out_size=[2,2];
out_interval=reshape([[1,1];[2,3];[3,4];[4,5]],cat(2,out_size,2));
test_out=rand_interval(out_interval,out_size)

%% check that the test for min greater than max works (should produce a error)
out_size=[2,2];
out_interval=reshape([[2,1];[2,3];[3,4];[4,5]],cat(2,out_size,2));
test_out=rand_interval(out_interval,out_size)


