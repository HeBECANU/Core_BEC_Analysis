function deriv=num_grad(fun,x_list,delta)
%finds the 2 point derivative for a n dimensional scalar function
%code is vecorized to handle simultanious evaluation
%delta=delta*(0.5+rand(1));

dim=size(x_list,2);
num_pts=size(x_list,1);
deriv=x_list*NaN;
for i=1:dim
    deriv(:,i)= (fun(x_list+repmat([zeros(1,i-1),delta,zeros(1,dim-i)],num_pts,1))...
               -fun(x_list-repmat([zeros(1,i-1),delta,zeros(1,dim-i)],num_pts,1)))/(2*delta);
end

end