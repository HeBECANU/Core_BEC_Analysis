function deriv=num_grad_vecs(fun,x_mat,delta,vectorized_logical)
%finds the 2 point derivative of the output of a function with repsect to each input
% of a function which takes a row vector and returns a scalar
% the rows of x_mat specifies the position at which to evaluate the derivatives in each input
% the sperate inputs are stacked in the first dimension
% each row is therefore independent to one another
% if you sepcify vectorized logical the function is assumed to operate independently on rows of the the input matrix X
%       which can be usefull to vectorize the evaluation
% eg test_fun=@(x) sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);

% sometimes this can be usefull to prevent false minima
%delta=delta*(0.5+rand(1));

if nargin<4 || isnan(vectorized_logical) || isempty(vectorized_logical)
    vectorized_logical=0;
end

dim=size(x_mat,2);
num_pts=size(x_mat,1);
deriv=x_mat*NaN;
if vectorized_logical
    for ii=1:dim
        deriv(:,ii)= (fun(x_mat+repmat([zeros(1,ii-1),delta,zeros(1,dim-ii)],num_pts,1))...
                   -fun(x_mat-repmat([zeros(1,ii-1),delta,zeros(1,dim-ii)],num_pts,1)))/(2*delta);
    end
else
    for ii=1:dim
        for jj=1:num_pts
            deriv(jj,ii)= (fun(x_mat(jj,:)+[zeros(1,ii-1),delta,zeros(1,dim-ii)])...
                       -fun(x_mat(jj,:)-[zeros(1,ii-1),delta,zeros(1,dim-ii)]))/(2*delta);
        end
    end
end

end