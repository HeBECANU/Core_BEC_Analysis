%test_taylor_series
xvals=linspace(-1,1,1e3)';
yvals_tay=taylor_series(xvals,ones(1,100),0);
yvals_direct=exp(xvals);
diff=yvals_tay-yvals_direct;
plot(xvals,diff)
max(diff)

%% test taking the derivative

derivs=randn(10,1);
offset=1;
num_deriv=@(x,fun,h) (fun(x + h) - fun(x - h)) ./ (2*h);
xvals=col_vec(linspace(-1,1,1e3));
%pick the optimal derivative size step
%https://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating-point_arithmetic
num_2deriv=@(x,fun,h) num_deriv(x,@(y) num_deriv(y,fun,h),h);
dxdxzero=num_2deriv(0,@(x) taylor_series(x,derivs,offset),1e-6);
optimal_deriv_step=2*sqrt(eps*abs(taylor_series(0,derivs,offset)/dxdxzero));
optimal_deriv_step=1e-5;
num_deriv_vals=num_deriv(xvals,@(x) taylor_series(x,derivs,offset),optimal_deriv_step);
anlytic_deriv_vals=deriv_taylor_series(xvals,derivs,offset,1);

stfig('test analytic derivatives')
clf
plot(xvals,num_deriv_vals,'k')
hold on
plot(xvals,anlytic_deriv_vals,'r')
hold off
max(abs(num_deriv_vals-anlytic_deriv_vals))


%% test taking the second derivative
optimal_deriv_step=1e-4;
num_deriv_vals=num_2deriv(xvals,@(x) taylor_series(x,derivs,offset),optimal_deriv_step);
anlytic_deriv_vals=deriv_taylor_series(xvals,derivs,offset,2);

stfig('test analytic derivatives')
clf
plot(xvals,num_deriv_vals,'k')
hold on
plot(xvals,anlytic_deriv_vals,'r')
hold off
max(abs(num_deriv_vals-anlytic_deriv_vals))

