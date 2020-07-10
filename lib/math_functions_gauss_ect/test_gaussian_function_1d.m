%test_gaussian_function_1d

%% test derivatives
deriv_order=1;
xsamp_conv=linspace(-10,10,1e3+1);
xsamp_conv=col_vec(xsamp_conv);
ngrad_v=derivest(@(x) gaussian_function_1d(x,w_gauss,center),...
    xsamp_conv,'DerivativeOrder',deriv_order,'MaxStep',1);
agrad_v=gaussian_function_1d(xsamp_conv,w_gauss,center,'derivative',deriv_order);

stfig('derivative compare');
subplot(2,1,1)
plot(xsamp_conv,ngrad_v)
hold on
plot(xsamp_conv,agrad_v)
hold off
legend('num','anal')
xlim([-1,1]*2*approx_vfwhm+center)

subplot(2,1,2)
plot(xsamp_conv,ngrad_v-agrad_v)
xlim([-1,1]*2*approx_vfwhm+center)


%% hack to test the 5th derivative
% super slow to evaluate
deriv_order=5;
xsamp_conv=linspace(-10,10,1e3+1);
xsamp_conv=col_vec(xsamp_conv);
ngrad3=@(y) derivest(@(x) gaussian_function_1d(x,w_gauss,center),...
    y,'DerivativeOrder',3,'MaxStep',1);
ngrad_v=derivest(@(x) ngrad3(x),...
    xsamp_conv,'DerivativeOrder',2,'MaxStep',1);
agrad_v=gaussian_function_1d(xsamp_conv,w_gauss,center,'derivative',deriv_order);


stfig('1st derivative');
subplot(2,1,1)
plot(xsamp_conv,ngrad_v)
hold on
plot(xsamp_conv,agrad_v)
hold off
legend('num','anal')
xlim([-1,1]*2*approx_vfwhm+center)

subplot(2,1,2)
plot(xsamp_conv,ngrad_v-agrad_v)
xlim([-1,1]*2*approx_vfwhm+center)



stfig('1st derivative');
subplot(2,1,1)
plot(xsamp_conv,ngrad_v)
hold on
plot(xsamp_conv,agrad_v)
hold off
legend('num','anal')
xlim([-1,1]*2*approx_vfwhm+center)

subplot(2,1,2)
plot(xsamp_conv,ngrad_v-agrad_v)
xlim([-1,1]*2*approx_vfwhm+center)