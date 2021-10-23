% voigt git play

w_gauss=0.1;
w_lorz=0.2;
center=0.0;

xsamp_comp=linspace(-10,10,1e4+1); %odd number of bins allows better results from conv
xsamp_comp=col_vec(xsamp_comp);
dx_comp=abs(xsamp_comp(2)-xsamp_comp(1));

gauss_component=gaussian_function_1d(xsamp_comp,w_gauss,center,'norm','int');
trapz(dx_comp,gauss_component)
lorentzian_component=lorentzian_function_1d(xsamp_comp,w_lorz,center,'norm','int');
trapz(dx_comp,lorentzian_component)

stfig('plot componets');
clf
subplot(2,1,1)
plot(xsamp_comp,gauss_component)
hold on
plot(xsamp_comp,lorentzian_component)

conv_comp=conv(gauss_component,lorentzian_component,'same')*dx_comp;

xsamp_conv=xsamp_comp(1)+dx_comp*(0:(numel(conv_comp)-1));
xsamp_conv=col_vec(xsamp_conv)

plot(xsamp_conv,conv_comp)

tic
vo_approx=voigt_function_1d(xsamp_conv,w_gauss,w_lorz,center,'method','approx');
toc
plot(xsamp_conv,vo_approx,'--')


vo_fadeva=voigt_function_1d(xsamp_conv,w_gauss,w_lorz,center,'method','fadd');

plot(xsamp_conv,vo_fadeva,'-.')

legend('gauss','lorz','conv','approx',' fadeva')
hold off

approx_vfwhm=voigt_approx_fwhm(w_gauss,w_lorz);

xlim([-1,1]*1*approx_vfwhm+center)


subplot(2,1,2)
plot(xsamp_conv,conv_comp-vo_approx)
hold on
plot(xsamp_conv,conv_comp-vo_fadeva,'--')
hold off
xlim([-1,1]*1*approx_vfwhm+center)


%% test derivatives
xsamp_conv=linspace(-10,10,1e3+1);
xsamp_conv=col_vec(xsamp_conv);
ngrad_v=derivest(@(x) voigt_function_1d(x,w_gauss,w_lorz,center,'method','fadd'),...
    xsamp_conv,'DerivativeOrder',1)
agrad_v=voigt_function_1d(xsamp_conv,w_gauss,w_lorz,center,'method','fadd','derivative',1);

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



%% 2nd order deriv
%ngrad2_v=gradient(ngrad_v,spacing);
xsamp_conv=linspace(-10,10,1e3+1);
xsamp_conv=col_vec(xsamp_conv);
ngrad2_v=derivest(@(x) voigt_function_1d(x,w_gauss,w_lorz,center,'method','fadd'),...
    xsamp_conv,'DerivativeOrder',2)
agrad2_v=voigt_function_1d(xsamp_conv,w_gauss,w_lorz,center,'method','fadd','derivative',2);

stfig('2nd derivative');
subplot(2,1,1)
plot(xsamp_conv,ngrad2_v)
hold on
plot(xsamp_conv,agrad2_v)
hold off
legend('num','anal')
xlim([-1,1]*2*approx_vfwhm+center)

subplot(2,1,2)
plot(xsamp_conv,agrad2_v-ngrad2_v)
xlim([-1,1]*2*approx_vfwhm+center)



%% Test lorentzian function
xsamp=linspace(-10,10,1e3+1);
xsamp=col_vec(xsamp);
order=1;
nderiv=derivest(@(x) lorentzian_function_1d(x,w_lorz,center,'norm','int'),...
    xsamp,'DerivativeOrder',order)
aderiv=lorentzian_function_1d(xsamp,w_lorz,center,'norm','int','derivative',order);

stfig('2nd derivative');
subplot(2,1,1)
plot(xsamp,nderiv)
hold on
plot(xsamp,aderiv)
hold off
legend('num','anal')
xlim([-1,1]*2*approx_vfwhm+center)

subplot(2,1,2)
plot(xsamp,aderiv-nderiv)
xlim([-1,1]*2*approx_vfwhm+center)

