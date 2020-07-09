fsamp=linspace(1,1000,1e5)*1e12;
polz_au=aprox_he_polz(fsamp);
stfig('atomic polz');
clf
plot(fsamp*1e-12,polz_au)
ylim([-200,600])
yline(0)


% tune out freq
fminopt = optimset('TolX',1,'TolFun',1e-12);  
tune_out_freq=fminsearch(@(x) abs(aprox_he_polz(x)),600e12,fminopt);  

hold on
plot(tune_out_freq*1e-12,aprox_he_polz(tune_out_freq),'xr')
hold off
fprintf('found tune out at %.6f THz, %.6f nm, polz %g au \n',tune_out_freq*1e-12,f2wl(tune_out_freq)*1e9,aprox_he_polz(tune_out_freq))

%%
derivest(@(x) aprox_he_polz(x),tune_out_freq,'DerivativeOrder',1)

%%

derivest(@(x) aprox_he_polz(x),tune_out_freq,'DerivativeOrder',4)
%%
fminopt = optimset('TolX',1,'TolFun',1e-6);  
tune_out_freq=fzero(@(x) aprox_he_polz(x*1e9),f2wl(413e-9)*1e-9,fminopt);  
tune_out_freq=tune_out_freq*1e-9;
au_polz_at_to=aprox_he_polz(tune_out_freq)
f2wl(tune_out_freq)
hold on
plot(tune_out_freq,au_polz,'xr')
hold off


%% try fiting a linear about the tune out

fsamp=col_vec(linspace(-1,1,1e5))*4e10+tune_out_freq;
[~,polz_si]=aprox_he_polz(fsamp);
stfig('atomic polz');
clf
xscale=1e-9;
yscale=1e44;
plot((fsamp-tune_out_freq)*xscale,polz_si*yscale)
xlabel('$\omega_{\mathrm{TO}}-\omega$ ($2\pi$ GHz)')
ylabel(sprintf('$\\alpha$ ($10^{%g} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',-log10(yscale)))
yline(0)


predictor=(fsamp-tune_out_freq)*1e-9;
response=polz_si*yscale;

fit_in=cat(2,predictor,response,response*nan);
meth_lin_fit=fit_poly_with_int(fit_in,1,0,0);
fprintf('fit intercept %f MHz \n',(meth_lin_fit.x_intercept.val/xscale)*1e-6)
