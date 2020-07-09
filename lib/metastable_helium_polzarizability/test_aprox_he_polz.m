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


