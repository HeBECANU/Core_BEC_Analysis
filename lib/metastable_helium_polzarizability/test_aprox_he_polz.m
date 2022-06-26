fsamp=linspace(276.4,277,1e3)*1e12;


[polz_au,polz_si,higher_terms]=aprox_he_polz(fsamp);
stfig('polz au');
clf
plot(fsamp*1e-12,polz_au)
yrange=1e6;
ylim([-1,1]*yrange)
yline(0)
xlim([fsamp(1),fsamp(end)]*1e-12)
xlabel('Frequency, $f$ (THz)')
ylabel('Polarizability, $\alpha$ (atom. u.)')

stfig('polz si');
plot(fsamp*1e-12,polz_si)
hold on
plot(fsamp*1e-12,higher_terms.si_polz.scalar)
hold off
xlim([fsamp(1),fsamp(end)]*1e-12)
yrange=1e-34;
%ylim([-1,1]*yrange)
xlabel('Frequency, $f$ (THz)')
ylabel('Polarizability, $\alpha$ (SI)')

stfig('polz_vec si');
plot(fsamp*1e-12,higher_terms.si_polz.scalar)
hold on
plot(fsamp*1e-12,higher_terms.si_polz.vector)
plot(fsamp*1e-12,higher_terms.si_polz.tensor)
hold off
xlim([fsamp(1),fsamp(end)]*1e-12)
yrange=1e-34;
%ylim([-1,1]*yrange)
xlabel('Frequency, $f$ (THz)')
ylabel('Polarizability, $\alpha$ (SI)')

%%

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
plot(tune_out_freq,polz_au,'xr')
hold off


