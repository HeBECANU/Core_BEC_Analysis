%% whats the  bandwitdh of a gaussian filter
% im sure i could do this analyticaly but i aint got time foo dat
% Bryce Henson 2022-02-20
% turns out sigma*BW(-3db) =  0.1876

filt_sigma=1;
fit_bw_guess=1/filt_sigma;
num_samp_cycles=1e3;
tmax=num_samp_cycles*(1/fit_bw_guess);
xsamp=linspace(0,tmax,num_samp_cycles*100);
filt_bw_pow=fzero(@(f) rms(gaussfilt(xsamp,sin(xsamp*f*2*pi),filt_sigma))-0.5*1/sqrt(2),fit_bw_guess);
filt_bw_pow
%the product of the time sigma and the 3db baddwidth
sigma_bw_product=filt_bw_pow*filt_sigma

%%
stfig('filtered');
clf
plot(xsamp,gaussfilt(xsamp,sin(xsamp*filt_bw_pow*2*pi),filt_sigma))
hold on
plot(xsamp,sin(xsamp*filt_bw_pow*2*pi))
hold off