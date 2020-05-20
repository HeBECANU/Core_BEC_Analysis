function fwhm_v=voigt_approx_fwhm(sigma,gamma)
% using approx from 
% https://doi.org/10.1016%2F0022-4073%2877%2990161-3
% Bryce Henson 2020-04-22

fwhm_l=2*gamma;
fwhm_g=2*sigma*sqrt(2*log(2));
fwhm_v=0.5346*fwhm_l+sqrt(0.2166*(fwhm_l^2)+fwhm_g^2);

end

