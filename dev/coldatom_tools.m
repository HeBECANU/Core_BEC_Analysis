% Some useful functions for cold atom calcs

doppler_limit(2*pi*1.6e6)
% doppler_broadening(doppler_limit(2*pi*1.6e6))*f2wl(1083e-9)
% doppler_broadening(doppler_limit(2*pi*1.6e6))*f2wl(410e-9)

T = 100e-6;
doppler_broadening(T)*f2wl(1083e-9)
doppler_broadening(T)*f2wl(410e-9)
% 7.44e14*(cos(15/pi)*recoil_impulse(f2wl(1083.331e-9),6.64e-27)/299792458)

function dv = recoil_impulse(frequency,mass)
    % accepts gamma in rad/sec
    hbar = 1.054571817e-34;
    c = 299792458;
    k = 2*pi*frequency/c;
    dv= hbar*k/mass;
end





% function recoil_limi