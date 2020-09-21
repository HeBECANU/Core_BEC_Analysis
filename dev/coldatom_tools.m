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

function T_dop = doppler_limit(gamma)
    % accepts gamma in rad/sec
    hbar = 1.054571817e-34;
    kB = 1.380649e-23;
    T_dop = hbar*gamma/(2*kB);
end

function FWHM = doppler_broadening(T,varargin)
    if nargin < 2
        m = 6.6433e-27;
    else
        m = varargin{1};
    end
    kB = 1.380649e-23;
    c = 299792458;
    % the multiplicative factor to obtain FWHM of a delta-function distribution 
    % at frequency f0 for atoms of mass m at temperature T.
    % eg doppler_width = doppler_broadening(T,m)*f0
    FWHM = sqrt(8*kB*T*log(2)/(m*c^2));
end

% function recoil_limi