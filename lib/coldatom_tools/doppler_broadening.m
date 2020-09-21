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