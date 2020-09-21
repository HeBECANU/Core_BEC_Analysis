function lambda = deBroglie_wavelength(temperature, varargin)
%  Returns the de Broglie wavelength at temperature T (in kelvin). 
%  Returns inf if passed T=0.
% If only one argument is passed, assumes 4He mass. Eg
% ```
%  lambda_db(300) == 5.0367e-11
% ```
% Accepts particle mass as second argument, for example  3He
% ```
%  lambda_db(300,4.9825e-27) == 5.8158e-11

% ```
    if nargin>1
        mass = varargin{1};
    else
        mass = 6.64332e-27;
    end
    lambda = sqrt((2*pi*(1.0540e-34)^2) ./ (mass * 1.3806e-23 *temperature));

end