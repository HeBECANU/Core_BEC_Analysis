function k = r2k(varargin)
    % Rough - correct in horizontal ballistic picture (does not account for
    % gravity) for estimating only, generally. 
    % Use txy_to_vel for more precise conversion.
    
    % Converts a distance (in m) to a wavenumber (m^-1) after a falltime
    % (default 417.6 ms) for particle of mass m (default 6.64e-27 kg)
    % 
    % Syntax: k = r2k(r)
    %         k = r2k(k,m)
    % Test: k2r(r2k(1)) == 1
    %            r2k(4e-2) == 6.0275e+06
    
    t = 0.4176; %fall time 
    r = varargin{1};
    h = 6.63e-34;
    hbar = h/(2*pi);
    if nargin == 1 % just passed distance, no mass
        m = 6.64e-27; %kg
    else
        m = varargin{2};
    end
   
    k = m*r/(hbar*t);
    
end
    