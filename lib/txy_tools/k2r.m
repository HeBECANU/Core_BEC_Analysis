function [r,T] = k2r(varargin)
    % Rough - correct in horizontal ballistic picture (does not account for
    % gravity) for estimating only, generally. 
    % Use txy_to_vel for more precise conversion.
    
    % Converts a distance (in m) to a wavenumber (m^-1) after a falltime
    % (default 417.6 ms) for particle of mass m (default 6.64e-27 kg)
    % Also returns the time separation in the case that particles are
    % moving purely horizontally
    % 
    % Syntax: [r,t] = k2r(k)
    %         r = k2r(k,m)
    % Test: k2r(r2k(1)) == 1
    %            k2r(1) == 1.6069e8
    
    t = 0.4176; %fall time 
    k = varargin{1};
    h = 6.63e-34;
    v_c = 4.9*t^2; % center of mass velocity
    hbar = h/(2*pi);
    if nargin == 1 % just passed distance, no mass
        m = 6.64e-27; %kg
    else
        m = varargin{2};
    end
   
    r = hbar*k*t/m;
    T = r/v_c;
end
    