function T_dop = doppler_limit(gamma)
    % accepts gamma in Hz
    hbar = 1.054571817e-34;
    kB = 1.380649e-23;
    T_dop = 2*pi*hbar*gamma/(2*kB);
end