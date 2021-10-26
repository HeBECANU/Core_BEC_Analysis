function unc_est=analy_err_in_fit_sine(in_st)
    % calculate the error in a sine wave fit as derived in
    % "A derivation of the errors for least squares fitting to time series data"
    % http://adsabs.harvard.edu/full/1999DSSN...13...28M
    % https://www.researchgate.net/publication/234204877_A_derivation_of_the_errors_for_least_squares_fitting_to_time_series_data
    
    amp=in_st.amp; % amplitude of the sine wave
    num_samp=in_st.samp_num; % number of samples
    time_samp=in_st.samp_time; %duration over which the oscillation was sampled
    sigma_obs=in_st.sigma_obs; %uncertianty in the observerd variable
    
    sigma_amp=sqrt(2./num_samp).*sigma_obs;
    sigma_phi=sqrt(2./num_samp).*sigma_obs/amp;
    sigma_freq=sqrt(6./num_samp).*(1/(pi*time_samp)).*sigma_obs/amp;

    unc_est=[];
    unc_est.amp=sigma_amp;
    unc_est.freq=sigma_freq;
    unc_est.phase=sigma_phi;
end