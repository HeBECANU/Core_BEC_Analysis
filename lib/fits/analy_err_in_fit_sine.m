function unc_est=analy_err_in_fit_sine(in_st)
    % calculate the error in a sine wave fit as derived in
    % "A derivation of the errors for least squares fitting to time series data"
    % http://adsabs.harvard.edu/full/1999DSSN...13...28M
    % https://www.researchgate.net/publication/234204877_A_derivation_of_the_errors_for_least_squares_fitting_to_time_series_data
    
    amp=in_st.amp; % amplitude of the sine wave
    num_samp=in_st.samp_num; % number of samples
    time_samp=in_st.samp_time; %duration over which the oscillation was sampled
    sigma_obs=in_st.sigma_obs; %uncertianty in the observerd variable
    if ~isfield(in_st,'damp_rate')
        in_st.damp_rate=[];
    end
    damp_rate=in_st.damp_rate;

    if ~(isempty(damp_rate)) 
        % we just find the average of 1/amp from t=0 to t=samp_time
        % \frac{1}{T} \int_0^T \frac{1}{\exp (\lambda  (-t))} \, dt
        % = \frac{e^{\lambda  T}-1}{\lambda  T}

        % i have explored other approaches such as weighting the mean, but they all give worse results
        % and generaly overcompensate
        %mean_inv_amp=(1/amp)*(exp(damp_rate.*time_samp)-1)/(damp_rate.*time_samp);
        %inv_amp=mean_inv_amp;
        %wmean_inv_amp=2*(1-exp(-damp_rate.*time_samp))/(-1+exp(2*damp_rate.*time_samp));
        mean_amp_factor=(1-exp(-damp_rate.*time_samp))/(damp_rate.*time_samp);
        %fprintf('amp correction fac %.3f \n',mean_amp_factor)
    else
        mean_amp_factor=1;
    end

    
    sigma_amp=sqrt(2./num_samp).*sigma_obs;
    sigma_freq=sqrt(6./num_samp).*(1/(pi*time_samp)).*sigma_obs*1/(amp*mean_amp_factor);
    sigma_phi=sqrt(2./num_samp).*sigma_obs*1/(amp*mean_amp_factor);

    unc_est=[];
    unc_est.amp=sigma_amp;
    unc_est.freq=sigma_freq;
    unc_est.phase=sigma_phi;
end