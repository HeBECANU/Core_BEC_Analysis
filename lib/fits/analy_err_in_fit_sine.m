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
    damp_threshold=1e-3;
    if  (~isempty(damp_rate) && damp_rate*time_samp> damp_threshold)
        %sigma_amp=sigma_obs*sqrt(2/num_samp)*...
        %    sqrt(time_samp*damp_rate.*(1+coth(time_samp*damp_rate)));

        sigma_amp=sigma_obs*sqrt(2/num_samp)*...
            sqrt(time_samp*damp_rate.*(1/16+coth(time_samp*damp_rate)));
        sigma_amp=sigma_amp*2;  % emp factor

        sigma_freq=(2*sigma_obs)./(amp*sqrt(pi*num_samp)).*    ...
        sqrt(...
            (time_samp*damp_rate^3).*(-1+exp(2*time_samp*damp_rate))/...
            (-1-2*(time_samp*damp_rate)^2+cosh(2*time_samp*damp_rate)) ...
            );
        sigma_freq=sigma_freq*(2/5); % emp factor this seems to be exact

        % sigma_freq=1.0921*(sigma_obs*damp_rate)./(amp*sqrt(num_samp));


        sigma_phi=sqrt(4./num_samp).*sigma_obs/(amp).*...
            sqrt(...
            (  (time_samp*damp_rate).*(-1+exp(2*time_samp*damp_rate)...
                            -2.*time_samp*damp_rate.*(1+time_samp*damp_rate)))./...
            (-1-2.*(time_samp*damp_rate)^2+cosh(2.*time_samp*damp_rate)) ...
            );


    else
        % the normal formula

        sigma_amp=sqrt(2./num_samp).*sigma_obs;
        sigma_amp=sigma_amp*2;  % emp factor
        sigma_freq=sqrt(6/num_samp).*(1/(pi*time_samp))...
            .*(sigma_obs/amp);

        sigma_phi=sqrt(2./num_samp).*sigma_obs/(amp);
    end

    

    unc_est=[];
    unc_est.amp=sigma_amp;
    unc_est.freq=sigma_freq;
    unc_est.phase=sigma_phi;
end