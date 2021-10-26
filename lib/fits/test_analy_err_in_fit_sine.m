% test_analy_err_in_fit_sine
% evaluate this approxmiation for a number of sine waves in different conditions

% define some conditions for the test
osc_param=[];
osc_param.amp=3;
osc_param.freq=1;
osc_param.phase=1;
osc_param.samp_num=1e2;
osc_param.samp_duration=50;
osc_param.samp_err=0.1;

num_rep_samp=100; % number of times to repeat the sample and fit
%modelfun=@(b,x) b(1).*sin( (b(2)*2*pi).*x  +b(3)  )+b(4)+b(5).*x;
%names={'amp','freq','phase','offset','grad'};

modelfun=@(b,x) b(1).*sin( (b(2)*2*pi).*x  +b(3)  );
names={'amp','freq','phase'};
beta0=[osc_param.amp,osc_param.freq,osc_param.phase];
t_samp=col_vec(linspace(0,osc_param.samp_duration,osc_param.samp_num));

fit_params_reps.vals=repmat(beta0*nan,num_rep_samp,1);
fit_params_reps.unc=repmat(beta0*nan,num_rep_samp,1);

for ii=1:num_rep_samp
    x_samp=modelfun(beta0,t_samp);
    x_samp=x_samp+normrnd(0,osc_param.samp_err,numel(x_samp),1);
    mdl = fitnlm(t_samp,x_samp,modelfun,beta0,'CoefficientNames',names);
    fit_params_reps.vals(ii,:)=mdl.Coefficients.Estimate;
    fit_params_reps.unc(ii,:)=mdl.Coefficients.SE;
end

param_std_num=std(fit_params_reps.vals,[],1);

in_st=[];
in_st.amp=osc_param.amp;
in_st.samp_num=osc_param.samp_num;
in_st.samp_time=osc_param.samp_duration;
in_st.sigma_obs=osc_param.samp_err;
unc_est=analy_err_in_fit_sine(in_st);

fprintf('amp error   pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.amp,param_std_num(1), unc_est.amp/param_std_num(1))
fprintf('freq error  pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.freq,param_std_num(3), unc_est.freq/param_std_num(3))
fprintf('phase error pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.phase,param_std_num(4), unc_est.phase/param_std_num(4))

% i want to pass if the error is better than a factor of 3 in either direction

ok_factor =2;
if ~is_a_within_factor_of_b(unc_est.amp,param_std_num(1),ok_factor)
    error('predicted amp error was too different to measured')
end
if ~is_a_within_factor_of_b(unc_est.freq,param_std_num(2),ok_factor)
    error('predicted freq error was too different to measured')
end
ok_factor=4;
if ~is_a_within_factor_of_b(unc_est.phase,param_std_num(3),ok_factor)
    error('predicted phase error was too different to measured')
end


%% Now for a damped sine wave
osc_param.amp=10;
osc_param.freq=420;
osc_param.phase=1;
osc_param.damp_rate=1/0.7;
osc_param.samp_num=150;
osc_param.samp_duration=3;
osc_param.samp_err=0.1;

num_rep_samp=100; % number of times to repeat the sample and fit
%modelfun=@(b,x) b(1).*sin( (b(2)*2*pi).*x  +b(3)  )+b(4)+b(5).*x;
%names={'amp','freq','phase','offset','grad'};

modelfun=@(b,x) b(1).*exp(-x.*b(2)).*sin( (b(3)*2*pi).*x  +b(4)  );
names={'amp','damp','freq','phase'};
beta0=[osc_param.amp,osc_param.damp_rate,osc_param.freq,osc_param.phase];
t_samp=col_vec(linspace(0,osc_param.samp_duration,osc_param.samp_num));

fit_params_reps.vals=repmat(beta0*nan,num_rep_samp,1);
fit_params_reps.unc=repmat(beta0*nan,num_rep_samp,1);

for ii=1:num_rep_samp
    x_samp=modelfun(beta0,t_samp);
    x_samp=x_samp+normrnd(0,osc_param.samp_err,numel(x_samp),1);
    mdl = fitnlm(t_samp,x_samp,modelfun,beta0,'CoefficientNames',names);
    fit_params_reps.vals(ii,:)=mdl.Coefficients.Estimate;
    fit_params_reps.unc(ii,:)=mdl.Coefficients.SE;
end

param_std_num=std(fit_params_reps.vals,[],1);

in_st=[];
in_st.amp=osc_param.amp;
in_st.samp_num=osc_param.samp_num;
in_st.samp_time=osc_param.samp_duration;
in_st.sigma_obs=osc_param.samp_err;
unc_est=analy_err_in_fit_sine(in_st);

fprintf('with no damping correction \n')
fprintf('amp error   pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.amp,param_std_num(1), unc_est.amp/param_std_num(1))
fprintf('freq error  pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.freq,param_std_num(3), unc_est.freq/param_std_num(3))
fprintf('phase error pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.phase,param_std_num(4), unc_est.phase/param_std_num(4))


in_st=[];
in_st.amp=osc_param.amp;
in_st.samp_num=osc_param.samp_num;
in_st.samp_time=osc_param.samp_duration;
in_st.sigma_obs=osc_param.samp_err;
in_st.damp_rate=osc_param.damp_rate;
unc_est=analy_err_in_fit_sine(in_st);
fprintf('with  damping correction \n')
fprintf('amp error   pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.amp,param_std_num(1), unc_est.amp/param_std_num(1))
fprintf('freq error  pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.freq,param_std_num(3), unc_est.freq/param_std_num(3))
fprintf('phase error pred:%.3e meas:%.3e    ratio:%.3e \n',unc_est.phase,param_std_num(4), unc_est.phase/param_std_num(4))

% i want to pass if the error is better than a factor of 3 in either direction

ok_factor =2.5;
if ~is_a_within_factor_of_b(unc_est.amp,param_std_num(1),ok_factor)
    error('predicted amp error was too different to measured')
end
if ~is_a_within_factor_of_b(unc_est.freq,param_std_num(2),ok_factor)
    error('predicted freq error was too different to measured')
end
ok_factor=4;
if ~is_a_within_factor_of_b(unc_est.phase,param_std_num(3),ok_factor)
    error('predicted phase error was too different to measured')
end



%% Plot dependence with sampling time


samp_duations=linspace(1,6,30);
iimax=numel(samp_duations);
num_rep_samp=50; % number of times to repeat the sample and fit
scaling_st=[];
scaling_duration=samp_duations;
scaling_uncorrected_ratio_afp=nan(numel(samp_duations),3);
scaling_corrected_ratio_afp=nan(numel(samp_duations),3);

parfor_progress_imp(iimax)
parfor ii=1:iimax
    

    osc_param=[];
    osc_param.amp=10;
    osc_param.freq=420;
    osc_param.phase=1;
    osc_param.damp_rate=1/0.7;
    %osc_param.samp_duration=10;
    osc_param.samp_err=0.1;
    osc_param.samp_duration=samp_duations(ii);
    %samp_rate=10;
    %osc_param.samp_num=samp_rate*osc_param.samp_duration;
    osc_param.samp_num=1e3;

    modelfun=@(b,x) b(1).*exp(-x.*b(2)).*sin( (b(3)*2*pi).*x  +b(4)  );
    names={'amp','damp','freq','phase'};
    beta0=[osc_param.amp,osc_param.damp_rate,osc_param.freq,osc_param.phase];
    t_samp=col_vec(linspace(0,osc_param.samp_duration,osc_param.samp_num));
    
    fit_params_reps=[];
    fit_params_reps.vals=repmat(beta0*nan,num_rep_samp,1);
    fit_params_reps.unc=repmat(beta0*nan,num_rep_samp,1);
    
    for jj=1:num_rep_samp
        x_samp=modelfun(beta0,t_samp);
        x_samp=x_samp+normrnd(0,osc_param.samp_err,numel(x_samp),1);
        mdl = fitnlm(t_samp,x_samp,modelfun,beta0,'CoefficientNames',names);
        fit_params_reps.vals(jj,:)=mdl.Coefficients.Estimate;
        fit_params_reps.unc(jj,:)=mdl.Coefficients.SE;
    end

    %mask=isoutlier(fit_params_reps.vals,"mean",2);
    % set outliers to nan
    %fit_params_reps.vals(mask)=nan;

    param_std_num=nanstd(fit_params_reps.vals,[],1);
    
    in_st=[];
    in_st.amp=osc_param.amp;
    in_st.samp_num=osc_param.samp_num;
    in_st.samp_time=osc_param.samp_duration;
    in_st.sigma_obs=osc_param.samp_err;
    unc_est=analy_err_in_fit_sine(in_st);
    
    scaling_uncorrected_ratio_afp(ii,:)=[unc_est.amp/param_std_num(1),...
                                            unc_est.freq/param_std_num(3),...
                                            unc_est.phase/param_std_num(4)]

    in_st=[];
    in_st.amp=osc_param.amp;
    in_st.samp_num=osc_param.samp_num;
    in_st.samp_time=osc_param.samp_duration;
    in_st.sigma_obs=osc_param.samp_err;
    in_st.damp_rate=osc_param.damp_rate;
    unc_est=analy_err_in_fit_sine(in_st);

    scaling_corrected_ratio_afp(ii,:)=[unc_est.amp/param_std_num(1),...
                                            unc_est.freq/param_std_num(3),...
                                            unc_est.phase/param_std_num(4)]



    parfor_progress_imp;
end
%%
parfor_progress_imp(0)
stfig('deacy correction')
idx=2
subplot(2,1,1)
plot(scaling_duration,scaling_corrected_ratio_afp(:,idx))
xlabel('sample duration')
ylabel('predicted (with corr)/ measured $\sigma(f)$')
%set(gca, 'XScale', 'log')
subplot(2,1,2)
plot(scaling_duration,scaling_uncorrected_ratio_afp(:,idx))
%set(gca, 'XScale', 'log')
xlabel('sample duration')
ylabel('predicted (w/o corr)/ measured $\sigma(f)$')