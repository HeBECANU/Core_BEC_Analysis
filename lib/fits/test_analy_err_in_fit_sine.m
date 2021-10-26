% test_analy_err_in_fit_sine
% evaluate this approxmiation for a number of sine waves in different conditions

% define some conditions for the test
osc_param=[];
osc_param.amp=3;
osc_param.freq=1;
osc_param.phase=1;
osc_param.samp_num=1e1;
osc_param.samp_duration=100;
osc_param.samp_err=1;

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

fprintf('amp error pred:%.3e  meas:%.3e  ratio:%.3e \n',param_std_num(1),unc_est.amp, param_std_num(1)/unc_est.amp)
fprintf('amp freq  pred:%.3e  meas:%.3e  ratio:%.3e \n',param_std_num(2),unc_est.freq, param_std_num(2)/unc_est.freq)
fprintf('amp phase pred:%.3e  meas:%.3e  ratio:%.3e \n',param_std_num(3),unc_est.phase, param_std_num(3)/unc_est.phase)





