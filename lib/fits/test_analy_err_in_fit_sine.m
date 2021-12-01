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
 osc_param=[];
osc_param.amp=10;
osc_param.freq=10;
osc_param.phase=1;
osc_param.damp_rate=1/10;
%osc_param.samp_duration=10;
osc_param.samp_err=0.5;
osc_param.samp_duration=samp_durations(ii);
osc_param.samp_rate=10*osc_param.freq;
%osc_param.samp_num=1e4;

%samp_durations=col_vec(linspace(0.05,4,60));
samp_durations=col_vec(logspace(log10(2/osc_param.freq),log10(1000),1e3));
num_rep_samp=20; % number of times to repeat the sample and fit
do_fit_plots=false;


iimax=numel(samp_durations);
% store the unc values as a amplitude, frequency, phase vector
sim_unc_afp=nan(numel(samp_durations),4); % this incudes damping uncert
sim_unc_unc_afp=nan(numel(samp_durations),4);
uncorrected_unc_afp=nan(numel(samp_durations),3);
corrected_unc_afp=nan(numel(samp_durations),3);
samp_numbers=nan(numel(samp_durations),1);

parfor_progress_imp(iimax)
parfor ii=1:iimax
 
    osc_param_loc=osc_param;
    osc_param_loc.samp_duration=samp_durations(ii);
    if isfield(osc_param_loc,'samp_rate') && isfield(osc_param_loc,'samp_num')
        error('spec samp rate or number not both')
    end
    if isfield(osc_param_loc,'samp_rate')
        osc_param_loc.samp_num=osc_param_loc.samp_rate*osc_param_loc.samp_duration;
    end
    
    modelfun=@(b,x) b(1).*exp(-x.*b(4)).*sin( (b(2)*2*pi).*x  +b(3)  );
    names={'amp','freq','phase','damp'};
    beta0=[osc_param_loc.amp,osc_param_loc.freq,osc_param_loc.phase,osc_param_loc.damp_rate];
    t_samp=col_vec(linspace(0,osc_param_loc.samp_duration,osc_param_loc.samp_num));
    
    fit_params_reps=[];
    fit_params_reps.vals=repmat(beta0*nan,num_rep_samp,1);
    fit_params_reps.unc=repmat(beta0*nan,num_rep_samp,1);
    
    for jj=1:num_rep_samp
        x_samp=modelfun(beta0,t_samp);
        x_unc=osc_param_loc.samp_err;
        x_samp=x_samp+normrnd(0,x_unc,numel(x_samp),1);
        mdl = fitnlm(t_samp,x_samp,modelfun,beta0,'CoefficientNames',names);
        fit_params_reps.vals(jj,:)=mdl.Coefficients.Estimate;
        fit_params_reps.unc(jj,:)=mdl.Coefficients.SE;

        if do_fit_plots
            colors_main=[[255,0,0];[33,188,44];[0,0,0]]./255;
            lch=colorspace('RGB->LCH',colors_main(:,:));
            lch(:,1)=lch(:,1)+20;
            colors_detail=colorspace('LCH->RGB',lch);
            %would prefer to use srgb_2_Jab here
            color_shaded=colorspace('RGB->LCH',colors_main(3,:));
            color_shaded(1)=50;
            color_shaded=colorspace('LCH->RGB',color_shaded);
    
            stfig('fit example');
            clf
            hold on
            predictor_val=col_vec(linspace(min(t_samp),max(t_samp),1e5));
            [pred_val,pred_ci]=predict(mdl,predictor_val,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
            plot(predictor_val,pred_ci(:,1),'-','LineWidth',1.5,'Color',color_shaded)
            plot(predictor_val,pred_ci(:,2),'-','LineWidth',1.5,'Color',color_shaded)
            plot(predictor_val,pred_val,'-','LineWidth',1.0,'Color',colors_main(3,:))
            x_unc=(x_samp*0+1)*x_unc;
            errorbar(t_samp,x_samp,x_unc,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5)
            hold off
            xlabel('time')
            ylabel('val')
            pause
        end
    end

    mask=isoutlier(fit_params_reps.vals,"mean",2);
    % set outliers to nan
    fit_params_reps.vals(mask)=nan;

    sim_unc_afp(ii,:)=nanstd(fit_params_reps.vals,[],1);
    % calculate the 
    sim_unc_unc_afp(ii,:)=std_err_of_sample_std(sum(~mask,1),sim_unc_afp(ii,:));



    parfor_progress_imp;
end
parfor_progress_imp(0)

%% calculate the predicted in a seprate loop so i dont have to rerun the simulations when
% trying new corrections out

for ii=1:iimax

    osc_param_loc=osc_param;
    osc_param_loc.samp_duration=samp_durations(ii);
    if isfield(osc_param_loc,'samp_rate') && isfield(osc_param_loc,'samp_num')
        error('spec samp rate or number not both')
    end
    if isfield(osc_param_loc,'samp_rate')
        osc_param_loc.samp_num=osc_param_loc.samp_rate*osc_param_loc.samp_duration;
    end
    samp_numbers(ii)=osc_param_loc.samp_num;

    in_st=[];
    in_st.amp=osc_param_loc.amp;
    in_st.samp_num=osc_param_loc.samp_num;
    in_st.samp_time=osc_param_loc.samp_duration;
    in_st.sigma_obs=osc_param_loc.samp_err;
    unc_est=analy_err_in_fit_sine(in_st);
    
    uncorrected_unc_afp(ii,:)=[unc_est.amp,unc_est.freq,unc_est.phase];
    
    in_st=[];
    in_st.amp=osc_param_loc.amp;
    in_st.samp_num=osc_param_loc.samp_num;
    in_st.samp_time=osc_param_loc.samp_duration;
    in_st.sigma_obs=osc_param_loc.samp_err;
    in_st.damp_rate=osc_param_loc.damp_rate;
    unc_est=analy_err_in_fit_sine(in_st);
    
    corrected_unc_afp(ii,:)=[unc_est.amp,unc_est.freq,unc_est.phase];

end


ratio_uncorrected_unc_afp=uncorrected_unc_afp./sim_unc_afp(:,1:3);
ratio_corrected_unc_afp=corrected_unc_afp./sim_unc_afp(:,1:3);


% predicted sampled time optima
if ~isfield(osc_param_loc,'samp_num')
    uncfopt_time=2.0175/osc_param.damp_rate;
    uncfopt_val=(1.092)*osc_param.damp_rate*osc_param.samp_err...
                /(osc_param.amp*sqrt(osc_param.samp_num));
else
    uncfopt_time=nan;
    uncfopt_val=nan;
end


colors_main=[[255,0,0];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);

stfig('fit unc dep');
idx=1;
subplot(3,1,idx)
plot(samp_durations,sim_unc_afp(:,idx),'k')
%errorbar(samp_durations,sim_unc_afp(:,idx),sim_unc_unc_afp(:,idx),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5)
hold on
plot(samp_durations,uncorrected_unc_afp(:,idx),'r')
plot(samp_durations,corrected_unc_afp(:,idx),'b')
hold off
legend('sim','uncorr','corr')
xlabel('sample duration')
ylabel('fit amp unc')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
idx=2;
subplot(3,1,idx)
plot(samp_durations,sim_unc_afp(:,idx),'k')
%errorbar(samp_durations,sim_unc_afp(:,idx),sim_unc_unc_afp(:,idx),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5)
hold on
plot(samp_durations,uncorrected_unc_afp(:,idx),'r')
plot(samp_durations,corrected_unc_afp(:,idx),'b')

plot(uncfopt_time,uncfopt_val,'o')
hold off
legend('sim','uncorr','corr')
xlabel('sample duration')
ylabel('fit freq unc (Hz)')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
idx=3;
subplot(3,1,idx)
plot(samp_durations,sim_unc_afp(:,idx),'k')
%errorbar(samp_durations,sim_unc_afp(:,idx),sim_unc_unc_afp(:,idx),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5)
hold on
plot(samp_durations,uncorrected_unc_afp(:,idx),'r')
plot(samp_durations,corrected_unc_afp(:,idx),'b')
hold off
legend('sim','uncorr','corr')
xlabel('sample duration')
ylabel('fit phase unc (rad)')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

mean_corr_ratio=mean(ratio_corrected_unc_afp);
fprintf('mean ratios \n')
fprintf('%g, ', mean_corr_ratio);
fprintf('\n')
stfig('fit ratios');
clf;
idx=1;
subplot(3,1,idx)
hold on
plot(samp_durations,ratio_uncorrected_unc_afp(:,idx),'r')
plot(samp_durations,ratio_corrected_unc_afp(:,idx),'b')
yline(1,'k')
hold off
legend('uncorr','corr')
xlabel('sample duration')
ylabel('fit amp unc ratio')
set(gca, 'XScale', 'log')
idx=2;
subplot(3,1,idx)
hold on
plot(samp_durations,ratio_uncorrected_unc_afp(:,idx),'r')
plot(samp_durations,ratio_corrected_unc_afp(:,idx),'b')
yline(1,'k')
hold off
legend('uncorr','corr')
xlabel('sample duration')
ylabel('fit freq unc ratio')
set(gca, 'XScale', 'log')
idx=3;
subplot(3,1,idx)
hold on
plot(samp_durations,ratio_uncorrected_unc_afp(:,idx),'r')
plot(samp_durations,ratio_corrected_unc_afp(:,idx),'b')
yline(1,'k')
hold off
legend('uncorr','corr')
xlabel('sample duration')
ylabel('fit phase unc ratio')
set(gca, 'XScale', 'log')

%% lets calculate a empirical correction to the formula so that it is equal

predictor=samp_durations;
response=sim_unc_afp(:,2);

%corr_model=@(b,x) (b(1)+b(2)./(b(3)+x)).*(corrected_unc_afp(:,2));
damp_rate=osc_param.damp_rate;
mean_exp=@(time,lambda) (1-exp(-time.*lambda))./(time.*lambda) ;
corr_model=@(b,x) (1./mean_exp(x,damp_rate))...
                    .*(sqrt(samp_numbers)./sqrt((samp_numbers./x).*(1./damp_rate)))...
                    .*(uncorrected_unc_afp(:,2));
names={'offset','grad','quad','thing2'};
beta0=[1,-3,1,1];
mdl = fitnlm(predictor,response,corr_model,beta0,'CoefficientNames',names)
%mdl.Coefficients.Estimate;
%mdl.Coefficients.SE;


colors_main=[[255,0,0];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);

stfig('fit example');
clf
hold on
predictor_val=predictor;
[pred_val,pred_ci]=predict(mdl,predictor_val,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
plot(predictor_val,pred_ci(:,1),'-','LineWidth',1.5,'Color',color_shaded)
plot(predictor_val,pred_ci(:,2),'-','LineWidth',1.5,'Color',color_shaded)
plot(predictor_val,pred_val,'-','LineWidth',1.0,'Color',colors_main(3,:))
response_unc=sim_unc_unc_afp(:,2);
errorbar(predictor,response,response_unc,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5)
hold off
xlabel('time')
ylabel('val')
%set(gca,'XScale','log')
%set(gca,'YScale','log')   

