%test_std_err_of_sample_std


n_samp=30;
data_sd=2;
test_itteration=1000;
normdat=normrnd(0,data_sd,n_samp,test_itteration);

stds_est=std(normdat,[],1);
sim_se_sd=std(stds_est);

pred_se_sd=std_err_of_sample_std(n_samp,data_sd);
fprintf('predicted se in sd %.3f\n',pred_se_sd)
fprintf('sim se in sd       %.3f\n',sim_se_sd)
if abs(frac_diff(sim_se_sd,sim_se_sd))>0.2
    error('did not give the right prediction')
end