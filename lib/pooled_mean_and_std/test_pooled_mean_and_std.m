%test_pooled_mean_and_std



sub_samples=5;
sigma_sub=rand(sub_samples,1)*10;
mean_sub=rand(sub_samples,1)*1;
num_obs_sub=randi([100,10000],sub_samples,1);
observations=cell(sub_samples,1);

for ii=1:sub_samples
    observations{ii}=normrnd(mean_sub(ii),sigma_sub(ii),num_obs_sub(ii),1);
end

sub_mean=cellfun(@mean, observations);
sub_std=cellfun(@std, observations);
sub_obs=cellfun(@numel, observations);

[pooled_mean,pooled_std]=pooled_mean_and_std(sub_mean,sub_std,sub_obs);

combined_obs=cat(1,observations{:});
combined_mean=mean(combined_obs);
combined_std=std(combined_obs);


logic_str = {'FAIL','pass'};
test_tol=1e-6;

fprintf('INFO: mean combined  = %.6f\n',combined_mean)
fprintf('INFO: mean pooled    = %.6f\n',pooled_mean)
is_mean_right= frac_diff(pooled_mean,combined_mean)<test_tol;
fprintf('TEST: combined & pooled agree  ±%.0f%% :%s \n',test_tol*1e2,logic_str{is_mean_right+1})
fprintf('INFO: std combined   = %.6f\n',combined_std)
fprintf('INFO: std pooled     = %.6f\n',pooled_std)
is_mean_right= frac_diff(pooled_std,combined_std)<test_tol;
fprintf('TEST: combined & pooled agree  ±%.0fe-6 :%s \n',test_tol*1e6,logic_str{is_mean_right+1})


stfig('combined data');
hist(combined_obs,numel(combined_obs)/10)


%%
%test_pooled_mean_and_std

sub_samples=5;
sigma_sub=rand(sub_samples,1)*10;
mean_sub=rand(sub_samples,1)*1;
num_obs_sub=randi([100,10000],sub_samples,1);
observations=cell(sub_samples,1);

for ii=1:sub_samples
    observations{ii}=mean_sub(ii)+sigma_sub(ii)*rand(num_obs_sub(ii),1);
end

sub_mean=cellfun(@mean, observations);
sub_std=cellfun(@std, observations);
sub_obs=cellfun(@numel, observations);

[pooled_mean,pooled_std]=pooled_mean_and_std(sub_mean,sub_std,sub_obs);

combined_obs=cat(1,observations{:});
combined_mean=mean(combined_obs);
combined_std=std(combined_obs);


logic_str = {'FAIL','pass'};
test_tol=1e-6;

fprintf('INFO: mean combined  = %.6f\n',combined_mean)
fprintf('INFO: mean pooled    = %.6f\n',pooled_mean)
is_mean_right= frac_diff(pooled_mean,combined_mean)<test_tol;
fprintf('TEST: combined & pooled agree  ±%.0fe-6 :%s \n',test_tol*1e6,logic_str{is_mean_right+1})
fprintf('INFO: std combined   = %.6f\n',combined_std)
fprintf('INFO: std pooled     = %.6f\n',pooled_std)
is_mean_right= frac_diff(pooled_std,combined_std)<test_tol;
fprintf('TEST: combined & pooled agree  ±%.0fe-6 :%s \n',test_tol*1e6,logic_str{is_mean_right+1})


stfig('combined data');
hist(combined_obs,numel(combined_obs)/10)


%%
