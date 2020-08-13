function [out, boot_out] = corr_err(corr_func,corr_opts,counts)
warning('off','all')
est_fun = @(c,a,opts) corr_wrap(corr_func,c,a,opts);
if ~isfield(corr_opts,'samp_frac_lims')
    corr_opts.samp_frac_lims = [0.25,0.5];
end
if ~isfield(corr_opts,'num_samp_frac')
    corr_opts.num_samp_frac = 5;
end
if ~isfield(corr_opts,'num_samp_rep')
    corr_opts.num_samp_rep = 3e1;
end

boot_out=bootstrap_se(est_fun,counts,...
    'samp_frac_lims',corr_opts.samp_frac_lims,...
    'num_samp_frac',corr_opts.num_samp_frac,...
    'num_samp_rep',corr_opts.num_samp_rep,...
    'replace',false,...
    'use_frac_size',true,...
    'opp_arguments',{corr_opts},...
    'verbose',0);

grid_size = size(corr_opts.redges,2) - 1;
out.rawcorr.val = nansum(boot_out.sampling.mean(:,1:grid_size).*boot_out.sampling.sample_size)./nansum(boot_out.sampling.sample_size);
out.normcorr.val = nansum(boot_out.sampling.mean(:,(grid_size+1):(2*grid_size)).*boot_out.sampling.sample_size)./nansum(boot_out.sampling.sample_size);
out.corr_density.val = nansum(boot_out.sampling.mean(:,(2*grid_size+1):(3*grid_size)).*boot_out.sampling.sample_size)./nansum(boot_out.sampling.sample_size);
out.fit.val = nansum(boot_out.sampling.mean(:,end-1:end).*boot_out.sampling.sample_size)./nansum(boot_out.sampling.sample_size);

out.rawcorr.unc = boot_out.results.se_fun_whole_unweighted(1,1:grid_size);
out.normcorr.unc =boot_out.results.se_fun_whole_unweighted(1,(grid_size+1:(2*grid_size)));
out.corr_density.unc = boot_out.results.se_fun_whole_unweighted(1,(2*grid_size+1):(3*grid_size));
out.fit.unc = boot_out.results.se_fun_whole_unweighted(1,end-1:end);

warning('on','all')
end

function out = corr_wrap(corr_func,counts,atten,corr_opts)
corr_opts.attenuate_counts = atten;
corr_opts.plots = false;
corr_opts.print_update = false;
rawcorr = corr_func(corr_opts,counts);

norm_sort_dir=corr_opts.sorted_dir;
if size(counts,1)>1
    %set the number of chunks to be at least as many as the heighest
    %count number out of all halos
    chunk_num = max([cellfun(@(x)size(x,1),counts(1,:)),cellfun(@(x)size(x,1),counts(2,:))]);
    counts_chunked(1,:)=chunk_data(counts(1,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
    counts_chunked(2,:)=chunk_data(counts(2,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
else
    counts_chunked=chunk_data(counts,corr_opts.norm_samp_factor,norm_sort_dir);
end
normcorr=corr_func(corr_opts,counts_chunked);



if isfield(rawcorr,'one_d_corr_density')
    density='one_d_corr_density';
    centers='x_centers';
elseif isfield(rawcorr,'rad_corr_density')
    density='rad_corr_density';
    centers='rad_centers';
end

corr_density = rawcorr.(density)./normcorr.(density);

fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);

is_data_flat = isdataflat(corr_density,0.2);
if is_data_flat
    fun1d =  @(b,x)(b(1)+1).*ones(size(x));
    inital_guess=nanmean(corr_density)-1;
else
    fun1d =  @(b,x) b(1).*exp(-((x).^2)./(2*b(2).^2))+1;
    [~,sigmaHat] = normfit(rawcorr.(centers)',0.01,zeros(size(rawcorr.(centers)')),abs(corr_density-1).^2');
    inital_guess=[max(corr_density)-1,sigmaHat];
end

fit=fitnlm(rawcorr.(centers),corr_density,...
    fun1d,...
    inital_guess,...
    'Options',fo);

if fit.Rsquared.Adjusted>=0
    if is_data_flat
        b = [fit.Coefficients.Estimate; nan]; %fitted parameters
    else
        b = fit.Coefficients.Estimate; %fitted parameters
    end
else
    b = [nan; nan];
end

b(2) = abs(b(2));

out = [rawcorr.(density);normcorr.(density);corr_density;b];

end