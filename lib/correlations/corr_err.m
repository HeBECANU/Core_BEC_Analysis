function out = corr_err(corr_func,corr_opts,counts) 
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

out=bootstrap_se(est_fun,counts,...
    'samp_frac_lims',corr_opts.samp_frac_lims,...
    'num_samp_frac',corr_opts.num_samp_frac,...
    'num_samp_rep',corr_opts.num_samp_rep,...
    'replace',false,...
    'use_frac_size',true,...
    'opp_arguments',{corr_opts},...
    'verbose',0);
end

function corr_density = corr_wrap(corr_func,counts,atten,corr_opts)
corr_opts.attenuate_counts = atten;
corr_opts.plots = false;
corr_opts.print_update = false;
out = corr_func(corr_opts,counts);
if isfield(out,'one_d_corr_density')
    corr_density = out.one_d_corr_density;
elseif isfield(out,'rad_corr_density')
    corr_density = out.rad_corr_density;
end
end