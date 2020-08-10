function out=calc_any_g2_type(corr_opts,counts)
%caucluates normalized g2 functions for a lot of different cases
%uses a lot of correlators from https://github.com/spicydonkey/correlation-funs (with some modifications)
%input
%counts_txy- cell array (n_shots*[1 or 2]) of txy data (n counts * 3)
%            txy must be sorted in time for prewindowing option
%            cell array can be n_shots*2 for two port xcorr
%corr_opts.type- what type of g2 fo you want to calculate
%               '1d_cart_cl' %window out the other 2 axes and then do a 1d correlation in the remaining
%               '1d_cart_bb' %differences are cacluated as the vector sum instead of difference
%               'radial_cl'  %Euclidean distance/L2 norm of the differences
%               'radial_bb'  %differences are calculates as the vector sum
%               '3d_cart_cl'
%               '3d_cart_bb'

%output
%norm g2 amplitude
%norm g2 vector with cordinates

%improvements
%- fit g2 amplitude
%- implement all the options
%- think more carefully about how the coicidences should be normalized
%  - results should be invariant under
%  - QE change
%  - Change in how the shot is cut up
% - implement bootstraping
%   - most usefull at the shot level (more practical to implement & more informative)
% - different smoothing for G2(in-shot) G2(between-shot)

%find the num of counts in each shot
num_counts=cellfun(@(x)size(x,1),counts);

if ~isfield(corr_opts,'verbose') || ~islogical(corr_opts.verbose)
    corr_opts.verbose = true;
end

if isfield(corr_opts,'timer') && corr_opts.timer
    tic
end

if isfield(corr_opts,'calc_err') && corr_opts.calc_err
    corr_opts.fit=true;
elseif ~isfield(corr_opts,'fit')
    corr_opts.fit=false;
end

if ~isfield(corr_opts,'calc_err')
    corr_opts.calc_err = false;
end

if ~isfield(corr_opts,'direction_labels') %the label for each direction
    direction_label={'x','y','z'};
else
    direction_label=corr_opts.direction_labels;
end

% set the number of updates in a smart way
% should input check everything used here
if ~isfield(corr_opts,'progress_updates') || isnan(corr_opts.progress_updates)
    update_time=2;
    pairs_per_sec=5e7*(1+corr_opts.do_pre_mask*10);
    dyn_updates=round(update_time*size(counts,2)*...
        mean((corr_opts.attenuate_counts*num_counts(1,:)).^2)/(pairs_per_sec));
    corr_opts.progress_updates=min(100,max(5,dyn_updates));
end
if ~isfield(corr_opts,'sampling_method')
    corr_opts.sampling_method = 'basic';
end
if isfield(corr_opts,'norm_samp_factor')
    if corr_opts.norm_samp_factor<0.01 || corr_opts.norm_samp_factor>2e7
        error('corr_opts.norm_samp_factor exceeds limits');
    end
else
    corr_opts.norm_samp_factor=3;
end

if isfield(corr_opts,'sample_proportion')
    if corr_opts.sample_proportion<=0 || corr_opts.sample_proportion>1
        error('Sample proportion must be between 0 and 1');
    end
else
    corr_opts.sample_proportion=0.15;
end

if ~isfield(corr_opts,'sort_norm')
    corr_opts.sort_norm=false;
end

if ~isfield(corr_opts,'fig')
    corr_opts.fig='corr. output';
end

if ~isfield(corr_opts,'param_num')
    corr_opt.param_num = 2;
end

if corr_opt.param_num == 4 %full freedom gaussian fit
    fun1d =  @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+b(4);
    inital_guess=[2.5,0,0.01,1];
elseif corr_opt.param_num == 3 %gaussian fit with fixed offset
    fun1d =  @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+1;
    inital_guess=[2.5,0,0.01];
elseif corr_opt.param_num == 2 %centered gaussian fit with fixed offset
    fun1d =  @(b,x) b(1).*exp(-((x).^2)./(2*b(2).^2))+1;
    inital_guess=[2.5,0.01];
else
    warning('Invalid number of fit parameters, using 2 instead')
    fun1d =  @(b,x) b(1).*exp(-((x).^2)./(2*b(2).^2))+1;
    inital_guess=[2.5,0.01];
end
fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);
if isequal(corr_opts.type(end-2:end),'_cl')
    corr_opts.cl_or_bb=false;
elseif isequal(corr_opts.type(end-2:end),'_bb')
    corr_opts.cl_or_bb=true;
end
if isequal(corr_opts.type,'1d_cart_cl')  || isequal(corr_opts.type,'1d_cart_bb')
    corr_func=@corr_1d_cart;
    corr_density='one_d_corr_density';
    centers='x_centers';
    corr_unc='corr_unc';
    bin_centers=(corr_opts.redges(2:end)+corr_opts.redges(1:end-1))/2;
elseif isequal(corr_opts.type,'radial_cl')  || isequal(corr_opts.type,'radial_bb')
    corr_func=@corr_radial;
    corr_density='rad_corr_density';
    centers='rad_centers';
    corr_unc='rad_corr_unc';
    bin_centers=sqrt((corr_opts.redges(2:end).^3-corr_opts.redges(1:end-1).^3)./(3.*(corr_opts.redges(2:end)-corr_opts.redges(1:end-1))));
elseif isequal(corr_opts.type,'3d_cart_cl')  || isequal(corr_opts.type,'3d_cart_bb')
    warning('3d corrs temporarily removed')
end

% if isequal(corr_opts.type,'1d_cart_cl')  || isequal(corr_opts.type,'1d_cart_bb')
%     direction_label=direction_label{corr_opts.one_d_dimension};
%
%     if isequal(corr_opts.type,'1d_cart_cl')
%         corr_opts.cl_or_bb=false;
%     elseif isequal(corr_opts.type,'1d_cart_bb')
%         corr_opts.cl_or_bb=true;
%     end
%     if corr_opts.verbose
%         fprintf('calculating intra-shot correlations \n')
%     end
%     shotscorr=corr_1d_cart(corr_opts,counts);
%     if corr_opts.calc_err
%         shotscorr_err=corr_err(@corr_1d_cart,corr_opts,counts);
%     end
%     if corr_opts.plots
%         stfig(corr_opts.fig,'add_stack',1);
%         clf
%         set(gcf,'color','w');
%         subplot(1,3,1)
%         if corr_opts.calc_err
%             errorbar(shotscorr.x_centers,shotscorr.one_d_corr_density,shotscorr_err.results.se_fun_whole_unweighted,'kx-')
%         else
%             plot(shotscorr.x_centers,shotscorr.one_d_corr_density,'.k-','MarkerSize',10)
%         end
%         title('In Shot X Dist (windowed)')
%         ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
%         xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
%         pause(1e-6);
%     end
%     norm_sort_dir=corr_opts.sorted_dir;
%     if ~corr_opts.sort_norm,norm_sort_dir=nan; end
%     if size(counts,1)>1
%         chunk_num = max([cellfun(@(x)size(x,1),counts(1,:)),cellfun(@(x)size(x,1),counts(2,:))]);
%         counts_chunked(1,:)=chunk_data(counts(1,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
%         counts_chunked(2,:)=chunk_data(counts(2,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
%         %         counts_chunked=counts_chunked(:,1:end-1);
%     else
%         counts_chunked=chunk_data(counts,corr_opts.norm_samp_factor,norm_sort_dir);
%     end
%     corr_opts.do_pre_mask=corr_opts.sort_norm; %can only do premask if data is sorted
%     if corr_opts.verbose
%         fprintf('calculating inter-shot correlations \n')
%     end
%     normcorr=corr_1d_cart(corr_opts,counts_chunked);
%     if corr_opts.calc_err
%         normcorr_err=corr_err(@corr_1d_cart,corr_opts,counts_chunked);
%     end
%     if corr_opts.plots
%         subplot(1,3,2)
%         if corr_opts.calc_err
%             errorbar(normcorr.x_centers,normcorr.one_d_corr_density,normcorr_err.results.se_fun_whole_unweighted,'kx-')
%         else
%             plot(normcorr.x_centers,normcorr.one_d_corr_density,'.k-','MarkerSize',10)
%         end
%         title('Between Shot X Dist (windowed)')
%         ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
%         xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
%         subplot(1,3,3)
%         xg2=shotscorr.one_d_corr_density./normcorr.one_d_corr_density;
%         if corr_opts.calc_err
%             xg2_err=xg2.*sqrt((shotscorr_err.results.se_fun_whole_unweighted'./shotscorr.one_d_corr_density).^2+(normcorr_err.results.se_fun_whole_unweighted'./normcorr.one_d_corr_density).^2);
%             errorbar(shotscorr.x_centers,xg2,xg2_err,'kx-')
%         else
%             plot(shotscorr.x_centers,xg2,'.k-','MarkerSize',10)
%         end
%         title('Norm. Corr.')
%         ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
%         ylabel(sprintf('$g^{(2)}(\\Delta %s)$',direction_label))
%         xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
%         pause(1e-6);
%     end
%     out.in_shot_corr.x_centers=shotscorr.x_centers;
%     out.in_shot_corr.one_d_corr_density=shotscorr.one_d_corr_density;
%     out.between_shot_corr.x_centers=normcorr.x_centers;
%     out.between_shot_corr.one_d_corr_density=normcorr.one_d_corr_density;
%     out.norm_g2.x_centers=shotscorr.x_centers;
%     out.norm_g2.g2_amp=xg2;
%     if corr_opts.fit
%         if ~corr_opts.calc_err
%             fit=fitnlm(shotscorr.x_centers,xg2,...
%                 fun1d,...
%                 inital_guess,...
%                 'Options',fo);
%         else
%             w = 1./xg2_err.^2;
%             fit=fitnlm(shotscorr.x_centers,xg2,...
%                 fun1d,...
%                 inital_guess,...
%                 'Weight',w,...
%                 'Options',fo);
%         end
%         out.fit = fit;
%         b = fit.Coefficients.Estimate; %fitted parameters
%         b_unc = fit.Coefficients.SE; %uncertainty in parameters
%         out.norm_g2.fitted_g2peak = b(1)+1;
%         out.norm_g2.fitted_g2peak_unc = b_unc(1);
%         if corr_opts.plots
%             hold on
%             xx = linspace(min(shotscorr.x_centers),max(shotscorr.x_centers),3e3)';
%             [ypred,ypredci] = predict(fit,xx,'Simultaneous',true);
%             plot(xx,ypred,'b-', xx,ypredci,'r-');
%             text(0.003,2.0,sprintf('fitted g2(0) amplitude:%s\n',string_value_with_unc(b(1)+1,b_unc(1),'type','b','separator',0)));
%             text(0.003,1.8,sprintf('fitted g2 width:%s\n',string_value_with_unc(abs(b(2)),b_unc(2),'type','b','separator',0)));
%         end
%     end
%     [g2peak,indx] = max(xg2);
%     if corr_opts.calc_err
%         out.in_shot_corr.corr_unc = shotscorr_err;
%         out.between_shot_corr.corr_unc = normcorr_err;
%         out.norm_g2.g2_unc = xg2_err;
%         g2peak_unc = xg2_err(indx);
%         out.norm_g2.g2peak_unc = g2peak_unc;
%         if corr_opts.verbose
%             fprintf('g2 peak amplitude         %s\n',string_value_with_unc(g2peak,g2peak_unc,'type','b','separator',0))
%         end
%     else
%         if corr_opts.verbose
%             fprintf('g2 peak amplitude         %4.2f \n',g2peak)
%         end
%     end
%     if corr_opts.fit && corr_opts.verbose
%         fprintf('fitted g2(0) amplitude         %s\n',string_value_with_unc(b(1)+1,b_unc(1),'type','b','separator',0))
%         fprintf('fitted g2 width         %s\n',string_value_with_unc(abs(b(2)),b_unc(2),'type','b','separator',0))
%     end
%     out.norm_g2.g2peak=g2peak;



% elseif isequal(corr_opts.type,'radial_cl')  || isequal(corr_opts.type,'radial_bb')
%     if isequal(corr_opts.type,'radial_cl')
%         corr_opts.cl_or_bb=false;
%     elseif isequal(corr_opts.type,'radial_bb')
%         corr_opts.cl_or_bb=true;
%     end


if ~corr_opts.calc_err
    if corr_opts.verbose
        cli_format_text('Calculating Correlations','c',3)
        fprintf('calculating intra-shot correlations \n')
    end
    shotscorr=corr_func(corr_opts,counts);
    %%
    norm_sort_dir=corr_opts.sorted_dir;
    if ~corr_opts.sort_norm,norm_sort_dir=nan; end
    if strcmp(corr_opts.sampling_method,'basic')
        if size(counts,1)>1
            %set the number of chunks to be at least as many as the heighest
            %count number
            chunk_num = max([cellfun(@(x)size(x,1),counts(1,:)),cellfun(@(x)size(x,1),counts(2,:))]);
            counts_chunked(1,:)=chunk_data(counts(1,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
            counts_chunked(2,:)=chunk_data(counts(2,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
        else
            counts_chunked=chunk_data(counts,corr_opts.norm_samp_factor,norm_sort_dir);
        end
    elseif strcmp(corr_opts.sampling_method,'complete')
        if size(counts,1)>1
            chunk_num = min([sum(cellfun(@(x)size(x,1),counts(1,:)))-size(counts{1,1},1),sum(cellfun(@(x)size(x,1),counts(2,:)))-size(counts{2,1},1)]);
            counts_chunked(1,:)=chunk_data_complete(counts(1,:),corr_opts.sample_proportion,norm_sort_dir,chunk_num);
            counts_chunked(2,:)=chunk_data_complete(counts(2,:),corr_opts.sample_proportion,norm_sort_dir,chunk_num);
        else
            counts_chunked=chunk_data_complete(counts,corr_opts.sample_proportion,norm_sort_dir);
        end
    end
    corr_opts.do_pre_mask=corr_opts.sort_norm; %can only do premask if data is sorted
    if corr_opts.verbose
        fprintf('calculating inter-shot correlations \n')
    end
    normcorr=corr_func(corr_opts,counts_chunked);
    %%
    xg2=shotscorr.(corr_density)./normcorr.(corr_density);
    [g2peak,~] = max(xg2);
else
    corrs=corr_err(corr_func,corr_opts,counts);
    shotscorr.(centers) = bin_centers;%
    shotscorr.(corr_density) = corrs.rawcorr.val; %weighted average
    shotscorr.(corr_unc) = corrs.rawcorr.unc;
    normcorr.(centers) = shotscorr.(centers);
    normcorr.(corr_density) = corrs.rawcorr.val; %weighted average
    normcorr.(corr_unc) = corrs.normcorr.unc;
    xg2=corrs.corr_density.val;
    xg2_err=corrs.corr_density.unc;
    
    %         out.in_shot_corr.corr_unc = shotscorr_err;
    %         out.between_shot_corr.corr_unc = normcorr_err;
    [g2peak,indx] = max(xg2);
    out.norm_g2.g2_unc = xg2_err;
    g2peak_unc = xg2_err(indx);
    out.norm_g2.g2peak_unc = g2peak_unc;
end

out.in_shot_corr.(centers)=shotscorr.(centers);
out.in_shot_corr.(corr_density)=shotscorr.(corr_density);
out.between_shot_corr.(centers)=normcorr.(centers);
out.between_shot_corr.(corr_density)=normcorr.(corr_density);
out.norm_g2.(centers)=shotscorr.(centers);
out.norm_g2.g2_amp=xg2;
is_data_flat = isdataflat(xg2,0.2);%0.2
if corr_opts.fit
    if ~corr_opts.calc_err && ~is_data_flat
        fit=fitnlm(shotscorr.(centers),xg2,...
            fun1d,...
            inital_guess,...
            'Options',fo);
        out.fit = fit;
        b = fit.Coefficients.Estimate; %fitted parameters
        b_unc = fit.Coefficients.SE; %uncertainty in parameters
    elseif ~corr_opts.calc_err
        fun1d = @(b,x) (b(1)+1).*ones(size(x));
        inital_guess = 0;
        fit=fitnlm(shotscorr.(centers),xg2,...
            fun1d,...
            inital_guess,...
            'Options',fo);
        out.fit = fit;
        b = [fit.Coefficients.Estimate, 0]; %fitted parameters
        b_unc = [fit.Coefficients.SE, 0]; %uncertainty in parameters
    else
        out.fit = corrs.fit;
        b = corrs.fit.val; %fitted parameters
        b(2) = abs(b(2));
        b_unc = corrs.fit.unc; %uncertainty in parameters
        if is_data_flat
            fun1d = @(b,x) (b(1)+1).*ones(size(x));
        else
            fun1d =  @(b,x) b(1).*exp(-((x).^2)./(2*b(2).^2))+1;
        end
    end
    out.norm_g2.fitted_g2peak = b(1)+1;
    out.norm_g2.fitted_g2peak_unc = b_unc(1);
end

if corr_opts.plots
    stfig(corr_opts.fig);
    clf
    set(gcf,'color','w');
    subplot(1,3,1)
    if corr_opts.calc_err
        errorbar(shotscorr.(centers),shotscorr.rad_corr_density,shotscorr.(corr_unc),'kx-')
    else
        plot(shotscorr.(centers),shotscorr.rad_corr_density,'.k-','MarkerSize',10)
    end
    title('In Shot X Dist (windowed)')
    ylabel('$G^{(2)}(\Delta r)$ coincedence density')
    xlabel('$\delta r$ Seperation')
    subplot(1,3,2)
    if corr_opts.calc_err
        errorbar(normcorr.(centers),normcorr.rad_corr_density,normcorr.(corr_unc),'kx-')
    else
        plot(normcorr.(centers),normcorr.rad_corr_density,'.k-','MarkerSize',10)
    end
    title('Between Shot X Dist (windowed)')
    ylabel('$G^{(2)}(\Delta r)$ coincedence density')
    xlabel('$\delta r$ Seperation')
    subplot(1,3,3)
    if corr_opts.calc_err
        errorbar(shotscorr.(centers),xg2,xg2_err,'kx-')
    else
        plot(shotscorr.(centers),xg2,'.k-','MarkerSize',10)
    end
    title('Norm. Corr.')
    ylabel('$g^{(2)} (\Delta r)$')
    xlabel('$\delta r$ Seperation')
    if corr_opts.plots && corr_opts.fit
        hold on
        xx = linspace(0,max(shotscorr.(centers)),3e3)';
        if ~corr_opts.calc_err
            [ypred,ypredci] = predict(fit,xx,'Simultaneous',true);
            plot(xx,ypred,'b-', xx,ypredci,'r-');
        else
            plot(xx,fun1d(b,xx),'b-')
            plot(xx,fun1d(b-b_unc,xx),'r-')
            plot(xx,fun1d(b+b_unc,xx),'r-')
        end
    end
    if g2peak < 1.5
        ylim([0.8 2.1])
    end
    pause(1e-6);
end

if corr_opts.verbose
    if corr_opts.calc_err
        fprintf('g2 peak amplitude         %s\n',string_value_with_unc(g2peak,g2peak_unc,'type','b','separator',0))
    else
        fprintf('g2 peak amplitude         %4.2f \n',g2peak)
    end
    if corr_opts.fit
        fprintf('fitted g2(0) amplitude         %s\n',string_value_with_unc(b(1)+1,b_unc(1),'type','b','separator',0))
        fprintf('fitted g2 width         %s\n',string_value_with_unc(abs(b(2)),b_unc(2),'type','b','separator',0))
    end
end
out.norm_g2.g2peak=g2peak;



%     out = cell(1,3);
%     for dimension = 1:3
%         if corr_opts.calc_err
%             lin_style = {'kx-','bx-','rx-'};
%         else
%             lin_style = {'.k-','.b-','.r-'};
%         end
%         corr_opts.one_d_dimension = dimension;
%         if isequal(corr_opts.type,'3d_cart_cl')
%             corr_opts.cl_or_bb=false;
%         elseif isequal(corr_opts.type,'3d_cart_bb')
%             corr_opts.cl_or_bb=true;
%         end
%         shotscorr=corr_1d_cart(corr_opts,counts);
%         if corr_opts.calc_err
%             shotscorr_err=corr_err(@corr_1d_cart,corr_opts,counts);
%         end
%         %shotscorr_high=corr_1d_cart_high_mem(corr_opts,counts);
%         if corr_opts.plots
%             if dimension == 1
%                 stfig(corr_opts.fig,'add_stack',1);
%                 clf
%                 set(gcf,'color','w');
%             end
%             subplot(1,3,1)
%             hold on
%             if corr_opts.calc_err
%                 errorbar(shotscorr.x_centers,shotscorr.one_d_corr_density,shotscorr_err.results.se_fun_whole_unweighted,lin_style{dimension})
%             else
%                 plot(shotscorr.x_centers,shotscorr.one_d_corr_density,lin_style{dimension},'MarkerSize',10)
%             end
%             title('In Shot X Dist (windowed)')
%             ylabel('$G^{(2)}(\Delta)$ coincedence density')
%             xlabel('$\Delta$ Seperation')
%             pause(1e-6);
%         end
%         norm_sort_dir=corr_opts.sorted_dir;
%         if ~corr_opts.sort_norm,norm_sort_dir=nan; end
%         if size(counts,1)>1
%             chunk_num = max([cellfun(@(x)size(x,1),counts(1,:)),cellfun(@(x)size(x,1),counts(2,:))]);
%             counts_chunked(1,:)=chunk_data(counts(1,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
%             counts_chunked(2,:)=chunk_data(counts(2,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
%             %         counts_chunked=counts_chunked(:,1:end-1);
%         else
%             counts_chunked=chunk_data(counts,corr_opts.norm_samp_factor,norm_sort_dir);
%         end
%         corr_opts.do_pre_mask=corr_opts.sort_norm; %can only do premask if data is sorted
%         fprintf('calculating inter-shot correlations \n')
%         normcorr=corr_1d_cart(corr_opts,counts_chunked);
%         if corr_opts.calc_err
%             normcorr_err=corr_err(@corr_1d_cart,corr_opts,counts_chunked);
%         end
%         if corr_opts.plots
%             subplot(1,3,2)
%             hold on
%             if corr_opts.calc_err
%                 errorbar(normcorr.x_centers,normcorr.one_d_corr_density,normcorr_err.results.se_fun_whole_unweighted,lin_style{dimension})
%             else
%                 plot(normcorr.x_centers,normcorr.one_d_corr_density,lin_style{dimension},'MarkerSize',10)
%             end
%             title('Between Shot X Dist (windowed)')
%             ylabel('$G^{(2)}(\Delta)$ coincedence density')
%             xlabel('$\Delta$ Seperation')
%             subplot(1,3,3)
%             hold on
%             xg2=shotscorr.one_d_corr_density./normcorr.one_d_corr_density;
%             if corr_opts.calc_err
%                 xg2_err=xg2.*sqrt((shotscorr_err.results.se_fun_whole_unweighted'./shotscorr.one_d_corr_density).^2+(normcorr_err.results.se_fun_whole_unweighted'./normcorr.one_d_corr_density).^2);
%                 errorbar(shotscorr.x_centers,xg2,xg2_err,lin_style{dimension})
%             else
%                 plot(shotscorr.x_centers,xg2,lin_style{dimension},'MarkerSize',10)
%             end
%             title('Norm. Corr.')
%             ylabel('$g^{(2)}(\Delta)$')
%             xlabel('$\Delta$ Seperation')
%             pause(1e-6);
%         end
%         out{dimension}.in_shot_corr.x_centers=shotscorr.x_centers;
%         out{dimension}.in_shot_corr.one_d_corr_density=shotscorr.one_d_corr_density;
%         out{dimension}.between_shot_corr.x_centers=normcorr.x_centers;
%         out{dimension}.between_shot_corr.one_d_corr_density=normcorr.one_d_corr_density;
%         out{dimension}.norm_g2.x_centers=shotscorr.x_centers;
%         out{dimension}.norm_g2.g2_amp=xg2;
%         [g2peak,indx] = max(xg2);
%         if corr_opts.fit
%             if ~corr_opts.calc_err
%                 fit=fitnlm(shotscorr.x_centers,xg2,...
%                     fun1d,...
%                     inital_guess,...
%                     'Options',fo);
%             else
%                 w = 1./xg2_err.^2;
%                 fit=fitnlm(shotscorr.x_centers,xg2,...
%                     fun1d,...
%                     inital_guess,...
%                     'Weight',w,...
%                     'Options',fo);
%             end
%             out{dimension}.fit = fit;
%             if corr_opts.plots
%                 hold on
%                 xx = linspace(0,max(shotscorr.x_centers),3e3)';
%                 [ypred,ypredci] = predict(fit,xx,'Simultaneous',true);
%                 plot(xx,ypred,'b-', xx,ypredci,'r-');
%             end
%         end
%         if corr_opts.calc_err
%             out{dimension}.in_shot_corr.corr_unc = shotscorr_err;
%             out{dimension}.between_shot_corr.corr_unc = normcorr_err;
%             out{dimension}.norm_g2.g2_unc = xg2_err;
%             g2peak_unc = xg2_err(indx);
%             out{dimension}.norm_g2.g2peak_unc = g2peak_unc;
%             fprintf('g2 peak amplitude         %s\n',string_value_with_unc(g2peak,g2peak_unc,'type','b','separator',0))
%         else
%             fprintf('g2 peak amplitude         %4.2f \n',g2peak)
%         end
%         out{dimension}.norm_g2.g2peak=g2peak;
%     end
%     if corr_opts.plots
%         subplot(1,3,2)
%         legend(sprintf('$\\Delta %s$ Seperation',direction_label{1}),...
%             sprintf('$\\Delta %s$ Seperation',direction_label{2}),...
%             sprintf('$\\Delta %s$ Seperation',direction_label{3}))
%     end
% end

if isfield(corr_opts,'timer') && corr_opts.timer
    toc
end

end

