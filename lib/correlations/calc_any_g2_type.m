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

%% OPTIONS CHECKS

if ~isfield(corr_opts,'full_error')
    corr_opts.full_error = false;
end

if ~isfield(corr_opts,'verbose') || ~islogical(corr_opts.verbose)
    corr_opts.verbose = true;
end

if isfield(corr_opts,'timer') && corr_opts.timer
    tic
end

% if isfield(corr_opts,'calc_err') && corr_opts.calc_err
%     corr_opts.fit=true;
if ~isfield(corr_opts,'fit')
    corr_opts.fit=false;
end

if ~isfield(corr_opts,'calc_err')
    corr_opts.calc_err = false;
end

if ~isfield(corr_opts,'direction_labels') %the label for each direction
    direction_labels={'x','y','z'};
else
    direction_labels=corr_opts.direction_labels;
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
if ~isfield(corr_opts,'g2_norm') %whether or not to calculate g2 normalisation factor
    corr_opts.g2_norm = true;
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

if ~isfield(corr_opts,'one_d_smoothing')
    corr_opts.one_d_smoothing=nan;
end

if ~isfield(corr_opts,'param_num')
    corr_opts.param_num = 2;
end

if corr_opts.param_num == 4 %full freedom gaussian fit
    fun1d =  @(b,x) b(1).*exp(-((x-b(3)).^2)./(2*b(2).^2))+b(4);
elseif corr_opts.param_num == 3 %gaussian fit with fixed offset
    fun1d =  @(b,x) b(1).*exp(-((x-b(3)).^2)./(2*b(2).^2))+1;
elseif corr_opts.param_num == 2 %centered gaussian fit with fixed offset
    fun1d =  @(b,x) b(1).*exp(-((x).^2)./(2*b(2).^2))+1;
else
    warning('Invalid number of fit parameters, using 2 instead')
    fun1d =  @(b,x) b(1).*exp(-((x).^2)./(2*b(2).^2))+1;
end

fun2d =  @(b,x) b(1).*exp(-(x(:,:,1).^2./(2*b(2).^2)+x(:,:,2).^2./(2*b(3).^2)))+1;

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
    direction_label=direction_labels{corr_opts.one_d_dimension};
    bin_centers=(corr_opts.one_d_edges(2:end)+corr_opts.one_d_edges(1:end-1))/2;
elseif isequal(corr_opts.type,'radial_cl')  || isequal(corr_opts.type,'radial_bb')
    corr_func=@corr_radial;
    corr_density='rad_corr_density';
    centers='rad_centers';
    corr_unc='rad_corr_unc';
    direction_label='r';
    bin_centers=sqrt((corr_opts.redges(2:end).^3-corr_opts.redges(1:end-1).^3)./(3.*(corr_opts.redges(2:end)-corr_opts.redges(1:end-1))));
elseif isequal(corr_opts.type,'3d_cart_cl')  || isequal(corr_opts.type,'3d_cart_bb')
    warning('3d corrs temporarily removed')
elseif isequal(corr_opts.type,'angle')
    corr_func=@corr_angle;
    corr_density='ang_corr_density';
    centers='centers';
    corr_unc='corr_unc';
    direction_label='\cos(\theta)';
    bin_centers=(corr_opts.edges(2:end)+corr_opts.edges(1:end-1))/2;
elseif isequal(corr_opts.type,'2d_cart_cl')  || isequal(corr_opts.type,'2d_cart_bb')
    corr_func=@corr_2d_cart;
    corr_density='two_d_corr_density';
    centers='centers';
    corr_unc='corr_unc';
    direction_label=direction_labels(corr_opts.two_d_dimensions);
    bin_centers=(corr_opts.edges(2:end)+corr_opts.edges(1:end-1))/2;
elseif isequal(corr_opts.type,'1d_vol_cl')  || isequal(corr_opts.type,'1d_vol_bb')
    corr_func=@corr_vol_cart;
    corr_density='one_d_corr_density';
    centers='x_centers';
    corr_unc='corr_unc';
    direction_label='';%direction_labels(corr_opts.one_d_dimension);
    bin_centers=(corr_opts.one_d_edges(2:end)+corr_opts.one_d_edges(1:end-1))/2;
end

if ~isfield(corr_opts,'print_update')
    corr_opts.print_update = true;
end

if ~isfield(corr_opts,'plots')
    corr_opts.plots = true;
end

%% CALCULATE THE CORRELATIONS

%normal method without error calculations
if corr_opts.verbose
    cli_format_text('Calculating Correlations','c',3)
    fprintf('calculating intra-shot correlations \n')
end
shotscorr=corr_func(corr_opts,counts);
%%
if corr_opts.g2_norm
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
        corr_opts.normalisation_factor = (size(counts,2)/nanmean(cellfun(@(x)size(x,1),counts),'all'))^2;
    elseif strcmp(corr_opts.sampling_method,'complete')
        counts_chunked=chunk_data_complete(counts,corr_opts.sample_proportion,norm_sort_dir);
    end
    corr_opts.do_pre_mask=corr_opts.sort_norm; %can only do premask if data is sorted
    if corr_opts.verbose
        fprintf('calculating inter-shot correlations \n')
    end
    normcorr=corr_func(corr_opts,counts_chunked);
    %%
    xg2=shotscorr.(corr_density)./normcorr.(corr_density);
    [g2peak,~] = max(xg2);
end
if corr_opts.calc_err %bootstrap your errors
    if corr_opts.full_error && corr_opts.g2_norm
        corrs=corr_err(corr_func,corr_opts,counts);
        normcorr.(corr_unc) = corrs.normcorr.unc;
    else
        corrs_out = corr_err_raw(corr_func,corr_opts,counts)';
        shotscorr.(corr_unc) = corrs_out.range;
        shotscorr.corr_ci = corrs_out.ci;
        out.in_shot_corr.(corr_unc) = shotscorr.(corr_unc); 
        out.in_shot_corr.corr_ci = shotscorr.corr_ci; 
    end
    %     shotscorr.(centers) = bin_centers;%
    %     shotscorr.(corr_density) = corrs.rawcorr.val; %weighted average
    %     shotscorr.(corr_unc) = corrs.rawcorr.unc;
    %     normcorr.(centers) = shotscorr.(centers);
    %     normcorr.(corr_density) = corrs.rawcorr.val; %weighted average
    %     xg2=corrs.corr_density.val;
    if corr_opts.g2_norm
        xg2_err = shotscorr.(corr_unc).'./normcorr.(corr_density);
    [g2peak,indx] = max(xg2);
    out.norm_g2.g2_unc = xg2_err;
    g2peak_unc = xg2_err(indx);
    out.norm_g2.g2peak_unc = g2peak_unc;
    end
end

out.in_shot_corr.(centers)=shotscorr.(centers);
out.in_shot_corr.(corr_density)=shotscorr.(corr_density);
if isequal(corr_opts.type,'1d_vol_cl')  || isequal(corr_opts.type,'1d_vol_bb')
out.in_shot_corr.shots = shotscorr.relevant_shots;
end
if corr_opts.g2_norm
    out.between_shot_corr.(centers)=normcorr.(centers);
    out.between_shot_corr.(corr_density)=normcorr.(corr_density);
    out.norm_g2.(centers)=shotscorr.(centers);
    out.norm_g2.g2_amp=xg2;
else
    xg2 = shotscorr.(corr_density);
end

%% FIT CORRELATIONS

if corr_opts.fit
    % Check if data is flat
    %         is_data_flat = isdataflat(xg2,1.0);%0.2
    if isequal(corr_opts.type,'2d_cart_cl')  || isequal(corr_opts.type,'2d_cart_bb')
        %         [muHat,sigmaHat] = normfit(shotscorr.(centers),0.01,zeros(size(shotscorr.(centers))),abs(xg2-1).^2);
        inital_guess=[max(max(xg2))-1,0.001,0.001];
        xy=shotscorr.(centers);
        [X(:,:,1),X(:,:,2)]=meshgrid(xy{1},xy{2});
        b = lsqcurvefit(fun2d, inital_guess, X, xg2.', [0  -10  -10],  [100  10  10]);
        b_unc = b.*0.1;
        out.fit = b;
        %         fit=fitnlm(X,xg2,...
        %             fun2d,...
        %             inital_guess,...
        %             'Options',fo);
        %     elseif ~corr_opts.calc_err
    else%if ~corr_opts.calc_err || ~corr_opts.full_error %&& ~is_data_flat
        [muHat,sigmaHat] = normfit(shotscorr.(centers),0.01,zeros(size(shotscorr.(centers))),abs(xg2-1).^2);
        inital_guess=[max(xg2)-1,max(sigmaHat,0.001)];
        if corr_opts.param_num == 4 %full freedom gaussian fit
            inital_guess = [inital_guess,muHat,1];
        elseif corr_opts.param_num == 3 %gaussian fit with fixed offset
            inital_guess = [inital_guess,muHat];
        end
        fit=fitnlm(shotscorr.(centers),xg2,...
            fun1d,...
            inital_guess,...
            'Options',fo);
        %     elseif ~corr_opts.calc_err
        fun1d_flat = @(b,x) (b(1)+1).*ones(size(x));
        inital_guess = nanmean(xg2)-1;
        fit_flat=fitnlm(shotscorr.(centers),xg2,...
            fun1d_flat,...
            inital_guess,...
            'Options',fo);
        if fit.Rsquared.Adjusted>0.2 || (isfield(corr_opts,'gaussian_fit') && corr_opts.gaussian_fit)
            out.fit = fit;
            b = fit.Coefficients.Estimate; %fitted parameters
            b_unc = fit.Coefficients.SE; %uncertainty in parameters
            is_data_flat = false;
        else
            out.fit = fit_flat;
            fit = fit_flat;
            b = [fit_flat.Coefficients.Estimate, 0]; %fitted parameters
            b_unc = [fit_flat.Coefficients.SE, 0]; %uncertainty in parameters
            is_data_flat = true;
        end
        %     else
        %         out.fit = corrs.fit;
        %         b = corrs.fit.val; %fitted parameters
        %         b(2) = abs(b(2));
        %         b_unc = corrs.fit.unc; %uncertainty in parameters
        %         is_data_flat = isdataflat(xg2,1.0);%0.2
        %         if is_data_flat
        %             fun1d = @(b,x) (b(1)+1).*ones(size(x));
        %         else
        %             fun1d =  @(b,x) b(1).*exp(-((x).^2)./(2*b(2).^2))+1;
        %         end
    end
    if corr_opts.g2_norm
        out.norm_g2.fitted_g2peak = b(1)+1;
        out.norm_g2.fitted_g2peak_unc = b_unc(1);
    end
end

%% PLOT CORRELATIONS

if corr_opts.plots
    stfig(corr_opts.fig);
    clf
    set(gcf,'color','w');
    if isequal(corr_opts.type,'2d_cart_cl')  || isequal(corr_opts.type,'2d_cart_bb')
        xy=shotscorr.(centers);
        [X(:,:,1),X(:,:,2)]=meshgrid(xy{1},xy{2});
        surf(X(:,:,1),X(:,:,2),xg2.')
        zlabel(sprintf('$g^{(2)}(\\Delta %s, \\Delta %s)$ coincedence density',direction_label{1},direction_label{2}))
        ylabel(sprintf('$\\Delta %s$ Seperation',direction_label{2}))
        xlabel(sprintf('$\\Delta %s$ Seperation',direction_label{1}))
    else
        if corr_opts.g2_norm
            subplot(1,3,1)
        end
        %     hold on
        if corr_opts.calc_err
            errorbar(shotscorr.(centers),shotscorr.(corr_density),shotscorr.(corr_unc),'kx-')
        else
            plot(shotscorr.(centers),shotscorr.(corr_density),'.k-','MarkerSize',10)
        end
        title('In Shot X Dist (windowed)')
        ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
        xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
        if corr_opts.g2_norm
            subplot(1,3,2)
            %     hold on
            if corr_opts.calc_err && corr_opts.full_error
                errorbar(normcorr.(centers),normcorr.(corr_density),normcorr.(corr_unc),'kx-')
            else
                plot(normcorr.(centers),normcorr.(corr_density),'.k-','MarkerSize',10)
            end
            title('Between Shot X Dist (windowed)')
            ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
            xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
            subplot(1,3,3)
            %     hold on
            if corr_opts.calc_err
                errorbar(shotscorr.(centers),xg2,xg2_err,'kx-')
            else
                plot(shotscorr.(centers),xg2,'.k-','MarkerSize',10)
            end
            title('Norm. Corr.')
            ylabel(sprintf('$g^{(2)}(\\Delta %s)$',direction_label))
            xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
            if corr_opts.plots && corr_opts.fit
                hold on
                xx = linspace(min([min(shotscorr.(centers)),0]),max(shotscorr.(centers)),3e3)';
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
        end
    end
    pause(1e-6);
end

%% WRITE OUT RESULTS

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

if corr_opts.g2_norm
    out.norm_g2.g2peak=g2peak;
end

if isfield(corr_opts,'timer') && corr_opts.timer
    toc
end

end

