function [bec_centres,bec_widths,bec_counts,centre_OK] = find_dist_cen(data, opts_cent)
% Mandatory data inputs:
%     data.counts_txy
%     data.shot_num
% Mandatory Input options:
%     opts.cent.auto
%     opts.cent.crop
%     opts.cent.sigma
%     opts.cent.threshold
%     opts_cent.method
%     opts.partition.t_win
%     opts.visual
if isfield(opts_cent,'hotspot_mask') && opts_cent.hotspot_mask
    data=hotspot_mask(data); %mask out hot spots
    counts_txy=data.masked.counts_txy;
else
    counts_txy=data.counts_txy;
end

num_shots = length(data.shot_num);
bec_centres = zeros(num_shots, 3);
bec_widths= zeros(num_shots, 3);
bec_counts = zeros(num_shots, 1);
centre_OK = zeros(num_shots, 1);
axis_label = {'z','x','y'};
lims = opts_cent.crop;
if length(opts_cent.threshold) < 3
    opts_cent.threshold = opts_cent.threshold(1)*[1,1,1];
end
if length(opts_cent.sigma) < 3
    opts_cent.sigma = opts_cent.sigma(1)*[1,1,1];
end
if length(opts_cent.min_threshold) < 3
    opts_cent.min_threshold = opts_cent.min_threshold(1)*[1,1,1];
end
if ~iscell(opts_cent.method)
    opts_cent.method{2} = opts_cent.method{1};
    opts_cent.method{3} = opts_cent.method{1};
elseif length(opts_cent.method) < 3
    opts_cent.method{2} = opts_cent.method{1};
    opts_cent.method{3} = opts_cent.method{1};
end

fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);

if opts_cent.visual>1
    stfig('BEC centering dist');
    clf
end

for this_idx = 1:num_shots % Loop over all shots
    this_txy = counts_txy{this_idx};
    if size(this_txy,2) ~= 3 || isempty(this_txy)
        centre_OK(this_idx) = 0;
        bec_centres(this_idx, :) = nan;
        bec_widths(this_idx, :) = nan;
        continue
    end
    trim_txy = masktxy_square(this_txy, lims);
    if size(trim_txy,2) ~= 3 || isempty(trim_txy) || size(trim_txy,1)<5
        centre_OK(this_idx) = 0;
        bec_centres(this_idx, :) = nan;
        bec_widths(this_idx, :) = nan;
        continue
    end
    bec_counts(this_idx) = size(trim_txy,1);
    if opts_cent.visual && this_idx == 1
        clf;
    end
    % fancy centering
    centre_OK(this_idx) = 1; %starts off ok until told otherwise
    for axis = 1:3
        this_axis = trim_txy(:, axis);
        this_sigma = opts_cent.sigma(axis);
        this_method = opts_cent.method{axis};
        if isempty(this_axis)
            centre_OK(this_idx) = 0;
            bec_centres(this_idx, axis) = nan;
            bec_widths(this_idx, axis) = nan;
            %             continue
        end
        count_hist = smooth_hist(this_axis,'sigma',this_sigma);
        flux = count_hist.count_rate.smooth;
        bin_centres = count_hist.bin.centers;
        outer_cut = ~or((max(bin_centres)-bin_centres)./range(bin_centres)...
            <0.02,(-min(bin_centres)+bin_centres)./range(bin_centres)<0.02);
        flux = flux(outer_cut);
        bin_centres = bin_centres(outer_cut);
        mask_upper = (flux > opts_cent.threshold(axis));
        mask_lower = (flux < opts_cent.min_threshold(axis));
        mask = or(mask_upper,mask_lower);
        if all(mask) || sum(~mask)/length(mask)<0.1
            centre_OK(this_idx) = 0;
            bec_centres(this_idx, axis) = nan;
            bec_widths(this_idx, axis) = nan;
            %             continue
        elseif sum(mask_upper)<5 && strcmp(this_method,'margin')
            centre_OK(this_idx) = 0;
            bec_centres(this_idx, axis) = nan;
            bec_widths(this_idx, axis) = nan;
            %             continue
            %             this_method = 'average';
        else
            switch this_method
                case 'margin'
                    locs = find(mask_upper);
                    margins = [min(locs(2:end-1)), max(locs(2:end-1))];
                    t_margins = bin_centres(margins);
%                     if axis == 1
%                         bec_centres(this_idx, axis) = adjusted_t_mean(t_margins(1),t_margins(2));
%                     else
                        bec_centres(this_idx, axis) = mean(t_margins);
%                     end
                    bec_widths(this_idx, axis) = diff(t_margins);
                case 'average'
                    %                     bec_centres(this_idx,axis) = nansum(flux(~mask).*bin_centres(~mask))./nansum(flux(~mask));
                    flux_masked = flux;
                    flux_masked(mask) = 0;
                    flux_masked(~mask) = flux_masked(~mask)-opts_cent.min_threshold(axis);
                    bec_centres(this_idx,axis) = trapz(bin_centres,flux_masked.*bin_centres)./trapz(bin_centres,flux_masked);
                    bec_widths(this_idx,axis) = sqrt(nansum(flux(~mask).*(bin_centres(~mask)-bec_centres(this_idx,axis)).^2)./nansum(flux(~mask)));
                case 'gauss_fit'
                    gauss =  @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2));
                    center_g = nansum(flux(~mask).*bin_centres(~mask))./nansum(flux(~mask));
                    width_g = sqrt(nansum(flux(~mask).*(bin_centres(~mask)-bec_centres(this_idx,axis)).^2)./nansum(flux(~mask)));
                    inital_guess=[mean(flux(~mask)),center_g,width_g];
                    fit=fitnlm(bin_centres(~mask),flux(~mask),...
                        gauss,...
                        inital_guess,...
                        'Options',fo);
                    fit_params = fit.Coefficients.Estimate;
                    bec_centres(this_idx,axis) = fit_params(2);
                    bec_widths(this_idx,axis) = fit_params(3);
                case 'bimod_fit'
                    center_g = nansum(flux(~mask).*bin_centres(~mask))./nansum(flux(~mask));
                    width_g = sqrt(nansum(flux(~mask).*(bin_centres(~mask)-bec_centres(this_idx,axis)).^2)./nansum(flux(~mask)));
                    inital_guess=[2*mean(flux(~mask)),center_g,width_g,mean(flux(~mask)),width_g];
                    fit_both=fitnlm(bin_centres(~mask),flux(~mask),...
                        @bimod,...
                        inital_guess,...
                        'Options',fo);
                    fit_params_both = fit_both.Coefficients.Estimate;
                    bec_centres(this_idx,axis) = fit_params_both(2);
                    bec_widths(this_idx,axis) = fit_params_both(3);
                case 'bimodsc_fit' %bimodal fit with seperate centers
                    center_g = nansum(flux(~mask).*bin_centres(~mask))./nansum(flux(~mask));
                    width_g = sqrt(nansum(flux(~mask).*(bin_centres(~mask)-bec_centres(this_idx,axis)).^2)./nansum(flux(~mask)));
                    inital_guess=[2*mean(flux(~mask)),center_g,width_g,mean(flux(~mask)),width_g,center_g];
                    fit_both=fitnlm(bin_centres(~mask),flux(~mask),...
                        @bimod_sc,...
                        inital_guess,...
                        'Options',fo);
                    fit_params_both = fit_both.Coefficients.Estimate;
                    bec_centres(this_idx,axis) = fit_params_both(6);
                    bec_widths(this_idx,axis) = fit_params_both(3);
            end
        end
        if opts_cent.visual > 1 %&& mod(this_idx,100)==0
            subplot(3, 1, axis);
            plot(bin_centres-bec_centres(this_idx,axis), flux, 'k')
            hold on
            plot(bin_centres(mask_upper)-bec_centres(this_idx,axis), flux(mask_upper), 'r')
            plot(bin_centres(mask_lower)-bec_centres(this_idx,axis), flux(mask_lower), 'b')
        end
    end
end

% cut outliers if there is enough data
% if num_shots>9
%     for ii = 1:3
%         centre_outlier = isoutlier(bec_centres(:,ii));
%         width_outlier = isoutlier(bec_widths(:,ii));
%         centre_OK = centre_OK & ~centre_outlier & ~width_outlier;
%     end
% end
if opts_cent.visual > 1
    for splt = 1:3
        subplot(3,1,splt)
        ylabel('Counts')
        xlabel(axis_label{splt})
        grid on
    end
end

if opts_cent.visual
    f=stfig('BEC centering');
    centre_OK = logical(centre_OK);
    indexes = 1:length(centre_OK);
    clf
    subplot(2,1,1)
    plot(indexes(centre_OK),bec_centres(centre_OK,:)-nanmean(bec_centres(centre_OK,:)),'.')
    title('BEC centre deviation')
    xlabel('Shot number')
    ylabel('Centre offset (mm)')
    legend('z','x','y')
    subplot(2,1,2)
    plot(indexes(centre_OK),bec_widths(centre_OK,:))
    title('BEC width')
    xlabel('Shot number')
    ylabel('widths (mm)')
    legend('z','x','y')
    if opts_cent.savefigs
        cli_header(1,'Saving images...');
        saveas(f,fullfile(opts_cent.data_out,'centering.fig'));
        saveas(f,fullfile(opts_cent.data_out,'centering.svg'));
        cli_header(2,'Done.');
    end
end

cli_header(2, 'Done!');
end

function out=bimod(b,x)
%1 parabola height (unnorm)
%2 center
%3 TF rad
%4 gauss peak height
%5 gauss width
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos
parabola=(1-((x(:)-b(2))./b(3)).^2).^(3/2);
zerosformax=zeros(length(parabola),1);
therm=b(4)*exp(-1*((x(:)-b(2)).^2)./(2.*b(5).^2));
out=real(b(1).*max(zerosformax,parabola)+therm);
end

function out=bimod_sc(b,x)
%1 parabola height (unnorm)
%2 center
%3 TF rad
%4 gauss peak height
%5 gauss width
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos
parabola=(1-((x(:)-b(6))./b(3)).^2).^(3/2);
zerosformax=zeros(length(parabola),1);
therm=b(4)*exp(-1*((x(:)-b(2)).^2)./(2.*b(5).^2));
out=real(b(1).*max(zerosformax,parabola)+therm);
end