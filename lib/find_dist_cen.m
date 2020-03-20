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
data=hotspot_mask(data); %mask outhot spots
num_shots = length(data.shot_num);
bec_centres = zeros(num_shots, 3);
bec_widths= zeros(num_shots, 3);
bec_counts = zeros(num_shots, 1);
centre_OK = zeros(num_shots, 1);

lims = opts_cent.crop;
if length(opts_cent.threshold) < 3
    opts_cent.threshold = opts_cent.threshold(1)*[1,1,1];
end
if length(opts_cent.sigma) < 3
    opts_cent.sigma = opts_cent.sigma(1)*[1,1,1];
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

for this_idx = 1:num_shots % Loop over all shots
    this_txy = data.masked.counts_txy{this_idx};
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
    for axis = 1:3
        this_axis = trim_txy(:, axis);
        this_sigma = opts_cent.sigma(axis);
        this_method = opts_cent.method{axis};
        if isempty(this_axis)
            centre_OK(this_idx) = 0;
            bec_centres(this_idx, axis) = nan;
            bec_widths(this_idx, axis) = nan;
            continue
        end
        count_hist = smooth_hist(this_axis,'sigma',this_sigma);
        flux = count_hist.count_rate.smooth;
        bin_centres = count_hist.bin.centers;
        mask = flux > opts_cent.threshold(axis);
        if all(mask) || sum(~mask)/length(mask)<0.1
            centre_OK(this_idx) = 0;
            bec_centres(this_idx, axis) = nan;
            bec_widths(this_idx, axis) = nan;
            continue
        elseif all(~mask) && strcmp(this_method,'margin')
            this_method = 'average';
        else
            centre_OK(this_idx) = 1;
            switch this_method
                case 'margin'
                    locs = find(mask);
                    margins = [min(locs(2:end-1)), max(locs(2:end-1))];
                    t_margins = bin_centres(margins);
                    bec_centres(this_idx, axis) = mean(t_margins);
                    bec_widths(this_idx, axis) = diff(t_margins);
                case 'average'
                    bec_centres(this_idx,axis) = nansum(flux(~mask).*bin_centres(~mask))./nansum(flux(~mask));
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
        if opts_cent.visual > 1
            subplot(3, 1, axis);
            plot(bin_centres-bec_centres(this_idx,axis), flux, 'k')
            hold on
            plot(bin_centres(mask)-bec_centres(this_idx,axis), flux(mask), 'r')
        end
    end
end

if opts_cent.visual
    for splt = 1:3
        subplot(3,1,splt)
        ylabel('Counts')
    end
end

if opts_cent.visual
    f=stfig('BEC centering');
    clf
    subplot(2,1,1)
    plot(bec_centres-mean(bec_centres),'.')
    title('BEC centre deviation')
    xlabel('Shot number')
    ylabel('Centre offset (mm)')
    subplot(2,1,2)
    plot(bec_widths)
    title('BEC width')
    xlabel('Shot number')
    ylabel('widths')
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