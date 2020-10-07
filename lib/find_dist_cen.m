function [bec_centres,bec_widths,bec_counts,centre_OK] = find_dist_cen(data, opts_cent)
% Mandatory data inputs:
%     data.counts_txy
%     data.shot_num
% Mandatory Input options:
%     opts.cent.auto
%     opts.cent.crop
%     opts.cent.bin_size
%     opts.cent.threshold
%     opts.partition.t_win
%     opts.visual
cli_header(1,'Finding centres in %s',opts_cent.label);
num_shots = length(data.shot_num);
bec_centres = zeros(num_shots, 3);
bec_widths= zeros(num_shots, 3);
bec_counts = zeros(num_shots, 1);
centre_OK = zeros(num_shots, 1);
% For automagic centre guessing
% temp_txy = data.counts_txy{1};
% temp_t = temp_txy(:, 1);
% t_edges = min(temp_t):opts.partition.t_win:max(temp_t);
% t_cents = 0.5 * (t_edges(2:end) + t_edges(1:end-1));
% temp_profile = hist_adaptive_method(temp_t, t_edges', 1);
% temp_profile = temp_profile(2:end-1);
% [~, ploc] = max(temp_profile);
% t_peak = t_cents(ploc);

lims = opts_cent.crop;
if length(opts_cent.threshold) < 3
   opts_cent.threshold = opts_cent.threshold(1)*[1,1,1]; 
end
% Set the t limits by conversion to k
% lims(1, :) = t_peak + [-1, 1] * (opts.c.hbar * opts.spherify.k_max / opts.c.m_He) / (0.5 * opts.c.g0 * 0.4187);

for this_idx = 1:num_shots % Loop over all QD shots
    %         trim_counts(this_idx) = length(trim_txy);
    this_txy = data.counts_txy{this_idx};
    if size(this_txy,2) ~= 3
        dumy = 0;
    end
    trim_txy = masktxy_square(this_txy, lims);
    bec_counts(this_idx) = size(trim_txy,1);
    % find centres with the hist threshold method

    % fancy centering
    for axis = 1:3
        this_axis = trim_txy(:, axis);
        bin_edges = lims(axis, 1):opts_cent.bin_size(axis):lims(axis, 2);
        flux = hist_adaptive_method(this_axis, bin_edges', 0);
        bin_centres = 0.5 * (bin_edges(2:end) + bin_edges(1:end-1));
        flux = flux(2:end-1);
        mask = flux > opts_cent.threshold(axis);
        locs = find(mask);
        margins = [min(locs(2:end-1)), max(locs(2:end-1))];
        if isempty(margins)
            centre_OK(this_idx) = 0;
        else
            t_margins = bin_centres(margins);
            bec_centres(this_idx, axis) = mean(t_margins);
            bec_widths(this_idx, axis) = diff(t_margins);
            centre_OK(this_idx) = 1;
        end
        if opts_cent.visual > 1
            subplot(3, 1, axis);
            plot(bin_centres, flux, 'k')
            hold on
            plot(bin_centres(mask), flux(mask), 'r')
        end
    end
end
% if opts_cent.visual
%    for splt = 1:3
%       subplot(3,1,splt)
%       ylabel('Counts')
%    end
% end

if opts_cent.visual
   f=stfig('BEC centering');
%     f=sfigure(2);
   clf
    xmin = -3e-3;
    xmax = 3e-3;
    c_edges = 0.5*linspace(xmin,xmax,30);
    w_edges = 1.5*linspace(xmin,xmax,30);

   subplot(3,2,1)
   hold on
   plot(bec_centres(:,1)-mean(bec_centres(:,1)),'r.')
   plot(bec_centres(:,2)-mean(bec_centres(:,2)),'b.')
   plot(bec_centres(:,3)-mean(bec_centres(:,3)),'k.')
   title('BEC centre deviation')
   xlabel('Shot number')
   ylabel('Centre offset from mean')
   legend('T','X','Y')
   subplot(3,2,2)
   hold on
   histogram(bec_centres(:,1)-mean(bec_centres(:,1)),c_edges,'FaceColor','r','FaceAlpha',0.2)
    histogram(bec_centres(:,2)-mean(bec_centres(:,2)),c_edges,'FaceColor','b','FaceAlpha',0.2)
    histogram(bec_centres(:,3)-mean(bec_centres(:,3)),c_edges,'FaceColor','k','FaceAlpha',0.2)
    xlim([min(c_edges),max(c_edges)])
    legend('T','X','Y')
    title('BEC centre deviation')
   ylabel('Centre offset from mean')
   ylabel('Number of shots')
   
   subplot(3,2,3)
   hold on
   plot(bec_widths(:,1)-mean(bec_widths(:,1)),'r.')
   plot(bec_widths(:,2)-mean(bec_widths(:,2)),'b.')
   plot(bec_widths(:,3)-mean(bec_widths(:,3)),'k.')
   legend('T','X','Y')
   title('BEC width variation')
   xlabel('Shot number')
   ylabel('Width variation from mean')
   subplot(3,2,4)
   hold on
    histogram(bec_widths(:,1)-mean(bec_widths(:,1)),w_edges,'FaceColor','r','FaceAlpha',0.2)
    histogram(bec_widths(:,2)-mean(bec_widths(:,2)),w_edges,'FaceColor','b','FaceAlpha',0.2)
    histogram(bec_widths(:,3)-mean(bec_widths(:,3)),w_edges,'FaceColor','k','FaceAlpha',0.2)
    xlim([min(w_edges),max(w_edges)])
    legend('T','X','Y')
    title('BEC width variation')
   ylabel('Width offset from mean')
   ylabel('Number of shots')
   
   %    plot(bec_widths)
   if opts_cent.savefigs
    cli_header(1,'Saving images...');
    saveas(f,fullfile(opts_cent.data_out,'centering.fig'));
    saveas(f,fullfile(opts_cent.data_out,'centering.eps'));
    saveas(f,fullfile(opts_cent.data_out,'centering.svg'));
    cli_header(2,'Done.');
   end
end

cli_header(2, 'Done!');
end