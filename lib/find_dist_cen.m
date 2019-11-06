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
    trim_txy = masktxy_square(this_txy, lims);
    bec_counts(this_idx) = size(trim_txy,1);
    % find centres with the hist threshold method
    if opts_cent.visual && this_idx == 1
        clf;
    end
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
if opts_cent.visual
   for splt = 1:3
      subplot(3,1,splt)
      ylabel('Counts')
   end
end

if opts_cent.visual
   f=stfig('BEC centering');
%     f=sfigure(2);
%    clf
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
   if opts.savefigs
    cli_header(1,'Saving images...');
    saveas(f,fullfile(opts_cent.data_out,'centering.fig'));
    saveas(f,fullfile(opts_cent.data_out,'centering.svg'));
    cli_header(2,'Done.');
   end
end

cli_header(2, 'Done!');
end