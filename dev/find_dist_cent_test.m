function dumy = find_dist_cent_test(data, opts_cent)
% function [bec_centres,bec_widths,bec_counts,centre_OK] = find_dist_cen(data, opts_cent)
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
 data=hotspot_mask(data); %mask outhot spots
num_thr = 250;
num_shots = num_thr;%length(data.shot_num);
bec_centres = zeros(num_shots, 3);
bec_widths= zeros(num_shots, 3);
bec_centres_avg = zeros(num_shots, 3);
bec_centres_fit = zeros(num_shots, 3);
bec_centres_fit_bimod = zeros(num_shots, 3);
bec_centres_avg_3d = zeros(num_shots, 3);
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

threshold = [linspace(21,600,num_thr).*1e3;linspace(15,130,num_thr).*1e3;linspace(15,130,num_thr).*1e3]; %now in set in Hz
min_threshod = [0,4,0].*1e3;%[20,11,10.2].*1e3;%[12,11,10.2].*1e3;%[0,0,0];%[1e4,3e3,3e3];%[2e4,5e3,9e3];
threshold3d = linspace(0,1.5,num_thr).*1e11; %now in Hz/m^2
bin_size  = [20e-5,3e-3,3e-3];%3e-5 * [1, 10, 10];%
sigma = [3e-5,12e-5,12e-5];
bin_edges_3d = [];
for axis = 1:3
    this_edge = (lims(axis, 1):bin_size(axis):lims(axis, 2)).';
    bin_edges_3d{axis} = this_edge;
    bin_centers_3d{axis} =(this_edge(1:end-1)+this_edge(2:end))./2;
end
figure(1)
clf
for this_idx = 1:num_thr%num_shots % Loop over all QD shots
    %         trim_counts(this_idx) = length(trim_txy);
    this_txy = data.masked.counts_txy{45};
    if size(this_txy,2) ~= 3
        dumy = 0;
    end
    trim_txy = masktxy_square(this_txy, lims);
    bec_counts(this_idx) = size(trim_txy,1);
%     count3d = histcn(trim_txy,bin_edges_3d{1},bin_edges_3d{2},bin_edges_3d{3});
%     flux3d = count2flux(count3d,bin_edges_3d);
%     mask3d = flux3d > threshold3d(this_idx);
%     flux3d(mask3d) = nan;
%     [X, T, Y] = meshgrid(bin_centers_3d{2},bin_centers_3d{1},bin_centers_3d{3});
%     bec_centres_avg_3d(this_idx, 1) = nansum(flux3d.*T,'all')/nansum(flux3d,'all');
%     bec_centres_avg_3d(this_idx, 2) = nansum(flux3d.*X,'all')/nansum(flux3d,'all');
%     bec_centres_avg_3d(this_idx, 3) = nansum(flux3d.*Y,'all')/nansum(flux3d,'all');
%     temp1 = bin_edges_3d{1};
%     temp2 = bin_edges_3d{2};
%     figure(19)
%     surf(temp2(1:end-1),temp1(1:end-1),nanmean(flux3d(:,:,:),3))
%     shading flat
%     find centres with the hist threshold method
%     if opts_cent.visual && this_idx == 1
%         clf;
%     end
    % fancy centering
    for axis = 1:3
        this_axis = trim_txy(:, axis);
        this_sigma = sigma(axis);
%         bin_edges = lims(axis, 1):bin_size(axis):lims(axis, 2);
        count_hist = smooth_hist(this_axis,'sigma',this_sigma);
        flux = count_hist.count_rate.smooth;%hist_adaptive_method(this_axis, bin_edges', 0);
        bin_centres = count_hist.bin.centers;%0.5 * (bin_edges(2:end) + bin_edges(1:end-1));
        outer_cut = ~or((max(bin_centres)-bin_centres)./range(bin_centres)<0.05,(-min(bin_centres)+bin_centres)./range(bin_centres)<0.05);
        flux = flux(outer_cut);
        bin_centres = bin_centres(outer_cut);
%         flux = flux(2:end-1);
        mask = flux > threshold(axis,this_idx);
        locs = find(mask);
        margins = [min(locs(2:end-1)), max(locs(2:end-1))];
        if isempty(margins)
            centre_OK(this_idx) = 0;
            bec_centres(this_idx, axis) = nan;
            bec_widths(this_idx, axis) = nan;
%             mask = ~mask;
        else
            t_margins = bin_centres(margins);
            bec_centres(this_idx, axis) = mean(t_margins);
            bec_widths(this_idx, axis) = diff(t_margins);
        end
        mask = flux > threshold(axis,this_idx) | flux < min_threshod(axis);
%             bec_centres_avg(this_idx,axis) = nansum(flux(~mask).*bin_centres(~mask))./nansum(flux(~mask));
            bec_centres_avg(this_idx,axis) = trapz(bin_centres(~mask),flux(~mask).*bin_centres(~mask))./trapz(bin_centres(~mask),flux(~mask));
            bec_widths_avg(this_idx,axis) = sqrt(nansum(flux(~mask).*(bin_centres(~mask)-bec_centres_avg(this_idx,axis)).^2)./nansum(flux(~mask)));
            
            centre_OK(this_idx) = 1;
%             parab = @(b,x) if
            fun1d =  @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2));
            if axis>1
                inital_guess=[2.5e5,bec_centres_avg(this_idx,axis),0.01];
                inital_guess_b=[2.5e5,bec_centres_avg(this_idx,axis),0.01,2.5e5,0.01,bec_centres_avg(this_idx,axis)];
            else
                inital_guess=[12e5,bec_centres_avg(this_idx,axis),0.001];
                inital_guess_b=[12e5,bec_centres_avg(this_idx,axis),0.001,12e5,0.001,bec_centres_avg(this_idx,axis)];
            end
            fo = statset('TolFun',10^-6,...
                'TolX',1e-4,...
                'MaxIter',1e4,...
                'UseParallel',1);
            fit=fitnlm(bin_centres(~mask),flux(~mask),...
                fun1d,...
                inital_guess,...
                'Options',fo);
            fit_both=fitnlm(bin_centres(~mask),flux(~mask),...
                @bimod,...
                inital_guess_b,...
                'Options',fo);
            fit_params = fit.Coefficients.Estimate;
            fit_params_both = fit_both.Coefficients.Estimate;
            bec_centres_fit(this_idx,axis) = fit_params(2);
            bec_centres_fit_bimod(this_idx,axis) = fit_params_both(6);
            bec_widths_fit(this_idx,axis) = fit_params(3);
            bec_widths_fit_bimod(this_idx,axis) = fit_params_both(3);
%         end
        if opts_cent.visual > 1 && this_idx == 100
            figure(1)
            subplot(3, 1, axis);
            plot(bin_centres, flux, 'k')
            hold on
            plot(bin_centres(~mask), flux(~mask), 'b')
        end
    end
end

stfig('centering test x');
clf
plot(threshold(2,:)./1e3,bec_centres(:,2))
hold on
plot(threshold(2,:)./1e3,bec_centres_avg(:,2))
plot(threshold(2,:)./1e3,bec_centres_fit(:,2))
% plot(threshold(2,:)./1e3,bec_centres_fit_bimod(:,2))
% plot(threshold(2,:)./1e3,bec_centres_avg_3d(:,2))
legend('margin x','weighted x','fit x','fit bimod x')
xlabel('threshold value (\(mm^{-1}\))')
ylabel('x center (m)')
% xlim([0,500])
ylim([-4e-3,-1e-3])

stfig('centering test y');
clf
plot(threshold(3,:)./1e3,bec_centres(:,3))
hold on
plot(threshold(3,:)./1e3,bec_centres_avg(:,3))
plot(threshold(3,:)./1e3,bec_centres_fit(:,3))
% plot(threshold(3,:)./1e3,bec_centres_fit_bimod(:,3))
% plot(threshold(3,:)./1e3,bec_centres_avg_3d(:,3))
legend('margin y','weighted y','fit y','fit bimod y')
xlabel('threshold value (\(mm^{-1}\))')
ylabel('y center (m)')
% xlim([0,500])
ylim([3e-3,7e-3])

stfig('centering test z');
clf
plot(threshold(1,:)./1e3,bec_centres(:,1))
hold on
plot(threshold(1,:)./1e3,bec_centres_avg(:,1))
plot(threshold(1,:)./1e3,bec_centres_fit(:,1))
% plot(threshold(1,:)./1e3,bec_centres_fit_bimod(:,1))
% plot(threshold(1,:)./1e3,bec_centres_avg_3d(:,1))
legend('margin z','weighted z','fit z','fit bimod z')
xlabel('threshold value (kHz)')
ylabel('z center (s)')
ylim([3.8503 3.8513])

% stfig('widths')
% clf
% subplot(3,1,1)
% plot(threshold(2,:)./1e3,bec_widths(:,2))
% hold on
% plot(threshold(2,:)./1e3,bec_widths_avg(:,2))
% plot(threshold(2,:)./1e3,bec_widths_fit(:,2))
% % plot(threshold(2,:)./1e3,bec_widths_fit_bimod(:,2))
% % plot(threshold(2,:)./1e3,bec_centres_avg_3d(:,2))
% legend('margin x','weighted x','fit x','fit bimod x')
% xlabel('threshold value (\(mm^{-1}\))')
% ylabel('x width (m)')
% xlim([0,500])
% subplot(3,1,2)


% stfig('cent test x/y 3d');
% clf
% plot(threshold3d./1e3,bec_centres_avg_3d(:,2:3))
% legend('x','y')
% xlabel('threshold value (kHz/\(m^2\))')
% ylabel('x/y center (m)')
% 
% stfig('cent test z 3d');
% clf
% plot(threshold3d./1e3,bec_centres_avg_3d(:,1))
% xlabel('threshold value (kHz/\(m^2\))')
% ylabel('z center (s)')

% if opts_cent.visual
%    for splt = 1:3
%       subplot(3,1,splt)
%       ylabel('Counts')
%    end
% end

% if opts_cent.visual
%    f=stfig('BEC centering');
% %     f=sfigure(2);
% %    clf
%    subplot(2,1,1)
%    plot(bec_centres-mean(bec_centres),'.')
%    title('BEC centre deviation')
%    xlabel('Shot number')
%    ylabel('Centre offset (mm)')
%    subplot(2,1,2)
%    plot(bec_widths)
%    title('BEC width')
%    xlabel('Shot number')
%    ylabel('widths')
%    if opts_cent.savefigs
%     cli_header(1,'Saving images...');
%     saveas(f,fullfile(opts_cent.data_out,'centering.fig'));
%     saveas(f,fullfile(opts_cent.data_out,'centering.svg'));
%     cli_header(2,'Done.');
%    end
% end

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
parabola=(1-((x(:)-b(6))./b(3)).^2).^(3/2);
zerosformax=zeros(length(parabola),1);
therm=b(4)*exp(-1*((x(:)-b(2)).^2)./(2.*b(5).^2));
out=real(b(1).*max(zerosformax,parabola)+therm);
end