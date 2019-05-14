function [components,details]=dominant_freq_components(tdat,xdat,options)
% options.components_diff_freq
% options.num_components  in multiples of the inverse time range of the input
% options.components_min_amp in multiples of the input std

% TODO
% detrend x data 

tdat=tdat(:);
xdat=xdat(:);

if ~isfield(options,'num_components') || isnan(options.num_components)
    options.num_components=1e9;
end

if ~isfield(options,'components_diff_freq') || isnan(options.components_diff_freq)
    options.components_diff_freq=5;
end

if ~isfield(options,'freq_limits') || sum(isnan(options.freq_limits))~=0
    options.freq_limits=[0,inf];
end

if ~isfield(options,'components_min_amp') || isnan(options.components_min_amp)
    options.components_min_amp=1e-4;
end

mean_xdat=mean(xdat);
std_xdat=std(xdat);
fft_dat=fft_tx(tdat,xdat-mean_xdat,'window','chebyshev','win_param',{300},'padding',100);

%mask the fft data to lie within freq_lims
fft_idx_lims=fast_sorted_mask(fft_dat(1,:),options.freq_limits(1),options.freq_limits(2));
fft_dat=fft_dat(:,fft_idx_lims(1):fft_idx_lims(2));


% find peaks that are min_peak_factor*xstd and seperated by at least a few times the time resolution
min_pk_sep=options.components_diff_freq/diff(tdat([1,end])); %peak sep in hz
min_pk_sep_idx=round(min_pk_sep/diff(fft_dat(1,1:2))); %peak sep in fft bins
[pks_unsorted,pks_idx] = findpeaks(abs(fft_dat(2,:)),...
    'MinPeakHeight',std_xdat*options.components_min_amp,...
    'MinPeakDistance',min_pk_sep_idx,...
    'NPeaks',options.num_components);
pks_freq=fft_dat(1,pks_idx);
%amplitude sort the peaks from highest to smallest amplitude
[~,sort_order]=sort(pks_unsorted,'descend');
components=[];
components.amp=pks_unsorted(sort_order);
components.freq=pks_freq(sort_order);
components.phase=angle(fft_dat(2,pks_idx))+pi/2;
%make into col vec
components.freq=components.freq(:);
components.amp=components.amp(:);
components.phase=components.phase(:);

details=[];
details.fft_dat=fft_dat;
details.mean_xdat=mean_xdat;
details.std_xdat=std_xdat;

end