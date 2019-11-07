function out=corr_multi_low_mem(corr_opts,counts_txy)
%a low memory implementation of the correlation finder for radial,one D, and three D
%Inputs
%   counts_txy- as a cell array of zxy counts size=[num_counts,3]
%               must be ordered in z for pre window optimzation
%   corr_opts.one_d_window=[[tmin,tmax];[xmin,xmax];[ymin,ymax]];
%pass in a cell array of counts and then calculate the correlations using a
%bin on the fly method
%for the the gxyz we require the difference in the other axes to be less
%than a threshold
%improvements
% implement proper windows (same format as txy mask)
% implement 3d cartesian correlator
% implement r,thet,phi
% implement pre delta windowing for time sorted data using binary search


updates=10;

%find un-norm g2
shots =size(counts_txy,2);
update_interval=ceil(shots/updates);
parfor_progress_imp(ceil(shots/update_interval));



[rad_bins,corr_opts.redges]=histcounts([],corr_opts.redges);
[one_d_bins,corr_opts.one_d_edges]=histcounts([],corr_opts.one_d_edges);
one_d_bins=zeros(shots,size(one_d_bins,2));

%three_d_bins=zeros([shots,size(three_d_bins)]);
if corr_opts.do_three_d_corr
    three_d_bins=histcn([nan,nan,nan],corr_opts.three_d_edges{:});
    size_single_three_d_bins=size(three_d_bins);
    three_d_bins=zeros([size(three_d_bins),shots]);
end
rad_bins=zeros(shots,size(rad_bins,2));
pairs_count=zeros(1,shots);
for shotnum=1:shots
    shot_txy=counts_txy{shotnum};
    num_counts_shot=size(shot_txy,1);
    % full number of pairs in the shot
    pairs_count(shotnum)=num_counts_shot^2 -num_counts_shot;
    for ii=1:num_counts_shot
        %if each shot is sorted in time then this will only produce counts that are one direction in time
        %pre mask optimzation here
        delta=shot_txy(ii,:)-shot_txy(ii+1:num_counts_shot,:);
        
        %one dimensional cartesian correlations
        %mask for the 1d correlation, this is pretty messy..
        switch corr_opts.one_d_dimension
            case 1
                one_d_mask_pos=delta(:,2)>corr_opts.one_d_window(2,1) & delta(:,2)<corr_opts.one_d_window(2,2) &...
                delta(:,3)>corr_opts.one_d_window(3,1) & delta(:,3)<corr_opts.one_d_window(3,2);
                one_d_mask_neg=-delta(:,2)>corr_opts.one_d_window(2,1) & -delta(:,2)<corr_opts.one_d_window(2,2) &...
                -delta(:,3)>corr_opts.one_d_window(3,1) & -delta(:,3)<corr_opts.one_d_window(3,2);
            case 2 
                one_d_mask_pos=delta(:,1)>corr_opts.one_d_window(1,1) & delta(:,1)<corr_opts.one_d_window(1,2) &...
                delta(:,3)>corr_opts.one_d_window(3,1) & delta(:,3)<corr_opts.one_d_window(3,2);
                one_d_mask_neg=-delta(:,1)>corr_opts.one_d_window(1,1) & -delta(:,1)<corr_opts.one_d_window(1,2) &...
                -delta(:,3)>corr_opts.one_d_window(3,1) & -delta(:,3)<corr_opts.one_d_window(3,2);
            case 3
                one_d_mask_pos=delta(:,1)>corr_opts.one_d_window(1,1) & delta(:,1)<corr_opts.one_d_window(1,2) &...
                delta(:,2)>corr_opts.one_d_window(2,1) & delta(:,2)<corr_opts.one_d_window(2,2);
                one_d_mask_neg=-delta(:,1)>corr_opts.one_d_window(1,1) & -delta(:,1)<corr_opts.one_d_window(1,2) &...
                -delta(:,2)>corr_opts.one_d_window(2,1) & -delta(:,2)<corr_opts.one_d_window(2,2);
        end
        one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(delta(one_d_mask_pos,corr_opts.one_d_dimension),corr_opts.one_d_edges);
        one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(-delta(one_d_mask_neg,corr_opts.one_d_dimension),corr_opts.one_d_edges);
        
        %three dimensional cartesian correlations
        %this is very slow, and is dominated by some constant time overhead
        %it would be better to have this only called once (with all(or a lot of) the deltas)
        if corr_opts.do_three_d_corr & size(delta,1)~=0
            three_d_increment=histcn(delta,corr_opts.three_d_edges{:});
            three_d_bins(:,:,:,shotnum)=three_d_bins(:,:,:,shotnum)+three_d_increment;
            three_d_increment=histcn(-delta,corr_opts.three_d_edges{:});
            three_d_bins(shotnum,:,:,:)=three_d_bins(shotnum,:,:,:)+reshape(three_d_increment,[1,size_single_three_d_bins]);
        end

        %%
        
        %radial correlations
        %because the radial correlation is inherently symsetric we can just multiply by 2
        rad_increment=histcounts(sqrt(sum(delta.^2,2)),corr_opts.redges);
        %sum(rad_increment)
        rad_bins(shotnum,:)=rad_bins(shotnum,:)+rad_increment;
    end
    %pairs_count_bluk(shotnum)
    %pairs_count(shotnum)
    if mod(shotnum,update_interval)==0
        parfor_progress_imp;
    end
end
parfor_progress_imp(0);

out.rad_centers=(corr_opts.redges(2:end)+corr_opts.redges(1:end-1))/2;
out.x_centers=(corr_opts.one_d_edges(2:end)+corr_opts.one_d_edges(1:end-1))/2;
out.pairs=sum(pairs_count);

sub_index=[1,2,3];
sub_index(corr_opts.one_d_dimension)=[];
one_d_volume=corr_opts.one_d_window(sub_index,1)-corr_opts.one_d_window(sub_index,2);
one_d_volume=prod(one_d_volume)*(corr_opts.one_d_edges(2:end)-corr_opts.one_d_edges(1:end-1));
out.one_d_corr_density=sum(one_d_bins,1)./(one_d_volume.*out.pairs);
if ~(isnan(corr_opts.one_d_smoothing) || corr_opts.one_d_smoothing==0)
    out.one_d_corr_density=gaussfilt(out.x_centers,out.one_d_corr_density,corr_opts.one_d_smoothing);
end


rad_volume=(4/3)*pi*(corr_opts.redges(2:end).^3-corr_opts.redges(1:end-1).^3);
out.rad_corr_density=sum(rad_bins,1)./(rad_volume.*out.pairs);
if ~(isnan(corr_opts.rad_smoothing) || corr_opts.rad_smoothing==0)
    out.rad_corr_density=gaussfilt(out.rad_centers,out.rad_corr_density,corr_opts.rad_smoothing);
end

end