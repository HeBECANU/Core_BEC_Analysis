function out=corr_radial(corr_opts,counts)
%a low/high memory implementation of the one D correlation
%Inputs
%   counts_txy- as a cell array of zxy counts size=[num_counts,3]
%               must be ordered in z for pre window optimzation
%   corr_opts.one_d_window=[[tmin,tmax];[xmin,xmax];[ymin,ymax]]; %only used for the premask
%   corr_opts.cl_or_bb
%   corr_opts.one_d_edges
%   corr_opts.one_d_dimension
%   corr_opts.attenuate_counts- value 0-1 the ammount of data 
%   corr_opts.do_pre_mask- logical, use fast sorted mask to window to a range of 
%                           one_d_window in the sorted axis arround each count
%       corr_opts.sorted_dir - must be passed when corr_opts.do_pre_mask used
%pass in a cell array of counts and then calculate the correlations using a
%bin on the fly method
%for the the gxyz we require the difference in the other axes to be less
%than a threshold

%improvements
% [x] implement pre delta windowing for time sorted data using binary search
% input checking
% [x] kill counts, to speed things up when there is a lot of data remove some fraction

if size(counts,1)==2
    corr_opts.between_sets = true;
else
    corr_opts.between_sets = false;
end

if ~isfield(corr_opts,'attenuate_counts')
    corr_opts.attenuate_counts=1;
else
    if corr_opts.attenuate_counts>1 || corr_opts.attenuate_counts<0
        error('invalid attenuation value, pass a value between 0 and 1\n')
    end
end

if ~isfield(corr_opts,'cl_or_bb')
    error('corr_opts.cl_or_bb not specified!do you want co-linear or back-to-back? ')
end

if isfield(corr_opts,'do_pre_mask') && ~isfield(corr_opts,'sorted_dir')
        error('you must pass corr_opts.sorted_dir with corr_opts.do_pre_mask')
elseif ~isfield(corr_opts,'do_pre_mask')
    corr_opts.do_pre_mask=false;
end

if ~isfield(corr_opts,'redges')
    error('no radial edges specified')
else
    if min(corr_opts.redges)<0
        error('redges should only be positive')
    end
end

if ~isfield(corr_opts,'progress_updates')
    corr_opts.progress_updates=50;
end

num_counts=cellfun(@(x)size(x,1),counts);
if ~isfield(corr_opts,'low_mem') || isnan(corr_opts.low_mem)
    mem_temp=memory;
    max_arr_size=floor(mem_temp.MaxPossibleArrayBytes/(8*2)); %divide by 8 for double array size and 2 as a saftey factor

    if corr_opts.do_pre_mask
        premask_factor=1;  %premasking currently uses the same amount of memory 
        %premasking could in principle dramaticaly reduce the number of pairs that are stored in memory
        %this might be a factor of ~1/100 to be implemnted in future
    else
        premask_factor=1;
    end
    %use the corr_opts.norm_samp_factor so that both corr/uncorr parts are calculated with the same function
    corr_opts.low_mem=((max(num_counts)*corr_opts.norm_samp_factor*premask_factor)^2)>max_arr_size;
    if corr_opts.print_update
    if corr_opts.low_mem,fprintf('auto detection:using the low memory option\n'), end
    if ~corr_opts.low_mem,fprintf('auto detection:using the high memory option\n'), end
    end
end

% convert edges to col vec
corr_opts.redges=col_vec(corr_opts.redges);

%set the prewindow to be the max radial edge
prewindow=[-1,1]*max(corr_opts.redges);

updates=corr_opts.progress_updates; %number of updates in the progress bar to give the user, can slow thigs down if too high
shots =size(counts,2);
update_interval=ceil(shots/updates);
if corr_opts.print_update
    parfor_progress_imp(ceil(shots/update_interval));
end

pairs_count=zeros(1,shots);
delta_multiplier=[1,-1];
delta_multiplier=delta_multiplier(1+corr_opts.cl_or_bb); %gives 1 when cl and -1 when bb
       
[rad_bins,corr_opts.redges]=histcounts([],corr_opts.redges);
out.rad_centers=(corr_opts.redges(2:end)+corr_opts.redges(1:end-1))/2; %this is slightly wrong it should be r^3 weighted

if corr_opts.low_mem %calculate with the low memory mode
    rad_bins=col_vec(rad_bins);
    for shotnum=1:shots
        if corr_opts.between_sets
            shot_txy_1=counts{1,shotnum};
            shot_txy_2=counts{2,shotnum};
            num_counts_shot_1=num_counts(1,shotnum);
            num_counts_shot_2=num_counts(2,shotnum);
        else
            shot_txy=counts{shotnum};
            num_counts_shot=num_counts(shotnum);
        end
        if corr_opts.attenuate_counts~=1 %randomly keep corr_opts.attenuate_counts fraction of the data
            if corr_opts.between_sets
                mask1=rand(num_counts_shot_1,1)<corr_opts.attenuate_counts;
                mask2=rand(num_counts_shot_2,1)<corr_opts.attenuate_counts;
                shot_txy_1=shot_txy_1(mask1,:);
                shot_txy_2=shot_txy_2(mask2,:);
                num_counts_shot_1=size(shot_txy_1,1);
                num_counts_shot_2=size(shot_txy_2,1);
            else
                mask=rand(num_counts_shot,1)<corr_opts.attenuate_counts;
                shot_txy=shot_txy(mask,:);
                num_counts_shot=size(shot_txy,1); %recalulate the number of counts
            end
        end
        if corr_opts.between_sets
            pairs_count(shotnum)=num_counts_shot_1*num_counts_shot_2;
            num_counts_shot = num_counts_shot_1+1;
        else
            % full number of pairs in the shot
            pairs_count(shotnum)=num_counts_shot^2 -num_counts_shot;
        end
        if num_counts_shot<2
            warning('%u counts input\n',num_counts_shot)
        end
        for ii=1:num_counts_shot-1
            %if each shot is sorted in time then this will only produce counts that are one direction in time
            if corr_opts.do_pre_mask && ~corr_opts.between_sets %pre mask optimzation using sortd mado_pre_masksking
                temp_1d_diff=shot_txy(ii+1:num_counts_shot,corr_opts.sorted_dir)...
                    -delta_multiplier*shot_txy(ii,corr_opts.sorted_dir);
                mask_idx=fast_sorted_mask(...
                    temp_1d_diff,...
                    prewindow(1),...
                    prewindow(2));
                delta=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+mask_idx(1):ii+mask_idx(2),:);
            else
                if corr_opts.between_sets
                    delta=shot_txy_1(ii,:)-delta_multiplier*shot_txy_2(:,:);
                else
                    delta=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:);
                end
            end

            %Radial correlations
            rad_delta=sqrt(sum(delta.^2,2));
            % using fast histogram gives speedup of 
            rad_increment=hist_adaptive_method(rad_delta,corr_opts.redges);
            %rad_increment=histcounts(rad_delta,corr_opts.redges)';
            rad_bins=rad_bins+rad_increment*2; 
        end
        if mod(shotnum,update_interval)==0 && corr_opts.print_update
            parfor_progress_imp;
        end
    end%loop over shots
    
else%calculate with the high memory mode
    rad_bins=zeros(shots,size(rad_bins,2));
    parfor shotnum=1:shots
        shot_txy=counts{shotnum};
        num_counts_shot=num_counts(shotnum);
        if corr_opts.attenuate_counts~=1 %randomly keep corr_opts.attenuate_counts fraction of the data
            mask=rand(num_counts_shot,1)<corr_opts.attenuate_counts;
            shot_txy=shot_txy(mask,:);
            num_counts_shot=size(shot_txy,1); %recalulate the number of counts
        end
        if num_counts_shot<2
            warning('%u counts input\n',num_counts_shot)
        end
        % full number of pairs in the shot
        pairs_count(shotnum)=num_counts_shot^2 -num_counts_shot;
        upper_triangle_size=pairs_count(shotnum)/2;
        delta_idx=1;
        delta=nan(upper_triangle_size,3);
        for ii=1:num_counts_shot-1
            %if each shot is sorted in time then this will only produce counts that are one direction in time
            %pre mask optimzation here
            if corr_opts.do_pre_mask  %pre mask optimzation using sortd masking
                    temp_1d_diff=shot_txy(ii+1:num_counts_shot,corr_opts.sorted_dir)...
                        -delta_multiplier*shot_txy(ii,corr_opts.sorted_dir);
                    mask_idx=fast_sorted_mask(...
                        temp_1d_diff,...
                        prewindow(1),...
                        prewindow(2));
                    delta_inc=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+mask_idx(1):ii+mask_idx(2),:);
                else
                    delta_inc=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:);
            end
            delta_inc_size=size(delta_inc,1);
            delta(delta_idx:delta_idx+delta_inc_size-1,:)=delta_inc;
            delta_idx=delta_idx+delta_inc_size; %increment the index
        end   
        %Radial correlations
        rad_delta=sqrt(sum(delta.^2,2));
        rad_increment=hist_adaptive_method(rad_delta,corr_opts.redges);
        %rad_increment=histcounts(rad_delta,corr_opts.redges);
        rad_bins(shotnum,:)=rad_increment*2 ;
        %debug
%         sfigure(2)
%         plot(out.rad_centers,rad_increment)
%         pause(0.1)
%         rad_bins(shotnum,:)=rad_increment*2;
%         size(rad_delta)
%         size(delta)
        
        if mod(shotnum,update_interval)==0 && corr_opts.print_update
            parfor_progress_imp;
        end
    end%loop over shots
    
    % sum up the output of parfor
    rad_bins=col_vec(sum(rad_bins,1));
end %done calculating with either high or low mem

if corr_opts.print_update
    parfor_progress_imp(0);
end


out.pairs=sum(pairs_count);

rad_volume=(4/3)*pi*(corr_opts.redges(2:end).^3-corr_opts.redges(1:end-1).^3);
out.rad_corr_density=rad_bins./(rad_volume.*out.pairs);
if ~(isnan(corr_opts.rad_smoothing) || corr_opts.rad_smoothing==0)
    out.rad_corr_density=gaussfilt(out.rad_centers,out.rad_corr_density,corr_opts.rad_smoothing);
end

end

% if sum(sum(delta,2)==0)>0
%     fprintf('warning duplicate count')
% end
            
% %%you can check that the sorted mask this is equivelent to the below brute approach
% delta_brute=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:); 
% min_lim=corr_opts.one_d_window(corr_opts.sorted_dir,1);
% max_lim=corr_opts.one_d_window(corr_opts.sorted_dir,2);
% mask=delta_brute(:,corr_opts.sorted_dir)>min_lim & delta_brute(:,corr_opts.sorted_dir)<max_lim;
% delta_brute =delta_brute(mask,:);
% if ~isequal(delta_brute,delta)
%     fprintf('warning mask is not giving correct answer')
% end