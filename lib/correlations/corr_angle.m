function out=corr_angle(corr_opts,counts)
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
    corr_opts.low_mem=((max(num_counts,[],'all')*corr_opts.norm_samp_factor*premask_factor)^2)>max_arr_size;
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

[ang_bins,corr_opts.edges]=histcounts([],corr_opts.edges);
out.centers=(corr_opts.edges(2:end)+corr_opts.edges(1:end-1))./2;

if corr_opts.low_mem %calculate with the low memory mode
    % the low memory mode is serial and is a bit easier on mem requirements
    ang_bins=col_vec(ang_bins);
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
                dist = randperm(num_counts_shot_1);
                dist = dist(1:ceil(corr_opts.attenuate_counts*num_counts_shot_1));
                mask1 = zeros(num_counts_shot_1,1);
                mask1(dist) = 1;
                mask1 = logical(mask1);
                
                dist = randperm(num_counts_shot_2);
                dist = dist(1:ceil(corr_opts.attenuate_counts*num_counts_shot_2));
                mask2 = zeros(num_counts_shot_2,1);
                mask2(dist) = 1;
                mask2 = logical(mask2);
                
                shot_txy_1=shot_txy_1(mask1,:);
                shot_txy_2=shot_txy_2(mask2,:);
                num_counts_shot_1=size(shot_txy_1,1);
                num_counts_shot_2=size(shot_txy_2,1);
            else
                dist = randperm(num_counts_shot);
                dist = dist(1:ceil(corr_opts.attenuate_counts*num_counts_shot));
                mask = zeros(num_counts_shot,1);
                mask(dist) = 1;
                mask = logical(mask);
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
        %         if num_counts_shot<2
        %             warning('%u counts input\n',num_counts_shot)
        %         end
        for ii=1:num_counts_shot-1
            %if each shot is sorted in time then this will only produce counts that are one direction in time
            if corr_opts.do_pre_mask && ~corr_opts.between_sets %pre mask optimzation using sortd mado_pre_masksking
                temp_1d_diff=shot_txy(ii+1:num_counts_shot,corr_opts.sorted_dir)...
                    -delta_multiplier*shot_txy(ii,corr_opts.sorted_dir);
                mask_idx=fast_sorted_mask(...
                    temp_1d_diff,...
                    prewindow(1),...
                    prewindow(2));
                delta=cos_theta(shot_txy(ii,:),shot_txy(ii+mask_idx(1):ii+mask_idx(2),:));
            else
                if corr_opts.between_sets
                    if isequal(shot_txy_1,shot_txy_2)
                        if ii == num_counts_shot-1
                            continue
                        end
                        delta=cos_theta(shot_txy_1(ii,:),shot_txy_2(ii+1:(num_counts_shot-1),:));
                    else
                        delta=cos_theta(shot_txy_1(ii,:),shot_txy_2(:,:));
                    end
                else
                    delta=cos_theta(shot_txy(ii,:),shot_txy(ii+1:num_counts_shot,:));
                end
            end
            
            %Angular correlations
            ang_delta=delta;%
            % using fast histogram gives speedup of
            ang_increment=hist_adaptive_method(ang_delta,corr_opts.edges);
            %rad_increment=histcounts(rad_delta,corr_opts.redges)';
            if corr_opts.between_sets && ~isequal(shot_txy_1,shot_txy_2)
                ang_bins=ang_bins+ang_increment;
            else
                ang_bins=ang_bins+ang_increment*2;
            end
        end
        if mod(shotnum,update_interval)==0 && corr_opts.print_update
            parfor_progress_imp;
        end
    end%loop over shots
%     toc
else%calculate with the high memory mode
%     tic
%     rad_bins=zeros(shots,size(rad_bins,2));
    error('High mem not available yet')
end %done calculating with either high or low mem

if corr_opts.print_update
    parfor_progress_imp(0);
end


% out.pairs=sum(pairs_count);

rad_volume=(4/3)*pi*(corr_opts.redges(2:end).^3-corr_opts.redges(1:end-1).^3);
% out.rad_corr_density=rad_bins./(rad_volume.*out.pairs);
if isfield(corr_opts,'normalisation_factor')
    norm_factor = corr_opts.normalisation_factor;
else
    norm_factor = 1;
end

out.ang_corr_density=ang_bins./(shots*norm_factor);
if ~(isnan(corr_opts.rad_smoothing) || corr_opts.rad_smoothing==0)
    out.ang_corr_density=gaussfilt(out.centers,out.ang_corr_density,corr_opts.ang_smoothing);
end

end

function angs = cos_theta(P1,P2)
angs = zeros(size(P2,1),1);
for ii = 1:size(P2,1)
    angs(ii,1) = dot(P1,P2(ii,:))/(norm(P1)*norm(P2(ii,:)));
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