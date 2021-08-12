function h_data = show_txy_raw(data_input,varargin)
% Returns a structure with fields:
% h_data.flux3 = bin_flux;
% h_data.counts3 = voxel_counts;
% h_data.volumes = bin_volumes;
% h_data.edges{axis} = bin_edges{axis};
% h_data.centres{axis} = bin_cents{axis};
% h_data.flux_1d{axis} = col_vec(squeeze(sum(voxel_counts,ax_excl)));%./col_vec(diff(bin_edges{axis}));
% h_data.flux_mean{axis} = squeeze(nanmean(bin_flux,ax_excl));
% h_data.edges_2d{excl_axis}{1} = h_data.edges{ax_idxs(1)};
% h_data.edges_2d{excl_axis}{2} = h_data.edges{ax_idxs(2)};
% h_data.counts_2d{excl_axis} = squeeze(sum(voxel_counts,excl_axis));
% default display options

    p = inputParser;
    addParameter(p,'num_bins',[100,100,100]);
    addParameter(p,'blur_size',1);
    addParameter(p,'draw_plots',true);
    addParameter(p,'log_plot',false);
    addParameter(p,'log_hist',false);
    addParameter(p,'centre_bec',0);
    addParameter(p,'gamma',1);
    addParameter(p,'min_counts',5);
    addParameter(p,'verbose',true);
    addParameter(p,'lims',[nan,nan
                        -3e-2,3e-2
                        -3e-2,3e-2]);
    addParameter(p,'label','');
    addParameter(p,'txy2v',false);
    addParameter(p,'txy',false);
    addParameter(p,'t_COM',nan);
    addParameter(p,'auto_align',false);
    addParameter(p,'v_mask',nan);
    addParameter(p,'keep_3d',false);
    addParameter(p,'RotationAngle',0);
    addParameter(p,'FontSize',18);
    
    parse(p,varargin{:});
    num_bins = p.Results.num_bins;
    blur_size = p.Results.blur_size;
    centre_bec = p.Results.centre_bec;
    draw_plots = p.Results.draw_plots;
    log_plot= p.Results.log_plot;
    log_hist = p.Results.log_hist;
    verbose = p.Results.verbose;
    label = p.Results.label;
    min_counts= p.Results.min_counts;
    gamma = p.Results.gamma;
    lims = p.Results.lims;
    auto_align = p.Results.auto_align;
    txy2v = p.Results.txy2v;
    t_COM = p.Results.t_COM;
    v_mask = p.Results.v_mask;
    keep_3d = p.Results.keep_3d;
    keep_txy = p.Results.txy;
    fontsize =  p.Results.FontSize;
    theta =  p.Results.RotationAngle;
%     data_in = varargin{1};
%     if nargin > 1
%         opts = varargin{2};
%     end
% 
%     if ~isfield(opts,'min_counts')
%         opts.min_counts = 5;
%     end
% 
%     if ~isfield(opts,'centre_bec')
%         opts.centre_bec = false;
%     end
    const = hebec_constants();
    fall_time = sqrt(2*const.fall_distance/const.g0);
    if isstruct(data_input) %passed a bundle of shots
        shots_in = length(data_input.shot_num);
        OK_num = ~isnan(data_input.shot_num);
        OK_counts = ~isnan(data_input.num_counts);
        OK_txy = cellfun(@(x) numel(x) > min_counts, data_input.counts_txy);
        mask = OK_num & OK_counts & OK_txy;
        data_work = struct_mask(data_input,mask);
        num_shots = sum(mask);
        all_together = vertcat(data_work.counts_txy{:});
    elseif ismatrix(data_input) && size(data_input,2) == 3 % passed a raw TXY array
        shots_in = 1;
        num_shots = 1;
        data_work.num_counts = size(data_input,1);
        data_work.shot_num = 1;
        all_together = data_input;
    end
        
    if isnan(lims(1,2))
%         lims(1,2) = max(all_together(:,1));
        lims = getlims(all_together);
    end
    txy = masktxy_square(all_together,lims);
    R = [cos(theta),sin(theta);
        -sin(theta),cos(theta)];
    txy(:,[2,3]) = (R*txy(:,[2,3])')';
    if iscell(txy)
        txy = cell2mat(txy);
    end
    
   
    v_com_xy = [0,0];
    if txy2v
        if isnan(t_COM)
            warning('Release time not specified, will make a guess')
            t_COM = mean(txy(:,1)); %this will be incorrect
        end
        txy = txy2vzxy(txy,'t_COM',t_COM);
        v_com_xy = mean(txy(:,2:3));
        txy(:,2:3) = txy(:,2:3) - v_com_xy;
        h_data.pulse_cen = mean(txy);
        h_data.pulse_std = std(txy);
        if ~isnan(v_mask)
           txy = mask_square(txy,v_mask); 
        end
    end
    bounds = getlims(txy);
    bounds(2:3,:) = bounds(2:3,:) + v_com_xy';
    if auto_align
        C = cov(txy(:,2:3));
        [evec,~] = eig(C);
        [~,y_index] = max(abs(evec*[0;1]));
        if y_index == 1 % first evec is in y dir
            txy(:,2:3) = txy(:,[3,2])*evec;
        elseif y_index == 2
            txy(:,2:3) = txy(:,2:3)*evec;
        end
    end
    if verbose
        fprintf('DISPLAYING TXY DATA:\n')
        fprintf('%u shots passed\n',shots_in)
        fprintf('%u total counts\n',sum(data_work.num_counts))
        fprintf('Mean %.2f std %.2f stderr %.2f\n',mean(data_work.num_counts),std(data_work.num_counts),std(data_work.num_counts)/length(data_work.shot_num))
        fprintf('%u windowed counts\n',size(txy,1))
%         fprintf('t_COM (%.3f,%.3f,%.3f)\n',h_data.pulse_cen');
%         fprintf('VAR (%.3f,%.3f,%.3f)\n',h_data.pulse_std');
    end %verbose 
    
    
    bounds = lims;
    t_edges = linspace(bounds(1,1),bounds(1,2)*1.01,num_bins(1)+1);
    x_edges = linspace(bounds(2,1),bounds(2,2)*1.01,num_bins(2)+1);
    y_edges = linspace(bounds(3,1),bounds(3,2)*1.01,num_bins(3)+1);
    [voxel_counts,bin_edges,bin_cents,~] = histcn(txy,t_edges,x_edges,y_edges);
    bin_volumes = (OuterProduct(diff(t_edges)'.*diff(x_edges),diff(y_edges')));
    voxel_counts = voxel_counts/num_shots; % av counts per voxel
    bin_flux = voxel_counts./bin_volumes;
    
    h_data.flux3 = bin_flux;
    h_data.counts3 = voxel_counts;
    h_data.volumes = bin_volumes;
    h_data.mean_counts = sum(h_data.counts3,'all');
    
    if txy2v
        profile_labels = {'Z','X','Y'};    
        axis_labels = {'$v_z$ (m/s)','$v_x$ (m/s)','$v_y$ (m/s)'};
    else
         profile_labels = {'T','X','Y'};    
         axis_labels = {'T (s)','X (cm)','Y (cm)'};
    end
    
    
    h_data.edges = cell(3,1);
    h_data.counts_1d = cell(3,1);   
    h_data.centres = cell(3,1);
    h_data.flux_1d = cell(3,1);
    h_data.num_shots = num_shots;
    if size(txy,1) < min_counts
        warning('Counts too low');
    end
    
    ax_order = 1:3; % to align the X and Y axes in the figure
    hist_ax_order = [1 3 2];
    
    for axis = 1:3
        ax_excl = 1:3;
        ax_1d = ax_order(axis);
        ax_excl(axis) = [];
        h_data.edges{axis} = bin_edges{axis};
        h_data.centres{axis} = bin_cents{axis};
        h_data.counts_1d{axis} = histcounts(txy(:,axis),bin_edges{axis});
%         h_data.counts_1d{axis} = col_vec(squeeze(sum(voxel_counts,ax_excl)));
        h_data.mean_1d = h_data.counts_1d{axis}/num_shots;
        h_data.flux_1d{axis} = h_data.counts_1d{axis}./diff(h_data.edges{axis});
        if axis ~= 1
            h_data.flux_1d{axis} = h_data.flux_1d{axis}/diff(bounds(1,:));
        end
    end
    
    
    h_data.X_c = [0,0,0];
    thr = .00;
%     for ax=1:3
%         mask = rescale(h_data.flux_1d{ax}) > thr;
%         X=h_data.centres{ax};
%         h_data.X_c(ax) = median(X(mask));
%     end
    
    h_data.counts_2d = cell(3,1);
    for txy_count = 1:3
        h_data.counts_2d{txy_count} = squeeze(sum(voxel_counts,txy_count));        
    end
    h_data.txy = txy;
    %% PLOTTING
    grid_height = 2;
    grid_width = 3;
    
    subplot_indices = {[10],[2,3],[4,7],[5,6,8,9],[11],[12]};
    subplot_counter = 0;
    
    if draw_plots
        stfig(sprintf('BEC TXY display %s',label));
        clf
%         tiledlayout(grid_height,grid_width)
        for axis_count = 1:3
            subplot_counter = subplot_counter + 1;
            subplot(4,3,subplot_indices{subplot_counter})
            axis = ax_order(axis_count);
%             nexttile
            
            hold on
            if axis ~=1
                fact = 1e2;
            else
                fact = 1;
            end
            plot(fact*h_data.centres{axis},(h_data.flux_1d{axis})/1e3)
%             title(sprintf('%s profile',profile_labels{axis}))
            if axis_count == 1
                ylabel('Flux (Hz)')
            end
            if axis_count ~= 1
                ylabel('Flux (m$^2$ s)$^{-1}$')
            end
%             yticks([])
            
            xlim(fact*[min(lims(axis,:)),max(lims(axis,:))])
            if log_plot
                set(gca,'Yscale','log')
            end
            if axis_count == 3
                camroll(90)
            end
%             if axis_count ~=2
%                 set(gca,'XAxisLocation','top')
%             end
            if axis_count == 1
                xticks([])
                camroll(90)
            else
                xlabel(sprintf('%s',axis_labels{axis}))
            end
            set(gca,'FontSize',fontsize)
        end
%         h_data.flux_1d = h_data.flux_1d{1,3,2};
        
        
        axis_labels = {'T(s)','X (cm)','Y (cm)'};
        for txy_count = 1:3
            kept_axes = ax_order;
            ax_sel = hist_ax_order(txy_count);
            kept_axes(txy_count) = [];
            
%             nexttile
            subplot_counter = subplot_counter + 1;
            subplot(4,3,subplot_indices{subplot_counter})
            if blur_size == 1
                V = h_data.counts_2d{ax_sel};
            else
                V = imgaussfilt(h_data.counts_2d{excl_axis},blur_size);
            end
            if txy_count == 1
                V = V';
                kept_axes = fliplr(kept_axes);
                ax_labels = {axis_labels{2},axis_labels{3}};
                edges_2d{1} = h_data.edges{kept_axes(1)};
                edges_2d{2} = h_data.edges{kept_axes(2)};
                fact = 1;
            else
                fact = 1e2;
                ax_labels = {axis_labels{txy_count},axis_labels{1}};
                kept_axes(2) = txy_count;
                edges_2d{1} = h_data.edges{kept_axes(1)};
                edges_2d{2} = fact*h_data.edges{txy_count};
            end
            if log_hist
                V = log(V);
            end


            imagesc(edges_2d{2},edges_2d{1},V)
            
            if txy_count == 1
                xticks([])
            else
                xlabel(sprintf('%s',ax_labels{1}))
                
%                 if txy_count ~= 1
%                     ylabel(sprintf('%s',ax_labels{2}))
%                 end
            end
            if txy_count ~=3
                yticks([])
            else
                set(gca,'YAxisLocation','right')
                ylabel('T (s)')
            end
            xlim(fact*lims(kept_axes(2),:))
            ylim(lims(kept_axes(1),:))
            set(gca,'YDir','normal')
            
            set(gca,'FontSize',fontsize)
        end
        
             
%         if plot_3d > 1
%             stfig('3D BEC display');
%             clf
%             nexttile
%         end
    if ~keep_3d
        h_data = rmfield(h_data,'flux3');
        h_data = rmfield(h_data,'counts3');
        h_data = rmfield(h_data,'volumes');
    end
    
    if ~keep_txy
        h_data = rmfield(h_data,'txy');
    end
        colormap(viridis)
        suptitle(label)
        
    end
    
end