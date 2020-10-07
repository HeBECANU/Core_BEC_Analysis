function h_data = show_txy_raw(data_in,varargin)
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
    addParameter(p,'lims',[0,nan
                        -3e-2,3e-2
                        -3e-2,3e-2]);
    addParameter(p,'label','');
    
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
    
    if isstruct(data_in) %passed a bundle of shots
        shots_in = length(data_in.shot_num);
        OK_num = ~isnan(data_in.shot_num);
        OK_counts = ~isnan(data_in.num_counts);
        OK_txy = cellfun(@(x) numel(x) > min_counts, data_in.counts_txy);
        mask = OK_num & OK_counts & OK_txy;
        data_in = struct_mask(data_in,mask);
        num_shots = sum(mask);
        all_together = vertcat(data_in.counts_txy{:});
    elseif ismatrix(data_in) && size(data_in,2) == 3 % passed a raw TXY array
        shots_in = 1;
        num_shots = 1;
        data_in.num_counts = size(data_in,1);
        all_together = data_in;
    end
        
    if isnan(lims(1,2))
        lims(1,2) = max(all_together(:,1));
    end
    txy = masktxy_square(all_together,lims);
    if iscell(txy)
        txy = cell2mat(txy);
    end
    
   
    
    h_data.pulse_cen = mean(txy);
    
    h_data.pulse_std = std(txy);
    if centre_bec && numel(centre_bec) == 1
        % centre coords automatically (less accurate)
        txy = txy - h_data.pulse_cen;
    elseif centre_bec && numel(centre_bec) == 3 
        % Centre coords set manually
        txy = txy - centre_bec;
    end
    
    if verbose
        fprintf('DISPLAYING TXY DATA:\n')
        fprintf('%u shots passed\n',shots_in)
        fprintf('%u total counts\n',sum(data_in.num_counts))
        fprintf('Mean %.2f std %.2f stderr %.2f\n',mean(data_in.num_counts),std(data_in.num_counts),std(data_in.num_counts)/length(data_in.shot_num))
        fprintf('%u windowed counts\n',size(txy,1))
        fprintf('COM (%.3f,%.3f,%.3f)\n',h_data.pulse_cen');
        fprintf('VAR (%.3f,%.3f,%.3f)\n',h_data.pulse_std');
    end %verbose 
    
    
    bounds = getlims(txy);
    t_edges = linspace(bounds(1,1),bounds(1,2)*1.01,num_bins(1)+1);
    x_edges = linspace(bounds(2,1),bounds(2,2)*1.01,num_bins(2)+1);
    y_edges = linspace(bounds(3,1),bounds(3,2)*1.01,num_bins(3)+1);
    [voxel_counts,bin_edges,bin_cents,~] = histcn(txy,t_edges,x_edges,y_edges);
    bin_volumes = (OuterProduct(diff(t_edges)'.*diff(x_edges),diff(y_edges')));
    voxel_counts = voxel_counts/num_shots;
    bin_flux = voxel_counts./bin_volumes;
    
    h_data.flux3 = bin_flux;
    h_data.counts3 = voxel_counts;
    h_data.volumes = bin_volumes;
    h_data.mean_counts = sum(h_data.counts3,'all');
    
    profile_labels = {'T','X','Y'};    
    axis_labels = {'T (s)','X (m)','Y (m)'};
    h_data.edges = cell(3,1);
    h_data.counts_1d = cell(3,1);   
    h_data.centres = cell(3,1);
    h_data.flux_1d = cell(3,1);
    h_data.num_shots = num_shots;
    if size(txy,1) < min_counts
        warning('Counts too low');
    end
    
    ax_order = [1,3,2]; % to align the X and Y axes in the figure
    
    for axis = 1:3
        ax_excl = 1:3;
        ax_1d = ax_order(axis);
        ax_excl(axis) = [];
        
        h_data.edges{axis} = bin_edges{axis};
        h_data.centres{axis} = bin_cents{axis};
        h_data.counts_1d{axis} = col_vec(squeeze(sum(voxel_counts,ax_excl)));
        h_data.mean_1d = sum(voxel_counts,ax_excl);
        h_data.flux_1d{axis} = h_data.counts_1d{axis}./(col_vec(diff(bin_edges{axis}))); % normalized by 1d bin areas so will give different peak heights
    end
    
    h_data.counts_2d = cell(3,1);
    for txy_count = 1:3
        h_data.counts_2d{txy_count} = squeeze(sum(voxel_counts,txy_count));        
    end
    h_data.txy = txy;
    %% PLOTTING
    grid_height = 2 + sum([log_plot,log_hist]);
    grid_width = 3;
    
    if draw_plots
        stfig(sprintf('BEC TXY display %s',label));
        clf
        tiledlayout(grid_height,grid_width)
        for axis_count = 1:3
            axis = ax_order(axis_count);
%             kept_axes(axis) = [];
            nexttile
            hold on
            plot(h_data.centres{axis},(h_data.flux_1d{axis}))
            title(sprintf('Mean %s flux',profile_labels{axis}))
            ylabel('Flux ')
            xlabel(sprintf('%s',axis_labels{axis}))
            xlim([min(h_data.centres{axis}),max(h_data.centres{axis})])
        end
%         h_data.flux_1d = h_data.flux_1d{1,3,2};
        
        if log_plot
            for axis = 1:3
                nexttile
                title(sprintf('Max %s profile',profile_labels{axis}))            
                hold on
                plot(h_data.centres{axis},h_data.flux_1d{axis})
                ylabel('Peak Flux (Hz m$^{-2}$)')
                xlabel(sprintf('%s',axis_labels{axis}))
                xlim([min(h_data.centres{axis}),max(h_data.centres{axis})])
            end
        end
        
        
        for txy_count = 1:3
            kept_axes = ax_order;
            kept_axes(txy_count) = [];
            edges_2d{1} = h_data.edges{kept_axes(1)};
            edges_2d{2} = h_data.edges{kept_axes(2)};
            nexttile
            if blur_size == 1
                V = h_data.counts_2d{txy_count};
            else
                V = imgaussfilt(h_data.counts_2d{excl_axis},blur_size);
            end
            imagesc(edges_2d{2},edges_2d{1},V)
            title(sprintf('%s%s projection',profile_labels{kept_axes(1)},profile_labels{kept_axes(2)}));
            set(gca,'YDir','normal')
            xlabel(sprintf('%s',axis_labels{kept_axes(2)}))
            ylabel(sprintf('%s',axis_labels{kept_axes(1)}))
        end
        
        
        if log_hist
            for txy_count = 1:3
                kept_axes = ax_order;
                kept_axes(txy_count) = [];
                edges_2d{1} = h_data.edges{kept_axes(1)};
                edges_2d{2} = h_data.edges{kept_axes(2)};
                nexttile
                if blur_size == 1
                    V = h_data.counts_2d{txy_count};
                else
                    V = imgaussfilt(h_data.counts_2d{excl_axis},blur_size);
                end
                imagesc(edges_2d{2},edges_2d{1},log(V))
                title(sprintf('%s%s projection',profile_labels{kept_axes(1)},profile_labels{kept_axes(2)}));
                set(gca,'YDir','normal')
                xlabel(sprintf('%s',axis_labels{kept_axes(2)}))
                ylabel(sprintf('%s',axis_labels{kept_axes(1)}))
            end
            colormap(viridis)
        end
        
%         if plot_3d > 1
%             stfig('3D BEC display');
%             clf
%             nexttile
%         end
        colormap(viridis)
        suptitle(label)
        
    end
    
end