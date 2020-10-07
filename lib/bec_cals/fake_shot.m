function shot_data = fake_shot(total_number,trap_freqs,const,varargin)

% fake_shot(total_number,trap_freqs,const,varargin) generates data that 
% mimics the k-space density of thermal and quantum depleted fractions, 
% plus a spatially uniform dark count rate. Returns a structure with fields:
%     k_cart        Nx3 k-vector where N = num counts generated
%     k             Lists all k-vector magnitudes (k.detections) and the
%                   maximum index of the thermal, depleted, and background
%                   counts (k.partitions)
%     stats         Various properties of the generated data
% Mandatory inputs:
%     total_number
%     trap_freqs
%     const         struct - table of physical constants (consistent with
%                   hebec_constants)
% Example:
% ```
% const = hebec_constants()
% shot_data = fake_shot(5e5,[420,420,42],const,'temperature',1e-7,'QE',0.1);
% ```

% Parse the inputs
    p = inputParser;
    %     defaultFontSize = 12; 
    %     defaultFont = 'times'; 
    addParameter(p,'temperature',0);
    addParameter(p,'background_rate',0);
    addParameter(p,'QE',1);
    addParameter(p,'contact_scale',1);
    addParameter(p,'k_max',1e7);
    addParameter(p,'k_min',1e6);
    addParameter(p,'verbose',0);
    addParameter(p,'visual',0);
    addParameter(p,'phi_range',nan);
    parse(p,varargin{:});
    temperature = p.Results.temperature;
    background_rate = p.Results.background_rate;
    contact_scale = p.Results.contact_scale;
    QE = p.Results.QE;
    k_max = p.Results.k_max;
    k_min = p.Results.k_min;
    verbose = p.Results.verbose;
    visual= p.Results.visual;
    phi_range= p.Results.phi_range;

% Calculate some useful objects
    w_HO = 2*pi*geomean(trap_freqs);
    total_number = ceil(total_number);
    critical_temperature = (0.94*const.hbar*w_HO*total_number^(1/3))/const.kb;
    if temperature > critical_temperature
        error('T > T_c, many calculations will fail');
    end
    condensed_fraction = 1-(temperature/critical_temperature)^3;
    condensed_number = ceil(condensed_fraction*total_number);
    
    bec_prop = bec_properties(trap_freqs,total_number,'temperature',temperature);
    interparticle_distance = (1/bec_prop.density_peak)^(1/3);
    n0 = bec_prop.density_peak;
    contact = contact_scale*(64*pi^2*const.ahe_scat^2)*n0*condensed_number/7;
    lambda_db = deBroglie_wavelength(temperature);

    % factor of 2 comes from integrating n(k)
    num_depletion_before_cut = ceil(1*(contact/(2*pi^2*k_min)));
    num_depletion_after_cut = 1*(contact/(2*pi^2*k_min) - contact/(2*pi^2*k_max));
    k_thermal = 10/lambda_db;
    num_depletion_outside_thermal = 1*(contact/(2*pi^2*k_thermal) - contact/(2*pi^2*k_max));
    condensed_number = condensed_number - num_depletion_before_cut;
    thermal_number = total_number - condensed_number; % should be non-negative
    detector_volume = (2*k_max)^3;
    n_background = round(background_rate*detector_volume);

    % Generate the samples
    uniform_seed = rand(ceil(num_depletion_before_cut),1);
    alpha = 2; %for the power-law transform - note that this produces the *counts* which are then converted to density
    depletion_kvec = k_min*(1 - uniform_seed).^(-1/(alpha - 1));
    depletion_dirs = randn(ceil(num_depletion_before_cut),3);
    depletion_dirs = depletion_dirs./vecnorm(depletion_dirs')';
    depletion_counts = depletion_kvec.*depletion_dirs;

    sigma_thermal_v = sqrt(const.kb*temperature/const.mhe); 
    thermal_counts = const.mhe*sigma_thermal_v*randn(thermal_number,3)/const.hbar;

    background_counts = 2*k_max*(rand(n_background,3)-0.5); % 3D uniform distribution in a box

    % Account for detector efficiency
    QE_samples = 1:ceil(QE*length(depletion_counts));
    depletion_counts = depletion_counts(QE_samples,:);
    thermal_counts = thermal_counts(1:ceil(QE*length(thermal_counts)),:);
    
    depletion_kvec = depletion_kvec(QE_samples);
    thermal_kvec = vecnorm(thermal_counts,2,2);
    background_kvec = vecnorm(background_counts,2,2);%
    
    % trim to detector
    depletion_mask = (depletion_kvec > k_min & depletion_kvec < k_max);
    thermal_mask = (thermal_kvec > k_min & thermal_kvec < k_max);
    background_mask = (background_kvec > k_min & background_kvec < k_max);
    
            % write some output for feedback
    if verbose
        cli_header('Setting up sims with N=%u, T=%.2f nK w_HO = 2pi*%.1f Hz',total_number,1e9*temperature, w_HO/(2*pi));
        if verbose>1
        cli_header(1,'Wavevector range %.2e, %.2e',k_min,k_max);
        cli_header(1,'Background density %.2e produces %.2e counts', background_rate,n_background);
        cli_header(1,'Cond frac of %.2f produces %.2f counts', condensed_fraction,condensed_number);
        cli_header(1,' With peak density %.2e', n0);
        cli_header(1,' and thermal population is %.2f counts', thermal_number);
        cli_header(1,'The deBroglie wavelength is %.2f micron', 1e6*lambda_db);
        cli_header(1,'Interparticle spacing at peak is %.2e m', interparticle_distance);
        cli_header(1,'The total contact is %.2e (%.1fx Tan theory)',contact,contact_scale);
        cli_header(1,'The Contact region will be from k >> (%.2e,%.2e)',1/interparticle_distance,1/lambda_db);
        cli_header(1,'Av num in %.2e<k will be %.1f (%.2e pc.)',k_min,num_depletion_before_cut,100*num_depletion_before_cut/total_number );
        cli_header(1,'Av num in %.2e<k<%.2e will %.1f (%.2e pc.)',k_min, k_max,num_depletion_after_cut,100*num_depletion_after_cut/total_number );
        cli_header(1,'Av num in %.2e<k<%.2e will %.1f (%.2e pc.)',10/lambda_db, k_max,num_depletion_outside_thermal,100*num_depletion_outside_thermal/total_number );
        end
        cli_header(1,'%u total atoms (%u therm, %u dep, %u backg)',sum([thermal_number,num_depletion_before_cut,n_background]),thermal_number,num_depletion_before_cut,n_background);
        cli_header(1,'Detection efficiency %.2f ',QE);
        cli_header(1,'Total detections (%u therm, %u dep, %u backg)',length(thermal_kvec),length(depletion_counts),length(background_kvec));
        if verbose>1
            cli_header(1,'%u depletion outside thermal',sum(depletion_mask));
        end
    end
    
    if ~isnan(phi_range)
        qd_phi = atan(depletion_counts(:,3)./sqrt(depletion_counts(:,1).^2+depletion_counts(:,2).^2));
        th_phi = atan(thermal_counts(:,3)./sqrt(thermal_counts(:,1).^2+thermal_counts(:,2).^2));
        bg_phi = atan(background_counts(:,3)./sqrt(background_counts(:,1).^2+background_counts(:,2).^2));
        depletion_mask_phi =  abs(qd_phi) > phi_range(1) & abs(qd_phi) < phi_range(2);
        thermal_mask_phi =  abs(th_phi) > phi_range(1) & abs(th_phi) < phi_range(2);
        background_mask_phi =  abs(bg_phi) > phi_range(1) & abs(bg_phi) < phi_range(2);

        depletion_mask = depletion_mask_phi & depletion_mask;
        thermal_mask = thermal_mask_phi & thermal_mask;
        background_mask = background_mask_phi & background_mask;
        if verbose>1
        cli_header('Masking out counts outside %.2fdg<|phi|<%.2fdg, leaving %u/%u counts',rad2deg(phi_range(1)),rad2deg(phi_range(2)),...
                sum(depletion_mask_phi)+sum(thermal_mask_phi)+sum(background_mask_phi),length(depletion_mask)+length(thermal_mask)+length(background_mask));
        cli_header(1,'Final detections (%u therm, %u dep, %u backg)',sum(thermal_mask),sum(depletion_mask),sum(background_mask));
        end
    end
    
    depletion_counts = depletion_counts(depletion_mask,:);
    thermal_counts = thermal_counts(thermal_mask,:);
    background_counts = background_counts(background_mask,:);    

    depletion_kvec = depletion_kvec(depletion_mask);
    thermal_kvec = thermal_kvec(thermal_mask);
    background_kvec = background_kvec(background_mask);   
    
    k_cloud = [thermal_counts;depletion_counts];
    k_cart = [k_cloud;background_counts];    
    k_detections = vecnorm(k_cart,2,2);


        %Make some plots for display

        k_plot = linspace(k_min,k_max,3e2);
        p_plot = const.hb*k_plot;
        v_plot = p_plot/const.mhe;

        g_tol = 1e-5;
        density_scale = (2*pi)^-3;
        % constants out front to convert from p to k
        thermal_BoseEinstein = QE*(const.hbar)^3 * (1/(lambda_db*const.mhe*w_HO)^3)*g_bose(exp(-(p_plot).^2/(2*const.mhe*const.kb*temperature)),g_tol);
        thermal_MaxwellBoltzmann = QE *thermal_number*(const.hbar/const.mhe)^3 * 1/sqrt(2*pi)* (const.mhe/(const.kb*temperature)) * exp(-(1/2)*((p_plot).^2)/(2*const.mhe*const.kb*temperature));
        thermal_chang = QE*thermal_number*g_bose(exp(-((k_plot*lambda_db).^2)/(4*pi)),g_tol)/(1.202*(2*pi/lambda_db)^3);
        
        depletion_profile = QE*density_scale * contact./(k_plot.^4);


       % K cartesian coordinates
        shot_data.k_cart = k_cart;
        % k radial vectors
        shot_data.k.detections = k_detections;
        % Indices to re-partition the counts according to their sources
        shot_data.k.partitions = [size(thermal_counts,1),size(depletion_counts,1),size(background_counts,1)];
        % plotting stuff
%         shot_data.hist.edges = log_edges;
%         shot_data.hist.centres = log_bin_centres;
%         shot_data.hist.volumes = log_bin_volumes;
%         shot_data.hist.counts = hist_allcounts;

        shot_data.stats.N=total_number;
        shot_data.stats.T=temperature;
        shot_data.stats.data_domain = [k_min,k_max];
        shot_data.stats.background_density = background_rate;
        shot_data.stats.background_counts = n_background;
        shot_data.stats.cond_pop = condensed_number;
        shot_data.stats.peak_density = n0;
        shot_data.stats.thermal_pop = thermal_number;
        shot_data.stats.deBroglie_wavelength = lambda_db;
        shot_data.stats.contact = contact;
        shot_data.stats.contact_scale = contact_scale;
        shot_data.stats.k_ranges(1,:) = [k_min,inf,num_depletion_before_cut,100*num_depletion_before_cut/total_number];
        shot_data.stats.k_ranges(2,:) = [k_min,k_max,num_depletion_after_cut,100*num_depletion_after_cut/total_number];
        shot_data.stats.k_ranges(3,:) = [10/lambda_db,k_max,num_depletion_outside_thermal,100*num_depletion_outside_thermal/total_number];
        shot_data.stats.QE = QE;
        shot_data.stats.total_atoms = [thermal_number,num_depletion_before_cut,n_background;
                                      length(thermal_kvec),length(depletion_counts),length(background_kvec)];
        shot_data.stats.QD_detect_region = sum(depletion_mask);
        shot_data.stats.noise_detect_region = sum(background_mask);
        shot_data.properties = bec_prop;
        shot_data.detection_density = nan;
        
        if visual
%             xplt = linspace(1,100,101);
            log_edges = logspace(min(log10(k_detections)),max(log10(k_detections)),100);
            hist_signal = histcounts(depletion_kvec,log_edges);
            hist_noise = histcounts(background_kvec,log_edges);
            hist_thermal = histcounts(thermal_kvec,log_edges);
            [hist_allcounts,hcs] = histcounts(k_detections,log_edges);
            % these are in kspace so let's integrate in spherical coordinates
            log_bin_volumes= (4*pi/3)*(log_edges(2:end).^3-log_edges(1:end-1).^3);
            log_bin_centres= 0.5*(log_edges(2:end)+log_edges(1:end-1));
            noise_density = (hist_noise./log_bin_volumes);
            detection_density = hist_allcounts./log_bin_volumes;
%             shot_data.hist.density = detection_density;
%             shot_data.hist.edges = log_edges;
            plot_profile = thermal_chang + depletion_profile + background_rate;

            stfig('Data generation');
            clf

%             subplot(2,2,1)
%             title('Detection counts')
%             hold on
%             plot((log_bin_centres),(hist_signal),'--')
%             plot((log_bin_centres),(hist_noise),':')
%             plot((log_bin_centres),(hist_thermal),'-.')
%             plot((log_bin_centres),(hist_allcounts))
% 
%             legend('Depletion','Noise floor','Thermal','Combined','Location','Best')
%             set(gca,'Xscale','log')
%             set(gca,'Yscale','log')
%             box off
            subplot(2,1,1)
            title('Simulated detection densities')
            hold on
            plot((log_bin_centres),(hist_thermal./log_bin_volumes),'-.','LineWidth',2)
            plot((log_bin_centres),(hist_signal./log_bin_volumes),'-.','LineWidth',2)
            plot((log_bin_centres),noise_density,'-.','LineWidth',2)
            plot((log_bin_centres),(detection_density),'LineWidth',2)
            plot(k_plot,depletion_profile,':','LineWidth',2)
%             plot(k_plot,thermal_MaxwellBoltzmann,':','LineWidth',2)
%             plot(k_plot,thermal_BoseEinstein,':','LineWidth',2)
            plot(k_plot,thermal_chang,':','LineWidth',2) 
            plot(k_plot,plot_profile,'LineWidth',2)
            legend('Thermal','Depletion','Noise floor','All counts','QD thry','Thermal thry','Full thry',...
                'Location','Best')
            ylim([.3*min(plot_profile),10*max(plot_profile)])
            xlabel('k ($ m^{-1}$)')
            ylabel('n(k) ($m^{3}$)')
            set(gca,'Xscale','log')
            set(gca,'Yscale','log')
            set(gca,'FontSize',16)
%             set(gca,'ColorOrder',cmap(randsample(1:10,10),:));
            box off

%             subplot(2,2,3)
%             hold on
%             plot((log_bin_centres),(hist_thermal./log_bin_volumes),'-.')
%             plot((log_bin_centres),(hist_signal./log_bin_volumes),'-.')
%             plot((log_bin_centres),noise_density,'-.')
%             plot((log_bin_centres),(detection_density))
%             plot(k_plot,depletion_profile,':')
%             plot(k_plot,thermal_MaxwellBoltzmann,':')
%             plot(k_plot,thermal_BoseEinstein,':')
%             plot(k_plot,plot_profile)
% %             legend('Thermal','Depletion','Noise floor','All counts','QD thry','MB thry','BE thry','Full thry',...
% %                 'Location','Best')
%             ylim([.3e-20,1.5*max(plot_profile)])
%             xlabel('k ($ m^{-1}$)')
%             ylabel('n(k) ($ m^{3}$)')
%             set(gca,'FontSize',16)
% %             set(gca,'ColorOrder',cmap(randsample(1:10,10),:));
%             box off

            subplot(2,1,2)
            hold on
            plot((log_bin_centres),log_bin_centres.^4.*(hist_thermal./log_bin_volumes),'-.','LineWidth',2)
            plot((log_bin_centres),log_bin_centres.^4.*(hist_signal./log_bin_volumes),'-.','LineWidth',2)
            plot((log_bin_centres),log_bin_centres.^4.*noise_density,'-.','LineWidth',2)
            plot((log_bin_centres),log_bin_centres.^4.*(detection_density),'LineWidth',2)
            plot(k_plot,k_plot.^4.*depletion_profile,':','LineWidth',2)
%             plot(k_plot,k_plot.^4.*thermal_MaxwellBoltzmann,':')
            plot(k_plot,k_plot.^4.*thermal_chang,':','LineWidth',2)
            plot(k_plot,k_plot.^4.*plot_profile,'LineWidth',2)
%             legend('Thermal','Depletion','Noise floor','All counts','QD thry','MB thry','BE thry','Full thry',...
%                 'Location','Best')
            ylim([.5*min(k_plot.^4.*plot_profile),10*max(k_plot.^4.*plot_profile)])
            xlabel('k ($ m^{-1}$)')
            ylabel('n(k) ($m^{3}$)')
            set(gca,'Xscale','log')
            set(gca,'Yscale','log')
            set(gca,'FontSize',16)
            

%             subplot(2,2,3)
%             title('Scaled densities')
%             hold on
%             plot((log_bin_centres),(hist_signal./log_bin_volumes).*log_bin_centres.^4,'--')
%             plot((log_bin_centres),noise_density.*log_bin_centres.^4,'--')
%             plot((log_bin_centres),(hist_thermal./log_bin_volumes).*log_bin_centres.^4,'--')
%             plot((log_bin_centres),(hist_allcounts./log_bin_volumes).*log_bin_centres.^4)
%             plot(k_plot,contact*ones(size(k_plot))/(2*pi)^3)
%             legend('Depletion','Background','Thermal','Detections','Theory',...
%                 'Location','Best')
%             set(gca,'Xscale','log')
%             set(gca,'Yscale','log')
%             box off
% 
%             subplot(2,2,4)
%             hold on
%             plot(k_plot,intCDF(k_plot,k_detections))
%             plot(k_plot,intCDF(k_plot,thermal_kvec))
%             plot(k_plot,intCDF(k_plot,background_kvec))
%             plot(k_plot,(intCDF(k_plot,k_detections)-intCDF(k_plot,background_kvec)))
%             plot(k_plot,(intCDF(k_plot,depletion_kvec)))
% 
%             legend('Total','Thermal','Background','Total - Background','Depletion','Location','Best')
%             ylabel('N(k<K)')
%             xlabel('K')
            drawnow
        end
%     end
end