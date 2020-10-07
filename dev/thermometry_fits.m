function shot_data = thermometry_fits(data,opts)

    cache_opts=[];
    if isfield(opts.temp,'cache_opts'), cache_opts=opts.temp.cache_opts; end
    cache_opts.verbose=1;

    opts.shot_range = nan;
    
    
    opts.shotnums = data.tdc.shot_num;
    if ~isnan(opts.temp.shot_nums)
        shot_nums = opts.temp.shot_nums;
        shot_nums = shot_nums(shot_nums < max(opts.shotnums) & shot_nums >= 1);
        opts.shotnums = shot_nums;        
    end
    cli_header(1,'Fitting thermal profiles');
    shot_data = thermo_fit_core(data,opts);
    

    %%
    % for these temperatures, k_B*T ~ 5e-30, and hbar*omega_trap ~ 1.5e-31
    % so this is right about the crossover

    cli_header(2,'Done.');
end

function shot_data = thermo_fit_core(data,opts)
    mcp_data = data.tdc;
    mcp_data.all_ok = ones(size(mcp_data.N_atoms));
    num_shots = length(opts.shotnums);
    
    opts.c.fall_time = sqrt(2*opts.c.fall_distance/opts.c.g0);
    opts.c.fall_velocity = 9.796*opts.c.fall_time;
    
    % Options for binning and show_TXY
    opts.t0 = 0.4145; % Centre of first pulse (sec)
    opts.pulsedt = .008; %time between pulses (sec)
    opts.num_bins = [50,50,50];
    opts.hist_2d = false;
    opts.draw_plots = false;
    opts.centre_bec = true;
    opts.verbose = false;
    opts.xylim = [-0.03,.03;-.03,.03];
    opts.pulse_twindow = .003; 
    
%     min_flux = 1e4;
%     num_pulses = 75;
    single_shot_plot = opts.single_shot_plot;
    do_refits = true;
    min_shot_counts = 1e2;

    x_refit_lvl = 0.15;
    y_refit_lvl = 0.15;

    const.g0 = 9.796;
    const.mhe = 6.6433e-27;
    const.kb=1.380e-23;
    
    A = 5e3;
    T_guess=250; %nK
    x_guess.bose = [A,0,T_guess,0];
    x_guess.gauss = x_guess.bose;

%     gfit_guess(3) = sqrt(const.mhe*gfit_guess(3)/const.kb);
    g_tol = 1e-4;
    n_p = @(p,v) p(1)*g_bose(exp(-const.mhe*abs((v-p(2))).^2/(2*const.kb*p(3)*1e-9)),g_tol)+p(4);
    gfun = @(b,x) b(1).*exp(-const.mhe*((x-b(2)).^2)./(2*const.kb*b(3)*1e-9))+b(4); 
    fit_options = statset('TolFun',1e-6);
    
    shot_data=[];
    shot_data.pulse_data = cell(num_shots,1);
    shot_data.ok_mask = zeros(num_shots,1);
    shot_data.shot_num= ones(num_shots,1);
    shot_data.num_counts= ones(num_shots,1);
    shot_data.T = nan(num_shots,2); %val, SE
    
    for idx=1:num_shots
        i = opts.shotnums(idx);
        warning('off')
        if mod(idx,10)==0
            cli_header('Working on shot %u/%u (%u)...',idx,num_shots,i);
        end
        try
            mask = zeros(size(data.tdc.N_atoms));
            mask(i) = 1;
%             this_txy = data.tdc.counts_txy{i};
            tshot = struct_mask(data.tdc,logical(mask));
            topts.draw_plots = false;
            topts.num_bins = [100,100,100];
            topts.verbose = 0;
            pulse_data = show_txy_raw(tshot,topts);
            if tshot.N_atoms < min_shot_counts
                warning('on')
                warning('Insufficient counts in shot %u',i)
                warning('off')
                shot_data.ok_mask(i,:) = 0;
            else
                
%                 pulse_data.x_cens = zeros(opts.num_bins(2));
%                 pulse_data.y_cens= zeros(opts.num_bins(3));
%                 pulse_data.x_flux= zeros(num_pulses,opts.num_bins(2));
%                 pulse_data.y_flux= zeros(num_pulses,opts.num_bins(3));
%                 pulse_data.txy = cell(num_pulses,1);
%                 pulse_data.num= zeros(num_pulses,1);
%                 pulse_data.cen= zeros(num_pulses,3);
%                 pulse_data.std= zeros(num_pulses,3);

%                 for pulse_num = 1:num_pulses
%                     if pulse_num == 80
%                         fprintf('Hi');
%                     end
%                     opts.lims = [opts.t0+(pulse_num-1)*opts.pulsedt,opts.t0+(pulse_num)*opts.pulsedt;
%                                 opts.xylim(1,:);
%                                 opts.xylim(2,:)];
%                     
% %                     h_data = show_txy_raw(tshot,opts);
%                     pulse_data.x_cens = h_data.centres{2};
%                     pulse_data.y_cens = h_data.centres{3};
%                     pulse_data.x_flux(pulse_num,:) = h_data.flux_1d{2};
%                     pulse_data.y_flux(pulse_num,:) = h_data.flux_1d{3};
%                     pulse_data.txy{pulse_num} = h_data.txy;
%                     pulse_data.num(pulse_num) = size(h_data.txy,1);
%                     pulse_data.cen(pulse_num,:) = h_data.pulse_cen;
%                     pulse_data.std(pulse_num,:) = h_data.pulse_std;
%                 end  % loop over pulses

                X = pulse_data.centres{2};
                x_mean_flux = pulse_data.flux_1d{2};
                v_x = X/opts.c.fall_time;
                [max_flux,max_loc] = max(x_mean_flux);
                flux_cen_v = v_x(max_loc);
                x_mask = x_mean_flux < .05*max_flux;
                x_guess.bose = [A,0,T_guess,0];
                x_guess.gauss = x_guess.bose;
                x_guess.bose(2) =  flux_cen_v;
                x_guess.gauss(2) = flux_cen_v;

                gauss_guess_fun = gfun(x_guess.gauss,v_x(x_mask));
                bose_guess_fun = n_p(x_guess.bose,v_x(x_mask));
                
                bose_rescale = mean(bose_guess_fun'./x_mean_flux(x_mask));
                gauss_rescale = mean(gauss_guess_fun'./x_mean_flux(x_mask));
                x_guess.bose(1) =  x_guess.bose(1)/bose_rescale;
                x_guess.gauss(1) = x_guess.gauss(1)/gauss_rescale;
                gauss_guess_fun = gfun(x_guess.gauss,v_x(x_mask));
                bose_guess_fun = n_p(x_guess.bose,v_x(x_mask));
                
                if sum(x_mask) < 10
                    warning('on')
                    warning('Too few points to fix in shot %u X',i)
                    shot_data.ok_mask(i) = 0;
                else
                    if opts.do_fits
                        
%                         x_guess.bose=fit_guess;
                        x_mdl = fitnlm(v_x(x_mask),x_mean_flux(x_mask),n_p,x_guess.bose,...
                                    'Options',fit_options);
                        x_coef.bose= x_mdl.Coefficients.Estimate;
                        x_SE = x_mdl.Coefficients.SE;
                        bose_resid_x =(col_vec(n_p(x_coef.bose,v_x))-col_vec(x_mean_flux));
                        x.bose.rel_err = col_vec(bose_resid_x)./col_vec(x_mean_flux);

                         % Fit Gaussian functions
%                         x_guess.gauss=gfit_guess;
                        
                        x_mdlfit.gauss = fitnlm(v_x(x_mask),x_mean_flux(x_mask),gfun,x_guess.gauss,...
                                    'Options',fit_options);
                        x_coef.gauss = x_mdlfit.gauss.Coefficients.Estimate;
                        gauss_resid_x = (col_vec(gfun(x_coef.gauss,v_x))-col_vec(x_mean_flux));
                        x.gauss.rel_err = col_vec(gauss_resid_x)./col_vec(x_mean_flux); 

    %                     re-fit
%                         x_refit_mask = abs(x.gauss.rel_err) < x_refit_lvl | abs(x.bose.rel_err) < x_refit_lvl;
                        x_refit_mask = abs(x.bose.rel_err) < x_refit_lvl;
                        if sum(x_refit_mask) > sum(x_mask) && do_refits
                            % New things to fit
                            x_mask = x_refit_mask & abs(v_x) > .01;
                            x_mdl = fitnlm(v_x(x_mask),x_mean_flux(x_mask),n_p,x_guess.bose,...
                                        'Options',fit_options);
                            x_coef.bose= x_mdl.Coefficients.Estimate;
                            x_SE = x_mdl.Coefficients.SE;
                            bose_resid_x =(n_p(x_coef.bose,v_x)-x_mean_flux);
                            x.bose.rel_err = bose_resid_x./x_mean_flux;
                        end
                    temp_x =(abs(x_coef.gauss(3)));% *const.mhe/const.kb;
                    temp_x_SE =(abs(x_mdlfit.gauss.Coefficients.SE(3))) *const.mhe/const.kb;
                    shot_data.T(idx,:) = 1e-9*[temp_x,temp_x_SE];
                    end
                    shot_data.ok_mask(i) = 1;
                end



                [~, msgid] = lastwarn; %were there warnings in that attempt?
                if strcmp(msgid,'stats:nlinfit:IllConditionedJacobian')
%                     cli_header(2,'Shot %u Gaussian poorly fitted\n',i);
        %             shot_data.ok_mask(i) = 0;
                    warning('Restting warning...');
                end
                
            shot_data.shot_num(i) = i;
            shot_data.num_counts(i) = mcp_data.N_atoms(i);
%             shot_data.pulse_counts(idx,:) = pulse_data.num;    
            if opts.do_fits
                pulse_data.gauss.residual.x = gauss_resid_x;
                pulse_data.bose.residual.x = bose_resid_x;
                pulse_data.x_mask = x_mask;
            end

                
                if single_shot_plot   
                    cli_header(2,'Plotting...');
                    stfig('Pulse data');
                    clf
                    subplot(2,1,1)
                    hold on
                    plot(v_x,x_mean_flux,'k:')
                    plot(v_x(x_mask),x_mean_flux(x_mask),'k*')
                    plot(v_x(x_mask),gauss_guess_fun,'go')
                    plot(v_x(x_mask),bose_guess_fun,'bo')
%                     plot(v_x(x_mask),x_mean_flux(x_mask),n_p,fit_guess)
                    if opts.do_fits
%                         plot(v_x(x_refit_mask),x_mean_flux(x_refit_mask),'ko')
                        plot(v_x,n_p(x_coef.bose,(v_x)),'k-.')
                        plot(v_x,n_p(x_coef.gauss,(v_x)),'r-.')
                    end
%                     ylim([100,max(2*x_mean_flux)])
                    legend('Flux','Fit domain 1','GFit guess','BFit guess','Polylog','Gaussian')
                    plot(v_x,gfun(x_guess.gauss,v_x),'g:')
                    plot(v_x,n_p(x_guess.bose,v_x),'b:')
                    xlabel('$v_x$')
                    ylabel('Flux')
                    set(gca,'Yscale','log')
                    title('Mean X profile')

                    subplot(2,1,2)
                    hold on
                    if opts.do_fits
                        plot(v_x(x_mask),(n_p(x_coef.bose,abs(v_x(x_mask)))'-x_mean_flux(x_mask))./x_mean_flux(x_mask),'ro')
%                         suptitle(sprintf('Shot %u: T$_P$ =(%.2e,%.2e), T$_G$ =(%.2e,%.2e), N=%u',i,x_coef.bose(3),y_coef.bose(3),gauss_temp_x,gauss_temp_y,tshot.N_atoms));
                    end
                    legend('Polylog','Gaussian','Location','NorthEast')
                    ylim([-1,1])
                    xlabel('$v_x$')
                    ylabel('$\Delta/y$')
                    title('Normalised residuals')
                    
                    
                    drawnow
                    
                end % single shot plot
                
            end %if enough counts

        catch
            warning('on')
            warning('Error in shot %u',i);
            shot_data.ok_mask(i) = 0;
            
        end %try shot
    end %loop over shots
%     if opts.visual
%         stfig('PAL results');
%         clf
%         hold on
%         errorbar(1:length(shot_data.T),1e9*shot_data.T(:,1),1e9*shot_data.T(:,2),'kx')
%         xlabel('Shot number')
%         ylabel('T (nK)')
%         set(gca,'FontSize',16)
%     end
%     title('$v_y$ residuals')
    
end