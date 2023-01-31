%% Initializing path
clear all;
tic
% close all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder));
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')
hebec_constants
combined_struct = @(S,T) cell2struct(cellfun(@vertcat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);

%% Import directory
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
%  opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
% opts.data_root = 'C:\Users\BEC Machine\Documents\DATA_BACKUP\';

% data_dir = '20220406_TF_vs_evap';
data_dir = '20220729_he3_4_vs_hold_time';%'20220819_he3_he4_data_4\20220819_evap_865_detune_fracN_368.1';%
% data_dir = '20220706_large_he3_clouds_detuning';
% data_dir = '20220328_TF_dep_trial_5_fracNwizard_detuning';
% data_folder = {'20220324_he3_and_4_seperated_clouds'};
% Find all the data folders in the data directory

f_winfreak = 10767.5;%in MHz
% expected detuning 33.574 GHz

full_data_dir = fullfile(opts.data_root, data_dir);
% check if there's a log file
files = dir(full_data_dir);
name_file = {files.name};
log_mask = cellfun(@(x) contains(x,'log_'),name_file);
log_check = any(log_mask);

if ~log_check
    if exist('data_dir','var')
        f_indx = 1;
        dir_folders = dir(full_data_dir);
        for current_folder = {dir_folders.name}
            folder_str=(current_folder{1});
            loc = strfind(current_folder,'evap_');
            loc2 = strfind(current_folder,'_detune');
            loc3 = strfind(current_folder,'_fracN');
            loc = loc{1};
            loc2 = loc2{1};
            loc3 = loc3{1};
            if isempty(loc3)
                offset = 0;
            else
                offset = 6;
            end
            if isempty(loc) || isempty(loc2)
                continue
            end
            evap_indx = (loc+5):(loc2-1);
            detuning_indx = (loc2+8+offset):length(folder_str);
            evap(f_indx) = str2num(folder_str(evap_indx));
            detuning(f_indx) = -3.*f_winfreak-str2num(folder_str(detuning_indx))+33574;
            data_folder{f_indx} = fullfile(data_dir,current_folder{1});
            f_indx = f_indx +1;
        end
    end

    % check if detuning is in directory title
    % detuning_cehck
    detuning_cehck = contains(data_dir,'detuning');
    if detuning_cehck
        variable = detuning;
        var_str = '$f_{He^3}-f_{He^4}+33,574$ MHz';%detuning';
    else
        variable = evap-825;
        var_str = 'Evaporator height (kHz)';
    end
    dirs_list = fullfile(opts.data_root, data_folder);
    if iscell(dirs_list)
        num_dirs = length(dirs_list);
    else
        num_dirs = 1;
    end
    % data_folder = {'20220324_he3_and_4_seperated_clouds'};
else
    log_file = name_file(log_mask);
    log_file = log_file{1};
    log_table = readtable(fullfile(full_data_dir,log_file));
    variable_vec = log_table{:,4};
    variable = unique(variable_vec);
    var_str = 'shot type';%'hold time';
    num_dirs = size(variable);
    indx_list = 1:size(variable_vec,1);
    dirs_list = variable;
    if isstr(variable(1))
    for ii = 1:num_dirs
        shot_list{ii} = indx_list(strcmp(variable_vec,'mix_shot'));%variable{ii}
    end
    else
    for ii = 1:num_dirs
        shot_list{ii} = indx_list(variable_vec==variable(ii));
    end
    end
end

x = linspace(0.001,200,1e3);
y = real(1./(-polylog(3,-x).*6).^(1/3));

opts.import.force_reimport = true;
opts.import.force_cache_load = ~opts.import.force_reimport;
% opts.import.shot_num = 26:50;
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
%% Analysis Params

% for the hold time data
% he4_time = [0.4392,0.46;0.4387,0.46;0.4378,0.46;0.437,0.46;0.435,0.46;0.432,0.46];%[0.435,0.46];%[1.05,1.09];
% he3_time = [0.41,0.4392;0.41,0.4387;0.41,0.4378;0.41,0.437;0.41,0.435;0.41,0.432];%[0.41,0.435];%[1.01,1.05];

he4_time = [0.437,0.46];
he3_time = [0.41,0.437];

he4_cen = [-3.7796,-4.16263].*1e-3;
he3_cen = [-3.827,-6.5065].*1e-3;

he3_lim = [5,1000];%[500,1000];

min_num = 1e2;

ax = 2; %which axis do we analyse [t,x,y or radial]

axis_labe = {'t','x','y','r'};

do_bimod = 0;
plot_all = 0;
show_fit_goodness = 0;

plot_model_comparison = 0;

%% separate He 3 and He 4 distrabutions

N_he4_avg = nan.*ones(size(dirs_list,2),1);
N_he3_avg = nan.*ones(size(dirs_list,2),1);
if ax == 1
    num_vars = 3;
else
    num_vars = 4;
end
he4_fits = nan.*ones(size(dirs_list,2),num_vars);
he4_fits_bimodal = nan.*ones(size(dirs_list,2),num_vars+2);
he3_fits = nan.*ones(size(dirs_list,2),num_vars);
TF = nan.*ones(size(dirs_list,2),1);
TC = nan.*ones(size(dirs_list,2),1);
cahce_exists = isfile(fullfile(full_data_dir,['clean_data_',num2str(ax),'.mat']));
if isstr(variable(1))
    num_dirs = 1;
end
if ~cahce_exists
    for jj = 1:num_dirs
        %% Import data

        % import using different directories
        if ~log_check
            if iscell(dirs_list)
                opts.import.dir = dirs_list{jj};
                opts.import.cache_save_dir = fullfile(dirs_list{jj}, 'cache', 'import\');
                [data, ~] = import_mcp_tdc_data(opts.import);
            else
                opts.import.dir = fullfile(opts.data_root, data_folder);
                [data, ~] = import_mcp_tdc_data(opts.import);
            end
            % import using log file
        else
            opts.import.dir = full_data_dir;
            opts.import.shot_num = shot_list{jj};
            opts.import.force_reimport = true;
            opts.import.force_cache_load = ~opts.import.force_reimport;
            [data, ~] = import_mcp_tdc_data(opts.import);
            data=struct_mask(data,~cellfun('isempty',data.counts_txy));
        end
        %% remove any hotspts

        data_ht_spot=hotspot_mask(data);
        data.counts_txy=data_ht_spot.masked.counts_txy;
        data.num_counts=data_ht_spot.masked.num_counts;

        clear flux_he4 flux_he3 N_he3 N_he4 he4_txy_c he3_txy_c
        kk = 1;
        for ii = 1:length(data.counts_txy)
            lims_4= [he4_time(jj,:); -0.03, 0.03; -0.03, 0.03];
            he4_txy = masktxy_square(data.counts_txy{ii}, lims_4);
            lims_3 = [he3_time(jj,:); -0.03, 0.03; -0.03, 0.03];
            he3_txy = masktxy_square(data.counts_txy{ii}, lims_3);

            t_he4 = linspace(he4_time(jj,1),he4_time(jj,2),5e3).';
            t_he3 = linspace(he3_time(jj,1),he3_time(jj,2),5e3).';

            space_bins = linspace(-30e-3,30e-3,5e3).';
            rad_bins = linspace(0.1e-3,30e-3,5e3).';

            if ax == 1
                bin_cen_4 = t_he4;
                bin_cen_3 = t_he3;
            elseif ax ~= 4
                bin_cen_4 = space_bins;
                bin_cen_3 = space_bins;
            else
                bin_cen_4 = rad_bins;
                bin_cen_3 = rad_bins;
            end

            %% run good shot checks
            is_shot_good = (size(he4_txy,1)+size(he3_txy,1))>min_num;
            he3_num_check = size(he3_txy(:,1),1)>he3_lim(1) && size(he3_txy(:,1),1)<he3_lim(2);
            is_shot_good= is_shot_good&he3_num_check;
            shot_check(ii) = is_shot_good;
            if is_shot_good
                %% histogram in time
                sigma = 0.5e-4;
                if ax == 4
                    count_hist_he4 = smooth_hist(sqrt((he4_txy(:,2)-he4_cen(1)).^2+(he4_txy(:,3)-he4_cen(2)).^2),'sigma',sigma,'edges',bin_cen_4);
                    count_hist_he4.count_rate.smooth = count_hist_he4.count_rate.smooth./(2.*pi*count_hist_he4.bin.centers);
                    count_hist_he3 = smooth_hist(sqrt((he3_txy(:,2)-he3_cen(1)).^2+(he3_txy(:,3)-he3_cen(2)).^2),'sigma',sigma,'edges',bin_cen_3);
                    count_hist_he3.count_rate.smooth = count_hist_he3.count_rate.smooth./(2.*pi*count_hist_he3.bin.centers);
                else
                    count_hist_he4 = smooth_hist(he4_txy(:,ax),'sigma',sigma,'edges',bin_cen_4);
                    count_hist_he3 = smooth_hist(he3_txy(:,ax),'sigma',sigma,'edges',bin_cen_3);
                    
                    %take average over radial direction
                    sigmar = 0.8e-4;
                    count_hist_he4_r = smooth_hist(sqrt((he4_txy(:,2)-he4_cen(1)).^2+(he4_txy(:,3)-he4_cen(2)).^2),'sigma',sigmar,'edges',rad_bins);
                    count_hist_he3_r = smooth_hist(sqrt((he3_txy(:,2)-he3_cen(1)).^2+(he3_txy(:,3)-he3_cen(2)).^2),'sigma',sigmar,'edges',rad_bins);
                    
                    bin_centres_r_he4 = count_hist_he4_r.bin.centers;
                    bin_centres_r_he3 = count_hist_he3_r.bin.centers;
                    dbin_r = bin_centres_r_he4(2)-bin_centres_r_he4(1);

                    count_hist_he4_r.count_rate.smooth = count_hist_he4_r.counts.smooth./(2.*pi*(bin_centres_r_he4*dbin_r+dbin_r^2/2));
                    count_hist_he3_r.count_rate.smooth = count_hist_he3_r.counts.smooth./(2.*pi*(bin_centres_r_he3*dbin_r+dbin_r^2/2));

                    
                    
                    flux_r_he4{kk} = count_hist_he4_r.count_rate.smooth;
                    flux_r_he3{kk} = count_hist_he3_r.count_rate.smooth;
                    

                end
                flux_he4{kk} = count_hist_he4.count_rate.smooth;
                flux_he3{kk} = count_hist_he3.count_rate.smooth;
                bin_centres_he4 = count_hist_he4.bin.centers;
                bin_centres_he3 = count_hist_he3.bin.centers;
                %record atom number
                N_he4(kk) = size(he4_txy(:,1),1);
                N_he3(kk) = size(he3_txy(:,1),1);

                he4_txy_c{kk} = he4_txy;
                he3_txy_c{kk} = he3_txy;

                kk = kk+1;
            end
        end
        flux_he3_all{jj} = flux_he3; %flux per shot per detuning
        flux_he4_all{jj} = flux_he4;
        N_he3_all{jj} = N_he3; %atom numbers per detuning
        N_he4_all{jj} = N_he4;
        he4_txy_all{jj} = he4_txy_c;
        he3_txy_all{jj} = he3_txy_c;
        if ax ~= 4
            flux_he3_r_all{jj} = flux_r_he3; %flux per shot per detuning
            flux_he4_r_all{jj} = flux_r_he4;
        end

    end

    %% cleaning shots & saving
    [he3_flux_clean, he4_flux_clean,N_he3_clean,N_he4_clean] = clean_data(flux_he3_all,flux_he4_all,bin_centres_he4,bin_centres_he3,N_he3_all ,N_he4_all,ax);
    save(fullfile(full_data_dir,['clean_data_',num2str(ax),'.mat']),"he4_flux_clean","he3_flux_clean","N_he4_clean","N_he3_clean","N_he3_all","N_he4_all","flux_he3_all","flux_he4_all",'bin_centres_he4',"bin_centres_he3","he4_txy_all","he3_txy_all","flux_he3_r_all","flux_he4_r_all","bin_centres_r_he3","bin_centres_r_he4")
else
    load(fullfile(full_data_dir,['clean_data_',num2str(ax),'.mat']))
end


%% Analysing data
for kk = 1: num_dirs
    he4_flux_avg = mean(cell2mat(flux_he4_all{kk}),2); %average of all shots for current detuning
    he3_flux_avg = mean(cell2mat(flux_he3_all{kk}),2);

    he4_flux_avg_r = mean(cell2mat(flux_he4_r_all{kk}),2); %average of all shots for current detuning
    he3_flux_avg_r = mean(cell2mat(flux_he3_r_all{kk}),2);

    he4_flux_avg_all(kk,:)= mean(cell2mat(flux_he4_all{kk}),2); % average of all shots for all detuning
    he3_flux_avg_all(kk,:) = mean(cell2mat(flux_he3_all{kk}),2);
    N_he3_avg(kk,1) = nanmean(N_he3_clean{kk});
    N_he4_avg(kk,1) = nanmean(N_he4_clean{kk});
    N_he4_unc(kk,1) = nanstd(N_he4_clean{kk})./sqrt(length(N_he4_clean{kk}));
    N_he3_unc(kk,1) = nanstd(N_he3_clean{kk})./sqrt(length(N_he3_clean{kk}));
    hebec_constants

    omega = 2.*pi.*[605 60 600.50]; %z,x,y
    omega_unc = [1 5 4e-3];
    omega_bar_he4 = geomean(omega);%average trap freq
    omega_bar_he3 = sqrt(4/3).*omega_bar_he4;

    %     R_TF = ;% expected thomas fermi radius


    %% Fit clouds

    %He4 fit

    %     threshold = (100:5:400).*1e3;
    if ax == 1
        threshold = 0.005;%linspace(0.001,0.006,150);%0.008
        initials = [1e-7 0.65 1e3];
        initials_he3 = [1e-7 0.627 max(he3_flux_avg)./1e2];
    else
        threshold = linspace(0.004,0.014,14);%linspace(0.008,0.012,12);%x-axis limitslinspace(10e-3,15e-3,6);%y-axis lims  %linspace(0.008,0.01557,12);% 0.008;%10.8e-3;%linspace(0.005,0.01557,12);%linspace(0.0044,0.0144,8);%linspace(0.001,0.007,6);%0.011;%0.0022;%0.003;%%linspace(0.0025,0.012,15);%linspace(0.003,0.010,8);%
        initials = [1e-7 0.65 1e3 -0.012];
        initials_he3 = [1e-7 0.65 1e3 -0.012];%[1e-7 0.627 max(he3_flux_avg)./1e2 -0.012];
    end
    [he4_fit_avg,new_he4_flux,mask_he4,~,he4_fit_unc, he4_fit_std] = fit_model(bin_centres_he4,he4_flux_avg,initials,threshold,N_he4_avg(kk),'he4','gauss',ax);
    if do_bimod
        [he4_fit_avg_bimodal,new_he4_flux_bimodal,mask_he4_bimodal,] = fit_model(bin_centres_he4,he4_flux_avg,initials,threshold,N_he4_avg(kk),'he4','bimodal',ax);
        new_he4_flux_bimodal_all(kk,:) = new_he4_flux_bimodal;
        he4_fits_bimodal(kk,:) = he4_fit_avg_bimodal;
    end
    clear initials
    %He3 fit
    threshold_he3 = -0.001;
    [he3_fit,new_he3_flux,mask_c_he3,~,he3_fit_unc,he3_fit_std] = fit_model(bin_centres_he3,he3_flux_avg,initials_he3,threshold_he3,N_he3_avg(kk),'he3','gauss',ax);

    if do_bimod
        k=const.kb;
        qe = 0.1;%0.08;
        TF_c = const.hb*omega_bar_he3/k.*(N_he3_avg(kk,1)./qe*6).^(1/3);%fermi temp
        T_on_TF = he4_fit_avg(1)/TF_c;
        [val,indx] = min(abs(y-T_on_TF));
        muf(kk) = x(indx);
        he3_thr_r = 0.2e-3;%
        % fit with fixed temperature
        initials_he3 = [he4_fit_avg(1),-3e-3 1e3,1];
%         [he3_fit_bimodal_T,T_fix_fit,~,he3_fit_bimodal_T_unc,he3_fit_bimodal_T_std] = fit_model(bin_centres_r_he3,he3_flux_avg_r,initials_he3,he3_thr_r,N_he3_avg(kk),'he3','bimodal_T_fix',4);
        [he3_fit_bimodal_T,T_fix_fit,~,he3_fit_bimodal_T_unc,he3_fit_bimodal_T_std] = fit_model(bin_centres_he3,he3_flux_avg,initials_he3,threshold_he3,N_he3_avg(kk),'he3','bimodal_T_fix',ax);

        % fix with fit muf
        initials_he3 = [he3_fit(1),-3e-3,1e3,muf(kk)];
%         [he3_fit_bimodal_mu,~,mask_he3_r,he3_fit_bimodal_mu_unc,he3_fit_bimodal_mu_std] = fit_model(bin_centres_r_he3,he3_flux_avg_r,initials_he3,he3_thr_r,N_he3_avg(kk),'he3','bimodal_mu_fix',4);

        [he3_fit_bimodal_mu,~,mask_he3_r,he3_fit_bimodal_mu_unc,he3_fit_bimodal_mu_std] = fit_model(bin_centres_he3,he3_flux_avg,initials_he3,threshold_he3,N_he3_avg(kk),'he3','bimodal_mu_fix',ax);
 
        he3_fits_T_fix(kk,:) = he3_fit_bimodal_T;
        he3_fits_mu_fix(kk,:) = he3_fit_bimodal_mu;
    end


    he4_fits(kk,:) = he4_fit_avg;
    he3_fits(kk,:) = he3_fit;


    %         he3_fits_full(kk,:) = he3_fit_bimodal;
    new_he3_flux_all(kk,:) = new_he3_flux;
    new_he4_flux_all(kk,:) = new_he4_flux;

    he4_masks(kk,:) = mask_he4;


    if plot_all
        stfig(['flux fit ',var_str,' = ',num2str(variable(kk))]);
        clf
        hold on
        plot(bin_centres_he4(mask_he4),he4_flux_avg(mask_he4),'kx',...
            bin_centres_he4(~mask_he4),he4_flux_avg(~mask_he4),'rx',...
            bin_centres_he4,new_he4_flux)
        %     plot(bin_centres_he4,he4_flux_avg,'kx',...
        %         bin_centres_he4,new_he4_flux)
        plot(bin_centres_he3,he3_flux_avg,'bx',...
            bin_centres_he3,new_he3_flux)
        xlabel('time')
        ylabel('flux')
        ylim([0 max(max(he4_flux_avg),max(he3_flux_avg)).*1.2])
        if do_bimod
            stfig(['bimodal flux fit ',var_str,' = ',num2str(variable(kk))]);
            clf
            hold on
            plot(bin_centres_he4,he4_flux_avg,'kx',...
                bin_centres_he4,new_he4_flux_bimodal)
            %         plot(bin_centres_he4,he4_flux_avg,'kx',...
            %             bin_centres_he4,new_he4_flux)
            %     plot(bin_centres_he3,he3_flux_avg,'bx',...
            %         bin_centres_he3,new_he3_flux_bimodal)
            xlabel('time')
            ylabel('flux')
        end
    end

    %% Temperature & fits

    T_he4(kk,1) = he4_fit_avg(1,1);
    T_he3(kk,1) = he3_fit(1,1);

    T_he4_unc(kk,1) = sqrt(he4_fit_unc(1,1).^2+he4_fit_std(1,1).^2+(10e-9)^2);
    T_he3_unc(kk,1) = sqrt(he3_fit_unc(1,1).^2+he3_fit_std(1,1).^2+(10e-9)^2);
    if do_bimod
        T_he3_mu_unc(kk,1) = sqrt(he3_fit_bimodal_mu_unc(1,1).^2+he3_fit_bimodal_mu_std(1,1).^2);
        mu_he3_unc(kk,1) = sqrt(he3_fit_bimodal_T_unc(1,3).^2+he3_fit_bimodal_T_std(1,3).^2);
    end

    k=const.kb;
    qe = 0.08;%0.08;
    TF(kk,1) = const.hb*omega_bar_he3/k.*(N_he3_avg(kk,1)./qe*6).^(1/3);%fermi temp
    TC(kk,1) = const.hb*omega_bar_he4/k.*(N_he4_avg(kk,1)./qe./zeta(3)).^(1/3);%condensate temp

    tf = sqrt(2.*const.fall_distance./const.g0);
    % position space radii
    aho = sqrt(const.hb./(const.mhe.*omega));
    R_TF = aho.*(15.*N_he4_avg.*const.ahe_scat./aho).^(1/5);
    % momentum space values
    aho_k = 1./aho;
    R_TF_k = 3.83./R_TF;
    % time-of flight values
    c_tof = [1/const.g0,tf,tf].*const.hb./const.mhe;
    aho_tof = aho_k.*c_tof;
    R_TF_tof = R_TF_k.*c_tof;

    fit_he4_NMSE(kk,1) = goodnessOfFit(new_he4_flux(mask_he4),he4_flux_avg(mask_he4),'NMSE');
    fit_he3_NMSE(kk,1) = goodnessOfFit(new_he3_flux,he3_flux_avg,'NMSE');
    fit_he4_NRMSE(kk,1) = goodnessOfFit(new_he4_flux(mask_he4),he4_flux_avg(mask_he4),'NRMSE');
    fit_he3_NRMSE(kk,1) = goodnessOfFit(new_he3_flux,he3_flux_avg,'NRMSE');

end

%% calculate T/TC from bimodal fit
m4 = const.mhe;
for kk = 1:num_dirs
    if plot_model_comparison && do_bimod
        figure(kk+100)
        clf
        plot(bin_centres_he4,thomas_fermi_dist(he4_fits_bimodal(kk,:).',bin_centres_he4,ax))
        hold on
        plot(bin_centres_he4,thermal_dist(he4_fits_bimodal(kk,:).',bin_centres_he4,m4,ax))
        xlabel('bin center')
        ylabel('flux')
    end

    Nt(kk) = trapz(bin_centres_he4,thermal_dist(he4_fits(kk,:).',bin_centres_he4,m4,ax));
    if do_bimod
        Nc(kk) = trapz(bin_centres_he4,thomas_fermi_dist(he4_fits_bimodal(kk,:).',bin_centres_he4,ax));
        T_on_Tc(kk) = (1-Nc(kk)./(Nc(kk)+Nt(kk))).^(1/3);
    end
end


%% Plots

l=const.fall_distance;%0.847;
m4 = const.mhe;
g=const.g0;
k= const.kb;
hb = const.hb;
m3 = 5.008237293822000e-27;
qe = 0.08;%0.08;
t0 = sqrt(2*l/g);

stfig(['number vs ',var_str]);
clf
box on
subplot(2,1,1)
title('raw number')
errorbar(variable,N_he3_avg./qe,N_he3_unc./qe,'o','CapSize',0,'MarkerSize',10,'LineWidth',3.5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:))
hold on
errorbar(variable,N_he4_avg./qe,N_he4_unc./qe,'o','CapSize',0,'MarkerSize',10,'LineWidth',3.5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_main(2,:))
xlabel(var_str)
ylabel('Atom Number')
legend('He$^3$','He$^4$')
set(gca,'FontSize',30)
subplot(2,1,2)
title('Ratio of He$^3$ to He$^4$');
errorbar(variable,N_he3_avg./N_he4_avg,N_he3_avg./N_he4_avg.*sqrt((N_he4_unc./N_he4_avg).^2+(N_he3_unc./N_he3_avg)),'o','CapSize',0,'MarkerSize',10,'LineWidth',2.5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:))
xlabel(var_str)
ylabel('Num He$^3$/Num He$^4$')
set(gca,'FontSize',30)

stfig('he3-he4 temp 3');
% clf
hold on
box on
errorbar(variable,T_he3.*1e9,T_he3_unc.*1e9,'o','CapSize',0,'MarkerSize',10,'LineWidth',3.5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:))
errorbar(variable,T_he4.*1e9,T_he4_unc.*1e9,'o','CapSize',0,'MarkerSize',10,'LineWidth',3.5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_main(2,:))
xlabel(var_str)
ylabel('Fitted Temperature (nK)')
legend('He$^3$','He$^4$')
set(gca,'FontSize',30)

stfig('T/TF and T/TC');
clf
hold on
box on
% scatter(detuning,T_he3./TF,'x')
errorbar(variable,T_he4./TF,T_he4_unc./TF,'o','CapSize',0,'MarkerSize',10,'LineWidth',3.5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:))
errorbar(variable,T_he4./TC,T_he4_unc./TC,'o','CapSize',0,'MarkerSize',10,'LineWidth',3.5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_main(2,:))
xlabel(var_str)
if do_bimod
    scatter(variable,T_on_Tc,17,'o')
end
xlabel(var_str)
ylabel('Critical Temperature Ratio')
legend('T/T$_F$','T/T$_C$')
set(gca,'FontSize',30)


%% Goodness of fit
if show_fit_goodness
    stfig('he3 fit algo');
    clf
    hold on
    scatter(variable,fit_he3_NMSE)
    scatter(variable,fit_he3_NRMSE)
    xlabel(var_str)
    ylabel('goodness of fit')
    legend('NMSE','NRMSE')

    stfig('he4 fit algo');
    clf
    hold on
    scatter(variable,fit_he4_NMSE)
    scatter(variable,fit_he4_NRMSE)
    xlabel(var_str)
    ylabel('goodness of fit')
    legend('NMSE','NRMSE')
end

%% Plot of average flux per detuning
plt_indx = 4;%1;%
shots_indx = 1:49;%1:35;%
if ax == 1
    initials = [1e-7 0.04 1e5];
    initials_he3 = [1e-7 0.42 5e2];%[1e-7 0.005 1e5];
else
    initials = [1e-7 -8e-3 1e3 -0.012];
    initials_he3 = initials;
end
% initials_he3 = [2e-7,-3e-3,1e2,muf(plt_indx)];
% initials_he3 = [2.23e-7,-8e-3,1e3,0.6];
% figure
% plot(polylog(2.5,x))
variable(plt_indx)
he4_fit_plt = he4_fits(plt_indx,:);%
he3_fit_plt = he3_fits(plt_indx,:);%flux_he4_all{jj}
new_he3_flux = mean(cell2mat(flux_he3_all{plt_indx}(shots_indx)),2).';%new_he3_flux_all(plt_indx,:).';
new_he4_flux = mean(cell2mat(flux_he4_all{plt_indx}(shots_indx)),2).';%new_he4_flux_all(plt_indx,:).';
[val,ind] = max(new_he4_flux);
mask_he4 = ~(abs(bin_centres_he4-bin_centres_he4(ind))<0.0023);%logical(he4_masks(plt_indx,:)).';
threshold_plt = -17e-3;%16.8e-3;%0.0015;%-0.022;%0.0022;%0.004;
threshold_plt_he3 = 11e-3;%0.0015;%-0.022;%0.0022;%0.004;
[he4_fit_avg,fit_flux,mask_he4,~,he4_fit_unc, he4_fit_std] = fit_model(bin_centres_he4,new_he4_flux.',initials,threshold_plt,N_he4_avg(plt_indx),'he4','gauss',ax);
[he3_fit_avg,fit_flux_3,mask_he3,temp,he3_fit_unc, he3_fit_std] = fit_model(bin_centres_he3,new_he3_flux.',initials_he3,threshold_plt_he3,N_he3_avg(plt_indx),'he3','gauss',ax);


stfig('flux fit per indx');
clf
hold on
box on
scale_factor = 1;
if ax ~= 1
    scale_factor = 1e3;
    subplot(2,1,1)

    % plot(bin_centres_he4(mask_he4).*scale_factor,he4_flux_avg_all(plt_indx,mask_he4),'kx',...
    %             bin_centres_he4(~mask_he4).*scale_factor,he4_flux_avg_all(plt_indx,~mask_he4),'rx',...
    %             bin_centres_he4.*scale_factor,new_he4_flux,'linewidth',2)
    plot(bin_centres_he4(mask_he4).*scale_factor,new_he4_flux(:,mask_he4),'kx',...
        bin_centres_he4(~mask_he4).*scale_factor,new_he4_flux(:,~mask_he4),'rx')
    hold on
    plot(bin_centres_he4.*scale_factor,fit_flux,'linewidth',2)
else
    plot(bin_centres_he4(mask_he4).*scale_factor,new_he4_flux(:,mask_he4),'kx',...
        bin_centres_he4(~mask_he4).*scale_factor,new_he4_flux(:,~mask_he4),'rx')
    hold on
    plot(bin_centres_he4.*scale_factor,fit_flux,'linewidth',2)
end
if ax == 2
    xlabel('$x$-axis position (mm)')
    ylabel('Flux (1/m)')
elseif ax == 3
    xlabel('$y$-axis position (mm)')
    ylabel('Flux (1/m)')
elseif ax == 4
    xlabel('$r$-axis position (mm)')
    ylabel('Flux (1/m)')
else
    xlabel('Fall time (s)')
    ylabel('Flux (Hz/Shot)')
end

set(gca,'FontSize',19)
if ax == 1
    xlim([he3_time(1),he4_time(2)])
else
    xlim([-30 30])
end
if do_bimod
    new_he4_flux_bimodal = new_he4_flux_bimodal_all(plt_indx,:).';
    plot(bin_centres_he4.*scale_factor,new_he4_flux_bimodal,'r--','linewidth',2)
end
if ax ~= 1
    subplot(2,1,2)
    box on
end
% plot(bin_centres_he3.*scale_factor,he3_flux_avg_all(plt_indx,:).','bx')
% plot(bin_centres_he3.*scale_factor,new_he3_flux,'bx')
mask_he3= logical(mask_he3);
plot(bin_centres_he3(mask_he3).*scale_factor,new_he3_flux(:,mask_he3),'kx',...
    bin_centres_he3(~mask_he3).*scale_factor,new_he3_flux(:,~mask_he3),'bx')
hold on
plot(bin_centres_he3.*scale_factor,fit_flux_3,'linewidth',2)

if ax == 2
    xlabel('$x$-axis position (mm)')
    ylabel('Flux (1/m)')
elseif ax == 3
    xlabel('$y$-axis position (mm)')
    ylabel('Flux (1/m)')
else
    xlabel('Fall time (s)')
    ylabel('Flux (Hz/Shot)')
end

set(gca,'FontSize',19)
if ax == 1
    xlim([he3_time(plt_indx,1),he4_time(plt_indx,2)])
else
    xlim([-30 30])
end


he3_fit_avg

%% temperature vs threshold
threshold =linspace(0.000,0.018,16);% linspace(0.0025,0.010,16);%0.008
for plt_indx = 1:6
    he3_flux_avg = mean(cell2mat(flux_he3_all{plt_indx}(:,:)),2);
    he4_flux_avg = mean(cell2mat(flux_he4_all{plt_indx}(:,:)),2);
    if ax == 1
        initials = [1e-7 0.44 1e3];
    else
        initials = [1e-7 0.65 1e3 -0.012];
    end
    [he4_fit_avg,new_he4_flux,mask_he4,he_fit_threshold] = fit_model(bin_centres_he4,he4_flux_avg,initials,threshold,N_he4_avg(plt_indx),'he4','gauss',ax);
    [he3_fit_avg,new_he3_flux,mask_he3,he_fit_threshold_3] = fit_model(bin_centres_he3,he3_flux_avg,initials,threshold,N_he3_avg(plt_indx),'he3','gauss',ax);

    stfig('Temperature vs threshold he3');
    % clf
    % scatter(threshold,he_fit_threshold(:,1))
    hold on
    scatter(threshold.*1e3,he_fit_threshold_3(:,1).*1e9)
    legend(['he4 ',num2str(variable(plt_indx))])
    xlabel('threshold')
    ylabel('Fit temperature (nK)')
    stfig('Temperature diff vs threshold');

    stfig('Temperature vs threshold he4');
    % clf
    scatter(threshold.*1e3,he_fit_threshold(:,1).*1e9)
    hold on
    % scatter(threshold,he_fit_threshold_3(:,1))
    legend(['he4 ',num2str(variable(plt_indx))])
    xlabel('threshold')
    ylabel('Fit temperature (nK)')
    stfig('Temperature diff vs threshold');
    % clf
    scatter(threshold,1e9.*(he_fit_threshold(:,1)-he_fit_threshold_3(:,1)))
    hold on
    grid on
    legend(['he4 -he3 ',num2str(variable(plt_indx))])
    xlabel('threshold')
    ylabel('Fit temperature diff (nK)')
end
%%
threshold = linspace(0.000,0.004,20);%linspace(0.000,0.021,20);%linspace(0.0025,0.010,16);%linspace(0.000,0.018,16);%0.008
plt_indx = 4;
he3_flux_avg = mean(cell2mat(flux_he3_all{plt_indx}(:,:)),2);
he4_flux_avg = mean(cell2mat(flux_he4_all{plt_indx}(:,:)),2);
if ax == 1
    initials = [1e-7 0.04 1e5];
    initials_he3 = [1e-7 0.005 5e2];
else
    threshold = linspace(0.000,0.017,30);%
    initials = [1e-7 0.65 1e3 -0.012];
    initials_he3 = [1e-7 0.65 1e3 -0.012];
end
[he4_fit_avg,new_he4_flux,mask_he4,he_fit_threshold,avg_unc,avg_std,temp_unc] = fit_model(bin_centres_he4,he4_flux_avg,initials,threshold,N_he4_avg(plt_indx),'he4','gauss',ax);
[he3_fit_avg,new_he3_flux,mask_he3,he_fit_threshold_3] = fit_model(bin_centres_he3,he3_flux_avg,initials_he3,threshold,N_he3_avg(plt_indx),'he3','gauss',ax);

temp_unc = sqrt(temp_unc.^2+16e-9^2);
measure_unc = sqrt(avg_unc(1).^2+avg_std(1).^2+16e-9^2);
%nice plot
stfig(['Temperature vs threshold nice ',axis_labe{ax},'-axis']);
clf
% scatter(threshold,he_fit_threshold(:,1))
% scatter(threshold.*1e3,he_fit_threshold_3(:,1).*1e9,'kx')
hold on
box on

curve1 = ones(size(threshold)).*(T_he4(plt_indx)+measure_unc).*1e9';
curve2 = ones(size(threshold)).*(T_he4(plt_indx)-measure_unc).*1e9';
x1 = threshold.*1.1.*1e3';
x1(1) = -1;
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2,inBetween, 'g');
h.FaceColor = [0.31 0.31 0.32].*2;
h.FaceAlpha = 0.5;

plot(threshold.*1.1.*1e3,ones(size(threshold)).*T_he4(plt_indx).*1e9,'k-','linewidth',2)
plot(threshold.*1.1.*1e3,ones(size(threshold)).*(T_he4(plt_indx)+measure_unc).*1e9,'color',[1,1,1].*0.5,'linewidth',0.8)
plot(threshold.*1.1.*1e3,ones(size(threshold)).*(T_he4(plt_indx)-measure_unc).*1e9,'color',[1,1,1].*0.5,'linewidth',0.8)

plot(ones(size(threshold)).*2,linspace(-50,300,30),'k--','linewidth',2)
plot(ones(size(threshold)).*11,linspace(-50,300,30),'b:','linewidth',2)

errorbar(threshold.*1e3,he_fit_threshold(:,1).*1e9,temp_unc(:,1).*1e9,'r.')
%     legend(['he4 ',num2str(variable(plt_indx))])
xlabel('Radius removed (mm)')
ylabel('Fitted Temperature (nK)')
set(gca,'FontSize',17)
grid on
xlim([0 18])
ylim([0 250])

% legend('He$^3$','He$^4$')
% toc
%% Fit for each shot in detuning
% plt_index = 1;
% he3_flux_all_shot = flux_he3_all{plt_index};
% he4_flux_all_shot = flux_he4_all{plt_index};
% N_he3_all_shot = N_he3_all{plt_index};
% N_he4_all_shot = N_he4_all{plt_index};
% [r,c] = size(he3_flux_all_shot);
% clear he4_temp he3_temp initials he3_flux_pershot he4_flux_pershot N_he3_pershot N_he4_pershot threshold he4_num
% stfig(['flux fit per shot per ',var_str]);
% clf
% j = 1;
% for i=1:c
%     he3_flux_pershot = he3_flux_all_shot{i};
%     he4_flux_pershot = he4_flux_all_shot{i};
%     N_he3_pershot = N_he3_all_shot(i);
%     N_he4_pershot = N_he4_all_shot(i);
%     if isempty(he3_flux_pershot) || isempty(he4_flux_pershot)
%         continue
%     end
%     %He4 fit
%     if ax == 1
%         initials = [1e-7 0.65 1e3];
%     else
%         initials = [1e-7 0.65 1e3 0];
%     end
%     threshold = 0.01;%linspace(0.001,0.006,150);%0.008
%     [he4_fit_avg_pershot,new_he4_flux_pershot,mask_c_he4] = fit_model(bin_centres_he4,he4_flux_pershot,initials,threshold,N_he4_pershot,'he4','gauss',ax);
%     if abs(he4_fit_avg_pershot(1,1)) > 1e-5 || imag(he4_fit_avg_pershot(1,1))>0
%         continue
%     end
%     he4_temp(j,1) = he4_fit_avg_pershot(1,1);
%     he4_num(j,1) = N_he4_all{plt_index}(i);
%
%     %     clear initials
%
%     %He3 fit
%     if ax == 1
%         initials = [1e-7 0.627 max(he3_flux_pershot)./1e2];
%     else
%         initials = [1e-7 0.627 max(he3_flux_pershot)./1e2 0];
%     end
%     [he3_fit_pershot,new_he3_flux_pershot,mask_c_he3] = fit_model(bin_centres_he3,he3_flux_pershot,initials,threshold,N_he3_pershot,'he3','gauss',ax);
%
%     he3_temp(j,1)= he3_fit_pershot(1,1);
%     he3_num(j,1) = N_he3_all{plt_index}(i);
%
%     j = j+1;
%     hold on
%     plot(bin_centres_he4,new_he4_flux_pershot)
%     plot(bin_centres_he3,new_he3_flux_pershot)
%     xlabel('time')
%     ylabel('flux')
% end
%
% Temp_N = Nt(plt_index).^(1/3)./(k./(hb.*omega_bar_he4));%temperature from thermal number
%
% stfig('temp vs avg');
% clf
% subplot(2,1,1)
% scatter(he4_num,he4_temp)
% hold on
% plot([min(he4_num) max(he4_num)],[T_he4(plt_index) T_he4(plt_index)])
% plot([min(he4_num) max(he4_num)],[nanmean(he4_temp) nanmean(he4_temp)])
% plot([min(he4_num) max(he4_num)],[Temp_N Temp_N])
% xlabel('var')
% ylabel('Temperature')
% legend('individual shots','fit of avg','avg of shots')
% subplot(2,1,2)
% scatter(he3_num,he3_temp)
% hold on
% plot([min(he3_num) max(he3_num)],[T_he3(plt_index) T_he3(plt_index)])
% plot([min(he3_num) max(he3_num)],[nanmean(he3_temp) nanmean(he3_temp)])
% xlabel('var')
% ylabel('Temperature')
% legend('individual shots','fit of avg','avg of shots')

%% fitting function
function [he_fit_avg,new_he_flux,mask,he_fit,he_fit_unc_avg,he_fit_std,he_fit_unc] = fit_model(bin_centre,he_flux_avg,initials,threshold,N_he,ele,model,ax)
global const
g=const.g0;
m4 = const.mhe;
m3 = 5.008237293822000e-27;
if strcmp(ele,'he4')
    t0_cen = 0.65;
    mass = m4;
else
    t0_cen = 0.627;
    mass = m3;
end

warning('off','all')

thermal = @(b,x,m) min(ones(size(x)).*1e10,thermal_dist(b,x,m,ax));
thomas_fermi = @(b,t) thomas_fermi_dist(b,t,ax);

full_model = @(b,x) min(ones(size(x)).*1e10,thermal(b,x,mass));
modelfun_he4=@(b,t) thomas_fermi(b,t) + full_model(b,t);
modelfun_he3=@(b,x) min(ones(size(x)).*1e10,fermi_dirac_dist(b,x,mass,ax));

t_cen = trapz(bin_centre,bin_centre.*he_flux_avg)./trapz(bin_centre,he_flux_avg);
fo = statset('TolFun',10^-8,...
    'TolX',10^-10,...
    'MaxIter',10^10,...
    'UseParallel',1);

if N_he > 300%2000
    for threshold_c = threshold
        if ax == 4
            mask = abs(bin_centre)>threshold_c;
        else
            if strcmp(ele,'he4')
            mask = abs(t_cen-bin_centre)>threshold_c;
            else
                mask = abs(t_cen-bin_centre)>threshold_c;%(bin_centre>threshold_c);% | (bin_centre>-20e-3 & bin_centre<-15e-3);
            end
        end
        try
            [he_fit_c,~,J,CoV] = nlinfit(bin_centre(mask),he_flux_avg(mask),full_model,initials);
        catch
            try
                if ax == 1
                    initials = [5e-7 t0_cen max(he_flux_avg)./1e2];
                else
                    initials = [5e-7 t0_cen max(he_flux_avg)./1e2,-0.012];
                end
                [he_fit_c,~,J,CoV] = nlinfit(bin_centre,he_flux_avg,full_model,initials);
            catch
                he_fit_c = nan.*initials;
                CoV = he_fit_c;
            end
        end
        if threshold_c == threshold(1)
            he_fit = he_fit_c;
            he_fit_unc = diag(CoV.^(0.5)).';
        else
            he_fit = [he_fit;he_fit_c];
            he_fit_unc = [he_fit_unc;diag(CoV.^(0.5)).'];
        end

        if size(he_fit,1) == 1
            he_fit_std = he_fit.*0;
        else
            he_fit_std = nanstd(he_fit,1)./sqrt(size(he_fit,1));
        end
        he_fit_unc_avg = nanmean(he_fit_unc,1);
    end
    he_fit_avg = nanmean(he_fit,1);
elseif N_he > 200
    mask = ones(size(he_flux_avg));
    [he_fit,~,J,CoV] = nlinfit(bin_centre,he_flux_avg,full_model,initials);
    he_fit_avg = he_fit;
    he_fit_unc_avg = diag(CoV.^(0.5)).';
    he_fit_std = he_fit_unc_avg;
else
    mask = ones(size(he_flux_avg));
    he_fit = initials.*nan;
    he_fit_avg = he_fit;
    he_fit_unc_avg = he_fit;
    he_fit_std = he_fit_unc_avg;
end

new_he_flux = full_model(he_fit_avg,bin_centre);

if contains(model,'bimodal')
    if strcmp(ele,'he4')
        mask = ones(size(he_flux_avg));
        if sum(isnan(he_fit_avg))>1
            clear he_fit_avg
            T = 1e-6;
            sig_guess = 5e-3;
            mu_guess = t_cen;
            amp_guess = max(he_flux_avg)./1e2;
            if ax == 1
                initials_bimod = [T,mu_guess,amp_guess/2,sig_guess./4,max(he_flux_avg)/2];
            else
                initials_bimod = [T,0.65,amp_guess/2,mu_guess,sig_guess,max(he_flux_avg)/2];
            end
            fitobject=fitnlm(bin_centre,he_flux_avg,modelfun_he4,...
                initials_bimod,'Options',fo);
            he_fit_avg = fitobject.Coefficients.Estimate;
            new_he_flux = modelfun_he4(fitobject.Coefficients.Estimate,bin_centre);
        else
            sig_guess = 1;%sqrt(he_fit_avg(1)*c/m4)./10;
            initials_bimod = [sig_guess,max(he_flux_avg)];
            mdl_func = @(b,t) modelfun_he4( [he_fit_avg,b],t);
            fitobject=fitnlm(bin_centre,he_flux_avg,mdl_func,...
                initials_bimod,'Options',fo);
            he_fit_avg = [he_fit_avg,fitobject.Coefficients.Estimate.'];
            new_he_flux = modelfun_he4(he_fit_avg,bin_centre);
        end
    elseif strcmp(ele,'he3') % TO DO: add proper DFG function
        for threshold_c = threshold
            if ax == 4
                mask = abs(bin_centre)>threshold_c;
            else
                mask = abs(t_cen-bin_centre)>threshold_c;%(bin_centre>threshold_c);% | (bin_centre>-20e-3 & bin_centre<-15e-3);%abs(t_cen-bin_centre)>threshold_c;
            end
            if sum(isnan(he_fit))>1
                T = 1e-6;
                sig_guess = 5e-3;
                mu_guess = t_cen;
                amp_guess = max(he_flux_avg)./1e2;
                initials_bimod = [T,mu_guess,amp_guess/2,max(he_flux_avg)/2];
                fitobject=fitnlm(bin_centre(mask),he_flux_avg(mask),modelfun_he3,...
                    initials_bimod,'Options',fo);
                he_fit_avg = fitobject.Coefficients.Estimate;
                new_he_flux = modelfun_he3(fitobject.Coefficients.Estimate,bin_centre);
            else
                T = he_fit_avg(1);
                %                     sig_guess = sqrt(he_fit(1)*c/m3);
                mu_guess = he_fit_avg(2);
                amp_guess = he_fit_avg(3);
                %
                %             initials_bimod = [T,-7e-3,amp_guess,2];%[200e-9,0.0,4e2,3e-29];%[T,mu_guess,amp_guess/2];

                %             modelfun_he3_temp = @(b,x) modelfun_he3([1.826458958617504e-07,b(1),b(2),b(3)],x);

                if strcmp(model,'bimodal')
                    initials_bimod = [T,initials(2),amp_guess,1];
                    modelfun_he3_temp = modelfun_he3;
                    fitobject=fitnlm(bin_centre.',he_flux_avg,modelfun_he3,...
                        initials_bimod,'Options',fo);
                elseif strcmp(model,'bimodal_T_fix')
                    modelfun_he3_temp = @(b,x) modelfun_he3([initials(1),b(1),b(2),b(3)],x);
                    fitobject=fitnlm(bin_centre(mask).',he_flux_avg(mask),modelfun_he3_temp,...
                        [initials(2),amp_guess,1],'Options',fo);%[-7e-3,amp_guess,2]
                else
                    if ax == 4
                        modelfun_he3_temp = @(b,x) modelfun_he3([b(1),0,b(2),initials(4)],x);
                        fitobject=fitnlm(bin_centre(mask).',he_flux_avg(mask),modelfun_he3_temp,...
                        [T,amp_guess],'Options',fo);%[-7e-3,amp_guess,2]
                    else
                        modelfun_he3_temp = @(b,x) modelfun_he3([b(1),b(2),b(3),initials(4)],x);
                        fitobject=fitnlm(bin_centre(mask).',he_flux_avg(mask),modelfun_he3_temp,...
                        [T,initials(2),amp_guess],'Options',fo);%[-7e-3,amp_guess,2]
                    end
                    
                end
                he_fit_c = fitobject.Coefficients.Estimate.';
                if threshold_c == threshold(1)
                    he_fit = he_fit_c;
                    %             he_fit_unc = diag(CoV.^(0.5)).';
                else
                    he_fit = [he_fit;he_fit_c];
                    %             he_fit_unc = [he_fit_unc;diag(CoV.^(0.5)).'];
                end
            end
            he_fit_avg = nanmean(he_fit,1);
            new_he_flux = modelfun_he3_temp(he_fit_avg,bin_centre);
            %             new_he_flux = modelfun_he3(fitobject.Coefficients.Estimate,bin_centre);
        end
    end
end

warning('on','all')
end


%% cleaning shots function

function [he3_flux_clean, he4_flux_clean,N_he3_clean,N_he4_clean] = clean_data(flux_he3_all,flux_he4_all,bin_centres_he4,bin_centres_he3,N_he3_all ,N_he4_all,ax )
[~,c1] = size(flux_he3_all);
num_lim = 10;
for jj= 1:c1
    he3_flux_all_shot = flux_he3_all{jj};
    he4_flux_all_shot = flux_he4_all{jj};
    N_he3_all_shot = N_he3_all{jj};
    N_he4_all_shot = N_he4_all{jj};
    N_he4_mean = nanmean(N_he4_all_shot);
    N_he3_mean = nanmean(N_he3_all_shot);
    N_he4_std = nanstd(N_he4_all_shot);
    N_he3_std = nanstd(N_he3_all_shot);
    [~,c] = size(he3_flux_all_shot);
    j = 1;
    is_shot_good = [];
    %     clear he3_flux_perd he4_flux_perd initials he3_flux_pershot he4_flux_pershot N_he3_pershot N_he4_pershot threshold
    for i=1:c
        he3_flux_pershot = he3_flux_all_shot{i};
        he4_flux_pershot = he4_flux_all_shot{i};
        N_he3_pershot = N_he3_all_shot(i);
        N_he4_pershot = N_he4_all_shot(i);

        empty_check = isempty(he3_flux_pershot) || isempty(he4_flux_pershot);
        num_check = N_he3_pershot>num_lim && N_he4_pershot>num_lim;
        deviation_check = (N_he4_mean-N_he4_pershot)/N_he4_std<4;

        is_shot_good(i) = ~empty_check && num_check && deviation_check;

        if is_shot_good(i)
            he4_flux_perd{j} = he4_flux_pershot;
            N_he4_perd(j) = N_he4_pershot;

            he3_flux_perd{j} = he3_flux_pershot;
            N_he3_perd(j) = N_he3_pershot;

            j = j+1;
        end
    end
    is_shot_good_full{jj} = is_shot_good;
    he3_flux{jj} = he3_flux_perd;
    N_he3{jj} = N_he3_perd;
    he4_flux{jj} = he4_flux_perd;
    N_he4{jj} = N_he4_perd;


    he3_flux_clean = he3_flux;
    he4_flux_clean = he4_flux;
    N_he3_clean = N_he3;
    N_he4_clean = N_he4;
end
end
%% functions
% full_model_space = thermal_space(b,x,m)
% full_model_time = thermal_time(b,t,m)
function thermal = thermal_dist(b,x,m,ax)
global const
l=const.fall_distance;%0.847;
g=const.g0;
k= const.kb;
if m == 5.008237293822000e-27
    t0 = 0.37;
else
    t0 = 0.388;
end
% t0 = sqrt(2*l/g);

A = @(b,m) (m./(2*pi.*k.*b(1))).^(3/2); % A
v0 = @(b,m) sqrt((2.*k.*b(1)./m)) ; % v

part_1 = @(b,t)  (0.5*g.*(t-b(2)).^2 +l)./((t-b(2)).^2);
part_2 = @(b,t,m)  exp(-((0.5*g.*(t-b(2)).^2 -l).^2)./(v0(b,m).^2.*(t-b(2)).^2));
part_3 = @(b,x,t,m)  exp(-((x-b(4)).^2)./(v0(b,m).^2.*(t-b(2)).^2));

if ax == 1 % time axis
    thermal = b(3).*A(b,m).*v0(b,m).^2.*pi.*part_1(b,x).*part_2(b,x,m);
else % spatial axes
    %     tau = linspace(0.1.*t0,2.*t0,1e3);
    %     thermal = b(3)*A(b,m).*v0(b,m).*sqrt(pi).*trapz(tau,(0.5.*g.*tau.^2+l)./tau.^3.*part_2(b,tau,m).*part_3(b,x,tau,m),2);%for a spatial dimension (either x or y)
    thermal = b(3)*A(b,m).*v0(b,m).*pi.* exp(-((x-b(4)).^2)./(v0(b,m).^2.*t0.^2));
end
end
%
function fermi_dirac = fermi_dirac_dist(b,x,m,ax)
global const
% m = const.mhe*3/4;
hb = const.hb;
kb = const.kb;
l=const.fall_distance;%0.847;
g=const.g0;
if m == 5.008237293822000e-27
    t0 = 0.37;
else
    t0 = 0.388;
end
% t0 = sqrt(2*l/g);

omega = [60 600 600].*2*pi; %trapping frequency
omegar = 600*2*pi;
lambda = 1/10;
a0 = sqrt(hb./(m.*omega));
omega_bar = geomean(omega);
N = 1e5;%b(4);%atom number
T = b(1);%temperature

EF = hb.*omega_bar.*(6*N).^(1/3);
KF = (2*m*EF./hb^2).^0.5;%fermi wavevector
TF = EF/kb;
mu = b(4);%find_mu(omega_bar,T,N,EF);

H = @(k,x,y,z) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omega(1).^2.*x.^2+omega(2).^2.*y.^2+omega(3).^2.*z.^2);
Hr = @(k,r) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omegar.^2.*r.^2);

pr = @(k,r) 4*pi.*1./lambda.*1/(2*pi)^3.*r.^2./(exp((Hr(k,r)-mu)./(kb.*T))+1);

% set up grid
r_vec = linspace(0,a0(1).*10,1.5e3).';
r_vec_2 = linspace(0,30e-3,2.5e1);%30e-3


% Fermi at finite T (with reduced integral
vx = @(x,t) x./t;
vy = @(y,t) y./t;
vr = @(r,t) r./t;
vz = @(t) (0.5.*g.*t.^2-l)./t;
% t_vec = linspace(0.4,0.46,2e2);

n_F_r = @(k) int_loop_r(pr,k,r_vec);

n_tof = @(t,r) r.*n_F_r(m.*sqrt(vz(t).^2+vr(r,t).^2)./hb).*(0.5*g*t.^2+l)./t.^4;
%
% fermi_dist = @(t) 2.*pi.*int_loop_r(n_tof,t,r_vec_2);
v0 = @(b,m) sqrt((2.*kb.*b(1)./m)) ; % v
if ax == 1
    fermi_dist = @(b,t) -polylog(2.5,-abs(b(4)).*exp(-t.^2./(v0(b,m).^2.*t0.^2)));
elseif ax~= 4
    fermi_dist = @(b,t) -real(polylog(2.5,-abs(b(4)).*exp(-t.^2./(v0(b,m).^2.*t0.^2)))) - 494.*sqrt(b(1)).*real(polylog(3.5,-abs(b(4)).*exp(-t.^2./(v0(b,m).^2.*t0.^2))));%2.5;
else
    fermi_dist = @(b,x) -real(polylog(2,-abs(b(4)).*exp(-x.^2./(v0(b,m).^2.*t0.^2))));% -b(5).*real(polylog(3,-abs(b(4)).*exp(-x.^2./(v0(b,m).^2.*t0.^2))));
end

fermi_dirac = fermi_dist(b,x-b(2));
fermi_dirac = b(3).*fermi_dirac./trapz(x,fermi_dirac);
% plot(t_vec,tshermal_dist([T,0,0.001],t_vec,m,1))

end

% function out = fermi_function(b,t,m)
% global const
% hb = const.hb;
% kb = const.kb;
% l=const.fall_distance;%0.847;
% g=const.g0;
% if m == 5.008237293822000e-27
%     t0 = 0.37;
% else
%     t0 = 0.388;
% end
% % t0 = sqrt(2*l/g);
% 
% omega = [60 600 600].*2*pi; %trapping frequency
% omegar = 600*2*pi;
% lambda = 1/10;
% a0 = sqrt(hb./(m.*omega));
% omega_bar = geomean(omega);
% N = 1e5;%b(4);%atom number
% T = b(1);%temperature
% 
% EF = hb.*omega_bar.*(6*N).^(1/3);
% KF = (2*m*EF./hb^2).^0.5;%fermi wavevector
% TF = EF/kb;
% mu = b(4);%find_mu(omega_bar,T,N,EF);
% v0 = @(b,m) sqrt((2.*kb.*b(1)./m)) ; % v
% try
%     out = -real(polylog(2.5,-abs(b(4)).*exp(-t.^2./(v0(b,m).^2.*t0.^2)))) - 494.*sqrt(b(1)).*real(polylog(3.5,-abs(b(4)).*exp(-t.^2./(v0(b,m).^2.*t0.^2))));%2.5;
% catch
%     out = t.*0;
% end
% end

function mu = find_mu(omega_bar,T,N,EF)
global const
kb = const.kb;
hb = const.hb;

mu_temp = linspace(-4.*EF,EF,5e3).';
N_E = @(E,mu_vec) E.^2./(2*(hb.*omega_bar).^3).*1./(exp((E-mu_vec)./(kb.*T))+1);

E = linspace(0,20.*EF,5e3);

N_temp = abs(trapz(E,N_E(E,mu_temp),2)-N);
[val,indx] = min(N_temp);
mu = mu_temp(indx);

end


function n_vec = int_loop_r(p,k,r_vec)
for ii = 1:length(k)
    n_vec(ii) = trapz(r_vec,p(k(ii),r_vec));
end
end


function thomas_fermi = thomas_fermi_dist(b,t,ax)
global const
l=const.fall_distance;%0.847;
g=const.g0;
t0 = sqrt(2*l/g);

if ax == 1
    thomas_fermi = min(ones(size(t)).*1e10,b(5).*real(max(zeros(size(t)),(1-((t-b(2)-t0)./b(4)).^2)).^(3/2)));
else
    thomas_fermi = min(ones(size(t)).*1e10,b(6).*real(max(zeros(size(t)),(1-((t-b(4))./b(5)).^2)).^(3/2)));
end
end

function [y errors] = polylog(n,z)
%% polylog - Computes the n-based polylogarithm of z: Li_n(z)
% Approximate closed form expressions for the Polylogarithm aka de
% Jonquiere's function are used. Computes reasonably faster than direct
% calculation given by SUM_{k=1 to Inf}[z^k / k^n] = z + z^2/2^n + ...
%
% Usage:   [y errors] = PolyLog(n,z)
%
% Input:   z < 1   : real/complex number or array
%          n > -4  : base of polylogarithm
%
% Output: y       ... value of polylogarithm
%         errors  ... number of errors
%
% Approximation should be correct up to at least 5 digits for |z| > 0.55
% and on the order of 10 digits for |z| <= 0.55!
%
% Please Note: z vector input is possible but not recommended as precision
% might drop for big ranged z inputs (unresolved Matlab issue unknown to
% the author).
%
% following V. Bhagat, et al., On the evaluation of generalized
% Bose胞instein and Fermi縫irac integrals, Computer Physics Communications,
% Vol. 155, p.7, 2003
%
% v3 20120616
% -------------------------------------------------------------------------
% Copyright (c) 2012, Maximilian Kuhnert
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------
if nargin~=2
    errors=1;
    error('[Error in: polylog function] Inappropriate number of input arguments!')
end
if (isreal(z) && sum(z(:)>=1)>0) % check that real z is not bigger than 1
    errors=1;
    error('[Error in: polylog function] |z| > 1 is not allowed')
elseif isreal(z)~=1 && sum(abs(z(:))>1)>0 % check that imaginary z is defined on unit circle
    errors=1;
    error('[Error in: polylog function] |z| > 1 is not allowed')
elseif n<=-4 % check that n is not too largly negative (see paper)
    errors=1;
    error('[Error in: polylog function] n < -4 might be inaccurate')
end
% display more digits in Matlab terminal:
%format long
alpha = -log(z); % see page 12
% if |z| > 0.55 use Eq. (27) else use Eq. (21):
if abs(z) > 0.55
    preterm = gamma(1-n)./alpha.^(1-n);
    nominator = b(0) + ...
        - alpha.*( b(1) - 4*b(0)*b(4)/7/b(3) ) + ...
        + alpha.^2.*( b(2)/2 + b(0)*b(4)/7/b(2) - 4*b(1)*b(4)/7/b(3) ) + ...
        - alpha.^3.*( b(3)/6 - 2*b(0)*b(4)/105/b(1) + b(1)*b(4)/7/b(2) - 2*b(2)*b(4)/7/b(3) );
    denominator = 1 + alpha.*4*b(4)/7/b(3) +...
        + alpha.^2.*b(4)/7/b(2) +...
        + alpha.^3.*2*b(4)/105/b(1) +...
        + alpha.^4.*b(4)/840/b(0);
    y = preterm + nominator ./ denominator;
else
    nominator = 6435*9^n.*S(n,z,8) - 27456*8^n*z.*S(n,z,7) + ...
        + 48048*7^n*z.^2.*S(n,z,6) - 44352*6^n*z.^3.*S(n,z,5) + ...
        + 23100*5^n*z.^4.*S(n,z,4) - 6720*4^n.*z.^5.*S(n,z,3) + ...
        + 1008*3^n*z.^6.*S(n,z,2) - 64*2^n*z.^7.*S(n,z,1);
    denominator = 6435*9^n - 27456*8^n*z + ...
        + 48048*7^n*z.^2 - 44352*6^n*z.^3 + ...
        + 23100*5^n*z.^4 - 6720*4^n*z.^5 + ...
        + 1008*3^n*z.^6 - 64*2^n*z.^7 + ...
        + z.^8;
    y = nominator ./ denominator;
end
% define b:
    function out = b(i)
        out = zeta(n-i);
    end
% define S as partial sums of Eq. 12:
    function out = S(n,z,j)
        out =0;
        for i=1:j
            out = out + z.^i./i^n;
        end
    end
    function [out] = zeta(x)
        %% Zeta Function
        % Eq. 18
        % following V. Bhagat, et al., On the evaluation of generalized
        % Bose胞instein and Fermi縫irac integrals, Computer Physics Communications,
        % Vol. 155, p.7, 2003
        %
        % Usage: [out] = zeta(x)
        % with argument x and summation from 1 to j
        %
        % MK 20120615
        prefactor = 2^(x-1) / ( 2^(x-1)-1 );
        numerator = 1 + 36*2^x*eta(x,2) + 315*3^x*eta(x,3) + 1120*4^x*eta(x,4) +...
            + 1890*5^x*eta(x,5) + 1512*6^x*eta(x,6) + 462*7^x*eta(x,7);
        denominator = 1 + 36*2^x + 315*3^x + 1120*4^x + 1890*5^x + 1512*6^x +...
            + 462*7^x;
        out = prefactor * numerator / denominator;
        function [out] = eta(x,j)
            %% Eta Function
            % Eq. 17 (partial sums)
            % following V. Bhagat, et al., On the evaluation of generalized
            % Bose-Einstein and Fermi-Dirac integrals, Computer Physics Communications,
            % Vol. 155, p.7, 2003
            %
            % Usage: [out] = eta(x,j)
            % with argument x and summation from 1 to j
            %
            % MK 20120615

            out=0;
            for k=1:j
                out = out + (-1)^(k+1) ./ k.^x;
            end
        end
    end
end

