%% Initializing path
clear all;
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
data_dir = '20220328_TF_dep_trial_5_fracNwizard_detuning';
% data_folder = {'20220324_he3_and_4_seperated_clouds'};
% Find all the data folders in the data directory

f_winfreak = 10767.5;%in MHz
% expected detuning 33.574 GHz

if exist('data_dir','var')
    f_indx = 1;
    full_data_dir = fullfile(opts.data_root, data_dir);
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

% data_folder = {'20220324_he3_and_4_seperated_clouds'};


opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;
% opts.import.shot_num = 26:50;
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
%% Analysis Params

he4_time = [1.05,1.09];
he3_time = [1.01,1.05];

min_num = 1e2;

ax = 2; %which axis do we analyse

do_bimod = 0;
plot_all = 0;
show_fit_goodness = 0;

plot_model_comparison = 0;

%% separate He 3 and He 4 distrabutions
dirs_list = fullfile(opts.data_root, data_folder);
if iscell(dirs_list)
    num_dirs = length(dirs_list);
else
    num_dirs = 1;
end
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
if ~cahce_exists
    for jj = 1:num_dirs
        %% Import data
        if iscell(dirs_list)
            opts.import.dir = dirs_list{jj};
            opts.import.cache_save_dir = fullfile(dirs_list{jj}, 'cache', 'import\');
            [data, ~] = import_mcp_tdc_data(opts.import);
        else
            opts.import.dir = fullfile(opts.data_root, data_folder);
            [data, ~] = import_mcp_tdc_data(opts.import);
        end
        %% remove any hotspts

        data_ht_spot=hotspot_mask(data);
        data.counts_txy=data_ht_spot.masked.counts_txy;
        data.num_counts=data_ht_spot.masked.num_counts;

        clear flux_he4 flux_he3 N_he3 N_he4
        for ii = 1:length(data.counts_txy)
            lims_4= [he4_time; -0.03, 0.03; -0.03, 0.03];
            he4_txy = masktxy_square(data.counts_txy{ii}, lims_4);
            lims_3 = [he3_time; -0.03, 0.03; -0.03, 0.03];
            he3_txy = masktxy_square(data.counts_txy{ii}, lims_3);

            t_he4 = linspace(he4_time(1),he4_time(2),5e3).';
            t_he3 = linspace(he3_time(1),he3_time(2),5e3).';

            space_bins = linspace(-30e-3,30e-3,5e3).';

            if ax == 1
                bin_cen_4 = t_he4;
                bin_cen_3 = t_he3;
            else
                bin_cen_4 = space_bins;
                bin_cen_3 = space_bins;
            end

            %% run good shot checks
            is_shot_good = (size(he4_txy,1)+size(he3_txy,1))>min_num;
            shot_check(ii) = is_shot_good;
            if is_shot_good
                %% histogram in time
                sigma = 0.5e-4;
                count_hist_he4 = smooth_hist(he4_txy(:,ax),'sigma',sigma,'edges',bin_cen_4);
                count_hist_he3 = smooth_hist(he3_txy(:,ax),'sigma',sigma,'edges',bin_cen_3);
                flux_he4{ii} = count_hist_he4.count_rate.smooth;
                flux_he3{ii} = count_hist_he3.count_rate.smooth;
                bin_centres_he4 = count_hist_he4.bin.centers;
                bin_centres_he3 = count_hist_he3.bin.centers;
                %record atom number
                N_he4(ii) = size(he4_txy(:,1),1);
                N_he3(ii) = size(he3_txy(:,1),1);
            end
        end
        flux_he3_all{jj} = flux_he3; %flux per shot per detuning
        flux_he4_all{jj} = flux_he4;
        N_he3_all{jj} = N_he3; %atom numbers per detuning
        N_he4_all{jj} = N_he4;
    end

    %% cleaning shots & saving
    [he3_flux_clean, he4_flux_clean,N_he3_clean,N_he4_clean] = clean_data(flux_he3_all,flux_he4_all,bin_centres_he4,bin_centres_he3,N_he3_all ,N_he4_all,ax);
    save(fullfile(full_data_dir,['clean_data_',num2str(ax),'.mat']),"he4_flux_clean","he3_flux_clean","N_he4_clean","N_he3_clean","N_he3_all","N_he4_all","flux_he3_all","flux_he4_all",'bin_centres_he4',"bin_centres_he3")
else
    load(fullfile(full_data_dir,['clean_data_',num2str(ax),'.mat']))
end


%% Analysing data
for kk = 1: num_dirs
    he4_flux_avg= mean(cell2mat(he4_flux_clean{kk}),2); %average of all shots for current detuning
    he3_flux_avg = mean(cell2mat(he3_flux_clean{kk}),2);
    he4_flux_avg_all(kk,:)= mean(cell2mat(he4_flux_clean{kk}),2); % average of all shots for all detuning
    he3_flux_avg_all(kk,:) = mean(cell2mat(he3_flux_clean{kk}),2);
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
        threshold = linspace(0.0025,0.012,15);%linspace(0.003,0.010,8);%0.008;
        initials = [1e-7 0.65 1e3 -0.012];
        initials_he3 = [1e-7 0.627 max(he3_flux_avg)./1e2 -0.012];
    end
    [he4_fit_avg,new_he4_flux,mask_he4,~,he4_fit_unc, he4_fit_std] = fit_model(bin_centres_he4,he4_flux_avg,initials,threshold,N_he4_avg(kk),'he4','gauss',ax);
    if do_bimod
        [he4_fit_avg_bimodal,new_he4_flux_bimodal,mask_he4_bimodal,] = fit_model(bin_centres_he4,he4_flux_avg,initials,threshold,N_he4_avg(kk),'he4','bimodal',ax);
        new_he4_flux_bimodal_all(kk,:) = new_he4_flux_bimodal;
        he4_fits_bimodal(kk,:) = he4_fit_avg_bimodal;
    end
    clear initials
    %He3 fit

    [he3_fit,new_he3_flux,mask_c_he3,~,he3_fit_unc,he3_fit_std] = fit_model(bin_centres_he3,he3_flux_avg,initials_he3,threshold,N_he3_avg(kk),'he3','gauss',ax);
    %     [he3_fit_bimodal,new_he3_flux_bimodal,mask_c_he3_bimodal] = fit_model(bin_centres_he3,he3_flux_avg,initials,threshold,N_he3_avg(kk),'he3','bimodal',ax);



    he4_fits(kk,:) = he4_fit_avg;
    he3_fits(kk,:) = he3_fit;
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

    T_he4_unc(kk,1) = sqrt(he4_fit_unc(1,1).^2+he4_fit_std(1,1).^2);
    T_he3_unc(kk,1) = sqrt(he3_fit_unc(1,1).^2+he3_fit_std(1,1).^2);

    k=const.kb;
    qe = 0.02;%0.08;
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
qe = 0.02;%0.08;
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

stfig('he3-he4 temp');
clf
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
plt_indx = 5;
variable(plt_indx)
he4_fit_plt = he4_fits(plt_indx,:);
he3_fit_plt = he3_fits(plt_indx,:);
new_he3_flux = new_he3_flux_all(plt_indx,:).';
new_he4_flux = new_he4_flux_all(plt_indx,:).';
mask_he4 = logical(he4_masks(plt_indx,:));

stfig('flux fit per indx');
clf
hold on
box on
scale_factor = 1;
if ax ~= 1
    scale_factor = 1e3;
    subplot(2,1,1)

plot(bin_centres_he4(mask_he4).*scale_factor,he4_flux_avg_all(plt_indx,mask_he4),'kx',...
            bin_centres_he4(~mask_he4).*scale_factor,he4_flux_avg_all(plt_indx,~mask_he4),'rx',...
            bin_centres_he4.*scale_factor,new_he4_flux,'linewidth',2)
else

plot(bin_centres_he4.*scale_factor,he4_flux_avg_all(plt_indx,:),'kx')

end
    if ax == 2
        xlabel('$x$-axis position (mm)')
    elseif ax == 3
        xlabel('$y$-axis position (mm)')
    else
        xlabel('Fall time (s)')
    end
    ylabel('Flux (Hz/Shot)')
    set(gca,'FontSize',19)
    if ax == 1
        xlim([1.02,1.09])
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
plot(bin_centres_he3.*scale_factor,he3_flux_avg_all(plt_indx,:).','bx')
    %bin_centres_he3.*scale_factor,new_he3_flux)
if ax == 2
    xlabel('$x$-axis position (mm)')
elseif ax == 3
    xlabel('$y$-axis position (mm)')
else
xlabel('Fall time (s)')
end
ylabel('Flux (Hz/Shot)')
set(gca,'FontSize',19)
if ax == 1
        xlim([1.02,1.09])
    else
    xlim([-30 30])
    end
%% temperature vs threshold
threshold = linspace(0.001,0.015,5);%0.008
if ax == 1
    initials = [1e-7 0.65 1e3];
else
    initials = [1e-7 0.65 1e3 -0.012];
end
[he4_fit_avg,new_he4_flux,mask_he4,he_fit_threshold] = fit_model(bin_centres_he4,he4_flux_avg,initials,threshold,N_he4_avg(plt_indx),'he4','gauss',ax);

stfig('Temperature vs threshold');
scatter(threshold,he_fit_threshold(:,1))
xlabel('threshold')
ylabel('fit temp')
%% Fit for each shot in detuning
plt_index = 1;
he3_flux_all_shot = flux_he3_all{plt_index};
he4_flux_all_shot = flux_he4_all{plt_index};
N_he3_all_shot = N_he3_all{plt_index};
N_he4_all_shot = N_he4_all{plt_index};
[r,c] = size(he3_flux_all_shot);
clear he4_temp he3_temp initials he3_flux_pershot he4_flux_pershot N_he3_pershot N_he4_pershot threshold he4_num
stfig(['flux fit per shot per ',var_str]);
clf
j = 1;
for i=1:c
    he3_flux_pershot = he3_flux_all_shot{i};
    he4_flux_pershot = he4_flux_all_shot{i};
    N_he3_pershot = N_he3_all_shot(i);
    N_he4_pershot = N_he4_all_shot(i);
    if isempty(he3_flux_pershot) || isempty(he4_flux_pershot)
        continue
    end
    %He4 fit
    if ax == 1
        initials = [1e-7 0.65 1e3];
    else
        initials = [1e-7 0.65 1e3 0];
    end
    threshold = 0.008;%linspace(0.001,0.006,150);%0.008
    [he4_fit_avg_pershot,new_he4_flux_pershot,mask_c_he4] = fit_model(bin_centres_he4,he4_flux_pershot,initials,threshold,N_he4_pershot,'he4','gauss',ax);
    if abs(he4_fit_avg_pershot(1,1)) > 1e-5 || imag(he4_fit_avg_pershot(1,1))>0
        continue
    end
    he4_temp(j,1) = he4_fit_avg_pershot(1,1);
    he4_num(j,1) = N_he4_all{plt_index}(i);

    %     clear initials

    %He3 fit
    if ax == 1
        initials = [1e-7 0.627 max(he3_flux_pershot)./1e2];
    else
        initials = [1e-7 0.627 max(he3_flux_pershot)./1e2 0];
    end
    [he3_fit_pershot,new_he3_flux_pershot,mask_c_he3] = fit_model(bin_centres_he3,he3_flux_pershot,initials,threshold,N_he3_pershot,'he3','gauss',ax);

    he3_temp(j,1)= he3_fit_pershot(1,1);
    he3_num(j,1) = N_he3_all{plt_index}(i);

    j = j+1;
    hold on
    plot(bin_centres_he4,new_he4_flux_pershot)
    plot(bin_centres_he3,new_he3_flux_pershot)
    xlabel('time')
    ylabel('flux')
end

Temp_N = Nt(plt_index).^(1/3)./(k./(hb.*omega_bar_he4));%temperature from thermal number

stfig('temp vs avg');
clf
subplot(2,1,1)
scatter(he4_num,he4_temp)
hold on
plot([min(he4_num) max(he4_num)],[T_he4(plt_index) T_he4(plt_index)])
plot([min(he4_num) max(he4_num)],[nanmean(he4_temp) nanmean(he4_temp)])
plot([min(he4_num) max(he4_num)],[Temp_N Temp_N])
xlabel('var')
ylabel('Temperature')
legend('individual shots','fit of avg','avg of shots')
subplot(2,1,2)
scatter(he3_num,he3_temp)
hold on
plot([min(he3_num) max(he3_num)],[T_he3(plt_index) T_he3(plt_index)])
plot([min(he3_num) max(he3_num)],[nanmean(he3_temp) nanmean(he3_temp)])
xlabel('var')
ylabel('Temperature')
legend('individual shots','fit of avg','avg of shots')

%% fitting function
function [he_fit_avg,new_he_flux,mask,he_fit,he_fit_unc_avg,he_fit_std] = fit_model(bin_centre,he_flux_avg,initials,threshold,N_he,ele,model,ax)
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
modelfun_he3=@(b,t) thomas_fermi(b,t) + full_model(b,t);

t_cen = trapz(bin_centre,bin_centre.*he_flux_avg)./trapz(bin_centre,he_flux_avg);
fo = statset('TolFun',10^-8,...
    'TolX',10^-10,...
    'MaxIter',10^10,...
    'UseParallel',1);

if N_he > 2000
    for threshold_c = threshold
        mask = abs(t_cen-bin_centre)>threshold_c;
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
        he_fit_avg = nanmean(he_fit,1);
        he_fit_std = nanstd(he_fit,1)./sqrt(size(he_fit,1));
        he_fit_unc_avg = nanmean(he_fit_unc,1);
    end
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

if strcmp(model,'bimodal')
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
            sig_guess = sqrt(he_fit_avg(1)*c/m4)./10;
            initials_bimod = [sig_guess,max(he_flux_avg)];
            mdl_func = @(b,t) modelfun_he4( [he_fit_avg,b],t);
            fitobject=fitnlm(bin_centre,he_flux_avg,mdl_func,...
                initials_bimod,'Options',fo);
            he_fit_avg = [he_fit_avg,fitobject.Coefficients.Estimate.'];
            new_he_flux = modelfun_he4(he_fit_avg,bin_centre);
        end
    elseif strcmp(ele,'he3') % TO DO: add proper DFG function
        %         if sum(isnan(he3_fit))>1
        %             clear he3_fit
        %             T = 1e-6;
        %             sig_guess = 5e-3;
        %             mu_guess = t_cen;
        %             amp_guess = max(he_flux_avg)./1e2;
        %             initials_bimod = [T,mu_guess,amp_guess/2,sig_guess./4,max(he_flux_avg)/2];
        %             fitobject=fitnlm(bin_centre,he_flux_avg,modelfun_he3,...
        %                 initials_bimod,'Options',fo);
        %             he_fit_avg = fitobject.Coefficients.Estimate;
        %             new_he_flux = modelfun_he3(fitobject.Coefficients.Estimate,bin_centre);
        %         else
        %             T = he3_fit(1);
        %             sig_guess = sqrt(he3_fit(1)*c/m3);
        %             mu_guess = he3_fit(2);
        %             amp_guess = he3_fit(3);
        % %             initials_bimod = [T,mu_guess,amp_guess,sig_guess./4,max(he_flux_avg)/10];
        %             initials_bimod = [he_fit_avg,sig_guess,amp_guess/2];
        %             fitobject=fitnlm(bin_centre,he_flux_avg,modelfun_he3,...
        %                 initials_bimod,'Options',fo);
        %             he_fit_avg = fitobject.Coefficients.Estimate;
        %             new_he_flux = modelfun_he3(fitobject.Coefficients.Estimate,bin_centre);
        %         end
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
t0 = sqrt(2*l/g);

A = @(b,m) (m./(2*pi.*k.*b(1))).^(3/2); % A
v0 = @(b,m) sqrt((2.*k.*b(1)./m)) ; % v

part_1 = @(b,t)  (0.5*g.*(t-b(2)).^2 +l)./((t-b(2)).^2);
part_2 = @(b,t,m)  exp(-((0.5*g.*(t-b(2)).^2 -l).^2)./(v0(b,m).^2.*(t-b(2)).^2));
part_3 = @(b,x,t,m)  exp(-((x-b(4)).^2)./(v0(b,m).^2.*(t-b(2)).^2));

if ax == 1 % time axis
    thermal = b(3).*A(b,m).*v0(b,m).^2.*pi.*part_1(b,x).*part_2(b,x,m);
else % spatial axes
    tau = linspace(0.1.*t0,2.*t0,1e3);
    thermal = b(3)*A(b,m).*v0(b,m).*sqrt(pi).*trapz(tau,(0.5.*g.*tau.^2+l)./tau.^3.*part_2(b,tau,m).*part_3(b,x,tau,m),2);%for a spatial dimension (either x or y)
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


