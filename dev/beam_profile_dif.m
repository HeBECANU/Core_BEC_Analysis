% beam profile difference


%fit probe beam
% user param

%im_raw=imread('./data/beam_profile/20180724T211918.841.png');
%im_raw=imread('C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\output_of_fiber_6.png');

% testfiledir='C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\input_slide_rail_mon_20220228\';
% filedir_1='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\beam_images\on';
% filedir_2='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\beam_images\off';
fildedir_dark = '';
filedir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\beam_images\f=125_acr_193mm_wvplt_165_field_3.0';

%need Ilight, I atoms, and Idark
hebec_constants
%physical parameters
%we use the 23S1 to 23P2 transition
Gamma=1.0216e+07;%Hz transition linewidth
det=0;%detuning
dipole_moment = 3;
omega0=2*pi*2.767321865671393e+14;
%general case for Isat and sigma0
Isat=const.c*const.epsilon0*Gamma^2*const.hb^2/(4*dipole_moment);
sigma0=const.hb*omega0*Gamma/(2*Isat);
%sigma+ light and two level system
Isat=omega0^3*Gamma*const.hb/(12*pi*const.c^2);
lambda0=2*pi*const.c/omega0;
sigma0=3*lambda0^2/(2*pi);

I0=0.5.*Isat;

sigma=sigma0/(1+4.*(det/Gamma).^2+I0/Isat);

OD_sat = 1.58;%3;%

% x=[];
% y=[];
% dist = [310,320,330,350,400,450,500,550];
% dist = [310,320,330,380,400,420,450,500,580];
% dist = 1:10;
% dist = [330, 340, 350, 400, 450, 500, 550, 560, 570];
% dist = [0,100,150,200,250,300,350,400,450,500,550,580,50];
wavelen = 1083.3e-9;
pix_size=20e-6;%3e-6;
% close('all')
% files = {'C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\input_slide_rail_0.png'};

files_total = dir(fullfile(filedir, '*.png'));
files_total = {files_total.name};
num_shots = size(files_total,2);

for j=1:2
    files = files_total((0+j):2:(num_shots-2+j));
    
    for i=1:length(files)
        im_raw=imread(fullfile(filedir,files{i}));
        %     im_raw=imread(fullfile(files{i}));


        %%
        beam_image=double(im_raw)/255;
        % beam_image=sum(beam_image,3);
        if i == 1
            total_image{j} = beam_image;
        else
            total_image{j} = total_image{j}+beam_image;
        end



        %     % find the mean of the image excluding this subframe
        %     im_background=sum(beam_image(:))-sum(subframe(:));
        %     im_background=im_background/(numel(beam_image)-numel(subframe));
        %     subframe=subframe-im_background;
        %     subframe=subframe./max(subframe(:));
    end
    total_image{j} = total_image{j}./length(files);
    %     total_image{j}=total_image{j}./max(max(total_image{j}));
    sub_size=size(total_image{j});
    sub_size_m=pix_size*sub_size;
    xvals_sub=linspace(-sub_size(2)/2,sub_size(2)/2,sub_size(2))*pix_size;
    yvals_sub=linspace(-sub_size(1)/2,sub_size(1)/2,sub_size(1))*pix_size;


    %%

    x_lim=[-1,1]*30e-6;
    y_lim=x_lim;
    xy_factor=1e6;

    font_name='cmr10';
    linewidth=1.5;
    font_size=12;

    stfig(['camera dir',num2str(j)]);
    clf
    set(gca, 'FontName', font_name)
    set(gca, 'FontSize', font_size)



    sh=surf(xvals_sub*xy_factor,yvals_sub*xy_factor,total_image{j},...
        'FaceAlpha',0.4);
    ifh=gca;
    shading interp
    sh.EdgeColor='k';
    if j ==1
        colormap(gca,plasma())
    end
    %     caxis([0 0.8])
    %pbaspect([1,1,1])
    % xlim(x_lim*xy_factor)
    % ylim(y_lim*xy_factor)
    ifh.View= [7.8   14.04];
    box on
    xlabel('y ($\mathrm{\mu}$m)')
    ylabel('z ($\mathrm{\mu}$m)')
    zlabel('Intensity (arb. u.)')
    shading interp

    stfig(['integrated profile dir:',num2str(j)]);
    clf
    subplot(2,1,1)
    plot(xvals_sub*xy_factor,trapz(yvals_sub,total_image{j}));
    xlabel('y ($\mathrm{\mu}$m)')
    ylabel('Integrated Intensity (arb. u.)')
    subplot(2,1,2)
    plot(yvals_sub*xy_factor,trapz(xvals_sub,total_image{j},2));
    xlabel('z ($\mathrm{\mu}$m)')
    ylabel('Integrated Intensity (arb. u.)')

    %%

    fig_name=['beam_profile_img_',files{i}];
    %fig_name='pal_mean_v_acc_dyn_static';
    fig_dir='C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\figs';
    %     export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
    %     export_fig(fullfile(fig_dir,strcat(fig_name,'.png')))



    %%


    %colors_main=[[98,136,63];[95,109,187];[180,72,117]]./255;
    colors_main=[[40,136,40];[95,109,187];[218, 129, 57]]./255;

    hsv=colorspace('RGB->HSV',colors_main(:,:));
    hsv(:,2)=hsv(:,2);
    colors_shaded=colorspace('HSV->RGB',hsv);

end
stfig('optical density 1');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

I_light = total_image{2};
I_atoms = total_image{1};

OD_meas = log(I_light./I_atoms);%log((I_light-I_dark)./(I_atoms-I_dark));measured optical density
%basic thresholding
% A(abs(A)>3)=nan;
% power htesholding
I_min = 1;
OD_meas(((I_atoms<I_min)|(I_light<I_min))) = nan;
% spike and averaging and interpolation

sh=surf(xvals_sub*xy_factor,yvals_sub*xy_factor,OD_meas,...
    'FaceAlpha',0.75);
ifh=gca;
shading interp
sh.EdgeColor='k';
%     caxis([0 0.8])
%pbaspect([1,1,1])
% xlim(x_lim*xy_factor)
% ylim(y_lim*xy_factor)
ifh.View= [7.8   14.04];
box on
xlabel('y ($\mathrm{\mu}$m)')
ylabel('z ($\mathrm{\mu}$m)')
zlabel('Optical density')
shading interp

maximum = max(max(OD_meas));
[x_indx,y_indx]=find(OD_meas==maximum);

[xmesh,ymesh]=meshgrid(xvals_sub,yvals_sub);
    pos_val = zeros(size(xmesh,1),size(ymesh,2),3);
    pos_val(:,:,1) = xmesh;
    pos_val(:,:,2) = ymesh;
    pos_val(:,:,3) = OD_meas; 
    pos_val_rows=reshape(pos_val,[],3,1);
    % stfig('conversion check')
 opt = statset('TolFun',1e-10,'TolX',1e-10,...
        'MaxIter',1e4,... %1e4
        'UseParallel',1);
    
    
    cof_names={'sigma x','sigma y', 'cen x','cen y', 'theta', 'offset', 'amp'};
    
    modelfun = @(b,x) gauss_2d_rot_function(x,[b(1),b(2)],[b(3),b(4)],b(5),b(6),b(7),'sigma');
    
    fit_params_guess = [std(pos_val_rows(:,1),max(pos_val_rows(:,3),0))/5,...
                        std(pos_val_rows(:,2),max(pos_val_rows(:,3),0))/5,...
                        xvals_sub(y_indx),yvals_sub(x_indx),0,-0.2,max(OD_meas(:))]; %Inital guess parameters
    fitobject=fitnlm(pos_val_rows(:,1:2),pos_val_rows(:,3),modelfun,fit_params_guess,...
        'options',opt,...
        'CoefficientNames',cof_names)
    fitparam=fitobject.Coefficients;
%     x{i} = fitparam{1,1};
%     dx{i} = fitparam{1,2};
%     y{i} = fitparam{2,1};
    hold on
    %
    %subplot(1,4,3)

    fit_pos_val_rows=pos_val_rows;
    %fit_pos_val_rows(:,3)=modelfun(fit_params_guess,pos_val_rows(:,1:2));
    fit_pos_val_rows(:,3)=predict(fitobject,pos_val_rows(:,1:2));
    fit_pos_vals=reshape(fit_pos_val_rows,size(pos_val));
    %imagesc(fit_pos_vals(:,:,3))
    sh=surf(fit_pos_vals(:,:,1)*xy_factor,fit_pos_vals(:,:,2)*xy_factor,fit_pos_vals(:,:,3),...
        'FaceAlpha',0.8,'FaceColor',[100,200,53]./255,'edgecolor','none');

stfig('OD slice');
hold on
OD_avg = mean(OD_meas(:,abs(xvals_sub*xy_factor-550)<50),2);
plot(yvals_sub*xy_factor,OD_avg)
xlabel('z ($\mathrm{\mu}$m)')
ylabel('optical density')

OD_mod = log((1-exp(-OD_sat))./(exp(-OD_meas)-exp(-OD_sat)));
OD_actual = OD_mod+(1-exp(-OD_mod)).*I0./Isat;
OD_actual(isnan(OD_actual)) = 0;



stfig('difference 4');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)


sh=surf(xvals_sub*xy_factor,yvals_sub*xy_factor,total_image{1}-total_image{2},...
    'FaceAlpha',0.75);
ifh=gca;
shading interp
sh.EdgeColor='k';
%     caxis([0 0.8])
%pbaspect([1,1,1])
% xlim(x_lim*xy_factor)
% ylim(y_lim*xy_factor)
ifh.View= [7.8   14.04];
box on
xlabel('y ($\mathrm{\mu}$m)')
ylabel('z ($\mathrm{\mu}$m)')
zlabel('Intensity (arb. u.)')
shading interp
stfig('integrated profile diff');
% clf

E = @(b,x) 1./b(2).*besselj(1,b(1).*(x-b(3)))./(x-b(3))./b(1);
fit_params_guess = [1,0.06,-200];
fitobject=fitnlm(xvals_sub*xy_factor,trapz(yvals_sub,total_image{1}-total_image{2}),E,fit_params_guess);
fitparam=fitobject.Coefficients;

subplot(2,1,1)
hold on
plot(xvals_sub*xy_factor,trapz(yvals_sub,total_image{1}-total_image{2}));
xlabel('y ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity (arb. u.)')
subplot(2,1,2)
hold on
dif_z = trapz(xvals_sub,total_image{1}-total_image{2},2);
plot(yvals_sub*xy_factor,dif_z);
xlabel('z ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity (arb. u.)')
[v1,id1]=max(dif_z(yvals_sub>0));
[v2,id2]=max(dif_z(yvals_sub<0));

(yvals_sub(id1+128)-yvals_sub(id2))*xy_factor;

[X,Y]= meshgrid(xvals_sub,yvals_sub);
X_w = X.*OD_meas;
Y_w = Y.*OD_meas;
X_OD = nansum(OD_meas,2)./size(yvals_sub,1);
Y_OD = nansum(OD_meas)./size(xvals_sub,2);
stfig('integrated OD profile');
subplot(2,1,1)
hold on
plot(xvals_sub*xy_factor,Y_OD);
xlabel('y ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity OD')
subplot(2,1,2)
hold on
plot(yvals_sub*xy_factor,X_OD);
xlabel('z ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity OD')

num = 1/sigma.*trapz(yvals_sub,trapz(xvals_sub,OD_actual,2))


% end

%%

