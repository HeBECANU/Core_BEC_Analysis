

uniform_region=10;
g1_sigma=0.4;
corr_sigma=0.01;
headroom=1.2;
peak_value=headroom*1/(g1_sigma*sqrt(2*pi));
% initalize data structure;
thermal_data=[];
thermal_data.xcord=[];
%thermal_data.therm_amp=[];
thermal_data.xcord=0;
plot_handle=stfig('conditional prop plot');
max_order=4;
clf

%% plot the thermal pdf

stfig('conditional prop plot');
xsamp=linspace(-5,5,1e3);
plot(xsamp,g1(xsamp,g1_sigma))

%%
fprintf('\ncounts %05u',numel(thermal_data.xcord))
while true
trying_to_sample=true;
tries=1;
while trying_to_sample 
    sample_x_cord=uniform_region*(rand(1)-0.5)*2;
    sample_amp=rand(1)*peak_value;
    pdf_amp=conditional_amp(sample_x_cord,thermal_data,g1_sigma,corr_sigma,max_order);
    if sample_amp<pdf_amp
        trying_to_sample=false;
    end
    tries=tries+1;
end

thermal_data.xcord=cat(1,thermal_data.xcord,sample_x_cord);
%thermal_data.therm_amp=cat(1,thermal_data.therm_amp,pdf_amp);

% find the conditional probability

xsamp=col_vec(linspace(-uniform_region,uniform_region,round(uniform_region*1/corr_sigma)));
ysamp=conditional_amp(xsamp,thermal_data,g1_sigma,corr_sigma,max_order);
peak_value=max(ysamp)*headroom;
if mod(numel(thermal_data.xcord),10)==0 || numel(thermal_data.xcord)==1
    stfig(plot_handle);
    subplot(2,1,1)
    plot(xsamp,ysamp)
    hold on
    % plot the contributions from the multiple correlation amplitudes
%     ysamp=conditional_amp(xsamp,thermal_data,g1_sigma,corr_sigma,1);
%     plot(xsamp,ysamp)
%     ysamp=conditional_amp(xsamp,thermal_data,g1_sigma,corr_sigma,2);
%     plot(xsamp,ysamp)
%     ysamp=conditional_amp(xsamp,thermal_data,g1_sigma,corr_sigma,3);
%     plot(xsamp,ysamp)
    plot(thermal_data.xcord,conditional_amp(thermal_data.xcord,thermal_data,g1_sigma,corr_sigma,max_order),'x')
    plot(xsamp,peak_value*unit_amp_gauss(xsamp,0,g1_sigma),'-')
    line([min(xsamp),max(xsamp)],[1,1]*peak_value)
    hold off
    subplot(2,1,2)
    shout=smooth_hist(thermal_data.xcord,'sigma',0.5,'lim',[-1,1]*uniform_region);
    plot(shout.bin.centers,shout.count_rate.smooth/(diff(shout.bin.centers(1:2))*sum(shout.count_rate.smooth)),'b')
    hold on
    shout=smooth_hist(thermal_data.xcord,'sigma',1/sqrt(numel(thermal_data.xcord)),'lim',[-1,1]*uniform_region);
    plot(shout.bin.centers,shout.count_rate.smooth/(diff(shout.bin.centers(1:2))*sum(shout.count_rate.smooth)),'k')
    txt = {sprintf('cen=%.4f',mean(thermal_data.xcord)),...
        sprintf('std=%.4f',std(thermal_data.xcord))};
    text(0.8,0.5,txt,'FontSize',20,'Units','normalized')
    hold off
    xlim([-1,1]*uniform_region)
    pause(1e-6)
    fprintf('\b\b\b\b\b%05u',numel(thermal_data.xcord))
end


end
fprintf('\n')

%% now we find the max value of the prob dist which is easy bc there is only one particle
%peak_value=conditional_amp(thermal_data.xcord,thermal_data,g1_sigma,corr_sigma);




%%
conditional_amp(thermal_data.xcord,thermal_data,g1_sigma,corr_sigma)
%%

function amp_out=g1(x_in,g1_sigma)
   amp_out=normpdf(x_in,0,g1_sigma);
end

function amp_out=big_g2_uncorr(dx_in,n,g1_sigma)
 amp_out=normpdf(dx_in,0,sqrt(n)*g1_sigma);
end

function amp_out=gn_amp(delt_x_in,n,corr_sigma)
   amp_out=1+unit_amp_gauss(delt_x_in,0,corr_sigma*sqrt(n-1))*(factorial(n)-1);
end

function y=unit_amp_gauss(x,mu,sigma)
    y = exp(-0.5 * ((x - mu)./sigma).^2);
end

function amp_out=conditional_amp(x_in,thermal_data,g1_sigma,corr_sigma,ordermax)
amp_corr=zeros(size(x_in,1),1);
amp_corr=amp_corr+g1(x_in,g1_sigma);
if numel(thermal_data.xcord)>0 && ordermax>1
    delta_mat = bsxfun(@minus, x_in, thermal_data.xcord');
    % if numel(thermal_data.xcord)<2
    %     thermal_data.xcord=thermal_data.xcord(2:end);
    % else
    %     delta_mat=triu(delta_mat,1);
    % end
    for order=2:ordermax
        amp_corr_tmp=gn_amp(delta_mat,order,corr_sigma).*big_g2_uncorr(delta_mat,order,g1_sigma); %xin col, xcounts rows
        amp_corr=amp_corr+sum(amp_corr_tmp,2);
    end
else
    amp_corr=0;
end
amp_out=amp_corr;
end
