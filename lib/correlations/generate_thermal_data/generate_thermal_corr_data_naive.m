
% THIS ALGO DOES NOT WORK!
% the fairly direct approach to the conditional probability ends up narrowing up with a width about the correlation
% length and does not produce a thermal density distribution

% %% start with some seed data (does not work)
% num_seed=1e2;
% thermal_data.xcord=normrnd(0,g1_sigma,num_seed,1);
% thermal_data.therm_amp=g1(thermal_data.xcord,g1_sigma);


uniform_region=5;
g1_sigma=1;
corr_sigma=0.2;
peak_value=1/(g1_sigma*sqrt(2*pi));
headroom=1.1;
% initalize data structure;
thermal_data=[];
thermal_data.xcord=[];
thermal_data.therm_amp=[];

plot_handle=stfig('conditional prop plot');
clf

%% plot the thermal pdf

stfig('conditional prop plot');
xsamp=linspace(-5,5,1e4);
plot(xsamp,g1(xsamp,g1_sigma))


%%
while true

trying_to_sample=true;
tries=1;
while trying_to_sample 
    sample_x_cord=uniform_region*(rand(1)-0.5)*2;
    sample_amp=rand(1)*headroom*peak_value;
    pdf_amp=conditional_amp(sample_x_cord,thermal_data,g1_sigma,corr_sigma);
    if sample_amp<pdf_amp
        trying_to_sample=false;
    end
    tries=tries+1;
end
%fprintf('tries %u\n',tries)
thermal_data.xcord=cat(1,thermal_data.xcord,sample_x_cord);
thermal_data.therm_amp=cat(1,thermal_data.therm_amp,pdf_amp);

% find the conditional probability

xsamp=col_vec(linspace(-uniform_region,uniform_region,round(uniform_region*4/corr_sigma)));
ysamp=conditional_amp(xsamp,thermal_data,g1_sigma,corr_sigma);
peak_value=max(ysamp);
stfig(plot_handle);
subplot(2,1,1)
plot(xsamp,ysamp)
hold on
plot(thermal_data.xcord,conditional_amp(thermal_data.xcord,thermal_data,g1_sigma,corr_sigma),'x')
hold off
subplot(2,1,2)
shout=smooth_hist(thermal_data.xcord,'sigma',0.01,'lim',[-1,1]*uniform_region);
plot(shout.bin.centers,shout.count_rate.smooth/sum(shout.count_rate.smooth))
xlim([-1,1]*uniform_region)

pause(0.1)

end

%% now we find the max value of the prob dist which is easy bc there is only one particle
%peak_value=conditional_amp(thermal_data.xcord,thermal_data,g1_sigma,corr_sigma);




%%
conditional_amp(thermal_data.xcord,thermal_data,g1_sigma,corr_sigma)
%%

function amp_out=g1(x_in,g1_sigma)
   amp_out=normpdf(x_in,0,g1_sigma);
end

function amp_out=gn_amp(x_in,n,corr_sigma)
   amp_out=unit_amp_gauss(x_in,0,corr_sigma*sqrt(n));
end

function y=unit_amp_gauss(x,mu,sigma)
    y = exp(-0.5 * ((x - mu)./sigma).^2);
end

function amp_out=conditional_amp(x_in,thermal_data,g1_sigma,corr_sigma)
%
%%amp_corr=min(cat(2,2*thermal_amp,sum(amp_corr,2)),[],2);

thermal_amp=g1(x_in,g1_sigma);
if numel(thermal_data.xcord)>0
    delta_mat = bsxfun(@minus, x_in, thermal_data.xcord');
    % if numel(thermal_data.xcord)<2
    %     thermal_data.xcord=thermal_data.xcord(2:end);
    % else
    %     delta_mat=triu(delta_mat,1);
    % end
    amp_corr=gn_amp(delta_mat,2,corr_sigma); %xin col, xcounts rows
    amp_corr=amp_corr.*repmat(thermal_data.therm_amp',[size(x_in,1),1]);
    amp_corr=sum(amp_corr,2);
else
    amp_corr=0;
end
amp_out=amp_corr+thermal_amp;
end