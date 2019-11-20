function [cm_out]=nonlinear_colormap(cm_in,type,params,do_scaling_plot)

% logistic
%   -  peak gradient
%   -  x pos of steepest part of curve
%   -  logistic height
%   -  logistic offset
%   -  add a linear to make the total height (logistic + linear)=total height, 0 for off

%power
%   - power
%   - scaling
%   - offset



if ~(isstring(type) || ischar(type))
    error('type is not valid')
end




hsv=srgb_to_Jab(cm_in);
original_maping = linspace(0,1,size(cm_in,1));
if strcmp(type,'logistic')
    if numel(params)<1
        error('need params')
    end
    default_params=[10,0.5,1,0,0];
    if sum(~isnan(params))<sum(~isnan(default_params))
        tmp_params=nan*default_params;
        tmp_params(~isnan(params))=params(~isnan(params));
        params=tmp_params;
        params(isnan(params))=default_params(isnan(params));
    end
    logistic_k=params(1);
    logistic_x0=params(2);
    logistic_height=params(3);
    logistic_offset=params(4);
    logistic_linear=params(5);
    query_points=logistic_offset+ logistic_height./(1+exp(-logistic_k*(original_maping-logistic_x0)));
    if logistic_linear~=0
        query_points=query_points+original_maping*(logistic_linear-logistic_height);
    end
elseif strcmp(type,'power')
    default_params=[0.5,1,0];
    if sum(~isnan(params))<sum(~isnan(default_params))
        tmp_params=nan*default_params;
        tmp_params(~isnan(params))=params(~isnan(params));
        params=tmp_params;
        params(isnan(params))=default_params(isnan(params));
    end
    pow_pow=params(1);
    pow_scaling=params(2);
    pow_offset=params(3);
    
    query_points=pow_offset+pow_scaling*(original_maping).^pow_pow;
elseif strcmp(type,'linear')
   
else     
    original_maping=original_maping;
end

cm_data=interp1(original_maping,hsv,query_points);

%%
if nargin>3 && do_scaling_plot==1
    figure(3)
    title('colormap rescaling')
    plot(original_maping,query_points)
    xlabel('input')
    ylabel('output')
    xlim([0,1])
    ylim([0,1])
end

%%


cm_out=Jab_to_srgb(cm_data);



end
    
    