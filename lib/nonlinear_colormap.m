function [cm_out]=nonlinear_colormap(cm_in,dyn_rng_cmp)

hsv=rgb2hsv(cm_in);
query_points=linspace(0,1,size(cm_in,1));
if isstring(dyn_rng_cmp) && stringcmp(dyn_rng_cmp,'log10')
     scaled_den=log10(query_points);
elseif isnumeric(dyn_rng_cmp) && dyn_rng_cmp~=0
     scaled_den=(query_points).^dyn_rng_cmp;
elseif (isstring(dyn_rng_cmp) && stringcmp(dyn_rng_cmp,'none')) || (isnumeric(dyn_rng_cmp) && dyn_rng_cmp==0)
    scaled_den=query_points;
end

cm_data=interp1(query_points,hsv,linspace(0,1,m));
cm_out=hsv2rgb(cm_data);


end
    
    