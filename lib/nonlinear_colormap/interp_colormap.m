function colors_out=interp_colormap(cm_in,query_points,range)


if nargin<3
    range=[0,1];
end


hsv=srgb_to_Jab(cm_in);
original_maping = linspace(range(1),range(2),size(cm_in,1));

cm_data=interp1(original_maping,hsv,query_points);

%%

colors_out=Jab_to_srgb(cm_data);

end
    
    