function colors_out=make_simple_colormap(colors_in,values,out_pts)



hsv=srgb_to_Jab(colors_in);
query_points = linspace(0,1,out_pts);

cm_data=interp1(values,hsv,query_points);

%%

colors_out=Jab_to_srgb(cm_data);

end
    
    