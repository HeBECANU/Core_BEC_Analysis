function rgb_out=array_to_rgb(array,colormap,range)

    if isempty(range)
        range=[0,max(counts(:))*1.1];
    end

    num_colors=2^24-1;
    % now we scale array such that range(0) is mapped to zero and range(1)
    % maps to unit32max
    % following https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
    % scale the values in array in range to between 0 and 1
    norm_counts=(array-range(1))/(range(2)-range(1));
    norm_counts=uint32(norm_counts*num_colors);
    %norm_counts=transpose(norm_counts);
    %norm_counts=flipud(norm_counts);
    cmap_num_colors=interp_colormap(colormap,linspace(0,1,num_colors));
    rgb_out= ind2rgb(norm_counts, cmap_num_colors);
  

end