function [val,idx] = sat_pulse_centre(T,Y,minval,upperval)
    %[val,idx] = sat_pulse_centre(T,Y,minval,thr)    
    % A function to find the centre of a saturated peak.
    % Works by capping the flux at some value below the saturated level, 
    % then finding the centre of mass of the thresholded density.
    
    Y = col_vec(Y);
    T = col_vec(T);
    m = rescale(Y)>upperval;
    Y(m) = upperval*max(Y);
    Y(rescale(Y)<minval) = 0;
    val = col_vec(sum(T.*Y)/sum(Y));
    [~,idx] = closest_value(T,val);
    idx = col_vec(idx);
end