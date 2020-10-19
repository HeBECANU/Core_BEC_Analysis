function [val,idx] = sat_pulse_centre(T,Y,thr)
    % A function to find the centre of a saturated peak.
    % Works by capping the flux at some value below the saturated level, 
    % then finding the centre of mass of the thresholded density.
    % Inputs: Equally-sized vectors of X and Y (density values) and a
    % cutoff thr
    m = rescale(Y)>thr;
    Y(m) = thr*max(Y);
    val = col_vec(mean(T.*Y)/sum(Y));
    [~,idx] = closest_value(T,val);
    idx = col_vec(idx);
end