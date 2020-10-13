function v_out = txy2vel(txy_in,varargin)
% A function that converts txy data to velocity-space.
% Takes either a standard tdc struct and appends a field vel_zxy
% OR accepts an Nx3 array of txy data and returns one of the same size.
% Required inputs:
% tdc_struct - A struct with field (at least) counts_txy, an Nx1 cell of n_i x 3 arrays
% Optional inputs:
% COM - centre of mass of cloud. Good idea to specify if pulses are saturated.
% defaults to the mean position of all txy points (shotwise)
% gravity - defaults to hebec_constant value
    c = hebec_constants();
    p = inputParser;
    addParameter(p,'gravity',c.g0);
    addParameter(p,'fall_distance',c.fall_distance);
    addParameter(p,'release_time',0);
    addParameter(p,'t_COM',nan);
    parse(p,varargin{:});
    g = p.Results.gravity;
    s = p.Results.fall_distance;
    t0 = p.Results.release_time;
    t_COM = p.Results.t_COM;
    
%     if isstruct(txy_in)
%         mode = 'struct';
%     elseif isvector(txy_in)
%         mode = 'time';
%     elseif ismatrix(txy_in)
%         mode = 'array';
%     else
%         error('TXY data type not recognized');
%     end
%     if strcmP(mode,'array')
    if isnan(t_COM) %if nothing passed
        t_COM = mean(txy_in(:,1));
    else
        warning('t_COM argument not understood, using mean arrival time')
        t_COM = mean(txy_in(:,1));
    end
%     end
    
    fall_times = txy_in(:,1) - t0;
    vz = -g*(t_COM^2-fall_times.^2)./(2*fall_times);
    % Correct to within machine precision
    v_out = [vz, txy_in(:,2:3)./fall_times];

% function vzxy_out=txy_to_vel(txy_in,out_time,gravity,fall_distance)

end

