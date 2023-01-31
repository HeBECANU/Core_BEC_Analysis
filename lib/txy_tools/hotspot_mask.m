function [data,anal_opts]=hotspot_mask(data,anal_opts)
% this code deletes the hot spots from our detector using manualy defined region deletion
% this makes a flat-feild image look like swiss cheese !
% code sniped from forbidden_heating_method
% TODO
% - automate the mask devlopment process
%   - 

% gracefully handle if the input is a structure
if ~isfield(data,'mcp_tdc')
    out_data = data;
else
    out_data = data.mcp_tdc;
end

if nargin<2 || isempty(anal_args)
    % default hot spot deletion
    % good as of 2019-11-20
    tlim=[-inf,inf]; % somtimes we shift t=0 so this should accept all times
    tmp_xlim=[-50e-3, 50e-3];    
    tmp_ylim=[-50e-3, 50e-3];
    anal_opts.hotspot_mask=[];
    anal_opts.hotspot_mask.square_mask=[tlim;tmp_xlim;tmp_ylim];
    % mask updated 2020-06-14
    % using some tune out data
    % countrate with hotspot mask on 30.7 Hz
    % countrate without hotspot mask on 474 Hz
    % pass fraction relative to 40 mm circle 0.8208 , or 0.9095 rel to 38 mm

%     anal_opts.hotspot_mask.circ_mask=[[0,-0.0e-3,38e-3,1];
%                                    % start from north then go clockwise
%                                   [5.788e-3,35.61e-3,4e-3,0];
%                                   [13.99e-3,31.79e-3,2.5e-3,0];
%                                   [22.69e-3,28.41e-3,3e-3,0];
%                                   [29.82e-3,20.82e-3,3e-3,0];
%                                    
%                                   [30.5e-3,5.0e-3  ,4.5e-3,0];
%                                   [35.6e-3,5.6e-3,4.5e-3,0];
%                                   
%                                   [30.5e-3,-15.71e-3,6e-3,0];                                   
%                                   [25.36e-3,-19.8e-3,3e-3,0];
%                                   [27.79e-3,-23.24e-3,3.5e-3,0];
%                                   
%                                   [23.81e-3,-26.49e-3,3e-3,0];
%                                   
%                                   [19.16e-3,-29.46e-3,4e-3,0];
%                                   
%                                   [2.48e-3,-34.8e-3,2e-3,0];
%                                   
%                                   [26.96e-3,21.46e-3,1e-3,0];
%                                     ];

   anal_opts.hotspot_mask.circ_mask=[[3.44e-3,18.96e-3,5.5e-3,0]];
                                
%                                  [28.4e-3,-13.68e-3,2e-3,0];
%                                 [32.08e-3,-12.24e-3,2e-3,0];
%                                 [29.54e-3,-17.09e-3,1.5e-3,0];              

%     anal_opts.hotspot_mask.circ_mask=[[0,0,38e-3,1];
%                                       [32e-3,3.6e-3,3e-3,0];
%                                       [35.6e-3,5.6e-3,4e-3,0];
%                                       [25.36e-3,-19.8e-3,3e-3,0];
%                                       [13.68e-3,31.6e-3,2e-3,0];
%                                       [28.4e-3,-13.68e-3,2e-3,0];
%                                       [32.08e-3,-12.24e-3,2e-3,0];
%                                       [28.88e-3,-21.84e-3,2e-3,0];
%                                       [18.8e-3,-28.72e-3,2e-3,0];
%                                       [2.48e-3,-34.8e-3,2e-3,0];
%                                   ];

end


num_shots=numel(out_data.counts_txy);
empty_shots=cellfun(@isempty,out_data.counts_txy);
empty_shot_indx = [];
out_data.masked.counts_txy={};
out_data.masked.num_counts=nan(size(out_data.counts_txy));
fprintf('hotspot masking shots %04u:%04u',num_shots,0)
for ii=1:num_shots
    txy_shot=out_data.counts_txy{ii};
    if ~isempty(txy_shot)
        txy_shot=masktxy_square(txy_shot,anal_opts.hotspot_mask.square_mask);
        txy_shot=masktxy_2d_circle(txy_shot,anal_opts.hotspot_mask.circ_mask);
        out_data.masked.num_counts(ii)=size(txy_shot,1);
        out_data.masked.counts_txy{ii}=txy_shot;
    else
        out_data.masked.num_counts(ii)=numel(txy_shot);
        out_data.masked.counts_txy{ii}=txy_shot;
        empty_shot_indx = [empty_shot_indx ii]; %warning('empty shot')
    end
    if mod(ii,10)==0,fprintf('\b\b\b\b%04u',ii),end 
end
if ~isempty(empty_shot_indx)
    fprintf('\nWarning: empty shot(s) %g \n',empty_shot_indx)
end
if ~isfield(data,'mcp_tdc')
    data = out_data;
else
    data.mcp_tdc = out_data;
end
    
end