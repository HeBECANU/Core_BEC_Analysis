function data=hotspot_mask(data,anal_opts)
% this code deletes the hot spots from our detector using manualy defined region deletion
% this makes a flat-feild image look like swiss cheese
% code sniped from forbidden_heating_method
% TODO
% - automate the mask devlopment process
%   - 
%

if nargin<2 || isempty(anal_args)
    % default hot spot deletion
    % good as of 2019-11-20
    tlim=[0,inf];
    tmp_xlim=[-50e-3, 50e-3];    
    tmp_ylim=[-50e-3, 50e-3];
    anal_opts=[];
    anal_opts.hotspot_mask.square_mask=[tlim;tmp_xlim;tmp_ylim];
    anal_opts.hotspot_mask.circ_mask=[[0,0,38e-3,1];
                                      [32e-3,3.6e-3,3e-3,0];
                                      [35.6e-3,5.6e-3,4e-3,0];
                                      [25.36e-3,-19.8e-3,3e-3,0];
                                      [13.68e-3,31.6e-3,2e-3,0];
                                      [28.4e-3,-13.68e-3,2e-3,0];
                                      [32.08e-3,-12.24e-3,2e-3,0];
                                      [28.88e-3,-21.84e-3,2e-3,0];
                                      [18.8e-3,-28.72e-3,2e-3,0];
                                      [2.48e-3,-34.8e-3,2e-3,0];
                                  ];
end


num_shots=numel(data.mcp_tdc.counts_txy);
empty_shots=cellfun(@isempty,data.mcp_tdc.counts_txy);
data.mcp_tdc.masked.counts_txy={};
data.mcp_tdc.masked.num_counts=data.mcp_tdc.num_counts*nan;
fprintf('hotspot masking shots %04u:%04u',num_shots,0)
for ii=1:num_shots
    txy_shot=data.mcp_tdc.counts_txy{ii};
    if ~isempty(txy_shot)
        txy_shot=masktxy_square(txy_shot,anal_opts.hotspot_mask.square_mask);
        txy_shot=masktxy_2d_circle(txy_shot,anal_opts.hotspot_mask.circ_mask);
        data.mcp_tdc.masked.num_counts(ii)=numel(txy_shot);
        data.mcp_tdc.masked.counts_txy{ii}=txy_shot;
    else
        warning('empty shot')
    end
    if mod(ii,10)==0,fprintf('\b\b\b\b%04u',ii),end 
end
    
    
end