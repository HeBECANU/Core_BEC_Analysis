function zxy_dat = centre_BEC_COM(configs,input_txy)
%%%%
%Finding BEC centres. bec_boxlim (in txy) set in config file, encloses condensate.
%Boxcull then trims the data down to include only these points, and the for-loop below
%converts them to zxy. For each bec, the centre is stored in bec_cent, and
%then the average is taken. txy_0 is the offset-corrected condensate, where
%the condensate is centred on the mean position.
%%%%

    %set constants and intialize arrays
    tof=configs.const.tof;
    vz=9.8*tof; 
    nShot = numel(input_txy);
    txy_dat = cell(1,numel(nShot));
    zxy_dat = cell(1,numel(nShot));
    bec_cent=zeros(nShot,3);    %initialize arrays to store condensate centres
    num_txy_raw=zeros(1,nShot); %and for number of counts in loaded shot


    % Build BEC locator box
    bec_boxlim=cell(1,3);
    for i=1:3
        bec_boxlim{i}=configs.bec.txy_pos(i)+configs.bec.box_fwidth(i)*[-0.5,0.5];
    end
    
    %find the mean posn of each BEC and centre them, then convert to zxy
    for i=1:nShot
        [txy_dat{i},~,~,bec_cent(i,:)]=boxcull(input_txy{i},bec_boxlim);
        num_txy_crop(i)=size(txy_dat{i},1);    % number of counts in this box
        txy_dat{i}=txy_dat{i}-repmat(bec_cent(i,:),[num_txy_crop(i),1]);   % centre around self-average
            % T-Z conversion
        zxy_dat{i}=txy_dat{i};
        zxy_dat{i}(:,1)=zxy_dat{i}(:,1)*vz;
    end
end

nShot = numel(input_txy);
txy_dat = cell(1,numel(nShot));
zxy_dat = cell(1,numel(nShot));
bec_cent=zeros(nShot,3);    %initialize arrays to store condensate centres
num_txy_raw=zeros(1,nShot); %and for number of counts in loaded shot

for i=1:nShot
    [~,~,~,bec_cent(i,:)]=boxcull(txy_raw{i},bec_boxlim);
    num_txy_raw(i)=size(txy_raw{i},1);    % number of counts in this shot
    txy_dat{i}=txy_raw{i}-repmat(bec_cent(i,:),[num_txy_raw(i),1]);   % centre around self-average
        % T-Z conversion
    zxy_dat{i}=txy_dat{i};
    zxy_dat{i}(:,1)=zxy_dat{i}(:,1)*vz;
end