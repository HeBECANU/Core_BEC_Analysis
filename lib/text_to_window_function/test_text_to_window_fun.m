
%test if the text to windowing function works



%% fake dataset
%fake_data=fake_cl_corr(100,1);
%fake_data=fake_therm_cl_corr_small_hot(5000,1)
%fake_data=fake_therm_cl_corr_small_cold(5000,1);
%fake_data=fake_super_cl_corr(1000,1);
%fake_data=fake_cl_corr_medium_hot(5000,0.1);
data_faker=@() fake_super_cl_corr(5000,1);

fake_data=data_faker();

total_counts=sum(fake_data.num_counts);
%% Spliting data
fprintf('spliting data\n')
data_up_down={};
for ii=1:size(fake_data.counts_txy,2)
    shot_txy=fake_data.counts_txy{ii};
    mask=rand(size(shot_txy,1),1)>0.5;
    data_up_down.up_counts{ii}=shot_txy(mask,:);
    data_up_down.down_counts{ii}=shot_txy(~mask,:);
end

%% make windowing text function
fprintf('calculating text mask\n')
windowing_function=text_to_window_fun([],0);
%%
fprintf('applying text mask\n')
for ii=1:size(data_up_down.up_counts,2)
    down_counts=data_up_down.down_counts{ii};
    mask=windowing_function(down_counts);
    %simple mask
    %mask=(down_counts(:,2)<0 & down_counts(:,3)<0) | (down_counts(:,2)>0 & down_counts(:,3)>0 & down_counts(:,3)<0.5)  ;
    data_up_down.down_counts_mask{ii}=down_counts(mask,:);
end


%% Verify the mask has been applied
x_edges=linspace(-1,1,50);
y_edges=linspace(-1,1,50);

all_mask_counts=vertcat(data_up_down.down_counts_mask{:});
bins=histcn(all_mask_counts(:,2:3),x_edges,y_edges);
sfigure(1)
clf
imagesc(bins')
set(gca,'YDir','normal')
title('Mask port all counts')
colormap(viridis)
pause(1e-3)