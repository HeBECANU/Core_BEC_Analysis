%test_smooth_hist

num_counts=1000;

xdata=cat(1,col_vec(rand(1,num_counts)),col_vec(randn(1,num_counts)));

%% old way
bin_lims=[-3,3];
bin_width=1e-1;


x_bin_num=round(range(bin_lims)/bin_width);
x_bin_num=2*floor( x_bin_num/2)+1; %round to an odd number
x_bin_width=range(bin_lims)/x_bin_num; %change the value to the rounded one
%set(handles.status_handle ,'String',['No Counts in 2d Window']);
%set(handles.time_binsize,'String',num2str(T_bin_width)); %upate the
%used odd bins time window
%change the max time so that whole bins can fit
x_bin_edge=linspace(min(bin_lims),min(bin_lims)+x_bin_width*x_bin_num,x_bin_num);
[T1d_counts,edges]=histcounts(xdata,x_bin_edge);
t_centers=mean([edges(1:end-1);edges(2:end)]);
T1d_counts=T1d_counts/(bin_width);   

stfig('old way of making a histogram')
plot(t_centers,T1d_counts)



%% new way with smoth_hist 

output=smooth_hist(xdata,'lims',bin_lims,'sigma',bin_width)

stfig('new way of making a histogram')
plot( output.bin.centers,output.count_rate.smooth)

%% handle single input
smooth_hist(xdata,'sigma',bin_width)
stfig('new way of making a histogram')
plot( output.bin.centers,output.count_rate.smooth)

%% handle no inputs
smooth_hist(xdata)
stfig('new way of making a histogram')
plot( output.bin.centers,output.count_rate.smooth)
