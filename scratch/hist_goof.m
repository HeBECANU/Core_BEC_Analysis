cli_header('Testing begins')
clear all
num_bins = 3*[100,100,100]; %in T,X,Y
% bec_win = [0.5,0.54; %T
%             -.03,.03; %X
%             -.03,.03]; %Y

% TXT = data.qd.counts_txy{2};
% bec_win = [0.5,0.55;-.03,.03;-.03,.03];
TXY = [.1,2,1].*randn(1e5,3);
bec_win = 5*[-1,1;-1,1;-1,1];

% function hist_master(
T = TXY(:,1);
X = TXY(:,2);
Y = TXY(:,3);

TXY_counts = zeros(num_bins);
num_dims = size(TXY,2);


t_win = T > bec_win(1,1) & T < bec_win(1,2);
x_win = X > bec_win(2,1) & X < bec_win(2,2);
y_win = Y > bec_win(3,1) & Y < bec_win(3,2);

atoms = TXY(t_win& x_win & y_win,:);

t_edges = linspace(bec_win(1,1),bec_win(1,2),num_bins(1)+1)';
x_edges = linspace(bec_win(2,1),bec_win(2,2),num_bins(2)+1)';
y_edges = linspace(bec_win(3,1),bec_win(3,2),num_bins(3)+1)';

t_cents = 0.5*(t_edges(1:end-1)+t_edges(2:end));
x_cents = 0.5*(y_edges(1:end-1)+x_edges(2:end));
y_cents = 0.5*(x_edges(1:end-1)+y_edges(2:end));

%Need to modify to get T partition
T_sort = sort(T);
[t_counts,t_cuts]= mod_hist_adaptive_method(T_sort,t_edges,1); 
t_cuts = t_cuts(2:end-1);
t_counts = t_counts(2:end-1); %Collecting within hist bounds

[yt_slice_counts,y_cuts] = rec_hist(Y,t_cuts,y_edges);
[xt_slice_counts,x_cuts] = rec_hist(X,t_cuts,x_edges);
% voxel_counts = rec_hist3(this_shot,t_cuts,bec_win,num_bins);
% % Moving on lazily to get the layout sorted
X_sort = sort(X);
[xy_counts,xy_cuts] = mod_hist_adaptive_method(X_sort,x_edges,1); 
xy_counts = xy_counts(2:end-1);
xy_cuts = xy_cuts(2:end-1);
[xy_slice_counts,Y_cuts] = rec_hist(Y,xy_cuts,y_edges);

disp('done')

T_slice = sum(xt_slice_counts,2);
X_slice = sum(xt_slice_counts);
Y_slice = sum(xy_slice_counts);


sfigure(1);
clf
subplot(2,3,1)
imagesc(xt_slice_counts);
xlabel('X')
ylabel('T')

subplot(2,3,2)
plot(t_cents,T_slice,'.')
xlabel('T')

subplot(2,3,3)
imagesc(yt_slice_counts);
xlabel('Y')
ylabel('T')

subplot(2,3,4)
plot(x_cents,X_slice,'.')
xlabel('X')

subplot(2,3,5)
imagesc(xy_slice_counts);
xlabel('X')
ylabel('Y')

subplot(2,3,6)
plot(y_cents,Y_slice,'.')
xlabel('Y')

function [slice_counts,cuts] = rec_hist(data_in,lims,edges)
% A recursive algorithm for fast N-dimensional histogramming?
% Inputs: 
%         data_in     = NxD list of N data points in D dimension
%         lims        = 2xD list of range extrema in each dimension
%         dims        = 1xD list of number of bins in each dimension
    
    dims = size(edges)-1;
    dims = dims(dims>0);
%     num_bins = length(edges)-1;
    slice_counts = zeros(dims);
    cuts = zeros(dims);
    for idx = 1:dims(1)
           if idx == 1
               l_idx = 1;
               u_idx = lims(idx);
           end
           if lims(idx)>0
               l_idx = u_idx;
               u_idx = min(lims(idx),length(data_in));
           else
               continue
           end
           if u_idx>l_idx
                slice =  sort(data_in(l_idx+1:u_idx));
                [counts,cut_idxs]= mod_hist_adaptive_method(slice,edges,1); 
                counts = counts(2:end-1); %Collecting within hist bounds
                cuts(idx,:) = cut_idxs(2:end-1); % Be doubly sure this is right...
                slice_counts(idx,:) = counts;
           end
    end
end

function [slice_counts ,cut_idxs]= rec_hist3(data_in,lims,M,dims)
% A recursive algorithm for fast N-dimensional histogramming?
% Inputs: 
%         data_in     = NxD list of N data points in D dimension
%         lims        = list of bin edges taken by histogramming the sorted
                %         first input
%                 M     = extrema
%         dims        = 1xD list of number of bins in each dimension
% Each layer: Sort and histogram on first dimension
% ASSUME FIRST INPUT SORTED
% Return D-lvl dimensional histogram
    slice_counts = zeros(dims);
    if length(dims) == 1 %base case
        slice_counts = zeros(dims,1);
        edges = linspace(M(1),M(2),dims+1)';
        [counts,cut_idxs]= mod_hist_adaptive_method(data_in,edges,1); 
        slice_counts =  counts(2:end-1);
    else   
        for idx = 1:dims(1) % loop over slices of first index
           if idx == 1 %initialize for first slice
               l_idx = 1;
               u_idx = lims(idx);
           end
           if lims(idx)>0 %Later, update slices based on input limits
               l_idx = u_idx;
               u_idx = min(lims(idx),length(data_in));
           else
               continue
           end
           if u_idx>l_idx
               % Probably must be modified in higher dimensions... 
               slice =  sort(data_in(l_idx+1:u_idx)); % Pick the data in this slice
                edges = linspace(M(2,1),M(2,2),dims(2)+1)';
                [counts,cut_idxs]= mod_hist_adaptive_method(slice,edges,1); 
                slice_counts(idx,:) =  counts(2:end-1);
           end
        end
    end

end


function [bin_count,cut_idxs]=mod_hist_adaptive_method(x_dat,edges,is_x_sorted)
%adaptive_hist_method - a function that tries to chose the fastest histogram method
% uses a very simple case strucute. Future implmentations may use something like a SVM if it can be made fast
% enough. see scaling_tests to make a plot of the relative method speed and adjust the thresholds for your own
% computer.
%
% Syntax:         bin_counts=adaptive_hist_method(sorted_data,edges,1)
% Equivelent to:  bin_counts=histcounts(sorted_data,[-inf;edges;inf])
% Inputs:
%    x_dat           - column vector of data/counts, 
%    edges           - column vector of bin edges, MUST BE ORDERED!
%    is_x_sorted     - is the x_dat sorted, if you use issorted(x_dat) you will lose a lot of the potential
%    speedup

% Outputs:
%    bin_count - column vector, with length numel(edges)+1,  the first(last) element 
%                are the number of counts below(above) the first(last) edge
%
% Example: 
%     data=rand(1e5,1);
%     data=sort(data);
%     edges=linspace(0.1,1.1,1e6)';
%     out1=adaptive_hist_method(data,edges);
%     out2=histcounts(data,[-inf;edges;inf])';
%     isequal(out1,out2)
%
% Other m-files required: none
% Also See: scaling_tests,test_search_based_hist,adaptive_hist_method,compare_method_speeds
% Subfunctions: binary_search_first_elm
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%  - try some kind of ML eg SVM to predict which method to use.
%    - problem is the single inference of the svm can be up to 70ms which is a lot of overhead.
%  - a more complex (manual) inference procedure.
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2019-05-13

%------------- BEGIN CODE --------------

dat_size=numel(x_dat);
edges_size=numel(edges);

if is_x_sorted
    if edges_size<2e3
        [bin_count,cut_idxs]=mod_hist_count_search(x_dat,edges);
    else
        if dat_size>1e3 && log10(dat_size)>log10(edges_size)-1
            %call histcounts with edges
            bin_count=histcounts(x_dat,[-inf;edges;inf])';
        else
            bin_count=mod_hist_bin_search(x_dat,edges);
        end
    end
else
    if dat_size>1e3 && log10(dat_size)>log10(edges_size)-1
        %call histcounts with edges
        bin_count=histcounts(x_dat,[-inf;edges;inf])';
    else
         bin_count=mod_hist_bin_search(x_dat,edges);
    end
end



end


function [bin_count,cut_idxs]=mod_hist_count_search(data,edges)
%bin_search_hist - a histogram algorithm based on binary search of counts
% for each edge in the edge vector this code performs a binary search of
% the ordered data to find the count index of this edge
% !!!!!!!!!!!!!!! REQUIRES ORDERED DATA !!!!!!!!!!!!!!!!!!!!!!!!
% gives asymptotic speedup O(m·log(n)) over convertional histograming O(n·m) for dense histograms
% (many more counts than bins)
% Optimizations
%   - pre search for first(last) edge. Two inital searches for the hist limits eliminates
%     counts from search if they are not in the hisogram.
%   - moving search domain. Uses the previous edge as a lower limit to the
%     search domain.
%   - Simple sparse optimization. checks if there is any counts in the hist
%     bin before conducting search  (costs one compare but saves log(n)·unfilled bins)
%
% Syntax:         bin_counts=count_search_hist(data,edges)
% Equivelent to:  bin_counts=histcounts(data,[-inf;edges;inf])
% Designed to replicate histcounts(X,edges) "The value X(i)
%is in the kth bin if edges(k) ? X(i) < edges(k+1)" 
% Inputs:
%    data            - column vector of data/counts , MUST BE ORDERED!
%    edges           - column vector of bin edges, MUST BE ORDERED!

%
% Outputs:
%    bin_count - column vector, with length numel(edges)+1,  the first(last) element 
%                are the number of counts below(above) the first(last) edge
%



if ~iscolumn(data) || ~iscolumn(edges)
    error('inputs must be column vectors')
end

num_edges=size(edges,1);
num_bins=num_edges-1 +2;
num_counts=size(data,1);
%initalize output
bin_count=zeros(num_bins,1);
cut_idxs = zeros(num_bins,1); % the upper index bound for each bin 
%find the lowest edge
edge_lowest=edges(1);
idx_lowest=binary_search_first_elm(data,edge_lowest,1,num_counts);
val_lowest=data(idx_lowest);
if val_lowest<edge_lowest
    idx_lowest=idx_lowest+1;
end
%idx_lowest is now the first count in the first bin

% handle the case where the first edge is above the highest count
if idx_lowest==num_counts
    idx_lowest=idx_lowest+1;
end

%fprintf('lowest edge data idx %d\n',idx_lowest)
idx_u=idx_lowest;
idx_l=idx_lowest;
val_u=val_lowest; %for sparse opt
%the number of counts below the first edge is then 
bin_count(1)=idx_lowest-1; 
%check that there are still counts
rem_counts=num_counts-bin_count(1);

%do a binary search for the count index of the last edge to speed up the main loop (set search lims)
if rem_counts~=0
    edge_highest=edges(end);
    idx_max=binary_search_first_elm(data,edge_highest,idx_lowest,num_counts);
    
    %handle the case where the highest edge is below the first count
    if idx_max~=1
        idx_max=idx_max+1;
    end
    val_highest=data(idx_max);
    if val_highest<edge_highest
        idx_max=idx_max+1;
    end
    %idx_max is now the first count after the last bin
    %fprintf('hihgest edge data idx %d\n',idx_max)
    bin_count(end)=num_counts-idx_max+1;
    rem_counts=rem_counts-bin_count(1);
end

if rem_counts~=0
    ii=2;
    while idx_u<idx_max
        %fprintf('edge %d\n',ii)
        upper_bin_edge=edges(ii);
        %sparse opt, if the upper edge is smaller than the first count after the last bin
        if upper_bin_edge>val_u
            %fprintf('upper edge value %f\n',upper_bin_edge)
            %fprintf('upper edge count idx %d\n',idx_u)
            idx_u=binary_search_first_elm(data,upper_bin_edge,idx_u,idx_max);  
            val_u=data(idx_u);
            if val_u<upper_bin_edge 
                    idx_u=idx_u+1;
            end
            %idx_u is now the first count not in this bin
            %fprintf('upper edge count idx %d\n',idx_u)
            %fprintf('lower edge count idx %d\n',idx_l)
            if idx_u~=idx_l
                 bin_count(ii)=idx_u-idx_l;
                 cut_idxs(ii) = idx_u;
            end
            idx_l=idx_u;
        end
        ii=ii+1;
    end
end

end


function bin_count=mod_hist_bin_search(data,edges)
%bin_search_hist - a histogram algorithm based on binary search of bins
% for each count in the data vector this code performs a binary search of
% the edges to find the apropriate histogram bin to increment 
% gives asymptotic speedup O(n·log(m)) over convertional hitograming O(n·m) for sparse histograms
% (many more bins than counts)
%
% Syntax:         bin_counts=bin_search_hist(data,edges)
% Equivelent to:  bin_counts=histcounts(data,[-inf;edges;inf])
% Designed to replicate histcounts(X,edges) "The value X(i)
%is in the kth bin if edges(k) ? X(i) < edges(k+1)" 
% Inputs:
%    data            - column vector of data/counts , no ordering requirement
%    edges           - column vector of bin edges, MUST BE ORDERED!
%
% Outputs:
%    bin_count - column vector, with length numel(edges)+1,  the first(last) element 
%                are the number of counts below(above) the first(last) edge
% Example: 
%     data=rand(1e5,1);
%     data=sort(data);
%     edges=linspace(0.1,1.1,1e6)';
%     out1=bin_search_hist(data,edges);
%     out2=histcounts(data,[-inf;edges;inf])';
%     isequal(out1,out2)
% Other m-files required: none
% Also See: scaling_tests,test_search_based_hist,adaptive_hist_method,compare_method_speeds
% Subfunctions: binary_search_first_elm
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%  - try basic search reduction ,
%    - compare count with last value to search only edges above or below that.
%    - might give ~5% improvement, got me thinking about pre search lookup
%    tables 
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2019-05-13

%------------- BEGIN CODE --------------

if ~iscolumn(data) || ~iscolumn(edges)
    error('inputs must be column vectors')
end

% number of bins is edges-1 with 2 extra for below lowest and above highest
num_edges=size(edges,1);
num_bins=num_edges-1 +2;
bin_count=zeros(num_bins,1);
num_data=size(data,1);
for ii=1:num_data
    data_val=data(ii);
    closest_idx=binary_search_first_elm(edges,data_val,1,num_edges);
    closest_idx=closest_idx+1;
    %if closest is on the edge check if it should go up or down
     if closest_idx==2
         if data_val<edges(1)
             closest_idx=closest_idx-1;
         end
     elseif closest_idx==num_edges
         if data_val>edges(num_edges)
             closest_idx=closest_idx+1;
         end
     end
        
    bin_count(closest_idx)=bin_count(closest_idx)+1; 
end

end









%modified from mathworks submission by Benjamin Bernard 
%from https://au.mathworks.com/matlabcentral/fileexchange/37915-binary-search-for-closest-value-in-an-array
function idx_closest = binary_search_first_elm(vec, val,lower_idx,upper_idx)
% Returns index of vec that is closest to val, searching between min_idx start_idx . 
%If several entries
% are equally close, return the first. Works fine up to machine error (e.g.
% [v, i] = closest_value([4.8, 5], 4.9) will return [5, 2], since in float
% representation 4.9 is strictly closer to 5 than 4.8).
% ===============
% Parameter list:
% ===============
% arr : increasingly ordered array
% val : scalar in R
% use for debug in loop %fprintf('%i, %i, %i\n',btm,top,mid)

top = upper_idx(1);
btm = lower_idx(1);

% Binary search for index
while top > btm + 1
    mid = floor((top + btm)/2);
    % Replace >= here with > to obtain the last index instead of the first.
    if vec(mid) <= val %modified to work to suit histogram
        btm = mid;
    else
        top = mid;
    end
end

% Replace < here with <= to obtain the last index instead of the first.
%if top - btm == 1 && abs(arr(top) - val) < abs(arr(btm) - val)
%    btm = top;
%end  

idx_closest=btm;
end




