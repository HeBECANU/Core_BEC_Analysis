function out = fast_hist3(data,varargin)
% A function that quickly produces 3D histograms from data.
% INPUTS
%    Essential:
%    xyz  = real-valued triples (X,Y,Z), sorted by X.
%    
%    Optional:
%    bins = triple of number of bins to use in X, Y, and Z.
%    OR
%    edges = cell array of demarcations of the bin edges, {X_edge,Y_edge,Z_edge}
% OUTPUTS
%    centres: PxQxR array of midpoints of each cell
%    values: PxQxR array of counts in each histogram bin
%    
% Needs fast_search_histogram https://github.com/brycehenson/fast_search_histogram
%    For now, assumes dense histogramming (n>>m)


% Verify first element is sorted
if ~issorted(data(:,1))
   [X,Xidx] = sort(data(:,1));
    Y = data(Xidx,2);
    Z = data(Xidx,3);
else
    X = data(:,1);
    Y = data(:,2);
    Z = data(:,3);
end


% Histogram the first dimension. This returns a list of count indices
% within each X-bin, which separates the (YZ) pairs for the next steps.
% Takes O(mx log(n))

% For each X-bin, sort the Y coordinate. This partitions the Z-values as
% above. Now each data point has an X- and a Y-value.
% Takes O(mx my log(n))

% Finally, sort and histogram in Z. 


   
% Cheeky method:
% Assuming even bin size (wx,wy,wz) with minimum values (mx,my,mz);
% (X,Y,Z) = ((x,y,z)-(mx,my,mz)) mod (wx,wy,wz)
% Loop over points:
% hist(X,Y,Z)++
   
   
end