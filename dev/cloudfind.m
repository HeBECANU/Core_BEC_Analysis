function out = cloudfind(tdc_data,varargin)
% A script that accepts a tdc data struct, produces a 3D histogram with
% show_txy_raw, and identifies distinct clouds by masking with respect to a
% threshold such that the number of connected components is equal to the
% median of all such possible bipartitions. 
% The points are then partitioned by (say) k-means clustering or a support
% vector model. 
% The COM, the point of peak density, and some scale parameters 
% of each pulse are calculated and returned as a struct

% This function would be useful in identifying PAL pulses without manually
% setting cloud separations, detecting clouds after a BEC Stern-Gerlach
% measurement, etc.
%  Could also try finding conn comp for a single amp cutoff - eg the median
%  value? But then misses many small peaks. Consider some measure from the
%  rank-value plot?

% Test:
% cloudfind(PAL_TDC)
% cloudfind(PAL_SG)


h_data = show_txy_raw(tdc_data,'draw_plots','false');
max_density = max(data.density_3d,'all');
% consider using fewer bins as mem scales like h*w*l 
h_data.density_3d = rescale(h_data.density_3d);
amplitude_vals = 1:10; %at 0, the entire image lattice is one connected component
lb_integral = zeros(size(amplitude_vals));
c = 0;
while c <= length(lb_integral)
   c = c + 1;
   Im = h_data.density_3d;
   pix_on = Im<av;
   lb_integral(av(c)) = sum(pix_on);
   
%    Calculate the connected components. Some ways:
%     Dimension of the null space of the adjacency matrix
    %     For SVM: What happens when k>num_blobs?
%         Many of the vectors will be close to each other
    
% Adajcency matrix method
%     All of the values above the threshold are considered 'on'. Two pixels
%     i and j are connected iff they are both on and adjacent.
%     The elements $(Im^n)_{i,j})$ of powers of Im give the number of paths 
%     of length n from i to j.
%     The exponential of the matrix is a weighted sum of path-connectedness
%     with factorially decaying weights.
%     Therefore, elements of the exponential exp(Im) are 0 iff two nodes do
%     not share any connection.
%     The number of zero eigenvalues of Im is equal to the number of
%     connected components. The null space of Im span (???) and have (???)
%     interpretation in the graph.

%     We can compute the adjacency matrix (Are there memory savings by
%     dynamic programming) by 
    % for each pix, subtract from neighbour (vert and horiz) and they are
    % connected iff sum(im(i,j),im(i,j+1)) == 2 for j<J-1
%      Ah what we probably want is the connectivity matrix, which has the
%      list of pairs? This will be cheaper. Is this what a sparse matrix
%      does?
%  Straightforward way is to
npixel = numel(Im);
% A = zeros(size(Im)); 
on_conn = false(numel(Im)); % Should be sparse?
%  Would be memory save to compute in 2d then in each dimension to find
% [X,Y] = meshgrid(1:size(pix_on,1),1:size(pix_on,1));
for i = 1:npixel
    X(i) = ind2sub(i,size(pix_on));
    for j = 1:npixel
        X(j) = ind2sub(j,size(pix_on));
        A(i,j) = pix_on(i,j)*(sum(X(i) - X(J) == 1) == 1); % up & nearest xy neigbours 
% Can also use find in diff (X-X').*pix_on
        
    IJ = reshape(Im,[numel(Im),1]); % Expensive - unnecssary?
    IJ = IJ - IJ'; % square matrix 
    
    lr_diff = Im + [Im(:,2:end),zeros(size(Im,1),1)];
    A = 

end





end

end