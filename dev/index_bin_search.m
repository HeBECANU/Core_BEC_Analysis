function bin_idx=index_bin_search(data,edges)
if ~iscolumn(data) || ~iscolumn(edges)
    error('inputs must be column vectors')
end

% number of bins is edges-1 +  2 extra for below lowest and above highest
num_edges=size(edges,1);
num_data=size(data,1);
bin_idx=zeros(num_data,1);
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
	bin_idx(ii)=closest_idx; 
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