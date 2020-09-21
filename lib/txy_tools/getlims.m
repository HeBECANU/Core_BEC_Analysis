function lims = getlims(data)
% getlims(x) returns the [min,max] along the columns of x
% Motivation: Where x is an [N x D] array of N samples in D dimensions,
% return the bounds of the interval where the data points are found
% eg
% '''
% d1 = [1 0 0];
% d2 = [0,-1,3];
% d3 = [2,-.5,1];
% x = [d1;d2;d3];
% isequal(getlims(x),[0,2;-1,0;0,3]);
% ```
    naxis = size(data,2);
    lims = zeros(naxis,2);
    for axis = 1:naxis
       axmin = min(data(:,axis)); 
       axmax = max(data(:,axis)); 
       lims(axis,:) = [axmin,axmax];
    end    
end