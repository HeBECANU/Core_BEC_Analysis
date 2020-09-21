function lims = getlims(data)
% A really quick function to return the [min,max] in each axis of data
% where the COLUMNS of data are the data coordinates
    naxis = size(data,2);
    lims = zeros(naxis,2);
    for axis = 1:naxis
       axmin = min(data(:,axis)); 
       axmax = max(data(:,axis)); 
       lims(axis,:) = [axmin,axmax];
    end    
end