function nandata = nandata(data,window)
% An inline helper function to replace data(window) with nans
%```
% data = [1,-5,2,-12];
% window = data<0;
% nan_data = nandata(data,window);
% all(isnan(nan_data(window)))&all(nan_data(~window)==data(~window))
%```
    nandata = data;
    nandata(window) = nan;
end