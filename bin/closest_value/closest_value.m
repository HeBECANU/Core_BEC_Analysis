function [v, idx] = closest_value(arr, val)
% Returns value and index of arr that is closest to val. If several entries
% are equally close, return the first. Works fine up to machine error (e.g.
% [v, i] = closest_value([4.8, 5], 4.9) will return [5, 2], since in float
% representation 4.9 is strictly closer to 5 than 4.8).
% ===============
% Parameter list:
% ===============
% arr : increasingly ordered array
% val : scalar in R

%BMH 20190507 com: this function needs a test script

%BMH 20190507 change: var inf to idx to prevent variable clash


len = length(arr);
idx = 1;
sup = len;

% Binary search for index
while sup - idx > 1
    med = floor((sup + idx)/2);
    
    % Replace >= here with > to obtain the last index instead of the first.
    if arr(med) >= val 
        sup = med;
    else
        idx = med;
    end
end

% Replace < here with <= to obtain the last index instead of the first.
if sup - idx == 1 && abs(arr(sup) - val) < abs(arr(idx) - val)
    idx = sup;
end  

v = arr(idx);
end
