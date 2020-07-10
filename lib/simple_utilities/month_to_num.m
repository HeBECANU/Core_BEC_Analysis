function num=month_to_num(str_in)
% convert a string of the month to a number of the month
% allows partial matches
% can handle arb size cell inputs
% examples
% month_to_num('jan')
% month_to_num('NOV')
% month_to_num('dEc')
% month_to_num({'dEc','nov' ; 'apr', 'june'})
% also see
% test_month_to_num

    month_list = {'January','February','March','April','May','June','July','August','September','October','November','December'}';
    month_list=lower(month_list);
    % if a cell call this function with each element of the cell
    if iscell(str_in)
        num=nan(size(str_in));
        for ii=1:numel(num)
            num(ii)=month_to_num(str_in{ii});
        end
    else
        str_in=lower(str_in);
        if numel(str_in)<3
            error('string must be 3 chars or more')
        end
        num=find(cellfun(@(x)contains(x,str_in), month_list));
        if isempty(num) || numel(num)>1
            num=nan;
        end
    end
end