function out = struct_mask(s,m,maintain_orr)
% Assuming s is a structure with equal sizes in each field, mask the
% fields according to mask m
% Improvements
% More robust testing for multidimensional arrays
% Be more strict about row/column input (perhaps faster)
% And definitely safer if a square matrix is ever input

fieldnames = fields(s);
nfields = numel(fieldnames);
out = [];
% flip=0;
% num_entries = length(m);
for ii=1:nfields
    this_field = fieldnames{ii};
    this_data = s.(this_field);
    if size(this_data,1) == 1
        this_data = this_data';
        if nargin>2 && maintain_orr
            if size(this_data,2) > 1
                out.(this_field) =this_data(m,:)';
            else
                out.(this_field) =this_data(m)';
            end
        else
            if size(this_data,2) > 1
                out.(this_field) =this_data(m,:);
            else
                out.(this_field) =this_data(m);
            end
        end
    else
        if size(this_data,2) > 1
            out.(this_field) =this_data(m,:);
        else
            out.(this_field) =this_data(m);
        end
    end
end

end



% Tests to run:
%     Variable table sizes