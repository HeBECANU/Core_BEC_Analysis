function out = struct_mask(s,m)
% Assuming s is a structure with equal sizes in each field, mask the
% fields according to mask m
% Improvements
    % More robust testing for multidimensional arrays
    % Be more strict about row/column input (perhaps faster)
        % And definitely safer if a square matrix is ever input
        
fieldnames = fields(s);
nfields = numel(fieldnames);
out = [];
for ii=1:nfields
   this_field = fieldnames{ii};
   this_data = getfield(s,this_field);
   out = setfield(out,this_field,this_data(m));
end

end