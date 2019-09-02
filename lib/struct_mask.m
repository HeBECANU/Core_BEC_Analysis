function out = struct_mask(s,m)
% Assuming s is a structure with equal sizes in each field, mask the
% fields according to m
fieldnames = fields(s);
nfields = numel(fieldnames);
out = [];
for ii=1:nfields
   this_field = fieldnames{ii};
   this_data = s.(this_field);
   if size(this_data,1) == 1
       this_data = this_data';
   end
   if size(this_data,2) > 1
       out.(this_field) =this_data(m,:);
   else
       out.(this_field) =this_data(m);
   end
end

end