function s=rename_struct_field(s,old_field,new_field)
[s.(new_field)] = s.(old_field);
s = rmfield(s,old_field);
end