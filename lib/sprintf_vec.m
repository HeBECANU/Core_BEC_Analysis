function string=sprintf_vec(format_single_in,mat_in,delim)

if nargin<3 || isempty(delim)
    delim=',';
end
% format an array using 
string_elem=arrayfun(@(x) sprintf(format_single_in,x),mat_in,'UniformOutput',false);
string=strjoin(string_elem,delim);

end