function vec_in=col_vec(vec_in)
if ~isvector(vec_in)
    error('input is not vector')
end
vec_in=reshape(vec_in,[],1);
end