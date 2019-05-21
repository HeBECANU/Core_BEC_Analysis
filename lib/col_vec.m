function in_vec=col_vec(in_vec)
% make the input vector a col vector
if ~isvector(in_vec)
    error('thats not a vector or scalar')
end
if isrow(in_vec)
    in_vec=transpose(in_vec);
end

end