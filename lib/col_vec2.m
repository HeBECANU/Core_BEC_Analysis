function in_vec=col_vec2(in_vec)
% make the input vector a col vector
if ~isvector(in_vec)
    error('thats not a vector or scalar')
end
if isrow(in_vec)
    if size(in_vec,2)<1e5
        in_vec=in_vec(:);
    else
        in_vec=transpose(in_vec);
    end
end

end