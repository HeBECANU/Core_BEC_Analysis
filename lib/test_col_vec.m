% test col_vec2






vec_size=1e4;

%test with col input
in_vec=rand(1,vec_size)';
tic
out_vec=col_vec2(in_vec);
time_col_vec2_col=toc;
if ~iscolumn(out_vec) 
    error('not col vec')
end

in_vec=rand(1,vec_size)';
tic
out_vec=in_vec(:);
time_inbuilt_col=toc;
if ~iscolumn(out_vec) 
    error('not col vec')
end

%test with row input

in_vec=rand(1,vec_size);
tic
out_vec=col_vec2(in_vec);
time_col_vec2_row=toc;
if ~iscolumn(out_vec) 
    error('not col vec')
end

in_vec=rand(1,vec_size);
tic
out_vec=col_vec(in_vec);
time_col_vec_row=toc;
if ~iscolumn(out_vec) 
    error('not col vec')
end

in_vec=rand(1,vec_size);
tic
out_vec=in_vec(:);
time_inbuilt_row=toc;
if ~iscolumn(out_vec) 
    error('not col vec')
end
fprintf('vector size 10^%u\n',log10(vec_size))
fprintf('speedup for col vec %f \n',time_inbuilt_col/time_col_vec2_col)
fprintf('speedup for row vec %f \n',time_inbuilt_row/time_col_vec2_row)
fprintf('speedup for row vec %f \n',time_inbuilt_row/time_col_vec_row)


%% timing plot
sizes




%%
col_vec2(magic(5))


%% do timing tests

