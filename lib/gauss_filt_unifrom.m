function out_vec=gauss_filt_unifrom(in_vec,sigma)
in_vec=col_vec(in_vec);
out_vec=gaussfilt(1:numel(in_vec),in_vec,sigma);
end