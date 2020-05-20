function out=isdouble_integer(in)
% see if a double is (close) to an integer
% the tolerance defaults at 2 eps

match_tol=eps*2;
out=abs((floor(in)-in))<match_tol;
nan_mask=isinf(in) | isnan(in) ;
out(nan_mask)=false;
 

end