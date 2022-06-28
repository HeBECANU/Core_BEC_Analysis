function out=is_a_within_factor_of_b(a,b,factor)
% find if a is within some factor of b
% i find this usefull for error checking often i care about the relative difference
    out= a<b*factor &  a>b/factor;

end