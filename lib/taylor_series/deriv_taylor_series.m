function out=deriv_taylor_series(x,derivs,offset,deriv_order)
%taking the derivative of a taylor series just shifts the coeffecients down by one
derivs=derivs(:);
if numel(derivs)<1+deriv_order
    warning('taylor series not defined to that derivative order')
    out=0;
else
derivs=derivs(1+deriv_order:end);
out=taylor_series(x,derivs,offset);
end
end