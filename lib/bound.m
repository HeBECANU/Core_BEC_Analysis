function y = bound(x,lower,upper)
    % limits the retuned value to lower or upper
  y=min(max(x,lower),upper);
end