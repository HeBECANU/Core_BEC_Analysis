function y = bound(x,lower,upper)
     % TODO
     % - use match_tensor_sizes to match mismatching inputs
    
    if ~isequal(size(x),size(lower),size(upper))
        [x,lower,upper]=match_tensor_sizes_multi_in(x,lower,upper);
        %error('input sizes must be equal')
    end
    % limits the retuned value to lower or upper
    y=min(max(x,lower),upper);
end