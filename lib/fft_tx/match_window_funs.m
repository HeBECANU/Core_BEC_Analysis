function [opt_match_val,cost_val]=match_window_funs(handle_funref,handle_funadj,len,opt_start,norm_order,verbose)
% given one window function find the parameters of another that match in the time domain
% example usage
% match_window_funs(@(n) gausswin(n,5), @(n,p) chebwin(n,p),1e3,100,[],[])


if isempty(verbose) || isnan(verbose)
    verbose=2;
end
%check that they return the same size
ref_data=handle_funref(len);
match_data=handle_funadj(len,opt_start);
if ~isequal(size(match_data),size(ref_data))
    error('functions do not return the same size')
end

if isempty(norm_order) || isnan(norm_order)
    norm_order=2;
end
    

cost_fun=@(x) norm(ref_data-handle_funadj(len,x),norm_order);
[opt_match_val,cost_val]=fminsearch(cost_fun,opt_start);

fun_handle_str={func2str(handle_funref),func2str(handle_funadj)};

if verbose>1
    sfigure(1);
    set(gcf,'color','w')
    subplot(2,1,1)
    plot(handle_funref(len))
    hold on
    plot(handle_funadj(len,opt_match_val));
    hold off
    xlabel('index')
    ylabel('amplitude')
    legend(fun_handle_str)
    subplot(2,1,2)
    plot(handle_funref(len)-handle_funadj(len,opt_match_val))
    ylabel('window difference')
end

end