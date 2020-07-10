function [deriv,mean_val]=num_grad_tensors(fun_in,x_in,delta_val)
% this function finds the numerical derivative of the output of fun_in with respect to the inputs
% the inputs may be arbitrary sized
% the inputs must be stacked in the first dimension (this is because matlab automaticaly removes trailing singleton
% dimensions)
% the output of the function can be arbitrary size
% the output has the format [seprate input index, [input tensor indexs], [output tensor indexes]] 
% todo
% alow for more complex numerical derivative methods
% Bryce Henson 2020-02-23

num_inputs=size(x_in,1);
x_in_size=size(x_in);
x_in_dims=ndims(x_in);
all_other_dims_input=repmat({':'},[1,x_in_dims-1]);
num_elements=prod(x_in_size(2:end));

if nargout>1 % the mean of the finite difference can be usefull for some indication of the nonlinearity of the function
    do_mean=true;
    mean_val=x_in*nan;
else
    do_mean=false;
end

% find the dimensionality of the output of the function
x_in_single=x_in(1,all_other_dims_input{:});
fun_out_tmp=fun_in(x_in_single);
fun_out_size=size(fun_out_tmp);
if isequal(fun_out_size,[1,1])
    fun_out_dims={};
    deriv=nan(size(x_in)); %intialize output
else
    fun_out_dims=repmat({':'},[1,size(fun_out_size,2)]);
    deriv=nan([size(x_in),fun_out_size]);
end



delta_tensor_template=zeros(x_in_size(2:end));
delta_tensor_template=permute(delta_tensor_template,[x_in_dims,1:x_in_dims-1]);
sub_tmp=cell([1,x_in_dims-1]);
for ii=1:num_inputs
    x_in_single=x_in(ii,all_other_dims_input{:});
    for jj=1:num_elements
        [sub_tmp{:}]=ind2sub(x_in_size(2:end),jj);
        delta_tensor=delta_tensor_template;
        delta_tensor(1,sub_tmp{:})=delta_val;
        fdp=fun_in(x_in_single+delta_tensor);
        fdn=fun_in(x_in_single-delta_tensor);
        deriv(ii,sub_tmp{:},fun_out_dims{:})= (fdp-fdn)/(2*delta_val);
        if do_mean 
            mean_val(ii,sub_tmp{:},fun_out_dims{:})= (fdp+fdn)/(2*delta_val);
        end
    end
end


end