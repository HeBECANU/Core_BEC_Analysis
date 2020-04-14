function tensor_out=checkerboard_tensor(size_in)
% makes an alternating pattern of logicals
% Inputs
% size_in - size of the output tensor (arb. dim. array)
% Output
% tensor_out - tensor of locicals with alternating true false
% Bryce Henson 2020-04-14
% TODO
% - test cases
% - phase shift
% - arb period


size_in=row_vec(size_in);

tensor_out=zeros(size_in);
for ii=1:size(size_in,2)
    this_dim_vec= 1 : size_in(ii);
    this_dim_vec=(mod(this_dim_vec, 2)==0);
    % rotate the vector to point in the right direction
    tmp_reshape_size=0*size_in +1;
    tmp_reshape_size(ii)=size_in(ii);
    this_dim_vec=reshape(this_dim_vec,tmp_reshape_size);
    
    % then repeat it in all the other directions
    tmp_repmat_vec=size_in;
    tmp_repmat_vec(ii)=1;
    repeated_vec = repmat(this_dim_vec,  tmp_repmat_vec);   
    % then sum it with the cumulative tensor
    tensor_out=tensor_out+repeated_vec;
end

% turn this into a logivcal by finding if the modulus is zero
tensor_out=(mod(tensor_out,2)==0);

end