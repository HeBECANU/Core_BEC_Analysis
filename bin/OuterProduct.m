% A = magic(3);
% B = [1,2,3]';
% C = OuterProduct(A,B)

function C = OuterProduct(A, B)  % version 4
%     https://au.mathworks.com/matlabcentral/answers/445798-outer-product-of-two-rectangular-matrices
C = A .* reshape(B, [ones(1, ndims(A)), size(B)]);
% Matlab < R2016b:
% C = bsxfun(@times, A, reshape(B, [ones(1, ndims(A)), size(B)]))
end