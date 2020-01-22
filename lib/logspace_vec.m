function mat_out=logspace_vec(x1,x2,num_steps)
% a simple extenstion to linspace that steps between vectors
% linspace_vec([1,0],[0,1],10)

x1=col_vec(x1);
x2=col_vec(x2);

if size(x1)~=size(x2)
    error('start and end vectors must be the same size')
end

if nargin<3
    num_steps=100;
end

mat_out=zeros(size(x1,1),num_steps);

for ii=1:size(x1,1)
    mat_out(ii,:)=logspace(x1(ii),x2(ii),num_steps);
end


end