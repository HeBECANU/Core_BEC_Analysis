function out=inlidx(matrix,index)
% return a single index of an array, very handy for making inline functions
% a=magic(5)
% inlidx(a,[1,2])
% inlidx(a,{':',2})

if nargin>1

if isvector(index) && ~iscell(index)
    index=num2cell(index);
elseif  ~iscell(index)
    error('index wrong')
end

if size(size(matrix))~=numel(index)
    error('cell vector must be the same size as matrix dimensions')
end

%inline index
out=matrix(index{:});
else
    out=matrix;
end

end