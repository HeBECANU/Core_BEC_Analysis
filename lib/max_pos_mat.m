function [pos,max_num]=max_pos_mat(mat_in)
[max_num,max_idx] = max(mat_in(:));
pos=cell(size(size(mat_in)));
[pos{:}]=ind2sub(size(mat_in),max_idx);
pos=cell2mat(pos);
end