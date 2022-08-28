function col_sum = q_sum(q_vec)
% sum up columns of a quaternion matrix, or all entries of a vector
% Hai 08/25/22

if size(q_vec,1) > 1 && size(q_vec,2) > 1
  col_sum_parts = NaN(size(q_vec,1),4);
  for j=1:size(col_sum_parts,1)
    col_sum_parts(j,:) = sum(compact(q_vec(j,:))); % 
  end
else
  col_sum_parts = sum(compact(q_vec)); % 1 by 4 
end
col_sum = quaternion(col_sum_parts); 
 
end