function col_sum = q_sum(q_vec)
% sum up columns of a quaternion matrix, or all entries of a vector
% Hai 08/25/22

if size(q_vec,1) > 1 && size(q_vec,2) > 1
%   col_sum_parts = NaN(size(q_vec,1),4);
%   for j=1:size(col_sum_parts,1)
%     col_sum_parts(j,:) = sum(compact(q_vec(j,:))); % 
%   end
  tmp = compact(q_vec); % all parts, matlab reshapes to dim = []*4
  col_sum0 = sum(reshape(tmp(:,1),size(q_vec)),2);
  col_sum1 = sum(reshape(tmp(:,2),size(q_vec)),2);
  col_sum2 = sum(reshape(tmp(:,3),size(q_vec)),2);
  col_sum3 = sum(reshape(tmp(:,4),size(q_vec)),2);
  col_sum = quaternion(col_sum0,col_sum1,col_sum2,col_sum3); 
else
  col_sum_parts = sum(compact(q_vec)); % 1 by 4 
  col_sum = quaternion(col_sum_parts); 
end
 
end