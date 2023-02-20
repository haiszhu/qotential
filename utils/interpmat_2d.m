function [L,V,R] = interpmat_2d(t,s,idx)
% interpolation matrix from nodes s in 2D to target nodes
%
% size(s) = 2 * p^2
% size(t) = 2 * q
% 
% p is the order in x and y direction
% q is the number of target points...
%
% Hai 11/12/20

p = sqrt(numel(s)/2);
% matrix to solve for coefficients
Vx = ones(p^2,p); for j=2:p, Vx(:,j) = Vx(:,j-1).*s(1,:)'; end   % polyval matrix on nodes for x-coord
Vy = ones(p^2,p); for j=2:p, Vy(:,j) = Vy(:,j-1).*s(2,:)'; end   % polyval matrix on nodes for y-coord
V = ones(p^2,p^2); for j=1:p^2, tmp = Vx(j,:)'*Vy(j,:); V(j,:) = tmp(:)'; end
% matrix to apply coefficients to targets
q = numel(t(1,:));
if nargin < 3
    Rx = ones(q,p); for j=2:p, Rx(:,j) = Rx(:,j-1).*t(1,:)'; end   % polyval matrix on targs for x-coord
    Ry = ones(q,p); for j=2:p, Ry(:,j) = Ry(:,j-1).*t(2,:)'; end   % polyval matrix on targs for y-coord
    R = ones(q,p^2); for j=1:q, tmp = Rx(j,:)'*Ry(j,:); R(j,:) = tmp(:)'; end
else
    Rx = ones(q,p); for j=2:p, Rx(:,j) = Rx(:,j-1).*t(1,:)'; end   % polyval matrix on targs for x-coord
    Rx(idx,10:end) = 0;
    Ry = ones(q,p); for j=2:p, Ry(:,j) = Ry(:,j-1).*t(2,:)'; end   % polyval matrix on targs for y-coord
    Ry(idx,10:end) = 0;
    R = ones(q,p^2); for j=1:q, tmp = Rx(j,:)'*Ry(j,:); R(j,:) = tmp(:)'; end
end
% together
L = (V'\R')';

end