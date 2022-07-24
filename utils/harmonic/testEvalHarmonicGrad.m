% check gradient of harmonic polynomials
% universal eval vs explicit formula up to order 10
%
% Hai 07/24/22

fsurf = @(x,y) 1/2*(sin(x)+cos(y)+sin(2*x.*y)+x.*y-y.^2+x.^3+x+1/2*y+1/4*x); % preferably a flat one
ftest = @(x,y) sin(3*x)+cos(2*y)+x.^2+y.^3+exp(x.*y); %+sin(3*x.*(y+1/3)+1/2);

order = 11; % <= 21 for Vioreanu nodes
ordert = 41; % targets for verification

% setup approximation nodes
addpath ../
scale = 1/2;
[uvs,wts]=get_vioreanu_nodes(order-1); 
uvs = scale*2*(uvs-1/2);

figure(1),clf,
%plot3(uvs(1,:),uvs(2,:),fsurf(uvs(1,:),uvs(2,:)),'.'); axis equal, hold on
scatter3(uvs(1,:),uvs(2,:),fsurf(uvs(1,:),uvs(2,:)),16*ones(size(uvs(1,:))),ftest(uvs(1,:),uvs(2,:)),'filled')
axis equal, hold on, colorbar
scatter3(-uvs(1,:),-uvs(2,:),fsurf(-uvs(1,:),-uvs(2,:)),16*ones(size(uvs(1,:))),ftest(-uvs(1,:),-uvs(2,:)),'filled')

% prepare 1: transformation for lower left triangular patch
rl = [uvs(1,:);uvs(2,:);fsurf(uvs(1,:),uvs(2,:))];
mul = ftest(uvs(1,:),uvs(2,:));
% three vertices of the triangle
rl1 = [scale*[-1;-1];fsurf(scale*(-1),scale*(-1))];
rl2 = [scale*[ 1;-1];fsurf(scale*( 1),scale*(-1))];
rl3 = [scale*[-1; 1];fsurf(scale*(-1),scale*( 1))];
% xc (center of the path), nc (nornal direction, i.e. z-axis direction), nx (x-axis direction)
xcl = 1/3*(rl1+rl2+rl3); origin = 0;
ncl = normal(rl1,rl2,rl3); % normal to flat triangle (three vertices r1,r2,r3)
nxl = (rl2-rl1)/norm(rl2-rl1); % lower left point is the x-direction
Sxl = transcoord(xcl,ncl,nxl,origin,rl);

figure(2), clf,
plot3(Sxl(1,:),Sxl(2,:),Sxl(3,:),'.r'); axis equal, hold on

% targets to check approximation performance
xt = linspace(-1,1,ordert); xt = scale*xt;
[xt1 xt2] = meshgrid(xt); xxt = [xt1(:)';xt2(:)']; 
rt = [xxt;fsurf(xxt(1,:),xxt(2,:))]; mut = ftest(xxt(1,:),xxt(2,:));
idxlt = tril(ones(ordert)); idxut = triu(ones(ordert)); 
idxlt = idxlt(end:-1:1,:); idxut = idxut(end:-1:1,:);   % logical value for lower left and upper right tri
idxt = 1:ordert^2; idxlt = idxt(logical(idxlt)); idxut = idxt(logical(idxut));  % idx for lower left and upper right tri
rlt = rt(:,idxlt); rut = rt(:,idxut);     % coordinates for two triangles
Sxlt = transcoord(xcl,ncl,nxl,origin,rlt); Sxut = transcoord(xcl,ncl,nxl,origin,rut); 

%%  lower left
% setup matrix
n = order; [fx,fy,fz] = evalHarmonicGrad(Sxl,n);    % harmonic gradient
[fx2,fy2,fz2] = harmonic_basis_gradient(n);     % explicit formula

Err = [];
for m=1:10
    for n=1:m
        diffx = fx{m}{n} - fx2{m}{n}(Sxl(1,:),Sxl(2,:),Sxl(3,:));
        diffy = fy{m}{n} - fy2{m}{n}(Sxl(1,:),Sxl(2,:),Sxl(3,:));
        diffz = fz{m}{n} - fz2{m}{n}(Sxl(1,:),Sxl(2,:),Sxl(3,:));
        Err = max([Err abs([diffx diffy diffz])]);
    end
end

Err

keyboard

