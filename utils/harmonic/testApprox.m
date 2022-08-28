% switch from Legendre nodes to vioreanu nodes
% test approximation on a triangular patch...
%
% constant together with high order...
% pick a reference point that is different from quadrature nodes
%
%
% Hai 07/22/22

fsurf = @(x,y) 1/2*(sin(x)+cos(y)+sin(2*x.*y)+x.*y-y.^2+x.^3+x+1/2*y+1/4*x); % preferably a flat one
ftest = @(x,y) sin(3*x)+cos(2*y)+x.^2+y.^3+exp(x.*y); %+sin(3*x.*(y+1/3)+1/2);

Scale = [1/2,1/4,1/8,1/12,1/16,1/24,1/32,1/64];
Err = zeros(size(Scale));
order = 12; % <= 21 for Vioreanu nodes
ordert = 41; % targets for verification

% setup approximation nodes
addpath ../

figure(1),clf,
figure(2),clf,
figure(3),clf,
for j_scale = 1:numel(Scale) 
    scale = Scale(j_scale);
    [uvs,wts]=get_vioreanu_nodes(order-1); 
    uvs = scale*2*(uvs-1/2);
    
    figure(1), subplot(2,ceil(numel(Scale)/2),j_scale)
    %plot3(uvs(1,:),uvs(2,:),fsurf(uvs(1,:),uvs(2,:)),'.'); axis equal, hold on
    scatter3(uvs(1,:),uvs(2,:),fsurf(uvs(1,:),uvs(2,:)),16*ones(size(uvs(1,:))),ftest(uvs(1,:),uvs(2,:)),'filled')
    axis equal, hold on, colorbar
    scatter3(-uvs(1,:),-uvs(2,:),fsurf(-uvs(1,:),-uvs(2,:)),16*ones(size(uvs(1,:))),ftest(-uvs(1,:),-uvs(2,:)),'filled')
    sgtitle(['surface and density'])

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
    
    figure(2), subplot(2,ceil(numel(Scale)/2),j_scale),
    plot3(Sxl(1,:),Sxl(2,:),Sxl(3,:),'.r'); axis equal, hold on
    sgtitle(['surface after transformation'])

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
    n = order; [fx,fy,fz,gradF] = evalHarmonicGrad(Sxl,n,1);   % harmonic gradient
    rhs = [mul, zeros(1,3*order*(order+1)/2)]';
    F0 = zeros(order*(order+1)/2);
    F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3;
    Mmatrix = [[ F0 -F1 -F2 -F3];...
               [ F1  F0 -F3  F2];...
               [ F2  F3  F0 -F1];...
               [ F3 -F2  F1  F0]];  % eq(27)
    soln = Mmatrix\rhs;
    
    % test approximation at a given target
    testpt1 =  Sxlt;
    [fxpt,fypt,fzpt,gradFt] = evalHarmonicGrad(testpt1,n,1);
    fxpt{1}{1} = 0*fxpt{1}{1}; fypt{1}{1} = 0*fypt{1}{1}; fzpt{1}{1} = ones(size(fzpt{1}{1}));
    F0t = zeros(size(testpt1,2),order*(order+1)/2); 
    F1t = gradFt.F1; F2t = gradFt.F2; F3t = gradFt.F3;
    Mu_parts = [[ F0t -F1t -F2t -F3t];...
                [ F1t  F0t -F3t  F2t];...
                [ F2t  F3t  F0t -F1t];...
                [ F3t -F2t  F1t  F0t]]*soln;
    diff  = Mu_parts(1:end/4)-mut(idxlt)';
    mui = Mu_parts(end/4+1:end/2); 
    muj = Mu_parts(end/2+1:3*end/4); 
    muk = Mu_parts(3*end/4+1:end);

    err = NaN(size(xt1));
    err(idxlt) = diff;
    figure(3), subplot(2,ceil(numel(Scale)/2),j_scale),imagesc(log10(abs(err)))
    colorbar, colormap(jet)
    caxis([-15 0])
    sgtitle(['order = ', num2str(order),''])
    
    max(abs(err(:)))
    Err(j_scale) = max(abs(err(:)));

end

figure(4), 
loglog(Scale,Err,'o-')
hold on
errfit = @(x) Err(1)/Scale(1)^order*x.^order;
loglog(Scale,errfit(Scale),'--r')

keyboard


