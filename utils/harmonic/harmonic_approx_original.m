function [Err,Tx1,Tx2]=harmonic_approx_original(fsurf,ftest,ref,order)
% approximation on [-1/2,1/2]^2 using harmonic_basis_gradient_original.m 
% 

origin = 0;
[x0 w0] = gauss(order);
X0c=linspace(-1+1/ref,1-1/ref,ref); 
x=[]; w=[]; for i=1:ref, x=[x;x0/ref+X0c(i)]; w=[w;w0/ref]; end
x=1/2*x; w=1/2*w;
[x1 x2] = meshgrid(x); xx = [x1(:)';x2(:)'];  
gridsize = 1/ref;

[fc,fi,fj,fk,idx1st,idx2nd] = harmonic_basis_gradient_original(order); % explicit formula
err = [];
Err = []; Tx1 = []; Tx2 = []; maxvel = 0;

M = ceil(5*16*3*gridsize);
tij = linspace(-1,1,M)/2*gridsize; [t1ij t2ij] = meshgrid(tij);
xij = x0/2*gridsize; [x1ij x2ij] = meshgrid(xij);
midx_l = tril(ones(M)); midx_u = triu(ones(M));
midx_l = midx_l(end:-1:1,:); midx_u = midx_u(end:-1:1,:);
for i=1:numel(X0c)
    Erri=[]; Tx1i = []; Tx2i = [];
    for j=1:numel(X0c)
        Errij=NaN(size(t1ij));
        
        x0c = 1/2*[X0c(i);X0c(j)];
        xxij = [x0c(1)+x1ij(:)';x0c(2)+x2ij(:)']; 
        temp.x = [xxij;fsurf(xxij(1,:),xxij(2,:))];
        temp.hx = ftest(xxij(1,:),xxij(2,:));
        % two triangles
        A = temp.x(:,1); B = temp.x(:,order^2-order+1);
        C = temp.x(:,order^2); D = temp.x(:,order);
        xc1 = A; nc1 = normal(A,B,D); nx1 = (B-xc1)/norm(B-xc1);
        xc2 = C; nc2 = normal(B,C,D); nx2 = (D-xc2)/norm(D-xc2);
        % lower left
        sx_1 = transcoord(xc1,nc1,nx1,origin,temp.x(:,idx1st));
        sigma_1 = temp.hx(:,idx1st) - temp.hx(:,1);
        Mmatrix = zeros(2*(order+2)*(order-1)); rhs = zeros(2*(order+2)*(order-1),1);
        Mmatrix(1:4:end,:) = fc(sx_1(1,:)',sx_1(2,:)',sx_1(3,:)');
        Mmatrix(2:4:end,:) = fi(sx_1(1,:)',sx_1(2,:)',sx_1(3,:)');
        Mmatrix(3:4:end,:) = fj(sx_1(1,:)',sx_1(2,:)',sx_1(3,:)');
        Mmatrix(4:4:end,:) = fk(sx_1(1,:)',sx_1(2,:)',sx_1(3,:)');
        rhs(1:4:end) = sigma_1;
        soln = Mmatrix\rhs;
        
        temp1_l=t1ij(logical(midx_l)); temp2_l=t2ij(logical(midx_l));
        ttij_l = [x0c(1)+temp1_l(:)';x0c(2)+temp2_l(:)']; 
        tx_temp = [ttij_l;fsurf(ttij_l(1,:),ttij_l(2,:))]; 
        tx_1_exa = transcoord(xc1,nc1,nx1,origin,tx_temp);
        tsigma_1_exa = ftest(tx_temp(1,:),tx_temp(2,:))- temp.hx(:,1);
        
        maxvel = max(maxvel,max(abs(tx_temp(3,:))));
        maxvel = max(maxvel,max(abs(tsigma_1_exa+temp.hx(:,1))));
        
        vecc1 = fc(tx_1_exa(1,:)',tx_1_exa(2,:)',tx_1_exa(3,:)')*soln-tsigma_1_exa(:);
        veci1 = fi(tx_1_exa(1,:)',tx_1_exa(2,:)',tx_1_exa(3,:)')*soln;
        vecj1 = fj(tx_1_exa(1,:)',tx_1_exa(2,:)',tx_1_exa(3,:)')*soln;
        veck1 = fk(tx_1_exa(1,:)',tx_1_exa(2,:)',tx_1_exa(3,:)')*soln;
        
        Errij(logical(midx_l)) = max([abs(vecc1),abs(veci1),abs(vecj1),abs(veck1)],[],2);
        
        % upper right
        sx_2 = transcoord(xc2,nc2,nx2,origin,temp.x(:,idx2nd));
        sigma_2 = temp.hx(:,idx2nd) - temp.hx(:,order^2);
        Mmatrix = zeros(2*(order+2)*(order-1)); rhs = zeros(2*(order+2)*(order-1),1);
        Mmatrix(1:4:end,:) = fc(sx_2(1,:)',sx_2(2,:)',sx_2(3,:)');
        Mmatrix(2:4:end,:) = fi(sx_2(1,:)',sx_2(2,:)',sx_2(3,:)');
        Mmatrix(3:4:end,:) = fj(sx_2(1,:)',sx_2(2,:)',sx_2(3,:)');
        Mmatrix(4:4:end,:) = fk(sx_2(1,:)',sx_2(2,:)',sx_2(3,:)');
        rhs(1:4:end) = sigma_2;
        soln = Mmatrix\rhs;
        
        temp1_u=t1ij(logical(midx_u)); temp2_u=t2ij(logical(midx_u));
        ttij_u = [x0c(1)+temp1_u(:)';x0c(2)+temp2_u(:)']; 
        tx_temp = [ttij_u;fsurf(ttij_u(1,:),ttij_u(2,:))];
        tx_2_exa = transcoord(xc2,nc2,nx2,origin,tx_temp);
        tsigma_2_exa = ftest(tx_temp(1,:),tx_temp(2,:))- temp.hx(:,order^2);
        
        maxvel = max(maxvel,max(abs(tx_temp(3,:))));
        maxvel = max(maxvel,max(abs(tsigma_2_exa+temp.hx(:,order^2))));
        
        vecc2 = fc(tx_2_exa(1,:)',tx_2_exa(2,:)',tx_2_exa(3,:)')*soln-tsigma_2_exa(:);
        veci2 = fi(tx_2_exa(1,:)',tx_2_exa(2,:)',tx_2_exa(3,:)')*soln;
        vecj2 = fj(tx_2_exa(1,:)',tx_2_exa(2,:)',tx_2_exa(3,:)')*soln;
        veck2 = fk(tx_2_exa(1,:)',tx_2_exa(2,:)',tx_2_exa(3,:)')*soln;
        
        Errij(logical(midx_u)) = max([abs(vecc2),abs(veci2),abs(vecj2),abs(veck2)],[],2);
        err = [err;max(abs([vecc1;veci1;vecj1;veck1;vecc2;veci2;vecj2;veck2]))];
%         keyboard
        Erri=[Erri;Errij];
        Tx1ij=x0c(1)+t1ij; Tx1i=[Tx1i;Tx1ij];
        Tx2ij=x0c(2)+t2ij; Tx2i=[Tx2i;Tx2ij];
    end
    Err = [Err, Erri];  Tx1 = [Tx1, Tx1i]; Tx2 = [Tx2, Tx2i];
end
Err = Err/maxvel;
end
