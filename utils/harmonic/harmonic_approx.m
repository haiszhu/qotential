function [Err,Tx1,Tx2]=harmonic_approx(fsurf,ftest,ref,order)
% approximation on [-1/2,1/2]^2 using harmonic_basis_gradient.m 
% 
% Hai 07/25/22

% fsurf = @(x,y) 1/2*(sin(x)+cos(y)+sin(2*x.*y)+x.*y-y.^2+x.^3+x+1/2*y+1/4*x); % preferably a flat one
% ftest = @(x,y) sin(3*x)+cos(2*y)+x.^2+y.^3+exp(x.*y); %+sin(3*x.*(y+1/3)+1/2);
% ref = 4; order = 16;

Err = []; Tx1 = []; Tx2 = []; maxvel = 0;

% vioreanu nodes
origin = 0;
X0c=linspace(-1+1/ref,1-1/ref,ref); 
[uvs,wts] = get_vioreanu_nodes(order-1); 
gridsize = 1/ref; uvs = (uvs-1/2)*gridsize;

% target points
M = ceil(5*16*3*gridsize);
tij = linspace(-1,1,M)/2*gridsize; [t1ij t2ij] = meshgrid(tij);
midx_l = tril(ones(M)); midx_u = triu(ones(M));
midx_l = midx_l(end:-1:1,:); midx_u = midx_u(end:-1:1,:);


for i=1:numel(X0c)
    Erri=[]; Tx1i = []; Tx2i = [];
    for j=1:numel(X0c)
        
        x0c = 1/2*[X0c(i);X0c(j)];
        Errij=NaN(size(t1ij));

        %%  lower left
        % prepare: transformation for lower left triangular patch
        xxij = [x0c(1)+uvs(1,:);x0c(2)+uvs(2,:)]; 
        rl = [xxij;fsurf(xxij(1,:),xxij(2,:))];
        mul = ftest(xxij(1,:),xxij(2,:));
        % three vertices of the triangle
        rl1 = [x0c+gridsize/2*[-1;-1];fsurf(x0c(1)+gridsize/2*(-1),x0c(2)+gridsize/2*(-1))];
        rl2 = [x0c+gridsize/2*[ 1;-1];fsurf(x0c(1)+gridsize/2*( 1),x0c(2)+gridsize/2*(-1))];
        rl3 = [x0c+gridsize/2*[-1; 1];fsurf(x0c(1)+gridsize/2*(-1),x0c(2)+gridsize/2*( 1))];
        % xc (center of the path), nc (nornal direction, i.e. z-axis direction), nx (x-axis direction)
        xcl = 1/3*(rl1+rl2+rl3); 
        ncl = normal(rl1,rl2,rl3); % normal to flat triangle (three vertices r1,r2,r3)
        nxl = (rl2-rl1)/norm(rl2-rl1); % lower left point is the x-direction
        Sxl = transcoord(xcl,ncl,nxl,origin,rl);

        % targets to check approximation performance
        temp1_l=t1ij(logical(midx_l)); temp2_l=t2ij(logical(midx_l));
        ttij_l = [x0c(1)+temp1_l(:)';x0c(2)+temp2_l(:)']; 
        rlt = [ttij_l;fsurf(ttij_l(1,:),ttij_l(2,:))]; mult = ftest(ttij_l(1,:),ttij_l(2,:));
        Sxlt = transcoord(xcl,ncl,nxl,origin,rlt); 

        % setup matrix
        n = order; [fx,fy,fz] = evalHarmonicGrad(Sxl,n);   % harmonic gradient
        fx{1}{1} = 0*fx{1}{1}; fy{1}{1} = 0*fy{1}{1}; fz{1}{1} = ones(size(fz{1}{1})); % technically speaking, we can do grad(x), grad(y), or grad(z), grad(z) seems to be more convenient...
        rhs = zeros(4*order*(order+1)/2,1); rhs(1:4:end) = mul';
        Fc = []; Fi = []; Fj = []; Fk = [];
        for k=n:-1:1
            for j=1:k
                Fc = [Fc, zeros(size(Sxl(1,:)')),             -fx{k}{j}',             -fy{k}{j}',             -fz{k}{j}'];
                Fi = [Fi,              fx{k}{j}', zeros(size(Sxl(1,:)')),             -fz{k}{j}',              fy{k}{j}'];
                Fj = [Fj,              fy{k}{j}',              fz{k}{j}', zeros(size(Sxl(1,:)')),             -fx{k}{j}'];
                Fk = [Fk,              fz{k}{j}',             -fy{k}{j}',              fx{k}{j}', zeros(size(Sxl(1,:)'))];
            end
        end
        Mmatrix = zeros(4*order*(order+1)/2);
        Mmatrix(1:4:end,:)=Fc; Mmatrix(2:4:end,:)=Fi; Mmatrix(3:4:end,:)=Fj; Mmatrix(4:4:end,:)=Fk;
        soln = Mmatrix\rhs;

        % test approximation at a given target
        testptl =  Sxlt;
        [fxpt,fypt,fzpt] = evalHarmonicGrad(testptl,n);
        fxpt{1}{1} = 0*fxpt{1}{1}; fypt{1}{1} = 0*fypt{1}{1}; fzpt{1}{1} = ones(size(fzpt{1}{1}));
        Fcpt = []; Fipt = []; Fjpt = []; Fkpt = [];
        for k=n:-1:1
            for j=1:k
                Fcpt = [Fcpt, zeros(size(testptl(1,:)')), -fxpt{k}{j}', -fypt{k}{j}', -fzpt{k}{j}'];
                Fipt = [Fipt, fxpt{k}{j}', zeros(size(testptl(1,:)')),-fzpt{k}{j}', fypt{k}{j}'];
                Fjpt = [Fjpt, fypt{k}{j}', fzpt{k}{j}', zeros(size(testptl(1,:)')),-fxpt{k}{j}'];
                Fkpt = [Fkpt, fzpt{k}{j}',-fypt{k}{j}', fxpt{k}{j}', zeros(size(testptl(1,:)'))];
            end
        end
        mu0 = Fcpt*soln; 
        diff = mu0-mult'; mui = Fipt*soln; muj = Fjpt*soln; muk = Fkpt*soln;
        Errij(logical(midx_l)) = max([abs(diff),abs(mui),abs(muj),abs(muk)],[],2);

        maxvel = max(maxvel,max(abs(mult)));
        maxvel = max(maxvel,max(abs(rlt(3,:))));

        %%  upper right
        % prepare: transformation for lower left triangular patch
        xxij = [x0c(1)-uvs(1,:);x0c(2)-uvs(2,:)]; 
        ru = [xxij;fsurf(xxij(1,:),xxij(2,:))];
        muu = ftest(xxij(1,:),xxij(2,:));
        % three vertices of the triangle
        ru1 = [x0c+gridsize/2*[ 1; 1];fsurf(x0c(1)+gridsize/2*( 1),x0c(2)+gridsize/2*( 1))];
        ru2 = [x0c+gridsize/2*[-1; 1];fsurf(x0c(1)+gridsize/2*(-1),x0c(2)+gridsize/2*( 1))];
        ru3 = [x0c+gridsize/2*[ 1;-1];fsurf(x0c(1)+gridsize/2*( 1),x0c(2)+gridsize/2*(-1))];
        % xc (center of the path), nc (nornal direction, i.e. z-axis direction), nx (x-axis direction)
        xcu = 1/3*(ru1+ru2+ru3); origin = 0;
        ncu = normal(ru1,ru2,ru3); % normal to flat triangle (three vertices r1,r2,r3)
        nxu = (ru2-ru1)/norm(ru2-ru1); % lower left point is the x-direction
        Sxu = transcoord(xcu,ncu,nxu,origin,ru); 

        % targets to check approximation performance
        temp1_u=t1ij(logical(midx_u)); temp2_u=t2ij(logical(midx_u));
        ttij_u = [x0c(1)+temp1_u(:)';x0c(2)+temp2_u(:)']; 
        rut = [ttij_u;fsurf(ttij_u(1,:),ttij_u(2,:))]; muut = ftest(ttij_u(1,:),ttij_u(2,:));
        Sxut = transcoord(xcu,ncu,nxu,origin,rut); 

        % setup matrix
        n = order; [fx,fy,fz] = evalHarmonicGrad(Sxu,n);   % harmonic gradient
        fx{1}{1} = 0*fx{1}{1}; fy{1}{1} = 0*fy{1}{1}; fz{1}{1} = ones(size(fz{1}{1})); % technically speaking, we can do grad(x), grad(y), or grad(z), grad(z) seems to be more convenient...
        rhs = zeros(4*order*(order+1)/2,1); rhs(1:4:end) = muu';
        Fc = []; Fi = []; Fj = []; Fk = [];
        for k=n:-1:1
            for j=1:k
                Fc = [Fc, zeros(size(Sxu(1,:)')),             -fx{k}{j}',             -fy{k}{j}',             -fz{k}{j}'];
                Fi = [Fi,              fx{k}{j}', zeros(size(Sxu(1,:)')),             -fz{k}{j}',              fy{k}{j}'];
                Fj = [Fj,              fy{k}{j}',              fz{k}{j}', zeros(size(Sxu(1,:)')),             -fx{k}{j}'];
                Fk = [Fk,              fz{k}{j}',             -fy{k}{j}',              fx{k}{j}', zeros(size(Sxu(1,:)'))];
            end
        end
        Mmatrix = zeros(4*order*(order+1)/2);
        Mmatrix(1:4:end,:)=Fc; Mmatrix(2:4:end,:)=Fi; Mmatrix(3:4:end,:)=Fj; Mmatrix(4:4:end,:)=Fk;
        soln = Mmatrix\rhs;

        % test approximation at a given target
        testptu =  Sxut;
        [fxpt,fypt,fzpt] = evalHarmonicGrad(testptu,n);
        fxpt{1}{1} = 0*fxpt{1}{1}; fypt{1}{1} = 0*fypt{1}{1}; fzpt{1}{1} = ones(size(fzpt{1}{1}));
        Fcpt = []; Fipt = []; Fjpt = []; Fkpt = [];
        for k=n:-1:1
            for j=1:k
                Fcpt = [Fcpt, zeros(size(testptu(1,:)')), -fxpt{k}{j}', -fypt{k}{j}', -fzpt{k}{j}'];
                Fipt = [Fipt, fxpt{k}{j}', zeros(size(testptu(1,:)')),-fzpt{k}{j}', fypt{k}{j}'];
                Fjpt = [Fjpt, fypt{k}{j}', fzpt{k}{j}', zeros(size(testptu(1,:)')),-fxpt{k}{j}'];
                Fkpt = [Fkpt, fzpt{k}{j}',-fypt{k}{j}', fxpt{k}{j}', zeros(size(testptu(1,:)'))];
            end
        end
        mu0 = Fcpt*soln; 
        diff = mu0-muut'; mui = Fipt*soln; muj = Fjpt*soln; muk = Fkpt*soln;
        Errij(logical(midx_u)) = max([abs(diff),abs(mui),abs(muj),abs(muk)],[],2);

        maxvel = max(maxvel,max(abs(muut)));
        maxvel = max(maxvel,max(abs(rut(3,:))));

        %
        Erri=[Erri;Errij];
        Tx1ij=x0c(1)+t1ij; Tx1i=[Tx1i;Tx1ij];
        Tx2ij=x0c(2)+t2ij; Tx2i=[Tx2i;Tx2ij];
        % keyboard

    end
    Err = [Err, Erri];  Tx1 = [Tx1, Tx1i]; Tx2 = [Tx2, Tx2i];
end
Err = Err/maxvel;

end
 