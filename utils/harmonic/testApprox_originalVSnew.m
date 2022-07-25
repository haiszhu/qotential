% verify that approximations are equivalent
% may take a minute or two to see the result
%
% this is similar to the setup of Fig3 from "High-order close evaluation of Laplace layer potentials: A differential geometry approach" paper 
% in that case, fsurf & ftest are all admissible pairs of tensor-product monic polynomials
% here is just 1 specific example...
% 
% order for original basis should be less than 7
% order for new basis should be less than 20
%
% Hai 07/25/22


fsurf = @(x,y) 1/2*(sin(x)+cos(3*y)+sin(2*x.*y)+x.*y-y.^2+x.^3+x+1/2*y+1/4*x); % preferably a flat one
ftest = @(x,y) sin(3*x)+cos(2*y)+x.^2+y.^3+exp(x.*y); %+sin(3*x.*(y+1/3)+1/2);

Ref = [2,4,8,12,16,24,32,48];
Erro = zeros(size(Ref));
Err  = zeros(size(Ref));

ordero = 6;     % approximation order for original basis, <= 7
order = ordero; % same order for new basis rho^p*cos(k\theta)*P^k_p(cos(\theta)), which gives about the same performance
order = 12;     % or increase the approximation order to see faster convergence rate

figure(1),clf, 
warning('off','MATLAB:nearlySingularMatrix')
for k=1:numel(Ref)
    ref = Ref(k);  
    
    [erro,Tx1o,Tx2o]=harmonic_approx_original(fsurf,ftest,ref,ordero);
    Erro(k) = max(erro(:));

    
    
    [err,Tx1,Tx2]=harmonic_approx(fsurf,ftest,ref,order);
    Err(k) = max(err(:));

    subplot(1,3,1)
    imagesc(Tx1o(1,:),Tx2o(:,1),log10(erro)); hold on, axis equal tight,
    colorbar, caxis([-15 0]);
    title([' original basis error, p = ', num2str(ordero)])
    drawnow

    subplot(1,3,2)
    imagesc(Tx1(1,:),Tx2(:,1),log10(err)); hold on, axis equal tight,
    colorbar, caxis([-15 0])
    title([' basis \rho^p cos(k\theta)P^k_p(cos(\theta)) error, p = ', num2str(order)])
    drawnow

    pause(0.5)

end
warning('on','MATLAB:nearlySingularMatrix')

subplot(1,3,3)
loglog(Ref,Erro,'o-r'), hold on
errfito = @(x) Erro(1)*Ref(1)^ordero*x.^(-ordero);
loglog(Ref,errfito(Ref),'--r')
loglog(Ref,Err,'o-b'), hold on
errfit  = @(x) Err(1)*Ref(1)^order*x.^(-order);
loglog(Ref,errfit(Ref),'--b')
title([' convergence'])


keyboard

