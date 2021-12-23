function [fx,fy,fz] = evalHarmonicGrad(r,n)
% numerical evaluation on scaled gradient of harmonic polynoimals, no
% dependent on target...
%
% Hai 07/29/20

x = r(1,:); y = r(2,:); z = r(3,:);
rho2 = x.^2+y.^2+z.^2;  % rho square

c = legendre_coefficients(n);
Cp = rescaled_derivative_coefficients(c);
Hxy = real(bsxfun(@power,x+1i*y,(1:n)')); 
Hxy_x = bsxfun(@times,real(bsxfun(@power,x+1i*y,(0:n-1)')),(1:n)'); 
Hxy_y = bsxfun(@times,-imag(bsxfun(@power,x+1i*y,(0:n-1)')),(1:n)');

for k=1:n   % loop over different order of derivative of n-th degree Legendre polynomial, since this is a mix of (x,y) and z, hard to precompute
    if mod(n,2)==0
        cpnk = Cp{k}(end,(2-mod(k-1,2)):2:end);                                         % k-th order derivative coefficient of 10th degree Legendre polynomial
        Z = bsxfun(@power,z,(mod(k,2):2:n-k)');                                         % various power of z
        if (mod(k,2)-1)==0
            Zp = bsxfun(@times,bsxfun(@power,z,(mod(k,2)-1:2:n-k-1)'),(mod(k,2):2:n-k)');   % derivative of various power of z
        else
            Zp = bsxfun(@times,bsxfun(@power,z,(mod(k,2)-1:2:n-k-1)'),(mod(k,2):2:n-k)');   % derivative of various power of z
            Zp(:,z==0) = 0;
        end
    else
        cpnk = Cp{k}(end,(2-mod(k,2)):2:end); 
        Z = bsxfun(@power,z,(mod(k-1,2):2:n-k)');
        if (mod(k-1,2)-1)==0
            Zp = bsxfun(@times,bsxfun(@power,z,(mod(k-1,2)-1:2:n-k-1)'),(mod(k-1,2):2:n-k)');
        else
            Zp = bsxfun(@times,bsxfun(@power,z,(mod(k-1,2)-1:2:n-k-1)'),(mod(k-1,2):2:n-k)');
            Zp(:,z==0) = 0;     % this is really annoying, hopefully won't need to come back and fix more...
        end
    end
    nl = numel(cpnk);                                                              % num of nonzero coefficients
    Rho2 = bsxfun(@power,rho2,(0:nl-1)');                                          % various power of rho square
    Rho2_x = bsxfun(@times,2*x,(0:nl-1)').*[zeros(size(x));Rho2(1:end-1,:)];       % partial derivative of Rho2 w.r.t x
    Rho2_y = bsxfun(@times,2*y,(0:nl-1)').*[zeros(size(y));Rho2(1:end-1,:)];       % partial derivative of Rho2 w.r.t y
    Rho2_z = bsxfun(@times,2*z,(0:nl-1)').*[zeros(size(z));Rho2(1:end-1,:)];       % partial derivative of Rho2 w.r.t z
    ZRho2 = Z.*Rho2(end:-1:1,:);                                                    % z^(j-k+1)*rho^(n-j-1)
    ZpRho2 = Zp.*Rho2(end:-1:1,:);                                                  % (j-k+1)*z^(j-k)*rho^(n-j-1)
    ZRho2_x = Z.*Rho2_x(end:-1:1,:);                                                % z^(j-k+1)*d(rho^(n-j-1))/dx
    ZRho2_y = Z.*Rho2_y(end:-1:1,:);                                                % z^(j-k+1)*d(rho^(n-j-1))/dy
    ZRho2_z = Z.*Rho2_z(end:-1:1,:);                                                % z^(j-k+1)*d(rho^(n-j-1))/dy

    for j=k:-1:1
        m = n-(k-j);    % approx order
        if mod(m,2)==0
            cpmj = Cp{j}(end-(k-j),(2-mod(j-1,2)):2:(end-(k-j)));
        else
            cpmj = Cp{j}(end-(k-j),(2-mod(j,2)):2:(end-(k-j))); 
        end
        
        lj = cpmj*ZRho2;
        lj_x = cpmj*ZRho2_x;
        lj_y = cpmj*ZRho2_y;
        lj_z = cpmj*ZRho2_z + cpmj*ZpRho2;
        
        hj = Hxy(j,:);
        hj_x = Hxy_x(j,:);
        hj_y = Hxy_y(j,:);
        
        fx{m}{j} = hj_x.*lj+hj.*lj_x; 
        fy{m}{j} = hj_y.*lj+hj.*lj_y; 
        fz{m}{j} = hj.*lj_z;
        
    end
end

end