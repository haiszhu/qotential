function [fx,fy,fz] = evalHarmonicGrad(r,n)
% numerical evaluation on scaled gradient of harmonic polynoimals, no
% dependency on target...
%
% Hai 07/29/20
% arbitrarily high order, Hai 07/23/22

x = r(1,:); y = r(2,:); z = r(3,:);
rho2 = x.^2+y.^2+z.^2;  % rho square

c = legendre_coefficients(n);   % lower triangular matrix of dim (n+1)*(n+1), coefficient of Legendre polynomial of degree less or equal than n, i.e. P_k(x), k<=n 
                                % row i--degree i-1, column j--coeff of x^(j-1) 
                                % for example, 1st  row (0th  degree): 1   
                                %              2nd  row (1st  degree): x
                                %              3rd  row (2nd  degree): 1/2*(-1+3*x^2) ...
                                %              11th row (10th degree): 1/256*(-63+3465*x^2-30030*x^4+90090*x^6-109395*x^8+46189*x^10)
                                % non-zero entry pattern: odd entries --> even entries --> odd entries
Cp = derivative_coefficients(c);    % cell array of dim n, 1st, 2nd, ... n-th derivative of Legendre polynomial
                                    % each cell is of dim (n+1 minus derivative order)^2
                                    % each entry is rescaled by 1/leading order to avoid coefficient blow up
                                    % for example, 
                                    %        ** 1st cell, Cp{1} is rescaled 1st derivative of P_k(x), k<=n, therefore of 1 dim less than the dim of c
                                    %          1st  row (rescaled derivative of 1st  degree): 1  
                                    %          2nd  row (rescaled derivative of 2nd  degree): 1/2*1/2*(3*2*x)
                                    %          3rd  row (rescaled derivative of 3rd  degree): 1/3*1/2*(-3+5*3*x^2)
                                    %              :
                                    %          10th row (rescaled derivative of 10th degree): 1/10*1/256*(3465*2*x-30030*4*x^3+90090*6*x^5-109395*8*x^7+46189*10*x^9) 
                                    %
                                    %        ** 2nd cell, Cp{2} is rescaled 2nd derivative of P_k(x), k<=n, therefore of 2 dim less than the dim of c
                                    %          1st  row (rescaled 2nd derivative of 2nd  degree): 1/2*1/2*(3*2)
                                    %          2nd  row (rescaled 2nd derivative of 3rd  degree): 1/2*1/3*1/2*(5*3*2*x)
                                    %          3rd  row (rescaled 2nd derivative of 4th  degree): 1/3*1/4*1/8*(-30*2*1+35*4*3*x^2)
                                    %              :
                                    %          9th  row (rescaled 2nd derivative of 4th  degree): 1/9*1/10*1/256*(3465*2-30030*4*3*x^2+90090*6*5*x^4-109395*8*7*x^6+46189*10*9*x^8) 
                                    %
                                    %        ** last cell, Cp{n} is rescaled n-nd derivative of P_n(x)
                                    %          only 1 entry: c(end,end)
% Cp = rescaled_derivative_coefficients(c);

% partial derivative of Spherical harmonics: rho^n*cos(k*phi)*sin^k(theta)*P^(k)_n(cos(theta)), k=1,...,n
%                                            when written in terms of x,y,z, they share the format 
%                                            (2d harmonics) * (derivative of Legendre(z, then (x^2+y^2+z^2) to some power to make all term of the same order))
% compute rho^k*cos(k*phi)*sin^k(theta), and its partial derivative. only depends on x & y. Actually, they are 2d harmonics, or real part of complex monomial
Hxy = real(bsxfun(@power,x+1i*y,(1:n)'));                                   % compute x, x^2-y^2, x^3-3*x*y^2, x^4+y^4-6*x^2*y^2, ...
Hxy_x = bsxfun(@times,real(bsxfun(@power,x+1i*y,(0:n-1)')),(1:n)');         % compute partial derivative of Hxy w.r.t x
Hxy_y = bsxfun(@times,-imag(bsxfun(@power,x+1i*y,(0:n-1)')),(1:n)');        % compute partial derivative of Hxy w.r.t. y

for k=1:n 
    % partial derivative of spherical harmonics w.r.t. x y z will be stored in fx{m}{j}, fy{m}{j}, fz{m}{j}, m is the m-th order basis, j=1,...,m linearly independent ones
    % (2d harmonics) * (derivative of Legendre(z, then (x^2+y^2+z^2) to some power to make all term of the same order))
    % k to keep track of length of Legendre part, later nl = numel(cpnk). 
    % (n-k)th off diagonal, from lower left most entry fx{n}{1} to diagonal... 

    % Legendre(z) part from the 2nd factor, i.e. (derivative of Legendre(z, then (x^2+y^2+z^2) to some power to make all term of the same order))
    % each term of Legendre(z) will be multiplied by (x^2+y^2+z^2)^some power to be of the same degree in x y z
    if mod(n,2)==0 % when Legendre poly order is even 
        cpnk = Cp{k}(end,(2-mod(k-1,2)):2:end);                                         % does not matter, num of nonzero elements matters
        Z = bsxfun(@power,z,(mod(k,2):2:n-k)');                                         % various power of z, either z^0, z^2, ... z^(n-k), or z^1, z^3, ... z^(n-k), depending on mod(k,2)
        Zp = bsxfun(@times,bsxfun(@power,z,(mod(k,2)-1:2:n-k-1)'),(mod(k,2):2:n-k)');   % derivative of various power of z, i.e. Z with 1 less power, multiplied by [0,2,...] or [1,3,...]
        if   mod(k,2)==0, Zp(:,z==0) = 0; end                                           % if mod(k,2)=0, mod(k,2)-1 in Zp computation needs to be modified to 0 for z=0 entries 
    else % when Legendre poly order is odd
        cpnk = Cp{k}(end,(2-mod(k,2)):2:end); 
        Z = bsxfun(@power,z,(mod(k-1,2):2:n-k)');
        Zp = bsxfun(@times,bsxfun(@power,z,(mod(k-1,2)-1:2:n-k-1)'),(mod(k-1,2):2:n-k)');
        if mod(k-1,2)==0, Zp(:,z==0) = 0; end
    end
    % each term of Legendre(z) and (x^2+y^2+z^2)^some power are multiplied in the following way
    % each row is one term z^alpha * (x^2+y^2+z^2)^beta, such that alpha+2*beta = n-k, also their partial derivatives
    nl = numel(cpnk);                                                              % num of nonzero coefficients, this num is "monotonically decreasing", as a result of the design of the loop
    Rho2 = bsxfun(@power,rho2,(0:nl-1)');                                          % various power of rho square
    Rho2_x = bsxfun(@times,2*x,(0:nl-1)').*[zeros(size(x));Rho2(1:end-1,:)];       % partial derivative of Rho2 w.r.t x
    Rho2_y = bsxfun(@times,2*y,(0:nl-1)').*[zeros(size(y));Rho2(1:end-1,:)];       % partial derivative of Rho2 w.r.t y
    Rho2_z = bsxfun(@times,2*z,(0:nl-1)').*[zeros(size(z));Rho2(1:end-1,:)];       % partial derivative of Rho2 w.r.t z
    ZRho2 = Z.*Rho2(end:-1:1,:);                                                    % for z^(j-k+1)*rho^(n-j-1)
    ZpRho2 = Zp.*Rho2(end:-1:1,:);                                                  % for (j-k+1)*z^(j-k)*rho^(n-j-1)
    ZRho2_x = Z.*Rho2_x(end:-1:1,:);                                                % for z^(j-k+1)*d(rho^(n-j-1))/dx
    ZRho2_y = Z.*Rho2_y(end:-1:1,:);                                                % for z^(j-k+1)*d(rho^(n-j-1))/dy
    ZRho2_z = Z.*Rho2_z(end:-1:1,:);                                                % for z^(j-k+1)*d(rho^(n-j-1))/dy
    % for example, say, P_5(z)=1/8*(15*x-70*x^3+63*x^5), then we need to know P^(1)_5, P^(2)_5, P^(3)_5, P^(4)_5, P^(5)_5
    %              j=1, then (degree 1 2d harmonics) * (degree 5-1 factor P^(1)_5=[z^0*(x^2+y^2+z^2)^2, z^2*(x^2+y^2+z^2), z^4]*some coefficient )
    %              j=2, then (degree 2 2d harmonics) * (degree 5-2 factor P^(2)_5=[z*(x^2+y^2+z^2), z^3]*some coefficient )

    % for given k, j = k, k-1, ..., 1,  
    % Cp{j}(end-(k-j),:) contains coefficient of j-th order derivative of all Legendre P_(x) of degree >= k 
    for j=k:-1:1
        m = n-(k-j);    % approx order
        if mod(m,2)==0
            cpmj = Cp{j}(end-(k-j),(2-mod(j-1,2)):2:(end-(k-j)));
        else
            cpmj = Cp{j}(end-(k-j),(2-mod(j,2)):2:(end-(k-j))); 
        end
        % Legendre part
        lj = cpmj*ZRho2;
        lj_x = cpmj*ZRho2_x;
        lj_y = cpmj*ZRho2_y;
        lj_z = cpmj*ZRho2_z + cpmj*ZpRho2;
        % 2d harmonic part
        hj = Hxy(j,:);
        hj_x = Hxy_x(j,:);
        hj_y = Hxy_y(j,:);
        % combine together
        fx{m}{j} = hj_x.*lj+hj.*lj_x; 
        fy{m}{j} = hj_y.*lj+hj.*lj_y; 
        fz{m}{j} = hj.*lj_z;
        
    end
end

end