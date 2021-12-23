function c = legendre_coefficients(n)
% return legendre polynomial coefficients of all order up to n, based on
% P(i,x) = (2*i-1)/i*x*P(i-1,x)-(i-1)/i*P(i-2,x)
% 
% 07/22/20, Hai

% initialize coefficient matrix
c = zeros(n+1); c(1,1) = 1.0; c(2,2) = 1.0;

% ith row corresponds to (i-1)th degree polynomial coefficient (i coefficients including constant)
for i=2:n
    c(i+1,1:i-1) = -(i-1) * c(i-1,1:i-1)/i; % -(i-1)/i*P(i-2,x) part
    c(i+1,2:i+1) = c(i+1,2:i+1) + (2*i-1) * c(i,1:i)/i; % (2*i-1)/i*x*P(i-1,x) part
end
 
end