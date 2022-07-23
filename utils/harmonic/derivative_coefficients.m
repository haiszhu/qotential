function Cp = derivative_coefficients(c0)
% derivative coefficients with rescaling factor 1/n!
%
% without rescaling, for n=16, coefficient can blow up tp more than 1e+16
%
% Hai 07/21/22

[n,~] = size(c0); c=c0;
k = 1; % derivative order
while n>1

    % coefficient of the derivative of legendre polynomial
    % cp = bsxfun(@times,c,linspace(0,n-1,n)); cp = cp(2:end,2:end);

    % compare to above, seems more efficient to do it elementwise 
    cp = zeros(n,n); 
    for i=2:n
        % power becomes coefficient
        % recale back according to largest coefficient being multiplied...
        cp(i,1:i-1) = ((1:i-1)/(i-1)).*c(i,2:i);
    end
    cp = cp(2:end,1:end-1);
    
    % store and continue to next order derivative
    n = n-1;
    Cp{k} = cp;
    k = k+1;
    c = cp;
end

end