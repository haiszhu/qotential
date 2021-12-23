function Cp = rescaled_derivative_coefficients(c0)
% rescaled version of derivative coefficients
%
% Hai 07/22/20

[n,~] = size(c0); c=c0;
k = 1; % derivative order
while n>1
    cp = zeros(n,n); 
    for i=2:n
        cp(i,1:i-1) = (1:i-1).*c(i,2:i);
    end
    cp = cp(2:end,1:end-1);
    n = n-1;
    Cp{k} = cp;
    k = k+1;
    c = cp;
end

[n,~] = size(c0);
k = 1; % derivative order
while n>1
    cp = Cp{k};
    ce = rescaled_coeff(k);
    for i=1:n-1
        cp(i,1:i) = ce(i).*cp(i,1:i);
    end
    Cp{k} = cp;
    n = n-1;
    k = k+1;
end

% keyboard

end

function ce = rescaled_coeff(k)

if k==1
    ce = [ 1,1/3,2/3,2/5,8/15,8/21,16/7,16/9,128/45,128/55];
end
if k==2
    ce = [ 1/3,1/15,2/15,2/105,8/105,8/63,16/315,16/495,128/495];
end
if k==3
    ce = [ 1/15,1/105,2/105,2/315,8/315,8/3465,16/3465,16/6435];
end
if k==4
    ce = [ 1/105,1/945,2/945,2/315,8/10395,8/135135,16/45045];
end
if k==5
    ce = [ 1/945,1/10395,2/945,2/135135,8/135135,8/135135];
end
if k==6
    ce = [ 1/10395,1/135135,2/135135,2/675675,8/675675];
end
if k==7
    ce = [ 1/135135,1/2027025,2/2027025,2/11486475];
end
if k==8
    ce = [ 1/2027025,1/34459425,2/34459425];
end
if k==9
    ce = [ 1/34459425,1/654729075];
end
if k==10
    ce = [ 1/654729075];
end

end