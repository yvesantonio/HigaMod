function B=bern(i,n,u)
    B=factorial(n)/factorial(i)/factorial(n-i)*u.^i.*(1-u).^(n-i);
end
