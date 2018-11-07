function [A,b]=DirichletInflowREL(A,b,dir_fourier,m,fedof)

gamma=1e30*max( max( max( abs(A) ), max(abs(b)) ) );

for k=1:m
    if(k<=2)
        A(k,k)  = gamma;
        b(k)                  = gamma*dir_fourier(k);
    else
        A( (k-1-2)*fedof +2+ 1, (k-1-2)*fedof +2+1)  = gamma;
        b( (k-1-2)*fedof +2+ 1 )                  = gamma*dir_fourier(k);
    end
end