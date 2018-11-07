function [A,b]=DirichletInflow(A,b,dir_fourier,m,fedof)

gamma=1e30*max( max( max( abs(A) ), max(abs(b)) ) );

for k=1:m
    A( (k-1)*fedof + 1, (k-1)*fedof + 1)  = gamma;
    b( (k-1)*fedof + 1 )                  = gamma*dir_fourier(k);
end

% for k=1:m
%      A( (k-1)*fedof + 1 , :)=zeros(size(A( k*fedof , :)));
%      A( (k-1)*fedof + 1 , (k-1)*fedof + 1 )  = 1;
%      b( (k-1)*fedof + 1  )           = dir_fourier(k);
% end