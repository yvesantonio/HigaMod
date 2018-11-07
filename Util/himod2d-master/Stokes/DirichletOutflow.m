function [A,b]=DirichletOutflow(A,b,dir_fourier,m,fedof)

gamma=1e30*max( max( max( abs(A) ), max(abs(b)) ) );

for k=1:m
    A( k*fedof , k*fedof )  = gamma;
    b( k*fedof  )           = gamma*dir_fourier(k);
end

% for k=1:m
%      A( k*fedof , :)=zeros(size(A( k*fedof , :)));
%      A( k*fedof , k*fedof )  = 1;
%      b( k*fedof  )           = dir_fourier(k);
% end