function fourier_coeff = projection(basis,func,weight)

m=size(basis,2);
fourier_coeff=zeros(m,1);
for i = 1 : m
        fourier_coeff(i) = integrate( (func.*basis(:,i))', weight );
end
