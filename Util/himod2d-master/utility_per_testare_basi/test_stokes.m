addpath ../Core
x=linspace(-1,1,200);
m=16
BS=new_modal_basis_legendre(m-2,x,'rob','rob')
[BV,BVD]=new_modal_basis_legendre(m,x,'dir','dir')

delta=zeros(m,m);
for i=1:m-2
	for k=1:m
			aus = BS(:,i).*BVD(:,k);
			delta(i,k) = integrami(aus,x);				
	end
end
delta