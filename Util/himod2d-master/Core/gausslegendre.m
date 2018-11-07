function [nqn,qnodes,w] = gausslegendre(n) 
%
%GAUSSLEGENDRE    One-dimensional Gauss-Legendre quadrature rule.
%
%    [NQN,QNODES,W]=GAUSSLEGENDRE(NP) returns the
%    NP quadrature nodes QNODES and the NP weights
%    of the Gauss-Legendre quadrature rule in (0,1).
 
%   Gervasio Paola, 1997.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $

nqn = n;
if n <= 1    
   x=1/2;w=1;    
   return 
end 
jac = zeros(n);   k=[1:n-1]; 
v   = (k)./(sqrt(4*(k.^2)-1)); 
jac = jac+diag(v,1)+diag(v,-1); 
[w,x]=eig(jac); 
norm2 = sqrt(diag(w'*w));
w = (2*w(1,:)'.^2)./norm2;
x = diag(x); 
 
for k=1:n-1 
   l=k;  a=x(k);    
   for j=k+1:n       
      if x(j)<a          
         l=j;          
         a=x(j);          
      end       
   end    
   if l ~= k       
      t = x(k);  x(k) = x(l);  x(l) = t;       
      t = w(k);  w(k) = w(l);  w(l) = t;       
   end    
end 


qnodes = 1/2*(x+1);
w = w/2;
return
