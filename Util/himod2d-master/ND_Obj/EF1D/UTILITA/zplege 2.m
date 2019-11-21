function [x,w] = zplege(n)
% [x,w] = zplege(n)
%
% Calcola nodi e pesi di quadratura della formula di Gauss-Legendre a
% n nodi sull'intervallo di riferimento [-1,1], n>=1.
%
% Tratto da "A. Quarteroni, R. Sacco, F. Saleri, Matematica Numerica,
% Springer Verlag, 2008"; si veda anche "W. Gautschi, Orthogonal
% Polynomials: Computation and Approximation, Oxford University Press,
% 2004".

 if(n==1), 
   x = 0;
   w = 2;
   return; 
 end
 
 [a,b]=coeflege(n); 
 JacM=diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);
 [w,x]=eig(JacM); 
 x=diag(x); 
 scal=2; 
 w=w(1,:)'.^2*scal;
 [x,ind]=sort(x); 
 w=w(ind);

return
