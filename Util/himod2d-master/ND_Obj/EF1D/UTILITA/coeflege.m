function [a,b] = coeflege(n)
% [a,b] = coeflege(n)
%
% Calcola i coefficienti del polinomio di Legendre di grado n, n>=2.
%
% Tratto da "A. Quarteroni, R. Sacco, F. Saleri, Matematica Numerica,
% Springer Verlag, 2008"; si veda anche "W. Gautschi, Orthogonal
% Polynomials: Computation and Approximation, Oxford University Press,
% 2004".

 if (n <= 1), 
   error(' n deve essere > 1 ');
 end 
 
 a=zeros(n,1);
 b=a;
 b(1)=2;
 for k=2:n, 
   b(k)=1/(4-1/(k-1)^2); 
 end

return
