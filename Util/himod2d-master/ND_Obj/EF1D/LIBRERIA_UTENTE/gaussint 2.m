function Int = gaussint(fun,a,b,N)
% Int = gaussint(fun,a,b,N)
%
% Calcolo dell'integrale di fun su [a,b] usando il metodo di Gauss a
% due punti con N suddivisioni.
%
% fun: stringa con il nome della funzione integranda
% a,b: estremi dell'intervallo di integrazione
% N: numero di suddivisioni

 % nodi e pesi di quadratura sull'intervallo [-1,1]
 csiA = -1/sqrt(3);
 csiB =  1/sqrt(3);
 wA = 1;
 wB = 1;

 % ampiezza dei sottointervalli
 h = (b-a)/N;
 % estremi dei sottointervalli
 xnod = [ a : h : b ];

 Int = 0; % inizializzazione
 for j=1:N % ciclo sui sottointervalli
   % nodi di quadratura mappati sul j-esimo sottointervallo: la
   % trasformazione  [-1,1] -> [x1,x2]  è
   %   x = h/2*csi + (x1+x2)/2
   xA = h/2*csiA + (xnod(j)+xnod(j+1))/2;
   xB = h/2*csiB + (xnod(j)+xnod(j+1))/2;
   % valutazione della funzione integranda nei nodi di quadratura
   fA = feval(fun,xA);
   fB = feval(fun,xB);
   % calcolo dell'integrale: include il fattore di scala h/2
   Int = Int + ( wA*fA + wB*fB )*h/2;
 end

end

