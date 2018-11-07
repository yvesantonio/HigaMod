function Int = trapint(fun,a,b,N)
% Int = trapint(fun,a,b,N)
%
% Calcolo dell'integrale di fun su [a,b] usando il metodo dei trapezi
% con N suddivisioni.
%
% fun: stringa con il nome della funzione integranda
% a,b: estremi dell'intervallo di integrazione
% N: numero di suddivisioni

 % ampiezza dei sottointervalli
 h = (b-a)/N;
 % estremi dei sottointervalli
 xnod = [ a : h : b ];

 Int = 0; % inizializzazione
 for j=1:N % ciclo sui sottointervalli
   % estremi del j-esimo sottointervallo
   xA = xnod(j);
   xB = xnod(j+1);
   % valutazione della funzione integranda
   fA = feval(fun,xA);
   fB = feval(fun,xB);
   % calcolo dell'integrale: area del trapezio
   Int = Int + (fA+fB)*h/2;
 end

end

