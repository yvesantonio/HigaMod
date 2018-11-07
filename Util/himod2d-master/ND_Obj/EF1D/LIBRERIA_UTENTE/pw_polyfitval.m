function Pf = pw_polyfitval(Xnod,Ynod,x,n)
% Pf = pw_polyfitval(Xnod,Ynod,x,n)
%
% Interpolazione polinomiale composita
% 
% Xnod: nodi di interpolazione
% Ynod: valori della funzione da interpolare nei nodi Xnod
% x:    punti usati per la valutazione dell'interpolante
% n:    grado polinomiale
% 
% Pf: valori del polinomio interpolante nei punti x
%
% Note
% 1) I vettory Xnod e x devono essere in ordine crescete.
% 2) Il numero di nodi in Xnod deve essere consistente con il grado
% polinomiale. Per esempio, se n=2, Xnod puo' contenere 3, 5, 7, ...
% nodi, mentre non sarebbe possibile usare 2 o 4 nodi. In generale, il
% numero di nodi deve essere  1+n*k  per un qualunque intero k.

 % Numero di intervalli
 Nint = (size(Xnod,2)-1)/n;
 % Il numero di intervalli deve essere intero
 if(floor(Nint)~=Nint)
   error([ ...
 'Numero di nodi non compatibile con il grado polinomiale richiesto.']);
 end
 
 Pf = zeros(1,size(x,2));

 for i=1:Nint % ciclo sui sottointervalli

   % inidici iniziali e finali per l'intervallo i-esimo
   is = 1 + (i-1)*n;
   ie = is + n;
   % nodi dell'intervallo i-esimo
   X = Xnod(is:ie);
   Y = Ynod(is:ie);
   % individuare i punti di x che appartengono all'intervallo i-esimo
   idx = find(and(X(1)<=x,x<=X(end)));
   if(~isempty(idx))
     xx = x(idx);

     % interpolazione polinomiale sull'intervallo
     [P,s,mu] = polyfit(X,Y,n);
     Pf(idx) = polyval(P,xx,[],mu);
   end

 end

return
