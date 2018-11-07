function base_ef = valuta_base_ef(griglia,base,x)
% base_ef = valuta_base_ef(griglia,base,x)
%
% Costruzione della base lagrangiana dello spazio di elementi definito
% da base e griglia.
% 
% griglia: rappresentazione della griglia di calcolo secondo lo schema
%       descritto in struttura_griglia
% base: base dello spazio di elementi finiti riferita all'elemento
%       campione (funzioni di forma): si veda struttura_base
% x: vettore di punti nei quali è richiesta la valutazione delle
%    funzioni di base
%
% base_ef: valutazioni delle funzioni di base nei punti x; più
%          precisamente, base_ef ha dimensioni [griglia.dim,lebgth(x)]
%          e ciascuna riga corrisponde ad una funzione di base.
%
% NOTA: Il grado dello spazio di elementi finiti indicato in base deve
% essere coerente con il numero di nodi per ciescun elemento idicato
% in griglia.

 % assicuriamoci che x sia un vettore riga
 if(size(x,1)>size(x,2))
   x = x';
 end
 base_ef = zeros(griglia.dim,length(x));

 for ie=1:griglia.ne % ciclo sugli elementi

   % estremi dell'elemento
   XX = griglia.xnodi(griglia.lel(ie,:));
   X1 = min(XX);
   X2 = max(XX);

   % individuare i punti di x che appartengono all'intervallo i-esimo
   idx = find((X1<=x).*(x<=X2));
   if(~isempty(idx))
     xx = x(idx);

     % valutazione delle funzioni di base dell'elemento
     for i=1:base.N_r
       base_ef(griglia.lel(ie,:),idx) = LagrPoli(XX,xx);
     end
   end

 end

return
