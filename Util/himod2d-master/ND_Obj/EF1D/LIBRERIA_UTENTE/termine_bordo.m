function [R,r] = termine_bordo(griglia,base,dati_bordo)
% [R,r] = termine_bordo(griglia,base,dati_bordo)
%
% Funzioni per il calcolo dei termini associati alle condizioni
% al contorno di 2a e 3a specie.
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% dati_bordo: struttura con i seguenti campi
%    bc: vettore di due elementi corrispondenti agli estremi
%        dell'intervallo di definizione del problema. Ciascun elemento
%        vale 1 se la corrispondente condizione è di Dirichlet, 2 se
%        di Neumann, 3 per il caso Robin.
%    gamma: vettore di due elementi con i coefficienti gamma della
%        condizione di Robin ai due estremi. Nota: è necessario
%        assegnare un valore anche in presenza condizioni di Dirichlet
%        o di Neumann; tale valore verrà poi ignorato.
%    r: vettore di due elementi con il valore del dato r per
%        condizioni al bordo di Neumann e Robin. Nota: è necessario
%        assegnare un valore anche in presenza condizioni di Diriclhet
%        o di Neumann; tale valore verrà poi ignorato.
%
% R: matrcie associata al termine di Robin, da sommare ai contributi
%    interni
% r: vettore dei contributi di bordo al termine noto del sistema lineare
%
% IMPORTANTE: per l'esatta definizione dei coefficienti gamma ed r si
% faccia riferimento al file "guida.pdf".

 % allocazioni
 R = sparse(griglia.dim,griglia.dim);
 r = zeros(griglia.dim,1);
 
 % primo nodo
 if(dati_bordo.bc(1)~=1)
   % condizione di tipo Neumann o Robin
   r(1) = -dati_bordo.r(1);
   if(dati_bordo.bc(1)~=2) % condizione di tipo Robin
     R(1,1) = dati_bordo.gamma(1);
   end
 end
 
 % secondo nodo
 if(dati_bordo.bc(2)~=1)
   % condizione di tipo Neumann o Robin
   r(end) = -dati_bordo.r(2);
   if(dati_bordo.bc(2)~=2) % condizione di tipo Robin
     R(end,end) = dati_bordo.gamma(2);
   end
 end
  
return

