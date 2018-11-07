function S = termine_reazione_lumping(griglia,base,sigma)
% S = termine_reazione_lumping(griglia,base,sigma)
%
% Calcolo della matrice associata al termine di reazione (matrice di
% massa), con condensazione (lumping).
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% sigma: stringa con il nome della funzione per la valutazione del
%        coefficiente di reazione. Si veda termine_reazione per i
%        dettagli.

  % matrice di massa in assenza di condensazione
  S = termine_reazione(griglia,base,sigma);
  
  % condensazione: somma per righe su S
  S = spdiags(sum(S,2),0,griglia.dim,griglia.dim);
  
return

