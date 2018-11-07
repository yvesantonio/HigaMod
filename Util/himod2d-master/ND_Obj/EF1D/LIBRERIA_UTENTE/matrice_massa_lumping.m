function M = matrice_massa_lumping(griglia,base,rho)
% M = matrice_massa_lumping(griglia,base,rho)
%
% Calcolo della matrice di massa con condensazione (lumping).
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% rho: stringa con il nome della funzione per la valutazione della
%      densità. Si veda la funzione matrice_massa per ulteriori
%      dettagli.

  % matrice di massa in assenza di condensazione
  M = matrice_massa(griglia,base,rho)
  
  % condensazione: somma per righe
  M = spdiags(sum(M,2),0,griglia.dim,griglia.dim);
  
return

