function M = matrice_massa(griglia,base,rho)
% M = matrice_massa(griglia,base,rho)
%
% Calcolo della matrice di massa.
%   M_ij = \int rho phi_i phi_j dx
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% rho: stringa con il nome della funzione per la valutazione della
%      densità. Tale funzione deve avere l'interfaccia seguente:
%         rho = rho_fun(x)
%      dove x è un vettore con le coordinate dei punti in cui è
%      richiesta la densità e rho è il vettore (della medesima
%      dimensione di x) che raccoglie i corrispondenti valori di
%      densità.
%      N.B. e' anche possibile utilizzare un puntatore a funzione.


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ciclo sugli elementi per assemblare M
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  lel = griglia.lel;
  N_r = base.N_r;
  
  % allocazione di M in formato sparso
  M = sparse(griglia.dim,griglia.dim);
  
  % ciclo sugli elementi
  for ie = 1:griglia.ne
     % calcolo del contributo elementare
     Mk = matrice_massa_loc(griglia.xnodi(lel(ie,:)), ...
                            base.PHI,base.gradPHI, ...
                            base.pesiGauss, ...
                            base.N_r,base.M, ...
                            rho);

     % inserimento in M del contributo elementare (assemblaggio)
     for i=1:N_r
       for j=1:N_r
         M(lel(ie,i),lel(ie,j)) = M(lel(ie,i),lel(ie,j)) + Mk(i,j);
       end
     end
  end
  
return

