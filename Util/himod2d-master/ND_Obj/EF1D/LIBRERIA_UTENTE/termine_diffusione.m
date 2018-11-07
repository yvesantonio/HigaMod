function K = termine_diffusione(griglia,base,mu)
% K = termine_diffusione(griglia,base,mu)
%
% Calcolo della matrice di rigidezza:
%   K_ij = \int mu grad(phi_i) grad(phi_j) dx
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% mu: stringa con il nome della funzione per la valutazione del
%     coefficiente di diffusione. Tale funzione deve avere
%     l'interfaccia seguente:
%        mu = mu_fun(x)
%     dove x è un vettore con le coordinate dei punti in cui è
%     richiesto il coefficiente e mu è il vettore (della medesima
%     dimensione di x) che raccoglie i corrispondenti valori del
%     coefficiente di diffusione.
%     N.B. e' anche possibile utilizzare un puntatore funzione.


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ciclo sugli elementi per assemblare K
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  lel = griglia.lel;
  N_r = base.N_r;
  
  % allocazione di K in formato sparso
  K = sparse(griglia.dim,griglia.dim);
  
  % ciclo sugli elementi
  for ie = 1:griglia.ne
     % calcolo del contributo elementare
     Kk = termine_diffusione_loc(griglia.xnodi(lel(ie,:)), ...
                                 base.PHI,base.gradPHI, ...
                                 base.pesiGauss, ...
                                 base.N_r,base.M, ...
                                 mu);

     % inserimento in K del contributo elementare (assemblaggio)
     for i=1:N_r
       for j=1:N_r
         K(lel(ie,i),lel(ie,j)) = K(lel(ie,i),lel(ie,j)) + Kk(i,j);
       end
     end
  end
  
return

