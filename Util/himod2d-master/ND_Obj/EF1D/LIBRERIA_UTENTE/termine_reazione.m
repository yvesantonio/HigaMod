function S = termine_reazione(griglia,base,sigma)
% S = termine_reazione(griglia,base,sigma)
%
% Calcolo della matrice associata al termine di reazione (matrice di
% massa)
%   S_ij = \int sigma phi_i phi_j dx
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% sigma: stringa con il nome della funzione per la valutazione del
%        coefficiente di reazione. Tale funzione deve avere
%        l'interfaccia seguente:
%           sigma = sigma_fun(x)
%        dove x è un vettore con le coordinate dei punti in cui è
%        richiesto il coefficiente e sigma è il vettore (della
%        medesima dimensione di x) che raccoglie i corrispondenti
%        valori del coefficiente di reazione.
%        N.B. e' anche possibile utilizzare un puntatore a funzione.


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ciclo sugli elementi per assemblare S
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  lel = griglia.lel;
  N_r = base.N_r;
  
  % allocazione di S in formato sparso
  S = sparse(griglia.dim,griglia.dim);
  
  % ciclo sugli elementi
  for ie = 1:griglia.ne
     % calcolo del contributo elementare
     Sk = termine_reazione_loc(griglia.xnodi(lel(ie,:)), ...
                               base.PHI,base.gradPHI, ...
                               base.pesiGauss, ...
                               base.N_r,base.M, ...
                               sigma);

     % inserimento in S del contributo elementare (assemblaggio)
     for i=1:N_r
       for j=1:N_r
         S(lel(ie,i),lel(ie,j)) = S(lel(ie,i),lel(ie,j)) + Sk(i,j);
       end
     end
  end
  
return

