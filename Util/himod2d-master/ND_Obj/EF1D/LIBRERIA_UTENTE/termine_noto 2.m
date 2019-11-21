function b = termine_noto(griglia,base,fvol)
% b = termine_noto(griglia,base,fvol)
%
% Calcolo del termine noto:
%   b_i = \int f phi_i dx
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% fvol: stringa con il nome della funzione per la valutazione del
%       termine noto. Tale funzione deve avere
%       l'interfaccia seguente:
%          fvol = fvol_fun(x)
%       dove x è un vettore con le coordinate dei punti in cui è
%       richiesto il coefficiente e fvol è il vettore (della medesima
%       dimensione di x) che raccoglie i corrispondenti valori della
%       forzante.
%       N.B. e' anche possibile utilizzare un puntatore a funzione.


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ciclo sugli elementi per assemblare b
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  lel = griglia.lel;
  N_r = base.N_r;
  
  % allocazione di b
  b = zeros(griglia.dim,1);
  
  % ciclo sugli elementi
  for ie = 1:griglia.ne
     % calcolo del contributo elementare
     bk = termine_noto_loc(griglia.xnodi(lel(ie,:)), ...
                           base.PHI,base.gradPHI, ...
                           base.pesiGauss, ...
                           base.N_r,base.M, ...
                           fvol);
     
     % inserimento in b del contributo elementare (assemblaggio)
     for i=1:N_r
        b(lel(ie,i)) = b(lel(ie,i)) + bk(i);
     end
  end
  
return

