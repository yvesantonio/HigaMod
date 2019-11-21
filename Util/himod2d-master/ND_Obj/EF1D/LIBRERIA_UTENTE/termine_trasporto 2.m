function C = termine_trasporto(griglia,base,a)
% C = termine_trasporto(griglia,base,a)
%
% Calcolo della matrice associata al trasporto (forma conservativa).
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta
%          con 'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% a: stringa con il nome della funzione che valuta il coefficiente
%    di trasporto. Tale funzione deve avere
%    l'interfaccia seguente:
%       a = a_fun(x)
%    dove x è un vettore con le coordinate dei punti in cui è
%    richiesto il coefficiente e a è il vettore (della medesima
%    dimensione di x) che raccoglie i corrispondenti valori del
%    coefficiente di trasporto.
%    N.B. e' anche possibile utilizzare un puntatore a funzione.


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ciclo sugli elementi per assemblare C
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  lel = griglia.lel;
  N_r = base.N_r;
  
  % allocazione di C in formato sparso
  C = sparse(griglia.dim,griglia.dim);
  
  % ciclo sugli elementi
  for ie = 1:griglia.ne
     % calcolo del contributo elementare
     Ck = termine_trasporto_loc(griglia.xnodi(lel(ie,:)), ...
                                base.PHI,base.gradPHI, ...
                                base.pesiGauss, ...
                                base.N_r,base.M, ...
                                a);
     
     % inserimento in C del contributo elementare (assemblaggio)
     for i=1:N_r
       for j=1:N_r
         C(lel(ie,i),lel(ie,j)) = C(lel(ie,i),lel(ie,j)) + Ck(i,j);
       end
     end
  end
  
  % introduzione dei termini di flusso in uscita: in corrispondenza 
  % della prima e dell'ultima fuzione di base
  a_estremi = feval(a,griglia.xnodi([1 griglia.dim]));
  apiu = max([0,a_estremi(1)*griglia.normali(1)]);
  C(1,1) = C(1,1) + apiu;
  apiu = max([0,a_estremi(2)*griglia.normali(2)]);
  C(end,end) = C(end,end) + apiu;
  
return

