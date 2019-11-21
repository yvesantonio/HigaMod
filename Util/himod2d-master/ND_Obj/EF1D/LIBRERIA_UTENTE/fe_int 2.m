function Int = fe_int(fun,griglia,base)
% Int = fe_int(fun,griglia,base)
%
% Calcolo dell'integrale di una funzione fun appartenente allo spazio
% di elementi finiti descritto da griglia e base.
%
% fun: vettore con i valori nodali della funzione integranda
% griglia: rappresentazione della griglia di calcolo secondo lo schema
%       descritto in struttura_griglia
% base: base dello spazio di elementi finiti riferita all'elemento
%       campione (funzioni di forma): si veda struttura_base

 Int = 0; % inizializzazione
 for ie=1:griglia.ne % ciclo sugli elementi

   % valutazione della funzione integranda nei nodi di quadratura
   fGauss = zeros(1,base.M);
   for i=1:base.N_r
     fGauss = fGauss + fun( griglia.lel(ie,i) )*base.PHI(i,:);
   end

   % trasformazione di coordinate isoparametrica
   xnodi = griglia.xnodi(griglia.lel(ie,:)); % nodi dell'elemento
   dx_dcsi = xnodi * base.gradPHI;
   % calcolo dell'integrale
   for l=1:base.M
     Int = Int + base.pesiGauss(l)*fGauss(l)*dx_dcsi(l);
   end

 end

return

