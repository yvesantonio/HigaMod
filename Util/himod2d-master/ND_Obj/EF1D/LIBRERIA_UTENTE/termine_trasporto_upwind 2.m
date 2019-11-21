function C = termine_trasporto_upwind(griglia,base,a,mu)
% C = termine_trasporto_upwind(griglia,base,a,mu)
%
% Calcolo della matrice associata al trasporto (forma conservativa)
% con stabilizzazione "upwind".
%
% Si noti che per applicare la stabilizzazione e' necessario conoscere
% sia il coefficiente di trasporto che quello di diffusione, necessari
% per determinare il numero di Peclet.
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta
%          con 'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% a, mu: funzioni per la valutazione dei coefficienti del problema, si
%        vedano le funzioni termine_trasporto e termine_diffusione per
%        maggiori dettagli.


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
     Kk = termine_diffnum_loc(griglia.xnodi(lel(ie,:)), ...
                                base.PHI,base.gradPHI, ...
                                base.pesiGauss, ...
				base.N_r,base.M, ...
                                a,mu,'upwind');
     
     % inserimento in C del contributo elementare (assemblaggio)
     for i=1:N_r
       for j=1:N_r
         C(lel(ie,i),lel(ie,j)) = C(lel(ie,i),lel(ie,j)) + Ck(i,j) ...
	                         + Kk(i,j);
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

