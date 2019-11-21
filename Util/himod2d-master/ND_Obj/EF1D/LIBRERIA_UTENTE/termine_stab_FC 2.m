function [FCA,FCb] = termine_stab_FC(griglia,base, ...
                                     mu,gradmu,a,diva,sigma,fvol, ...
			             d,rho)
% [FCA FCb] = termine_stab_FC(griglia,base,mu,gradmu,a,diva,sigma, ...
%                             fvol,d,rho)
%
% Calcolo dei termini correttivi per stabilizzazione fortemente
% consistente (matrice FCA e termine noto FCb). Si impiega la forma
% conservativa dell'equazione.
%
% I due termini FCA e FCb devono essere sommati rispettivamente alla
% matrice ed al termine noto corrispondenti alla discretizzazione
% senza stabilizzazione del problema.
%
% Si noti che per la stabilizzazione fortemente consistente è
% necessario specificare informazioni aggiuntive sulle derivate
% seconde delle funzioni di base e sulle derivate prime dei
% coefficienti.
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta
%          con 'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma.
%       IMPORTANTE: rispetto alla struttura fornita da
%       'struttura_base' è necessario aggiungere un campo 'lapPHI' che
%       contenga le derivate seconde delle funzioni di base (usare ad
%       esempio la funzione 'aggiunta_lapPHI').
% mu, a, sigma, fvol: valgono le indicazioni date in
%                     termine_diffusione, termine_trasporto ecc.
% gradmu: stringa con il nome della funzione che valuta il gradiente
%         del coefficiente di diffusione mu
%         N.B. e' anche possibile utilizzare un puntatore a funzione.
% diva: stringa con il nome della funzione che valuta la divergenza
%       del coefficiente di trasporto a
%       N.B. e' anche possibile utilizzare un puntatore a funzione.
% d, rho: parametri per la stabilizzazione fortemente consistente.


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ciclo sugli elementi per assemblare FCA, FCb
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  lel = griglia.lel;
  N_r = base.N_r;
  
  % allocazione di FCA e FCb in formato sparso
  FCA = sparse(griglia.dim,griglia.dim);
  FCb = sparse(griglia.dim,1);
  
  % ciclo sugli elementi
  for ie = 1:griglia.ne
     % calcolo del contributo elementare
     [FCAk,FCbk] = termine_stab_FC_loc(griglia.xnodi(lel(ie,:)), ...
                                  base.PHI,base.gradPHI,base.lapPHI, ...
                                  base.pesiGauss,base.N_r,base.M, ...
                                  mu,gradmu,a,diva,sigma,fvol, ...
				  d,rho);
     
     % inserimento in FCA e FCb dei contributi elementari (assemblaggio)
     for i=1:N_r
       for j=1:N_r
	  FCA(lel(ie,i),lel(ie,j)) = FCA(lel(ie,i),lel(ie,j)) + ...
	                             FCAk(i,j);
       end
       FCb(lel(ie,i)) = FCb(lel(ie,i)) + FCbk(i);
     end
  end
  
return
