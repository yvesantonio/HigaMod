function [errL2,errH1] = norme_errore(griglia,base,u,uexfun,graduexfun)
% [errL2,errH1] = norme_errore(griglia,base,u,uexfun,graduexfun)
%
% Funzione per il calcolo delle norme dell'errore.
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% u: vettore rappresentante la soluzione
% uexfun: stringa con il nome della funzione che valuta la soluzione
%         esatta
% graduexfun: stringa con il nome della funzione che valuta il gradiente
%             della soluzione esatta
%
% errL2: errore norma L2
% errH1: errore norma H1

% Per non intaccare l'ordine di convergenza del metodo degli elementi
% finiti, impieghiamo un ordine di esattezza pari a 2r+2, essendo r il
% grado degli elementi finiti. In pratica occorre considerare un nodo
% di quadratura in piu' rispetto a quanto fatto per il calcolo delle
% matrici. 

 % u vettore colonna
 if(size(u,1)<size(u,2))
   u = u';
 end
 
 Merr = base.M+1; % nuovo numero di nodi di quadratura
 [csiGauss, pesiGauss] = zplege(Merr);
 csiGauss = csiGauss';  % vettore riga
 pesiGauss = pesiGauss';
 
 % valuto le funzioni di base nei nuovi nodi di quadratura
 [PHI, gradPHI] = LagrPoli(base.csinodi,csiGauss);

 % valuto gli integrali elemento per elemento
 err_u = 0;
 err_gradu = 0;
 for ie = 1:griglia.ne
   [err_uk,err_graduk] = ...
        errori_elemento(griglia.xnodi(griglia.lel(ie,:)), ...
                                      u(griglia.lel(ie,:)), ...
                                      PHI,gradPHI, pesiGauss, ...
                                      uexfun,graduexfun);
   % accumulazione degli errori
   err_u = err_u+err_uk;
   err_gradu = err_gradu+err_graduk;
 end
 err_u = sqrt(err_u);
 err_gradu = sqrt(err_gradu);
 
 errL2 = err_u;
 errH1 = sqrt(err_u^2 + err_gradu^2);
 
return


function [err_uk,err_graduk] = errori_elemento(xnodi,u,PHI,gradPHI, ...
                                               pesi,uexfun,graduexfun);
 
 % per maggiori dettagli si vedano anche i commenti nel file
 % UTILITA/termine_diffusione_loc.m

 % soluzione elementi finiti
 xGauss = xnodi * PHI;
 uGauss = u' * PHI;
 dx_dcsi = xnodi * gradPHI;
 graduGauss = (u' * gradPHI)./dx_dcsi;
 
 % soluzione esatta
 uex = feval(uexfun,xGauss);
 grad_uex = feval(graduexfun,xGauss);
 
 % errori
 err_uk = sum(pesi.*((uex-uGauss).^2).*dx_dcsi);
 err_graduk = sum(pesi.*((grad_uex-graduGauss).^2).*dx_dcsi);
 
return
 
