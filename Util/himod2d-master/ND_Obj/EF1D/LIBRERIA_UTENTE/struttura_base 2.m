function [base,csiGauss] = struttura_base(r)
% [base,csiGauss] = struttura_base(r)
%
% Costruzione della struttura descrivente la base di elementi finiti
% in funzione del grado r. Tale struttura contiene tutte le
% informazioni relative alle funzioni di forma PHI definite
% sull'elemento di riferimento [-1 1]. In particolare sono presenti i
% seguenti campi:
%
% N_r: numero di funzioni di forma
% M: numero di nodi di Gauss (per ciascun elemento)
% PHI: valutazione delle funzioni di forma sui nodi di Gauss
%      nell'elemento di riferimento
% gradPHI: valutazione delle derivate delle funzioni di forma sui nodi 
%          di Gauss nell'elemento di riferimento
% csinodi: nodi equispaziati nell'elemento di riferimento
% csiGauss: nodi di integrazione gaussiana nell'elemento di riferimento
% pesiGauss: pesi di integrazione gaussiana
%
% N.B. oltre alla struttura base, vengono restituiti anche i punti di 
% Gauss csiGauss.

 base = struct('N_r',[], ...
               'M',[], ...
               'PHI',[], ...
               'gradPHI',[], ...
               'csinodi',[], ...
               'csiGauss',[], ...
               'pesiGauss',[]);
 
 base.N_r = r+1;  % una funzione di forma per ogni nodo
 base.M = r+1;  % numero di nodi
 N_r = base.N_r;
 M = base.M;
 
 % coordinate dei nodi sull'intervallo di riferimento
 base.csinodi = linspace(-1,1,N_r);
 
 % nodi e pesi per la qudratura gaussiana
 [base.csiGauss, base.pesiGauss] = zplege(base.M);
 base.csiGauss  = base.csiGauss';  % vettore riga
 base.pesiGauss = base.pesiGauss';
 
 % valutazione dei polinomi nodali
 [base.PHI, base.gradPHI] = LagrPoli(base.csinodi,base.csiGauss);
 
return

