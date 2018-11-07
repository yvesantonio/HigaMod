function griglia = struttura_griglia(xnodi,r)
% griglia = struttura_griglia(xnodi,r)
%
% Costruzione della struttura descrivente la griglia di calcolo. Tale
% struttura ha i campi seguenti:
%
% xnodi: coordinate dei nodi della griglia
% lel: per ciascun elemento riporta gli indici dei nodi appartenenti.
%      Gli indici vanno da 1 fino a size(xnodi,2). Gli indici dei nodi
%      dell'elemento i-esimo sono dati da:
%         i_nodi_el = lel(i,:)
% dim: dimensione della griglia, ovvero numero totale di nodi
% ne: numero di elementi
% normali: le normali ai lati di bordo: trattandosi di un problema
%          monodimensionale tale campo sarà semplicemente il vettore
%          [-1 1], essendo i "lati" di bordo gli estremi del dominio
%          considerato
%
% Parametri di ingresso:
% xnodi: vettore RIGA con le coordinate dei nodi della griglia (non
%        necessariamente equispaziati)
% r: grado degli elementi utilizzati: i nodi vengono ripartiti tra i
%    singoli elementi in funzione del grado r scelto. Naturalmente il 
%    numero di nodi deve essere compatibile con il grado r.
%
% IMPORTANTE: occorre fare attenzione alla distinzione tra vertici e
% nodi: i vertici sono sempre ne+1 (essendo ne il numero di elementi),
% poiché ciascun elemento è delimitato da due vertici, mentre il
% numero di nodi dipende dal grado dello spazio di elementi finiti
% considerato. Considerando che ciascun elemento di grado r ha r+1
% nodi, e che i nodi estremi per ciascun elemento coincidono con i
% suoi vertici e sono collegati a due elementi adiacenti, otteniamo:
%   ne elementi di grado r  ->  ne*r+1 nodi
 
 griglia = struct('xnodi',[], ...
                  'lel',[], ...
                  'dim',[], ...
                  'normali',[]);
 
 griglia.xnodi = xnodi;
 for i=1:r:size(xnodi,2)-r
   % numerazione dei nodi: da sinistra a destra
   griglia.lel = [griglia.lel ; [i:1:i+r]];
 end
 griglia.normali = [-1 1];
 
 if size(griglia.xnodi,2)~=max(max(griglia.lel))
   error(['Il numero dei nodi non e'' compatibile ', ...
          'con il grado della base!']);
 end
 
 griglia.dim = size(xnodi,2);
 griglia.ne = size(griglia.lel,1);
 
return
