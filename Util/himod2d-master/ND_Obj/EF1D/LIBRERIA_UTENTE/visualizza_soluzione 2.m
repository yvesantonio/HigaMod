function [x,y] = visualizza_soluzione(griglia,base,u,N,varargin)
% Funzione per la visualizzazione della soluzione numerica.
%
% La visualizzazione viene costruita considerando il valore della
% soluzione u_h nei nodi della griglia, nei punti di Gauss, e in
% ulteriori punti equispaziati specificati dal parametro N. Le ascisse
% e le ordinate usate per tracciare il grafico sono restituite come
% parametri x,y che possono poi essere usati per valutare altre
% quantità nei medesimi punti.
%
% [x,y] = visualizza_soluzione(griglia,base,u,N)
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta con
%          'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% u: vettore dei valori nodali della soluzione u_h
% N: numero di punti (in aggiunta a nodi e punti di Gauss) da
%    impiegare per la visualizzazione su ciascun elemento: la
%    rappresentazione sarà tanto più "liscia" quanto più N è elevato
%
% x,y: valori usati per il grafico: plot(x,y)
%
% [x,y] = visualizza_soluzione(griglia,base,u,N,fig)
%
% Come nel caso precedente, ma viene usata la figura fig senza crearne
% una nuova

 N_r = base.N_r;
 M   = base.M;
 PHI = base.PHI;

 % Gli N punti richiesti per la visualizzazione vengono generati
 % sull'elemento di riferimento [-1,1].
 csi_grafica = linspace(-1,1,N+2);
 % Gli estremi dell'intervallo -1 e 1 compaiono già tra i nodi
 % dell'elemento, e dunque possono essere eliminati.
 csi_grafica(end) = [];
 csi_grafica(1) = [];
 
 % Valutazione delle funzioni di forma nei punti csi_grafica
 PHI_grafica = LagrPoli(base.csinodi,csi_grafica); 
 
 % Ciascun elemento viene ora visualizzato utilizzando i punti appena
 % definiti, ai quali vengono aggiunti anche i nodi e i punti di
 % Gauss. Tutti i valori necessari per tracciare il grafico vengono
 % raccolti nei vettori x, y, XG, YG.
 x = [];
 y = [];
 XG = [];
 YG = [];
 for ie=1:griglia.ne

    % trasformazione di coordinate per elementi isoparametrici: si
    % veda UTILITA/termine_diffusione_loc.m per ulteriori dettagli
    xnodi = griglia.xnodi(griglia.lel(ie,:))';
    xGauss    =         PHI' * xnodi;
    x_grafica = PHI_grafica' * xnodi;
    
    % valutazione della soluzione u_h nei punti usati per tracciare il
    % grafico: vengono prima ottenuti i valori nodali ynodi, poi tali
    % valori sono utilizzati come coefficienti delle funzioni di forma
    ynodi = u(griglia.lel(ie,:)); % valori nodali per l'elemento ie
    yGauss    =         PHI' * ynodi; % u_h nei nodi di Gauss
    y_grafica = PHI_grafica' * ynodi; % u_h nei punti x_grafica
    
    % vettori globali per i soli punti di Gauss
    XG = [XG;xGauss];
    YG = [YG;yGauss];
    
    % vettori locali che includono tutti i punti usati per tracciare
    % il grafico: nodi, punti di Gauss e punti aggiuntivi
    x_loc = [xnodi;xGauss;x_grafica];
    y_loc = [ynodi;yGauss;y_grafica];
    
    % i vettori devono essere riordinati affinché la coordinata x sia
    % sempre crescente
    [x_loc,indici] = sort(x_loc);
    y_loc = y_loc(indici);
    
    % escludo l'ultimo elemento per evitare di ripetere i nodi
    % estremi, comuni a due elementi
    x_loc(end) = [];
    y_loc(end) = [];
    
    % vettori globali che includono tutti i punti
    x = [x;x_loc];
    y = [y;y_loc];
    
 end
 % Bisogna adesso reintrodurre l'ultimo nodo dell'ultimo elemento.
 x = [x;griglia.xnodi(end)];
 y = [y;      u(end)      ];
 
 % A questo punto ci sono tutti i dati per tracciare il grafico.
 indici_vertici = [1:N_r-1:griglia.dim];
 indici_nodi_interni = setdiff([1:griglia.dim],indici_vertici);
 if(nargin==4) % creare una nuova figura
   fig = figure;
 else % usare una figura preesistente
   fig = figure(varargin{1});
   hold off
 end
 hh = plot(griglia.xnodi(indici_vertici),u(indici_vertici),'sg');
 set(hh,'MarkerSize',10);
 hold on
 r = N_r-1; % grado degli elementi
 if(r>1)
   plot(griglia.xnodi(indici_nodi_interni), ...
                    u(indici_nodi_interni),'or');
   plot(XG,YG,'+k');
   legend('vertici','nodi','punti di Gauss');
 else
   plot(XG,YG,'+k');
   legend('vertici==nodi','punti di Gauss');
 end
 plot(x,y);
 
return

