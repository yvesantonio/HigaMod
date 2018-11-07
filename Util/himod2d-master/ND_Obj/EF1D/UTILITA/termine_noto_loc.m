function bk = termine_noto_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M, ...
                               fvolfun)
% bk = termine_noto_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M,fvolfun)
% 
% Funzione per il calcolo del termine noto locale.
%
% xnodi: coordinate dei nodi dell'elemento
% PHI, gradPHI: valutazioni delle funzioni di forma nei nodi di
%               quadratura dell'elemento
% pesiGauss: pesi di quadratura
% N_r: numero di funzioni di forma
% M: numero di nodi di quadratura
% fvolfun: funzione per la valutazione del termine noto f

 % PRINCIPALI VARIABILI LOCALI
 %
 % xGauss: nodi di quadratura nell'intervallo fisico
 % dx_dcsi: derivata della trasformazione di coordinate isoparametrica
 %          valutata in xGauss
 % f: termine noto valutato in xGauss
 
 % N.B: per chiarezza si impiegano numerosi cicli `for' che, al fine
 % di velocizzare il programma, potrebbero essere vettorializzati.

 % Per ulteriori dettagli, si vedano i commenti nel file
 % termine_diffusione_loc.m, che ha una struttura del tutto analoga.
 
 xGauss = xnodi * PHI;
 dx_dcsi = xnodi * gradPHI;
 f = feval(fvolfun,xGauss);
 
 bk = zeros(N_r);
 for l=1:M   % ciclo sui nodi di quadratura
   for i=1:N_r  % ciclo sulle funzioni test
     
      bk(i) = bk(i) + pesiGauss(l)*f(l)*PHI(i,l)*dx_dcsi(l);
		  
   end
 end
 
return

