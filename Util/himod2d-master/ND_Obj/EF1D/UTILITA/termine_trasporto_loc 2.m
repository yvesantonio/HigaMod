function Ck = termine_trasporto_loc(xnodi,PHI,gradPHI,pesiGauss, ...
                                    N_r,M,afun)
% Ck = termine_trasporto_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M,afun)
% 
% Funzione per il calcolo della matrice di trasporto locale.
%
% xnodi: coordinate dei nodi dell'elemento
% PHI, gradPHI: valutazioni delle funzioni di forma nei nodi di
%               quadratura dell'elemento
% pesiGauss: pesi di quadratura
% N_r: numero di funzioni di forma
% M: numero di nodi di quadratura
% afun: funzione per la valutazione del coefficiente di trasporto
%
% Ck: matrice di trasporto locale
 
 % PRINCIPALI VARIABILI LOCALI
 %
 % xGauss: nodi di quadratura sull'elemento reale
 % dx_dcsi: derivata della trasformazione di coordinate isoparametrica
 %          valutata in xGauss
 % a: coefficiente di trasporto valutato in xGauss
 
 % N.B: per chiarezza si impiegano numerosi cicli `for' che, al fine
 % di velocizzare il programma, potrebbero essere vettorializzati.
 
 % Per ulteriori dettagli, si vedano i commenti nel file
 % termine_diffusione_loc.m, che ha una struttura del tutto analoga.

 xGauss = xnodi * PHI;
 dx_dcsi = xnodi * gradPHI;
 a = feval(afun,xGauss);
 
 Ck = zeros(N_r,N_r);
 for l=1:M   % ciclo sui nodi di quadratura
   for i=1:N_r  % ciclo sulle funzioni test
     for j=1:N_r  % ciclo sulle funzioni shape
      
        Ck(i,j) = Ck(i,j) - pesiGauss(l)*a(l) * ...
	                    gradPHI(i,l)*PHI(j,l);
		  
     end
   end
 end
 
return

