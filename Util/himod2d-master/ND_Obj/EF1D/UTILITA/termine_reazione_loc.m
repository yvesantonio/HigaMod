function Sk = termine_reazione_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M, ...
                                   sigmafun)
% Sk = termine_reazione_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M,sigmafun)
% 
% Funzione per il calcolo della matrice del termine di reazione locale
% (matrice di massa).
%
% xnodi: coordinate dei nodi dell'elemento
% PHI, gradPHI: valutazioni delle funzioni di forma nei nodi di
%               quadratura dell'elemento
% pesiGauss: pesi di quadratura
% N_r: numero di funzioni di forma
% M: numero di nodi di quadratura
% sigmafun: funzione per la valutazione del coefficiente di reazione
%
% Sk: matrice del termine di reazione locale
 
 % PRINCIPALI VARIABILI LOCALI
 %
 % xGauss: nodi di quadratura sull'elemento reale
 % dx_dcsi: derivata della trasformazione di coordinate isoparametrica
 %          valutata in xGauss
 % sigma: coefficiente di reazione
 
 % N.B: per chiarezza si impiegano numerosi cicli `for' che, al fine
 % di velocizzare il programma, potrebbero essere vettorializzati.

 % Per ulteriori dettagli, si vedano i commenti nel file
 % termine_diffusione_loc.m, che ha una struttura del tutto analoga.
 
 xGauss = xnodi * PHI;
 dx_dcsi = xnodi * gradPHI;
 sigma = feval(sigmafun,xGauss);
 
 Sk = zeros(N_r,N_r);
 for l=1:M   % ciclo sui nodi di quadratura
   for i=1:N_r  % ciclo sulle funzioni test
     for j=1:N_r  % ciclo sulle funzioni shape
      
        Sk(i,j) = Sk(i,j) + pesiGauss(l)*sigma(l) * ...
	                    PHI(i,l)*PHI(j,l) * dx_dcsi(l);
		  
     end
   end
 end
 
return

