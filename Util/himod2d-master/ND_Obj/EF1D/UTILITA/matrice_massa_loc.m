function Mk = matrice_massa_loc(xnodi,PHI,gradPHI,pesiGauss, ...
                                N_r,M,rhofun)
% Mk = matrice_massa_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M,rhofun)
% 
% Funzione per il calcolo della matrice di massa locale.
%
% xnodi: coordinate dei nodi dell'elemento
% PHI, gradPHI: valutazioni delle funzioni di forma nei nodi di
%               quadratura dell'elemento
% pesiGauss: pesi di quadratura
% N_r: numero di funzioni di forma
% M: numero di nodi di quadratura
% rhofun: funzione per la valutazione della densità
%
% Mk: matrice di massa locale

 % PRINCIPALI VARIABILI LOCALI
 %
 % xGauss: nodi di quadratura sull'elemento reale
 % dx_dcsi: derivata della trasformazione di coordinate isoparametrica
 %          valutata in xGauss
 % rho: densità valutata in xGauss
 
 % N.B: per chiarezza si impiegano numerosi cicli `for' che, al fine
 % di velocizzare il programma, potrebbero essere vettorializzati.
 
 % Per ulteriori dettagli, si vedano i commenti nel file
 % termine_diffusione_loc.m, che ha una struttura del tutto analoga.

 xGauss = xnodi * PHI;
 dx_dcsi = xnodi * gradPHI;
 rho = feval(rhofun,xGauss);
 
 Mk = zeros(N_r,N_r);
 for l=1:M   % ciclo sui nodi di quadratura
   for i=1:N_r  % ciclo sulle funzioni test
     for j=1:N_r  % ciclo sulle funzioni shape
      
        Mk(i,j) = Mk(i,j) + pesiGauss(l)*rho(l) * ...
	                    PHI(i,l)*PHI(j,l) * dx_dcsi(l);
		  
     end
   end
 end
 
return

