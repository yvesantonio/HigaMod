function Kk = termine_diffusione_loc(xnodi,PHI,gradPHI,pesiGauss, ...
                                     N_r,M,mufun)
% Kk = termine_diffusione_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M,mufun)
%
% Funzione per il calcolo della matrice di rigidezza locale.
%
% xnodi: coordinate dei nodi dell'elemento
% PHI, gradPHI: valutazioni delle funzioni di forma nei nodi di
%               quadratura dell'elemento
% pesiGauss: pesi di quadratura
% N_r: numero di funzioni di forma
% M: numero di nodi di quadratura
% mufun: funzione per la valutazione del coefficiente di diffusione
%
% Kk: matrice di rigidezza locale
 
 % PRINCIPALI VARIABILI LOCALI
 %
 % xGauss: nodi di quadratura sull'elemento reale
 % dx_dcsi: derivata della trasformazione di coordinate isoparametrica
 %          valutata in xGauss
 % mu: coefficiente di diffusione valutato in xGauss
 
 % N.B: per chiarezza si impiegano numerosi cicli `for' che, al fine
 % di velocizzare il programma, potrebbero essere vettorializzati.
 
 % xGauss viene valutato utilizzando la trasformazione di coordinate
 % isoparametrica. Per un generco punto csi nell'elemento di
 % riferimento questa è
 %
 %   x(csi) = sum_i xn_i * phi_i(csi)
 % 
 % essendo xn_i la coordinata del nodo i-esimo e phi_i la i-esima
 % funzione di forma (si ricordi che, per definizione, in un elemento
 % finito isoparametrico la trasformazione di coordinate el. di
 % riferimento -> el. reale utilizza le stesse funzioni di forma che
 % approssimano l'incognita).
 %
 % Nel caso particolare dei nodi di quadratura, essendo PHI e gradPHI
 % le valutazioni delle funzioni di forma nei nodi di quadratura
 % dell'elelemto di riferimento, si ottiene:
 xGauss  = xnodi * PHI;
 dx_dcsi = xnodi * gradPHI;

 % A questo punto possiamo valutare il coefficiente di diffusione:
 mu = feval(mufun,xGauss);
 
 % Calcolo con quadratura gaussiana dei termini
 %
 % K_ij = integrale sull'elemento di (mu * dphi_i/dx * dphi_j/dx * dx)
 %
 % N.B. L'integrale viene riportato sull'elemento campione:
 %
 % K_ij = integrale di (mu * dphi_i/dcsi*dcsi/dx * ...
 %                      dphi_j/dcsi*dcsi/dx * dx/dcsi*dcsi) = 
 %        integrale di (mu * dcsi/dx * dphi_i/dcsi*dphi_j/dcsi * dcsi)
 %
 % N.B.: la quadratura gaussiana di una generica funzione f(x)
 % consiste in
 %
 % integrale su [-1 1] di (f(x) * dx) ~= ...
 %   somma sugli l nodi di (pesoGauss_l * f(nodo_l))
 
 Kk = zeros(N_r,N_r);
 for l=1:M   % ciclo sui nodi di quadratura
   for i=1:N_r  % ciclo sulle funzioni test
     for j=1:N_r  % ciclo sulle funzioni shape
        
        Kk(i,j) = Kk(i,j) + pesiGauss(l)*mu(l)/dx_dcsi(l) * ...
	                    gradPHI(i,l)*gradPHI(j,l);
		  
     end
   end
 end
 
return

