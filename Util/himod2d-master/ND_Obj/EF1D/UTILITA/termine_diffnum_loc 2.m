function Kk = termine_diffnum_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M, ...
                                  afun,mifun,stabilizzazione)
% Kk = termine_diffnum_loc(xnodi,PHI,gradPHI,pesiGauss,N_r,M, ...
%                          afun,mifun,stabilizzazione)
%
% Funzione per il calcolo del terminie di diffusione numerico per la
% stabilizzazione dei problemi di diffusione-trasporto. Si ricalca la
% struttura di termine_diffusione_loc, con la differenza che la
% viscosità è quella numerica, funzione del numero di Péclet locale.
%
% xnodi: coordinate dei nodi dell'elemento
% PHI, gradPHI: valutazioni delle funzioni di forma nei nodi di
%               quadratura dell'elemento
% pesiGauss: pesi di quadratura
% N_r: numero di funzioni di forma
% M: numero di nodi di quadratura
% afun, mifun: funzioni per la valutazione dei coefficienti di
%              trasporto e diffusione, rispettivamente
% stabilizzazione: stringa con valore 'upwind' o 'SG' (Scharfetter 
%                  Gummel)
 
 xGauss = xnodi * PHI;
 dx_dcsi = xnodi * gradPHI;

 a = feval(afun,xGauss);     % valutazione del coefficiente a
 mi = feval(mifun,xGauss);   % valutazione del coefficiente mi
 h = max(xnodi)-min(xnodi);  % passo di griglia (ampiezza dell'elemento)
 Pe = (mean(abs(a))*h) / (2*mean(mi)); % Péclet locale medio
 
 % Viene ora calcolata, in funzione del numero di Péclet e del metodo
 % di stabilizzazione scelto, la viscosità numerica
 switch stabilizzazione
  case{'upwind'}
   mi_num = Pe*mi;
  case{'SG'}
   if Pe>0
     B = 2*Pe/(exp(2*Pe)-1);
   else
     B = 1;
   end
   mi_num = (Pe-1+B)*mi;
  otherwise
   error(sprintf('stabilizzazione "%s" sconosciuta',stabilizzazione));
 end
 
 Kk = zeros(N_r,N_r);
 for l=1:M   % ciclo sui nodi di quadratura
   for i=1:N_r  % ciclo sulle funzioni test
     for j=1:N_r  % ciclo sulle funzioni shape
      
        Kk(i,j) = Kk(i,j) + pesiGauss(l)*mi_num(l)/dx_dcsi(l) * ...
	                    gradPHI(i,l)*gradPHI(j,l);
		  
     end
   end
 end
 
return

