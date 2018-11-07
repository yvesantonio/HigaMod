function [FCAk,FCbk] = termine_stab_FC_loc(xnodi,PHI,gradPHI,lapPHI, ...
    pesiGauss,N_r,M,mifun,gradmifun,afun,divafun,sigmafun,fvolfun,d,rho)
% [FCAk,FCbk] = termine_stab_FC_loc(xnodi,PHI,gradPHI,lapPHI, ...
%   pesiGauss,N_r,M,mifun,gradmifun,afun,divafun,sigmafun,fvolfun,d,rho)
%
% Funzione per il calcolo dei termini locali (matrice e termine noto)
% nel caso di stabilizzazione fortemente consistente del trasporto.
%
% xnodi: coordinate dei nodi dell'elemento
% PHI, gradPHI, lapPHI: valutazioni delle funzioni di forma e delle 
%                       loro derivate prime e seconde nei nodi di
%                       quadratura dell'elemento
% pesiGauss: pesi di quadratura
% N_r: numero di funzioni di forma
% M: numero di nodi di quadratura
% ***fun: dati del problema
% d,rho: parametri per la stabilizzazione fortemente consistente
 
 % PRINCIPALI VARIABILI LOCALI
 %
 % xGauss: nodi di quadratura nell'intervallo fisico
 % dx_dcsi, d2x_dcsi: derivate della trasformazione di coordinate
 %                    isoparametrica valutate in xGauss
 % dcsi_dx, d2csi_dx: derivate della trasformazione inversa
 % gradPHIx, lapPHIx: derivate rispetto alla coordinata fisica x
 % LsPHI, LssPHI, LPHI, rLPHI: operatori del problema differenziale 
 %                             valutati nei nodi di quadratura
 
 % N.B: per chiarezza si impiegano numerosi cicli `for' che, al fine
 % di velocizzare il programma, potrebbero essere vettorializzati.
 
 % Per ulteriori dettagli, si vedano i commenti nel file
 % termine_diffusione_loc.m, che ha una struttura del tutto analoga.

 xGauss    = xnodi * PHI;
 dx_dcsi   = xnodi * gradPHI;
 d2x_dcsi2 = xnodi * lapPHI;
 dcsi_dx   = dx_dcsi.^(-1);
 d2csi_dx2 = -(dx_dcsi.^(-3)).*d2x_dcsi2;
 
 gradPHIx = zeros(size(gradPHI));  % derivata delle PHI rispetto ad x
 lapPHIx = zeros(size(lapPHI)); % derivata seconda delle PHI rispetto 
                                % ad x
 for i=1:N_r
    gradPHIx(i,:) = gradPHI(i,:).*dcsi_dx;
    lapPHIx(i,:) = lapPHI(i,:).*(dcsi_dx.^2) + gradPHI(i,:).*d2csi_dx2;
 end
 
 hk = max(xnodi)-min(xnodi); % passo di griglia (ampiezza dell'elemento)
 
 mi = feval(mifun,xGauss);  % valutazione del coeff. mi
 gradmi = feval(gradmifun,xGauss);
 
 a = feval(afun,xGauss);  % valutazione del coeff. a
 diva = feval(divafun,xGauss);
 ak = mean(abs(a));
 
 sigma = feval(sigmafun,xGauss);  % valutazione del coeff. sigma
 
 f = feval(fvolfun,xGauss);  % valutazione di f
 
 % valutazione di LsPHI, LssPHI e LPHI
 LsPHI = zeros(size(PHI));
 LssPHI = LsPHI;
 LPHI = LsPHI;
 rLPHI = LsPHI;
 for i=1:N_r
    LsPHI(i,:) = -gradmi.*gradPHIx(i,:) - mi.*lapPHIx(i,:) + ...
                 (0.5*diva + sigma).*PHI(i,:);
    LssPHI(i,:) = a.*gradPHIx(i,:) + 0.5*diva.*PHI(i,:);
    LPHI(i,:) = LsPHI(i,:)+LssPHI(i,:);
    rLPHI(i,:) = LssPHI(i,:)+rho*LsPHI(i,:);
 end
 
 FCAk = zeros(N_r,N_r);
 FCbk = zeros(N_r,1);
 for l=1:M   % ciclo sui nodi di quadratura
   for i=1:N_r  % ciclo sulle funzioni test
     for j=1:N_r  % ciclo sulle funzioni shape
        
        FCAk(i,j) = FCAk(i,j) + d*hk/ak*pesiGauss(l) * ...
                                LPHI(j,l) * rLPHI(i,l) * ...
                                dx_dcsi(l);

     end
     
     FCbk(i) = FCbk(i) + d*hk/ak*pesiGauss(l) * ...
                         f(l) * rLPHI(i,l) * ...
                         dx_dcsi(l);
     
   end
 end
 
return

