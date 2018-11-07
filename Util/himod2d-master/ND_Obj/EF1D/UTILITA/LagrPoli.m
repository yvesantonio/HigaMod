function [PHI,gradPHI,lapPHI] = LagrPoli(xnodi,x)
% [PHI,gradPHI,lapPHI] = LagrPoli(xnodi,x)
%
% Valutazione dei polinomi nodali di Lagrange e delle loro derivate
% prime e seconde (gradiente e laplaciano).
%
% xnodi: nodi di interpolazione, non necessariamente equispaziati
%        (vettore riga o colonna, indifferentemenete)
% x: ascisse in cui valutare i polinomi (vettore RIGA)
%
% PHI, gradPHI, lapPHI: polinomi e derivate
%
% N.B. si impiega la formula di Lagrange. Il grado e' scelto in modo
% da interpolare tutti i nodi: r = length(xnodi)-1

 % Il numero di funzioni di base e il grado sono dedotti dal numero di
 % nodi, utilizzando la condizione di interpolazione.
 N_r = length(xnodi);
 r = N_r-1;
 % Numero di punti in cui valutare i polinomi.
 M = length(x);
 
 if(r<=0)
   error('In LagrPoli servono almeno due nodi di Lagrange');
 end
 
 % inizializzazioni
 PHI     =  ones(N_r,M);
 gradPHI = zeros(N_r,M);
 lapPHI  = zeros(N_r,M);

 % PHI
 for i=1:N_r  % ciclo sulle funzioni di forma
   for j=1:N_r  % ciclo sui nodi
     if j~=i  % polinomi caratteristici di Lagrange
       PHI(i,:) = PHI(i,:) .* (x-xnodi(j))/(xnodi(i)-xnodi(j));
     end
   end
 end

 % gradPHI
 for i=1:N_r  % ciclo sulle funzioni di forma
   for s=1:N_r  % ciclo sui nodi
     if s~=i
       prod = ones(1,M);  % inizializzazione
       for j=1:N_r  % ciclo interno sui nodi
         if and((j~=i),(j~=s))  % polinomi caratteristici di Lagrange
           prod = prod .* (x-xnodi(j))/(xnodi(i)-xnodi(j));
         end
       end
       gradPHI(i,:) = gradPHI(i,:) + 1/(xnodi(i)-xnodi(s))*prod;
     end
   end
 end
 
 % lapPHI
 if r~=1  % altrimenti lasciare 0  
   for i=1:N_r  % ciclo sulle funzioni di forma
     for s=1:N_r  % primo ciclo sui nodi
       if (s~=i)  
         sommal = zeros(1,M);
         for l=1:N_r
           if and((l~=i),(l~=s))
             prod = ones(1,M);  % inizializzazione
             for j=1:N_r  % ciclo interno sui nodi
               if and((j~=i),and((j~=s),(j~=l)))
                 % polinomi caratteristici di Lagrange
                 prod = prod .* (x-xnodi(j))/(xnodi(i)-xnodi(j));
               end
             end
             sommal = sommal + 1/(xnodi(i)-xnodi(l))*prod;
           end
         end
         lapPHI(i,:) = lapPHI(i,:) + 1/(xnodi(i)-xnodi(s))*sommal;
       end
     end
   end
 end
 
return

