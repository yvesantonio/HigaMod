function visualizza_evoluzione(griglia,base,u,N,varargin)
% Animazione della soluzione numerica per problemi evolutivi.
%
% visualizza_evoluzione(griglia,base,u,N)
%
% griglia, base, N: come in visualizza_soluzione
% u: matrice le cui colonne rappresentano la soluzione
%    ai diversi istanti temporali
%
% visualizza_evoluzione(griglia,base,u,N,tempi)
%
% Rispetto alla forma precendente, inserisce un titolo con il livello
% temporale di ciascuna immagine, utilizzando i valori raccolti nel
% vettore tempi.
% 
% griglia, base, u, N: come al caso precedente
% tempi: vettore con i tempi in cui è calcolata u
%
% N.B. per terminare premere "CTRL + C" !

 h = figure;
 xnodi = griglia.xnodi;
 limiti = [min(xnodi) max(xnodi) min(min(u)) max(max(u))];
 if nargin==5
   tempi = varargin{1};
 end
 disp('CTRL + C per terminare')
 for j=0:inf  % ripeti indefinitamente
   for i=1:size(u,2)  % ciclo sui livelli temporali
     visualizza_soluzione(griglia,base,u(:,i),N,h);
     axis(limiti);
     if nargin==5
       title(sprintf('Tempo: %f',tempi(i)));
     end
     pause(1)
   end
 end
 
return
     
