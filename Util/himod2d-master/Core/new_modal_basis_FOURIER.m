function [mb,mb_y]=new_modal_basis_FOURIER(m, nodi)
%
% function [mb,mb_y]=new_modal_basis_FOURIER(m, nodi)
%
% m         = numero delle funzioni di base (Dirichlet), m+1 nel caso Robin
% nodi        = punti in cui si intende valutare le funzioni di base, coincidono spesso con i nodi di quadratura
% bc_up     = natura delle boundary conditions in 1 'dir' oppure 'rob'
% bc_down   = natura delle boundary conditions in 0 'dir' oppure 'rob'
%
%
% Restituisce:
% mb   = matrice che contiene il valore delle funzioni di base (colonne) valutate nei punti (righe)
% mb_y = matrice che contiene il valore delle derivate delle funzioni di base (colonne) valutate nei punti (righe)
%



     mb   = zeros( length(nodi), m);
     mb_y = zeros( length(nodi), m);

     for i=1:length(nodi)
         mb  ( i , 1 ) = 1;
         mb_y( i , 1 ) = 0;
     end
     % Ciclo sulle funzioni base
     for n = 2:2:m
         % Ciclo nella sezione
         for i = 1:length(nodi)
             mb  ( i, n ) = sqrt(2)*sin(n*pi*nodi(i));
             mb_y( i, n ) = n*pi*sqrt(2)*cos(n*pi*nodi(i));
         end

     end
     for n= 3:2:m
         % Ciclo nella sezione
         for i = 1:length(nodi)
             mb  ( i, n ) = sqrt(2)*cos((n-1)*pi*nodi(i));
             mb_y( i, n ) = -(n-1)*pi*sqrt(2)*sin((n-1)*pi*nodi(i));
         end
     end
