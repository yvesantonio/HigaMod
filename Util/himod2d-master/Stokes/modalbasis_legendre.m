function [mb,mb_y]=modalbasis_legendre(m, nodi, bc_up, bc_down)
%
% function [mb,mb_y]=new_modal_basis_FOURIER(m, nodi)
%
% m         = numero delle funzioni di base (Dirichlet)
% nodi        = punti in cui si intende valutare le funzioni di base, coincidono spesso con i nodi di quadratura
%
% Restituisce:
% mb   = matrice che contiene il valore delle funzioni di base (colonne) valutate nei punti (righe)
% mb_y = matrice che contiene il valore delle derivate delle funzioni di base (colonne) valutate nei punti (righe)

if(strcmp(bc_up,'rob') && strcmp(bc_down,'rob'))
    mb   = zeros( length(nodi), m);
    mb_y = zeros( length(nodi), m);
    
    [P,Pd]=polynomial_legendre(m);
    % Ciclo sulle funzioni base
    for n = 1:m
        % Ciclo nella sezione
        mb(:,n)=polyval(P{n},nodi);
        mb_y(:,n)=polyval(Pd{n},nodi);
        %          for i = 1:length(nodi)
        %              mb  ( i, n ) = nodi(i).^(n-1);
        %              mb_y( i, n ) = (n-1)*nodi(i).^(n-2);
        %          end
    end
    
elseif(strcmp(bc_up,'dir') && strcmp(bc_down,'dir'))
    
    mb   = zeros( length(nodi), m);
    mb_y = zeros( length(nodi), m);
    
    % Ciclo sulle funzioni base
    for n = 1:m
        % Ciclo nella sezione
        for i = 1:length(nodi)
            mb  ( i, n ) = nodi(i).^(n-1)*(1-nodi(i).^2);
            mb_y( i, n ) = (n-1)*nodi(i).^(n-2)*(1-nodi(i).^2)+nodi(i).^(n-1)*(-2*nodi(i));
        end
    end
    
end
end