function [u,dir_l,dir_r,neu_l,neu_r,rob_l,rob_r]=extract_solution(size_mb,omega_xl,omega_xr,hx,BC_l, boundary_l, BC_r, boundary_r, ur, A, b,Coeff_forma,gamma)
%
%	Aggiunge i valori di Dirichlet (INFLOW) alla soluzione e calcola le quantita' di interfaccia per il ciclo DD
%
dir_l=0;
dir_r=0;
neu_l=0;
neu_r=0;
rob_l=0;
rob_r=0;

% ne numero di intervalli e nx estremi (nx = ne+1)
ne = round((omega_xr-omega_xl)/hx);
nx = ne+1; % Numero di nodi in x

% Mesh FEM in x: nodi equispaziati
meshx    = zeros( nx, 1);
meshx(1) = omega_xl;
for i = 2 : nx
    meshx(i) = meshx(i-1) + hx;
end

if ( strcmp(BC_l, 'dir') && strcmp(BC_r, 'dir'))
    
    u=[];
    
    for imb = 1 : size_mb
        u = [u; boundary_l( imb ); ur( (imb-1)*(nx-2) + 1 : imb*(nx-2) ); boundary_r( imb )];
    end
    
elseif ( strcmp(BC_l, 'dir') && strcmp(BC_r, 'neu') )
    
    u=[];
    
    for imb = 1 : size_mb
        u = [u; boundary_l( imb ); ur( (imb-1)*(nx-1) + 1 : imb*(nx-1) ) ];
    end
    
elseif ( strcmp(BC_l, 'neu') && strcmp(BC_r, 'dir'))
    
    u=[];
    
    for imb = 1 : size_mb
        u = [u; ur( (imb-1)*(nx-1) + 1 : imb*(nx-1) ); boundary_r( imb )];
    end
    
elseif ( strcmp(BC_l, 'dir') && strcmp(BC_r, 'rob') )
    
    u=[];
    
    for imb = 1 : size_mb
        u = [u; boundary_l( imb ); ur( (imb-1)*(nx-1) + 1 : imb*(nx-1) ) ];
    end
    
elseif (strcmp(BC_l, 'rob') && strcmp(BC_r, 'rob') )
    u=ur;
elseif (strcmp(BC_l, 'rob') && strcmp(BC_r, 'neu') )
    u=ur;
elseif (strcmp(BC_l, 'neu') && strcmp(BC_r, 'neu') )
    u=ur;
else 
	display('Error in extract_solution.m');
end

% Vecchio codice per l'algoritmo Dirichlet Neumann, non viene mai chiamato
% Viene lasciato se si volesse provare DN
if(~   ( strcmp(BC_l,'rob') || strcmp(BC_r,'rob') )    )
    
    % Estrazione della soluzione al bordo (Dirichlet e Neumann)
    dir_l = zeros( size_mb, 1 );
    dir_r = zeros( size_mb, 1 );
    neu_l = zeros( size_mb, 1 );
    neu_r = zeros( size_mb, 1 );
    
    for imb = 1 : size_mb
        dir_l( imb ) = u ( (imb-1)*nx + 1  );
        dir_r( imb ) = u ( (imb-1)*nx + nx );
    end
    
    residuo=b-A*u; %residuo globale
    for imb = 1 : size_mb
        neu_l( imb ) = - residuo( (imb-1)*nx + 1  );
        neu_r( imb ) = - residuo( (imb-1)*nx + nx );
    end
else %Caso Robin Robin!
    mu_valutato=Coeff_forma.mu(1,1);%attenzione va cambiato nel caso mu non costante, non ha nessun senso fisico così.
    for imb= 1 : size_mb
        rob_l(imb) = mu_valutato*(u( (imb-1)*nx + 2 )-u( (imb-1)*nx + 1))/hx+gamma.R*u( (imb-1)*nx + 1 );
        rob_r(imb) = -mu_valutato*(u(imb*nx)-u(imb*nx-1))/hx+gamma.L*u(imb*nx);
    end
end
end
