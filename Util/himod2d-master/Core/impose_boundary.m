function [Arid, brid]=impose_boundary(size_mb,omega_xl,omega_xr,hx,BC_l, boundary_l, BC_r, boundary_r,A,b)
%
% Impone tutte le condizioni di bordo relative ad interfacce, inflow e outflow.
% Se necessario riduce il problema eliminando i dof di Dirichlet.
%

% ne intervalli e nx estremi
ne = round((omega_xr-omega_xl)/hx);
nx = ne+1;
% Mesh FEM in x: nodi equispaziati
meshx = zeros(nx,1);
meshx(1) = omega_xl;
for i=2:nx
    meshx(i) = meshx(i-1)+hx;
end

% inizializzazione b ridotto: verranno aggiunti blocchi per ogni frequenza.
brid=[];

% Trattamento del bordo: ritaglia nella matrice globale quello che serve in base alle bc scelte
if (strcmp(BC_l,'dir') && strcmp(BC_r,'dir'))
    
    for imb = 1 : size_mb     % Per ogni frequenza
        
        buff = zeros( nx-2, 1); % Inizializzo il contributo al termine noto del nodo di Dirichlet lungo quanto il numero dei nodi interni
        
        for jmb = 1 : size_mb
            
            Arid( (imb-1)*(nx-2) + 1 : imb*(nx-2) , (jmb-1)*(nx-2) + 1 : jmb*(nx-2) ) = ...
                A( (imb-1)*nx + 2     : imb*nx - 1 , (jmb-1)*nx + 2     : jmb*nx - 1 );
            
            buff = buff + A( (imb-1)*nx + 2 : imb*nx-1, (jmb-1)*nx + 1 )*boundary_l( jmb ) ...
                + A( (imb-1)*nx + 2 : imb*nx-1, jmb*nx         )*boundary_r( jmb );
            
        end
        
        brid = [ brid; b( (imb-1)*nx + 2 : imb*nx - 1 ) - buff ];
        
    end
    
elseif ( strcmp(BC_l,'dir') && strcmp(BC_r,'neu') )
    for imb = 1 : size_mb
        
        buff = zeros( nx-1, 1);
        
        for jmb = 1 : size_mb
            Arid( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                A( (imb-1)*nx + 2     : imb*nx     , (jmb-1)*nx + 2     : jmb*nx     );
            
            buff = buff + A( (imb-1)*nx + 2 : imb*nx, (jmb-1)*nx + 1 )*boundary_l(jmb);
            
        end
        
        brid= [ brid; b( (imb-1)*nx + 2 : imb*nx ) - buff ];
        
        brid(end) = brid( end ) + boundary_r ( imb ); % Dato di Neumann
    end
    
elseif ( strcmp(BC_l,'neu') && strcmp(BC_r,'dir') )
    
    for imb = 1 : size_mb
        
        buff = zeros( nx-1, 1);
        
        for jmb = 1 : size_mb
            Arid( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                A( (imb-1)*nx + 1     : imb*nx-1   , (jmb-1)*nx + 1     : jmb*nx - 1 );
            
            buff = buff + A( (imb-1)*nx + 1 : imb*nx - 1 , jmb*nx )*boundary_r( jmb );
        end
        
        brid = [ brid ; b( (imb-1)*nx + 1 : imb*nx - 1 ) - buff ];
        
        brid(end-nx+2)= brid(end-nx+2) + boundary_l( imb );
        
    end
elseif ( strcmp(BC_l,'rob') && (strcmp(BC_r,'rob') || (strcmp(BC_r,'neu') ) ) ) %robrob e robneu condividono lo stesso codice.
    
    Arid=A;
    brid=b;
    for imb = 1 : size_mb
        brid((imb-1)*nx+1) = brid((imb-1)*nx+1) + boundary_l ( imb );
        brid(imb*nx) = brid(imb*nx) + boundary_r ( imb );
    end
elseif ( strcmp(BC_l,'neu') && (strcmp(BC_r,'neu') ) )
    Arid=A;
    brid=b;
    for imb = 1 : size_mb
        brid((imb-1)*nx+1) = brid((imb-1)*nx+1) + boundary_l ( imb );
        brid(imb*nx) = brid(imb*nx) + boundary_r ( imb );
    end
elseif ( strcmp(BC_l,'dir') && strcmp(BC_r,'rob') )
    
    for imb = 1 : size_mb
        
        buff = zeros( nx-1, 1);
        
        for jmb = 1 : size_mb
            Arid( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                A( (imb-1)*nx + 2     : imb*nx     , (jmb-1)*nx + 2     : jmb*nx     );
            buff = buff + A( (imb-1)*nx + 2 : imb*nx, (jmb-1)*nx + 1 )*boundary_l(jmb);
            
        end
        
        brid= [ brid; b( (imb-1)*nx + 2 : imb*nx ) - buff ];
        brid(end) = brid( end ) + boundary_r ( imb ); % Dato di Robin Robin
    end
end
