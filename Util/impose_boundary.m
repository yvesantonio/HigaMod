function [Arid, brid]=impose_boundary(size_mb,numbNodes,BC_l, boundary_l, BC_r, boundary_r,A,b)
%
% Impone tutte le condizioni di bordo relative ad interfacce, inflow e outflow.
% Se necessario riduce il problema eliminando i dof di Dirichlet.
%

nx = numbNodes;

% inizializzazione b ridotto: verranno aggiunti blocchi per ogni frequenza.
brid=[];

% Trattamento del bordo: ritaglia nella matrice globale quello che serve in base alle bc scelte
if (strcmp(BC_l,'dir') && strcmp(BC_r,'dir'))
    
    for imb = 1 : size_mb     % Per ogni frequenza
        
        % buff = zeros( nx-2, 1);
        % 
        % for jmb = 1 : size_mb
        %     
        %     Arid( (imb-1)*(nx-2) + 1 : imb*(nx-2) , (jmb-1)*(nx-2) + 1 : jmb*(nx-2) ) = ...
        %         A( (imb-1)*nx + 2     : imb*nx - 1 , (jmb-1)*nx + 2     : jmb*nx - 1 );
        %     
        %     buff = buff + A( (imb-1)*nx + 2 : imb*nx-1, (jmb-1)*nx + 1 )*boundary_l( jmb ) ...
        %         + A( (imb-1)*nx + 2 : imb*nx-1, jmb*nx         )*boundary_r( jmb );
        %     
        % end
        % 
        % brid = [ brid; b( (imb-1)*nx + 2 : imb*nx - 1 ) - buff ];
        
        % New
        % .................................................................
        
        buff = zeros( nx-2, 1);
        
        for jmb = 1 : size_mb
            
            Arid( (imb-1)*(nx-2) + 1 : imb*(nx-2) , (jmb-1)*(nx-2) + 1 : jmb*(nx-2) ) = ...
                A( (imb-1)*nx + 2     : imb*nx - 1 , (jmb-1)*nx + 2     : jmb*nx - 1 );
            
            Aaux((imb-1) * nx + 1 : imb * nx , (jmb-1) * nx + 1 : jmb * nx) = ...
                [1e10                                                                        , ...
                 zeros(1,nx-2)                                                            , ...
                 0                                                                        ; ...
                 zeros(nx-2,1)                                                            , ...
                 Arid( (imb-1)*(nx-2) + 1 : imb*(nx-2) , (jmb-1)*(nx-2) + 1 : jmb*(nx-2) ), ...
                 zeros(nx-2,1)                                                            ; ...
                 0                                                                        , ...
                 zeros(1,nx-2)                                                            , ...
                 1e10];
            
            buff = buff + A( (imb-1)*nx + 2 : imb*nx-1, (jmb-1)*nx + 1 )*boundary_l( jmb ) ...
                + A( (imb-1)*nx + 2 : imb*nx-1, jmb*nx         )*boundary_r( jmb );
            
        end
        
        Baux = [boundary_l( imb ) * 1e10; b( (imb-1)*nx + 2 : imb*nx - 1 ) - buff; boundary_r( imb ) * 1e10];
        brid = [ brid; Baux ];
        
        % .................................................................
        
    end
    
    Arid = Aaux;
    
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
