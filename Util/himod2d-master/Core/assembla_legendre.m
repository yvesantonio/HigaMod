function [Al,bl]=assembla_legendre(imb,kmb,nx, mesh_wx,wyq,mb_i,mb_yi,mb_k,mb_yk, FEM_basis, FEM_basis_x,bc_up, bc_down, Coeff_forma, Computed)

nqnx=length(mesh_wx)/(nx-1);

coeffrobin=Coeff_forma.coeffrobin;

forc  = ( Computed.jacforce)*( mb_i.*wyq);

%Calcolo dei coefficienti \hat r_{ik}^{st} k here refers to the solution
%muv       = integrate( Computed.mu_c          , mb_k .*mb_i .*wyq );
%muv_y     = integrate( Computed.mu_c          , mb_yk.*mb_yi.*wyq );
%betav     = integrate( Computed.beta1_c       , mb_k .*mb_i .*wyq );
%betav_y   = integrate( Computed.beta2_c       , mb_yk.*mb_i .*wyq );
%sigmav    = integrate( Computed.sigma_c       , mb_k .*mb_i .*wyq );


r11v = ( Computed.jacmu)*( mb_k .*mb_i .*wyq);
r01v = ( Computed.jacmuD)*(mb_yk .*mb_i .*wyq );
r10v = ( Computed.jacmuD)*(mb_k .*mb_yi .*wyq );
r10v = r10v + ( Computed.jacbeta1 )*(mb_k .*mb_i .*wyq );
r00v = ( Computed.jacmud2j2)*(mb_yk.*mb_yi.*wyq );
r00v = r00v + ( Computed.jacsigma      )*( mb_k .*mb_i .*wyq );
r00v = r00v + ( Computed.jacbeta1Dejacbeta2j)*(mb_yk .*mb_i .*wyq);

%determina le basi valutate solo nei xpunti della frontiera UP e DOWN
base_mod_aus = new_modal_basis_legendre( max(imb,kmb), [-1,1],bc_up,bc_down);

if(strcmp(bc_up,'rob'))
    r00v  = r00v + coeffrobin*base_mod_aus(2, imb)*base_mod_aus(2, kmb)*Computed.a_up;
end

if(strcmp(bc_down,'rob'))
    r00v  = r00v + coeffrobin*base_mod_aus(1, imb)*base_mod_aus(1, kmb)*Computed.a_down;
end

%Inizializzazione del blocco di matrice
%Al = sparse(nx,nx);
bl = zeros(nx,1);

%determina le basi valutate solo nei punti della frontiera UP e DOWN
% base_mod_aus = new_modal_basis_legendre( max(imb,kmb), [-1,1]);
%
% if(strcmp(bc_up,'rob'))
%     funzionale_robin_up   = coeffrobin*base_mod_aus(2, imb)*base_mod_aus(2, kmb);
% else
%     funzionale_robin_up   = 0;
% end
%
% if(strcmp(bc_down,'rob'))
%     funzionale_robin_down = coeffrobin*base_mod_aus(1, imb)*base_mod_aus(1, kmb);
% else
%     funzionale_robin_down = 0;
% end

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag = zeros( nx, 1);
diagsup  = zeros( nx, 1); % Spdiags taglierà il primo elemento
diaginf  = zeros( nx, 1); % Spdiags invece taglierà l'ultimo

for ie = 1 : nx-1  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    %m  = muv   ( (ie-1)*nqnx+1 : ie*nqnx)';
    %my = muv_y ( (ie-1)*nqnx+1 : ie*nqnx)';
    %b  = betav ( (ie-1)*nqnx+1 : ie*nqnx)';
    %by = betav_y ( (ie-1)*nqnx+1 : ie*nqnx)';
    %s  = sigmav( (ie-1)*nqnx+1 : ie*nqnx)';
    r11 = r11v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r01 = r01v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:2                  % Per ogni funzione di forma EF
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basis  ( (ie-1)*nqnx+1:ie*nqnx, i)';
        femb_xi = FEM_basis_x( (ie-1)*nqnx+1:ie*nqnx, i)';
        for j = 1:2              % Per ogni funzione di forma EF
            
            femb_j  = FEM_basis  ( (ie-1)*nqnx+1:ie*nqnx, j)';
            femb_xj = FEM_basis_x( (ie-1)*nqnx+1:ie*nqnx, j)';
            
            if(i==j)
                %                 maindiag(ie+i-1)=maindiag(ie+i-1)+integrate(...
                %                       r11.*femb_xi.*femb_xj ...
                %                     + r01.*femb_xi.*femb_j...
                %                     + r10.*femb_i.*femb_xj ...
                %                     + r00.*femb_i.*femb_j  ...
                %                     , w);
                maindiag(ie+i-1)=maindiag(ie+i-1)+...
                    (  r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i>j)
                diaginf(ie)=diaginf(ie)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            else
                diagsup(ie+1)=diagsup(ie+1)+( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            end
            % Versione compatta: più leggibile, ma meno efficiente si risparmia anche il 30% in fase di assemblaggio
            %                         Al(ie+i-1,ie+j-1) =  Al(ie+i-1,ie+j-1)+integrate( m.*femb_xi.*femb_xj+b.*femb_i.*femb_xj ...
            %                             + g.*femb_xi.*femb_j+s.*femb_i.*femb_j   ...
            %                             + funzionale_robin_up.*femb_i.*femb_j    ...
            %                             + funzionale_robin_down.*femb_i.*femb_j  ...
            %                             , w);
        end
        bl(ie+i-1) = bl(ie+i-1) + ( forc((ie-1)*nqnx+1:ie*nqnx)'.*femb_i )*(w );
    end
end
Al=spdiags( [diaginf,maindiag,diagsup], -1:1, nx, nx); % Assegno le diagonali alla matrice A
end