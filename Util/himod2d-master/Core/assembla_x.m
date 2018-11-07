function [Al,bl,a_ril,b_ril]=assembla_x(imb,kmb,size_fb,degree,n,nqnx, ...
    mesh_wx,wyq,mb_i,mb_yi,mb_k,mb_yk, FEM_basis, FEM_basis_x,L,bc_up,...
    bc_down,dato_uploc,dato_downloc,posizione_dominio,Coeff_forma,gamma,Computed,nd,coupling)

coeffrobin=Coeff_forma.coeffrobin;

%Calcolo del rilevamento
[a_ril,b_ril]=coeff_ril(bc_up,bc_down,dato_uploc,dato_downloc,Coeff_forma,L(1,1));
aus=@(x,y) a_ril*y+b_ril;
rilevamento_c=aus(0,Computed.yq)';

% ... e della forzante.
 forc  = integrate( L.*Computed.force_c - a_ril.*Computed.beta2_c - Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);

%Calcolo dei coefficienti r_{ik}^{st}
muv       = integrate( L.*Computed.mu_c          , mb_k .*mb_i .*wyq );
muv_y     = integrate( L.*Computed.mu_c          , mb_yk.*mb_yi.*wyq );
betav     = integrate( L.*Computed.beta1_c       , mb_k .*mb_i .*wyq );
betav_y   = integrate( L.*Computed.beta2_c       , mb_yk.*mb_i .*wyq );
sigmav    = integrate( L.*Computed.sigma_c       , mb_k .*mb_i .*wyq );


%Inizializzazione del blocco di matrice
Al = sparse(size_fb,size_fb);
bl = zeros(size_fb,1);

%determina le basi valutate solo nei punti della frontiera UP e DOWN
base_mod_aus = new_modal_basis( max(imb,kmb), [0,1], bc_up,bc_down,Coeff_forma,L(1,1));


if(strcmp(bc_up,'rob'))
    funzionale_robin_up   = coeffrobin*base_mod_aus(2, imb)*base_mod_aus(2, kmb);
else
    funzionale_robin_up   = 0;
end

if(strcmp(bc_down,'rob'))
    funzionale_robin_down = coeffrobin*base_mod_aus(1, imb)*base_mod_aus(1, kmb);
else
    funzionale_robin_down = 0;
end

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag = zeros( size_fb, 1);
diagsup  = zeros( size_fb, 1); % Spdiags taglier?? il primo elemento
diaginf  = zeros( size_fb, 1); % Spdiags invece taglier?? l'ultimo

for ie = 1 : n  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    m  = muv   ( (ie-1)*nqnx+1 : ie*nqnx)';
    my = muv_y ( (ie-1)*nqnx+1 : ie*nqnx)';
    b  = betav ( (ie-1)*nqnx+1 : ie*nqnx)';
    by = betav_y ( (ie-1)*nqnx+1 : ie*nqnx)';
    s  = sigmav( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:degree+1                  % Per ogni funzione di forma EF
        for j = 1:degree+1              % Per ogni funzione di forma EF
            
            %selezione dei pezzi corretti per le basi FEM
            femb_i  = FEM_basis  ( (ie-1)*nqnx+1:ie*nqnx, i)';
            femb_xi = FEM_basis_x( (ie-1)*nqnx+1:ie*nqnx, i)';
            femb_j  = FEM_basis  ( (ie-1)*nqnx+1:ie*nqnx, j)';
            femb_xj = FEM_basis_x( (ie-1)*nqnx+1:ie*nqnx, j)';
            
            if(i==j)
                maindiag(ie+i-1)=maindiag(ie+i-1)+integrate(...
                    m.*femb_xi.*femb_xj + my.*femb_i.*femb_j...
                    + b.*femb_i.*femb_xj + by.*femb_i.*femb_j...
                    + s.*femb_i.*femb_j   ...
                    + funzionale_robin_up.*femb_i.*femb_j    ...
                    + funzionale_robin_down.*femb_i.*femb_j  ...
                    , w);
                
            elseif(i>j)
                diaginf(ie)=diaginf(ie)+integrate(...
                    m.*femb_xi.*femb_xj + my.*femb_i.*femb_j...
                    + b.*femb_i.*femb_xj + by.*femb_i.*femb_j...
                    + s.*femb_i.*femb_j   ...
                    + funzionale_robin_up.*femb_i.*femb_j    ...
                    + funzionale_robin_down.*femb_i.*femb_j  ...
                    , w);
            else
                diagsup(ie+1)=diagsup(ie+1)+integrate( ...
                    m.*femb_xi.*femb_xj + my.*femb_i.*femb_j...
                    + b.*femb_i.*femb_xj + by.*femb_i.*femb_j...
                    + s.*femb_i.*femb_j   ...
                    + funzionale_robin_up.*femb_i.*femb_j    ...
                    + funzionale_robin_down.*femb_i.*femb_j  ...
                    , w);
            end
            % Versione compatta: pi?? leggibile, ma meno efficiente si risparmia anche il 30% in fase di assemblaggio
            %                         Al(ie+i-1,ie+j-1) =  Al(ie+i-1,ie+j-1)+integrate( m.*femb_xi.*femb_xj+b.*femb_i.*femb_xj ...
            %                             + g.*femb_xi.*femb_j+s.*femb_i.*femb_j   ...
            %                             + funzionale_robin_up.*femb_i.*femb_j    ...
            %                             + funzionale_robin_down.*femb_i.*femb_j  ...
            %                             , w);
        end
        bl(ie+i-1) = bl(ie+i-1) + integrate( forc((ie-1)*nqnx+1:ie*nqnx)'.*femb_i , w );
    end
end
Al=spdiags( [diaginf,maindiag,diagsup], -1:1, size_fb, size_fb); % Assegno le diagonali alla matrice A

if(strcmp(coupling,'RR'))
    if(nd>1)
        %Contributo alla forma bilineare dovuto ai termini del Robin-Robin
        mu=Computed.mu_c(1,1); %Attenzione se mu non ?? costante va cambiato, e serve un vettore che dipende da y valutato in X_left
        if(posizione_dominio~=1) %Il primo dominio non ha bordo sinistro, l'ultimo non ha bordo destro.
            Al(1,1)=Al(1,1) + integrate(L(1,1)*mu*gamma.L*mb_k'.*mb_i',wyq); % se si vuole usare con domini non larghi L controllare qui e sotto
        end
        %Sarebbe un mu diverso dalla riga sopra
        if(posizione_dominio~=-1)
            Al(end,end)=Al(end,end) + integrate(L(1,1)*mu*gamma.R*mb_k'.*mb_i',wyq);
        end
    end
end

end
