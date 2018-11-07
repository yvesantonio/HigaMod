function Al = AssembleVelXY(Map,wyq,DMbKX,MbRY,DMbRY,FEM_basis,FEM_basis_x,ne,...
    mesh_wx,nqnx,nu)

% per ogni x (vettorialmente)
% calcola r^st_kr
r01v = ( nu*Map.Jac.*Map.J)*( DMbKX.*MbRY .*wyq );
r00v = ( nu*Map.Jac.*Map.D.*Map.J ) *(DMbKX.*DMbRY.*wyq );

nx=ne+1;

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag =   zeros( ne + nx, 1);
diagsup  =   zeros( ne + nx, 1); % Spdiags taglierà il primo elemento
diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglierà l'ultimo
diagsuper  = zeros( ne + nx, 1); % Spdiags taglierà i primi due elementi
diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglierà gli ultimi due

for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r01 = r01v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:3                  % Per ogni funzione di forma EF test
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basis  (i,:);
        femb_xi = FEM_basis_x(i,:);
        for j = 1:3             % Per ogni funzione di forma EF sol
            
            femb_j  = FEM_basis  ( j,:);
            
            if(i==j)
                maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+...
                    (  r01.*femb_xi.*femb_j...
                    +r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j+1 && i==2)
                diaginf(2*ie-1)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+1 && i==3)
                diaginf(2*ie)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+2)
                diaginfer(2*ie-1)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-2)
                diagsuper(2*ie+1)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-1 && i==2)
                diagsup(2*ie+1)=( ...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j-1 && i==1)
                diagsup(2*ie)=( ...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            end
        end
    end
end
Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A
end
