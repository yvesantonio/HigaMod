function Al=AssembleBXREL(Map,wyq,mb_k,mb_s,mb_yk,ne,mesh_wx,nqnx,FEM_basis,FEM_basis_x,FEM_basisP,k)

r00v = -( Map.Jac.*Map.D) *(mb_yk.*mb_s.*wyq );
r10v = -( Map.Jac )*(mb_k.*mb_s.*wyq );

nx=ne+1;
if(k>2)
Al=zeros(nx,nx+ne);
for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for j = 1:3     % Per ogni funzione di forma EF vel
        femb_j  = FEM_basis  ( j, :);
        femb_xj = FEM_basis_x( j, :);
        for i = 1:2  % Per ogni funzione di forma EF press
            %selezione dei pezzi corretti per le basi FEM
            femb_i  = FEM_basisP  ( i, :);
            Al(ie-1+i,2*ie-2+j) =  Al(ie-1+i,2*ie-2+j)+(r00.*femb_i.*femb_j+r10.*femb_i.*femb_xj)*w;
        end
    end
end
else
    
Al=zeros(nx,1);
ie=1;
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    j=1;
        femb_j  = FEM_basis  ( j, :);
        femb_xj = FEM_basis_x( j, :);
        for i = 1:2  % Per ogni funzione di forma EF press
            %selezione dei pezzi corretti per le basi FEM
            femb_i  = FEM_basisP  ( i, :);
            Al(ie-1+i,2*ie-2+j) =  Al(ie-1+i,2*ie-2+j)+(r00.*femb_i.*femb_j+r10.*femb_i.*femb_xj)*w;
        end
end
