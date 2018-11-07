function Al=AssembleBY(Map,wyq,mb_s,mb_yk,ne,mesh_wx,nqnx,FEM_basis,FEM_basisP)

r00v = -( Map.Jac.*Map.J) *(mb_yk.*mb_s.*wyq );

nx=ne+1;
Al=zeros(nx,nx+ne);
for ie = 1 : ne  % Per ogni intervallo
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx( (ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    for i = 1:2                  % Per ogni funzione di forma EF pressione
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basisP  ( i, :);
        for j = 1:3             % Per ogni funzione di forma EF velocità
            femb_j  = FEM_basis  ( j, :);
            Al(ie-1+i,2*ie-2+j) = Al(ie-1+i,2*ie-2+j) + (r00.*femb_i.*femb_j)*w;
        end
    end
end
end
