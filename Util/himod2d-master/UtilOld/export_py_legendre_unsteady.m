function export_py_legendre_unsteady(size_mb,cutx,hx,u,nome_file,bc_up,bc_down,Map,t)

% ne intervalli e nx estremi
ne = round((cutx(2)-cutx(1))/hx);
nx = ne+1;

% Mesh FEM in x, nodi equispaziati
meshx    = zeros(3*ne,1);
meshx(1) = cutx(1);
meshx(2) = cutx(1)+hx/2;
meshx(3) = cutx(1)+hx;

for h=2:ne
    meshx( (h-1)*3 + 1) = meshx( (h-1)*3 );
    meshx( (h-1)*3 + 2) = meshx( (h-1)*3 ) + hx/2;
    meshx( (h-1)*3 + 3) = meshx( (h-1)*3 ) + hx;
end

M = 60;
meshy=linspace(-1,1,M)';
yvis=length(meshy);
[X,Yhat]=meshgrid(meshx,meshy);
Y=Map.invhaty(t,X,Yhat);

sol=zeros(3*nx,yvis);

mb = new_modal_basis_legendre(size_mb,meshy,bc_up,bc_down);

for h=1:ne
    for k=1:M
        for imb=1:size_mb
            sol(3*(h-1) + 1 , k)=sol(3*(h-1) + 1, k ) + u(h+(imb-1)*nx)*mb(k,imb);
            sol(3*(h-1) + 2 , k)=sol(3*(h-1) + 2, k ) + 1/2*( u(h+(imb-1)*nx) + u(h+1+(imb-1)*nx) ) *mb(k,imb);
            sol(3*(h-1) + 3 , k)=sol(3*(h-1) + 3, k ) + u(h+1+(imb-1)*nx)*mb(k,imb);
        end
    end
end

orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    nome_file ...
    );

nome_file=strcat('./',titolo,'.out');

fid =fopen(nome_file,'w');
fprintf(fid,'            %u\n',3*ne);
fprintf(fid,'             %u\n',M);

for k=1:ne
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',meshx( (k-1)*3 + 1),Y(j,(k-1)*3+1),sol((k-1)*3 + 1,j));
    end
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',meshx( (k-1)*3 + 2),Y(j,(k-1)*3+2),sol((k-1)*3 + 2,j));
    end
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',meshx( (k-1)*3 + 3),Y(j,(k-1)*3+3),sol((k-1)*3 + 3,j));
    end
end

fclose(fid);
return