function export_py(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,nome_file,L,a)
if nargin < 11
    L=1;
    a=0;
    warning('last argument not specified, are you sure?')
end
nd=size(size_mb,2);

orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    nome_file ...
    );

nome_file=strcat('./',titolo,'.out');
ne=0;
for i=1:nd
    ne=ne+round((cutx(i+1)-cutx(i))/hx(i));
end
M=60;

fid =fopen(nome_file,'w');
fprintf(fid,'            %u\n',3*ne);
fprintf(fid,'             %u\n',M);

for i=1:nd
% ne intervalli e nx estremi
ne = round((cutx(i+1)-cutx(i))/hx(i));
nx = ne+1;

% Mesh FEM in x, nodi equispaziati
meshx{i} = zeros(3*nx,1);
meshx{i}(1) = cutx(i);
meshx{i}(2) = cutx(i)+hx(i)/2;
meshx{i}(3) = cutx(i)+hx(i);

for h=2:ne
    meshx{i}( (h-1)*3 + 1) = meshx{i}( (h-1)*3 );
    meshx{i}( (h-1)*3 + 2) = meshx{i}( (h-1)*3 ) + hx(i)/2;
    meshx{i}( (h-1)*3 + 3) = meshx{i}( (h-1)*3 ) + hx(i);
end

meshyref=linspace(0,1,M)';
meshy=linspace(a,a+L,M)';
yvis=length(meshy);

sol{i}=zeros(3*nx,yvis);

mb = new_modal_basis(size_mb(i),meshyref,bc_up{i},bc_down{i},Coeff_forma,L);

for h=1:ne
    for k=1:M
        for imb=1:size_mb(i)
            sol{i}(3*(h-1) + 1 , k)=sol{i}(3*(h-1) + 1, k ) + u{i}(h+(imb-1)*nx)*mb(k,imb);
            sol{i}(3*(h-1) + 2 , k)=sol{i}(3*(h-1) + 2, k ) + 1/2*(u{i}(h+(imb-1)*nx) + u{i}(h+1+(imb-1)*nx) ) *mb(k,imb);
            sol{i}(3*(h-1) + 3 , k)=sol{i}(3*(h-1) + 3, k ) + u{i}(h+1+(imb-1)*nx)*mb(k,imb);
        end
        sol{i}(3*(h-1) + 1,k)=sol{i}(3*(h-1)+1,k)+ a_ril(i)*meshy(k)+b_ril(i);
        sol{i}(3*(h-1) + 2,k)=sol{i}(3*(h-1)+2,k)+ a_ril(i)*meshy(k)+b_ril(i);
        sol{i}(3*(h-1) + 3,k)=sol{i}(3*(h-1)+3,k)+ a_ril(i)*meshy(k)+b_ril(i);
    end
end


for k=1:ne
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',meshx{i}( (k-1)*3 + 1),meshy(j),sol{i}((k-1)*3 + 1,j));
    end
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',meshx{i}( (k-1)*3 + 2),meshy(j),sol{i}((k-1)*3 + 2,j));
    end
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',meshx{i}( (k-1)*3 + 3),meshy(j),sol{i}((k-1)*3 + 3,j));
    end
end
end
fclose(fid);
return