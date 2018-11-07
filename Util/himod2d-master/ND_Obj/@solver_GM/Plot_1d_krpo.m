function []=Plot_1d_krpo(obj,num)
%funzione ausiliaria per plottare la soluzione ottenuta con modello tipo
%Krpo

for i=1:length(obj.nodi)
    sol_display(i)=obj.Solution{num,2}(2*i-1);
end

Z=obj.profile{num}(:,2)*sol_display/sol_display(1);

if(obj.treD)
    surf(obj.nodi,obj.profile{num}(:,1),Z,'EdgeColor', 'none');
else
    contourf(obj.nodi,obj.profile{num}(:,1),Z,obj.levels);
end

% N=100;
% h=1./N;
% meshx = zeros(3*6*N,1);
% for i=1:6*N
% 	meshx[3*i]=i*h;
% 	meshx[3*i+1]=i*h+h/2.;
% 	meshx[3*i+2]=i*h+h;
% end

%griglia di export 100y 301x
[n,m] = size(Z)
f6 = fopen('matrice.out','w');
for i=1:n
    for j=1:m
        fprintf(f6,'%8.4f\t',Z(i,j));
    end
    fprintf(f6,'\n');
end
fclose(f6);

f3 = fopen('xcoordinate.out','w');
fprintf(f3,'%8.4f \n', obj.nodi);
fclose(f3);

ycord = zeros(size(obj.profile{num}(:,1)));
ycord = obj.profile{num}(:,1);
f4 = fopen('ycoordinate.out','w');
fprintf(f4,'%8.4f \n', ycord);
fclose(f4);
end