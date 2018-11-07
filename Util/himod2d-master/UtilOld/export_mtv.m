function export_mtv(nome_file,sol_mtv,x,meshy)
%   function export_mtv(nome_file,sol_mtv,x,meshy)
%
%   This function export the solution obtained by plot_solution_DD in a suitable "plotmtv" format.
%	The first argument: nome_file is only a part of the name of the file that will be saved in the local dir.
%   export_mtv works not only for single domains but also in multiple domains.
orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    nome_file ...
    );
nome_file=strcat('./',titolo,'.mtv');

fid =fopen(nome_file,'w');
fprintf(fid,'$ DATA=CONTCURVE\n');
fprintf(fid,'%% contstyle=3\n');
fprintf(fid,'%% toplabel="" \n');
fprintf(fid,'%% xlabel="" \n');
fprintf(fid,'%% ylabel="" \n');
fprintf(fid,'%% zlabel="" \n');
fprintf(fid,'%% axisguides=false\n');
fprintf(fid,'%% contfill\n');



for dominio=1:length(sol_mtv)
    
sizex=length(x{dominio});
sizey=length(meshy);

y=zeros(sizex,sizey);

for i=1:sizex
    y(i,:)=meshy;
end


[plm,pln]=size(sol_mtv{dominio});
N_couple_elem = (plm-1)*(pln-1);
for k=1:N_couple_elem
    ipl=mod(k-1,plm-1)+1;
    jpl=floor((k-1)/(plm-1))+1;
    fprintf(fid,'%-e  %-e %-e\n',x{dominio}(ipl),y(ipl,jpl),sol_mtv{dominio}(ipl,jpl));
    fprintf(fid,'%-e  %-e %-e\n',x{dominio}(ipl+1),y(ipl+1,jpl),sol_mtv{dominio}(ipl+1,jpl));
    fprintf(fid,'%-e  %-e %-e\n',x{dominio}(ipl),y(ipl,jpl+1),sol_mtv{dominio}(ipl,jpl+1));
    fprintf(fid,'\n');
    fprintf(fid,'%-e  %-e %-e\n',x{dominio}(ipl+1),y(ipl+1,jpl),sol_mtv{dominio}(ipl+1,jpl));
    fprintf(fid,'%-e  %-e %-e\n',x{dominio}(ipl+1),y(ipl+1,jpl+1),sol_mtv{dominio}(ipl+1,jpl+1));
    fprintf(fid,'%-e  %-e %-e\n',x{dominio}(ipl),y(ipl,jpl+1),sol_mtv{dominio}(ipl,jpl+1));
    fprintf(fid,'\n');
end
end
fprintf(fid,'$ END');
fclose(fid);
return
