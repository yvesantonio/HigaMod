function export_function_mtv(nome_file,f,x,y)
%
%	function export_function_mtv(nome_file,f,x,y)
%   
%	This function let you export a function handle in a format suitable for "plotmtv".
%	This function gets 4 four arguments: the (partial) name of the file where you are going to save the output.
%   The file will be saved in the local directory (where you are running the code). 
%   f is meant to be a function handle (@-matlabfunction) the will be evaluated in a grid which will be 
%   the product of a grid in x direction (x) and a grid in y direction (y).
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

for i=1:length(x)-1
    for j=1:length(y)-1
        fprintf(fid,'%-e  %-e %-e\n',x(i),y(j),  f(x(i),y(j)));
        fprintf(fid,'%-e  %-e %-e\n',x(i+1),y(j),f(x(i+1),y(j)));
        fprintf(fid,'%-e  %-e %-e\n',x(i),y(j+1),f(x(i),y(j+1)));
        fprintf(fid,'\n');
        fprintf(fid,'%-e  %-e %-e\n',x(i+1),y(j),  f(x(i+1),y(j)));
        fprintf(fid,'%-e  %-e %-e\n',x(i+1),y(j+1),f(x(i+1),y(j+1)));
        fprintf(fid,'%-e  %-e %-e\n',x(i),y(j+1),f(x(i),y(j+1)));
        fprintf(fid,'\n');
    end
end

fprintf(fid,'$ END');
fclose(fid);
