function export_py_fun(func,domain,xStep,fileName)
% exporting in a format suitable for py scripts
% the function FUNC define on DOMAIN
% on a grid with spacing xStep.
% The filenam is in the format
%    DATE_TIME_FILENAME
% TODO this way of exporing in python is a very hold one.
% several values are exported twice. CHANGE IT!


% creating a filename
% format
% DATA_ORA_FILENAME
orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    fileName ...
    );
fileName=strcat('./',titolo,'.out');


% visualizzation grid
ne=round( ( domain(2)-domain(1) ) / xStep );
nx = ne+1;
M=60;

% writing header information
fid =fopen(fileName,'w');
fprintf(fid,'            %u\n',3*ne);
fprintf(fid,'             %u\n',M);

% Mesh FEM in x, nodi equispaziati
meshx = zeros(3*nx,1);
% first element
meshx(1) = domain(1);
meshx(2) = domain(1)+xStep/2;
meshx(3) = domain(1)+xStep;
% all other elements
for h=2:ne
    meshx( (h-1)*3 + 1) = meshx( (h-1)*3 );
    meshx( (h-1)*3 + 2) = meshx( (h-1)*3 ) + xStep/2;
    meshx( (h-1)*3 + 3) = meshx( (h-1)*3 ) + xStep;
end
meshy=linspace(0,1,M)';

% starting writing phase:
% format for each line 
%        x      y      f(x,y)
for k=1:ne
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',...
            meshx( (k-1)*3 + 1),meshy(j),func(meshx( (k-1)*3 + 1),meshy(j)));
    end
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',...
            meshx( (k-1)*3 + 2),meshy(j),func(meshx( (k-1)*3 + 2),meshy(j)));
    end
    for j=1:M
        fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',...
            meshx( (k-1)*3 + 3),meshy(j),func(meshx( (k-1)*3 + 3),meshy(j)));
    end
end
fclose(fid);
return
