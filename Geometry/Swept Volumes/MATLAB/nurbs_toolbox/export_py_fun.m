%-+--------------------------------------------------------------------
% HiModlab is a general purpose Hierarchical model reduction library.
% Copyright (C) 2006-2017  by the HiMod authors (see authors.txt).
%
% This file is part of the HiModlab library.
%
% The HiModlab library is free software: you can use it, redistribute
% it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either
% version 2 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Acknowledgement and Disclaimer
%
% The development of HiModLab was supported by the US National Science
% Foundation (NSF Project, DMS 1419060 ). Any opinions, findings and conclusions
% or recomendations expressed in the source code and documentation are
% those of the authors and do not necessarily reflect the views of the
% National Science Foundation (NSF).
%
%
%
%-+--------------------------------------------------------------------
function export_py_fun(u,cutx,hx,nome_file)
nd=1;

orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    nome_file ...
    );
nome_file=strcat('./',titolo,'.out');

% griglia visualizzazione
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
    
    meshy=linspace(0,1,M)';
    yvis=length(meshy);
    
    for k=1:ne
        for j=1:M
            fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',...
                meshx{i}( (k-1)*3 + 1),meshy(j),u(meshx{i}( (k-1)*3 + 1),meshy(j)));
        end
        for j=1:M
            fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',...
                meshx{i}( (k-1)*3 + 2),meshy(j),u(meshx{i}( (k-1)*3 + 2),meshy(j)));
        end
        for j=1:M
            fprintf(fid,'       %1.7E       %1.7E       %1.7E\n',...
                meshx{i}( (k-1)*3 + 3),meshy(j),u(meshx{i}( (k-1)*3 + 3),meshy(j)));
        end
    end
end
fclose(fid);
return
