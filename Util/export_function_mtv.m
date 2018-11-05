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
