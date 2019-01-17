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
function write_msh(p,e,t,filename)

t(4,:) = t(4,:) - 1;
e = ordina_edge(e,t);
tmp = zeros(1,size(p,2));
tmp(e(1,:)) = e(5,:);
p = [p;tmp];

% write msh file
fid = fopen([filename,'.msh'],'w');
np = size(p,2);
nt = size(t,2);
ne = size(e,2);
fprintf(fid,'%d %d %d \n',[np nt ne]);
fprintf(fid,'%f %f %d \n',p);
fprintf(fid,'%d %d %d %d \n',t);
fprintf(fid,'%d %d %d \n',e([1:2,5],:));
fclose(fid);
