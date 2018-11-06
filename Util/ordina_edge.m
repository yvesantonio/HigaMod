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
function [e] = ordina_edge(e,t)

% riordino nodi sugli edge in modo che siano percorsi in  verso
% antiorario, secondo l'ordine imposto da t

nb = size(e,2);

for i = 1:nb
  % trova triangolo che si affaccia sull'edge i
  tri = find(sum(ismember(t(1:3,:),e(1:2,i))) == 2);
  % trova indici locali in tri dei due nodi dell'edge i
  i1 = find(e(1,i) == t(1:3,tri));
  i2 = find(e(2,i) == t(1:3,tri));
  % stabilisce se i due nodi dell'edge non sono antiorari
  % i nodi sono antiorari se sono 1,2 - 2,3 - 3,1  cioe' i2 = i1 + 1 mod 3
  if (i2 ~= (rem(i1,3) + 1))
    iold = e(1,i);
    e(1,i) = e(2,i);
    e(2,i) = iold;
  end
end
