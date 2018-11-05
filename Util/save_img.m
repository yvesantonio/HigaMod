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
function save_img(fig,name,varargin)
%
%  function save_img(fig,name,varargin)
%
%  senza extra-argomenti salva in png altrimenti salva nel formato dentro
%  name è parte del nome del file che viene salvato, il resto è data. 
% il tutto viene salvato dentro al cartella IMG
% fig è la figura da salvare
%  varargin, per ora c'è solo eps, ma funziona male oltre a png.
%
if(isempty(varargin))
    formato='png';
elseif(strcmp(varargin{1},'png'))
    formato='png';
elseif(strcmp(varargin{1},'eps'))
    formato='eps';
else
    display('in save_img: formato non riconosciuto')
    return
end
orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];

titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    name ...
    );

if(strcmp(formato,'eps'))
    filename=strcat('./',titolo,'.eps');
    print(fig, '-depsc2', '-r300', filename);
else
    filename=strcat('./',titolo,'.png');
    print(fig, '-dpng', '-r300', filename);
end
end
