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
% Demonstration of ruled surface construction. 
% 
 
%  D.M. Spink 
%  Copyright (c) 2000. 
 
% clear all recorded graphics and the current window 
title('Ruled surface construction from two NURBS curves.'); 
 
crv1 = nrbtestcrv; 
crv2 = nrbtform(nrbcirc(4,[4.5;0],pi,0.0),vectrans([0.0 4.0 -4.0])); 
srf = nrbruled(crv1,crv2); 
nrbplot(srf,[40 20]); 
title('Ruled surface construction from two NURBS curves.'); 
