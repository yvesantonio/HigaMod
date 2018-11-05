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
% Construction of a bilinearly blended Coons surface 
% 
 
%  D.M. Spink 
%  Copyright (c) 2000. 
 
% Boundary curve 1 
pnts = [ 0.0  3.0  4.5  6.5 8.0 10.0; 
         0.0  0.0  0.0  0.0 0.0  0.0;  
         2.0  2.0  7.0  4.0 7.0  9.0];    
crv1 = nrbmak(pnts, [0 0 0 1/3 0.5 2/3 1 1 1]); 
 
% Boundary curve 2 
pnts= [ 0.0  3.0  5.0  8.0 10.0; 
        10.0 10.0 10.0 10.0 10.0; 
        3.0  5.0  8.0  6.0 10.0]; 
crv2 = nrbmak(pnts, [0 0 0 1/3 2/3 1 1 1]); 
 
% Boundary curve 3 
pnts= [ 0.0 0.0 0.0 0.0; 
        0.0 3.0 8.0 10.0; 
        2.0 0.0 5.0 3.0]; 
crv3 = nrbmak(pnts, [0 0 0 0.5 1 1 1]); 
 
% Boundary curve 4 
pnts= [ 10.0 10.0 10.0 10.0 10.0; 
        0.0   3.0  5.0  8.0 10.0; 
        9.0   7.0  7.0 10.0 10.0]; 
crv4 = nrbmak(pnts, [0 0 0 0.25 0.75 1 1 1]); 
 
srf = nrbcoons(crv1, crv2, crv3, crv4); 
 
% Draw the surface 
nrbplot(srf,[20 20]); 
title('Construction of a bilinearly blended Coons surface.'); 
