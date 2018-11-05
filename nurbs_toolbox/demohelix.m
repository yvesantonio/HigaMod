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
function demohelix 
% Demonstration of a 3D helical curve 
% 
 
% D.M. Spink 
% Copyright (c) 2000 
 
coefs =[ 6.0  0.0  6.0  1; 
        -5.5  0.5  5.5  1; 
        -5.0  1.0 -5.0  1; 
         4.5  1.5 -4.5  1; 
         4.0  2.0  4.0  1; 
        -3.5  2.5  3.5  1; 
        -3.0  3.0 -3.0  1; 
         2.5  3.5 -2.5  1; 
         2.0  4.0  2.0  1; 
        -1.5  4.5  1.5  1; 
        -1.0  5.0 -1.0  1; 
         0.5  5.5 -0.5  1; 
         0.0  6.0  0.0  1]'; 
knots = [0 0 0 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1 1 1]; 
 
crv = nrbmak(coefs,knots); 
nrbplot(crv,100); 
title('3D helical curve.'); 
grid on; 
