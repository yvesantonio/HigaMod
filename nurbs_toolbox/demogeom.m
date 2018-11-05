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
function demogeom 
% Demonstration of how to construct a 2D geometric   
% shape from a piece-wise set of NURBSs 
% 
 
% D.M. Spink 
% Copyright (c) 2000. 
 
coefs = [0.0 7.5 15.0 25.0 35.0 30.0 27.5 30.0; 
         0.0 2.5  0.0 -5.0  5.0 15.0 22.5 30.0]; 
knots = [0.0 0.0 0.0 1/6 1/3 1/2 2/3 5/6 1.0 1.0 1.0]; 
 
% Geometric boundary curve 
geom = [ 
nrbmak(coefs,knots) 
nrbline([30.0 30.0],[20.0 30.0]) 
nrbline([20.0 30.0],[20.0 20.0]) 
nrbcirc(10.0,[10.0 20.0],1.5*pi,0.0) 
nrbline([10.0 10.0],[0.0 10.0]) 
nrbline([0.0 10.0],[0.0 0.0]) 
nrbcirc(5.0,[22.5 7.5]) 
]; 
 
ng = length(geom); 
for i = 1:ng 
  nrbplot(geom(i),50); 
  hold on; 
end 
hold off; 
axis equal; 
title('2D Geometry formed by a series of NURBS curves'); 
