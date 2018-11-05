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
function demodercrv 
% Demonstrates the construction of a general 
% curve and determine of the derivative. 
% 
 
% D.M. Spink 
% Copyright (c) 2000 
 
% make and draw nurbs test curve 
crv = nrbtestcrv; 
nrbplot(crv,48); 
title('First derivatives along a test curve.'); 
 
npts = 9; 
tt = linspace(0.0,1.0,npts); 
 
dcrv = nrbderiv(crv); 
 
% first derivative 
[p1, dp] = nrbdeval(crv,dcrv,tt); 
 
p2 = vecnorm(dp); 
 
hold on; 
plot(p1(1,:),p1(2,:),'ro'); 
h = quiver(p1(1,:),p1(2,:),p2(1,:),p2(2,:),0); 
set(h,'Color','black'); 
hold off; 
