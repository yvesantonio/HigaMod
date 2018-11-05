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
function demodersrf 
% Demonstrates the construction of a general 
% surface derivatives. 
% 
 
% D.M. Spink 
% Copyright (c) 2000 
 
% make and draw a test surface 
srf = nrbtestsrf; 
p = nrbeval(srf,{linspace(0.0,1.0,20) linspace(0.0,1.0,20)}); 
h = surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:))); 
set(h,'FaceColor','blue','EdgeColor','blue'); 
title('First derivatives over a test surface.'); 
 
npts = 5; 
tt = linspace(0.0,1.0,npts); 
 
dsrf = nrbderiv(srf); 
 
[p1, dp] = nrbdeval(srf, dsrf, {tt, tt}); 
 
up2 = vecnorm(dp{1}); 
vp2 = vecnorm(dp{2}); 
 
hold on; 
plot3(p1(1,:),p1(2,:),p1(3,:),'ro'); 
h1 = quiver3(p1(1,:),p1(2,:),p1(3,:),up2(1,:),up2(2,:),up2(3,:)); 
h2 = quiver3(p1(1,:),p1(2,:),p1(3,:),vp2(1,:),vp2(2,:),vp2(3,:)); 
set(h1,'Color','black'); 
set(h2,'Color','black'); 
 
hold off; 
 
