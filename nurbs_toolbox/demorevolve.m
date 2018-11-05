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
function demorevolve 
% Demonstration of surface construction by revolving a 
% profile curve. 
 
% D.M. Spink 
% Copyright (c) 2000 
 
% Construct a test profile in the x-z plane  
pnts = [3.0 5.5 5.5 1.5 1.5 4.0 4.5; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.5 1.5 4.5 3.0 7.5 6.0 8.5]; 
crv = nrbmak(pnts,[0 0 0 1/4 1/2 3/4 3/4 1 1 1]); 
 
% rotate and vectrans by some arbitrary amounts. 
xx = vecrotz(deg2rad(25))*vecroty(deg2rad(15))*vecrotx(deg2rad(20)); 
nrb = nrbtform(crv,vectrans([5 5])*xx); 
 
% define axes of rotation 
pnt = [5 5 0]'; 
vec = xx*[0 0 1 1]'; 
srf = nrbrevolve(nrb,pnt,vec(1:3)); 
 
% make and draw nurbs curve 
p = nrbeval(srf,{linspace(0.0,1.0,20) linspace(0.0,1.0,20)}); 
surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:))); 
title('Construct of a 3D surface by revolution of a curve.'); 
shading interp; 
colormap(copper); 
axis equal; 
 
