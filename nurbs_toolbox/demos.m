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
function tbxStruct=Demos 
% Demos   Demo List infomation for NURBS Toolbox 
 
% D.M. Spink 
% Copyright (c) 2000 
 
if nargout == 0; demo toolbox; return; end 
 
tbxStruct.Name='nurbs'; 
tbxStruct.Type='Toolbox'; 
 
tbxStruct.Help= ... 
        {' The NURBS Toolbox provides commands for the construction' 
         ' and use of Non-Uniform Rational B-Splines (NURBS).' 
         ' '}; 
 
tbxStruct.DemoList={ 
                '3D Line','demoline', 
                'Rectangle','demorect', 
                'Circle and Arc','democirc', 
                'Helix','demohelix', 
                'Ellipse', 'demoellip', 
                'Test Curve','democurve', 
                'Test Surface','demosurf', 
                'Bilinear Surface','demo4surf', 
                'Knot Insertion','demokntins', 
                'Degree Elevation','demodegelev', 
                'Curve derivatives','demodercrv', 
                'Surface derivatives','demodersrf' 
                'Cylinderical arc','democylind', 
                'Ruled Surface','demoruled', 
                'Coons Surface','democoons', 
                'Test Extrusion','demoextrude', 
                'Test Revolution - Cup','demorevolve', 
                'Test Revolution - Ball & Torus','demotorus', 
                '2D Geomtery','demogeom'}; 
 
