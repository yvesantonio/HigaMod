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
function deg = rad2deg(rad) 
%  
% Function Name: 
%  
%   rad2deg - Convert radians to degrees. 
%  
% Calling Sequence: 
%  
%   rad = rad2deg(deg); 
%  
% Parameters: 
%  
%   rad		: Angle in radians. 
% 
%   deg		: Angle in degrees. 
%  
% Description: 
%  
%   Convenient utility function for converting radians to degrees, which are 
%   often the required angular units for functions in the NURBS toolbox. 
%  
% Examples: 
%  
%   Convert 0.3 radians to degrees 
%  
%   rad = deg2rad(0.3); 
 
%  D.M. Spink 
%  Copyright (c) 2000. 
 
deg = 180.0*rad/pi; 
