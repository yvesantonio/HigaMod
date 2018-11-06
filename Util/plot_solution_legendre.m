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
function [aus_figure,sol,meshx,meshy] = plot_solution_legendre(size_mb,mesh_x, u,Map,bc_up,bc_down)

persistent calls current_figure

if isempty(calls)
    calls = 0;
    current_figure=figure('visible','off');
%     current_figure=figure('visible','on');
else
    current_figure=figure('visible','off');
%       current_figure=figure('visible','on');
end

colormap jet;
cmin = -0.04;
cmax = +0.1;

% ne intervalli e nx estremi
nx = length(mesh_x);
%ne=nx-1;

M = 30;
meshy=linspace(-1,1,M)';
yvis=length(meshy);

sol=zeros(nx,yvis);

mb = new_modal_basis_legendre(size_mb,meshy,bc_up,bc_down);

for h=1:nx
    for k=1:M
        for imb=1:size_mb
            sol(h,k)=sol(h,k)+ u(h+(imb-1)*nx)*mb(k,imb);
        end
    end
end

minimo=min(min(sol));
massimo=max(max(sol));

display(minimo)
display(massimo)

%contourLines = linspace(minimo, massimo, 20);
contourLines = linspace(cmin,cmax, 20);

[X,Yhat]=meshgrid(mesh_x,meshy);
Y=Map.invhaty(X,Yhat);
warning('OFF','MATLAB:contourf:ConstantData');
contourf(X,Y,sol(:,:)',contourLines);
grid on
warning('ON','MATLAB:contourf:ConstantData');
hold off
set(gca, ...
    'FontUnits','points', ...
    'FontSize',11);

axis equal
axis([mesh_x(1) mesh_x(end) -1.5 1.5]);

caxis([cmin,cmax]);
calls = calls + 1;
aus_figure=current_figure;
%export_fig(current_figure, 'filename',['himod_m=',num2str(size_mb)],'-eps','-transparent')
end
