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
close all
clear all
clc

ne=20;
nx=ne+1;
nqnx=4;
mesh_x=linspace(0,1,nx);
%Setting up a quadrature rule for the x direction.
[~, node_x, weights_x] = gausslegendre( nqnx );
quad_nodes_mesh_x = zeros( ne*nqnx, 1);
quad_weights_mesh_x = zeros( ne*nqnx, 1);
for i=1:ne
    [ ...
        quad_nodes_mesh_x( (i-1)*nqnx+1 : i*nqnx ), ...
        quad_weights_mesh_x( (i-1)*nqnx+1 : i*nqnx ) ...
        ] = quadrature_rule( mesh_x(i), mesh_x(i+1), node_x, weights_x);
end
mesh_wx=quad_weights_mesh_x;
[ FEM_basis , FEM_basis_x]=new_fem_basis(2,mesh_x(1:2),quad_nodes_mesh_x(1:4));
% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag =   zeros( ne + nx, 1);
diagsup  =   zeros( ne + nx, 1); % Spdiags taglierà il primo elemento
diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglierà l'ultimo
diagsuper  = zeros( ne + nx, 1); % Spdiags taglierà i primi due elementi
diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglierà gli ultimi due
bl         = zeros( ne + nx, 1);

for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r11 = 1;
    r10 = 0;
    r01 = 0;
    r00 = 0;
    
    for i = 1:3                  % Per ogni funzione di forma EF test
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basis  ( :, i)';
        femb_xi = FEM_basis_x( :, i)';
        for j = 1:3             % Per ogni funzione di forma EF sol
            
            femb_j  = FEM_basis  ( :, j)';
            femb_xj = FEM_basis_x( :, j)';
            
            if(i==j)
                maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+...
                    (  r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j+1 && i==2)
                diaginf(2*ie-1)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+1 && i==3)
                diaginf(2*ie)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+2)
                diaginfer(2*ie-1)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-2)
                diagsuper(2*ie+1)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-1 && i==2)
                diagsup(2*ie+1)=( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j-1 && i==1)
                diagsup(2*ie)=( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            end
        end
        bl(2*(ie-1)+i) = bl(2*(ie-1)+i) + ( 1.*femb_i )*w;
    end
end
Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A

Al(1,1)=1e30;
bl(1)=1e30;
Al(end,end)=1e30;
bl(end)=1e30;

spy(Al)
figure;
u=Al\bl;
h=(mesh_x(2)-mesh_x(1));
xd=0:h/2:1;
plot(xd,u);
hold on
plot(xd,0.5*xd.*(1-xd)+1,'dr')