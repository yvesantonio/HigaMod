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
%-+------------------------------------------------------------------------

function [L2, H1] = error_himod(size_mb,omega_xl,omega_xr,hx, u,a_ril,b_ril,bc_up, bc_down, Coeff_forma,true_sol,true_solx,true_soly)
%
% function [L2, H1] = error_himod(size_mb,omega_xl,omega_xr,hx, u,a_ril,b_ril,bc_up, bc_down, Coeff_forma,true_sol,true_solx,true_soly)
% 
% Calcola l'errore, nelle due norme, di una soluzione ottenuta su un singolo dominio con il metodo Hi-Mod basi istruite.

  nqnx = 80;
  nqny = 80;
  
% Numero di sezioni verticali: ne intervalli e nx estremi (nx = ne+1)
  ne = round((omega_xr-omega_xl)/hx);
  nx = ne+1;%nodi in x

% Mesh grossolana in x, nodi equispaziati
  meshx = zeros(nx,1);
  meshx(1) = omega_xl;
  for i=2:nx
     meshx(i) = meshx(i-1)+hx; 
  end;

% Nodi di intergrazione alla Gass-Legendre
  [nqnx,xq,wxq] = gausslegendre(nqnx);

% Collocazione dei nodi di quadratura nelle sezioni verticali
  [nqny,mesh_y,mesh_wy] = gausslegendre(nqny);

% Mesh fine in x, quella dei soli nodi di quadratura
  [mesh_xx, mesh_wx] = quadrature_rule(omega_xl, omega_xr, xq, wxq); 
  
  index = zeros(nqnx,2);
  j = 1;
  
  for i = 1:nqnx%ad ogni nodo associo l'indice j del nodo della mesh fem che si trova a sx
     while( mesh_xx(i) > meshx(j+1) && j < ne)
       j = j+1;
     end
     index(i,1) = j;% nodo a sx nodo a dx
     index(i,2) = j+1;% nodo che si trova a destra
  end
  
  % Visualizzazione
  sol=zeros(nqnx,nqny);
  solx=zeros(nqnx,nqny);
  soly=zeros(nqnx,nqny);

  [mod_b, mod_by] = new_modal_basis(size_mb, mesh_y, bc_up, bc_down, Coeff_forma);

 for i=1:nqnx
    for j=1:nqny
      for imb=1:size_mb
         sol(i,j)=sol(i,j)+(  u(index(i,1)+(imb-1)*nx)*(meshx(index(i,2))-mesh_xx(i))/hx ...
                            + u(index(i,2)+(imb-1)*nx)*(mesh_xx(i)-meshx(index(i,1)))/hx )*mod_b(j,imb);
         solx(i,j)=solx(i,j)+(  (u(index(i,2)+(imb-1)*nx) - u(index(i,1)+(imb-1)*nx))/hx )*mod_b(j,imb);
         soly(i,j)=soly(i,j)+(  u(index(i,1)+(imb-1)*nx)*(meshx(index(i,2))-mesh_xx(i))/(meshx(index(i,2))-meshx(index(i,1))) ...
                            +   u(index(i,2)+(imb-1)*nx)*(mesh_xx(i)-meshx(index(i,1)))/(meshx(index(i,2))-meshx(index(i,1))) )*mod_by(j,imb);
      end
	sol(i,j)  = sol(i,j)  + a_ril*mesh_y(j) + b_ril;
    soly(i,j) = soly(i,j) + a_ril;
    end
    
 end

 [X,Y] = meshgrid(mesh_xx, mesh_y);
 %true_sol=matlabFunction(true_sol);
 %true_solx=matlabFunction(true_solx);
 %true_soly=matlabFunction(true_soly);
 
 Tsol = true_sol(X,Y)';
 Tsolx= true_solx(X,Y)';
 Tsoly= true_soly(X,Y)';
 
 L2 = sqrt( (((sol-Tsol).^2)' * mesh_wx)' * mesh_wy);
 H1 = sqrt( (((sol-Tsol).^2 + (solx-Tsolx).^2 + (soly-Tsoly).^2)' * mesh_wx)' * mesh_wy);

 
 %subplot(1,3,1); surf(X,Y,sol');
 %subplot(1,3,2); surf(X,Y,Tsol');
 %subplot(1,3,3); surf(X,Y,abs(sol'-Tsol'));
 
 %pause;
end
