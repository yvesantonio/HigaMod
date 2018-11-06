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
function [L2, H1] = error_himod_P2(size_mb,omega_xl,omega_xr,hx, u,a_ril,b_ril,bc_up, bc_down, Coeff_forma,true_sol,true_solx,true_soly,dforma,fisx,fisy,caso,p,k,psJac2)
%
% function [L2, H1] = error_himod(size_mb,omega_xl,omega_xr,hx, u,a_ril,b_ril,bc_up, bc_down, Coeff_forma,true_sol,true_solx,true_soly)
% 
% Calcola l'errore, nelle due norme, di una soluzione ottenuta su un singolo dominio con il metodo Hi-Mod basi istruite.
%
%
% Riassunto dei parametri di input della function:
% - size_mb       = dimensione base modale;
% - omega_xl      = estremo sinistro dell'intervallo di integrazione;
% - omega_xr      = estremo destro dell'intervallo di integrazione;
% - hx            = passo di discretizzazione dell'intervallo;
% - u             = soluzione numerica del problema;
% - a_ril         = rilevamento;
% - b_ril         = rilevamento;
% - bc_up         = condizioni al contorno sopra (tipo 'dir'...);
% - bc_down       = condizioni al contorno sotto (tipo 'dir'...);
% - Coeff_forma   = coefficienti della forma bilineare;
% - true_sol      = soluzione analitica;
% - true_solx     = derivata dx analitica;
% - true_soly     = derivata dy analitica;
%


%% ATTENZIONE AI SQRT(2) CHE SONO DA SISTEMARE CHE SONO
%  MESSI APPOSTA PER IL CASO PARTICOLARE Y=X.
%%
  nqnx = 80;
  nqny = 80;
    
% Numero di sezioni verticali: ne intervalli e nx estremi (nx = ne+1)
  ne = round((omega_xr-omega_xl)/hx);
  nx = ne+1; %nodi in x
  nxIGA = 2*(ne*k+p+1-k);     % numero nodi con IGA
  hx=hx/p;              % distanza si dimezza tra un nodo e un altro
  
% Mesh grossolana in x, nodi equispaziati
  meshx = zeros(nx,1);
  meshx2 = zeros(nxIGA,1);
  meshx(1) = omega_xl;
  meshx2(1) = omega_xl;

  %Servono le varie L dei casi?
  
  for i=2:nx
      meshx(i) = meshx(i-1)+hx*2*sqrt(2);
  end
  % meshx è praticamente inutile
  psJac = psJac2;
  
  L = sum(psJac);
  disp('Lunghezza del dominio [error]:')
  disp(L)
  
  for i=2:nxIGA
     %meshx2(i) = meshx2(i-1)+hx*dforma(meshx2(i-1));   %MODIFICATI PER IL CASO 
     meshx2(i) = meshx2(i-1)+psJac(i-1);
  end;

  

% Nodi di intergrazione alla Gass-Legendre
  [nqnx,xq,wxq] = gausslegendre(nqnx);

% I nodi devono essere riscalati
%xq = L*xq;
%OK ?????????????

% - xq           = valore in x dei nodi di quadratura;
% - wxq          = peso dei vari nodi di quadratura;
 
% Collocazione dei nodi di quadratura nelle sezioni verticali
  [nqny,mesh_y,mesh_wy] = gausslegendre(nqny);
%mesh_wy = mesh_wy*sqrt(2);
% Mesh fine in x, quella dei soli nodi di quadratura
 

[mesh_xx, mesh_wx] = quadrature_rule(omega_xl, omega_xr*L, xq, wxq); 

%mesh_xx' è un vettore che contiene gli 80 punti di quadratura tra 0 e L
% Adesso creo un vettore di lunghezza nqnx*size_mb dove salvo quanto vale
% il polinomio polyfit per ogni mb nei nqnx punti
approx_quad = zeros(nqnx*size_mb);

for m = 1 : size_mb
    % Estraggo dal vettore della soluzione i 2ncp punti del modo m-simo
    interpol = u((m-1)*(nxIGA)+1:m*nxIGA);
    % Interpolo i dati della soluzione per il modo m-simo
    xxx =linspace(0,1*L,nxIGA)';
    coef = polyfit(xxx,interpol,7);
    % Defiizione del polinomio di interpolazione
    interpoly = @(x) coef(1)*x.^7 + coef(2)*x.^6 + coef(3)*x.^5 + coef(4)*x.^4+...
                     coef(5)*x.^3 + coef(6)*x.^2 + coef(7)*x + coef(8);
                 
   % interpoly = @(x) coef(1)*x.^10 + coef(2)*x.^9 + coef(3)*x.^8 + coef(4)*x.^7+...
   %                  coef(5)*x.^6  + coef(6)*x.^5 + coef(7)*x.^4 + coef(8)*x.^3+...
   %                  coef(9)*x.^2  + coef(10)*x.^1 + coef(11); 
    disp('***********************************')
    disp('Coefficienti della interpolazione:')
    disp(coef')
    % Valutazione nei nqnx punti
    intervalue = interpoly(mesh_xx);
    %plot(mesh_xx,intervalue)
    %hold on
    %ll = 0.5*(sqrt(5)+1/2*log(2+sqrt(5)));
    %dipx  = @(x) (1/3*x.^3-ll^2*x);
    %figure;
    %plot(mesh_xx,dipx(mesh_xx),'r')
    % Salvataggio nel vettore
    approx_quad((m-1)*nqnx+1:m*nqnx) = intervalue;
end
  
disp('Punto di controllo semifinale');

% Visualizzazione
sol  = zeros(nqnx,nqny);
solx = zeros(nqnx,nqny);
soly = zeros(nqnx,nqny);
  
  
[mod_b, mod_by] = new_modal_basis(size_mb, mesh_y, bc_up, bc_down, Coeff_forma);

  
 for i=1:nqnx
    for j=1:nqny
        for imb = 1 : size_mb
            solx(i,j) = solx(i,j) + approx_quad((imb-1)*nqnx+i)*mod_b(j,imb);
            % Per questo devo valutare la derivata MMH
            soly(i,j) = soly(i,j) + approx_quad((imb-1)*nqnx+i)*mod_by(j,imb);
            sol(i,j)  = sol(i,j)  + approx_quad((imb-1)*nqnx+i)*mod_b(j,imb);
        end
    sol(i,j)  = sol(i,j)  + a_ril*mesh_y(j)+ b_ril;
    soly(i,j) = soly(i,j) + a_ril;
    end
 end

 [X,Y] = meshgrid(mesh_xx, mesh_y-0.5);
 Tsol = true_sol(X,Y)';
 Tsolx= true_solx(X,Y)';
 Tsoly= true_soly(X,Y)';

 L2 = sqrt( (((sol-Tsol).^2)' * mesh_wx)' * mesh_wy);
 H1 = sqrt( (((sol-Tsol).^2 + (solx-Tsolx).^2 + (soly-Tsoly).^2)' * mesh_wx)' * mesh_wy);
 
 disp('Punto controllo finale')
 % Ho cancellato un *2 che c'era dentro mmm
 
 figure;
 subplot(1,3,1); surf(X,Y,sol');
 subplot(1,3,2); surf(X,Y,Tsol');
 subplot(1,3,3); surf(X,Y,abs(sol'-Tsol'));
 
  
% figure;
% plot(X(1,:),Tsolx(1,:)')
% hold on
% plot(X(1,:),solx(1,:)','r')
 %pause;
end
