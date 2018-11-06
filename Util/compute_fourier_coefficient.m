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
function [uk]=compute_fourier_coefficient(u,Coeff_forma,bc_up,bc_down)
%
% function [uk]=compute_fourier_coefficient(u,Coeff_forma,bc_up,bc_down)
% 
% teo.aletti@gmail.com
% Calcola la norma L2 dei coefficienti di fourier della funzione u, attenzione 
% controllare la compatibilità con i dati al bordo, sarebbe necessario aggiungere il rilevamento
% 
N=32;  % numero dei coefficienti di Fourier calcolato
n=160; % numero dei nodi di quadratura

uk=zeros(N,1);

[~,x,w] = gausslegendre(n);
y=x;

[mb,~]=new_modal_basis(N,x,bc_up,bc_down,Coeff_forma);

[x,y]=meshgrid(x,y);
u_valued=u(x,y);% le righe scorrono la coordinata x (colonne scorrono la y)
for i=1:N
	%fisso la x (cioè fisso una riga di u_valued)
	% e integro 
	coeff_pwise=integrate(u_valued.*(ones(n,1)*mb(:,i)'),w);%coeff_pwise rappresenta per ogni x fissato il coefficiente di fourier
	uk(i)=sqrt(integrate(coeff_pwise'.^2,w));   % norma L2
end
