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
