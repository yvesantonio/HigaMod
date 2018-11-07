%**********************************************%
% Coefficienti della forma bilineare           %
%**********************************************%
chi=1;
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) ( 20 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  2 + 0*x+0*y ); % Trasporto verticale
sigma = @(x,y) (  2 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);

addpath ../Core
x=linspace(0,1,2000);
m=5;

[mb,mb_y]=new_modal_basis(m,x,'dir','dir',Coeff_forma,2);

set(gca,...
    'FontUnits','points', ...
    'FontSize',20)
figure;
plot(x,mb,'linewidth',2)
xlabel('y');
ylabel('\phi')
title('Basi')
legend('0','1','2','3','4','Location','EastOutside');
%save_img('eps','basidirrob3')
grid on

figure;
plot(x,mb_y)
legend('0_x','1_x','2_x','3_x','4_x')
title('Derivate')
grid on

x=linspace(0,2,2000);
delta=zeros(m,m);
for i=1:m
	for k=1:m
			aus = mb(:,i).*mb(:,k);
			delta(i,k) = integrami(aus,x);				
	end
end
delta
