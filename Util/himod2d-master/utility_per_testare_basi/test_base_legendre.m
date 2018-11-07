addpath ../Core
x=linspace(-1,1,2000);
m=6;

[mb,mb_y]=new_modal_basis_legendre(m,x,'rob','rob');

set(gca,...
    'FontUnits','points', ...
    'FontSize',20)
figure(1);
plot(x,mb,'linewidth',2)
xlabel('y');
ylabel('\phi')
title('Basi')
legend('1','2','3','4','5','6','Location','EastOutside');
%save_img('eps','basidirrob3')
grid on

figure(2);
plot(x,mb_y)
legend('0_x','1_x','2_x','3_x','4_x')
title('Derivate')
grid on

delta=zeros(m,m);
for i=1:m
	for k=1:m
			aus = mb(:,i).*mb(:,k);
			delta(i,k) = integrami(aus,x);				
	end
end
delta