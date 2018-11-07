addpath ../
%% Caso Dir-Rob

x=linspace(0,1,200);
m=18;

global coeffrobin mu 
mu=@(x,y) 1;
coeffrobin = 2;

[mb,mb_y]=new_modal_basis(m,x,'dir','rob');

% figure(1);
% plot(x,mb)
% title('Basi')
% legend('0','1s','1c','2s','2c');
% grid on
% 
% figure(2);
% plot(x,mb_y)
% legend('0_x','1s_x','1c_x','2s_x','2c_x')
% title('Derivate')
% grid on

for i=1:m
	for k=1:m
			aus = mb(:,i).*mb(:,k);
			delta(i,k) = integrami(aus,x);				
	end
end


R_0 = @(x,y) -mu(0,0)*x+coeffrobin*y;

display('Ortonormale ? ')
delta
display('Rispetta Robin omogeneo?')
R_0(mb_y(1,:),mb(1,:))
display('Dirichlet omogeneo?')
mb(end,:)


%% Caso Rob-Rob

x=linspace(0,1,200);
m=18;

global coeffrobin mu 
mu=@(x,y) 1;
coeffrobin = 2;

[mb,mb_y]=new_modal_basis(m,x,'rob','rob');

% figure(1);
% plot(x,mb)
% title('Basi')
% legend('0','1s','1c','2s','2c');
% grid on
% 
% figure(2);
% plot(x,mb_y)
% legend('0_x','1s_x','1c_x','2s_x','2c_x')
% title('Derivate')
% grid on

for i=1:m
	for k=1:m
			aus = mb(:,i).*mb(:,k);
			delta(i,k) = integrami(aus,x);				
	end
end

R_0 = @(x,y) -mu(0,0)*x+coeffrobin*y;
R_1 = @(x,y) +mu(0,0)*x+coeffrobin*y;

display('Ortonormale ? ')
    delta
display('Rispetta Robin omogeneo in 0?')    
    R_0(mb_y(1,:),mb(1,:))
display('e in 1?')
    R_1(mb_y(end,:),mb(end,:))