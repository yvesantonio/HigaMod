addpath ../Core
addpath ../Util

close all
clear all
clc

numberofnodes_y=64;
[trash, node_y, weights_y] = gausslegendre( numberofnodes_y );
m=4;
projected=zeros(length(node_y),1);
% alpha=.01;
% xl=0.4-alpha;
% xr=0.6+alpha;
% if(alpha~=0)
%     f=@(x) (x>0.4).*(x<0.6) ...
%     +1/alpha*(x-xl).*(x<=0.4).*(x>=xl) ...
% + 1/alpha*(xr-x).*(x<=xr).*(x>=0.6);
% else
%     f=@(x) (x>0.4).*(x<0.6);
% end
f=@(y) y;
feval=f(node_y);
C=struct('mu',@(x,y) 1+0*x+0*y,'coeffrobin',1);
% basis=new_modal_basis(m, node_y, 'dir','dir',C);
basis=new_modal_basis_FOURIER(m,node_y);
fcoeff=projection(basis,feval,weights_y);
for i=1:m
    projected=projected+basis(:,i)*fcoeff(i);
end

x=linspace(0,1,1e6);
% plot(x,f(x),'k','linewidth',3);
% hold on
% plot(node_y,projected,'r--','linewidth',3);
% legend('\fontsize{20} f(x)','Projection m=8')
% set(gca,'fontsize',20)
% grid on

createfigure(x,f(x),node_y,projected)
