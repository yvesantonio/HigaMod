close all
clear all
clc
x=linspace(0,0.9);

f=@(x) sqrt(1-x.^2);
plot(x,f(x) ,'linewidth',2,'color','k')
xlim([-2,2])
xa=0.7;
xb=1.1;
hold on
r=@(x) f(xa)/xa*x;
plot([xa,xb],[f(xa),r(xb)],'linewidth',2,'color','k')
plot([xb,xb-0.2],[r(xb),r(xb)-0.1],'linewidth',2,'color','k')
plot([xb,xb-0.1],[r(xb),r(xb)-0.2],'linewidth',2,'color','k')
plot([xa,xb],[f(xa),f(xa)],'--k','linewidth',2)
text(xa+0.4,f(xa)+0.1,'\theta','FontSize',30)
r=0.2;
fa=@(x) f(xa) + sqrt(r^2-(x-xa).^2);
x=linspace(xa+0.15,xa+0.2);
plot(x,fa(x),'linewidth',2,'color','k')
text(xa+0.55,f(xa)+0.5,'n(x)','FontSize',30)
axis equal