clc
clear all
close all

% theta=pi/4;
% haty=sym('2*(y+1-x*sin(t)*tantheta)/(2-x*sin(t)*tantheta)-1');
% haty=subs(haty,'tantheta',tan(theta));

%haty=sym('y/(0.1*cos(12*x-pi*t)+1)');
%invhaty=sym('y+y*cos(x*6.283185307179586+1.570796326794897)*sin(t*2.0e1)*5.0e-1');
haty=sym('y/(0.5*sin(1e4*t)*cos(pi/3*6*x+pi/2)+1)');
% haty=finverse(invhaty,'y');
[Map]=provideMapUnsteady(haty);
xhat=linspace(0,1);
yhat=linspace(-1,1);
[Xhat,Yhat]=meshgrid(xhat,yhat);
X=Xhat;
dt=0.01;
for t=0:dt:10
    Y=Map.invhaty(t,Xhat,Yhat);
    
    fhat=@(t,xhat,yhat) 0*t+(xhat.^2-xhat)*(yhat.^2-1);
    fhatev=fhat(t,Xhat',Yhat');
    f=@(t,x,y) (x.^2-x)*(Map.haty(t,x,y).^2-1)+0*t;
    fev=f(t,X',Y');
    
    subplot(3,1,1);
    contourf(X,Y,fev)
    xlim([0,1])
    ylim([-2,2])
    title(['t = ',num2str(t)]);
    subplot(3,1,2);
    contourf(Xhat,Yhat,fhatev)
     xlim([0,1])
    ylim([-2,2])
    subplot(3,1,3);
    contourf(X,Y,fhatev)
     xlim([0,1])
    ylim([-2,2])
    pause
end