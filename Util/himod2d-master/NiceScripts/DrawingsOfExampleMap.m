addpath ../MapHandler/
close all
haty=sym('y/(0.1*cos(2*x)+1)');
[Map]=provideMap(haty);
x=linspace(30,40);
xref=linspace(0,10);
[X,Y]=meshgrid(x,[1,-1]);
y=Map.invhaty(X,Y);
LW=2;
figure 
hold on 
plot(x,y(1,:),'k','linewidth',LW)
plot(x,y(2,:),'k','linewidth',LW)
plot([x(end),x(end)],y(:,end),'k','linewidth',LW)
plot([x(1),x(1)],y(:,1),'k','linewidth',LW)
plot([x(1),x(end)],[0,0],'--k','linewidth',LW)
plot(xref,ones(size(x)),'k','linewidth',LW)
plot(xref,-ones(size(x)),'k','linewidth',LW)
plot([xref(end),xref(end)],[-1,1],'k','linewidth',LW)
plot([xref(1),xref(1)],[-1,1],'k','linewidth',LW)
plot([xref(1),xref(end)],[0,0],'--k','linewidth',LW)

xlim([-1,50])
ylim([-2,2])


set(gcf,'Units','normalized')
set(gca,'Units','normalized')

ax = axis;
ap = get(gca,'Position');

% xgamma=x(end)-2;
% ygamma=1.5;
% xo = [xgamma , x(1)+5.6];
% yo = [ygamma , 1.1];
% x1 = [xgamma , x(1)+6.7];
% y1 = [ygamma , -0.85];
% xp = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
% yp = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
% xp1 = (x1-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
% yp1 = (y1-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
% annotation('arrow',xp,yp,'linewidth',LW);
% annotation('arrow',xp1,yp1,'linewidth',LW);
% text(xgamma+0.1,ygamma+0.1,'\Gamma_{Lat}','FontSize',30)

% xo = [x(1) , x(end)+1];
% yo = [0 , 0];
% x1 = [x(1) , x(1)];
% y1 = [0 , 1.5];
% xp = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
% yp = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
% xp1 = (x1-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
% yp1 = (y1-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
% annotation('arrow',xp,yp,'linewidth',LW);
% annotation('arrow',xp1,yp1,'linewidth',LW);
% text(x(end)+1 +0.1 , 0  + 0.1,'x','FontSize',30)
% text(x(1)+0.1   ,1.5 + 0.1,'y','FontSize',30)

text(x(end)-2,-0.7,'\Omega','FontSize',30)
text(xref(end)-2,-0.7,'\Omega','FontSize',30)


ar=linspace(xref(end)+1,x(1)-1);
plot(ar,0.005*(ar-ar(1)).*(ar(end)-ar),'linewidth',LW,'color','k')
text(ar(floor(end/2)),0,'\psi','FontSize',30)
plot([ar(end),ar(end)-0.5],[0,0.2],'linewidth',LW,'color','k')
plot([ar(end),ar(end)-1.5],[0,0],'linewidth',LW,'color','k')