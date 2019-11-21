cinf=-0.01;
csup=0.1;
addpath ../../../../Core
addpath ../../../../Util
addpath ../../../../MapHandler
clines=linspace(cinf,csup,20);
[p,~,~]=importfilemesh('data/referencedomain.msh');
nx=601
ny=size(p,2)/nx
x=p(1,1 :nx);
y=p(2,1:nx:ny*nx);
dt=0.01;
[X,Yhat]=meshgrid(x,y);
umax=0.4;
tau=0.5;
haty=sym('(y-umax*sin(t/tau)/2)/(umax*sin(t/tau)/2+1)');
haty=subs(subs(haty,'umax',umax),'tau',tau);
[Map]=provideMapUnsteady(haty);
kmax=49;
t=0;
orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
nomefilmato=strcat(pwd,'/ffem',ora,'.avi');
filmatoavi = avifile(nomefilmato,'fps',5);
for k=0:kmax
    t=t+dt;
    %mesh=['mesh',num2str(k),'.msh'];
    soln=['movingffem',num2str(k),'.dat'];
    %read mesh
    %read solution
    [sol]=importfiledata(soln);
    sols=zeros(ny,nx);
    for i=1:ny
        sols(i,:)=sol((i-1)*nx+1:i*nx);
    end
    Y=Map.invhaty(t,X,Yhat);
    %plot
    f(k+1)=figure(k+1);
    warning('OFF','MATLAB:contourf:ConstantData');
    contourf(X,Y,sols,clines)
    warning('ON','MATLAB:contourf:ConstantData');
    %pdeplot(p,[],t,'xydata',sol,'colormap','jet','contour','on','levels',clines);%,'zdata',sol
    grid on
    axis equal
    xlim([0,6])
    ylim([-2,2])
    caxis([cinf,csup]);
    colorbar;
    title(['t=',num2str(t)])
    
    set(f(k+1),'visible','off')
    if(k>1)
        close(f(k-1));
    end
    F = export_fig(f(k+1), '-nocrop', '-a1','-transparent');
    filmatoavi=addframe(filmatoavi,F);
    %view(0,90)
    %zlim([cinf,csup])
    display(['status: ' num2str(k/kmax*100),'/%']);
end
filmatoavi=close(filmatoavi);