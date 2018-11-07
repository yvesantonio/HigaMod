cinf=-0.01;
csup=0.1;
addpath ../../../../Core
addpath ../../../../Util
addpath ../../../../MapHandler
clines=linspace(cinf,csup,20);
[p,~,~]=importfilemesh('referencedomain.msh');
nx=251
ny=size(p,2)/nx
x=p(1,1:nx);
y=p(2,1:nx:ny*nx);
dt=0.01;
[X,Yhat]=meshgrid(x,y);

haty=sym('y/(0.5*sin(50*t)*cos(pi/3*x+pi/2)+1)');
[Map]=provideMapUnsteady(haty);
kmax=99;
t=0;
nomefilmato=strcat(pwd,'/ffem_ref_h=',num2str(5/nx),'.avi');
filmatoavi = avifile(nomefilmato,'fps',5);
for k=0:1:kmax
    t=t+dt;
    %mesh=['mesh',num2str(k),'.msh'];
    soln=['wavy_unsteady_mappedback_ffem',num2str(k),'.dat'];
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
    set(f(k+1),'visible','off')
    warning('OFF','MATLAB:contourf:ConstantData');
    contourf(X,Y,sols,clines)
    warning('ON','MATLAB:contourf:ConstantData');
    %pdeplot(p,[],t,'xydata',sol,'colormap','jet','contour','on','levels',clines);%,'zdata',sol
    grid on
    axis equal
    xlim([0,5])
    ylim([-2,2])
    caxis([cinf,csup]);
    colorbar;
    title(['t=',num2str(t)])
    
    if(k>1)
        close(f(k-1));
    end
    if(k==59||k==66||k==93)
         export_fig(f(k+1), 'filename',['t=',num2str(t),'ffem_ref_h',num2str(5/nx)],'-png','-transparent')
    end
    F = export_fig(f(k+1), '-nocrop', '-a1','-transparent');
    filmatoavi=addframe(filmatoavi,F);
    %view(0,90)
    %zlim([cinf,csup])
    display(['status: ' num2str(k/kmax*100),'/%']);
end
filmatoavi=close(filmatoavi);