clc
close all

addpath ../../../Stokes
addpath ../../../Util
addpath ../../../MapHandler

Cx=struct('cinf',-0.1,'csup',1.2);
Cy=struct('cinf',-0.1,'csup',+0.5);
Cp=struct('cinf',-0.1,'csup',+12);

Contour         = struct('out',1,'video',0,'visible','on','exp',[1,1.4,1.5,1.6,2,2.3]);
Profile         = struct('out',0,'video',0,'visible','off','exp',[],'x',0.5);
ProfileVelocity = struct('out',0,'video',0,'visible','off','exp',[1,1.4,1.5,1.6,2,2.3],'x',0.5);
Pressure        = struct('out',0,'video',0,'visible','off','exp',[]);
Domain          = struct('out',0,'video',0,'visible','off','exp',[]);

out=struct(...
    'Contour',Contour,'Profile',Profile,'Pressure',Pressure,'Domain',Domain,'ProfileVelocity',ProfileVelocity,...
    'Cx',Cx,'Cy',Cy,'Cp',Cp, ...
    'ended',false);

I=[0,5];
nb_elem=10;
mesh=struct('I',I,'nb_elem',nb_elem);

m=8;
N=30;
dt=1/N;
ntimestep=N*1;
nu=1;
fx=@(x,y,t) 0+0*x+0*y+0*t;
fy=@(x,y,t) 0+0*x+0*y+0*t;
bx=@(x,y,t) 0+0*x.*(1-x)+0*y+0*t;

Pinflow=@(t) 10 + 0*t;
%Boundary conditions.
ZERO=@(t) 0*t;
INF_x='neu';
INF_y='dir';
OUT_x='neu';
OUT_y='dir';
BC=struct('INF_x',INF_x,'dataINFx',Pinflow,...
          'INF_y',INF_y,'dataINFy',zeros(m,1),...
          'OUT_x',OUT_x,'dataOUTx',ZERO,...
          'OUT_y',OUT_y,'dataOUTy',zeros(m,1));
%nel codice ?? imposta u_y = 0 all'inflow
Data=struct('nu',nu,'m',m,'dt',dt,'ntimestep',ntimestep, ...
'forcex',fx,'forcey',fy,'betax',bx,'BC',BC);
nb_quadnode_y=32; %quadratura su tutto l'intervallo -1,1
nb_quadnode_x=4;  %quadratura sul singolo elemento finito
quad=struct('nx',nb_quadnode_x,'ny',nb_quadnode_y);

umax=0.4;
tau=1/2/pi;
haty=sym('(y-umax*sin(t/tau)/2)/(umax*sin(t/tau)/2+1)');
haty=subs(subs(haty,'umax',umax),'tau',tau);
%haty=sym('y');
[Map]=provideMapUnsteady(haty);
sol=StokesSolver(mesh,quad,Data,out,Map,false);