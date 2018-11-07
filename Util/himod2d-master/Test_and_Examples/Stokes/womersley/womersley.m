clc
close all
clear all

addpath ../../../Stokes
addpath ../../../Util
addpath ../../../MapHandler

Cx=struct('cinf',-0.01,'csup',0.03);
Cy=struct('cinf',-0.1,'csup',+0.1);
Cp=struct('cinf',-0.3,'csup',+0.3);

Contour         = struct('out',0,'video',0,'visible','on','exp',[1,1.4,1.5,1.6,2,2.3]);
Profile         = struct('out',0,'video',0,'visible','off','exp',[],'x',0.5);
ProfileVelocity = struct('out',1,'video',0,'visible','off','exp',[1,1.4,1.5,1.6,2,2.3],'x',0.5);
Pressure        = struct('out',0,'video',0,'visible','off','exp',[]);
Domain          = struct('out',0,'video',0,'visible','off','exp',[]);

out=struct(...
    'Contour',Contour,'Profile',Profile,'Pressure',Pressure,'Domain',Domain,'ProfileVelocity',ProfileVelocity,...
    'Cx',Cx,'Cy',Cy,'Cp',Cp, ...
    'ended',false);
%%%%%%%%%%%%
I=[0,1];
nb_elem=6;
mesh=struct('I',I,'nb_elem',nb_elem);

m=3;
N=200;
dt=1/N
ntimestep=N*2.5;
global nu w
nu=1/10;
fx=@(x,y,t) 0+0*x+0*y+0*t;
fy=@(x,y,t) 0+0*x+0*y+0*t;
bx=@(x,y,t) 0+0*x.*(1-x)+0*y+0*t;
% Inflow data
essential_inflow=false;
w=pi;
W=sqrt(w/nu)
Pinflow=@(t) 2*nu*sin(w*t);

%Boundary conditions.
INF_x='neu';
INF_y='dir';
OUT_x='neu';
OUT_y='dir';
ZERO=@(t) 0*t;
BC=struct('INF_x',INF_x,'dataINFx',Pinflow,...
          'INF_y',INF_y,'dataINFy',zeros(m,1),...
          'OUT_x',OUT_x,'dataOUTx',ZERO,...
          'OUT_y',OUT_y,'dataOUTy',zeros(m,1));
%nel codice è imposta u_y = 0 all'inflow
Data=struct('nu',nu,'m',m,'dt',dt,'ntimestep',ntimestep, ...
'forcex',fx,'forcey',fy,'betax',bx,'BC',BC);

nb_quadnode_y=32; %quadratura su tutto l'intervallo -1,1
nb_quadnode_x=4;  %quadratura sul singolo elemento finito
quad=struct('nx',nb_quadnode_x,'ny',nb_quadnode_y);

haty=sym('2*y-1');
[Map]=provideMapUnsteady(haty);
sol=StokesSolver(mesh,quad,Data,out,Map,false);