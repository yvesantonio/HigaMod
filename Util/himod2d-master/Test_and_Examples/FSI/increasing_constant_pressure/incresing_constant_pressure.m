clc
close all
clear all

addpath ../../../Stokes
addpath ../../../Util
addpath ../../../MapHandler
addpath ../../../FSI

Cx=struct('cinf',-0.2,'csup',0.2);
Cy=struct('cinf',-0.2,'csup',+0.2);
Cp=struct('cinf',-5,'csup',+5);

Contour         = struct('out',0,'video',0,'visible','off','exp',[1,1.4,1.5,1.6,2,2.3]);
Profile         = struct('out',0,'video',0,'visible','off','exp',[],'x',0.5);
ProfileVelocity = struct('out',0,'video',0,'visible','off','exp',[1,1.4,1.5,1.6,2,2.3],'x',0.5);
Pressure        = struct('out',0,'video',0,'visible','off','exp',[]);
Domain          = struct('out',1,'video',0,'visible','on','exp',[1 2.5 5 10 20]);

out=struct(...
    'Contour',Contour,'Profile',Profile,'Pressure',Pressure,'Domain',Domain,'ProfileVelocity',ProfileVelocity,...
    'Cx',Cx,'Cy',Cy,'Cp',Cp, ...
    'ended',false);

I=[0,5];
nb_elem=15;
mesh=struct('I',I,'nb_elem',nb_elem);

toll=1e-6;
m=7;

N=4;
dt=1/N;
ntimestep=20*N;

nu=1;
K=10;

fx=@(x,y,t) 0+0*x+0*y+0*t;
fy=@(x,y,t) 0+0*x+0*y+0*t;

% Inflow data
essential_inflow=false;
Pinflow=@(t) 0.5*t*(t<5)+2.5*(t>=5);
Poutflow=Pinflow;


%Boundary conditions.
INF_x='neu';
OUT_x='neu';

INF_y='neu';
OUT_y='neu';

BC=struct('INF_x',INF_x,'dataINFx',Pinflow,...
          'INF_y',INF_y,'dataINFy',zeros(m,1),...
          'OUT_x',OUT_x,'dataOUTx',Poutflow,...
          'OUT_y',OUT_y,'dataOUTy',zeros(m,1));
      
Data=struct('nu',nu,'m',m,'dt',dt,'ntimestep',ntimestep, ...
'forcex',fx,'forcey',fy,'BC',BC,'toll',toll,'K',K);

nb_quadnode_y=32; %quadratura su tutto l'intervallo -1,1
nb_quadnode_x=4;  %quadratura sul singolo elemento finito
quad=struct('nx',nb_quadnode_x,'ny',nb_quadnode_y);

sol=FSISolver(mesh,quad,Data,out);