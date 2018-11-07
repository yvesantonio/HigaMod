clc
close all
clear all

addpath ../../../Stokes
addpath ../../../Util
addpath ../../../MapHandler/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing setting %
Cx=struct('cinf',-0.05,'csup',0.05);
Cy=struct('cinf',-0.01,'csup',+0.01);
Cp=struct('cinf',-0.4,'csup',0.4);

Contour         = struct('out',1,'video',0,'visible','off' ,'exp',[2]);
Profile         = struct('out',0,'video',0,'visible','off','exp',[100],'x',2);
ProfileVelocity = struct('out',0,'video',0,'visible','off','exp',[1,1.4,1.5,1.6,2,2.3],'x',4);
Pressure        = struct('out',0,'video',0,'visible','off','exp',[100]);
Domain          = struct('out',0,'video',0,'visible','off','exp',[100]);

out=struct(...
    'Contour',Contour,'Profile',Profile,'Pressure',Pressure,'Domain',Domain,'ProfileVelocity',ProfileVelocity,...
    'Cx',Cx,'Cy',Cy,'Cp',Cp, ...
    'ended',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadrature data
nb_quadnode_y=32; %quadratura su tutto l'intervallo -1,1
nb_quadnode_x=4;  %quadratura sul singolo elemento finito
quad=struct('nx',nb_quadnode_x,'ny',nb_quadnode_y);

% Map generation
haty=sym('(y-0.4*exp(-(x-3)^2/0.2))/(1-0.4*exp(-(x-3)^2/0.2))');
Map=provideMapUnsteady(haty);

% Omega 1D
I=[0,6];
number_elem=60;
mesh=struct('I',I,'nb_elem',number_elem);

m=16;

N=20;
dt=1/N;
ntimestep=N*2.1;

nu=0.01;
fx=@(x,y,t) 0+0*x+0*y+0*t;
fy=@(x,y,t) 0+0*x+0*y+0*t;
bx=@(x,y,t) 0+0*x.*(1-x)+0*y+0*t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Boundary conditions.
%dirdata=zeros(m,1);
w=pi;
Px=@(t) 2*nu*sin(w*t);
ZERO=@(t) 0*t;
INF_x='neu';
INF_y='dir';
OUT_x='neu';
OUT_y='dir';
BC=struct('INF_x',INF_x,'dataINFx',Px,...
          'INF_y',INF_y,'dataINFy',zeros(m,1),...
          'OUT_x',OUT_x,'dataOUTx',ZERO,...
          'OUT_y',OUT_y,'dataOUTy',zeros(m,1));
      
Data=struct('nu',nu,'m',m,'dt',dt,'ntimestep',ntimestep, ...
'forcex',fx,'forcey',fy,'betax',bx,'BC',BC,'toll',1e-4);

% sol=StokesSolver(mesh,quad,Data,out,Map,false);
 sol=FixedpointSolver(mesh,quad,Data,out,Map,false);
 load handel
 sound(y(1:2e4),Fs)