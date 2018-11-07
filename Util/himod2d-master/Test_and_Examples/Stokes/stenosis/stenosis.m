clc
close all
clear all

addpath ../../../Stokes
addpath ../../../Util
addpath ../../../MapHandler/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing setting %
Cx=struct('cinf',-0.05,'csup',0.5);
Cy=struct('cinf',-0.5,'csup',+0.5);
Cp=struct('cinf',-0.1,'csup',5.1);

Contour         = struct('out',1,'video',0,'visible','on' ,'exp',[1]);
Profile         = struct('out',1,'video',0,'visible','off','exp',[100],'x',2);
ProfileVelocity = struct('out',1,'video',0,'visible','off','exp',[100],'x',2);
Pressure        = struct('out',1,'video',0,'visible','off','exp',[100]);
Domain          = struct('out',1,'video',0,'visible','off','exp',[100]);

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
number_elem=120;
mesh=struct('I',I,'nb_elem',number_elem);

m=7;

dt=1;
ntimestep=100;

nu=1;
fx=@(x,y,t) 0+0*x+0*y+0*t;
fy=@(x,y,t) 0+0*x+0*y+0*t;
bx=@(x,y,t) 0+0*x.*(1-x)+0*y+0*t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Boundary conditions.
%dirdata=zeros(m,1);
Px=@(t) 5 + 0*t;
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
'forcex',fx,'forcey',fy,'betax',bx,'BC',BC);
steady=true;
sol=StokesSolver(mesh,quad,Data,out,Map,steady);