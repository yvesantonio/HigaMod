clc
close all
clear all

addpath ../../../Stokes
addpath ../../../Util
addpath ../../../FSI

Cx=struct('cinf',-0.1,'csup',2);
Cy=struct('cinf',-2,'csup',+2);
Cp=struct('cinf',-3,'csup',+7);

v=0.4:0.4:1.6;

Contour         = struct('out',0,'video',0,'visible','off','exp',[1,1.4,1.5,1.6,2,2.3]);
Profile         = struct('out',0,'video',0,'visible','off','exp',[],'x',0.5);
ProfileVelocity = struct('out',0,'video',0,'visible','off','exp',[1,1.4,1.5,1.6,2,2.3],'x',0.5);
Pressure        = struct('out',0,'video',0,'visible','off','exp',v);
Domain          = struct('out',1,'video',0,'visible','off','exp',v);

out=struct(...
    'Contour',Contour,'Profile',Profile,'Pressure',Pressure,'Domain',Domain,'ProfileVelocity',ProfileVelocity,...
    'Cx',Cx,'Cy',Cy,'Cp',Cp, ...
    'ended',false);

I=[0,10];
nb_elem=200;
mesh=struct('I',I,'nb_elem',nb_elem);

toll=1e-3;
m=7;
N=100;
T=2;
dt=T/N;
ntimestep=N;
nu=0.1;
K=50;
fx=@(x,y,t) 0+0*x+0*y+0*t;
fy=@(x,y,t) 0+0*x+0*y+0*t;

%Pinflow=@(t) 10*(t<T);
Pinflow=@(t) 5*((sin(2*pi/T*t)).^2*(t<T/4)+(t>=T/4));
%Pinflow=@(t) 5*((sin(2*pi/T*t)).^2)*( sin(2*pi/T*t)>0);
%Pinflow=@(t) 5+0*t;
Poutflow=@(t) 0*t;


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

sol=FSISolverNotSimGrad(mesh,quad,Data,out);
% load handel
% sound(y(1:20000),Fs)