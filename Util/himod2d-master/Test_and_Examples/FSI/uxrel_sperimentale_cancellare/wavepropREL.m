clc
close all
clear all

addpath ../../core
addpath ../../util
addpath ../../FSI
addpath ../../../MapHandler
addpath ../../../Util

Cx=struct('cinf',-0.1,'csup',5);
Cy=struct('cinf',-1,'csup',+1);
Cp=struct('cinf',-3,'csup',+15.);

nb_quiv=15;
out=struct('cont',1,'quiv',0,'profile',1,'pressure',1,'domain',1,   ...
    'videocont',1,'videoquiv',0,'videoprofile',1,'videopressure',1,'videodomain',1, ...
    'viscont','on','visquiver','off','visprofile','off','vispressure','off','visdomain','off',...
    'Cx',Cx,'Cy',Cy,'Cp',Cp,'nb_quiv',nb_quiv, ...
    'ended',false);

I=[0,2.5];
nb_elem=50;
mesh=struct('I',I,'nb_elem',nb_elem);

toll=1e-5;
m=6;
N=40;
T=2;
dt=T/N;
ntimestep=N;
nu=1;
K=10;
fx=@(x,y,t) 0+0*x+0*y+0*t;
fy=@(x,y,t) 0+0*x+0*y+0*t;
% Inflow data
essential_inflow=false;
%Pinflow=@(t) 10*(t<T);
Pinflow=@(t) 10*(sin(2*pi/T*t)).^2*(t<T/4)+10*(t>=T/4);
%Pinflow=@(t) 10+0*t;
Poutflow=@(t) 0*t;
dir_data=zeros(m);
%dir_data(3)=1;
%nel codice è imposta u_y = 0 all'inflow
Data=struct('nu',nu,'m',m,'dt',dt,'ntimestep',ntimestep,...
    'forcex',fx,'forcey',fy,...
    'essential_inflow',essential_inflow,'inflow_dirichlet',dir_data,...
    'Pinf',Pinflow,'Pout',Poutflow,'K',K,'toll',toll);

nb_quadnode_y=32; %quadratura su tutto l'intervallo -1,1
nb_quadnode_x=4;  %quadratura sul singolo elemento finito
quad=struct('nx',nb_quadnode_x,'ny',nb_quadnode_y);

spl=FSISolver_x_released(mesh,quad,Data,out);