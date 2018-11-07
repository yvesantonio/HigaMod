clc
close all
clear all

addpath ../../../Stokes
addpath ../../../FSI
addpath ../../../MapHandler
addpath ../../../Util

Cx=struct('cinf',-0.5,'csup',0.5);
Cy=struct('cinf',-0.5,'csup',+0.5);
Cp=struct('cinf',-5,'csup',+5.);

nb_quiv=30;
out=struct('cont',1,'quiv',1,'profile',0,'pressure',1,'domain',1,   ...
    'videocont',1,'videoquiv',1,'videoprofile',0,'videopressure',1,'videodomain',1, ...
    'viscont','off','visquiver','off','visprofile','off','vispressure','off','visdomain','off',...
    'Cx',Cx,'Cy',Cy,'Cp',Cp,'nb_quiv',nb_quiv, ...
    'ended',false);

I=[0,20];
nb_elem=100;
mesh=struct('I',I,'nb_elem',nb_elem);

toll=1e-4;
m=6;
N=15;
T=1;
dt=T/N; 
ntimestep=6*N;
nu=0.01;
K=20;
d=2;
fx=@(x,y,t) 10*(2*t.*(t<T/4)+(1-2*t).*(t>=T/4).*(t<T/2)).*(-(x-I(2)/2).*(x-I(2)/2+d).*(x<I(2)/2).*(x>I(2)/2-d)+(x-I(2)/2).*(x-I(2)/2-d).*(x>I(2)/2).*(x<I(2)/2+d)+0*y+0*t);
fy=@(x,y,t) 0+0*x+0*y+0*t;
% Inflow data
essential_inflow=false;
%Pinflow=@(t) 10*(t<T);
%Pinflow=@(t) 10*(sin(2*pi/T*t)).^2*(t<T/4)+10*(t>=T/4);
Pinflow=@(t)  0*t;
Poutflow=@(t) 0*t;
dir_data=zeros(m-2);
dir_data(1)=1;
%nel codice è imposta u_y = 0 all'inflow
Data=struct('nu',nu,'m',m,'dt',dt,'ntimestep',ntimestep,...
    'forcex',fx,'forcey',fy,...
    'essential_inflow',essential_inflow,'inflow_dirichlet',dir_data,...
    'Pinf',Pinflow,'Pout',Poutflow,'K',K,'toll',toll);

nb_quadnode_y=32; %quadratura su tutto l'intervallo -1,1
nb_quadnode_x=4;  %quadratura sul singolo elemento finito
quad=struct('nx',nb_quadnode_x,'ny',nb_quadnode_y);

sol=FSISolver(mesh,quad,Data,out);
load handel
sound(y(1:20000),Fs)
