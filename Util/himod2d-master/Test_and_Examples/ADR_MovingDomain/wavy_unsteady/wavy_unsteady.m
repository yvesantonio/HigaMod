%% HI-MOD MONODOMINIO

close all
clear all
clc

addpath ../../../Core
addpath ../../../Util
addpath ../../../MapHandler

%**********************************************%
% setting number of domains, mesh, dof & modes %
%**********************************************%
cutx    = [0,5];
m1       = 4;
size_mb  = m1;
nd       = length(size_mb);
hx       = 0.02*ones(size(size_mb));
mesh_x   = cutx(1):hx:cutx(2);
dt=0.01;
ntimestep=100;
out=[0.2,0.4,0.6,0.8,1]/dt;
chi = 2;
force = @(x,y,t) 10*((x-0.5).^2+.4*y.^2<0.01)+0*t;

mu    = @(x,y,t) (  1 + 0*x+0*y +0*t); % Diffusione
beta1 = @(x,y,t) (  10 + 0*x+0*y +0*t); % Trasporto orizzontale
beta2 = @(x,y,t) (  0 + 0*x+0*y +0*t); % Trasporto verticale
sigma = @(x,y,t) (  1 + 0*x+0*y +0*t); % Reazione
u0=@(x,y) 0*x+0*y;
Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);
haty=sym('y/(0.5*sin(50*t)*cos(pi/3*x+pi/2)+1)');
[Map]=provideMapUnsteady(haty);

solver_unsteady(mesh_x, Coeff_forma, size_mb, force,Map,dt,ntimestep,u0,'dir','dir',out);   