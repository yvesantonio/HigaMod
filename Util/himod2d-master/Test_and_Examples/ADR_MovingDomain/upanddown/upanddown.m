%% HI-MOD MONODOMINIO

close all
clear all
clc

addpath ../../../Core
addpath ../../../Util
addpath ../../../MapHandler

cutx     = [0,6];
m1       = 16;
size_mb  = m1;
hx       = 0.05;
mesh_x   = cutx(1):hx:cutx(2);
dt       = 0.1;
ntimestep= 50;
u0 = @(x,y) 0*x+0*y;

mu    = @(x,y,t) (  1 + 0*x+0*y +0*t); 
beta1 = @(x,y,t) (  10 + 0*x+0*y +0*t);
beta2 = @(x,y,t) (  0 + 0*x+0*y +0*t); 
sigma = @(x,y,t) (  1 + 0*x+0*y +0*t); 
force = @(x,y,t) 10*((x-0.5).^2+.4*y.^2<0.01)+0*t;

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',NaN);
umax=0.4;
tau=0.5;
haty=sym('(y-umax*sin(t/tau)/2)/(umax*sin(t/tau)/2+1)');
haty=subs(subs(haty,'umax',umax),'tau',tau);
[Map]=provideMapUnsteady(haty);

out=[10,20,30,40,50];

solver_unsteady(mesh_x, Coeff_forma, size_mb, force,Map,dt,ntimestep,u0,'dir','dir',out);