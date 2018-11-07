%% Identity.m
close all
clear all

addpath ../../../Core
addpath ../../../Util
addpath ../../../MapHandler

cutx    = [0,6];
m1       = 16;
size_mb  = m1;
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));
mesh_x   = cutx(1):hx:cutx(2);

force = @(x,y) 1.3*(0.5*(x-1).^2+.4*y.^2<0.05);

%**********************************************%
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) (  2 +  3*x+0*y); % Trasporto orizzontale
beta2 = @(x,y) (  2*cos(2*pi*x/3)+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',NaN);

haty=sym('y');
[Map]=provideMap(haty);
solver_leg(mesh_x, Coeff_forma, size_mb, force,Map,'dir','dir');