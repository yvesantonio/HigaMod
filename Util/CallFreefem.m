close all
clear all
clc

mu = @(x,y) 1;
beta1 = @(x,y) 0;
beta2 = @(x,y) 0;
sigma = @(x,y) 0;
chi = @(x,y) 1;

Coeff = struct('mu',mu,'beta1',beta1,'beta2',beta2,'sigma',sigma,'chi',chi);
                   
up = 1;
down = 0;
left = 0;
right = 2;

force = @(x,y) -(4*x.^3-1).*sin(2*pi*(y+0.5)) + (0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5));

Domain = struct('up',up,'down',down,'left',left,'right',right);

numbHorNodes = 100;
numbVerNodes = 100;

Mesh = struct('numbHorNodes',numbHorNodes,'numbVerNodes',numbVerNodes);

[L2,H1] = analysisFreeFem(Domain, Mesh,Coeff,force);

disp(L2)
disp(H1)