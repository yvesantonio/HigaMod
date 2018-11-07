%% HI-MOD MONODOMINIO

close all
clear all

addpath ../../../Core
addpath ../../../Util
addpath ../../../MapHandler

%**********************************************%
% setting number of domains, mesh, dof & modes %
%**********************************************%
cutx    = [0,6];
m1       = 4;
size_mb  = m1;
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));    
mesh_x   = cutx(1):hx:cutx(2);

%**********************************************%
% BC laterali:                                 %
%        * 'rob': mu du/dnu + chi u = dato     %
%        * 'dir': u=dato		       %
% Attenzione: funziona bene solo nel caso si   %
% utilizzino le stesse condizioni di bordo     %
% laterali su tutti i domini                   %
%**********************************************%
chi = 1;
%**********************************************%
% forzante:                                    %
%**********************************************%
force = @(x,y) 10*(  ((x-2).^2+.4*(y-0.75).^2<0.01) + ((x-1).^2+.4*(y+0.2).^2<0.01) );
%**********************************************%

%**********************************************%
% Grafico del dato in ingresso                 %
%**********************************************%
% fplot(dato_dir,[0,1]);
% title('Dato in ingresso');
% xlabel('y');
% ylabel('u_{in}');
%**********************************************%


%**********************************************%
% Coefficienti della forma bilineare           %
%**********************************************%
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) (  20 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  0 + 0*x+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);
%*******************%
% Map of the domain %
%*******************%
theta=pi/12;
haty=sym('2*(y+1-x*tantheta)/(2-x*tantheta)-1');
haty=subs(haty,'tantheta',tan(theta));
%haty=sym('y');
[Map]=provideMap(haty);
solver_leg(mesh_x, Coeff_forma, size_mb, force,Map,'rob','rob');
