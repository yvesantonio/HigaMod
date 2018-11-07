%% HI-MOD MONODOMINIO

close all
clear all
clc
addpath ../../../Core
addpath ../../../Util

cutx    = [0,1.5,5];
m1       = [4,3];
size_mb  = m1;
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));

dato_dir_up=0.0;
dato_dir_down=0.0;

BC_laterali_up='dir';
BC_laterali_down='dir';

bc_up={ BC_laterali_up BC_laterali_up};
dato_up={dato_dir_up dato_dir_up};
bc_down={ BC_laterali_down BC_laterali_down};
dato_down={dato_dir_up dato_dir_up};

dato_dir=@(y) 0.05*y.*(1-y);
force = @(x,y) 1.3*(0.5*(x-1).^2+.4*(y-0.5).^2<0.03);
               
Dati = struct('dato_dir',dato_dir,'force',force);


a     = @(x)     0.*x; % Equazione del lato inferiore del canale
L     = @(x) 1 + 0.*x; % Spessore del canale
psi_x = @(x)     0.*x; % Riscalatura
Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

%**********************************************%
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) (  2 +  3*x+0*y); % Trasporto orizzontale
beta2 = @(x,y) (  4*cos(2*pi*x/3)+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',NaN);
gammaL=25;
gammaR=-1;
gamma=struct ('L',gammaL,'R',gammaR);
solver_DD(cutx,size_mb,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,gamma,true,'RR');
axis equal