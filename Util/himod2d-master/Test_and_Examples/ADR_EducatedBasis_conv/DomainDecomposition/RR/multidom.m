%% HI-MOD MONODOMINIO

close all
clear all

addpath ../../../../Core
addpath ../../../../Util

cutx    = [0,3,6];
m1=5;
m2=3;

size_mb  = [m1 m2];
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));
chi      = 3;
dato_dir_up=chi*0.05;
dato_dir_down=chi*0.05;

mu=1;


BC_laterali_up='rob';
BC_laterali_down='rob';

bc_up={ BC_laterali_up BC_laterali_up BC_laterali_up BC_laterali_up };
dato_up={dato_dir_up dato_dir_up dato_dir_up dato_dir_up};
bc_down={ BC_laterali_down BC_laterali_down BC_laterali_down BC_laterali_down};
dato_down={dato_dir_down dato_dir_down dato_dir_down dato_dir_down};

%**********************************************%
% Creazione di un profilo di Dirichlet per il  %
% dato all'inflow                              %
% Questo dato è compatibile con le condizioni  %
%**********************************************%
b=0.01;

if(strcmp(BC_laterali_up,'dir')&&strcmp(BC_laterali_down,'dir'))
    C=dato_dir_down;
    A=dato_dir_up-b-C;
else
    if(strcmp(BC_laterali_up,'dir')&&strcmp(BC_laterali_down,'rob'))
        C=(dato_down{1}+mu*b)/chi;
        A=dato_dir_up-b-C;
    else
        if(strcmp(BC_laterali_up,'rob')&&strcmp(BC_laterali_down,'rob'))
            C=(dato_down{1}+mu*b)/chi;
            A=(dato_up{1}-mu*b-b-(dato_down{1}-mu*b)/chi )/(2*mu+chi);
        end
    end
end

dato_dir = @(x) A*x.^2+b*x+C;

force = @(x,y) 10*(  ((x-1.5).^2+.4*(y-0.25).^2<0.01) + ((x-1.5).^2+.4*(y-0.75).^2<0.01) );
Dati = struct('dato_dir',dato_dir,'force',force);

a     = @(x)     0.*x; % Equazione del lato inferiore del canale
L     = @(x) 1 + 0.*x; % Spessore del canale
psi_x = @(x)     0.*x; % Riscalatura
Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) ( 20 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  0 + 0*x+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);

gamma=struct ('L',NaN,'R',NaN);
solver_DD(cutx,size_mb,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,gamma,true,'ND');