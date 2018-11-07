%% HI-MOD MULTIDOMINIO

close all
clear all

addpath ../../../Core
addpath ../../../Util
%**********************************************%
%    FUNZIONA ANCHE CON UN SOLO DOMINIO        %
%**********************************************%
% setting number of domains, mesh, dof & modes %
%**********************************************%
cutx     = [0,1,3,5,6];
m1       = 7;
m2       = 7;
m3       = 7;
m4		 = 7;
size_mb  = [m1,m2,m3,m4];
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));

%**********************************************%
% BC laterali:                                 %
%        * 'rob': mu du/dnu + chi u = dato     %
%        * 'dir': u=dato		       %
% Attenzione: funziona bene solo nel caso si   %
% utilizzino le stesse condizioni di bordo     %
% laterali su tutti i domini                   %
%**********************************************%

chi      = 3;
dato_dir_up=0.05;
dato_dir_down=0.05;

mu=1;


BC_laterali_up='dir';
BC_laterali_down='dir';

bc_up={ BC_laterali_up BC_laterali_up BC_laterali_up BC_laterali_up};
dato_up={dato_dir_up dato_dir_up dato_dir_up dato_dir_up};
bc_down={ BC_laterali_down BC_laterali_down BC_laterali_down BC_laterali_down};
dato_down={dato_dir_down dato_dir_down dato_dir_down dato_dir_down};

%**********************************************%
% Creazione di un profilo di Dirichlet per il  %
% dato all'inflow                              %
% Questo dato è compatibile con le condizioni %
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
%**********************************************%
% forzante:                                    %
%**********************************************%
%DUE CAMINI
force = @(x,y) 10*(  ((x-3).^2+.4*(y-0.25).^2<0.01) + ((x-3).^2+.4*(y-0.75).^2<0.01) );
%UN CAMINO
%force = @(x,y) 10*((x-0.5).^2+.4*(y-0.5).^2<0.01);
%DUE CAMINI SPOSTATI
%force=@(x,y) 10*(  ((x-1).^2+.4*(y-0.25).^2<0.01) + ((x-1.5).^2+.4*(y-0.75).^2<0.01) );
%RETTANGOLO
%force=@(x,y) 3*((x<2).*(x>1).*(y<0.75).*(y>0.25));

%**********************************************%
Dati = struct('dato_dir',dato_dir,'force',force);

%**********************************************%
% Grafico del dato in ingresso                 %
%**********************************************%
% fplot(dato_dir,[0,1]);
% title('Dato in ingresso');
% xlabel('y');
% ylabel('u_{in}');
%**********************************************%

%**********************************************%
% Dati relativi alla forma del dominio         %
% ATTENZIONE: il codice funziona solo per L=1  %
% e psi_x=0 in realtà per L!=1 e psi_x!=0      %
% vanno ricontrollati tutti i coefficienti.    %
% Soprattutto in assembla_x                    %
%**********************************************%

a     = @(x)     0.*x; % Equazione del lato inferiore del canale
L     = @(x) 1 + 0.*x; % Spessore del canale
psi_x = @(x)     0.*x; % Riscalatura
Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);
%**********************************************%


%**********************************************%
% Coefficienti della forma bilineare           %
%**********************************************%

%larghezza=1;
%kappa=2;
%delta=0.5;
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
%beta1 = @(x,y) 20*(kappa+2)/(kappa+2*(1-delta))*(1-delta*((y-0.5)*2/larghezza).^kappa); 
beta1 = @(x,y) ( -20 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  0 + 0*x+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);
%**********************************************%

      	gammaL=25;
        gammaR=-1;   
		if(beta1(0,0)<0)
			gammaL=-1;
			gammaR=25;
		end
        gamma=struct ('L',gammaL,'R',gammaR);
        %***********************************************%
        %                    SOLVER                     %
        %***********************************************%
		% here you can make a comparison between the different type of domain decomposition coupling available.
        % solver_DD(cutx,size_mb,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,gamma,true,'DN');
        % solver_DD(cutx,size_mb,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,gamma,true,'ND');
        solver_DD(cutx,size_mb,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,gamma,false,'RR');
