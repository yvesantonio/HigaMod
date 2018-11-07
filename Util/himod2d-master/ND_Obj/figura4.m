%% full2d
clear all
close all
clc
solver_GM.SetPath;
solver_GM.Clean;
CASOTEST=solver_GM;

CASOTEST.SetTemp(0,0.01,0.4,2./3,false);
%t0,dt,tf,theta,timedipendent

CASOTEST.SetFisic(3,12,0.05,0.05,1); 
%chi,kappa,cin,cest,mu

CASOTEST.SetGeom(3,3,1,0.01);
%lund2d,lun1d,larghezza,spess

CASOTEST.SetDD(1e-4,100,0.8);
%toll,maxiter,thetaDD

CASOTEST.SetMetodo('full2d');
%(krpo)(naive)(full2d)

CASOTEST.SetCaseTest(1);
%STAZIONARIO
%(1) Due camini Bisogna cambiare a mano la forzante nel file ND2d_staz.edp
%(2) Un camino

CASOTEST.Solve;
CASOTEST.SetGraphics(0,false,45,15,0.0,0.16);
%treD,frontale,az,el,cinf,cs1up

CASOTEST.PlotSolution;

%%
!mv *.out ./figura4/input/
!mv mesh_full2d.msh ./figura4/input/mesh_full2d.msh
!rm *.msh
!mv u2dkmeno1 ./figura4/input/u2dkmeno1full2d
!mv dati_geometrici ./figura4/input/dati_geometricifull2d

%% KRPO
clear all
close all
clc
solver_GM.SetPath;
solver_GM.Clean;
CASOTEST=solver_GM;

CASOTEST.SetTemp(0,0.01,0.4,2./3,false);
%t0,dt,tf,theta,timedipendent

CASOTEST.SetFisic(3,12,0.05,0.05,1); 
%chi,kappa,cin,cest,mu

CASOTEST.SetGeom(3,3,1,0.01);
%lund2d,lun1d,larghezza,spess

CASOTEST.SetDD(1e-4,100,0.8);
%toll,maxiter,thetaDD

CASOTEST.SetMetodo('krpo');
%(krpo)(naive)(full2d)

CASOTEST.SetCaseTest(1);
%STAZIONARIO
%(1) Due camini Bisogna cambiare a mano la forzante nel file ND2d_staz.edp
%(2) Un camino

CASOTEST.Solve;
CASOTEST.SetGraphics(0,false,45,15,0.0,0.16);
%treD,frontale,az,el,cinf,cs1up

CASOTEST.PlotSolution;

!mv *.out ./figura4/input/
!mv mesh_1d.msh ./figura4/input/mesh_1d.msh
!mv mesh_2d.msh ./figura4/input/mesh_2d.msh
!rm *.msh
!mv u2dkmeno1 ./figura4/input/u2dkmeno1krpo
!mv dati_geometrici ./figura4/input/dati_geometricikrpo


%% HI-MOD MONODOMINIO

close all
clear all

addpath ../Core
addpath ../Util

cutx    = [0,3,6];
m1=3;
m2=1;

size_mb  = [m1 m2];
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));
chi      = 3;
dato_dir_up=chi*0.05;
dato_dir_down=chi*0.05;

mu=1;


BC_laterali_up='rob';
BC_laterali_down='rob';

bc_up={ BC_laterali_up BC_laterali_up};
dato_up={dato_dir_up dato_dir_up};
bc_down={ BC_laterali_down BC_laterali_down};
dato_down={dato_dir_down dato_dir_down};

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
!mv *1.out figura4/input/himo31.out

%
% HI-MOD MONODOMINIO

close all
clear all

addpath ../Core
addpath ../Util

cutx    = [0,3,6];
m1=5;
m2=1;

size_mb  = [m1 m2];
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));
chi      = 3;
dato_dir_up=chi*0.05;
dato_dir_down=chi*0.05;

mu=1;


BC_laterali_up='rob';
BC_laterali_down='rob';

bc_up={ BC_laterali_up BC_laterali_up};
dato_up={dato_dir_up dato_dir_up};
bc_down={ BC_laterali_down BC_laterali_down};
dato_down={dato_dir_down dato_dir_down};

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
!mv *1.out figura4/input/himo51.out
close all