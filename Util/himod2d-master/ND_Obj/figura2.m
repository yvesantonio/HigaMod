%% full2d
clear all
close all
clc
solver_GM.SetPath;
solver_GM.Clean;
CASOTEST=solver_GM;

CASOTEST.SetTemp(0,0.01,0.4,2./3,false);
%t0,dt,tf,theta,timedipendent

CASOTEST.SetFisic(1,12,0.05,0.02,1); 
%chi,kappa,cin,cest,mu

CASOTEST.SetGeom(5,5,1,0.01);
%lund2d,lun1d,larghezza,spess

CASOTEST.SetDD(1e-4,100,0.8);
%toll,maxiter,thetaDD

CASOTEST.SetMetodo('full2d');
%(krpo)(naive)(full2d)

CASOTEST.SetCaseTest(1);
%STAZIONARIO
%(1) Due camini Bisogna cambiare a mano la forzante nel file ND2d_staz.edp

CASOTEST.Solve;
CASOTEST.SetGraphics(0,false,45,15,0.0,0.11);
%treD,frontale,az,el,cinf,cs1up

CASOTEST.PlotSolution;

%%
!mv *.out ./figura2/input/
!mv mesh_full2d.msh ./figura2/input/mesh_full2d.msh
!rm *.msh
!mv u2dkmeno1 ./figura2/input/u2dkmeno1full2d
!mv dati_geometrici ./figura2/input/dati_geometricifull2d
!mv media* ./figura2/input/sigma
%% krpogood
clear all
close all
clc
solver_GM.SetPath;
solver_GM.Clean;
CASOTEST=solver_GM;

CASOTEST.SetTemp(0,0.01,0.4,2./3,false);
%t0,dt,tf,theta,timedipendent

CASOTEST.SetFisic(1,12,0.05,0.02,1); 
%chi,kappa,cin,cest,mu

CASOTEST.SetGeom(4,6,1,0.01);
%lund2d,lun1d,larghezza,spess

CASOTEST.SetDD(1e-4,100,0.8);
%toll,maxiter,thetaDD

CASOTEST.SetMetodo('krpo');
%(krpo)(naive)(full2d)

CASOTEST.SetCaseTest(1);
%STAZIONARIO
%(1) Due camini Bisogna cambiare a mano la forzante nel file ND2d_staz.edp

CASOTEST.Solve;
CASOTEST.SetGraphics(0,false,45,15,0.0,0.11);
%treD,frontale,az,el,cinf,cs1up

CASOTEST.PlotSolution;

!mv matrice.out ./figura2/input/matricegood.out
!mv xcoordinate.out ./figura2/input/xcoordinategood.out
!mv ycoordinate.out ./figura2/input/ycoordinategood.out
!mv mesh_2d.msh ./figura2/input/mesh_2dgood.msh
!mv mesh_1d.msh ./figura2/input/mesh_1dgood.msh
!rm *.msh
!mv u2dkmeno1 ./figura2/input/u2dkmeno1good
!mv dati_geometrici ./figura2/input/dati_geometricigood

%% krpobad
clear all
close all
clc
solver_GM.SetPath;
solver_GM.Clean;
CASOTEST=solver_GM;

CASOTEST.SetTemp(0,0.01,0.4,2./3,false);
%t0,dt,tf,theta,timedipendent

CASOTEST.SetFisic(1,12,0.05,0.02,1); 
%chi,kappa,cin,cest,mu

CASOTEST.SetGeom(2.5,7.5,1,0.01);
%lund2d,lun1d,larghezza,spess

CASOTEST.SetDD(1e-4,100,0.8);
%toll,maxiter,thetaDD

CASOTEST.SetMetodo('krpo');
%(krpo)(naive)(full2d)

CASOTEST.SetCaseTest(1);
%STAZIONARIO
%(1) Due camini Bisogna cambiare a mano la forzante nel file ND2d_staz.edp

CASOTEST.Solve;
CASOTEST.SetGraphics(0,false,45,15,0.0,0.11);
%treD,frontale,az,el,cinf,cs1up

CASOTEST.PlotSolution;

!mv matrice.out ./figura2/input/matricebad.out
!mv xcoordinate.out ./figura2/input/xcoordinatebad.out
!mv ycoordinate.out ./figura2/input/ycoordinatebad.out
!mv mesh_2d.msh ./figura2/input/mesh_2dbad.msh
!mv mesh_1d.msh ./figura2/input/mesh_1dbad.msh
!rm *.msh
!mv u2dkmeno1 ./figura2/input/u2dkmeno1bad
!mv dati_geometrici ./figura2/input/dati_geometricibad

%% HI-MOD MONODOMINIO

close all
clear all
clc

addpath ../Core
addpath ../Util

cutx    = [0,4,10];
m1=5;
m2=2;

size_mb  = [m1 m2];
nd       = length(size_mb);
hx       = 0.01*ones(size(size_mb));
chi      = 1;
dato_dir_up=chi*0.02;
dato_dir_down=chi*0.02;

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

force = @(x,y) 10*(((x-1.5).^2+.4*(y-0.5).^2<0.01));
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
!mv *2.out figura2/input/himo52.out
close all