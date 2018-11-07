%% GM-STAZIONARIO
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
%SET PER STAZIONARIO (3,12,0.05,0.05,1)

CASOTEST.SetGeom(3,1,1,0.01);
%lund2d,lun1d,larghezza,spess
%SET per stazionario (3,3,1,0.01)

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

%Per foto python
!cp *.out ./proveGM/input/
!cp *.msh ./proveGM/input/
!cp u2dkmeno1 ./proveGM/input/
!cp dati_geometrici ./proveGM/input/



% %% GM-NON STAZIONARIO
% clear all
% close all
% clc
% solver_GM.SetPath;
% solver_GM.Clean;
% CASOTEST=solver_GM;
% 
% CASOTEST.SetTemp(0,0.01,0.4,2./3,true);
% %t0,dt,tf,theta,timedipendent
% 
% CASOTEST.SetFisic(0.5,12,0.3,0.0,0.05); 
% %chi,kappa,cin,cest,mu
% %SET per non stazionario (0.5,12,0.3,0.0,0.05)
% 
% CASOTEST.SetGeom(3,2,1,0.01);
% %lund2d,lun1d,larghezza,spess
% %SET per non stazionario (3,2,1,0.01)
% 
% CASOTEST.SetDD(1e-4,100,0.8);
% %toll,maxiter,thetaDD
% 
% CASOTEST.SetMetodo('krpo');
% %(krpo)(naive)(full2d)
% 
% CASOTEST.SetCaseTest(3);
% %NON STAZIONARIO
% %(1) Scalino prima dell'ostruzione
% %(2) Scalino dopo l'ostruzione
% %(3) Macchia
% 
% CASOTEST.Solve;
% CASOTEST.SetGraphics(0,false,45,20,-0.01,0.32);
% %treD,frontale,az,el,cinf,csup
% 
% CASOTEST.PlotSolution;
