clear all
close all
clc
solver_GM.SetPath;
solver_GM.Clean;
SIM=solver_GM;
SIM.SetTemp(0,0.0025,0.4,2./3,false);
%t0,dt,tf,theta,timedipendent
SIM.SetFisic(3,12,0.05,0.05,1);
%chi,kappa,cin,cest,mu
%SET PER STAZIONARIO (3,12,0.05,0.05,1)
SIM.SetGeom(3,3,1,0.1);
%lund2d,lun1d,larghezza,spess
%SET PER STAZIONARIO (3,3,1,0.2)
SIM.SetDD(1e-4,100,1);
%toll,maxiter,thetaDD
SIM.SetMetodo('krpo');
%(krpo)
%(naive)    MOLTO BRUTTO!!!
%(full2d)
SIM.SetCaseTest(1);
%(1)
SIM.Solve;
SIM.SetGraphics(0,false,45,20,0,0.13);
%treD,frontale,az,el,cinffig:cfr:krpo,csup
SIM.PlotSolution;
SIM.MakeMoreImg(1,'png');
