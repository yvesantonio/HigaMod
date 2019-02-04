    close all
clear

%%%%%%%%
% Slab %
%%%%%%%%
        
% PARAMETERS

L = 1;
W = 1;
H = 5;

% CREATE NURBS VOLUME

lin = nrbline([0 0 0], [ L 0 0 ]);
srf = nrbextrude(lin, [ 0 W 0 ]);
VOL = nrbextrude(srf, [ 0 0 H ]);

% CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE

geometry = geo_load(VOL);

% EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
% AND ITS JACOBIAN AND HESSIAN MATRICES

map = @(x,y,z) geometry.map([x;y;z]);
Jac = @(x,y,z) geometry.map_der([x;y;z]);
Hes = @(x,y,z) geometry.map_der2([x;y;z]);

geometricInfo.geometry = geometry;
geometricInfo.map = @(x,y,z) map(x,y,z);
geometricInfo.Jac = @(x,y,z) Jac(x,y,z);
geometricInfo.Hes = @(x,y,z) Hes(x,y,z);
geometricInfo.Type = 'Slab_3';
geometricInfo.W = W;
geometricInfo.L = L;
geometricInfo.H = H;
        
%%%%%%%%%%%%%%%%%%%%
% CREATE .STL FILE %
%%%%%%%%%%%%%%%%%%%%

type = geometricInfo.Type;

Nx = 1e2;
Ny = 2e1;
Nz = 2e1;

evalX = linspace(0,1,Nx);
evalY = linspace(0,1,Ny);
evalZ = linspace(0,1,Nz);

tic;
[X,Y,Z] = mapOut3D(evalX,evalY,evalZ,map,type);
t1 = toc;

[evalJac,structPhi,evalDetJac] = jacOut3D(evalX,evalY,evalZ,Jac,type);

auxX = Nx;
auxY = Ny;
auxZ = Nz;
Nx = auxZ;
Ny = auxY;
Nz = auxX;

evalJ11 = zeros(Nx,Ny,Nz);
evalJ12 = zeros(Nx,Ny,Nz);
evalJ13 = zeros(Nx,Ny,Nz);
evalJ21 = zeros(Nx,Ny,Nz);
evalJ22 = zeros(Nx,Ny,Nz);
evalJ23 = zeros(Nx,Ny,Nz);
evalJ31 = zeros(Nx,Ny,Nz);
evalJ32 = zeros(Nx,Ny,Nz);
evalJ33 = zeros(Nx,Ny,Nz);

for ii = 1:Nx
    for jj = 1:Ny
        for kk = 1:Nz
            evalJ11(ii,jj,kk) = evalJac{ii,jj,kk}(1,1);
            evalJ12(ii,jj,kk) = evalJac{ii,jj,kk}(1,2);
            evalJ13(ii,jj,kk) = evalJac{ii,jj,kk}(1,3);
            evalJ21(ii,jj,kk) = evalJac{ii,jj,kk}(2,1);
            evalJ22(ii,jj,kk) = evalJac{ii,jj,kk}(2,2);
            evalJ23(ii,jj,kk) = evalJac{ii,jj,kk}(2,3);
            evalJ31(ii,jj,kk) = evalJac{ii,jj,kk}(3,1);
            evalJ32(ii,jj,kk) = evalJac{ii,jj,kk}(3,2);
            evalJ33(ii,jj,kk) = evalJac{ii,jj,kk}(3,3);
        end
    end
end

disp('FINISHED BASE MAPPING')

[XBase,YBase,ZBase] = ndgrid(evalZ,evalY,evalX);

tic;
Phi1 = griddedInterpolant(XBase,YBase,ZBase,X,'spline');
Phi2 = griddedInterpolant(XBase,YBase,ZBase,Y,'spline');
Phi3 = griddedInterpolant(XBase,YBase,ZBase,Z,'spline');

J11 = griddedInterpolant(XBase,YBase,ZBase,evalJ11,'spline');
J12 = griddedInterpolant(XBase,YBase,ZBase,evalJ12,'spline');
J13 = griddedInterpolant(XBase,YBase,ZBase,evalJ13,'spline');
J21 = griddedInterpolant(XBase,YBase,ZBase,evalJ21,'spline');
J22 = griddedInterpolant(XBase,YBase,ZBase,evalJ22,'spline');
J23 = griddedInterpolant(XBase,YBase,ZBase,evalJ23,'spline');
J31 = griddedInterpolant(XBase,YBase,ZBase,evalJ31,'spline');
J32 = griddedInterpolant(XBase,YBase,ZBase,evalJ32,'spline');
J33 = griddedInterpolant(XBase,YBase,ZBase,evalJ33,'spline');

detJ = griddedInterpolant(XBase,YBase,ZBase,evalDetJac,'spline');

intrpPhi1_dx = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi1_dx,'spline');
intrpPhi1_dy = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi1_dy,'spline');
intrpPhi1_dz = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi1_dz,'spline');
intrpPhi2_dx = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi2_dx,'spline');
intrpPhi2_dy = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi2_dy,'spline');
intrpPhi2_dz = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi2_dz,'spline');
intrpPhi3_dx = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi3_dx,'spline');
intrpPhi3_dy = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi3_dy,'spline');
intrpPhi3_dz = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi3_dz,'spline');

t2 = toc;

geometricInfo.Phi1 = Phi1;
geometricInfo.Phi2 = Phi2;
geometricInfo.Phi3 = Phi3;
geometricInfo.J11  = J11;
geometricInfo.J12  = J12;
geometricInfo.J13  = J13;
geometricInfo.J21  = J21;
geometricInfo.J22  = J22;
geometricInfo.J23  = J23;
geometricInfo.J31  = J31;
geometricInfo.J32  = J32;
geometricInfo.J33  = J33;
geometricInfo.detJ  = detJ;

geometricInfo.Phi1_dx = intrpPhi1_dx;
geometricInfo.Phi1_dy = intrpPhi1_dy;
geometricInfo.Phi1_dz = intrpPhi1_dz;
geometricInfo.Phi2_dx = intrpPhi2_dx;
geometricInfo.Phi2_dy = intrpPhi2_dy;
geometricInfo.Phi2_dz = intrpPhi2_dz;
geometricInfo.Phi3_dx = intrpPhi3_dx;
geometricInfo.Phi3_dy = intrpPhi3_dy;
geometricInfo.Phi3_dz = intrpPhi3_dz;

disp('FINISHED CREATING THE MAP INTERPOLANT')

% DISPLAY SIMULATION TIMES

disp(['Number of D.O.F.s    : ',num2str(Nx * Ny * Nz)])
disp(['Exec. Time Base Map  : ',num2str(t1),' [sec]'])
disp(['Exec. Time Interpol  : ',num2str(t2),' [sec]'])

% EXPORT GEOMETRY DATA

filename = 'Slab3DGeo_3.mat';
geo = VOL;
save(filename,'geo')

filename = 'Slab3DMap_3.mat';
save(filename,'geometricInfo')

% PLOT THE NURBS VOLUME

figure;
nrbplot(VOL,[10 10 10]);