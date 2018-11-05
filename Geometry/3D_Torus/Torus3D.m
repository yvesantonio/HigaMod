close all
clear

%%%%%%%%%
% Torus %
%%%%%%%%%
        
% PARAMETERS

ang  = 2 * pi;
ang2 = 4 * pi / 4;
R    = 3;
L    = 10;

% CREATE NURBS VOLUME

knots = {[0 0 0 1 1 1], [0 0 0 1 1 1]};
coefs = zeros (4, 3, 3);
coefs(:, 1, 1) = [1, 0, 0, 1];
coefs(:, 2, 1) = [1/sqrt(2), -1/sqrt(2), 0 , 1/sqrt(2)];
coefs(:, 3, 1) = [0, -1, 0, 1];
coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 2, 2) = [0, 0, 0, sqrt(2)-1];
coefs(:, 3, 2) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 1, 3) = [0, 1, 0, 1];
coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 3, 3) = [-1, 0, 0, 1];

srf = nrbmak (coefs, knots);
srf = nrbtform (srf, vecscale ([R R 1]));

vol = nrbrevolve(srf, [ 0 L 0 ], [ 1 0 0 ], ang2);

% CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE

geometry = geo_load(vol);

% EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
% AND ITS JACOBIAN AND HESSIAN MATRICES

map = @(x,y,z) geometry.map([x;y;z]);
Jac = @(x,y,z) geometry.map_der([x;y;z]);
Hes = @(x,y,z) geometry.map_der2([x;y;z]);

geometricInfo.geometry = geometry;
geometricInfo.map = @(x,y,z) map(x,y,z);
geometricInfo.Jac = @(x,y,z) Jac(x,y,z);
geometricInfo.Hes = @(x,y,z) Hes(x,y,z);
geometricInfo.Type = 'Torus';
geometricInfo.ang = ang;
geometricInfo.ang2 = ang2;
geometricInfo.R = R;
geometricInfo.L = L;
        
% EXPORT GEOMETRY DATA

filename = 'Torus3DGeo.mat';
geo = vol;
save(filename,'geo')

filename = 'Torus3DMap.mat';
save(filename,'geometricInfo')

%%%%%%%
% PLOTS
%%%%%%%

% PLOT THE NURBS VOLUME


% PLOT THE FINAL VOLUME

figure
nrbplot(vol,[30 30 30]);