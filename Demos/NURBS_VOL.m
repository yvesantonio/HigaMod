%%%%%%%%%%%%
% CYLINDER %
%%%%%%%%%%%%

close all
% 
% %% CREATE A SIMPLE CYLINDER
% 
% % PARAMETERS
% 
% ang = 2*pi;
% R   = 10;
% L   = 3;
% 
% % CREEATE NURBS VOLUME
% 
% crv = nrbline([ 0 0 0 ],[ R 0 0 ]);
% srf = nrbrevolve(crv, [ 0 0 0 ], [ 0 0 1 ], ang);
% vol = nrbextrude(srf, [ 0 0 L ]);
% 
% % PLOT THE FINAL VOLUME
% 
% nrbplot(vol,[30 20 10]);
% 
% %% TEST THE MAP FUNCTION EXTRACTED
% 
% % EXTRACT THE GEOMETRY INFORMATION FROM THE GENERATED VOLUME
% 
% geometry = geo_load(vol);
% 
% map = @(x,y,z) geometry.map([x;y;z]);
% Jac = @(x,y,z) geometry.map_der([x;y;z]);
% Hes = @(x,y,z) geometry.map_der2([x;y;z]);
% 
% % Coordinates defined in the refenrece domain
% 
% phi = linspace(0,1,10);
% theta = linspace(0,1,5);
% radius = linspace(0,1,5);
% 
% [RR,TT,EE] = meshgrid(phi,theta,radius);
% 
% % Map into the physical domain
% 
% XX = mapOut3D(phi,theta,radius,map,1,'Torus');
% YY = mapOut3D(phi,theta,radius,map,2,'Torus');
% ZZ = mapOut3D(phi,theta,radius,map,3,'Torus');
% 
% [m,n,o] = size(XX);
% mySize  = m * n * o;
% 
% X = reshape(XX,[1,mySize]);
% Y = reshape(YY,[1,mySize]);
% Z = reshape(ZZ,[1,mySize]);
% 
% % Plot the obtained points in the physical domain
% 
% % % figure;
% % % plot3(X,Y,Z,'o');
% 
% %% CREATE A QUARTER OF A TORUS
% 
% % PARAMETERS
% 
% ang  = 2 * pi;
% ang2 = 2 * pi / 4;
% R    = 3;
% L    = 10;
% 
% % CREATE NURBS VOLUME
% 
% knots = {[0 0 0 1 1 1], [0 0 0 1 1 1]};
% coefs = zeros (4, 3, 3);
% coefs(:, 1, 1) = [1, 0, 0, 1];
% coefs(:, 2, 1) = [1/sqrt(2), -1/sqrt(2), 0 , 1/sqrt(2)];
% coefs(:, 3, 1) = [0, -1, 0, 1];
% coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 2, 2) = [0, 0, 0, sqrt(2)-1];
% coefs(:, 3, 2) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 1, 3) = [0, 1, 0, 1];
% coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 3, 3) = [-1, 0, 0, 1];
% 
% srf = nrbmak (coefs, knots);
% srf = nrbtform (srf, vecscale ([R R 1]));
% 
% vol = nrbrevolve(srf, [ 0 L 0 ], [ 1 0 0 ], ang2);
% 
% % PLOT THE FINAL VOLUME
% 
% fig1 = figure;
% nrbplot(vol,[30 30 30]);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])
% figure.XRuler.Axle.LineStyle = 'none';  
% axis off
% 
% fig2 = figure;
% nrbkntplot(vol);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])
% 
% fig3 = figure;
% nrbctrlplot(vol);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])
% 
% % Refine the volume
% 
% newKnots = linspace(0,1,3);
% newKnots = newKnots(2:end-1);
% volRefined = nrbkntins(vol, {newKnots, newKnots, newKnots});
% 
% % Plot refined volume
% 
% fig4 = figure;
% nrbkntplot(volRefined);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])
% 
% fig5 = figure;
% nrbctrlplot(volRefined);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])

% Modify the control points of the refined volume 

% close all
% 
% index = [3 , 3 , 3];
% move = [2 , 2 , 2];
% nrb = nrbmodp (volRefined, move, index);
% nrb = nrbmodw (nrb, new_w, index);
% 
% fig4 = figure;
% nrbkntplot(nrb);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])
% 
% fig5 = figure;
% nrbctrlplot(nrb);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODIFIED CIRCLE PARAMETRIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 20;
tau = 0;
t = linspace(0 + tau * pi/4,pi - tau * pi/4,n);
p = 1;

x1Func = @(t) R .* cos(t);
y1Func = @(t) R .* sin(t);

x2Func = @(t) R .* cos(t);
y2Func = @(t) -R .* sin(t);

x1Coord = x1Func(t);
y1Coord = y1Func(t);
x2Coord = x2Func(t);
y2Coord = y2Func(t);

%%%%%%%%%%
% Weights 
%%%%%%%%%%

z1Coord = zeros(1,n);
z2Coord = zeros(1,n);

w1Coord = ones(1,n);
w2Coord = ones(1,n);
% wCoord(ceil(n/2) - range) = sqrt(2);
% wCoord(ceil(n/2) + range) = sqrt(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create control points and knot vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 1;

knot1 = linspace(0,1,n);
knot1 = repmat(knot1,1,p)';
knot1 = knot1(:)';
knot1 = sort(knot1);

knot2 = linspace(0,1,n);
knot2 = repmat(knot2,1,p)';
knot2 = knot2(:)';
knot2 = sort(knot2);

strKnot = 0 * ones(1,p+1);
endKnot = ones(1,p+1);
knot1    = [ strKnot knot1(2:end-1) endKnot ];
knot2    = [ strKnot knot2(2:end-1) endKnot ];

pnts1 = [x1Coord ; y1Coord ; z1Coord; w1Coord];
pnts2 = [x2Coord ; y2Coord ; z2Coord; w2Coord];

crv1 = nrbmak(pnts1,knot1);
crv2 = nrbmak(pnts2,knot2);

srf = nrbruled(crv1,crv2);

figure;
nrbplot(srf,[10 10]);

%%%%%%%%%%%%%%%%%%%
% ANEURYSM STENOSIS
%%%%%%%%%%%%%%%%%%%

n = 500;
t = linspace(0,1,n);
p = 1;

R1 = 10;
R2 = 15;
R = (R2 + R1)/2;

r1 = 2;     % Amplitude of lower stenosis
r2 = -2;    % Amplitude of lower aneurysm
r3 = -1;    % Amplitude of upper stenosis
r4 = 2;     % Amplitude of upper aneurysm

k1 = (pi/180) * 5;      % Range of the lower stenosis
k2 = (pi/180) * 12;     % Range of the lower aneurysm
k3 = (pi/180) * 6;      % Range of the upper stenosis
k4 = (pi/180) * 10;     % Range of the upper aneurysm

phi = pi;   % Ending angle of artery profile

theta1 = pi/5;      % Position of the center of lower stenosis
theta2 = 3*pi/4;    % Position of the center of lower aneurysm
theta3 = pi/5;      % Position of the center of upper stenosis
theta4 = 3*pi/4;    % Position of the center of upper aneurysm

A1 = 0.01; % Amplitude of the roughness function in lower stenosis
A2 = 0.01; % Amplitude of the roughness function in lower aneurysm
A3 = 0.01; % Amplitude of the roughness function in upper stenosis
A4 = 0.01; % Amplitude of the roughness function in upper aneurysm

M1 = (pi/180) * 90; % Range of the roughness function in lower stenosis
M2 = (pi/180) * 90; % Range of the roughness function in lower aneurysm
M3 = (pi/180) * 90; % Range of the roughness function in upper stenosis
M4 = (pi/180) * 90; % Range of the roughness function in upper aneurysm

w1 = 100 * phi; % Pace of roughness peaks in lower stenosis
w2 = 50 * phi; % Pace of roughness peaks in lower aneurysm
w3 = 40 * phi; % Pace of roughness peaks in upper stenosis
w4 = 50 * phi; % Pace of roughness peaks in upper aneurysm

rghtFunc1 = @(t) (1 + A1 * cos(w1 * t) .* exp(-( (phi * t - theta1)/(M1 * k1)).^2)) .* (1 + A2 * cos(w2*t) .* exp(-( (phi * t - theta2)/(M2 * k2)).^2));
rghtFunc2 = @(t) (1 + A3 * cos(w3 * t) .* exp(-( (phi * t - theta3)/(M3 * k3)).^2)) .* (1 + A4 * cos(w4*t) .* exp(-( (phi * t - theta4)/(M4 * k4)).^2));

stn1Func = @(t) R1 + r1 * exp(-((phi * t - theta1)/k1).^2) + r2*exp(-((phi * t - theta2)/k2).^2);
stn2Func = @(t) R2 + r3 * exp(-((phi * t - theta3)/k3).^2) + r4*exp(-((phi * t - theta4)/k4).^2);

x1Func = @(t) rghtFunc1(t) .* stn1Func(t) .* cos(phi * t);
y1Func = @(t) rghtFunc1(t) .* stn1Func(t) .* sin(phi * t);

x2Func = @(t) rghtFunc2(t) .* stn2Func(t) .* cos(phi * t);
y2Func = @(t) rghtFunc2(t) .* stn2Func(t) .* sin(phi * t);

x1Coord = x1Func(t);
y1Coord = y1Func(t);
x2Coord = x2Func(t);
y2Coord = y2Func(t);

%%%%%%%%%%
% Weights 
%%%%%%%%%%

z1Coord = zeros(1,n);
z2Coord = zeros(1,n);

w1Coord = ones(1,n);
w2Coord = ones(1,n);
% wCoord(ceil(n/2) - range) = sqrt(2);
% wCoord(ceil(n/2) + range) = sqrt(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create control points and knot vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

knot1 = linspace(0,1,n);
knot1 = repmat(knot1,1,p)';
knot1 = knot1(:)';
knot1 = sort(knot1);

knot2 = linspace(0,1,n);
knot2 = repmat(knot2,1,p)';
knot2 = knot2(:)';
knot2 = sort(knot2);

strKnot = 0 * ones(1,p+1);
endKnot = ones(1,p+1);
knot1    = [ strKnot knot1(2:end-1) endKnot ];
knot2    = [ strKnot knot2(2:end-1) endKnot ];

pnts1 = [x1Coord ; y1Coord ; z1Coord; w1Coord];
pnts2 = [x2Coord ; y2Coord ; z2Coord; w2Coord];

crv1 = nrbmak(pnts1,knot1);
crv2 = nrbmak(pnts2,knot2);

srf = nrbruled(crv1,crv2);
% 
% figure;
% nrbplot(srf,[200 30]);

% %% TEST THE MAP FUNCTION EXTRACTED
% 
% % EXTRACT THE GEOMETRY INFORMATION FROM THE GENERATED VOLUME
% 
% geometry = geo_load(vol);
% 
% map = @(x,y,z) geometry.map([x;y;z]);
% Jac = @(x,y,z) geometry.map_der([x;y;z]);
% Hes = @(x,y,z) geometry.map_der2([x;y;z]);
% 
% % Coordinates defined in the refenrece domain
% 
% phi = linspace(0,1,10);
% theta = linspace(0,1,10);
% radius = linspace(0,1,10);
% 
% [RR,TT,EE] = meshgrid(phi,theta,radius);
% 
% % Map into the physical domain
% 
% XX = mapOut3D(phi,theta,radius,map,1,'Torus');
% YY = mapOut3D(phi,theta,radius,map,2,'Torus');
% ZZ = mapOut3D(phi,theta,radius,map,3,'Torus');
% 
% [m,n,o] = size(XX);
% mySize  = m * n * o;
% 
% X = reshape(XX,[1,mySize]);
% Y = reshape(YY,[1,mySize]);
% Z = reshape(ZZ,[1,mySize]);
% 
% % Plot the obtained points in the physical domain
% 
% figure;
% plot3(X,Y,Z,'.');
% 
% w = 1;
% L = 10;
% n = 100;
% 
% % Coordinates in X
% 
% xCoord = linspace(0,L,n);
% 
% % Coordinates in X
% 
% yCoord = 0 * linspace(0,L,n);
% 
% % Coordinates in X
% 
% zCoord = 0 * linspace(0,L,n);
% zCoord(ceil(n/2)) = 1;
% 
% % Weights
% 
% wCoord = ones(1,n);
% 
% knot   = linspace(0,1,n);
% knot   = [ 0 0 0 knot(2:end-1) 1 1 1 ];
% 
% pnts = [xCoord ; yCoord ; zCoord; wCoord];
% crv = nrbmak(pnts,knot);
% 
% figure
% nrbplot(crv,[100]);
%  
% srf = nrbrevolve(crv,[0 0 -3],[1 0 0]);
% 
% figure
% p = nrbeval(srf,{linspace(0.0,1.0,100) linspace(0.0,1.0,100)});
% surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
% title('Construct of a 3D surface by revolution of a curve.');
% shading interp;
% colormap(copper);
% axis equal;
% hold off

% %% SECOND TEST
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Define initial coordinates for the lower boundary 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% w  = 1;
% R  = 10;
% R2 = 20;
% r  = 2;
% n  = 100;
% p  = 1;
% range = (n/2) * (1/5);
% theta = linspace( 0 , pi/2 , n);
% 
% xCoord = R * cos(theta);
% yCoord = R * sin(theta);
% zCoord = zeros(1,n);
% 
% a = theta(ceil(n/2) - range);
% b = theta(ceil(n/2) + range);
% amp = 1 + range * 2;
% int = (b - a);
% func = @(x) exp(-( x - (a + b)/2 ).^2);
% 
% xCoord(ceil(n/2) - range : ceil(n/2) + range) = (R + r * (sin(linspace(-0.5*pi,1.5*pi,amp)) + 1)) .* cos(theta(ceil(n/2) - range : ceil(n/2) + range));
% yCoord(ceil(n/2) - range : ceil(n/2) + range) = (R + r * (sin(linspace(-0.5*pi,1.5*pi,amp)) + 1)) .* sin(theta(ceil(n/2) - range : ceil(n/2) + range));
% 
% %%%%%%%%%%
% % Weights 
% %%%%%%%%%%
% 
% wCoord = ones(1,n);
% % wCoord(ceil(n/2) - range) = sqrt(2);
% % wCoord(ceil(n/2) + range) = sqrt(2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Create control points and knot vector
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% knot = linspace(0,1,n);
% knot = repmat(knot,1,p)';
% knot = knot(:)';
% knot = sort(knot);
% 
% strKnot = 0 * ones(1,p+1);
% endKnot = ones(1,p+1);
% knot    = [ strKnot knot(2:end-1) endKnot ];
% 
% pnts = [xCoord ; yCoord ; zCoord; wCoord];
% crv1 = nrbmak(pnts,knot);
% 
% crv2 = nrbcirc(R2,[ 0 0 0 ],0,pi/2);
% 
% srf = nrbruled(crv1,crv2);
% 
% figure
% nrbplot(crv1,[100]);
% hold on;
% nrbplot(crv2,[100]);
% 
% figure
% nrbplot(srf,[100,50]);
% 
% %% TEST THE MAP FUNCTION EXTRACTED
% 
% % EXTRACT THE GEOMETRY INFORMATION FROM THE GENERATED VOLUME
% 
% geometry = geo_load(srf);
% 
% map = @(x,y) geometry.map([x;y]);
% Jac = @(x,y) geometry.map_der([x;y]);
% Hes = @(x,y) geometry.map_der2([x;y]);
% 
% % Coordinates defined in the refenrece domain
% 
% theta = linspace(0,1,10);
% radius = linspace(0,1,10);
% 
% [TT,RR] = meshgrid(theta,radius);
% 
% % Map into the physical domain
% 
% XX = mapOut(theta,radius,map,1);
% YY = mapOut(theta,radius,map,2);
% 
% S = size(XX);
% mySize  = S(1) * S(2);
% 
% X = reshape(XX,[1,mySize]);
% Y = reshape(YY,[1,mySize]);
% 
% % Plot the obtained points in the physical domain
% 
% figure;
% plot(X,Y,'.');
%  
% srf = nrbrevolve(crv,[0 0 -3],[1 0 0]);
% 
% figure
% p = nrbeval(srf,{linspace(0.0,1.0,100) linspace(0.0,1.0,100)});
% surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
% title('Construct of a 3D surface by revolution of a curve.');
% shading interp;
% colormap(copper);
% axis equal;
% hold off

%% THIRD Test

% knots = {[0 0 0 1 1 1], [0 0 0 1 1 1]};
% coefs = zeros (4, 3, 3);
% coefs(:, 1, 1) = [1, 0, 0, 1];
% coefs(:, 2, 1) = [1/sqrt(2), -1/sqrt(2), 0 , 1/sqrt(2)];
% coefs(:, 3, 1) = [0, -1, 0, 1];
% coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 2, 2) = [0, 0, 0, sqrt(2)-1];
% coefs(:, 3, 2) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 1, 3) = [0, 1, 0, 1];
% coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 3, 3) = [-1, 0, 0, 1];
% 
% pnts = [];
% pnts(1,:) = [1, 0, 0, 1];
% pnts(2,:) = [1/sqrt(2), -1/sqrt(2), 0 , 1/sqrt(2)];
% pnts(3,:) = [0, -1, 0, 1];
% pnts(4,:) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% pnts(5,:) = [0, 0, 0, sqrt(2)-1];
% pnts(6,:) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
% pnts(7,:) = [0, 1, 0, 1];
% pnts(8,:) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% pnts(9,:) = [-1, 0, 0, 1];
% 
% srf = nrbmak (coefs, knots);
% 
% figure;
% nrbplot(srf,[10,10]);
% hold on;
% plot3(pnts(:,1),pnts(:,2),pnts(:,3),'.r', 'MarkerSize' , 30)
% view([0 0 1]);
% axis auto
% daspect([1 1 1])
% pbaspect([1 1 1])
% figure.XRuler.Axle.LineStyle = 'none';  
% axis off

%% My tests
%!test
% crv = nrbdegelev (nrbcirc (1, [], 0, pi/2), 2);
% crv = nrbunclamp (crv, 3);
% xx = linspace (0, 1, 20);
% crv1 = nrbclamp (crv);
% assert (crv1.knots, [0 0 0 0 0 1 1 1 1 1])
% assert (nrbeval(crv, xx), nrbeval(crv1, xx), 1e-14)
% crv1 = nrbclamp (crv, 2);
% assert (crv1.knots, [-3 -2 -1 0 0 1 1 2 3 4])
% assert (nrbeval(crv, xx), nrbeval(crv1, xx), 1e-14)