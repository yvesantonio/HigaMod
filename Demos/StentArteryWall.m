%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STENOSIS SMOOTH HEAVISIDE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 1;
n = 500;
t = linspace(0,10,n);
p = 1;

th1 = 4;
th2 = 6;
ep = 1e-3;

theta1 = pi/5;      % Position of the beggining stenosis
theta2 = 3*pi/4;    % Position of the ending stenosis

stentStruct = @(tt) 0.5 .* (1 + (2/pi).*atan((tt - th1)./ep)) - 0.5 .* (1 + (2/pi).*atan((tt - th2)./ep));

lowerBorderX = @(tt) linspace(0,1,length(tt));
lowerBorderY = @(tt) zeros(size(tt));
upperBorderX = @(tt) linspace(0,1,length(tt));
upperBorderY = @(tt) R + stentStruct(tt);

x1Coord = lowerBorderX(t);
y1Coord = lowerBorderY(t);
x2Coord = upperBorderX(t);
y2Coord = upperBorderY(t);
z1Coord = zeros(1,n);
z2Coord = zeros(1,n);

%%%%%%%%%%
% Weights 
%%%%%%%%%%

w1Coord = ones(1,n);
w2Coord = ones(1,n);

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

figure
nrbplot(srf,[50,100]);
grid off

axis auto
axis off
% xMin = 0;
% xMax = 1;
% yMin = 0;
% yMax = 2;
% axis([xMin xMax yMin yMax])

% CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE

geometry = geo_load(srf);

% EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
% AND ITS JACOBIAN AND HESSIAN MATRICES

map = @(x,y) geometry.map([x;y]);
Jac = @(x,y) geometry.map_der([x;y]);
Hes = @(x,y) geometry.map_der2([x;y]);

geometricInfo.geometry = geometry;
geometricInfo.map = @(x,y) map(x,y);
geometricInfo.Jac = @(x,y) Jac(x,y);
geometricInfo.Hes = @(x,y) Hes(x,y);
geometricInfo.Type = 'StenosisHeaviside';