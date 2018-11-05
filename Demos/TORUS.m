%% TORUS PARAMETRIZATION
close all

%%%%%%%%%%%%%%%%%%%%%
% SURFACE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS OF THE SURFACE OF THE LUMEN WALL

P = 1;
R = 10;
r = 3;

rHat1 = 1;
rHat2 = 0;
rHat3 = 0;
rHat4 = 0.3;

w1 = 2 * pi;
w2 = 2 * pi;
w3 = 2 * pi;
w4 = 2 * pi;

phi1 = 3.0 * pi/4;
phi2 = 2.0 * pi/4;
phi3 = 2.5 * pi/4;
phi4 = 4.5 * pi/4;

theta1 = 3.0 * pi/4;
theta2 = 2.0 * pi/4;
theta3 = 1.0 * pi/4;
theta4 = 0.5 * pi/4;

km1 = .3;
km2 = .2;
km3 = .4;
km4 = .5;

kn1 = .2;
kn2 = .2;
kn3 = .2;
kn4 = .3;

% PARAMETERS OF THE SURFACE ROUGHNESS OF THE LUMEN WALL

wR1 = 25;
wR2 = 25;
wR3 = 25;
wR4 = 25;

kR1 = .5;
kR2 = .3;
kR3 = .4;
kR4 = .5;

A1 = .1;
A2 = .2;
A3 = .1;
A4 = .1;

phiR1 = 2 * pi;
phiR2 = 2 * pi;
phiR3 = 2 * pi;
phiR4 = 2 * pi;

WR = .8;
AR = 0.1;

N = 100;

%%%%%%%%%%%%%%%%%%
% SURFACE EQUATION
%%%%%%%%%%%%%%%%%%

% ANEURYSMS EQUATIONS

funcAne1 = @(s,t) (1 + rHat1 * exp(-((w1 * s - phi1)/km1).^2) .* exp(-((w1 * t - theta1)/kn1).^2));
funcAne2 = @(s,t) (1 + rHat2 * exp(-((w2 * s - phi2)/km2).^2) .* exp(-((w2 * t - theta2)/kn2).^2));

% STENOSIS EQUATIONS

funcSte1 = @(s,t) (1 - rHat3 * exp(-((w3 * s - phi3)/km3).^2) .* exp(-((w3 * t - theta3)/kn3).^2));
funcSte2 = @(s,t) (1 - rHat4 * exp(-((w4 * s - phi4)/km4).^2) .* exp(-((w4 * t - theta4)/kn4).^2));

% ANEURYSIM ROUGHNESS EQUATIONS

aux1S = @(q) cos((wR1 * pi * q).^1) .* exp(-( (phiR1 * q - phi1)/kR1).^2);
aux1T = @(q) cos((wR1 * pi * q).^1) .* exp(-( (phiR1 * q - theta1)/kR1).^2);
funcRoughAne1 = @(s,t) (1 + A1 * aux1S(s) * aux1T(t));

aux2S = @(q) cos((wR2 * pi * q).^1) .* exp(-( (phiR2 * q - phi2)/kR2).^2);
aux2T = @(q) cos((wR2 * pi * q).^1) .* exp(-( (phiR2 * q - theta2)/kR2).^2);
funcRoughAne2 = @(s,t) (1 + A2 * aux2S(s) * aux2T(t));

% STENOSIS ROUGHNESS EQUATIONS

aux3S = @(q) cos((wR3 * pi * q).^1) .* exp(-( (phiR3 * q - phi3)/kR3).^2);
aux3T = @(q) cos((wR3 * pi * q).^1) .* exp(-( (phiR3 * q - theta3)/kR3).^2);
funcRoughSte1 = @(s,t) (1 + A3 * aux3S(s) * aux3T(t));

aux4S = @(q) cos((wR4 * pi * q).^1) .* exp(-( (phiR4 * q - phi4)/kR4).^2);
aux4T = @(q) cos((wR4 * pi * q).^1) .* exp(-( (phiR4 * q - theta4)/kR4).^2);
funcRoughSte2 = @(s,t) (1 + A4 * aux4S(s) * aux4T(t));

% COORDINATES OF THE POINTS IN THE SURFACE

rNew = @(s,t) r .* funcAne1(s,t) .* funcAne2(s,t) .* funcSte1(s,t) .* funcSte2(s,t);
rNew = @(s,t) rNew(s,t) .* funcRoughSte1(s,t) .* funcRoughSte2(s,t) .* funcRoughAne1(s,t) .* funcRoughAne2(s,t);
RNew = @(s,t) R;

xCoord = @(m,n) (RNew(m,n) + rNew(m,n) .* cos(2 * pi * m)) .* cos(2 * pi * n);
yCoord = @(m,n) (RNew(m,n) + rNew(m,n) .* cos(2 * pi * m)) .* sin(2 * pi * n);
zCoord = @(m,n) rNew(m,n) .* sin(2 * pi * m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTERNAL SURFACE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS OF THE EXTERNAL SURFACE OF THE ARTERIAL WALL

P = 1;
R = 10;
r = 3 + 1;

rHat1 = 1 - .25;
rHat2 = 0;
rHat3 = 0;
rHat4 = 0;

w1 = 2 * pi;
w2 = 2 * pi;
w3 = 2 * pi;
w4 = 2 * pi;

phi1 = 3.0 * pi/4;
phi2 = 2.0 * pi/4;
phi3 = 2.5 * pi/4;
phi4 = 4.5 * pi/4;

theta1 = 3.0 * pi/4;
theta2 = 2.0 * pi/4;
theta3 = 1.0 * pi/4;
theta4 = 0.5 * pi/4;

km1 = .3;
km2 = .2;
km3 = .4;
km4 = .5;

kn1 = .2;
kn2 = .2;
kn3 = .2;
kn4 = .3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTERNAL SURFACE EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ANEURYSMS EQUATIONS

funcAne1 = @(s,t) (1 + rHat1 * exp(-((w1 * s - phi1)/km1).^2) .* exp(-((w1 * t - theta1)/kn1).^2));
funcAne2 = @(s,t) (1 + rHat2 * exp(-((w2 * s - phi2)/km2).^2) .* exp(-((w2 * t - theta2)/kn2).^2));

% STENOSIS EQUATIONS

funcSte1 = @(s,t) (1 - rHat3 * exp(-((w3 * s - phi3)/km3).^2) .* exp(-((w3 * t - theta3)/kn3).^2));
funcSte2 = @(s,t) (1 - rHat4 * exp(-((w4 * s - phi4)/km4).^2) .* exp(-((w4 * t - theta4)/kn4).^2));

% COORDINATES OF THE POINTS IN THE SURFACE

rNewExt = @(s,t) r .* funcAne1(s,t) .* funcAne2(s,t) .* funcSte1(s,t) .* funcSte2(s,t);
RNewExt = @(s,t) R;

xCoordExt = @(m,n) (RNewExt(m,n) + rNewExt(m,n) .* cos(2 * pi * m)) .* cos(2 * pi * n);
yCoordExt = @(m,n) (RNewExt(m,n) + rNewExt(m,n) .* cos(2 * pi * m)) .* sin(2 * pi * n);
zCoordExt = @(m,n) rNewExt(m,n) .* sin(2 * pi * m);

%%%%%%%%%%%%%
% FIRST CURVE
%%%%%%%%%%%%%

m = linspace(0,0.25,N); % Internal radius
n = linspace(0,0.5,N); % External radius
k = linspace(0,1,100);
pnts1 = zeros(3,N,N);

for ii = 1:N
    pnts1(:,:,ii) = [ xCoord(m,n(ii)) ; yCoord(m,n(ii)) ; zCoord(m,n(ii)) ];
end

knot1 = linspace(0,1,N);
knot1 = repmat(knot1,1,P)';
knot1 = knot1(:)';
knot1 = sort(knot1);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot1    = [ strKnot knot1(2:end-1) endKnot ];

%%%%%%%%%%%%%%
% SECOND CURVE
%%%%%%%%%%%%%%

m = linspace(0.25,0.5,N); % Internal radius
n = linspace(0,0.5,N); % External radius
k = linspace(0,1,100);
pnts2 = zeros(3,N,N);

for ii = 1:N
    pnts2(:,:,ii) = [ xCoord(m,n(ii)) ; yCoord(m,n(ii)) ; zCoord(m,n(ii)) ];
end

knot2 = linspace(0,1,N);
knot2 = repmat(knot2,1,P)';
knot2 = knot2(:)';
knot2 = sort(knot2);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot2    = [ strKnot knot2(2:end-1) endKnot ];

%%%%%%%%%%%%%
% THIRD CURVE
%%%%%%%%%%%%%

m = linspace(0.5,0.75,N); % Internal radius
n = linspace(0,0.5,N); % External radius
k = linspace(0,1,100);
pnts3 = zeros(3,N,N);

for ii = 1:N
    pnt3 = [ xCoord(m,n(ii)) ; yCoord(m,n(ii)) ; zCoord(m,n(ii)) ];
    pnt3rev(1,:) = wrev(pnt3(1,:));
    pnt3rev(2,:) = wrev(pnt3(2,:));
    pnt3rev(3,:) = wrev(pnt3(3,:));
    pnts3(:,:,ii) = pnt3rev;
end

knot3 = linspace(0,1,N);
knot3 = repmat(knot3,1,P)';
knot3 = knot3(:)';
knot3 = sort(knot3);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot3    = [ strKnot knot3(2:end-1) endKnot ];

%%%%%%%%%%%%%%
% FOURTH CURVE
%%%%%%%%%%%%%%

m = linspace(0.75,1,N); % Internal radius
n = linspace(0,0.5,N); % External radius
k = linspace(0,1,100);
pnts4 = zeros(3,N,N);

for ii = 1:N
    pnt4 = [ xCoord(m,n(ii)) ; yCoord(m,n(ii)) ; zCoord(m,n(ii)) ];
    pnt4rev(1,:) = wrev(pnt4(1,:));
    pnt4rev(2,:) = wrev(pnt4(2,:));
    pnt4rev(3,:) = wrev(pnt4(3,:));
    pnts4(:,:,ii) = pnt4rev;
end

knot4 = linspace(0,1,N);
knot4 = repmat(knot4,1,P)';
knot4 = knot4(:)';
knot4 = sort(knot4);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot4    = [ strKnot knot4(2:end-1) endKnot ];

%%%%%%%%%%%%%%%%%%
% DESIGN SECTIONS
%%%%%%%%%%%%%%%%%%

srf = cell(N);
crv1 = cell(N);
crv2 = cell(N);
crv3 = cell(N);
crv4 = cell(N);
for ii = 1:N
    crv1{ii} = nrbmak(pnts1(:,:,ii),knot1);
    crv2{ii} = nrbmak(pnts2(:,:,ii),knot2);
    crv3{ii} = nrbmak(pnts3(:,:,ii),knot3);
    crv4{ii} = nrbmak(pnts4(:,:,ii),knot4);
    srf{ii} = nrbcoons(crv1{ii}, crv3{ii}, crv4{ii}, crv2{ii});
end

%%%%%%%%%%%%%%%%
% CREATE VOLUME
%%%%%%%%%%%%%%%%

% Extract the control points 

volCoefs = zeros(4,N,N,N);
for ii = 1:N
    cut = srf{ii};
    volCoefs(:,:,:,ii) = cut.coefs;
end

% Create knots for the third parametric component

knots = linspace(0,1,N);
knots = repmat(knots,1,P)';
knots = knots(:)';
knots = sort(knots);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knots    = [ strKnot knots(2:end-1) endKnot ];

% Create the volume

VOL = nrbmak(volCoefs,{knots,knots,knots});

%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%

figure;
nrbplot(VOL,[50 50 70]);
az = 0;
el = 40;
view(az, el);
ax = gca;
ax.Visible = 'off';
export_fig(sprintf([ 'Lumen_Volume' ]),'-pdf');

figure;
nrbplot(srf{50},[50 50]);hold on;
nrbplot(srf{60},[50 50]);
nrbplot(srf{70},[50 50]);
nrbplot(srf{80},[50 50]);
nrbplot(srf{90},[50 50]);
nrbplot(srf{100},[50 50]);
hold off;
az = -20;
el = 40;
view(az, el);
ax = gca;
ax.Visible = 'off';
export_fig(sprintf([ 'Lumen_Aneurysm' ]),'-pdf');

figure;
nrbplot(srf{1},[50 50]);hold on;
nrbplot(srf{5},[50 50]);
nrbplot(srf{10},[50 50]);
nrbplot(srf{15},[50 50]);
nrbplot(srf{20},[50 50]);
nrbplot(srf{25},[50 50]);
hold off;
az = -20;
el = 40;
view(az, el);
ax = gca;
ax.Visible = 'off';
export_fig(sprintf([ 'Lumen_Stenosis' ]),'-pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: THIS METHOD DOES NOT WORK BECAUSE THE BECAUSE THE RULED FUNCTION
% DOES NOT WORK FOR SURFACES, ONLY FOR CURVES. THIS WAY, WE HAVE TO CREATE
% EACH INDIVIDUAL CUT AND GLUE THEM TOEGHETER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE ARTERIAL WALL VOLUME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% FIRST CURVE
%%%%%%%%%%%%%%

m = linspace(0,1,N); % Internal radius
n = linspace(0,0.5,N); % External radius
k = linspace(0,1,100);
pnts1 = zeros(3,N,N);

for ii = 1:N
    pnts1(:,:,ii) = [ xCoord(m,n(ii)) ; yCoord(m,n(ii)) ; zCoord(m,n(ii)) ];
end

knot1 = linspace(0,1,N);
knot1 = repmat(knot1,1,P)';
knot1 = knot1(:)';
knot1 = sort(knot1);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot1    = [ strKnot knot1(2:end-1) endKnot ];

%%%%%%%%%%%%%%%
% SECOND CURVE
%%%%%%%%%%%%%%%

m = linspace(0,1,N); % Internal radius
n = linspace(0,0.5,N); % External radius
k = linspace(0,1,100);
pnts2 = zeros(3,N,N);

for ii = 1:N
    pnts2(:,:,ii) = [ xCoordExt(m,n(ii)) ; yCoordExt(m,n(ii)) ; zCoordExt(m,n(ii)) ];
end

knot2 = linspace(0,1,N);
knot2 = repmat(knot2,1,P)';
knot2 = knot2(:)';
knot2 = sort(knot2);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot2    = [ strKnot knot2(2:end-1) endKnot ];

%%%%%%%%%%%%%%%%%%
% DESIGN SECTIONS
%%%%%%%%%%%%%%%%%%

srf = cell(N);
crv1 = cell(N);
crv2 = cell(N);
for ii = 1:N
    crv1{ii} = nrbmak(pnts1(:,:,ii),knot1);
    crv2{ii} = nrbmak(pnts2(:,:,ii),knot2);
    srf{ii} = nrbruled(crv2{ii}, crv1{ii});
end

%%%%%%%%%%%%%%%%
% CREATE VOLUME
%%%%%%%%%%%%%%%%

% Extract the control points 

volCoefs = zeros(4,N,2,N);
for ii = 1:N
    cut = srf{ii};
    volCoefs(:,:,:,ii) = cut.coefs;
end

% Create knots for the third parametric component

knots = linspace(0,1,N);
knots = repmat(knots,1,P)';
knots = knots(:)';
knots = sort(knots);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knots    = [ strKnot knots(2:end-1) endKnot ];

% Create the volume

VOL = nrbmak(volCoefs,{knots,[0 0 1 1],knots});

%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%

figure;
nrbplot(VOL,[100 10 100]);
az = 0;
el = 40;
view(az, el);
ax = gca;
ax.Visible = 'off';
export_fig(sprintf([ 'Wall_Volume' ]),'-pdf');

figure;
nrbplot(srf{70},[100 10]);hold on;
nrbplot(srf{75},[100 10]);
nrbplot(srf{80},[100 10]);
nrbplot(srf{85},[100 10]);
nrbplot(srf{90},[100 10]);
nrbplot(srf{95},[100 10]);
hold off;
az = -20;
el = 40;
view(az, el);
ax = gca;
ax.Visible = 'off';
export_fig(sprintf([ 'Wall_Stenosis' ]),'-pdf');

figure;
nrbplot(srf{1},[100 10]);hold on;
nrbplot(srf{5},[100 10]);
nrbplot(srf{10},[100 10]);
nrbplot(srf{15},[100 10]);
nrbplot(srf{20},[100 10]);
nrbplot(srf{25},[100 10]);
hold off;
az = -20;
el = 40;
view(az, el);
ax = gca;
ax.Visible = 'off';
export_fig(sprintf([ 'Wall_Aneurysm' ]),'-pdf');

coeffVec = [1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];
for ii = 1:length(coeffVec)
    figure;
    nrbplot(srf{coeffVec(ii)},[150 10]);hold on;
    hold off;
    az = 180 * coeffVec(ii)/100;
    el = 0;
    view(az, el);
    ax = gca;
    ax.Visible = 'off';
    export_fig(sprintf([ 'Wall_Section_' num2str(az) '_Degrees' ]),'-pdf');
end