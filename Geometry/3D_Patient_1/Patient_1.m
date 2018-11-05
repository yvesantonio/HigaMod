close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANAGE PATIENT GEOMETRY
% INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ .CSV DATA FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = [pwd,'\morphology\centerlines.csv'];
DATA = csvread(filename);

minX = min(DATA(:,1));
minY = min(DATA(:,2));
minZ = min(DATA(:,3));
minR = min(DATA(:,4));

maxX = max(DATA(:,1));
maxY = max(DATA(:,2));
maxZ = max(DATA(:,3));
maxR = max(DATA(:,4));

disp('FINISHED READING DATA')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE JUMPS IN CENTERLINE INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the jumps in the dataset. If the jump is gretater than 0.1 mm,
% the precision of the centerline data, then we have another centerline.

len = size(DATA,1);

X     = DATA(1:end-1,1:3);
Xnext = DATA(2:end,1:3);
jump  = (Xnext - X).^2;
jump  = sum(jump,2);
jump  = jump.^(1/2);

I = find(jump > 0.5);

disp('FINISHED DETECTING CENTERLINES')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT CENTERLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The dataset contains a list of points of the various centerlines in the
% arterial network. We first have to extract each individual segment, to
% separate the simulation domains.

lineStruct = cell(size(I,1)+1,1);

for ii = 1:size(I,1)+1
    if(ii == 1)
        lineStruct{1} = DATA(1:I(1),:);
    elseif(ii == size(I,1)+1)
        lineStruct{end} = DATA(I(end) + 2:end,:);
    else
        lineStruct{ii} = DATA(I(ii-1)+2:I(ii),:);
    end
end

disp('FINISHED EXTRACTING CENTERLINES')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT OBTAINED CENTERLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
for ii = 1:size(lineStruct,1)
    plot3(lineStruct{ii}(:,1),lineStruct{ii}(:,2),lineStruct{ii}(:,3),'.');
    hold on
end
grid on

disp('FINISHED PLOTTING CENTERLINES')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATE CENTERLINE DATA SCATTERED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compute the intersection points and its plane properties in each
% bifurcation of the arterial network.

interpStruct = cell(size(lineStruct,1),4);

for ii = 1:size(lineStruct,1)
    
    param = linspace(0,1,length(lineStruct{ii}(:,1)));
    
    interpX = griddedInterpolant(param,lineStruct{ii}(:,1));
    interpY = griddedInterpolant(param,lineStruct{ii}(:,2));
    interpZ = griddedInterpolant(param,lineStruct{ii}(:,3));
    interpR = griddedInterpolant(param,lineStruct{ii}(:,4));
    
    interpStruct{ii,1} = interpX;
    interpStruct{ii,2} = interpY;
    interpStruct{ii,3} = interpZ;
    interpStruct{ii,4} = interpR;
    
end

disp('FINISHED INTERPOLATING CENTERLINE DATA')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE UNIFORM SET OF PARAMETRIC POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NSec = 2e1;
N = 2e2;
P = 1;

param = linspace(0,1,N);
    
xCoord = interpStruct{ii,1}(param);
yCoord = interpStruct{ii,2}(param);
zCoord = interpStruct{ii,3}(param);
rCoord = interpStruct{ii,4}(param);

disp('FINISHED EVALUATING CENTERLINE PARAMETRIC CURVES')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATE THE POINTS IN THE CENTERLINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

points = [xCoord,yCoord,zCoord];
points = points';
Cline = cscvn(points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE TANGENT, NORMAL AND BINORMAL
% VECTORS OF THE CENTERLINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Tn, Nm, Bm] = myfrenet(xCoord,yCoord,zCoord);
total = size(Tn,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THE CROSS SECTION EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

funcStruct = cell(1,total);
for ii = 1:total
    
    % CENTER OF THE CIRCLE
    
    Pc = [xCoord(ii),yCoord(ii),zCoord(ii)];
    
    % RADIUS OF THE CIRCLE
    
    R = rCoord(ii);
    
    % PARAMETRIZATION SUPPORT VECTORS
    
    uVect = Nm(ii,:);
    vVect = Bm(ii,:);
    
    % CIRCLE PARAMETRIC EQUATION
    
    funcStruct{ii} = @(t) (R .* cos(t) .* uVect) + (R .* sin(t) .* vVect) + Pc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFACE EQUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
% FIRST CURVE
%%%%%%%%%%%%%

theta = linspace(0,2*pi*0.25,NSec);
pnts1 = zeros(3,NSec,N);

for jj = 1:N
    
    % DEFINE THE CURRENT ANALYSED CROSS SECTION
    
    myFunc = funcStruct{jj};
    ptsFunc = zeros(NSec,3);
    
    % COMPUTE THE CROSS SECTION BOUNDARY
    
    for ii = 1:NSec
        ptsFunc(ii,:) = myFunc(theta(ii));
    end
    
    % STORE THE DATA FROM EACH CROSS SECTION
    
    pnts1(:,:,jj) = ptsFunc';
    
end

knot1 = linspace(0,1,NSec);
knot1 = repmat(knot1,1,P)';
knot1 = knot1(:)';
knot1 = sort(knot1);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot1    = [ strKnot knot1(2:end-1) endKnot ];

%%%%%%%%%%%%%%
% SECOND CURVE
%%%%%%%%%%%%%%

theta = linspace(2*pi*0.25,2*pi*0.5,NSec);
pnts2 = zeros(3,NSec,N);

for jj = 1:N
    
    % DEFINE THE CURRENT ANALYSED CROSS SECTION
    
    myFunc = funcStruct{jj};
    ptsFunc = zeros(NSec,3);
    
    % COMPUTE THE CROSS SECTION BOUNDARY
    
    for ii = 1:NSec
        ptsFunc(ii,:) = myFunc(theta(ii));
    end
    
    % STORE THE DATA FROM EACH CROSS SECTION
    
    pnts2(:,:,jj) = ptsFunc';
    
end

knot2 = linspace(0,1,NSec);
knot2 = repmat(knot2,1,P)';
knot2 = knot2(:)';
knot2 = sort(knot2);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot2    = [ strKnot knot2(2:end-1) endKnot ];

%%%%%%%%%%%%%
% THIRD CURVE
%%%%%%%%%%%%%

theta = linspace(2*pi*0.5,2*pi*0.75,NSec);
pnts3 = zeros(3,NSec,N);
pnt3rev = zeros(3,NSec);

for jj = 1:N
    
    % DEFINE THE CURRENT ANALYSED CROSS SECTION
    
    myFunc = funcStruct{jj};
    ptsFunc = zeros(NSec,3);
    
    % COMPUTE THE CROSS SECTION BOUNDARY
    
    for ii = 1:NSec
        ptsFunc(ii,:) = myFunc(theta(ii));
    end
    
    % REVERSE THE DIRECTION OF THE POINTS
    
    pnt3 = ptsFunc';
    pnt3rev(1,:) = wrev(pnt3(1,:));
    pnt3rev(2,:) = wrev(pnt3(2,:));
    pnt3rev(3,:) = wrev(pnt3(3,:));
    
    % STORE THE DATA FROM EACH CROSS SECTION
    
    pnts3(:,:,jj) = pnt3rev;
    
end

knot3 = linspace(0,1,NSec);
knot3 = repmat(knot3,1,P)';
knot3 = knot3(:)';
knot3 = sort(knot3);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot3    = [ strKnot knot3(2:end-1) endKnot ];

%%%%%%%%%%%%%%
% FOURTH CURVE
%%%%%%%%%%%%%%

theta = linspace(2*pi*0.75,2*pi,NSec);
pnts4 = zeros(3,NSec,N);
pnt4rev = zeros(3,NSec);

for jj = 1:N
    
    % DEFINE THE CURRENT ANALYSED CROSS SECTION
    
    myFunc = funcStruct{jj};
    ptsFunc = zeros(NSec,3);
    
    % COMPUTE THE CROSS SECTION BOUNDARY
    
    for ii = 1:NSec
        ptsFunc(ii,:) = myFunc(theta(ii));
    end
    
    % REVERSE THE DIRECTION OF THE POINTS
    
    pnt4 = ptsFunc';
    pnt4rev(1,:) = wrev(pnt4(1,:));
    pnt4rev(2,:) = wrev(pnt4(2,:));
    pnt4rev(3,:) = wrev(pnt4(3,:));
    
    % STORE THE DATA FROM EACH CROSS SECTION
    
    pnts4(:,:,jj) = pnt4rev;
    
end

knot4 = linspace(0,1,NSec);
knot4 = repmat(knot4,1,P)';
knot4 = knot4(:)';
knot4 = sort(knot4);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knot4    = [ strKnot knot4(2:end-1) endKnot ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESIGN SECTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

disp('FINISHED CREATING VOLUME SECTIONS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE VOLUME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the control points 

volCoefs = zeros(4,NSec,NSec,N);
for ii = 1:N
    cut = srf{ii};
    volCoefs(:,:,:,ii) = cut.coefs;
end

% Create knots for the third parametric component

knotsSurf = linspace(0,1,NSec);
knotsSurf = repmat(knotsSurf,1,P)';
knotsSurf = knotsSurf(:)';
knotsSurf = sort(knotsSurf);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knotsSurf    = [ strKnot knotsSurf(2:end-1) endKnot ];

knotsVol = linspace(0,1,N);
knotsVol = repmat(knotsVol,1,P)';
knotsVol = knotsVol(:)';
knotsVol = sort(knotsVol);
strKnot = 0 * ones(1,P+1);
endKnot = ones(1,P+1);
knotsVol    = [ strKnot knotsVol(2:end-1) endKnot ];

% Create the volume

VOL = nrbmak(volCoefs,{knotsSurf,knotsSurf,knotsVol});

disp('FINISHED CREATING VOLUME')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geometry = geo_load(VOL);

% EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
% AND ITS JACOBIAN AND HESSIAN MATRICES

map = @(x,y,z) geometry.map([x;y;z]);
Jac = @(x,y,z) geometry.map_der([x;y;z]);
Hes = @(x,y,z) geometry.map_der2([x;y;z]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE GEOMETRIC INFO STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/
% geometricInfo.geometry = geometry;
% geometricInfo.Type = 'Cylinder';
% geometricInfo.map = @(x,y,z) map(x,y,z);
% geometricInfo.Jac = @(x,y,z) Jac(x,y,z);
% geometricInfo.Hes = @(x,y,z) Hes(x,y,z);

% type = geometricInfo.Type;

% Nx = 50;
% Ny = 50;
% Nz = 200;

% evalX = linspace(0,1,Nx);
% evalY = linspace(0,1,Ny);
% evalZ = linspace(0,1,Nz);

% tic;
% [X,Y,Z] = mapOut3D(evalX,evalY,evalZ,map,type);
% [evalJac,structPhi,evalDetJac] = jacOut3D(evalX,evalY,evalZ,Jac,type);
% t1 = toc;

% [Nx,Ny,Nz] = size(evalJac);
% evalJ11 = zeros(Nx,Ny,Nz);
% evalJ12 = zeros(Nx,Ny,Nz);
% evalJ13 = zeros(Nx,Ny,Nz);
% evalJ21 = zeros(Nx,Ny,Nz);
% evalJ22 = zeros(Nx,Ny,Nz);
% evalJ23 = zeros(Nx,Ny,Nz);
% evalJ31 = zeros(Nx,Ny,Nz);
% evalJ32 = zeros(Nx,Ny,Nz);
% evalJ33 = zeros(Nx,Ny,Nz);

% for ii = 1:Nx
%     for jj = 1:Ny
%         for kk = 1:Nz
%             evalJ11(ii,jj,kk) = evalJac{ii,jj,kk}(1,1);
%             evalJ12(ii,jj,kk) = evalJac{ii,jj,kk}(1,2);
%             evalJ13(ii,jj,kk) = evalJac{ii,jj,kk}(1,3);
%             evalJ21(ii,jj,kk) = evalJac{ii,jj,kk}(2,1);
%             evalJ22(ii,jj,kk) = evalJac{ii,jj,kk}(2,2);
%             evalJ23(ii,jj,kk) = evalJac{ii,jj,kk}(2,3);
%             evalJ31(ii,jj,kk) = evalJac{ii,jj,kk}(3,1);
%             evalJ32(ii,jj,kk) = evalJac{ii,jj,kk}(3,2);
%             evalJ33(ii,jj,kk) = evalJac{ii,jj,kk}(3,3);
%         end
%     end
% end

% disp('FINISHED BASE MAPPING')

% evalX = linspace(0,1,Nx);
% evalY = linspace(0,1,Ny);
% evalZ = linspace(0,1,Nz);
% [XBase,YBase,ZBase] = ndgrid(evalX,evalY,evalZ);

% tic;
% Phi1 = griddedInterpolant(XBase,YBase,ZBase,X);
% Phi2 = griddedInterpolant(XBase,YBase,ZBase,Y);
% Phi3 = griddedInterpolant(XBase,YBase,ZBase,Z);

% J11 = griddedInterpolant(XBase,YBase,ZBase,evalJ11);
% J12 = griddedInterpolant(XBase,YBase,ZBase,evalJ12);
% J13 = griddedInterpolant(XBase,YBase,ZBase,evalJ13);
% J21 = griddedInterpolant(XBase,YBase,ZBase,evalJ21);
% J22 = griddedInterpolant(XBase,YBase,ZBase,evalJ22);
% J23 = griddedInterpolant(XBase,YBase,ZBase,evalJ23);
% J31 = griddedInterpolant(XBase,YBase,ZBase,evalJ31);
% J32 = griddedInterpolant(XBase,YBase,ZBase,evalJ32);
% J33 = griddedInterpolant(XBase,YBase,ZBase,evalJ33);

% detJ = griddedInterpolant(XBase,YBase,ZBase,evalDetJac);

% intrpPhi1_dx = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi1_dx);
% intrpPhi1_dy = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi1_dy);
% intrpPhi1_dz = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi1_dz);
% intrpPhi2_dx = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi2_dx);
% intrpPhi2_dy = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi2_dy);
% intrpPhi2_dz = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi2_dz);
% intrpPhi3_dx = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi3_dx);
% intrpPhi3_dy = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi3_dy);
% intrpPhi3_dz = griddedInterpolant(XBase,YBase,ZBase,structPhi.Phi3_dz);

% t2 = toc;

% geometricInfo.Phi1 = Phi1;
% geometricInfo.Phi2 = Phi2;
% geometricInfo.Phi3 = Phi3;
% geometricInfo.J11  = J11;
% geometricInfo.J12  = J12;
% geometricInfo.J13  = J13;
% geometricInfo.J21  = J21;
% geometricInfo.J22  = J22;
% geometricInfo.J23  = J23;
% geometricInfo.J31  = J31;
% geometricInfo.J32  = J32;
% geometricInfo.J33  = J33;
% geometricInfo.detJ  = detJ;

% geometricInfo.Phi1_dx = intrpPhi1_dx;
% geometricInfo.Phi1_dy = intrpPhi1_dy;
% geometricInfo.Phi1_dz = intrpPhi1_dz;
% geometricInfo.Phi2_dx = intrpPhi2_dx;
% geometricInfo.Phi2_dy = intrpPhi2_dy;
% geometricInfo.Phi2_dz = intrpPhi2_dz;
% geometricInfo.Phi3_dx = intrpPhi3_dx;
% geometricInfo.Phi3_dy = intrpPhi3_dy;
% geometricInfo.Phi3_dz = intrpPhi3_dz;


disp('FINISHED CREATING THE MAP INTERPOLANT')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SIMULATION TIMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['Number of D.O.F.s    : ',num2str(Nx * Ny * Nz)])
% disp(['Exec. Time Base Map  : ',num2str(t1),' [sec]'])
% disp(['Exec. Time Interpol  : ',num2str(t2),' [sec]'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT GEOMETRY DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename1 = '1Patient3DGeo.mat';
% geo = VOL;
% save(filename1,'geo')
% 
% disp('FINISHED EXPORTING GEOMETRY')
% 
% filename2 = '1Patient3DMap.mat';
% save(filename2,'geometricInfo','-v7.3')
% 
% disp('FINISHED EXPORTING MAPPING')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE NURBS VOLUME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
nrbplot(VOL,[10 10 500]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE BOUNDARY OF A CROSS SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% nrbplot(crv1{10},20); hold on;
% nrbplot(crv2{10},20);
% nrbplot(crv3{10},20);
% nrbplot(crv4{10},20); hold off;

% PLOT THE CROSS SECTION

% figure;
% nrbplot(srf{10},[20 20]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE .STL FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE UNIFORM SET OF PARAMETRIC POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interpStruct = cell(size(lineStruct,1),4);

for ii = 1:size(lineStruct,1)
    
    param = linspace(0,1,length(lineStruct{ii}(:,1)));
    
    interpX = griddedInterpolant(param,lineStruct{ii}(:,1));
    interpY = griddedInterpolant(param,lineStruct{ii}(:,2));
    interpZ = griddedInterpolant(param,lineStruct{ii}(:,3));
    interpR = griddedInterpolant(param,lineStruct{ii}(:,4));
    
    interpStruct{ii,1} = interpX;
    interpStruct{ii,2} = interpY;
    interpStruct{ii,3} = interpZ;
    interpStruct{ii,4} = interpR;
    
end

NSec = NSec * 4;
N = 2e2;
P = 1;

param = linspace(0,1,N);

xCoord = interpStruct{ii,1}(param);
yCoord = interpStruct{ii,2}(param);
zCoord = interpStruct{ii,3}(param);
rCoord = interpStruct{ii,4}(param);

points = [xCoord,yCoord,zCoord];
points = points';
Cline = cscvn(points);

[Tn, Nm, Bm] = myfrenet(xCoord,yCoord,zCoord);
total = size(Tn,1);

funcStruct = cell(1,total);
for ii = 1:total
    
    % CENTER OF THE CIRCLE
    
    Pc = [xCoord(ii),yCoord(ii),zCoord(ii)];
    
    % RADIUS OF THE CIRCLE
    
    R = rCoord(ii);
    
    % PARAMETRIZATION SUPPORT VECTORS
    
    uVect = Nm(ii,:);
    vVect = Bm(ii,:);
    
    % CIRCLE PARAMETRIC EQUATION
    
    funcStruct{ii} = @(t) (R .* cos(t) .* uVect) + (R .* sin(t) .* vVect) + Pc;
end

% CREATE .STL FILE FOR THE SURFACE

NSecStl = NSec;

theta = linspace(0,2*pi,NSecStl);
Ntotal = N*NSecStl;
srfPts = zeros(Ntotal,3);

for jj = 1:N
    
    % DEFINE THE CURRENT ANALYSED CROSS SECTION
    
    myFunc = funcStruct{jj};
    ptsFunc = zeros(NSec,3);
    
    % COMPUTE THE CROSS SECTION BOUNDARY
    
    for ii = 1:NSecStl
        ptsFunc(ii,:) = myFunc(theta(ii));
    end
    
    srfPts(1 + (jj - 1)*NSecStl:jj*NSecStl,:) = ptsFunc;
    
end

% aux = srfPts;
% srfPts = zeros(Ntotal+2,3);
% 
% srfPts(1,:) = [xCoord(1),yCoord(1),zCoord(1)];
% srfPts(2:end-1,:) = aux;
% srfPts(end,:) = [xCoord(end),yCoord(end),zCoord(end)];

FILE = 'srfPts.txt';
dlmwrite(FILE,srfPts,'delimiter','\t','precision',5)

% figure;
% scatter3(srfPts(:,1),srfPts(:,2),srfPts(:,3),'.');

p = srfPts;
% [t,tnorm]=MyRobustCrust(p);
t=MyCrustOpen(p);

% figure(1);
% set(gcf,'position',[0,0,1280,800]);
% subplot(1,2,1)
% hold on
% axis equal
% title('Points Cloud','fontsize',14)
% plot3(p(:,1),p(:,2),p(:,3),'g.')
% view(3);
% axis vis3d

% figure(1)
% subplot(1,2,2)
% hold on
% title('Output Triangulation','fontsize',14)
% axis equal
% trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
% view(3);
% axis vis3d

stlwrite('3D_Patient_1_PRE.stl',t,[p(:,1),p(:,2),p(:,3)]);

model = createpde(1);
importGeometry(model,'3D_Patient_1_PRE.stl');
mesh = generateMesh(model,'Hmax',.5,'GeometricOrder','linear');
[p,e,t] = meshToPet(mesh);

TR = triangulation(t(1:end-1,:)',p');
[F,Pbound] = freeBoundary(TR);
% figure
% trisurf(F,Pbound(:,1),Pbound(:,2),Pbound(:,3), ...
%        'FaceColor','cyan','FaceAlpha',0.8);

prefix = '3D_Patient_1';
FILE = [prefix,'.stl'];

stlwrite(FILE,F,Pbound);
% figure;
% pdeplot3D(model);

disp('FINISHED CREATING .STL FILE FOR THE SURFACE')
