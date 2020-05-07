%-+--------------------------------------------------------------------
% HiModlab is a general purpose Hierarchical model reduction library.
% Copyright (C) 2006-2017  by the HiMod authors (see authors.txt).
%

% This file is part of the HiModlab library.
%
% The HiModlab library is free software: you can use it, redistribute
% it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either
% version 2 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Acknowledgement and Disclaimer
%
% The development of HiModLab was supported by the US National Science
% Foundation (NSF Project, DMS 1419060 ). Any opinions, findings and conclusions
% or recomendations expressed in the source code and documentation are
% those of the authors and do not necessarily reflect the views of the
% National Science Foundation (NSF).
%
%
%
%-+--------------------------------------------------------------------
%% HI-MOD MONODOMAIN
% The following script allows the solution of an Advection - Diffusion -
% Reaction differential problem in 2D using the HIGAMod solution.

    clear
    close all

    disp('******************************************')
    disp('*           HIGAMod Simulation           *');
    disp('******************************************');

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
MU = 1e-5 * [0.5 0.25];
ww = [1 2 4];

strPlotStruct = cell(length(MU),length(ww));

for iii = 1:length(MU)
    for jjj = 1:length(ww)
        
    %% Simulation case   
    caso  = 9;      % Analysed Case
    
    %-------------------------------------------------------------------------%
    % Note;
    % The exciting force acting uppon the system is different depending on the
    % case we are analysing. In the current version of the code we are
    % considering the following exciting forces:
    %
    %  (1) :  Straight centerline;
    %  (2) :  Parabolic centerline;
    %  (3) :  Cubic centerline;
    %-------------------------------------------------------------------------%

    switch caso
        
    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
        minHor     = 0;
        maxHor     = 1;
        minVer     = 0;
        maxVer     = 1;
    end

    %-------------------------------------------------------------------------%
    % NOTE:
    % If one wants to change the values of min_x and max_x, the procedure can
    % be performed inside the class 'AssemblerIGA'.
    %-------------------------------------------------------------------------%

    %% Discrtization parameters
    
    domainLimit_inX      = [minHor,maxHor];  
    domainLimit_inY      = [minVer,maxVer];
    
    numbModes       = 30;
    nd              = length(numbModes); 
    stepHorMesh     = (maxHor-minHor)*0.005*ones(size(numbModes));
    numbKnots    = round((maxHor-minHor)/stepHorMesh);
    
    %% Isogeometric basis properties

    % Polynomial Degree of the B-Spline Base
    
    degreeSplineBasis    = 2;
    
    % Continuity of the Base 'C^(p-k)'
    
    continuityParameter  = 1;
    
    % Number of control points
    
    numbControlPoints = numbKnots*(degreeSplineBasis - continuityParameter) + 1 + continuityParameter;
    
    if (degreeSplineBasis<continuityParameter)
    error('Wrong Choice for the Parameters!');
    end
    
    %% Time discretization properties
    
    Leff    = 3.7e-3;
    Tp      = (3.7e-3)^2 /(1e-5);
    
    period          = ww(jjj) * max(Tp);
    initialTime     = 0.0;
    finalTime       = 2.0 * Tp;
    numbSteps       = 400;
    timeStep        = (finalTime - initialTime)/numbSteps;
    timeDomain      = initialTime:timeStep:finalTime;
    initialState    = zeros(numbControlPoints * numbModes,1);

    %% Boundary conditions
    
    %-------------------------------------------------------------------------%
    % LATERAL BOUNDARY CONDITIONS
    %
    % * 'rob': mu du/dn u + chi u = dato
    % * 'dir': u = dato
    %
    % ATTENTION:
    % The algorithm works fine only if the lateral boundary conditions are the
    % same in the whole domain.
    %-------------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%
    % CHANGE BC HERE!
    %%%%%%%%%%%%%%%%%
    
    igaBoundCond.BC_UP_TAG     = 'neu';
    igaBoundCond.BC_DOWN_TAG   = 'neu';
    igaBoundCond.BC_INF_TAG   = 'dir';
    igaBoundCond.BC_OUT_TAG  = 'neu';
    igaBoundCond.BC_UP_DATA    = 0;
    igaBoundCond.BC_DOWN_DATA  = 0;
    igaBoundCond.BC_INF_DATA  = @(rho) 0;%0 + 1.75 * rho + 1.75 * -rho.^2;
    igaBoundCond.BC_OUT_DATA = @(rho) 0 + 0 * rho + 0 * rho.^2;
    
    BC_laterali_up   = igaBoundCond.BC_UP_TAG;
    BC_laterali_down = igaBoundCond.BC_DOWN_TAG;
    bc_up={ BC_laterali_up };
    dato_up={0};
    bc_down={ BC_laterali_down};
    dato_down={0};
    chi = 1; 

    %% Physical domain
    %---------------------------------------------------------------------%
    % Note: Complete domain is defined using the nurbs functions, not only
    % the centreline. This way, we can automatically compute the map, its
    % first and second derivatives and the jacobian of the transformation
    % from the pysical domain to the reference domain, where the reduction
    % procedure is defined.
    %---------------------------------------------------------------------%

    switch caso
        
    case {1}
        
        %%%%%%%%%%%%%%%%%
        % CIRCULAR RING %
        %%%%%%%%%%%%%%%%%
        
        R1 = 0.5;
        R2 = 1.5;
        c = [0.0 0.0 0.0];
        alpha = pi;
        
        % IMPORT THE GEOMETRY FILES FROM "Demos/ScatterGeometry" FOLDER
        
        geo_name = nrbrevolve(nrbline([R1 0 0], [R2 0 0]),c, [0 0 1], alpha);
        
        % CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE
        
        geometry = geo_load(geo_name);
        
        % EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
        % AND ITS JACOBIAN AND HESSIAN MATRICES
        
        map = @(x,y) geometry.map([x;y]);
        Jac = @(x,y) geometry.map_der([x;y]);
        Hes = @(x,y) geometry.map_der2([x;y]);
        
        geometricInfo.geometry = geometry;
        geometricInfo.map = @(x,y) map(x,y);
        geometricInfo.Jac = @(x,y) Jac(x,y);
        geometricInfo.Hes = @(x,y) Hes(x,y);
        geometricInfo.Type = 'Ring';
        
    case{2}
        
        %%%%%%%%%%%%%
        % RECTANGLE %
        %%%%%%%%%%%%%
        
        minX = +0.0;
        maxX = +3.0;
        minY = +0.0;
        maxY = +1.0;
                
        % IMPORT THE GEOMETRY FILES FROM "Demos/ScatterGeometry" FOLDER
        
        line1 = nrbline([minX minY],[maxX minY]);
        line2 = nrbline([minX maxY],[maxX maxY]);
        geo_name = nrbruled (line1, line2);
        
        % CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE
        
        geometry = geo_load(geo_name);
        
        % EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
        % AND ITS JACOBIAN AND HESSIAN MATRICES
        
        map = @(x,y) geometry.map([x;y]);
        Jac = @(x,y) geometry.map_der([x;y]);
        Hes = @(x,y) geometry.map_der2([x;y]);
        
        geometricInfo.geometry = geometry;
        geometricInfo.map = @(x,y) map(x,y);
        geometricInfo.Jac = @(x,y) Jac(x,y);
        geometricInfo.Hes = @(x,y) Hes(x,y);
        geometricInfo.Type = 'Rect';
        
    case {3}
        
        %%%%%%%%%%%%
        % STENOSIS %
        %%%%%%%%%%%%
        
        w  = 1;
        R  = 10;
        R2 = 20;
        r  = 2;
        n  = 100;
        p  = 1;
        range = (n/2) * (1/5);
        theta = linspace( 0 , pi/2 , n);

        xCoord = R * cos(theta);
        yCoord = R * sin(theta);
        zCoord = zeros(1,n);

        a = theta(ceil(n/2) - range);
        b = theta(ceil(n/2) + range);
        amp = 1 + range * 2;
        int = (b - a);
        func = @(x) exp(-( x - (a + b)/2 ).^2);

        xCoord(ceil(n/2) - range : ceil(n/2) + range) = (R + r * (sin(linspace(-0.5*pi,1.5*pi,amp)) + 1)) .* cos(theta(ceil(n/2) - range : ceil(n/2) + range));
        yCoord(ceil(n/2) - range : ceil(n/2) + range) = (R + r * (sin(linspace(-0.5*pi,1.5*pi,amp)) + 1)) .* sin(theta(ceil(n/2) - range : ceil(n/2) + range));

        %%%%%%%%%%
        % Weights 
        %%%%%%%%%%

        wCoord = ones(1,n);
        % wCoord(ceil(n/2) - range) = sqrt(2);
        % wCoord(ceil(n/2) + range) = sqrt(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create control points and knot vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        knot1 = linspace(0,1,n);
        knot1 = repmat(knot1,1,p)';
        knot1 = knot1(:)';
        knot1 = sort(knot1);

        strKnot = 0 * ones(1,p+1);
        endKnot = ones(1,p+1);
        knot1    = [ strKnot knot1(2:end-1) endKnot ];

        pnts = [xCoord ; yCoord ; zCoord; wCoord];
        crv1 = nrbmak(pnts,knot1);

        crv2 = nrbcirc(R2,[ 0 0 0 ],0,pi/2);

        srf = nrbruled(crv1,crv2);
        
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
        geometricInfo.Type = 'Stenosis';
        
    case {4}
        
        %%%%%%%%%%%%
        % ANEURYSM %
        %%%%%%%%%%%%
        
        w  = 1;
        R  = 10;
        R2 = 20;
        r  = 2;
        n  = 100;
        p  = 1;
        range = (n/2) * (1/5);
        theta = linspace( 0 , pi/2 , n);

        xCoord = R * cos(theta);
        yCoord = R * sin(theta);
        zCoord = zeros(1,n);

        a = theta(ceil(n/2) - range);
        b = theta(ceil(n/2) + range);
        amp = 1 + range * 2;
        int = (b - a);
        func = @(x) exp(-( x - (a + b)/2 ).^2);

        xCoord(ceil(n/2) - range : ceil(n/2) + range) = (R - r * (sin(linspace(-0.5*pi,1.5*pi,amp)) + 1)) .* cos(theta(ceil(n/2) - range : ceil(n/2) + range));
        yCoord(ceil(n/2) - range : ceil(n/2) + range) = (R - r * (sin(linspace(-0.5*pi,1.5*pi,amp)) + 1)) .* sin(theta(ceil(n/2) - range : ceil(n/2) + range));

        %%%%%%%%%%
        % Weights 
        %%%%%%%%%%

        wCoord = ones(1,n);
        % wCoord(ceil(n/2) - range) = sqrt(2);
        % wCoord(ceil(n/2) + range) = sqrt(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create control points and knot vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        knot1 = linspace(0,1,n);
        knot1 = repmat(knot1,1,p)';
        knot1 = knot1(:)';
        knot1 = sort(knot1);

        strKnot = 0 * ones(1,p+1);
        endKnot = ones(1,p+1);
        knot1    = [ strKnot knot1(2:end-1) endKnot ];

        pnts = [xCoord ; yCoord ; zCoord; wCoord];
        crv1 = nrbmak(pnts,knot1);

        crv2 = nrbcirc(R2,[ 0 0 0 ],0,pi/2);

        srf = nrbruled(crv1,crv2);
        
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
        geometricInfo.Type = 'Stenosis';
        
    case {5}
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Aneurism + Stenosis %
        %%%%%%%%%%%%%%%%%%%%%%%
        
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
        geometricInfo.Type = 'AneurismStenosis';
        
    case {6}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STENOSIS SMOOTH HEAVISIDE %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        n = 500;
        t = linspace(0,1,n);
        p = 1;
        
        th1 = 55 * pi/180;
        th2 = 35 * pi/180;
        ep = 5e-2;
        r1 = 0.3;
        r2 = 0.3;
        R1 = 1;
        R2 = 2;
        
        theta1 = pi/5;      % Position of the beggining stenosis
        theta2 = 3*pi/4;    % Position of the ending stenosis

        stenLesion = @(tt) 0.5 .* (1 + (2/pi).*atan((tt - th1)./ep)) - 0.5 .* (1 + (2/pi).*atan((tt - th2)./ep));

        lowerBorderX = @(tt) (R1 - r1 * stenLesion(tt * pi/2)) .* sin(tt * pi/2);
        lowerBorderY = @(tt) (R1 - r1 * stenLesion(tt * pi/2)) .* cos(tt * pi/2);
        upperBorderX = @(tt) (R2 + r2 * stenLesion(tt * pi/2)) .* sin(tt * pi/2);
        upperBorderY = @(tt) (R2 + r2 * stenLesion(tt * pi/2)) .* cos(tt * pi/2);
        
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
        
    case {7}
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % DRUG ELUTING PROBLEM %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        minX = +0.0;
        maxX = +42.0;
        minY = +0.0;
        maxY = +20.0;
                
        % IMPORT THE GEOMETRY FILES FROM "Demos/ScatterGeometry" FOLDER
        
        line1 = nrbline([minX minY],[maxX minY]);
        line2 = nrbline([minX maxY],[maxX maxY]);
        geo_name = nrbruled (line1, line2);
        
        % CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE
        
        geometry = geo_load(geo_name);
        
        % EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
        % AND ITS JACOBIAN AND HESSIAN MATRICES
        
        map = @(x,y) geometry.map([x;y]);
        Jac = @(x,y) geometry.map_der([x;y]);
        Hes = @(x,y) geometry.map_der2([x;y]);
        
        geometricInfo.geometry = geometry;
        geometricInfo.map = @(x,y) map(x,y);
        geometricInfo.Jac = @(x,y) Jac(x,y);
        geometricInfo.Hes = @(x,y) Hes(x,y);
        geometricInfo.Type = 'Rect';
        
    case {8,9,10,11}
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % DRUG ELUTING PROBLEM %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        minX = +0.0;
        maxX = 2 * 2 * 3.7e-3;
        minY = +0.0;
        maxY = 2 * 3.7e-3;
                
        % IMPORT THE GEOMETRY FILES FROM "Demos/ScatterGeometry" FOLDER
        
        line1 = nrbline([minX minY],[maxX minY]);
        line2 = nrbline([minX maxY],[maxX maxY]);
        geo_name = nrbruled (line1, line2);
        
        % CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE
        
        geometry = geo_load(geo_name);
        
        % EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
        % AND ITS JACOBIAN AND HESSIAN MATRICES
        
        map = @(x,y) geometry.map([x;y]);
        Jac = @(x,y) geometry.map_der([x;y]);
        Hes = @(x,y) geometry.map_der2([x;y]);
        
        geometricInfo.geometry = geometry;
        geometricInfo.map = @(x,y) map(x,y);
        geometricInfo.Jac = @(x,y) Jac(x,y);
        geometricInfo.Hes = @(x,y) Hes(x,y);
        geometricInfo.Type = 'Rect';
        
    case {12,13,14,15}
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % WOMERSLEY SIMULATION %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        minX = +0.0;
        maxX = +4.0;
        minY = +0.0;
        maxY = +1.0;
                
        % IMPORT THE GEOMETRY FILES FROM "Demos/ScatterGeometry" FOLDER
        
        line1 = nrbline([minX minY],[maxX minY]);
        line2 = nrbline([minX maxY],[maxX maxY]);
        geo_name = nrbruled (line1, line2);
        
        % CONSTRUCT THE GEOMETRY FROM THE IMPORT FILE
        
        geometry = geo_load(geo_name);
        
        % EXTRACT THE MAP FROM THE PHYSICAL DOMAIN TO THE REFERENCE DOMAIN
        % AND ITS JACOBIAN AND HESSIAN MATRICES
        
        map = @(x,y) geometry.map([x;y]);
        Jac = @(x,y) geometry.map_der([x;y]);
        Hes = @(x,y) geometry.map_der2([x;y]);
        
        geometricInfo.geometry = geometry;
        geometricInfo.map = @(x,y) map(x,y);
        geometricInfo.Jac = @(x,y) Jac(x,y);
        geometricInfo.Hes = @(x,y) Hes(x,y);
        geometricInfo.Type = 'Rect';

    end
    
    %% Quadrature propeties
    %---------------------------------------------------------------------%
    % Specifies the number of quadrature nodes to be used on the horizontal
    % and vertical direction to computes the integrals inside the Build
    % class.
    %---------------------------------------------------------------------%
    
    % Horizontal direction
    
    numbHorNodes = 32;
    
    % Vertical direction
    
    numbVerNodes = 10 * numbModes;

    %% Coefficients of the bilinear form
    %-------------------------------------------------------------------------%

    switch caso
    case {1,2,3,4,5}
        
        mu    = @(x,y,t) (  1.00 + 0*x + 0*y + 0*t ); % Difusion
        beta1 = @(x,y,t) (  0.00 + 0*x + 0*y + 0*t ); % Horizontal Advection
        beta2 = @(x,y,t) (  0.00 + 0*x + 0*y + 0*t ); % Vertical Advection
        sigma = @(x,y,t) (  0.00 + 0*x + 0*y + 0*t ); % Reaction
        
    case {6} 
        
        %%%%%%%%%%%%%%%%%%%
        % Poiseuille Flow %
        %%%%%%%%%%%%%%%%%%%
        
        % Component of the convective flow over the tangent and radial
        % directions with respect to the centerline of the domain
        
        bt = @(x,y) 1000 + 0*x + 0*y;
        br = @(x,y) 0 + 0*x + 0*y;
        b  = @(x,y) sqrt(bt(x,y).^2 + br(x,y).^2);
        
        % Modulating function to give the flow a parabolic profile with
        % respect to the direction of flow
        
        fm = @(x,y) (sqrt(x.^2 + y.^2) - R1) .* (sqrt(x.^2 + y.^2) - R2) * (-4/((R1 + R2)^2));
        
        % Equations for the horizontal and vertical components of the
        % convective field
        
        poiFlow1 = @(x,y) fm(x,y) .* b(x,y) .* cos( -pi/2 + atan2(y,x));
        poiFlow2 = @(x,y) fm(x,y) .* b(x,y) .* sin( -pi/2 + atan2(y,x));
        
        % Definition of the bilinear coefficients
        
        mu    = @(x,y,t) (  1 + 0*x + 0*y  + 0*t); % Difusion
        beta1 = @(x,y,t) poiFlow1(x,y) + 0*t; % Horizontal Advection
        beta2 = @(x,y,t) poiFlow2(x,y) + 0*t; % Vertical Advection
        sigma = @(x,y,t) (  0 + 0*x + 0*y + 0*t ); % Reaction
        
    case {7} 
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % DRUG ELUTING PROBLEM %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        % Definition of the bilinear coefficients
        
        mu    = @(x,y,t) (  1 + 0*x + 0*y + 0*t ); % Difusion
        beta1 = @(x,y,t) (  5 + 0*x + 0*y + 0*t ); % Horizontal Advection
        beta2 = @(x,y,t) (  0 + 0*x + 0*y + 0*t ); % Vertical Advection
        sigma = @(x,y,t) (  0 + 0*x + 0*y + 0*t ); % Reaction
        
    case {8} 
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % DRUG ELUTING PROBLEM %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        % Definition of the bilinear coefficients
        
        mu    = @(x,y,t) (  1 + 0*x + 0*y + 0*t ); % Difusion
        beta1 = @(x,y,t) (  5 + 0*x + 0*y + 0*t ); % Horizontal Advection
        beta2 = @(x,y,t) (  0 + 0*x + 0*y + 0*t ); % Vertical Advection
        sigma = @(x,y,t) (  0 + 0*x + 0*y + 0*t ); % Reaction
        
            
    case {9,10}
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIFFUSION COEFFICIENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Lower Stent Branches
        
        lx1 = 0.25 * maxX;
        dx1 = 0.07 * maxX;
        ly1 = 0.03 * maxY;
        dy1 = 0.10 * maxY;

        lx2 = 0.35 * maxX;
        dx2 = 0.07 * maxX;
        ly2 = 0.03 * maxY;
        dy2 = 0.10 * maxY;

        lx3 = 0.45 * maxX;
        dx3 = 0.07 * maxX;
        ly3 = 0.03 * maxY;
        dy3 = 0.10 * maxY;

        lx7 = 0.15 * maxX;
        dx7 = 0.07 * maxX;
        ly7 = 0.03 * maxY;
        dy7 = 0.10 * maxY;

        lx8 = 0.05 * maxX;
        dx8 = 0.07 * maxX;
        ly8 = 0.03 * maxY;
        dy8 = 0.10 * maxY;
        
        % Upper Stent Branches
        
        lx4 = 0.25 * maxX;
        dx4 = 0.07 * maxX;
        ly4 = 0.87 * maxY;
        dy4 = 0.10 * maxY;

        lx5 = 0.35 * maxX;
        dx5 = 0.07 * maxX;
        ly5 = 0.87 * maxY;
        dy5 = 0.10 * maxY;

        lx6 = 0.45 * maxX;
        dx6 = 0.07 * maxX;
        ly6 = 0.87 * maxY;
        dy6 = 0.10 * maxY;

        lx9 = 0.15 * maxX;
        dx9 = 0.07 * maxX;
        ly9 = 0.87 * maxY;
        dy9 = 0.10 * maxY;

        lx10 = 0.05 * maxX;
        dx10 = 0.07 * maxX;
        ly10 = 0.87 * maxY;
        dy10 = 0.10 * maxY;
        
        % Heaviside functions
        
        Square1D = @(xx,l,d) heaviside(xx - l) - heaviside(xx - l - d);
        Square2D = @(xxx,yyy,lx,ly,dx,dy) Square1D(xxx,lx,dx) .* Square1D(yyy,ly,dy);
        
        % Definition of the bilinear coefficients
        
        stent = @(x,y) (Square2D(x,y,lx1,ly1,dx1,dy1) + Square2D(x,y,lx2,ly2,dx2,dy2) + ...
                        Square2D(x,y,lx3,ly3,dx3,dy3) + Square2D(x,y,lx4,ly4,dx4,dy4) + ...
                        Square2D(x,y,lx5,ly5,dx5,dy5) + Square2D(x,y,lx6,ly6,dx6,dy6) + ...
                        Square2D(x,y,lx7,ly7,dx7,dy7) + Square2D(x,y,lx8,ly8,dx8,dy8) + ...
                        Square2D(x,y,lx9,ly9,dx9,dy9) + Square2D(x,y,lx10,ly10,dx10,dy10));
        
        % Definition of the bilinear coefficients
        
	Tbeta = Tp / ww(jjj);
	W     = 2 * pi / Tbeta;

        mu    = @(x,y,t) 1e-7 * stent(x,y) + 1e-6 * ~stent(x,y);
        beta1 = @(x,y,t) 0.5 * (1 + sin(W * t)) .* 5e-2 * (4/maxY^2) * (-y.^2 + y * maxY);
        beta2 = @(x,y,t) 0;
        sigma = @(x,y,t) 0;
        
    case {11,12,13,14,15} 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIFFUSION COEFFICIENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Lower Stent Branches
        
        lx1 = 0.5 * 0.25 * maxX;
        dx1 = 0.5 * 0.07 * maxX;
        ly1 = 0.01 * maxY;
        dy1 = 0.05 * maxY;
        
        lx2 = 0.5 * 0.35 * maxX;
        dx2 = 0.5 * 0.07 * maxX;
        ly2 = 0.01 * maxY;
        dy2 = 0.05 * maxY;
        
        lx3 = 0.5 * 0.45 * maxX;
        dx3 = 0.5 * 0.07 * maxX;
        ly3 = 0.01 * maxY;
        dy3 = 0.05 * maxY;
        
        lx7 = 0.5 * 0.15 * maxX;
        dx7 = 0.5 * 0.07 * maxX;
        ly7 = 0.01 * maxY;
        dy7 = 0.05 * maxY;
        
        lx8 = 0.5 * 0.05 * maxX;
        dx8 = 0.5 * 0.07 * maxX;
        ly8 = 0.01 * maxY;
        dy8 = 0.05 * maxY;
        
        % Upper Stent Branches
        
        lx4 = 0.5 * 0.25 * maxX;
        dx4 = 0.5 * 0.07 * maxX;
        ly4 = 0.94 * maxY;
        dy4 = 0.05 * maxY;

        lx5 = 0.5 * 0.35 * maxX;
        dx5 = 0.5 * 0.07 * maxX;
        ly5 = 0.94 * maxY;
        dy5 = 0.05 * maxY;

        lx6 = 0.5 * 0.45 * maxX;
        dx6 = 0.5 * 0.07 * maxX;
        ly6 = 0.94 * maxY;
        dy6 = 0.05 * maxY;

        lx9 = 0.5 * 0.15 * maxX;
        dx9 = 0.5 * 0.07 * maxX;
        ly9 = 0.94 * maxY;
        dy9 = 0.05 * maxY;

        lx10 = 0.5 * 0.05 * maxX;
        dx10 = 0.5 * 0.07 * maxX;
        ly10 = 0.94 * maxY;
        dy10 = 0.05 * maxY;
        
        % Heaviside functions
        
        Square1D = @(xx,l,d) heaviside(xx - l) - heaviside(xx - l - d);
        Square2D = @(xxx,yyy,lx,ly,dx,dy) Square1D(xxx,lx,dx) .* Square1D(yyy,ly,dy);
        
        % Definition of the bilinear coefficients
        
        stent = @(x,y,t) (Square2D(x,y,lx1,ly1,dx1,dy1) + Square2D(x,y,lx2,ly2,dx2,dy2) + ...
                          Square2D(x,y,lx3,ly3,dx3,dy3) + Square2D(x,y,lx4,ly4,dx4,dy4) + ...
                          Square2D(x,y,lx5,ly5,dx5,dy5) + Square2D(x,y,lx6,ly6,dx6,dy6) + ...
                          Square2D(x,y,lx7,ly7,dx7,dy7) + Square2D(x,y,lx8,ly8,dx8,dy8) + ...
                          Square2D(x,y,lx9,ly9,dx9,dy9) + Square2D(x,y,lx10,ly10,dx10,dy10)) + 0*t;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLOW PROPERTIES 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        R         = 0.5 * maxY;
        numbNodes = 50;
        Y         = linspace(0,R,numbNodes);
        Ns        = 100;
        T         = period;
        w         = 2 * pi ./ T;
        rho       = 1060; 
        mu        = 3.5e-3;
        omega     = w;

        alpha = R * sqrt(omega * rho / mu);
        
        disp(['Womersley Number   : ', num2str(alpha)])
        disp(['Fluid Density      : ', num2str(rho)])
        disp(['Dynamic Viscosity  : ', num2str(mu)])
        disp(['Omega (Ang. Freq.) : ', num2str(w)])
        disp(['Period Fluid       : ', num2str(period)])
        disp(['Charac. Time Diff. : ', num2str(Tp)])
        disp(['Rel. Period Fluid  : ', num2str(period./Tp)])
        disp(['Diffusion Coeff.   : ', num2str(MU(iii))])
        disp(['Initial Time       : ', num2str(initialTime)])
        disp(['Final Time         : ', num2str(finalTime)])
        disp(['Numb. of Steps     : ', num2str(numbSteps)])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EXTRACT THE BLOOD FLOW INFORMATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        BloodFlow = dlmread('FFTBloodFlow.txt');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EXTRACT FFT COEFFICIENTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        BloodFlowCmplx = complex(BloodFlow(:,2),BloodFlow(:,3));
        FlowProfile = ifft(BloodFlowCmplx);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOT THE BLOOD FLOW INFORMATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % figure('visible','off')
        % plot(FlowProfile(1:length(FlowProfile) * 0.5));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE WAVE FORM INFORMATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PPP = FlowProfile(1:length(FlowProfile) * 0.5);% - mean(FlowProfile(1:length(FlowProfile) * 0.5));
        fftPPP = fft(PPP);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE FOURIER COEFFICIENTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % NORMAL FOURIER

        N   = 15;
        dt  = T/length(PPP);
        p   = PPP;
        [p0,pn,phi] = FourierSeries(p,dt,N,T);

        tt = 0:dt:T;
        pt = p0;
        for n=1:N
          pt=pt+pn(n)*cos(2*pi*tt*n/T-phi(n));
        end

        % figure('visible','on')
        % plot(tt,pt,'-k',0:dt:T-dt,p,'x')
        % grid on

        % figure('visible','on')
        % plot(p0,'k-'); hold on;
        % plot(pn,'o'); hold on;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE FOURIER TRANSFORM OF THE BLOOD FLOW PROFILE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        WoProp.R       = R;
        WoProp.rho     = rho;
        WoProp.mu      = mu;
        WoProp.omega   = omega;
        WoProp.p0      = p0;
        WoProp.pn      = pn;
        WoProp.phi     = phi;
        WoProp.numbNodes = numbNodes;
        
        % Definition of the bilinear coefficients
        
        mu    = @(x,y,t) MU(iii) .* ((x - R).^2 + (y - R).^2 < (.25 * R)^2) + 1e-1 * MU(iii) .* ~((x - R).^2 + (y - R).^2 < (.25 * R)^2);%  .* ~stent(x,y,t) + 1e-10 .* stent(x,y,t);
        beta1 = @(x,y,t) y .* (y - maxY) .* (2 * 1e-2 / maxY); % 2e2 * Woo(WoProp,y,t);
        beta2 = @(x,y,t) 0;
        sigma = @(x,y,t) 0;
        
    end

    Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);

    %% Exact solution

    switch caso
    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
        true_sol  = @(x,y) (0.2*x.^5-0.5*(maxHor^3)*x.^2).*sin(2*pi*((y/maxVer)+0.5));
        true_solx = @(x,y) (x.^4-(maxHor^3)*x).*sin(2*pi*((y/maxVer)+0.5));
        true_soly = @(x,y) (0.2*x.^5-0.5*(maxHor^3)*x.^2)*2*pi.*cos(2*pi*((y/maxVer)+0.5));
    end

    %% Force and Dirichlet profile for the inflow data
    %-------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % NOTE:
    % This Dirichlet profile must be compatible witht the boundary conditions
    % of the problem.
    %-------------------------------------------------------------------------%

    switch caso
    case {1,2,4,5,6}
        
        dato_dir = @(y) 0;
        force = @(x,y,t) 1 + 0 * x + 0 * y + 0 * t;
    
    case {3}
        
        dato_dir = @(y) 0;
        force = @(x,y,t) 1 + 0 * x + 0 * y + 0 * t;
        
    case {7}
        
        % Circular Heaviside function
        
        C = 10;
        lx1 = 6;
        dx1 = 8;
        ly1 = 1;
        dy1 = 3;
        lx2 = 6;
        dx2 = 8;
        ly2 = 16;
        dy2 = 3;
        
        Square1D = @(xx,l,d) heaviside(xx - l) - heaviside(xx - l - d);
        Square2D = @(xxx,yyy,lx,ly,dx,dy) Square1D(xxx,lx,dx) .* Square1D(yyy,ly,dy);
        
        % Definition of the bilinear coefficients
        
        force    = @(x,y,t) C * (Square2D(x,y,lx1,ly1,dx1,dy1) + Square2D(x,y,lx2,ly2,dx2,dy2)) + 0*t; % Difusion
        
        % Plot source term
        
        xx = linspace(0,42,100);
        yy = linspace(0,20,100);
        [XX,YY] = meshgrid(xx,yy);
        ZZ = force(XX,YY);
        
        % figure('visible','off')
        % [C,h] = contourf(XX,YY,ZZ,10);
        % set(h,'LineColor','none')
        
    case {8,9,10,11}
        
        % Circular Heaviside function
        
        C = 10;
        
        % Definition of the bilinear coefficients
        
        force    = @(x,y,t) C * stent(x,y); % Difusion
        
        % Plot source term
        
        xx = linspace(0,maxX,100);
        yy = linspace(0,maxY,100);
        [XX,YY] = meshgrid(xx,yy);
        ZZ = force(XX,YY,1);
        
        % figure('visible','on')
        % [C,h] = contourf(XX,YY,ZZ,10);
        % set(h,'LineColor','none')
        
    case {12,13,14,15}
        
        % Definition of the bilinear coefficients
        
        force    = @(x,y,t) 10; % Difusion
        
    end

    Dati = struct('igaBoundCond',igaBoundCond,'force',force);

    %-------------------------------------------------------------------------%
    % Note;
    % The following loop varies in order to change the coefficients of the
    % interface and to try different configurations at the same time.
    %-------------------------------------------------------------------------%
    
    %% Solver
    % Definition of the Object of the EvaluationHandler Class

    import Core.SolverHandler

    obj_solverIGA = SolverHandler();

    % Properties Assignment

    obj_solverIGA.domainLimit_inX = domainLimit_inX;
    obj_solverIGA.domainLimit_inY = domainLimit_inY;
    obj_solverIGA.dimModalBasis = numbModes;
    obj_solverIGA.stepMeshX = stepHorMesh;
    obj_solverIGA.label_upBoundDomain = bc_up;
    obj_solverIGA.label_downBoundDomain = bc_down;
    obj_solverIGA.data_upBoundDomain = dato_up;
    obj_solverIGA.data_downBoundDomain = dato_down;
    obj_solverIGA.dirCondFuncStruct = Dati;
    obj_solverIGA.coefficientForm = Coeff_forma;
    obj_solverIGA.geometricInfo = geometricInfo;
    obj_solverIGA.dataExportOption = true;
    obj_solverIGA.simulationCase = caso;
    obj_solverIGA.exactSolution = true_sol;
    obj_solverIGA.exactSolution_dX = true_solx;
    obj_solverIGA.exactSolution_dY = true_soly;
    obj_solverIGA.degreePolySplineBasis = degreeSplineBasis;
    obj_solverIGA.continuityParameter = continuityParameter;
    obj_solverIGA.numbHorQuadNodes = numbHorNodes;
    obj_solverIGA.numbVerQuadNodes = numbVerNodes;
    obj_solverIGA.timeDomain = timeDomain;
    obj_solverIGA.timeStep = timeStep;
    obj_solverIGA.initialState = initialState;
    obj_solverIGA.numbControlPoints = numbControlPoints;
    
    warning('off','all')
    
    % COMPUTE HIGAMOD SOLUTION
    
    tic;
    [u2,plotStruct] = solverIGAScatterStentTransient(obj_solverIGA);
    toc;    
    
    % STORE SOLUTION STRUCTURE
    
    strPlotStruct{iii,jjj} = plotStruct;
    
    % SAVE PLOT STRUCTURE INTO .MAT FILE
    
    % save(sprintf(['Sim_',num2str(iii),'_',num2str(jjj),'.mat']),'-struct','plotStruct','-v7.3')
    
    end
end

% save StentTransient.mat -v7.3

%% Plot Results

plotSolutionTransientParallel(strPlotStruct);
