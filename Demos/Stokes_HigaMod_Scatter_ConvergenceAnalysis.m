%----------------------------------------------------------------------
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

    clear all
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

%% SIMULATION CASE

caso  = 2;      % Analysed Case

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

mVectP = [1 3 5 7 9 11 13];
hVectP = [0.2 0.1 0.05 0.025 0.0125];

structError = cell(length(mVectP),length(hVectP));

for ii = 1:length(mVectP)
    for jj = 1:length(hVectP)

    %% Discrtization parameters

    discStruct = [];

    modesVelc = mVectP(ii) + 2;
    modesPres = mVectP(ii);

    % Step used to generate the knot vector

    discStruct.stepHorMesh     = hVectP(jj);

    % Number of knots/elements in the isogeometric space

    discStruct.numbElements    = 1/discStruct.stepHorMesh;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISCRETIZATION PARAMETERS FOR THE PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    discStruct.numbModesP       = modesPres;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISCRETIZATION PARAMETERS FOR THE X COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Number of transverse modes

    discStruct.numbModesUx       = modesVelc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISCRETIZATION PARAMETERS FOR THE Y COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Number of transverse modes

    discStruct.numbModesUy       = modesVelc;

    %% Isogeometric basis properties

    igaBasisStruct = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IGA BASIS FOR THE PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Polynomial Degree of the B-Spline Base

    igaBasisStruct.degreeSplineBasisP    = 1;

    % Continuity of the Base 'C^(p-k)'

    igaBasisStruct.continuityParameterP  = 0;

    % Number of control points

    igaBasisStruct.numbControlPtsP = discStruct.numbElements * ...
                                     (igaBasisStruct.degreeSplineBasisP - igaBasisStruct.continuityParameterP) + ...
                                     1 + igaBasisStruct.continuityParameterP;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IGA BASIS FOR THE X COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Polynomial Degree of the B-Spline Base

    igaBasisStruct.degreeSplineBasisUx    = 2;

    % Continuity of the Base 'C^(p-k)'

    igaBasisStruct.continuityParameterUx  = 1;

    % Number of control points

    igaBasisStruct.numbControlPtsUx = discStruct.numbElements * ...
                                      (igaBasisStruct.degreeSplineBasisUx - igaBasisStruct.continuityParameterUx) + ...
                                      1 + igaBasisStruct.continuityParameterUx;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IGA BASIS FOR THE Y COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Polynomial Degree of the B-Spline Base

    igaBasisStruct.degreeSplineBasisUy    = 2;

    % Continuity of the Base 'C^(p-k)'

    igaBasisStruct.continuityParameterUy  = 1;

    % Number of control points

    igaBasisStruct.numbControlPtsUy = discStruct.numbElements * ...
                                      (igaBasisStruct.degreeSplineBasisUy - igaBasisStruct.continuityParameterUy) + ...
                                      1 + igaBasisStruct.continuityParameterUy;

    %% Boundary conditions

    boundCondStruct = [];

    %%%%%%%%%%%%%
    % PARAMTERS %
    %%%%%%%%%%%%%

    % INFLOW PRESSURE

    Pinf = 2;

    % OUTFLOW PRESSURE

    Pout = 0;

    % PRESSURE DROP

    dP = Pinf - Pout;

    % DYNAMIC VISCOSITY

    nu = 1;

    % LENGTH OF THE CYLINDER

    L = 5;

    % RADIUS OF THE CYLINDER

    R = 0.5;

    % PARABOLIC PROFILE PARAMETERS

    a0 = 0;
    a1 = (dP * R) / (2 * nu * L);
    a2 = (-1 * dP) / (4 * nu * L);

    Umax = a1^2 / 4 * a2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFLOW, OUTFLOW AND LATERAL BOUNDARY CONDITIONS FOR THE PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    boundCondStruct.bc_up_tag_P    = 'neu';
    boundCondStruct.bc_down_tag_P  = 'neu';
    boundCondStruct.bc_inf_tag_P   = 'neu';
    boundCondStruct.bc_out_tag_P   = 'neu';
    boundCondStruct.bc_up_data_P   = 0;
    boundCondStruct.bc_down_data_P = 0;
    boundCondStruct.bc_inf_data_P  = @(rho) 0 + 0 * rho + 0 * rho.^2;
    boundCondStruct.bc_out_data_P  = @(rho) 0 + 0 * rho + 0 * rho.^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFLOW, OUTFLOW AND LATERAL BOUNDARY CONDITIONS FOR Ux %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    boundCondStruct.bc_up_tag_Ux    = 'dir';
    boundCondStruct.bc_down_tag_Ux  = 'dir';
    boundCondStruct.bc_inf_tag_Ux   = 'dir';
    boundCondStruct.bc_out_tag_Ux   = 'neu';
    boundCondStruct.bc_up_data_Ux   = 0;
    boundCondStruct.bc_down_data_Ux = 0;
    boundCondStruct.bc_inf_data_Ux  = @(rho) a0 + a1 * rho + a2 * rho.^2;
    boundCondStruct.bc_out_data_Ux  = @(rho) 0 + 0 * rho + 0 * rho.^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFLOW, OUTFLOW AND LATERAL BOUNDARY CONDITIONS FOR Uy %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    boundCondStruct.bc_up_tag_Uy    = 'dir';
    boundCondStruct.bc_down_tag_Uy  = 'dir';
    boundCondStruct.bc_inf_tag_Uy   = 'dir';
    boundCondStruct.bc_out_tag_Uy   = 'dir';
    boundCondStruct.bc_up_data_Uy   = 0;
    boundCondStruct.bc_down_data_Uy = 0;
    boundCondStruct.bc_inf_data_Uy  = @(rho) 0 + 0 * rho + 0 * rho.^2;
    boundCondStruct.bc_out_data_Uy  = @(rho) 0 + 0 * rho + 0 * rho.^2;

    %% Physical domain
    %---------------------------------------------------------------------%
    % Note: Complete domain is defined using the nurbs functions, not only
    % the centreline. This way, we can automatically compute the map, its
    % first and second derivatives and the jacobian of the transformation
    % from the pysical domain to the reference domain, where the reduction
    % procedure is defined.
    %---------------------------------------------------------------------%

    switch caso

    case {1,6,7,8,9,10}

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

        minX = 0.0;
        maxX = L;
        minY = 0.0;
        maxY = 2 * R;

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

    end

    %% Quadrature propeties
    %---------------------------------------------------------------------%
    % Specifies the number of quadrature nodes to be used on the horizontal
    % and vertical direction to computes the integrals inside the Build
    % class.
    %---------------------------------------------------------------------%

    quadProperties = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pressure quadrature nodes along X direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    quadProperties.numbHorNodesP = 16;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pressure quadrature nodes along Y direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    quadProperties.numbVerNodesP = 64;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected X quadrature nodes along X direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    quadProperties.numbHorNodesUx = 16;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected X quadrature nodes along Y direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    quadProperties.numbVerNodesUx = 64;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected Y quadrature nodes along X direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    quadProperties.numbHorNodesUy = 16;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected Y quadrature nodes along Y direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    quadProperties.numbVerNodesUy = 64;

    %% Problem parameters

    probParameters = [];

    switch caso
    case {1,2,3,4,5,6,7,8,9,10}

        % Pressure compensation

        probParameters.delta    = @(x,y) ( 0.0 + 0*x + 0*y );

        % Artificial reaction

        probParameters.sigma    = @(x,y) ( 0.0 + 0*x + 0*y );

        % Fluid kinetic viscosity

        probParameters.nu    = @(x,y) (  nu + 0*x + 0*y );

        % Forcing term acting on the fluid

        probParameters.force_x = @(x,y) (  0.00 + 0*x + 0*y );
        probParameters.force_y = @(x,y) (  0.00 + 0*x + 0*y );

        % Robin coefficient for the modal basis

        probParameters.coeffrobin = 1;

    end

    %% Solver

    % Definition of the Object of the EvaluationHandler Class

    import Core.SolverHandler

    obj_solverIGA = SolverHandler();

    % Properties Assignment

    obj_solverIGA.discStruct = discStruct;
    obj_solverIGA.boundCondStruct = boundCondStruct;
    obj_solverIGA.igaBasisStruct = igaBasisStruct;
    obj_solverIGA.probParameters = probParameters;
    obj_solverIGA.geometricInfo = geometricInfo;
    obj_solverIGA.quadProperties = quadProperties;

    tic;
    [errStruct] = solverIGAScatterStokes(obj_solverIGA);
    toc;

    %% Store Error Information

    structError{ii,jj} = errStruct;

    end
end
    
%% EXPORT CONVERGENCE ANALYSIS RESULTS

exportConvAnalysisStokes(structError,mVectP,hVectP)