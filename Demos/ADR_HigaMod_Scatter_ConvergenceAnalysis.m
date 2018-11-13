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

%     clear all
%     close all

    disp('******************************************')
    disp('*           HIGAMod Simulation           *');
    disp('******************************************');

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler

    %% Simulation case   
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

    switch caso
        
    case {1,2,3,4,5,6,7,8,9,10}
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
    
    mVect = [ 1 3 5 7 9 11 13 15];
    hVect = [ 0.05 0.025 0.0125 0.00625 0.003125];
    
    matErrL2 = zeros(length(mVect),length(hVect));
    matErrH1 = zeros(length(mVect),length(hVect));

    for ii = 1:length(mVect)
        for jj = 1:length(hVect)

            %% Discrtization parameters

            domainLimit_inX      = [minHor,maxHor];  
            domainLimit_inY      = [minVer,maxVer];

            numbModes       = mVect(ii);
            nd              = length(numbModes); 
            stepHorMesh     = (maxHor-minHor) * hVect(jj) * ones(size(numbModes));
            numbElements    = round((maxHor-minHor)/stepHorMesh);

            %% Isogeometric basis properties

            % Polynomial Degree of the B-Spline Base

            degreeSplineBasis    = 2;

            % Continuity of the Base 'C^(p-k)'

            continuityParameter  = 1;

            % Number of control points

            numbControlPts = numbElements * continuityParameter + degreeSplineBasis + 1 -continuityParameter;

            if (degreeSplineBasis<continuityParameter)
                error('Wrong Choice for the Parameters!');
            end

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

            dato_dir_up   = 0.0;
            dato_dir_down = 0.0;
            chi           = 1;
            mu            = 1;
            cest          = 1;

            BC_laterali_up   ='dir';
            BC_laterali_down ='dir';

            switch caso
            case {1,2,3,4,5,6,7,8,9,10}
                bc_up={ BC_laterali_up };
                dato_up={0};
                bc_down={ BC_laterali_down};
                dato_down={0};
            end

            % INFLOW AND OUTFLOW BOUNDARY CONDITIONS
            % Note: The current version of the code works only for constant
            % boundary conditions in inflow and outflow;

            % SELECT THE SIDES

            dirSides = [1 2];
            neuSides = [];
            robSides = [];

            % DEFINE DIRICHLET AND NEUMANN BOUNDARY CONDITIONS

            dir  = @(x,side) (side == 1) * 0 + (side == 2) * 0;
            neu  = @(x,side) (side == 1) * 0 + (side == 2) * 1;
            rob.value  = @(x,side) (side == 1) * 0 + (side == 2) * 1;
            rob.mu     = 1.0;
            rob.chi    = 1.0;

            % CREATE DATA STRUCTURE

            igaBoundCond.dirSides   = dirSides;
            igaBoundCond.neuSides   = neuSides;
            igaBoundCond.robSides   = robSides;
            igaBoundCond.dir  = @(x,side) dir(x,side);
            igaBoundCond.neu  = @(x,side) neu(x,side);
            igaBoundCond.rob.value  = @(x,side) rob.value;
            igaBoundCond.rob.mu     = rob.mu;
            igaBoundCond.rob.chi    = rob.chi;

            %%%%%%%%%%%%%%%%%
            % CHANGE BC HERE!
            %%%%%%%%%%%%%%%%%

            igaBoundCond.BC_UP_TAG     = 'dir';
            igaBoundCond.BC_DOWN_TAG   = 'dir';
            igaBoundCond.BC_INF_TAG   = 'dir';
            igaBoundCond.BC_OUT_TAG  = 'dir';
            igaBoundCond.BC_UP_DATA    = 0;
            igaBoundCond.BC_DOWN_DATA  = 0;
            igaBoundCond.BC_INF_DATA  = @(rho) 0 + 0 * rho + 0 * rho.^2;
            igaBoundCond.BC_OUT_DATA = @(rho) 0 + 0 * rho + 0 * rho.^2;

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

                minX = -1.0;
                maxX = +5.0;
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

            end

            %% Quadrature propeties
            %---------------------------------------------------------------------%
            % Specifies the number of quadrature nodes to be used on the horizontal
            % and vertical direction to computes the integrals inside the Build
            % class.
            %---------------------------------------------------------------------%

            % Horizontal direction

            numbHorNodes = 16;

            % Vertical direction

            numbVerNodes = 32;

            %% Coefficients of the bilinear form
            %-------------------------------------------------------------------------%

            switch caso
            case {1,2,3,4,5,6}

                mu    = @(x,y) (  0.24 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  -5.0 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0.00 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0.00 + 0*x + 0*y ); % Reaction

            case {7,8} 

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

                mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) poiFlow1(x,y); % Horizontal Advection
                beta2 = @(x,y) poiFlow2(x,y); % Vertical Advection
                sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction

            case {9,10}

                %%%%%%%%%%%%%%%%%%%
                % Heaviside Coeff %
                %%%%%%%%%%%%%%%%%%%

                % Circular Heaviside function

                c = [0 1];
                lag = .25;
                delta = 100;

                Hc = @(x,y,a,d) 1./(1 + exp(-2 .* d .* (sqrt((x - c(1)).^2 + (y - c(2)).^2) - a)));
                func = @(x,y,a,d) Hc(x,y,-a,d) - Hc(x,y,a,d);

                % Definition of the bilinear coefficients

                mu    = @(x,y) 1 + 100 * func(x,y,lag,delta); % Difusion
                beta1 = @(x,y) (  0 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0 + 0*x + 0*y ); %func(x,y,lag,delta); % Reaction

            end

            Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
            'sigma',sigma,'coeffrobin',chi);

            %% Exact solution

            switch caso
            case {1,2,3,4,5,6,7,8,9,10}
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
            case {1,2,4,5,6,7,8,9,10}

                dato_dir = @(y) 0;
                % force = @(x,y) 1 + 0 * x + 0 * y;

                force = @(x,y) 50.*( (x>=2.7).*(x<=3.).*(y>=.35).*(y<=.65) + (x>=3.6).*(x<=3.9).*(y>=.35).*(y<=.65) );

            case {3}

                dato_dir = @(y) 0;
                force = @(x,y) 1 + 0 * x + 0 * y;

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

            tic;
            [u2,a,b,L2_2,H1_2] = solverIGAScatter(obj_solverIGA);
            toc;

            disp('Maximum L2 Norm Error with Matlab');
            disp(max(max(L2_2)));

            disp('Maximum H1 Norm Error with Matlab');
            disp(max(max(H1_2)));
            
            %% Store Error Information
            
            matErrL2(ii,jj) = L2_2;
            matErrH1(ii,jj) = H1_2;
            
        end
    end
    
%% Convergence Analysis Plot

figL2 = figure;
for kk = 1:length(hVect)
    loglog(mVect,matErrL2(:,kk),'-o','LineWidth',3,'MarkerSize',3);
    set(gcf, 'Color', 'w');
    set(gca, 'FontSize', 14);
    xlim auto
    ylim auto
    grid on
    hold on
end

loglog(mVect,mVect.^-1,'--','LineWidth',2);
loglog(mVect,mVect.^-2,'--','LineWidth',2);

LegendTitles = cell(1,length(hVect));
for kk = 1:length(hVect)
    LegendTitles{kk} = ['h = ' num2str(hVect(kk))];
end
LegendTitles{kk + 1} = ['Order 1'];
LegendTitles{kk + 2} = ['Order 2'];
legend(LegendTitles,'Location','northeast')
set(figL2, 'Visible', 'on')
export_fig(sprintf(['ConvAnalysisL2']),'-pdf');

figH1 = figure;
for kk = 1:length(hVect)
    loglog(mVect,matErrH1(:,kk),'-o','LineWidth',3,'MarkerSize',3);
    set(gcf, 'Color', 'w');
    set(gca, 'FontSize', 14);
    xlim auto
    ylim auto
    grid on
    hold on
end

loglog(mVect,mVect.^-1,'--','LineWidth',2);
loglog(mVect,mVect.^-2,'--','LineWidth',2);

LegendTitles = cell(1,length(hVect));
for kk = 1:length(hVect)
    LegendTitles{kk} = ['h = ' num2str(hVect(kk))];
end
LegendTitles{kk + 1} = ['Order 1'];
LegendTitles{kk + 2} = ['Order 2'];
legend(LegendTitles,'Location','northeast')
set(figL2, 'Visible', 'on')
export_fig(sprintf(['ConvAnalysisH1']),'-pdf');