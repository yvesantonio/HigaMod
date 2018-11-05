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

    clc
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

    %% Simulation case   
    simulationCase  = 4;      % Analysed Case
    
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

    switch simulationCase
        case {1,3,4,5,6,7,10}
            minHor     = 0;
            maxHor     = 1;
            minVer     = 0;
            maxVer     = 1;
        
        case {2,8,9}
            
            minHor     = 0;
            maxHor     = 2;
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
    
    numbModes       = 5;
    nd              = length(numbModes); 
    stepHorMesh     = (maxHor-minHor)*0.1*ones(size(numbModes));
    numbElements    = round((maxHor-minHor)/stepHorMesh);

    timeStep   = .02;
    endTime    = .1 ;
    timeDomain = 0:timeStep:endTime;
   
    %% Isogeometric basis properties

    % Polynomial Degree of the B-Spline Base
    
    degreeSplineBasis    = 1;
    
    % Continuity of the Base 'C^(p-k)'
    
    continuityParameter  = 0;
    
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

    switch simulationCase
        
        case 1
            
            BC_laterali_up   ='dir';
            BC_laterali_down ='dir';

            bc_up={ BC_laterali_up };
            dato_up = 0;
            bc_down={ BC_laterali_down};
            dato_down = 0;
        
        case {2,3,4,5,6,7,8}
        
            BC_laterali_up   ='dir';
            BC_laterali_down ='dir';

            bc_up={ BC_laterali_up };
            dato_up = 0;
            bc_down={ BC_laterali_down};
            dato_down = 0;
            
        case 9
            
            BC_laterali_up   ='dir';
            BC_laterali_down ='dir';

            bc_up={ BC_laterali_up };
            dato_up= 0;
            bc_down={ BC_laterali_down};
            dato_down= 0;
            
    end
    
    %% Initial State
    
    % initialState = [0;0.350390794315583;0.424581602787768;0.445184966653434;0.449680813103781;0.450143306783941;0.449680813103781;0.445184966653434;0.424581602787768;0.350390794315582;0.108878893429393;0;-2.35727052384867e-16;-2.61059026450664e-16;-1.41974433446596e-16;-1.36289364330610e-16;-1.33828951861278e-16;-1.36289364330610e-16;-1.41974433446596e-16;-2.61059026450664e-16;-2.35727052384867e-16;-8.14996194947619e-17;0;-0.0467371606955923;0.0809713333396734;0.135491573848283;0.148626166437808;0.150008175278703;0.148626166437808;0.135491573848283;0.0809713333396733;-0.0467371606955923;-0.0610990572829693;0;2.23836998553038e-16;3.78902785682020e-16;5.32238696762621e-16;5.09924044465286e-16;5.03797018839029e-16;5.09924044465286e-16;5.32238696762620e-16;3.78902785682020e-16;2.23836998553038e-16;1.74893747746089e-16;0;-0.0585270401585809;-0.00333999882720612;0.0669121291438901;0.0876721943794969;0.0899574241439401;0.0876721943794969;0.0669121291438901;-0.00333999882720613;-0.0585270401585809;0.0119533599405041];
    initialState = zeros(numbControlPts * numbModes,1);
%     for ii = numbModes
%         inflowDir  = 10;
%         outflowDir = 10;
%         
%         initialState((ii - 1) * numbControlPts + 1) = inflowDir;
%         initialState( ii * numbControlPts - 1) = inflowDir;
%     end

    %% Physical domain

    %-------------------------------------------------------------------------%
    % DATA RELATIVE TO THE SHAPE OF THE DOMAIN
    %-------------------------------------------------------------------------%
    % NOTE:
    % The current algorithm works only for the following set of parameters:
    %
    % * L = 1
    % * psi_x = 0
    %
    % Actually, if we have different values for those parameters, we must
    % verify again the coefficients, especially in the assembling methods
    % inside the AssemblerADRHandler class.
    %-------------------------------------------------------------------------%

    switch simulationCase
        
    case 6
        a     = @(x) 10 * sin(pi*x/5);                   % Equation of the Inferior Bound of the Channel (Considered As Central)
        L     = @(x) (maxVer - minVer) + 0*x; % Thickness of the Channel
        D1    = 0;
        D2    = 1;                          %/(max_y - min_y);
        psi_x = @(x) 0*x;                   % Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) 10 * sin(pi * x/5);
        domainProfileDer = @(x) 2 * pi * cos(pi * x/5);

    case {5}
        
        a     = @(x) x.^2;                   % Equation of the Inferior Bound of the Channel (Considered As Central)
        L     = @(x) (maxVer - minVer) + 0*x; % Thickness of the Channel
        D1    = 0;
        D2    = 1;                          %/(max_y - min_y);
        psi_x = @(x) 0*x;                   % Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) x.^2;
        domainProfileDer = @(x) 2*x;

    case {1,2,3}
        a     = @(x) 0*x;                   % Equation of the Inferior Bound of the Channel (Considered As Central)
        L     = @(x) (maxVer - minVer) + 0*x; % Thickness of the Channel
        D1    = 0;
        D2    = 1;                          %/(max_y - min_y);
        psi_x = @(x) 0*x;                   % Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) 0*x;
        domainProfileDer = @(x) 0*x;
        
     case 4
        a     = @(x) 0*x;                   % Equation of the Inferior Bound of the Channel (Considered As Central)
        L     = @(x) (maxVer - minVer) + 0*x; % Thickness of the Channel
        D1    = 0;
        D2    = 1/(maxVer - minVer);
        psi_x = @(x) 0*x;                   % Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) 0*x;
        domainProfileDer = @(x) 0*x;
        
     case 7
         
        a     = @(x) (1 + sqrt(1 - x.^2)).*(x<1) + (1 + sqrt(1 - (x - 2).^2)).*(x>=1);
        L     = @(x) (maxVer - minVer) + 0*x; % Thickness of the Channel
        D1    = 0;
        D2    = 1/(maxVer - minVer);
        psi_x = @(x) 0*x;                   % Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) (1 - sqrt(1 - x.^2)).*(x<1) + (1 + sqrt(1 - (x - 2).^2)).*(x>=1);
        domainProfileDer = @(x) (0.5 * (2 * x)./sqrt(1 - x.^2)).*(x<1) + (-(x-2)./sqrt(1 - (x - 2).^2)).*(x>=1);
        
     case {8,9}
         
        a     = @(x) sqrt(1 - (x - 1).^2);
        L     = @(x) (maxVer - minVer) + 0*x; % Thickness of the Channel
        D1    = 0;
        D2    = 1/(maxVer - minVer);
        psi_x = @(x) 0*x;                   % Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) sqrt(1 - (x - 1).^2);
        domainProfileDer = @(x) -(x-1)./sqrt(1 - (x - 1).^2);


    end
    
    %% Quadrature propeties
    %---------------------------------------------------------------------%
    % Specifies the number of quadrature nodes to be used on the horizontal
    % and vertical direction to computes the integrals inside the Build
    % class.
    %---------------------------------------------------------------------%
    
    % Horizontal direction
    
    numbHorNodes = 64;
    
    % Vertical direction
    
    numbVerNodes = 64;
    
    %% Vertical Mesh

    mesh_fisy = [];

    %% Horizontal Mesh
    % For the X Axis

    knot      = augknt([minHor maxHor],degreeSplineBasis+1);
    h         = stepHorMesh;
    nel       = numbElements;
    
    internKnot = linspace(0,1,numbElements+1)';
    internKnot = internKnot(2:end-1);
    
    ins       = sort(reshape(internKnot*ones(1,continuityParameter),1,[]));      %continuityParameter*(nel-1)));    
    ins       = (maxHor-minHor)*ins+minHor;
    
    % Creates the control points and knots points of the desired
    % B-Spline curve

    [cp,knot] = bspkntins(degreeSplineBasis,minHor:(maxHor-minHor)*1/degreeSplineBasis:maxHor,knot,ins);

    mesh_fisx  = linspace(minHor,maxHor,numbElements+1);
    mesh_fisx2 = linspace(minHor,maxHor,2*numbControlPts);

    psJac  = [];
    psJac2 = [];
    num    = 8;

    for ii = 1:numbElements

        pos1 = ii;
        pos2 = ii+1;

        % 1. Take 8 Quadrature Nodes and Respective Weights

        import Core.IntegrateHandler

        obj_gaussLegendre_1 = IntegrateHandler();

        obj_gaussLegendre_1.numbQuadNodes = num;

        [numnodes,x_nodes,wei] = gaussLegendre(obj_gaussLegendre_1); 

        % 2. Rescale Them Depending on the Iteration

        obj_quadratureRule_1 = IntegrateHandler();

        obj_quadratureRule_1.leftBoundInterval = mesh_fisx(ii);
        obj_quadratureRule_1.rightBoundInterval = mesh_fisx(ii+1);
        obj_quadratureRule_1.inputNodes = x_nodes;
        obj_quadratureRule_1.inputWeights = wei;

        [nod,wei] = quadratureRule(obj_quadratureRule_1);

        % 3. Definition of the Jacobian Computed for Every Node

        jaco = [];

        for it = 1:num

            jaco(it) = sqrt(1+(domainProfileDer(nod(it))).^2); 

        end

        % 4. Computation of the Length of the Quadrature Rule

        % 5. Save the Information in psJac

        psJac = [psJac sum(jaco.*wei')];

    end

    %% Coefficients of the bilinear form
    %-------------------------------------------------------------------------%

    switch simulationCase
    case {1,2,3,4,5,6,7,8,9}

        mu    = @(x,y,t) (   + 0*x + 0*y + 0*t); % Difusion
        beta1 = @(x,y,t) (  0 + 0*x + 0*x + 0*t ); % Horizontal Advection
        beta2 = @(x,y,t) (  0 + 0*x + 0*y + 0*t ); % Vertical Advection
        sigma = @(x,y,t) (  0 + 0*x + 0*y + 0*t ); % Reaction

    end

    Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);

    %% Exact solution
    %-------------------------------------------------------------------------%
    
    switch simulationCase
        
    case {1} % Homogeneous Neumann B.C. in the Vertical direction + HDHN

        true_sol  = @(x,y,t) cos(pi*y/maxVer) .* sin(pi*x/(2*maxHor)) .* 10*sin(t);
        true_solx = @(x,y,t) cos(pi*y/maxVer) .* cos(pi*x/(2*maxHor)) .* (pi/(2*maxHor)) .* 10*sin(t);
        true_soly = @(x,y,t) sin(pi*y/maxVer) .* sin(pi*x/(2*maxHor)) .* (-pi/(maxVer)) .* 10*sin(t);
        
    case {2} % Homogeneous Dirichlet B.C. in the Vertical direction + HDHN

        true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*(exp(0.1*t) - 1);
        true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*(exp(0.1*t) - 1);
        true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*(exp(0.1*t) - 1);


    case {3,4,5,6,7,8,9}

        true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*(exp(0.1*t) - 1);
        true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*(exp(0.1*t) - 1);
        true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*(exp(0.1*t) - 1);

    end

    %% Force and Dirichlet profile for the inflow data
    %-------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % NOTE:
    % This Dirichlet profile must be compatible witht the boundary conditions
    % of the problem.
    %-------------------------------------------------------------------------%

    switch simulationCase
        
    case 1 % Homogeneous Neumann B.C. in the Vertical direction + HDHN

        dato_dir = @(y) 0;
        force = @(x,y,t) mu(x,y,t) .* (pi^2/(4*maxHor^2)) .* cos(pi*y/maxVer) .* sin(pi*x/(2*maxHor)) .* 10*sin(t) + ...
                       mu(x,y,t) .* (pi^2/(maxVer^2)) .* cos(pi*y/maxVer) .* sin(pi*x/(2*maxHor)) .* 10*sin(t) + ...
                       beta1(x,y,t) .* cos(pi*y/maxVer) .* cos(pi*x/(2*maxHor)) .* (pi/(2*maxHor)) .* 10*sin(t) + ...
                       beta2(x,y,t) .* sin(pi*y/maxVer) .* sin(pi*x/(2*maxHor)) .* (-pi/(maxVer)) .* 10*sin(t) + ...
                       sigma(x,y,t) .* cos(pi*y/maxVer) .* sin(pi*x/(2*maxHor)) .* 10*sin(t) + ...
                       10 * cos(t) .* cos(pi*y/maxVer) .* sin(pi*x/(2*maxHor));
                   
    case 2 % Homogeneous Dirichlet B.C. in the Vertical direction + HDHN

        dato_dir = @(y) 0;
        force = @(x,y,t) -(4*x.^3-1).*sin(2*pi*(y+0.5)).*(exp(0.1*t) - 1)...
                         +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5)).*(exp(0.1*t) - 1)...
                         +(0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t)*0.1;
    
    case {3,4,5,6,7,8}

        dato_dir = @(y) 0;
        force = @(x,y,t) 10 * ((x-.5).^2+(y-.5).^2<0.25); 
        
    case {9}
        
        dato_dir = @(y) 0;
        force = @(x,y,t) 0; 

    end

    Dati = struct('dato_dir',dato_dir,'force',force);

    %-------------------------------------------------------------------------%
    % Note;
    % The following loop varies in order to change the coefficients of the
    % interface and to try different configurations at the same time.
    %-------------------------------------------------------------------------%
    
    %% Solver
    % Definition of the Object of the EvaluationHandler Class

    import Core.SolverHandler

    obj_solverIGATransientFull = SolverHandler();

    % Properties Assignment

    obj_solverIGATransientFull.domainLimit_inX = domainLimit_inX;
    obj_solverIGATransientFull.domainLimit_inY = domainLimit_inY;
    obj_solverIGATransientFull.numbControlPoints = numbControlPts;
    obj_solverIGATransientFull.dimModalBasis = numbModes;
    obj_solverIGATransientFull.stepMeshX = stepHorMesh;
    obj_solverIGATransientFull.label_upBoundDomain = bc_up;
    obj_solverIGATransientFull.label_downBoundDomain = bc_down;
    obj_solverIGATransientFull.data_upBoundDomain = dato_up;
    obj_solverIGATransientFull.data_downBoundDomain = dato_down;
    obj_solverIGATransientFull.dirCondFuncStruct = Dati;
    obj_solverIGATransientFull.coefficientForm = Coeff_forma;
    obj_solverIGATransientFull.geometricInfo = geometricInfo;
    obj_solverIGATransientFull.dataExportOption = true;
    obj_solverIGATransientFull.simulationCase = simulationCase;
    obj_solverIGATransientFull.exactSolution = true_sol;
    obj_solverIGATransientFull.exactSolution_dX = true_solx;
    obj_solverIGATransientFull.exactSolution_dY = true_soly;
    obj_solverIGATransientFull.physicMesh_inX = mesh_fisx;
    obj_solverIGATransientFull.physicMesh_inY = mesh_fisy;
    obj_solverIGATransientFull.domainProfile = domainProfile;
    obj_solverIGATransientFull.domainProfileDer = domainProfileDer;
    obj_solverIGATransientFull.jacAtQuadNodes = psJac;
    obj_solverIGATransientFull.jacAtQuadNodes2 = psJac2; 
    obj_solverIGATransientFull.degreePolySplineBasis = degreeSplineBasis;
    obj_solverIGATransientFull.continuityParameter = continuityParameter;
    obj_solverIGATransientFull.numbHorQuadNodes = numbHorNodes;
    obj_solverIGATransientFull.numbVerQuadNodes = numbVerNodes;
    obj_solverIGATransientFull.D1 = D1;
    obj_solverIGATransientFull.D2 = D2;
    obj_solverIGATransientFull.timeStep = timeStep;
    obj_solverIGATransientFull.timeDomain = timeDomain;
    obj_solverIGATransientFull.initialState = initialState;

    tic;
    [solutionMatrix,liftCoeffA,liftCoeffB,errorNormL2,errorNormH1,stateMatrix, forceMatrix, forceHistory,force] = solverIGATransientFull(obj_solverIGATransientFull);
    toc;
    
    % Plot the error evolution
    
%     figure
%     plot(timeDomain(2:end),errorNormL2,'-o');
%     
%     figure
%     plot(timeDomain(2:end),errorNormH1,'-*');