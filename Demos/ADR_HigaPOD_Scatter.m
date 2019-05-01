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
    disp('*           HIGAPOD Simulation           *');
    disp('******************************************');

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.PODSolverHandler

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
        
    case {2}
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
    
    numbModes       = 5;
    nd              = length(numbModes); 
    stepHorMesh     = (maxHor-minHor)*0.1*ones(size(numbModes));
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
    case {2}
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
    case {2}
        
        mu    = {{@(param) 0.5*param(1), @(x,y) (  0.24 + 0*x + 0*y )},
                 {@(param) 0.5*param(1), @(x,y) (  0.24 + 0*x + 0*y )}}; % Diffusion
        beta1 = {{@(param) 1, @(x,y) (  -5.0 + 0*x + 0*y )}}; % Horizontal Advection
        beta2 = {{@(param) 1, @(x,y) (  0.00 + 0*x + 0*y )}}; % Vertical Advection
        sigma = {{@(param) 1, @(x,y) (  0.00 + 0*x + 0*y )}}; % Reaction
        
    end
    
    % Mind the extra {...} around parametrized quantities, per matlab documentation on
    % struct containing cell arrays
    Coeff_forma = struct('mu', {mu}, 'beta1', {beta1}, 'beta2', {beta2}, 'sigma', {sigma}, 'coeffrobin',chi);
    %% Force and Dirichlet profile for the inflow data
    %-------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % NOTE:
    % This Dirichlet profile must be compatible witht the boundary conditions
    % of the problem.
    %-------------------------------------------------------------------------%

    switch caso
    case {2}
        
        dato_dir = @(y) 0;
        % force = @(x,y) 1 + 0 * x + 0 * y;
        
        force = {{@(param) 1, @(x,y) 50.*( (x>=2.7).*(x<=3.).*(y>=.35).*(y<=.65) + (x>=3.6).*(x<=3.9).*(y>=.35).*(y<=.65) )}};
    
    end

    Dati = struct('igaBoundCond', igaBoundCond, 'force', {force});

    %-------------------------------------------------------------------------%
    % Note;
    % The following loop varies in order to change the coefficients of the
    % interface and to try different configurations at the same time.
    %-------------------------------------------------------------------------%
    
    %% Solver
    % Definition of the Object of the EvaluationHandler Class

    import Core.PODSolverHandler

    obj_solverHigaPOD = PODSolverHandler();

    % Properties Assignment

    obj_solverHigaPOD.domainLimit_inX = domainLimit_inX;
    obj_solverHigaPOD.domainLimit_inY = domainLimit_inY;
    obj_solverHigaPOD.dimModalBasis = numbModes;
    obj_solverHigaPOD.stepMeshX = stepHorMesh;
    obj_solverHigaPOD.label_upBoundDomain = bc_up;
    obj_solverHigaPOD.label_downBoundDomain = bc_down;
    obj_solverHigaPOD.data_upBoundDomain = dato_up;
    obj_solverHigaPOD.data_downBoundDomain = dato_down;
    obj_solverHigaPOD.dirCondFuncStructExpansion = Dati;
    obj_solverHigaPOD.coefficientFormExpansion = Coeff_forma;
    obj_solverHigaPOD.geometricInfo = geometricInfo;
    obj_solverHigaPOD.dataExportOption = true;
    obj_solverHigaPOD.simulationCase = caso;
    obj_solverHigaPOD.degreePolySplineBasis = degreeSplineBasis;
    obj_solverHigaPOD.continuityParameter = continuityParameter;
    obj_solverHigaPOD.numbHorQuadNodes = numbHorNodes;
    obj_solverHigaPOD.numbVerQuadNodes = numbVerNodes;
    obj_solverHigaPOD.parameterRange = [[1,2]];
    
    tic;
    eigenvalues = offline(obj_solverHigaPOD, 100, 10);
    toc;
    
    [errors, speedups] = online(obj_solverHigaPOD, 30);
    
    eigenvalues
    errors
    speedups
