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

MM = [ 5 ];
HH = [ 0.1 ];
    
for ii = 1:length(MM)
    for jj = 1:length(HH)
    
    %% Simulation case   
    caso  = 3;      % Analysed Case
    
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

    %% Discrtization parameters
    
    domainLimit_inX      = [minHor,maxHor];  
    domainLimit_inY      = [minVer,maxVer];
    
    numbModes       = MM(ii);
    nd              = length(numbModes); 
    stepHorMesh     = (maxHor-minHor)*HH(jj)*ones(size(numbModes));
    numbElements    = round((maxHor-minHor)/stepHorMesh);
    
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

    BC_laterali_up   ='dir';
    BC_laterali_down ='dir';

    switch caso
    case {1,2,3,4,5,6,7,8}
        bc_up={ BC_laterali_up };
        dato_up={0};
        bc_down={ BC_laterali_down};
        dato_down={0};
    end
    
    % INFLOW AND OUTFLOW BOUNDARY CONDITIONS
    % Note: The current version of the code works only for constant
    % boundary conditions in inflow and outflow;
    
    % BOUNDARY CONDITION LABELS
    
    infLabelBC = 'neu';
    outLabelBC = 'neu';
    
    % DEFINE DIRICHLET AND NEUMANN BOUNDARY CONDITIONS
    
    infBCdir = @(x,y,z) (  0 + 0*x + 0*y + 0*z ); % Inflow
    outBCdir = @(x,y,z) (  0 + 0*x + 0*y + 0*z ); % Outflow
    infBCneu = @(x,y,z) (  0 + 0*x + 0*y + 0*z ); % Inflow
    outBCneu = @(x,y,z) (  0 + 0*x + 0*y + 0*z ); % Outflow
    
    % CREATE DATA STRUCTURE
    
    BoundCond.inflowDATAdir  = infBCdir;
    BoundCond.inflowDATAneu  = infBCneu;
    BoundCond.outflowDATAdir = outBCdir;
    BoundCond.outflowDATAneu = outBCneu;
    BoundCond.inflowINFO     = infLabelBC;
    BoundCond.outflowINFO    = outLabelBC;

    %% Physical domain
    %---------------------------------------------------------------------%
    % Note: Complete domain is defined using the nurbs functions, not only
    % the centreline. This way, we can automatically compute the map, its
    % first and second derivatives and the jacobian of the transformation
    % from the pysical domain to the reference domain, where the reduction
    % procedure is defined.
    %---------------------------------------------------------------------%

    switch caso
        
    case {1,4,5,6,7,8,9,10}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Artery from PATIENT 1 %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        currentFolder = pwd;
        
        folder = [pwd,'/Geometry/3D_Patient_1/'];
        cd(folder);
        
        filename1 = '1Patient3DGeo.mat';
        filename2 = '1Patient3DMap.mat';
        
        load(filename1);
        load(filename2);
        
        cd ..
        cd ..
        
        Vol = geo;
        
        % figure
        % nrbplot(Vol,[15 15 500]);
        
    case{2}
        
        %%%%%%%%%%%%
        % Cylinder %
        %%%%%%%%%%%%
        
        currentFolder = pwd;
        
        folder = [pwd,'/Geometry/3D_Cylinder/'];
        cd(folder);
        
        filename1 = 'Cylinder3DGeo.mat';
        filename2 = 'Cylinder3DMap.mat';
        
        load(filename1);
        load(filename2);
        
        cd ..
        cd ..
        
        Vol = geo;
        
        % figure
        % nrbplot(Vol,[10 10 10]);
        
    case {3}
        
        %%%%%%%%
        % Slab %
        %%%%%%%%
        
        currentFolder = pwd;
        
        folder = [pwd,'/Geometry/3D_Slab/'];
        cd(folder);
        
        filename1 = 'Slab3DGeo.mat';
        filename2 = 'Slab3DMap.mat';
        
        load(filename1);
        load(filename2);
        
        cd ..
        cd ..
        
        Vol = geo;

    end
    
    %% Quadrature propeties
    %---------------------------------------------------------------------%
    % Specifies the number of quadrature nodes to be used on the horizontal
    % and vertical direction to computes the integrals inside the Build
    % class.
    %---------------------------------------------------------------------%
    
    % Horizontal direction
    
    numbHorNodes = 8;
    
    % Vertical direction
    
    numbVerNodes = numbModes * 3;

    %% Coefficients of the bilinear form
    %-------------------------------------------------------------------------%

    switch caso
    case {1,2,3,4,5,6,7,8,9,10}
        mu    = @(x,y,z) (  1 + 0*x + 0*y + 0*z ); % Difusion
        beta1 = @(x,y,z) (  1 + 0*x + 0*y + 0*z ); % 1st Component Convective Field
        beta2 = @(x,y,z) (  1 + 0*x + 0*y + 0*z ); % 2nd Component Convective Field
        beta3 = @(x,y,z) (  1 + 0*x + 0*y + 0*z ); % 3rd Component Convective Field
        sigma = @(x,y,z) (  0 + 0*x + 0*y + 0*z ); % Reaction
    end
    
    Coeff_forma.mu = mu;
    Coeff_forma.beta1 = beta1;
    Coeff_forma.beta2 = beta2;
    Coeff_forma.beta3 = beta3;
    Coeff_forma.sigma = sigma;
    Coeff_forma.coeffrobin = chi;

    %% Exact solution

    switch caso
    case {1,2,3,4,5,6,7,8,9,10}
        true_sol  = @(x,y,z) 0 + 0*x + 0*y + 0*z;
        true_solx = @(x,y,z) 0 + 0*x + 0*y + 0*z;
        true_soly = @(x,y,z) 0 + 0*x + 0*y + 0*z;
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
        force = @(x,y,z) 1 + 0 * x + 0 * y + 0*z;    
    case {3}
        dato_dir = @(y) 0;
        force = @(x,y,z) 1 + 0 * x + 0 * y + 0*z;
    end

    Dati.BoundCond = BoundCond;
    Dati.force = force;

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
    [u2,a,b,L2_2,H1_2] = solverIGAScatter3D(obj_solverIGA);
    toc;
    
    disp('Maximum L2 Norm Error with Matlab');
    disp(max(max(L2_2)));
    
    disp('Maximum H1 Norm Error with Matlab');
    disp(max(max(H1_2)));
    end
end