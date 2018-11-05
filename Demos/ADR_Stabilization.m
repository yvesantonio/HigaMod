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

    %% Simulation case   
    caso  = 1;      % Analysed Case
    
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
    case {1,2,3}
        min_x     = 0;
        max_x     = 1;
        min_y     = 0;
        max_y     = 1;
    end

    %-------------------------------------------------------------------------%
    % NOTE:
    % If one wants to change the values of min_x and max_x, the procedure can
    % be performed inside the class 'AssemblerIGA'.
    %-------------------------------------------------------------------------%

    %% Discrtization parameters
    
    domainLimit_inX      = [min_x,max_x];  
    domainLimit_inY      = [min_y,max_y];
    
    numbModes       = 6;
    nd              = length(numbModes); 
    stepHorMesh     = (max_x-min_x)*0.1*ones(size(numbModes));
    numbElements    = round((max_x-min_x)/stepHorMesh);
    
    %% Isogeometric basis properties

    % Polynomial Degree of the B-Spline Base
    
    degreeSplineBasis    = 1;
    
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
    case {1,2,3}
        bc_up={ BC_laterali_up };
        dato_up={0};
        bc_down={ BC_laterali_down};
        dato_down={0};
    end

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

    switch caso
    case 1
        a     = @(x)     0*x;     	% Equation of the Inferior Bound of the Channel (Considered As Central)
        L     = @(x)  1+ 0*x;     	% Thickness of the Channel
        psi_x = @(x)     0*x;    	% Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) 0*x;
        domainProfileDer = @(x) 0*x;

    case 2
        a     = @(x)     0*x;     	% Equation of the Inferior Bound of the Channel (Considered As Central)
        L     = @(x)  1+ 0*x;     	% Thickness of the Channel
        psi_x = @(x)     0*x;    	% Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) 0*x;
        domainProfileDer = @(x) 0*x;

    case 3
        a     = @(x)     0*x;     	% Equation of the Inferior Bound of the Channel (Considered As Central)
        L     = @(x)  1+ 0*x;     	% Thickness of the Channel
        psi_x = @(x)     0*x;    	% Rescaling

        geometricInfo = struct('a',a,'L',L,'psi_x',psi_x);

        % Assignment of the Profile and its Derivative:

        domainProfile  = @(x) 0*x;
        domainProfileDer = @(x) 0*x;
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

    knot      = augknt([min_x max_x],degreeSplineBasis+1);
    h         = stepHorMesh;
    nel       = numbElements;
    ins       = sort(reshape((h:h:1-h)'*ones(1,continuityParameter),1,continuityParameter*(nel-1)));
    ins       = (max_x-min_x)*ins+min_x;

    % Creates the control points and knots points of the desired
    % B-Spline curve

    [cp,knot] = bspkntins(degreeSplineBasis,min_x:(max_x-min_x)*1/degreeSplineBasis:max_x,knot,ins);

    mesh_fisx  = linspace(min_x,max_x,numbElements+1);
    mesh_fisx2 = linspace(min_x,max_x,2*numbControlPts);

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

    for ii = 1:2*numbControlPts-1

        pos1 = ii;
        pos2 = ii+1;

        % 1. Take 8 Quadrature Nodes and Respective Weights

        obj_gaussLegendre_2 = IntegrateHandler();

        obj_gaussLegendre_2.numbQuadNodes = num;

        [numnodes,x_nodes,wei] = gaussLegendre(obj_gaussLegendre_2); 

        %2. Rescale Them Depending on the Iteration

        obj_quadratureRule_2 = IntegrateHandler();

        obj_quadratureRule_2.leftBoundInterval = mesh_fisx2(ii);
        obj_quadratureRule_2.rightBoundInterval = mesh_fisx2(ii+1);
        obj_quadratureRule_2.inputNodes = x_nodes;
        obj_quadratureRule_2.inputWeights = wei;

        [nod,wei] = quadratureRule(obj_quadratureRule_2);

        %3. Definition of the Jacobian Computed for Every Node

        jaco = [];

        for it = 1:num

            jaco(it) = sqrt(1+(domainProfileDer(nod(it))).^2); 

        end

        %4. Computation of the Length of the Quadrature Rule

        %5. Save the Information in psJac

        psJac2 = [psJac2 sum(jaco.*wei')];

    end

    %% Coefficients of the bilinear form
    %-------------------------------------------------------------------------%

    switch caso
    case 1

        mu    = @(x,y) (  .1 + 0*x + 0*y ); % Difusion
        beta1 = @(x,y) (  5 + 0*x + 0*y ); % Horizontal Advection
        beta2 = @(x,y) (  0 + sin(6*x) + 0*y ); % Vertical Advection
        sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction

    case 2

        mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
        beta1 = @(x,y) (  50 + 0*x + 0*y ); % Horizontal Advection
        beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
        sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction

    case 3

        mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
        beta1 = @(x,y) (  0 + 0*x + 0*y ); % Horizontal Advection
        beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
        sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction

    end

    Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);

    %% Exact solution
    %-------------------------------------------------------------------------%
    L = 1;
    switch caso
    case 1

        true_sol  = @(x,y) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5));
        true_solx = @(x,y) (x.^4-x).*sin(2*pi*(y+0.5));
        true_soly = @(x,y) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5));

    case 2

        true_sol  = @(x,y) (1/4*x.^4-L^3*x).*sin(2*pi*(y+0.5));
        true_solx = @(x,y) (x.^3-L^3).*sin(2*pi*(y+0.5));
        true_soly = @(x,y) (1/4*x.^4-L^3*x)*2*pi.*cos(2*pi*(y+0.5)); 

    case 3

        true_sol  = @(x,y) (1/5*x.^5-L^4*x).*sin(2*pi*(y+0.5));
        true_solx = @(x,y) (x.^4-L^4).*sin(2*pi*(y+0.5));
        true_soly = @(x,y) (1/5*x.^5-L^4*x)*2*pi.*cos(2*pi*(y+0.5));     

    end

    %% Force and Dirichlet profile for the inflow data
    %-------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % NOTE:
    % This Dirichlet profile must be compatible witht the boundary conditions
    % of the problem.
    %-------------------------------------------------------------------------%

    switch caso
    case 1

        dato_dir = @(y) 0;
        
        force = @(x,y) 10*(( ((x-0.15).^2 + 4*(y-0.25).^2) < 0.01) + ( ((x-0.15).^2 + 4*(y+0.25).^2) < 0.01));
        
%         force = @(x,y) -(4*x.^3-1).*sin(2*pi*(y+0.5))...
%                        +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5))...
%                        +(x.^4-x).*sin(2*pi*(y+0.5));              

    case 2

        dato_dir = @(y) 0;
        len = L;
        
        force = @(x,y) 50*(( ((x-0.3).^2 + (y-0.25).^2) < 0.02) + ( ((x-0.3).^2 + (y+0.25).^2) < 0.01));
        
%         force = @(x,y)  (1/4*x.^4-len^3*x).*4*pi^2.*sin(2*pi*(y+0.5))+...
%                         (x.^3-len^3).*sin(2*pi*(y+0.5))...
%                         - 3*x.^2.*sin(2*pi*(y+0.5));     

    case 3

        dato_dir = @(y) 0;
        len = L;
        force = @(x,y) ((x-0.2).^2+(y).^2<0.01);
        force = @(x,y)  (1/5*x.^5-len^4*x).*4*pi^2.*sin(2*pi*(y+0.5))+...
                        (x.^4-len^4).*sin(2*pi*(y+0.5))...
                        - 4*x.^3.*sin(2*pi*(y+0.5));  

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
    obj_solverIGA.physicMesh_inX = mesh_fisx;
    obj_solverIGA.physicMesh_inY = mesh_fisy;
    obj_solverIGA.domainProfile = domainProfile;
    obj_solverIGA.domainProfileDer = domainProfileDer;
    obj_solverIGA.jacAtQuadNodes = psJac;
    obj_solverIGA.jacAtQuadNodes2 = psJac2; 
    obj_solverIGA.degreePolySplineBasis = degreeSplineBasis;
    obj_solverIGA.continuityParameter = continuityParameter;
    obj_solverIGA.numbHorQuadNodes = numbHorNodes;
    obj_solverIGA.numbVerQuadNodes = numbVerNodes;

    tic;
    [u2,a,b,L2_2,H1_2] = solverIGA(obj_solverIGA);
    toc;
    
    disp('Maximum L2 Norm Error');
    disp(max(max(L2_2)));
    
    disp('Maximum H1 Norm Error');
    disp(max(max(H1_2)));

%     % Error in the L2-Norm
% 
%     L2 = [L2 L2_2];
% 
%     % Error in the H1-Norm
% 
%     H1 = [H1 H1_2];
% 
%     ncpp = [ncpp ncp];
% 
%     display('The error in the desired norms are:')
%     display('L2 - Norm')
%     disp(L2')
%     display('H1 - Norm')
%     disp(H1')
% 
% 
%     ordL2=[];
%     ordH1=[];
% 
%     for kk = 1:nn-1
% 
%         % Convergence in the L2-Norm
% 
%         ordL2 = [ordL2 -log(L2(kk+1)/L2(kk))/log(2)];
% 
%         % Convergence in the H1-Norm
% 
%         ordH1 = [ordH1 -log(H1(kk+1)/H1(kk))/log(2)];
% 
%     end
% 
%     display('The orders of convergence are:')
%     display('Norma L2')
%     disp(ordL2')
%     display('Norma H1')
%     disp(ordH1')
% 
%     hh = [];
%     hh1 = [];
% 
%     for i=1:nn
% 
%         hh = [hh 1/(2^(i*(p+1)))];
%         hh1 = [hh1 1/(2^(i*(p)))];
% 
%     end

%     % Plots for the Order of Convergence
% 
%     figure;
% 
%     loglog(ncpp,hh,'r','LineWidth',2);
%     hold on
% 
%     loglog(ncpp,L2,'b-o','LineWidth',2);
%     xlabel('ndof','Fontsize',20);
%     ylabel('error L2','Fontsize',20);
%     legend Order3 Simulation
%     set(gca, 'FontSize', 20)
% 
%     figure;
% 
%     loglog(ncpp,hh1,'r','LineWidth',2);
%     hold on
% 
%     loglog(ncpp,H1,'b-o','LineWidth',2);
%     xlabel('ndof','Fontsize',20);
%     ylabel('error H1','Fontsize',20);
%     legend Order2 Simulation
%     set(gca, 'FontSize', 20)
% 
%     display('Tempi impiegati')
%     disp(times')