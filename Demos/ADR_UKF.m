%% HI-MOD MONODOMAIN
% The following script allows the solution of an Advection - Diffusion -
% Reaction differential problem in 2D using the HIGAMod solution.

    %% Demonstration Import Data
    
    close all;
    clear all;
    clc;

    disp('******************************************')
    disp('*           HIGAMod Simulation           *');
    disp('******************************************');

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    import Core.FilterADRHandler

    %% Definition of Domain Parameters
    
    simulationCase  = 9;
    filterTime = [];
    
% for ll = 1:4
%     for pp = 1:4
ll=2;
pp=3;

    m = 1;              % Number of Iterations for the Convergence Order

    errorNormL2 = [];
    errorNormH1 = [];
    ncpp = [];
    times = zeros(1,m);

    %-------------------------------------------------------------------------%
    % Note;
    % The exciting force acting uppon the system is different depending on the
    % case we are analysing. In the current version of the code we are
    % considering the following exciting forces:
    %
    %  (1) :              b = [1,0];      I = [0,1]; [NEU]
    %  (2) :  f = 1;      b = [0,0];      I = [0,1];
    %  (3) :  f = 0;      b = [1,0];      I = [0,1];
    %-------------------------------------------------------------------------%
    %

        %----------------------------------------------%
        %        ALSO WORKS FOR ONLY ONE DOMAIN        %
        %----------------------------------------------%
        % setting number of domains, mesh, dof & modes %
        %----------------------------------------------%

        switch simulationCase
            case {1}
                
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = 1;
                endTime    = 10;
                timeDomain = 0:timeStep:endTime;
                
                cutx      = [min_x,max_x];  
                cuty      = [min_y,max_y];
                m1        = [8 10];
                h1        = [ 0.2 0.1 0.05 0.025 0.0125];
                size_mb   = 3;
                nd        = length(size_mb); 
                hx        = (max_x-min_x)*0.05/(2^(m-1))*ones(size(size_mb));
                ne        = round((max_x-min_x)/hx);
        
                initialState = ones(size_mb * (ne + 1),1);
                
            case {2}
                
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = .01;
                endTime    = .1;
                timeDomain = 0:timeStep:endTime;
                
                cutx      = [min_x,max_x];  
                cuty      = [min_y,max_y];
                m1        = [8 10];
                h1        = [ 0.2 0.1 0.05 0.025 0.0125];
                size_mb   = 3;
                nd        = length(size_mb); 
                hx        = (max_x-min_x)*0.05/(2^(m-1))*ones(size(size_mb));
                ne        = round((max_x-min_x)/hx);
                
                initialState = ones(size_mb * (ne + 1),1);
                
            case {3,4}
                
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = .01;
                endTime    = .1;
                timeDomain = 0:timeStep:endTime;
                
                cutx      = [min_x,max_x];  
                cuty      = [min_y,max_y];
                m1        = [8 10];
                h1        = [ 0.2 0.1 0.05 0.025 0.0125];
                size_mb   = 3;
                nd        = length(size_mb); 
                hx        = (max_x-min_x)*0.05/(2^(m-1))*ones(size(size_mb));
                ne        = round((max_x-min_x)/hx);
                
                initialState = ones(size_mb * (ne + 1),1);
                
            case {8}
                
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = 1;
                endTime    = 10;
                timeDomain = 0:timeStep:endTime;
                
                cutx      = [min_x,max_x];  
                cuty      = [min_y,max_y];
                size_mb   = 3;
                nd        = length(size_mb); 
                hx        = (max_x-min_x)*0.1/(2^(m-1))*ones(size(size_mb));
                ne        = round((max_x-min_x)/hx);
                
                initialState = ones(size_mb * (ne + 1),1);
                
            case {9}
                
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = 1;
                endTime    = 10;
                timeDomain = 0:timeStep:endTime;
                
                cutx      = [min_x,max_x];  
                cuty      = [min_y,max_y];
                m1        = [2 3 4 5];
                h1        = [ 0.2 0.1 0.05 0.025 0.0125];
                size_mb   = m1(ll);
                nd        = length(size_mb); 
                hx        = (max_x-min_x)*h1(pp)/(2^(m-1))*ones(size(size_mb));
                ne        = round((max_x-min_x)/hx);
                
                initialState = 0 * ones(size_mb * (ne + 1),1);
                
            case {10}
                
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = 1;
                endTime    = 10;
                timeDomain = 0:timeStep:endTime;
                
                cutx      = [min_x,max_x];  
                cuty      = [min_y,max_y];
                m1        = [8 10];
                h1        = [ 0.2 0.1 0.05 0.025 0.0125];
                size_mb   = 5;
                nd        = length(size_mb); 
                hx        = (max_x-min_x)*0.05/(2^(m-1))*ones(size(size_mb));
                ne        = round((max_x-min_x)/hx);
                
                initialState = ones(size_mb * (ne + 1),1);
        end
        
    %% Definition of Measurement Parameters
    
        filterCase = 2;
    
        RSN_Measure = 10;
        RSN_Proccess = 1;
        
        measureCovOpt = [0 -1 -2 -3];
        
        outputMatrix = eye(size_mb * (ne + 1));
        
        proccessCovCoeff = -2;
        measureCovCoeff = measureCovOpt(filterCase);

    %% Definition of Isogeometric Parameters
        
        % IGA Parameters

        p         = 1;                     % Polynomial Degree of the B-Spline Base
        k         = 1;                     % Continuity of the Base 'C^(p-k)'
        ncp       = ne * k + p + 1 -k;     % Number of Control Points

        if (p<k)
            error('Wrong Choice for the Parameters!');
        end
        
    %% Definition of Boundary Conditions

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
        chi           = 3;
        mu            = 1;
        cest          = 0.05;

        switch simulationCase
            case {1,2,3,4,5,6,7,8}
                
                BC_laterali_up   ='rob';
                BC_laterali_down ='rob';
                
                bc_up={ BC_laterali_up };
                dato_up={0};
                bc_down={ BC_laterali_down};
                dato_down={0};
                
            case {9}
                
                BC_laterali_up   ='dir';
                BC_laterali_down ='dir';
                
                bc_up={ BC_laterali_up };
                dato_up={0};
                bc_down={ BC_laterali_down};
                dato_down={0};
                
            case {10}
                
                BC_laterali_up   ='rob';
                BC_laterali_down ='rob';
                
                bc_up={ BC_laterali_up };
                dato_up={0};
                bc_down={ BC_laterali_down};
                dato_down={0};
        end
        
    %% Definition of Channel Properties

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
            case {1,2,3,4,5,6,7}
                a     = @(x)     0*x;     	% Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x)  1+ 0*x;     	% Thickness of the Channel
                psi_x = @(x)     0.*x;    	% Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) 0*x;
                dforma = @(x) 0*x;

            case 9
                a     = @(x) sqrt(1 - x.^2); % Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x)  1+ 0*x;     	 % Thickness of the Channel
                psi_x = @(x)     0.*x;    	 % Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) sqrt(1 - x.^2);
                dforma = @(x) -x./sqrt(1 - x.^2);  

            case 8
                a     = @(x) sqrt(1 - (x - 1).^2);  % Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x) 1+ 0*x;                % Thickness of the Channel
                psi_x = @(x) 0.*x;                  % Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) sqrt(1 - (x - 1).^2);
                dforma = @(x) -(x-1)./sqrt(1 - (x - 1).^2);
                
            case 10
                a     = @(x) sqrt(1 - (x - 1).^2);  % Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x) 1+ 0*x;                % Thickness of the Channel
                psi_x = @(x) 0.*x;                  % Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) -sqrt(1 - (x - 1).^2);
                dforma = @(x) (x-1)./sqrt(1 - (x - 1).^2);
        end

    %% Definition of Horizontal Mesh

        mesh_fisy = [];
        knot      = augknt([min_x max_x],p+1);
        h         = hx;
        nel       = ne;
        ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
        ins       = (max_x-min_x)*ins+min_x;
        [cp,knot] = bspkntins(p,min_x:(max_x-min_x)*1/p:max_x,knot,ins);

        mesh_fisx  = linspace(min_x,max_x,ne+1);
        mesh_fisx2 = linspace(min_x,max_x,2*ncp);

        psJac  = [];
        psJac2 = [];
        num    = 8;

        for ii = 1:ne

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

                jaco(it) = sqrt(1+(dforma(nod(it))).^2); 

            end

            % 4. Computation of the Length of the Quadrature Rule

            % 5. Save the Information in psJac

            psJac = [psJac sum(jaco.*wei')];

        end

        for ii = 1:2*ncp-1

            pos1 = ii;
            pos2 = ii+1;

            % 1. Take 8 Quadrature Nodes and Respective Weights

            obj_gaussLegendre_2 = IntegrateHandler();

            obj_gaussLegendre_2.numbQuadNodes = num;

            [numnodes,x_nodes,wei] = gaussLegendre(obj_gaussLegendre_2); 

            %2. Li riscalo a seconda della iterazione

            obj_quadratureRule_2 = IntegrateHandler();

            obj_quadratureRule_2.leftBoundInterval = mesh_fisx2(ii);
            obj_quadratureRule_2.rightBoundInterval = mesh_fisx2(ii+1);
            obj_quadratureRule_2.inputNodes = x_nodes;
            obj_quadratureRule_2.inputWeights = wei;

            [nod,wei] = quadratureRule(obj_quadratureRule_2);

            %3. Definizione dello jacobiano calcolato nei num nodi

            jaco = [];

            for it = 1:num

                jaco(it) = sqrt(1+(dforma(nod(it))).^2); 

            end

            %4. Calcolo lunghezza con formula di quadratura

            %5. Salvo in psJac

            psJac2 = [psJac2 sum(jaco.*wei')];

        end

    %% Definition of Bilinear Form Coefficients
        
        %-------------------------------------------------------------------------%
        %                   COEFFICIENTS OF THE BILINEAR FORM
        %-------------------------------------------------------------------------%

        switch simulationCase
            case {1}

                mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  1 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction

            case {2}

                mu    = @(x,y) (  0.5 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  1 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  1 + 0*x + 0*y ); % Reaction

            case {3}

                mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  1 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction
                
            case {4}

                mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  0 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction
                
            case {8}
                
                mu    = @(x,y) (  0.5 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  1.0 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  1.5 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  2.0 + 0*x + 0*y ); % Reaction
                
            case {9}
                
                mu    = @(x,y) (  0.05 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  0.1 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0.0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0.0 + 0*x + 0*y ); % Reaction
                
            case {10}
                
                mu    = @(x,y) (  1.0 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  0.5 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  1.5 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  2.0 + 0*x + 0*y ); % Reaction

        end

        Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
            'sigma',sigma,'coeffrobin',chi);

    %% Definition of Exact Solution for Testing
        
        %-------------------------------------------------------------------------%
        %                   EXACT SOLUTION FOR THE CASES PROPOSED
        %-------------------------------------------------------------------------%

        switch simulationCase
            case 1

                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*exp(0.1*t);

            case 2

                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*0.05*t;
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*0.05*t;
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*0.05*t;
                
            case 3

                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);

            case 4

                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);
                
            case {5,6,7,8,9,10}
                
                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2);

        end

    %% Definition of Dirichlet Profile at Inflow
        
        %-------------------------------------------------------------------------%
        %           CREATION OF THE DIRICHLET PROFILE FOR THE INFLOW DATA
        %-------------------------------------------------------------------------%

        %-------------------------------------------------------------------------%
        % NOTE:
        % This Dirichlet profile must be compatible witht the boundary conditions
        % of the problem.
        %-------------------------------------------------------------------------%

        switch simulationCase
            case 1 

                dato_dir = @(y) 0;
                force = @(x,y,t) 10 * ((x - sqrt(2)/2).^2 + (y - (sqrt(2) + 1)/2).^2 < 0.01);
%                 force = @(x,y,t) -(4*x.^3-1).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
%                                  +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5)).*exp(0.1*t)...
%                                  +(x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
%                                  +(0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t)*0.1;

            case 2

                dato_dir = @(y) 0;
                
                force = @(x,y,t) 10;
                
%                 force = @(x,y,t) -(4*x.^3-1).*sin(2*pi*(y+0.5)).*0.05*t...
%                                  +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5)).*0.05*t...
%                                  +(x.^4-x).*sin(2*pi*(y+0.5)).*0.05*t...
%                                  +(0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*0.05;
 
            case 3

                dato_dir = @(y) 0;
                force = @(x,y,t) 10;
%                 force = @(x,y,t) -(4*x.^3-1).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2)...
%                                  +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2)...
%                                  +(x.^4-x).*sin(2*pi*(y+0.5)).*(sin(0.1*pi*t) + 2)...
%                                  +(0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*cos(0.1*pi*t)*0.1*pi;
                             
            case 4

                dato_dir = @(y) 0;
                force = @(x,y,t) 2;
                % force = @(x,y,t) 1000 * ((x-0.5).^2+(y).^2<0.01);
                
            case 8
                
                dato_dir = @(y) 0;
                force = @(x,y,t) 10;
                
            case 9
                
                dato_dir = @(y) 0;
                force = @(x,y,t) 1;
                
            case 10
                
                dato_dir = @(y) 0;
                force = @(x,y,t) x .* y;
        end

        Dati = struct('dato_dir',dato_dir,'force',force);
        
    %% Definition of Filter Parameters
        
    %% Equation Solver

        %-------------------------------------------------------------------------%
        % Note;
        % The following loop varies in order to change the coefficients of the
        % interface and to try different configurations at the same time.
        %-------------------------------------------------------------------------%

        %--------%
        % SOLVER %
        %--------%

        % Definition of the Object of the EvaluationHandler Class

        import Core.SolverHandler

        obj_solverUKF = SolverHandler();

        % Properties Assignment

        obj_solverUKF.domainLimit_inX = cutx;
        obj_solverUKF.domainLimit_inY = cuty;
        obj_solverUKF.dimModalBasis = size_mb;
        obj_solverUKF.stepMeshX = hx;
        obj_solverUKF.label_upBoundDomain = bc_up;
        obj_solverUKF.label_downBoundDomain = bc_down;
        obj_solverUKF.data_upBoundDomain = dato_up;
        obj_solverUKF.data_downBoundDomain = dato_down;
        obj_solverUKF.dirCondFuncStruct = Dati;
        obj_solverUKF.coefficientForm = Coeff_forma;
        obj_solverUKF.geometricInfo = Dati_geometrici;
        obj_solverUKF.dataExportOption = true;
        obj_solverUKF.couplingCond_DD = 'RR';
        obj_solverUKF.simulationCase = simulationCase;
        obj_solverUKF.exactSolution = true_sol;
        obj_solverUKF.exactSolution_dX = true_solx;
        obj_solverUKF.exactSolution_dY = true_soly;
        obj_solverUKF.physicMesh_inX = mesh_fisx;
        obj_solverUKF.physicMesh_inY = mesh_fisy;
        obj_solverUKF.domainProfile = forma;
        obj_solverUKF.domainProfileDer = dforma;
        obj_solverUKF.jacAtQuadNodes = psJac;
        obj_solverUKF.jacAtQuadNodes2 = psJac2; 
        obj_solverUKF.degreePolySplineBasis = p;
        obj_solverUKF.continuityParameter = k;
        obj_solverUKF.numbControlPoints = ncp; 
        obj_solverUKF.numbVerQuadNodes = 64;

        obj_solverUKF.timeStep = timeStep;
        obj_solverUKF.timeDomain = timeDomain;
        obj_solverUKF.initialState = initialState;

        obj_solverUKF.RSN_Measure = RSN_Measure;
        obj_solverUKF.RSN_Proccess = RSN_Proccess;
        obj_solverUKF.measureMatrix = outputMatrix;
        obj_solverUKF.proccessCov = proccessCovCoeff;
        obj_solverUKF.measureCov = measureCovCoeff;

        obj_solverUKF.caseSimul = simulationCase;
        obj_solverUKF.caseFilter = size_mb;               

        [solutionMatrix,noisyMeasure,estimateMatrix,fTime] = solverUKF(obj_solverUKF);
        
        filterTime(ll,pp) = fTime;
        
%     end
% end

fid = fopen('FilterTimeMH.txt','wt');
for ii = 1:size(filterTime,1)
    fprintf(fid,'%0.15f \t',filterTime(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
    
%     fid = fopen(['Cond',num2str(ll+3),num2str(pp),'.txt'],'wt');
%     for ii = 1:size(condVectMatrixObs,1)
%     fprintf(fid,'%g\t',condVectMatrixObs(ii,:));
%     fprintf(fid,'\n');
%     fclose(fid)
    
    %% Error Plots
    
%     figure;
%     for kk = 1:length(numbEnsembles)
%         plot(log(condVectMatrixState(kk,:)));
%         hold on
%     end
%     
%     figure;
%     for kk = 1:length(numbEnsembles)-1
%         loglog(numbEnsembles,condVectMatrixObs(:,kk));
%         hold on
%     end
%     
%     % Export Condition Number Data
%     dlmwrite('numbEnsembles.txt',numbEnsembles);
%     dlmwrite('condVectMatrixState.txt',condVectMatrixState);
%     dlmwrite('condVectMatrixObs.txt',condVectMatrixObs);
%     fclose(fileID);

%     figure;
%     
%     time = timeDomain(2:end);
%     stem(time,errorNormL2);hold on;
%     stem(time,errorNoiseL2);hold on;
%     stem(time,errorEstimateL2);hold off;
%     titleName = sprintf('Error L2-Norm Case: %d Filter: %d',simulationCase,mcase);
%     title(titleName)
%     
%     figure;
%     
%     stem(time,errorNormH1); hold on;
%     stem(time,errorNoiseH1); hold on;
%     stem(time,errorEstimateH1); hold off;
%     titleName = sprintf('Error H1-Norm Case: %d Filter: %d',simulationCase,mcase);
%     title(titleName)