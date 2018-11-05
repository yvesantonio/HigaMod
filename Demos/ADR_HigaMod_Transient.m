%% HI-MOD MONODOMAIN
% The following script allows the solution of an Advection - Diffusion -
% Reaction differential problem in 2D using the HIGAMod solution.

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
    
% for jj = 1:4

    nn    = 1;              % Number of Iterations for the Convergence Order
    simulationCase  = 4;    % Analysed Case

    errorNormL2 = [];
    errorNormH1 = [];
    ncpp = [];
    times = zeros(1,nn);

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

    for m = 1:nn

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
                
            case {2}
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = 0.031;
                endTime    = 1;
                timeDomain = 0:timeStep:endTime;
                
            case {3,4}
                min_x     = 0;
                max_x     = 1;
                min_y     = 0;
                max_y     = 1;
                timeStep   = 0.01;
                endTime    = 0.1;
                timeDomain = 0:timeStep:endTime;
        end

        %-------------------------------------------------------------------------%
        % NOTE:
        % If one wants to change the values of min_x and max_x, the procedure can
        % be performed inside the function 'assembla_IGA'.
        %-------------------------------------------------------------------------%

        cutx      = [min_x,max_x];  
        cuty      = [min_y,max_y];
        m1        = 5;
        size_mb   = m1;
        nd        = length(size_mb); 
        hx        = (max_x-min_x)*0.1/(2^(m-1))*ones(size(size_mb));
        ne        = round((max_x-min_x)/hx);

        % IGA Parameters

        p         = 1;                     % Polynomial Degree of the B-Spline Base
        k         = 1;                     % Continuity of the Base 'C^(p-k)'
        ncp       = ne * k + p + 1 -k;     % Number of Control Points
        
        initialState = zeros(m1 * ncp,1);

        if (p<k)
            error('Wrong Choice for the Parameters!');
        end

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

        BC_laterali_up   ='dir';
        BC_laterali_down ='dir';

        switch simulationCase
            case {1,2,3,4,5,6,7,8,9,10}
                bc_up={ BC_laterali_up };
                dato_up={0};
                bc_down={ BC_laterali_down};
                dato_down={0};
        end


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
            case {1,2,3,5,6}
                b = 0.5;
                n = 0.5;
                
                a     = @(x)     -x.^(1/2);     	% Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x)  1+ 0*x;     	% Thickness of the Channel
                psi_x = @(x)     0;    	% Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) -x.^(1/2);
                dforma = @(x) (1/2)*(-x.^(-1/2));

            case {7,10}
                cof   = 1;
                a     = @(x)     cof*x.^2;      % Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x)  1+ 0*x;           % Thickness of the Channel
                psi_x = @(x)     0.*x;          % Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) cof*x.^2;
                dforma = @(x) 2*cof*x;  

            case 8
                a     = @(x)     x.^3;     	% Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x)  1+ 0*x;     	% Thickness of the Channel
                psi_x = @(x)     0.*x;    	% Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) x.^3;
                dforma = @(x) 3*x.^2;  

            case 9
                a     = @(x)     3*x;     	% Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x)  1+ 0*x;     	% Thickness of the Channel
                psi_x = @(x)     0.*x;    	% Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) 3*x;
                dforma = @(x) 3;    
                
            case {4}
                b = 0.5;
                n = 0.5;
                
                a     = @(x)     0*x;     	% Equation of the Inferior Bound of the Channel (Considered As Central)
                L     = @(x)  1+ 0*x;     	% Thickness of the Channel
                psi_x = @(x)     0;    	% Rescaling

                Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);

                % Assignment of the Profile and its Derivative:

                forma  = @(x) 0*x;
                dforma = @(x) 0*x;
        end

        % For the X Axis

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

        %-------------------------------------------------------------------------%
        %                   COEFFICIENTS OF THE BILINEAR FORM
        %-------------------------------------------------------------------------%

        switch simulationCase
            case {1}

                mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  1 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction

            case {2,4}

                mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  1 + 0*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 0*x + 0*y ); % Vertical Advection
                sigma = @(x,y) (  0 + 0*x + 0*y ); % Reaction

            case {3}

                mu    = @(x,y) (  1 + 0*x + 0*y ); % Difusion
                beta1 = @(x,y) (  4 + 2*x + 0*y ); % Horizontal Advection
                beta2 = @(x,y) (  0 + 6*cos(2*pi*x/6) + cos(y) ); % Vertical Advection
                sigma = @(x,y) (  1 + 0*x + 0*y ); % Reaction

        end

        Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
            'sigma',sigma,'coeffrobin',chi);

        %-------------------------------------------------------------------------%
        %                   EXACT SOLUTION FOR THE CASES PROPOSED
        %-------------------------------------------------------------------------%

        switch simulationCase
            case 1

                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*exp(0.1*t);

            case 2

                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*exp(0.1*t);
                
            case {3,4}

                true_sol  = @(x,y,t) (0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_solx = @(x,y,t) (x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t);
                true_soly = @(x,y,t) (0.2*x.^5-0.5*x.^2)*2*pi.*cos(2*pi*(y+0.5)).*exp(0.1*t);

        end

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
                force = @(x,y,t) -(4*x.^3-1).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t)*0.1;

            case 2

                dato_dir = @(y) 0;
                
                force = @(x,y,t) -(4*x.^3-1).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t)*0.1;
 
            case 3

                dato_dir = @(y) 0;
                force = @(x,y,t) (1 + 5*t).*((x-0.2).^2+(y).^2<0.01);
                
            case 4

                dato_dir = @(y) 0;
                force = @(x,y,t) -(4*x.^3-1).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(0.2*x.^5-0.5*x.^2)*4*pi^2.*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(x.^4-x).*sin(2*pi*(y+0.5)).*exp(0.1*t)...
                                 +(0.2*x.^5-0.5*x.^2).*sin(2*pi*(y+0.5)).*exp(0.1*t)*0.1;
        end

        Dati = struct('dato_dir',dato_dir,'force',force);

        %-------------------------------------------------------------------------%
        % Note;
        % The following loop varies in order to change the coefficients of the
        % interface and to try different configurations at the same time.
        %-------------------------------------------------------------------------%

        for j = 25:5:25

            for i = -1:-2:-1

                gammaL = j;
                gammaR = i;        
                gamma = struct ('L',gammaL,'R',gammaR);

                %--------%
                % SOLVER %
                %--------%

                % Definition of the Object of the EvaluationHandler Class

                import Core.SolverHandler

                obj_solverTransientIGA = SolverHandler();

                % Properties Assignment

                obj_solverTransientIGA.domainLimit_inX = cutx;
                obj_solverTransientIGA.domainLimit_inY = cuty;
                obj_solverTransientIGA.dimModalBasis = size_mb;
                obj_solverTransientIGA.stepMeshX = hx;
                obj_solverTransientIGA.label_upBoundDomain = bc_up;
                obj_solverTransientIGA.label_downBoundDomain = bc_down;
                obj_solverTransientIGA.data_upBoundDomain = dato_up;
                obj_solverTransientIGA.data_downBoundDomain = dato_down;
                obj_solverTransientIGA.dirCondFuncStruct = Dati;
                obj_solverTransientIGA.coefficientForm = Coeff_forma;
                obj_solverTransientIGA.geometricInfo = Dati_geometrici;
                obj_solverTransientIGA.robinCondStruct = gamma;
                obj_solverTransientIGA.dataExportOption = true;
                obj_solverTransientIGA.couplingCond_DD = 'RR';
                obj_solverTransientIGA.simulationCase = simulationCase;
                obj_solverTransientIGA.exactSolution = true_sol;
                obj_solverTransientIGA.exactSolution_dX = true_solx;
                obj_solverTransientIGA.exactSolution_dY = true_soly;
                obj_solverTransientIGA.physicMesh_inX = mesh_fisx;
                obj_solverTransientIGA.physicMesh_inY = mesh_fisy;
                obj_solverTransientIGA.domainProfile = forma;
                obj_solverTransientIGA.domainProfileDer = dforma;
                obj_solverTransientIGA.jacAtQuadNodes = psJac;
                obj_solverTransientIGA.jacAtQuadNodes2 = psJac2; 
                obj_solverTransientIGA.degreePolySplineBasis = p;
                obj_solverTransientIGA.continuityParameter = k;  
                obj_solverTransientIGA.numbControlPoints = ncp;  
                
                obj_solverTransientIGA.timeStep = timeStep;
                obj_solverTransientIGA.timeDomain = timeDomain;
                obj_solverTransientIGA.initialState = initialState;
                

                tic;
                [u2,a,b,errL2,errH1,stateMatrix,inputMatrix,forceHistory] = solverIGATransient(obj_solverTransientIGA);
                times(m) = toc;

            end
        end
    end
% end

%     figure;
%     
%     title('Error Evolution on the L2-Norm')
%     plot(errorNormL2);
%     
%     figure;
%     
%     title('Error Evolution on the H1-Norm')
%     plot(errorNormH1);