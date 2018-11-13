classdef SolverHandler
    
    %% SOLVER HANDLER CLASS
    % The SolverHandler is a class that contains all of the scripts
    % responsible for the computation of the solution for the differential
    % problems proposed in the examples. It calls all of the other basic
    % classes and functions already incorporated.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% SOLVER HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the SolverHandler
        % encapsulate all of the variables needed to run the methods 
        % associated with the computation of the solution for the
        % differential problems proposed. They basically represent all of
        % the variables required to define complitily the equations, domain
        % and boundary conditions of the differential problem.
        
            %% Solver Properties

            domainLimit_inX;        % Vector Containing the Extremes of the Domains
                                    % in the X Direction
                          
            domainLimit_inY;        % Vector Containing the Extremes of the Domains
                                    % in the Y Direction
                          
            dimModalBasis;          % Dimension of the Modal Basis in Each Domain
            
            stepMeshX;              % Vector Containing the Step of the Finite
                                    % Element Mesh
                          
            label_upBoundDomain;    % Contains the Label Identifying the Nature of
                                    % the Boundary Conditions on the Upper Limit of
                                    % the Domain
                          
            label_downBoundDomain;  % Contains the Label Identifying the Nature of
                                    % the Boundary Conditions on the Lower Limit of
                                    % the Domain
                          
            data_upBoundDomain;     % Contains the Values of the Boundary Conditions
                                    % on the Upper Limit of the Domain
                          
            data_downBoundDomain;   % Contains the Values of the Boundary Conditions
                                    % on the Lower Limit of the Domain
                          
            coefficientForm;        % Data Strusture Containing All the @-Functions
                                    % and the Constants Relative to the Bilinear Form
                          
            dirCondFuncStruct;      % Data Structure Containing All the @-Functions
                                    % for the Dirichlet Conditions at the Inflow and
                                    % for the Exciting Force
                          
            geometricInfo;          % Data Structure Containing All the
                                    % Geometric Information regarding the
                                    % Domain. The current version of the code
                                    % works only for the specific condition of:
                                    % (L = 1, a = 0, psi_x = 0)
                                
            robinCondStruct;        % Data Structure Containing the Two Values of the
                                    % Coefficients (R, L) for the Robin Condition Used
                                    % in the Domain Decomposition
                          
            dataExportOption;       % Labels the Kind of Plot Function that Will Be
                                    % Used Once the Solution is Compluted
                          
            couplingCond_DD;        % Contains the Label Adressing the Coupling Condition
                                    % of the Problem
                          
            simulationCase;         % Specify What Problem Case is Being Solved
            
            exactSolution;          % Exact solution of the differential problem
            
            exactSolution_dX;       % Derivative of the exact solution in the X direction
            
            exactSolution_dY;       % Derivative of the exact solution in the Y direction
            
            physicMesh_inX;         % Vector Containing the Physical Mesh in the X
                                    % Direction
                          
            physicMesh_inY;         % Vector Containing the Physical Mesh in the Y
                                    % Direction
                          
            jacAtQuadNodes;         % Data Sturcture Used to Save the Value of the
                                    % Jacobians Computed in the Quadratures Nodes 
                          
            jacAtQuadNodes2;        % Data Sturcture Used to Save the Value of the
                                    % Jacobians Computed in the Quadratures Nodes 
                          
            degreePolySplineBasis;  % Degree of the Polynomial B-Spline Basis
            
            continuityParameter;    % Degree of Continuity of the Basis 'C^(p-k)'
            
            domainProfile;          % Symbolic Function Defining the Profile of the
                                    % Simulation Domain
                          
            domainProfileDer;       % Symbolic Function Defining the Derivative of
                                    % the Profile of the Simulation Domain
                          
            domainMesh_inX;         % Mesh of the Problem in the X Direction
            
            exciteForce;            % Force Exciting the System in the X Direction
            
            nodeMapping;            % Data Structure Containing the Map of the Nodes
            
            initialState            % Initial state of the solution for the ADR problem
            
            timeStep;               % Time step applied to the Theta Method in the 
                                    % time discretization
            
            timeDomain;             % Vector containing the time instant in which the
                                    % solution will be computed
                                    
            RSN_Measure;       % Signal to noise ratio considering the spectral
                               % power of both signals
                               
            RSN_Proccess;      % Signal to noise ratio considering the spectral
                               % power of both signals
            
            measureMatrix;
            
            proccessCov;
            
            measureCov;
            
            caseSimul;
            
            caseFilter;
            
            numbControlPoints;
            
            delta;                      %Parameter to tune for stabilization (in this case RTS)
            
            numbHorQuadNodes;           % Number of horizontal nodes to apply the quadrature
                                        % formula
                                    
            numbVerQuadNodes;           % Number of vertical nodes to apply the quadrature
                                        % formula
                                        
            D1; % First component of the Jacobian contribution to the
                % computation of the r** components for integration along
                % the vertical direction
                
            D2; % Second component of the Jacobian contribution to the
                % computation of the r** components for integration along
                % the vertical direction
                
            numbEnsembles   % Number of ensembles required to perform the ensemble Kalman
                            % Filter
                            
            numbParam       % Number of parameters to be identified
            
            centreline      % NURBS curve defining the centreline of 
                            % the physical domain
                
            stabMethod      % Stabilization method to be used in the solution
            
            stabDelta       % Scalling factor that may be used to control the
                            % amount of stabilization to be added in the
                            % problem
                            
            PeConfig        % Choice of the definition of the Pechlet number 
                            % to be used
                            
            selectedPlots   % Vector containing the plots to be saved at the 
                            % end of the simulation
                            
            refSolStruct    % Structure containing all the information of 
                            % refenrece solution
                            
            bcStruct
    end
    
    methods (Access = public)
        
        %% SOLVER HANDLER - CONSTRUCT METHOD
        % The construct method is the basic structure that defines a new
        % object uppon creation. It contains the basic properties that need
        % to be assigned when an object of the class AssemblerADRHandler is
        % created.
        % 
        % By default, MATLAB generates a default constructor with no input 
        % arguments. We chose to leave the default definition of the
        % constructor because in this way we can freely defines the
        % properties we need each time a method is used. 
        %
        % Since all the properties are public, they can me changed from the
        % outside of the class.
        
        % function obj = EvaluationHandler(varagin)
		% obj.options = gatherUserOptions(obj.options, varagin{:});		
        % end
        
        %% SOLVER HANDLER - FUNCTIONAL METHODS
        % The following functions correspond to the individual scripts
        % previously developed of the functions 'solverIGA_DD',
        % 'solver_DD' and 'solverLegendre'.
        
            %% Method 'solverFEM'
            
            function [u,liftCoeffA,liftCoeffB,errorNormL2,errorNormH1] = solverFEM(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%

                % Initialization of the Data Structure that will further contain the 
                % Labels relative to the Nature of the Boundary Conditions at the 
                % Interface

                BC_l = cell(2,1);
                BC_r = cell(2,1);
                
                BC_l{1}  = 'dir';   % Inflow Boundary Condition
                BC_r{2} = 'neu';   % Outflow Boundary Condition


                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%

                % Definition of the Object from the AssemblerADRHandler class

                build_IGA = AssemblerADRHandler();

                % Properties Assignment

                build_IGA.dimModalBasis = obj.dimModalBasis(1);
                build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                build_IGA.stepMeshX = obj.stepMeshX(1);
                build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                build_IGA.coefficientForm = obj.coefficientForm;
                build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                build_IGA.geometricInfo = obj.geometricInfo;
                build_IGA.robinCondStruct = obj.robinCondStruct;
                build_IGA.couplingCond_DD = obj.couplingCond_DD;
                build_IGA.physicMesh_inX = obj.physicMesh_inX;
                build_IGA.physicMesh_inY = obj.physicMesh_inY;
                build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                build_IGA.continuityParameter = obj.continuityParameter;
                build_IGA.domainProfile = obj.domainProfile;
                build_IGA.domainProfileDer = obj.domainProfileDer;
                build_IGA.numbHorQuadNodes = obj.numbHorQuadNodes;
                build_IGA.numbVerQuadNodes = obj.numbVerQuadNodes;

                % Call of the 'buildSystemIGA' Method

                [A,b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildSystemFEM(build_IGA); 

                disp('Finished BUILD SYSTEM IGA');                
                
                %-----------------------------------------------------------------%
                %                  	      PROBLEM SOLUTION                        %
                %-----------------------------------------------------------------%
                
                u = A\b;

                %-----------------------------------------------------------------%
                %                  	       SOLUTION PLOT                          %
                %-----------------------------------------------------------------%

                [errL2,errH1,~,~,~,~] = plot_solution_FEM( ...
                obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer,obj.jacAtQuadNodes, ...
                obj.jacAtQuadNodes2,obj.degreePolySplineBasis,obj.continuityParameter,obj.exactSolution_dX,obj.exactSolution_dY,obj.geometricInfo);

                % PREPARATION OF THE DATA FOR THE FUNCTION THAT COMPUTES THE ERROR

                errorNormH1 = errH1;
                errorNormL2 = errL2;

                disp('Finished Method SOLVER IGA');

            end % End function
            
            %% Method 'solverIGA'
            
            function [u,liftCoeffA,liftCoeffB,errorNormL2,errorNormH1] = solverIGA(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%

                % Initialization of the Data Structure that will further contain the 
                % Labels relative to the Nature of the Boundary Conditions at the 
                % Interface

                BC_l = cell(2,1);
                BC_r = cell(2,1);
                
                BC_l{1}  = 'dir';   % Inflow Boundary Condition
                BC_r{2} = 'neu';   % Outflow Boundary Condition


                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%

                % Definition of the Object from the AssemblerADRHandler class

                build_IGA = AssemblerADRHandler();

                % Properties Assignment

                build_IGA.dimModalBasis = obj.dimModalBasis(1);
                build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                build_IGA.stepMeshX = obj.stepMeshX(1);
                build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                build_IGA.coefficientForm = obj.coefficientForm;
                build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                build_IGA.geometricInfo = obj.geometricInfo;
                build_IGA.robinCondStruct = obj.robinCondStruct;
                build_IGA.couplingCond_DD = obj.couplingCond_DD;
                build_IGA.physicMesh_inX = obj.physicMesh_inX;
                build_IGA.physicMesh_inY = obj.physicMesh_inY;
                build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                build_IGA.continuityParameter = obj.continuityParameter;
                build_IGA.domainProfile = obj.domainProfile;
                build_IGA.domainProfileDer = obj.domainProfileDer;
                build_IGA.numbHorQuadNodes = obj.numbHorQuadNodes;
                build_IGA.numbVerQuadNodes = obj.numbVerQuadNodes;
                build_IGA.D1 = obj.D1;
                build_IGA.D2 = obj.D2;

                % Call of the 'buildSystemIGA' Method

                [A,b,~,liftCoeffA(1),liftCoeffB(1),Jac] = buildSystemIGA(build_IGA); 

                disp('Finished BUILD SYSTEM IGA');                
                
                %-----------------------------------------------------------------%
                %                  	      PROBLEM SOLUTION                        %
                %-----------------------------------------------------------------%
                
                % Debug 
                %---------------------------------------------------------------------%
                % plot(linspace(0,length(b),length(b)),b);
                %---------------------------------------------------------------------%

                u = A\b;
                
                disp('Finished SOLVE SYSTEM');

                %-----------------------------------------------------------------%
                %                  	       SOLUTION PLOT                          %
                %-----------------------------------------------------------------%

                [errL2,errH1,~,~,~] = plot_solution_IGA( ...
                obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer,obj.jacAtQuadNodes, ...
                obj.degreePolySplineBasis,obj.continuityParameter,obj.exactSolution_dX,obj.exactSolution_dY,...
                obj.domainLimit_inY);

                disp('Finished PLOT OPERATION / ERROR with EXACT SOLUTION')
            
%                 [errL2,errH1] = computeErrorIGA( ...
%                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
%                 obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer,obj.jacAtQuadNodes, ...
%                 obj.degreePolySplineBasis,obj.continuityParameter,obj.exactSolution_dX,obj.exactSolution_dY,...
%                 obj.domainLimit_inY);

                [errL2,errH1] = computeExactErrorIGA( ...
                obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer,obj.jacAtQuadNodes, ...
                obj.degreePolySplineBasis,obj.continuityParameter,obj.exactSolution_dX,obj.exactSolution_dY,...
                obj.domainLimit_inY);
             
                disp('Finished ERROR with FREEFEM++ SOLUTION')
            
%                 export_py(obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
%                 obj.coefficientForm,['matlab_solver'],1,0);

                % PREPARATION OF THE DATA FOR THE FUNCTION THAT COMPUTES THE ERROR

                errorNormH1 = errH1;
                errorNormL2 = errL2;

                disp('Finished Method SOLVER IGA');

            end % End function
            
            %% Method 'solverIGAScatter'
            
            function [u,liftCoeffA,liftCoeffB,errorNormL2,errorNormH1] = solverIGAScatter(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%

                % Definition of the Object from the AssemblerADRHandler class

                build_IGA = AssemblerADRHandler();

                % Properties Assignment

                build_IGA.dimModalBasis = obj.dimModalBasis(1);
                build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                build_IGA.stepMeshX = obj.stepMeshX(1);
                build_IGA.label_upBoundDomain = obj.dirCondFuncStruct.igaBoundCond.BC_UP_TAG;
                build_IGA.label_downBoundDomain = obj.dirCondFuncStruct.igaBoundCond.BC_DOWN_TAG;
                build_IGA.localdata_upBDomain = obj.dirCondFuncStruct.igaBoundCond.BC_UP_DATA;
                build_IGA.localdata_downBDomain = obj.dirCondFuncStruct.igaBoundCond.BC_DOWN_DATA;
                build_IGA.coefficientForm = obj.coefficientForm;
                build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                build_IGA.geometricInfo = obj.geometricInfo;
                build_IGA.robinCondStruct = obj.robinCondStruct;
                build_IGA.couplingCond_DD = obj.couplingCond_DD;
                build_IGA.physicMesh_inX = obj.physicMesh_inX;
                build_IGA.physicMesh_inY = obj.physicMesh_inY;
                build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                build_IGA.continuityParameter = obj.continuityParameter;
                build_IGA.domainProfile = obj.domainProfile;
                build_IGA.domainProfileDer = obj.domainProfileDer;
                build_IGA.numbHorQuadNodes = obj.numbHorQuadNodes;
                build_IGA.numbVerQuadNodes = obj.numbVerQuadNodes;
                build_IGA.igaBoundCond = obj.dirCondFuncStruct.igaBoundCond;
                
                % Call of the 'buildSystemIGA' Method

                [A,b,~,liftCoeffA(1),liftCoeffB(1),space,refDomain1D,boundStruct] = buildSystemIGAScatter(build_IGA); 

                disp('Finished BUILD SYSTEM IGA');                
                
                %-----------------------------------------------------------------%
                %                  	      PROBLEM SOLUTION                        %
                %-----------------------------------------------------------------%

                u = A\b;
                
                %------------------%
                % REBUILD SOLUTION %
                %------------------%
                
                buildBC = BoundaryConditionHandler();
                
                buildBC.bcStruct = boundStruct;
                buildBC.uRid = u;
                
                [uAug] = buildBoundCond(buildBC);
                u = uAug;
                
                disp('Finished SOLVE SYSTEM');

                %-----------------------------------------------------------------%
                %                  	       SOLUTION PLOT                          %
                %-----------------------------------------------------------------%

%                 [errL2,errH1] = plot_solution_IGA_scatter( ...
%                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
%                 obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map);
% 
%                 errorNormH1 = errL2;
%                 errorNormL2 = errH1;
                
                disp('Finished PLOT OPERATION / ERROR with EXACT SOLUTION')
            
                [errL2,errH1] = computeErrorIGA_scatter( ...
                obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map);
             
                disp('Finished ERROR with FREEFEM++ SOLUTION')

                errorNormH1 = errH1;
                errorNormL2 = errL2;

                disp('Finished Method SOLVER IGA')

            end % End function
            
            %% Method 'solverIGAScatterTransient'
            
            function [solutionMatrix,errL2,errH1] = solverIGAScatterTransient(obj)

            
            %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%

                % Initialization of the Data Structure that will further contain the 
                % Labels relative to the Nature of the Boundary Conditions at the 
                % Interface
                
                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%

                numbIterations = length(obj.timeDomain);
                numbStates = obj.numbControlPoints*obj.dimModalBasis;

                solutionMatrix = zeros(numbStates,numbIterations);
                
                solutionMatrix(:,1) = obj.initialState;
                
                for iteration = 1 : numbIterations - 1
                    
                    % Definition of the Object from the AssemblerADRHandler class

                    build_IGA = AssemblerADRHandler();

                    % Properties Assignment

                    build_IGA.dimModalBasis = obj.dimModalBasis(1);
                    build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                    build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                    build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                    build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                    build_IGA.stepMeshX = obj.stepMeshX(1);
                    build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                    build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                    build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                    build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                    build_IGA.coefficientForm = obj.coefficientForm;
                    build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                    build_IGA.geometricInfo = obj.geometricInfo;
                    build_IGA.robinCondStruct = obj.robinCondStruct;
                    build_IGA.couplingCond_DD = obj.couplingCond_DD;
                    build_IGA.physicMesh_inX = obj.physicMesh_inX;
                    build_IGA.physicMesh_inY = obj.physicMesh_inY;
                    build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                    build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                    build_IGA.continuityParameter = obj.continuityParameter;
                    build_IGA.domainProfile = obj.domainProfile;
                    build_IGA.domainProfileDer = obj.domainProfileDer;
                    build_IGA.numbHorQuadNodes = obj.numbHorQuadNodes;
                    build_IGA.numbVerQuadNodes = obj.numbVerQuadNodes;
                    build_IGA.timeInstant = iteration;
                    build_IGA.timeDomain = obj.timeDomain;


                    % Call of the 'buildSystemIGA' Method

                    [A,M,b,~,liftCoeffA(1),liftCoeffB(1),space,refDomain1D] = buildSystemIGAScatterTransient(build_IGA);            

                    % Computation of the State Matrix

                    stateMatrix = (M + obj.timeStep * A)\M;

                    % Computation of the Input Matrix

                    inputMatrix = M + obj.timeStep*A;
                    
                    %-----------------------------------------------------------------%
                    %                  	      PROBLEM SOLUTION                        %
                    %-----------------------------------------------------------------%
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\(obj.timeStep * b);

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                    display = ['Finished Time Iteration ',num2str(iteration)];
                    disp(display);
                    
                end

                %-----------------------------------------------------------------%
                %                  	       SOLUTION PLOT                          %
                %-----------------------------------------------------------------%
                
                numbIterations = length(obj.timeDomain) - 1;
                
                errorHistoryL2 = zeros(numbIterations,1);
                errorHistoryH1 = zeros(numbIterations,1);
                
                % Create the Freefem++ simulation folder
    
                for ii = 1:1000
        
                checkFolder = exist(['MatlabPlots',num2str(ii)]);
        
                    if (checkFolder == 7)
                        disp(['The folder MatlabPlots',num2str(ii),' already exists!'])
                    else
                        fileNameF = ['MatlabPlots',num2str(ii)];
                        mkdir(fileNameF)
                        break
                    end
                end
                
                % Print the solution
                
                for ii = 1:numbIterations

                    [errL2,errH1] = plot_solution_IGA_scatter_transient( ...
                                    obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                                    obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map,ii,obj.timeDomain);
            
                    errorHistoryL2(ii) = errL2;
                    errorHistoryH1(ii) = errH1;

                end

                disp('Finished Plot Operation');

                % Create simulation video

                imageNames = cell(1,numbIterations);

                for ii = 1:numbIterations
                    fileName = ['Plot_At_t=',num2str(ii),'.png'];
                    imageNames{ii} = fileName;
                end

                workingDir = [pwd,'/',fileNameF];

                outputVideo = VideoWriter(fullfile(workingDir,'TransientMatlabEvolution.avi'));
                outputVideo.FrameRate = 10;
                open(outputVideo)

                for ii = 1:length(imageNames)
                img = imread(fullfile(workingDir,imageNames{ii}));
                writeVideo(outputVideo,img)
                end

                close(outputVideo)

            end % End function
            
            %% Method 'solverIGAScatter3D'
            
            function [u,liftCoeffA,liftCoeffB,errorNormL2,errorNormH1] = solverIGAScatter3D(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%

                % Initialization of the Data Structure that will further contain the 
                % Labels relative to the Nature of the Boundary Conditions at the 
                % Interface

                BC_l = cell(2,1);
                BC_r = cell(2,1);
                
                BC_l{1}  = 'dir';   % Inflow Boundary Condition
                BC_r{2} = 'neu';   % Outflow Boundary Condition


                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%
                
                tic;
                
                % Definition of the Object from the AssemblerADRHandler class

                build_IGA = AssemblerADRHandler();

                % Properties Assignment

                build_IGA.dimModalBasis = obj.dimModalBasis(1);
                build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                build_IGA.stepMeshX = obj.stepMeshX(1);
                build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                build_IGA.coefficientForm = obj.coefficientForm;
                build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                build_IGA.geometricInfo = obj.geometricInfo;
                build_IGA.robinCondStruct = obj.robinCondStruct;
                build_IGA.couplingCond_DD = obj.couplingCond_DD;
                build_IGA.physicMesh_inX = obj.physicMesh_inX;
                build_IGA.physicMesh_inY = obj.physicMesh_inY;
                build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                build_IGA.continuityParameter = obj.continuityParameter;
                build_IGA.domainProfile = obj.domainProfile;
                build_IGA.domainProfileDer = obj.domainProfileDer;
                build_IGA.numbHorQuadNodes = obj.numbHorQuadNodes;
                build_IGA.numbVerQuadNodes = obj.numbVerQuadNodes;

                % Call of the 'buildSystemIGA' Method

                [A,b,modalBasisStruct,liftCoeffA(1),liftCoeffB(1),space,refDomain1D] = buildSystemIGAScatter3D(build_IGA); 

                disp('Finished BUILD SYSTEM IGA');                
                
                %-----------------------------------------------------------------%
                %                  	      PROBLEM SOLUTION                        %
                %-----------------------------------------------------------------%
                
                % Debug 
                %---------------------------------------------------------------------%
                % plot(linspace(0,length(b),length(b)),b);
                %---------------------------------------------------------------------%

                u = A\b;
                
                timeSolHigaMod = toc;
                
                % DEBUG
                
                disp(' ');
                disp(' ');
                disp(' ');
                disp('*************');
                disp('*** DEGUB ***');
                disp('*************');
                disp(' ');
                disp(['Size Sol. : ',num2str(size(u))]);
                disp(['Sol. Time : ',num2str(timeSolHigaMod)]);
                
                disp('Finished SOLVE SYSTEM');

                %-----------------------------------------------------------------%
                %                  	       SOLUTION PLOT                          %
                %-----------------------------------------------------------------%

                [errL2,errH1] = computeErrorIGA_scatter_3D( ...
                obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map,...
                obj.geometricInfo);

                errorNormH1 = errL2;
                errorNormL2 = errH1;
                
                disp('Finished PLOT OPERATION / ERROR with EXACT SOLUTION')

            end % End function
            
            %% Method 'solverIGAScatter3DTransient'
            
            function [solutionMatrix,liftCoeffA,liftCoeffB,errorHistoryL2,errorHistoryH1] = solverIGAScatter3DTransient(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%


                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%
                                
                numbIterations = length(obj.timeDomain);
                numbStates = obj.numbControlPoints*obj.dimModalBasis^2;
                
                solutionMatrix = zeros(numbStates,numbIterations);
                solutionMatrix(:,1) = zeros(numbStates,1);
                
                tic;
                
                for iteration = 1 : numbIterations - 1

                    % Definition of the Object from the AssemblerADRHandler class

                    build_IGA = AssemblerADRHandler();

                    % Properties Assignment

                    build_IGA.dimModalBasis = obj.dimModalBasis(1);
                    build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                    build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                    build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                    build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                    build_IGA.stepMeshX = obj.stepMeshX(1);
                    build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                    build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                    build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                    build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                    build_IGA.coefficientForm = obj.coefficientForm;
                    build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                    build_IGA.geometricInfo = obj.geometricInfo;
                    build_IGA.robinCondStruct = obj.robinCondStruct;
                    build_IGA.couplingCond_DD = obj.couplingCond_DD;
                    build_IGA.physicMesh_inX = obj.physicMesh_inX;
                    build_IGA.physicMesh_inY = obj.physicMesh_inY;
                    build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                    build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                    build_IGA.continuityParameter = obj.continuityParameter;
                    build_IGA.domainProfile = obj.domainProfile;
                    build_IGA.domainProfileDer = obj.domainProfileDer;
                    build_IGA.numbHorQuadNodes = obj.numbHorQuadNodes;
                    build_IGA.numbVerQuadNodes = obj.numbVerQuadNodes;
                    build_IGA.timeInstant = obj.timeDomain(iteration + 1);

                    % Call of the 'buildSystemIGA' Method

                    [A,M,b,modalBasisStruct,liftCoeffA(1),liftCoeffB(1),space,refDomain1D] = buildSystemIGAScatter3DTransient(build_IGA); 
                    
                    % Computation of the State Matrix

                    stateMatrix = (M + obj.timeStep * A)\M;

                    % Computation of the Input Matrix

                    inputMatrix = M + obj.timeStep*A;

                    % Save history of state and input matrices

                    stateMatHist{iteration} = stateMatrix;
                    inputMatHist{iteration} = inputMatrix;

                    % Save the force history

                    forceMatrix = obj.timeStep * inv(inputMatrix);
                    forceHistory(:,iteration) = b;
                    
                    % Compute state and source contribution
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\(obj.timeStep * b);
                    
                    % Compute and save the new solution

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                    disp(['Finished solving at iteration ',num2str(iteration)]);
                    
                end
                
                timeSolHigaMod = toc;
                
                % DEBUG
                
                disp(' ');
                disp(' ');
                disp(' ');
                disp('*************');
                disp('*** DEGUB ***');
                disp('*************');
                disp(' ');
                disp(['Size Sol. : ',num2str(size(solutionMatrix))]);
                disp(['Sol. Time : ',num2str(timeSolHigaMod)]);
                
                disp('Finished SOLVE SYSTEM');

                %-----------------------------------------------------------------%
                %                  	       SOLUTION PLOT                          %
                %-----------------------------------------------------------------%

                aux = [];
                numbIterations = length(obj.timeDomain) - 1;
                
                errorHistoryL2 = zeros(numbIterations,1);
                errorHistoryH1 = zeros(numbIterations,1);
                
                % Create the Freefem++ simulation folder
    
                for ii = 1:1000
        
                checkFolder = exist(['MatlabPlots3D_SimCase',num2str(ii)]);
        
                    if (checkFolder == 7)
                        disp(['The folder MatlabPlots3D_SimCase',num2str(ii),' already exists!'])
                    else
                        fileNameF = ['MatlabPlots3D_SimCase',num2str(ii)];
                        mkdir(fileNameF)
                        break
                    end
                end
                
                % Print the solution
                
                for ii = 1:numbIterations

                [errL2,errH1,MeshStruct] = computeErrorIGA_scatter_3D_Transient( ...
                obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map,...
                obj.geometricInfo,ii,obj.timeDomain,aux);
            
                aux = MeshStruct;

                errorHistoryL2(ii) = errL2;
                errorHistoryH1(ii) = errH1;

                disp(['FINISHED CREATING PLOT AT ITERATION ',num2str(ii)]);

                end
                
                disp('FINISHED PLOT OPERATION')
                
%                 %% Create simulation video
%                 
%                 imageNames = cell(numbIterations,1);
%                 for ii = 1:numbIterations
%                     fileName = ['Plot_At_it=',num2str(ii),'.jpg'];
%                     imageNames{ii} = fileName;
%                 end
% 
%                 workingDir = fileNameF;
% 
%                 outputVideo = VideoWriter(fullfile(workingDir,'HigaMod3D.avi'));
%                 outputVideo.FrameRate = 10;
%                 open(outputVideo)
% 
%                 for ii = 1:length(imageNames)
%                     img = imread(fullfile(workingDir,imageNames{ii}));
%                     writeVideo(outputVideo,img)
%                 end
% 
%                 close(outputVideo);     

            end % End function
            
            %% Method 'solverIGAScatterStabilized'
            
            function [u,liftCoeffA,liftCoeffB,errorNormL2,errorNormH1,CutCoord,refSolStructComp] = solverIGAScatterStabilized(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%

                % Initialization of the Data Structure that will further contain the 
                % Labels relative to the Nature of the Boundary Conditions at the 
                % Interface

                BC_l = cell(2,1);
                BC_r = cell(2,1);
                
                BC_l{1}  = 'dir';   % Inflow Boundary Condition
                BC_r{2} = 'neu';   % Outflow Boundary Condition


                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%

                % Definition of the Object from the AssemblerADRHandler class

                build_IGA = AssemblerADRHandler();

                % Properties Assignment

                build_IGA.dimModalBasis = obj.dimModalBasis(1);
                build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                build_IGA.stepMeshX = obj.stepMeshX(1);
                build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                build_IGA.coefficientForm = obj.coefficientForm;
                build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                build_IGA.geometricInfo = obj.geometricInfo;
                build_IGA.robinCondStruct = obj.robinCondStruct;
                build_IGA.couplingCond_DD = obj.couplingCond_DD;
                build_IGA.physicMesh_inX = obj.physicMesh_inX;
                build_IGA.physicMesh_inY = obj.physicMesh_inY;
                build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                build_IGA.continuityParameter = obj.continuityParameter;
                build_IGA.domainProfile = obj.domainProfile;
                build_IGA.domainProfileDer = obj.domainProfileDer;
                build_IGA.numbHorQuadNodes = obj.numbHorQuadNodes;
                build_IGA.numbVerQuadNodes = obj.numbVerQuadNodes;
                build_IGA.stabMethod = obj.stabMethod;
                build_IGA.stabDelta = obj.stabDelta;
                build_IGA.PeConfig = obj.PeConfig;

                % Call of the 'buildSystemIGA' Method

                [A,b,~,liftCoeffA(1),liftCoeffB(1),space,refDomain1D] = buildSystemIGAScatterStabilized(build_IGA); 

                disp('Finished BUILD SYSTEM IGA');                
                
                %-----------------------------------------------------------------%
                %                  	      PROBLEM SOLUTION                        %
                %-----------------------------------------------------------------%
                
                % Debug 
                %---------------------------------------------------------------------%
                % plot(linspace(0,length(b),length(b)),b);
                %---------------------------------------------------------------------%

                u = A\b;
                
                disp('Finished SOLVE SYSTEM');

                %-----------------------------------------------------------------%
                %                  	       SOLUTION PLOT                          %
                %-----------------------------------------------------------------%
                
                simulationInfo.H = obj.stepMeshX(1);
                simulationInfo.M = obj.dimModalBasis;
                simulationInfo.Pe = obj.PeConfig;
                simulationInfo.Delta = obj.stabDelta;
                simulationInfo.Method = obj.stabMethod;
                
                if (obj.selectedPlots(end) == 0)
                    
                    [errL2,errH1,CutCoord] = plot_solution_IGA_scatter_stab( ...
                    obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                    obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map,...
                    simulationInfo);

                    refSolStructComp = obj.refSolStruct;
                    errorNormH1 = errL2;
                    errorNormL2 = errH1;
                
                    disp('Finished PLOT OPERATION / ERROR with EXACT SOLUTION')
                
                else
            
                    [errL2,errH1,CutCoord,refSolStructComp] = computeErrorIGA_scatter_stab( ...
                    obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                    obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map, ...
                    simulationInfo,obj.selectedPlots,obj.refSolStruct);

                    disp('Finished ERROR with FREEFEM++ SOLUTION')

                    errorNormH1 = errH1;
                    errorNormL2 = errL2;

                end

                disp('Finished Method SOLVER IGA');

            end % End function
            
            %% Method 'solverIGANoisy'
            
            function [solutionMatrix,Z,X_POST,liftCoeffA,liftCoeffB,...
                      errorNormL2,errorNormH1,errorNoiseNormL2,errorNoiseNormH1,...
                      errorEstimateNormL2,errorEstimateNormH1,P_POST] = solverIGANoisy(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%

                % Initialization of the Data Structure that will further contain the 
                % Labels relative to the Nature of the Boundary Conditions at the 
                % Interface
                
                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%
                
                %% SOLUTION EVOLUTION

                numbIterations = length(obj.timeDomain);
                numbStates = (round(1/obj.stepMeshX(1)) + 1) * obj.dimModalBasis(1);
                
                solutionMatrix = zeros(numbStates,numbIterations);
                solutionMatrix(:,1) = obj.initialState;
                
                forceHistory = zeros(numbStates,numbIterations - 1);
                
                tic;
                
                for iteration = 1 : numbIterations - 1
                    
                    if (iteration == 1)
                
                        % Definition of the Object from the AssemblerADRHandler class

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = obj.coefficientForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildSystemIGATransient(build_IGA);             
                        
                        % Computation of the State Matrix
                        
                        Mdt = M * (obj.timeStep^-1);

                        stateMatrix = (Mdt + A)\Mdt;

                        % Computation of the Input Matrix

                        inputMatrix = (Mdt + A);
                        
                        forceHistory(:,1) = b;
                        
                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);    
                    
                    else
                        
                        % Definition of the Object from the AssemblerADRHandler class

                        build_force = AssemblerADRHandler();

                        % Properties Assignment

                        build_force.dimModalBasis = obj.dimModalBasis(1);
                        build_force.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_force.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_force.stepMeshX = obj.stepMeshX(1);
                        build_force.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_force.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_force.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_force.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_force.coefficientForm = obj.coefficientForm;
                        build_force.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_force.geometricInfo = obj.geometricInfo;
                        build_force.robinCondStruct = obj.robinCondStruct;
                        build_force.couplingCond_DD = obj.couplingCond_DD;
                        build_force.physicMesh_inX = obj.physicMesh_inX;
                        build_force.physicMesh_inY = obj.physicMesh_inY;
                        build_force.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_force.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_force.continuityParameter = obj.continuityParameter;
                        build_force.domainProfile = obj.domainProfile;
                        build_force.domainProfileDer = obj.domainProfileDer;
                        build_force.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildIGAForce(build_force); 
                        
                        forceHistory(:,iteration) = b;

                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);              

                    end
                    
                    %-----------------------------------------------------------------%
                    %                  	      PROBLEM SOLUTION                        %
                    %-----------------------------------------------------------------%
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\b;

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                end
                                
                solveTime = toc;
                
                display = sprintf('Finished Solving the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,solveTime);
                disp(display);
                
                %% MEASUREMENT PARAMETERS

                % Selection matrix
                                
                M = obj.numbVerQuadNodes;
                Nh = obj.numbControlPoints;
                
                % Hs = obj.measureMatrix;
                
                % SELECT EVERY POINT IN THE DOMAIN
                 Hs = eye(obj.numbControlPoints * obj.numbVerQuadNodes);
                
                % SELECT A LINE IN THE DOMAIN
                % Hs = zeros(Nh,M * Nh);
                
                % for ii = 1:Nh
                %     Hs(ii, ceil(M/2) * Nh + ii) = 1;
                %  end
                
                % IGA basis and modal basis matrices
                
                [Hb,Hm] = observStates(obj.dimModalBasis,obj.domainLimit_inX,obj.domainLimit_inY,...
                                       obj.stepMeshX,obj.label_upBoundDomain,obj.label_downBoundDomain,...
                                       obj.coefficientForm,obj.degreePolySplineBasis,obj.continuityParameter);
                                   
                % Observation matrix
                
                H = Hs * Hm * Hb;
                
                % Covariance matrices
                
                [numbStates,numbSamples] = size(solutionMatrix);
                
%                 Q1 = 10^-8.5;
%                 R1 = 10^-7;
                
                Q1 = 10^-2.5;
                R1 = 10^-1;
                
                Q = Q1 * eye(numbStates);
                R = R1 * eye(size(Hs,1));
                         
                %% FILTER ALGORITHM
                
                %-----------------------------------------------------------------%
                %           	       KALMAN FILTER ALGORITHM                    %
                %-----------------------------------------------------------------%
                
                tic;
                
                % Matrix Initialization

                Z = zeros(size(Hs,1),numbSamples);
                X_PRIOR = zeros(numbStates,numbSamples);
                X_POST = zeros(numbStates,numbSamples);
                RESIDUAL = zeros(size(Hs,1),numbSamples);
                P_PRIOR = cell(1,numbSamples);
                P_POST = cell(1,numbSamples);
                GAIN = cell(1,numbSamples);
                
                % Create measurements
                
                for jj = 1:numbSamples
                    measureStd = 0.5 * max(H * solutionMatrix(:,jj));
                    v = measureStd * randn(size(Hs,1),1);
                    Z(:,jj) = H * solutionMatrix(:,jj) + v;
                    
                end
                
                % Initial State

                X_PRIOR(:,1) = obj.initialState;
                X_POST(:,1) = obj.initialState;
                
                P_POST{1} = Q;

                for jj = 1:numbSamples-1
                    
                    % Prediction equations
                        
                    X_PRIOR(:,jj + 1) = stateMatrix * X_POST(:,jj);
                
                    P_PRIOR{jj + 1} = stateMatrix * (P_POST{jj} * stateMatrix') + Q;
                    
                    % Assimilation equations
                    
                    RESIDUAL(:,jj + 1) = Z(:,jj + 1) - H * X_PRIOR(:,jj + 1);
                    
                    IS = (R + H * P_PRIOR{jj + 1} * H');
                    GAIN{jj} = (P_PRIOR{jj + 1} * H')*inv(IS);
                
                    P_POST{jj+1} = (eye(size(GAIN{jj} * H)) - GAIN{jj} * H) * P_PRIOR{jj + 1};
                
                    X_POST(:,jj + 1) = X_PRIOR(:,jj + 1) + GAIN{jj} * RESIDUAL(:,jj + 1);

                   
                end

                filterTime = toc;
                
                display = sprintf('Finished Filtering the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,filterTime);
                disp(display);
                
                %% DATA EXPORT
                
                %-----------------------------------------------------------------%
                %                  	       DATA EXPORT                            %
                %-----------------------------------------------------------------%
                
                tic;
                
                %% CREATE EXPORT FOLDER
    
                for ii = 1:1000
        
                checkFolder = exist(['MatlabPlots',num2str(ii)]);
        
                    if (checkFolder == 7)
                        disp(['The folder MatlabPlots',num2str(ii),' already exists!'])
                    else
                        fileNameF = ['MatlabPlots',num2str(ii)];
                        mkdir(fileNameF)
                        break
                    end
                end
                
                %% PRINT SYSTEM EVOLUTION
                
                 for ii = 1:numbIterations
 
                 [errL2,errH1,~,~] = plot_KF( ...
                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.domainLimit_inY,obj.stepMeshX, ...
                 solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                 obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer, ...
                 obj.jacAtQuadNodes, obj.degreePolySplineBasis,obj.exactSolution_dX,obj.exactSolution_dY,obj.continuityParameter, ...
                 ii,obj.timeDomain,Z(:,ii),X_POST(:,ii),Hs,Hm,Hb);
             
                 errorHistoryL2(ii) = errL2;
                 errorHistoryH1(ii) = errH1;
             
                 end
                 
                 disp('Finished Plot Operation');
                 
                %% CREATE NOISY MEASUREMENT VIDEO
                
                for ii = 1:numbIterations
                    fileName = ['Plot_At_it=',num2str(ii),'.png'];
                    imageNames{ii} = fileName;
                end
                 
                workingDir = [pwd,'\',fileNameF];
                 
                outputVideo = VideoWriter(fullfile(workingDir,'Assimilated State.avi'));
                outputVideo.FrameRate = 10;
                open(outputVideo)

                for ii = 1:length(imageNames)
                    img = imread(fullfile(workingDir,imageNames{ii}));
                    writeVideo(outputVideo,img)
                end

                close(outputVideo)
                
                %%%%%%%%%
                % DEBUG %
                %%%%%%%%%
                
                errorNormH1 = 0;
                errorNormL2 = 0;
                errorNoiseNormH1 = 0;
                errorNoiseNormL2 = 0;
                errorEstimateNormH1 = 0;
                errorEstimateNormL2 = 0;

                disp('Finished Computing Error for the IGA Method');

            end % End function
            
            %% Method 'solverIGAForceID'
            
            function [solutionMatrix,Z,X_POST] = solverIGAForceID(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %---------------------------------------------------------------------%
                %                        PROBLEM INITIALIZATION                       %
                %---------------------------------------------------------------------%

                % Initialization of the Data Structure that will further contain the 
                % Labels relative to the Nature of the Boundary Conditions at the 
                % Interface
                
                %---------------------------------------------------------------------%
                % Note:
                % In the case of Physical Boundary, the Labels of the Boundary 
                % Conditions are known. At the inflow we have Dirichlet Boundary 
                % Conditions and at the outflow, Neumann (Homogeneous) Boundary 
                % Conditions.
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % The assembling of the block system is performed by the external
                % function 'buildSystem'. The assembling can be set differentling
                % depending on the discretization strategy used to solve the problem.
                % In the present function we use an Isogeometric (IGA) Approach to
                % solve the problem.
                %---------------------------------------------------------------------%
                
                %% Evolution of the solution

                numbIterations = length(obj.timeDomain);
                numbStates = (round(1/obj.stepMeshX(1)) + 1) * obj.dimModalBasis(1);
                
                solutionMatrix = zeros(numbStates,numbIterations);
                solutionMatrix(:,1) = obj.initialState;
                
                forceHistory = zeros(numbStates,numbIterations - 1);
                
                tic;
                
                for iteration = 1 : numbIterations - 1
                    
                    if (iteration == 1)
                
                        % Definition of the Object from the AssemblerADRHandler class

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = obj.coefficientForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildSystemIGATransient(build_IGA);             
                        
                        % Computation of the State Matrix
                        
                        Mdt = M * (obj.timeStep^-1);

                        stateMatrix = (Mdt + A)\Mdt;

                        % Computation of the Input Matrix

                        inputMatrix = (Mdt + A);
                        
                        forceHistory(:,1) = b;
                        
                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);    
                    
                    else
                        
                        % Definition of the Object from the AssemblerADRHandler class

                        build_force = AssemblerADRHandler();

                        % Properties Assignment

                        build_force.dimModalBasis = obj.dimModalBasis(1);
                        build_force.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_force.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_force.stepMeshX = obj.stepMeshX(1);
                        build_force.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_force.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_force.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_force.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_force.coefficientForm = obj.coefficientForm;
                        build_force.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_force.geometricInfo = obj.geometricInfo;
                        build_force.robinCondStruct = obj.robinCondStruct;
                        build_force.couplingCond_DD = obj.couplingCond_DD;
                        build_force.physicMesh_inX = obj.physicMesh_inX;
                        build_force.physicMesh_inY = obj.physicMesh_inY;
                        build_force.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_force.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_force.continuityParameter = obj.continuityParameter;
                        build_force.domainProfile = obj.domainProfile;
                        build_force.domainProfileDer = obj.domainProfileDer;
                        build_force.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildIGAForce(build_force); 
                        
                        forceHistory(:,iteration) = b;

                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);              

                    end
                    
                    %-----------------------------------------------------------------%
                    %                  	      PROBLEM SOLUTION                        %
                    %-----------------------------------------------------------------%
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\b;

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                end
                                
                solveTime = toc;
                
                display = sprintf('Finished Solving the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,solveTime);
                disp(display);
                
                %% Compute force components

                dato_dir = @(y) 0;
                force = @(x,y,t) 1;
                forceData = struct('dato_dir',dato_dir,'force',force);
                
                for iteration = 1 : numbIterations - 1

                    build_force = AssemblerADRHandler();

                    % Properties Assignment

                    build_force.dimModalBasis = obj.dimModalBasis(1);
                    build_force.leftBDomain_inX = obj.domainLimit_inX(1);
                    build_force.rightBDomain_inX = obj.domainLimit_inX(2);
                    build_force.stepMeshX = obj.stepMeshX(1);
                    build_force.label_upBoundDomain = obj.label_upBoundDomain{1};
                    build_force.label_downBoundDomain = obj.label_downBoundDomain{1};
                    build_force.localdata_upBDomain = obj.data_upBoundDomain{1};
                    build_force.localdata_downBDomain = obj.data_downBoundDomain{1};
                    build_force.coefficientForm = obj.coefficientForm;
                    build_force.dirCondFuncStruct = forceData;
                    build_force.geometricInfo = obj.geometricInfo;
                    build_force.robinCondStruct = obj.robinCondStruct;
                    build_force.couplingCond_DD = obj.couplingCond_DD;
                    build_force.physicMesh_inX = obj.physicMesh_inX;
                    build_force.physicMesh_inY = obj.physicMesh_inY;
                    build_force.jacAtQuadNodes = obj.jacAtQuadNodes;
                    build_force.degreePolySplineBasis = obj.degreePolySplineBasis;
                    build_force.continuityParameter = obj.continuityParameter;
                    build_force.domainProfile = obj.domainProfile;
                    build_force.domainProfileDer = obj.domainProfileDer;
                    build_force.timeInstant = iteration-1;

                    % Call of the 'buildSystemIGA' Method

                    [ref0RHS,~,~,~,~,~] = buildIGAForce(build_force); 

                    refHistory(:,iteration) = ref0RHS;
                    
                end
                
                %% Creation of noise matrices 

                [numbStates,numbSamples] = size(solutionMatrix);
                
                Q1 = 10^-8.5;
                R1 = 10^-7;
                
                Q = Q1 * eye(numbStates);
                R = R1 * eye(numbStates);
                
                %% Augmented system
                
                sysMatAug = [stateMatrix inv(inputMatrix); zeros(size(stateMatrix)) eye(size(inputMatrix))];
                meaMatAug = [obj.measureMatrix zeros(size(inputMatrix))];
                
                %-----------------------------------------------------------------%
                %           	       KALMAN FILTER ALGORITHM                    %
                %-----------------------------------------------------------------%
                
                tic;
                
                % Matrix Initialization

                Z = zeros(numbStates,numbSamples);              % Output Evolution
                X_PRIOR = zeros(2 * numbStates,numbSamples);
                X_POST = zeros(2 * numbStates,numbSamples);
                RESIDUAL = zeros(numbStates,numbSamples);
                P_PRIOR = cell(1,numbSamples);
                P_POST = cell(1,numbSamples);
                GAIN = cell(1,numbSamples);
                
                %% Create measurements
                
                for jj = 1:numbSamples
                    measureStd = 0.05 * max(solutionMatrix(:,jj));
                    v = measureStd * randn(size(solutionMatrix(:,jj)));
                    Z(:,jj) = solutionMatrix(:,jj) + v;
                    
                end
                
                %% Filtering algorithm
                
                % Initial State
                
                X_0_Aug = [ obj.initialState * 10; zeros(size(b)) ];

                X_PRIOR(:,1) = X_0_Aug;
                X_POST(:,1) = X_0_Aug;
                
                QAug = [ Q zeros(size(Q)); zeros(size(Q)) Q ];
                
                P_POST{1} = [ Q zeros(size(Q)); zeros(size(Q)) Q*1e3 ];

                
                for jj = 1:numbSamples-1
                    
                    %Predictor equations
                        
                    X_PRIOR(:,jj + 1) = sysMatAug * X_POST(:,jj);
                
                    P_PRIOR{jj + 1} = sysMatAug * (P_POST{jj} * sysMatAug') + QAug;
                    
                    % Corrector equations
                    
                    RESIDUAL(:,jj + 1) = Z(:,jj + 1) - meaMatAug * X_PRIOR(:,jj + 1);
                    
                    IS = (R + meaMatAug * P_PRIOR{jj + 1} * meaMatAug');
                    GAIN{jj} = (P_PRIOR{jj + 1} * meaMatAug')*inv(IS);
                
                    P_POST{jj+1} = (eye(2 * numbStates) - GAIN{jj} * meaMatAug) * P_PRIOR{jj + 1};
                
                    X_POST(:,jj + 1) = X_PRIOR(:,jj + 1) + GAIN{jj} * RESIDUAL(:,jj + 1);

                   
                end

                filterTime = toc;
                
                display = sprintf('Finished Filtering the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,filterTime);
                disp(display);
                
                %% Estimate the force coefficient
                
                coeffMean = zeros(numbSamples - 1,1);
                coeffVar  = zeros(numbSamples - 1,1);
                
                for jj = 1:numbSamples-1

                
                    estimatedCoeff = X_POST(numbStates + 1:numbStates + 11,jj + 1) ./ refHistory(1:11,jj);   

                    coeffMean(jj) = mean(estimatedCoeff(2:end));
                    coeffVar(jj)  = var(estimatedCoeff(2:end));
                    pause(1)
                   
                end

%                 disp(coeffMean);
%                 disp(coeffVar)
%                 
%                 figure;
%                 plot(linspace(0,1,length(coeffMean(2:end))),coeffMean(2:end));
%                 figure;
%                 plot(linspace(0,1,length(coeffVar(2:end))),coeffVar(2:end));


                coeffMeanAug = zeros(size(coeffMean,1) + 1,1);
                coeffMeanAug(2:end) = coeffMean;
                coeffVarAug = zeros(size(coeffVar,1) + 1,1);
                coeffVarAug(2:end) = coeffVar;

                figure;
                
                minX = 1;
                maxX = size(coeffMean,1) + 1;
                minY = 0;
                maxY = 2;
                axis([minX maxX minY maxY]);
                axis manual
                plot(2 * ones(length(coeffMeanAug)),'--k','LineWidth',2); hold on
                errorbar(1:1:length(coeffMeanAug),coeffMeanAug,coeffVarAug,'--o','MarkerSize',6,'MarkerEdgeColor','red','MarkerFaceColor','white','LineWidth',2)
                set(gca, 'FontSize', 10);
                hold off
                
%                 tic;
%                 
%                 % Matrix Initialization
%                 
%                 X = zeros(numbStates,numbSamples);              % State Evolution
%                 Z = zeros(numbStates,numbSamples);              % Output Evolution
%                 X_PRIOR = zeros(numbStates,numbSamples);
%                 X_POST = zeros(numbStates,numbSamples);
%                 RESIDUAL = zeros(numbStates,numbSamples);
%                 P_PRIOR = cell(1,numbSamples);
%                 P_POST = cell(1,numbSamples);
%                 GAIN = cell(1,numbSamples);
%                 
%                 % Create measurements
%                 
%                 for jj = 1:numbSamples
%                     measureStd = 0.05 * max(solutionMatrix(:,jj));
%                     v = measureStd * randn(size(solutionMatrix(:,jj)));
%                     Z(:,jj) = solutionMatrix(:,jj) + v;
%                     
%                     disp(measureStd^2)
%                     
%                 end
%                 
%                 % Initial State
%                 
%                 X_0 = obj.initialState;
% 
%                 X_PRIOR(:,1) = X_0;
%                 X_POST(:,1) = X_0;
%                 P_POST{1} = Q * 1e-2;
% 
%                 
%                 for jj = 1:numbSamples-1
%                     
%                     %Predictor equations
%                     
%                     stateInovation = stateMatrix * X_POST(:,jj);
%                     inputContrib = inputMatrix\(obj.timeStep * forceHistory(:,jj));
%                         
%                     X_PRIOR(:,jj + 1) = stateInovation + inputContrib;
%                 
%                     P_PRIOR{jj + 1} = stateMatrix * (P_POST{jj} * stateMatrix') + Q;
%                     
%                     % Corrector equations
%                     
%                     RESIDUAL(:,jj + 1) = Z(:,jj + 1) - (obj.measureMatrix) * X_PRIOR(:,jj + 1);
%                     
%                     IS = (R + (obj.measureMatrix) * P_PRIOR{jj + 1} * (obj.measureMatrix)');
%                     GAIN{jj} = (P_PRIOR{jj + 1} * obj.measureMatrix')*inv(IS);
%                 
%                     P_POST{jj+1} = (eye(numbStates) - GAIN{jj} * obj.measureMatrix) * P_PRIOR{jj + 1};
%                 
%                     X_POST(:,jj + 1) = X_PRIOR(:,jj + 1) + GAIN{jj} * RESIDUAL(:,jj + 1);
% 
%                    
%                 end
% 
%                 filterTime = toc;
%                 
%                 display = sprintf('Finished Filtering the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,filterTime);
%                 disp(display);
                
                %-----------------------------------------------------------------%
                %                  	       DATA EXPORT                            %
                %-----------------------------------------------------------------%
                
                tic;
                
%                 for ii = 1:1000
%         
%                 checkFolder = exist(['Paraview_',num2str(ii)]);
%         
%                     if (checkFolder == 7)
%                         disp(['The folder Paraview_',num2str(ii),' already exists!'])
%                     else
%                         directoryPath = ['Paraview_',num2str(ii)];
%                         mkdir(directoryPath)
%                         break
%                     end
%                 end
                
%                 oldFolder = cd(directoryPath);
%                 
%                 dataFileApprox = cell(numbSamples-1,1);
%                 dataFileNoisy = cell(numbSamples-1,1);
%                 dataFileEstimate = cell(numbSamples-1,1);
% 
%                 errorEvolH1 = zeros(numbSamples-1,1);
%                 errorEvolL2 = zeros(numbSamples-1,1);
%                 errorNoisyH1 = zeros(numbSamples-1,1);
%                 errorNoisyL2 = zeros(numbSamples-1,1);
%                 errorEstimatedH1 = zeros(numbSamples-1,1);
%                 errorEstimatedL2 = zeros(numbSamples-1,1);
                
%                 %% Create the Freefem++ simulation folder
%     
%                 for ii = 1:1000
%         
%                 checkFolder = exist(['MatlabPlots',num2str(ii)]);
%         
%                     if (checkFolder == 7)
%                         disp(['The folder MatlabPlots',num2str(ii),' already exists!'])
%                     else
%                         fileNameF = ['MatlabPlots',num2str(ii)];
%                         mkdir(fileNameF)
%                         break
%                     end
%                 end
%                 
%                 %% Print the solution
%                 
%                  for ii = 1:numbIterations
%  
%                  [errL2,errH1,~,~,~] = plot_KF( ...
%                  obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.domainLimit_inY,obj.stepMeshX, ...
%                  solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
%                  obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer, ...
%                  obj.jacAtQuadNodes, obj.degreePolySplineBasis,obj.exactSolution_dX,obj.exactSolution_dY,obj.continuityParameter, ...
%                  ii,obj.timeDomain,Z(1 : numbStates,ii),X_POST(1 : numbStates,ii),coeffMean,coeffVar);
%              
%                  errorHistoryL2(ii) = errL2;
%                  errorHistoryH1(ii) = errH1;
%              
%                  end
%                  
%                  disp('Finished Plot Operation');
%                  
%                 %% Create noisy measurements video
%                 
%                 for ii = 1:numbIterations
%                     fileName = ['Plot_At_it=',num2str(ii),'.png'];
%                     imageNames{ii} = fileName;
%                 end
%                  
%                 workingDir = ['/Users/YvesBarbosa/Documents/HigaModMat/ModifiedVersion/Demos/',fileNameF];
%                  
%                 outputVideo = VideoWriter(fullfile(workingDir,'Assimilated State.avi'));
%                 outputVideo.FrameRate = 10;
%                 open(outputVideo)
% 
%                 for ii = 1:length(imageNames)
%                     img = imread(fullfile(workingDir,imageNames{ii}));
%                     writeVideo(outputVideo,img)
%                 end
% 
%                 close(outputVideo)
% %                 
% % %                 jj = 10;
% % %                 
% % %                     [~,~,~,~,~,~] = plot_solution_transient_IGA( ...
% % %                     obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, solutionMatrix(:,jj+1),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
% % %                     obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer,obj.jacAtQuadNodes, ...
% % %                     obj.jacAtQuadNodes2,obj.degreePolySplineBasis,obj.continuityParameter,obj.exactSolution_dX,obj.exactSolution_dY,obj.timeDomain,jj);
% % %                 
% % %                 
% % %                 
% % %                     [~,~,~,~,~,~] = plot_solution_transient_IGA( ...
% % %                     obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, Z(:,jj + 1),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
% % %                     obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer,obj.jacAtQuadNodes, ...
% % %                     obj.jacAtQuadNodes2,obj.degreePolySplineBasis,obj.continuityParameter,obj.exactSolution_dX,obj.exactSolution_dY,obj.timeDomain,jj);
% % %                 
% % %                 
% % %                 
% % %                     [~,~,~,~,~,~] = plot_solution_transient_IGA( ...
% % %                     obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, X_POST(:,jj + 1),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
% % %                     obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer,obj.jacAtQuadNodes, ...
% % %                     obj.jacAtQuadNodes2,obj.degreePolySplineBasis,obj.continuityParameter,obj.exactSolution_dX,obj.exactSolution_dY,obj.timeDomain,jj);
% % %                 
% %              
% % %                 for ii = 1:numbSamples-1
% % % 
% % %                 [dataFileApprox{ii},errorEvolH1(ii),errorEvolL2(ii)] = export_solution_IGA( ...
% % %                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX, ...
% % %                 obj.stepMeshX, solutionMatrix(:,ii + 1),obj.label_upBoundDomain, ...
% % %                 obj.label_downBoundDomain,obj.coefficientForm,obj.simulationCase, ...
% % %                 obj.domainProfile,obj.domainProfileDer,obj.degreePolySplineBasis, ...
% % %                 obj.continuityParameter,ii,obj.jacAtQuadNodes,obj.exactSolution,  ...
% % %                 obj.exactSolution_dX,obj.exactSolution_dY,obj.jacAtQuadNodes2, obj.timeDomain);
% % %             
% % %                 [dataFileNoisy{ii},errorNoisyH1(ii),errorNoisyL2(ii)] = exportNoisyMeasure( ...
% % %                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX, ...
% % %                 obj.stepMeshX, Z(:,ii + 1),obj.label_upBoundDomain, ...
% % %                 obj.label_downBoundDomain,obj.coefficientForm,obj.simulationCase, ...
% % %                 obj.domainProfile,obj.domainProfileDer,obj.degreePolySplineBasis, ...
% % %                 obj.continuityParameter,ii,obj.jacAtQuadNodes,obj.exactSolution,  ...
% % %                 obj.exactSolution_dX,obj.exactSolution_dY,obj.jacAtQuadNodes2, obj.timeDomain);
% % %             
% % %                 [dataFileEstimate{ii},errorEstimatedH1(ii),errorEstimatedL2(ii)] = exportStateEstimate( ...
% % %                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX, ...
% % %                 obj.stepMeshX, X_POST(:,ii + 1),obj.label_upBoundDomain, ...
% % %                 obj.label_downBoundDomain,obj.coefficientForm,obj.simulationCase, ...
% % %                 obj.domainProfile,obj.domainProfileDer,obj.degreePolySplineBasis, ...
% % %                 obj.continuityParameter,ii,obj.jacAtQuadNodes,obj.exactSolution,  ...
% % %                 obj.exactSolution_dX,obj.exactSolution_dY,obj.jacAtQuadNodes2, obj.timeDomain,noisyMeasure(:,ii + 1));
% % %             
% % %                 end
% %                 
% % %                 errorHistoryL2 = errorEvolL2;
% % %                 errorHistoryH1 = errorEvolH1;
% % %                 
% % %                 errorNoiseHistoryL2 = errorNoisyL2;
% % %                 errorNoiseHistoryH1 = errorNoisyH1;
% % %                 
% % %                 errorEstimateHistoryL2 = errorEstimatedL2;
% % %                 errorEstimateHistoryH1 = errorEstimatedH1;
% % %                 
% % %                 disp('Finished Computing Paraview VTU Files');
% % %                 
% % %                 vtkseries(obj.timeDomain(2:end), dataFileApprox, 'ApproxEvolution.pvd');
% % %                 
% % %                 disp('Finished Approximated Surf VTK File');
% % %                 
% % %                 vtkseries(obj.timeDomain(2:end), dataFileNoisy, 'NoiseEvolution.pvd');
% % %                 
% % %                 disp('Finished Noisy Measurements Surf VTK File');
% % %                 
% % %                 vtkseries(obj.timeDomain(2:end), dataFileEstimate, 'StateEstimate.pvd');
% % %                 
% % %                 disp('Finished Noisy Measurements Surf VTK File');
% % %                 
% % %                 cd(oldFolder);
% % %                 
% % %                 exportTime = toc;
% % %                 
% % %                 display = sprintf('Finished Export the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,exportTime);
% % %                 disp(display);
% % %                 
%                 
%                 %% PREPARATION OF THE DATA FOR THE FUNCTION THAT COMPUTES THE ERROR
% 
%                 errorNormH1 = errorHistoryH1;
%                 errorNormL2 = errorHistoryL2;
%                 errorNoiseNormH1 = errorNoiseHistoryH1;
%                 errorNoiseNormL2 = errorNoiseHistoryL2;
%                 errorEstimateNormH1 = errorEstimateHistoryH1;
%                 errorEstimateNormL2 = errorEstimateHistoryL2;

                errorNormH1 = 0;
                errorNormL2 = 0;
                errorNoiseNormH1 = 0;
                errorNoiseNormL2 = 0;
                errorEstimateNormH1 = 0;
                errorEstimateNormL2 = 0;

                disp('Finished Computing Error for the IGA Method');

            end % End function
            
            %% Method 'solverEnKF'
            
            function [solutionMatrix,Z,assStateHistory, ...
                      condGain,filterTime] = solverEnKF(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                %% IMPORT CLASSES
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% SOLUTION EVOLUTION

                numbIterations = length(obj.timeDomain);
                numbStates = (round(1/obj.stepMeshX(1)) + 1) * obj.dimModalBasis(1);
                
                solutionMatrix = zeros(numbStates,numbIterations);
                solutionMatrix(:,1) = obj.initialState;
                
                forceHistory = zeros(numbStates,numbIterations - 1);
                
                tic;
                
                for iteration = 1 : numbIterations - 1
                    
                    if (iteration == 1)
                
                        % Definition of the Object from the AssemblerADRHandler class

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = obj.coefficientForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildSystemIGATransient(build_IGA);             
                        
                        % Computation of the State Matrix
                        
                        Mdt = M * (obj.timeStep^-1);

                        stateMatrix = (Mdt + A)\Mdt;

                        % Computation of the Input Matrix

                        inputMatrix = (Mdt + A);
                        
                        forceHistory(:,1) = b;
                        
                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);    
                    
                    else
                        
                        % Definition of the Object from the AssemblerADRHandler class

                        build_force = AssemblerADRHandler();

                        % Properties Assignment

                        build_force.dimModalBasis = obj.dimModalBasis(1);
                        build_force.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_force.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_force.stepMeshX = obj.stepMeshX(1);
                        build_force.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_force.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_force.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_force.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_force.coefficientForm = obj.coefficientForm;
                        build_force.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_force.geometricInfo = obj.geometricInfo;
                        build_force.robinCondStruct = obj.robinCondStruct;
                        build_force.couplingCond_DD = obj.couplingCond_DD;
                        build_force.physicMesh_inX = obj.physicMesh_inX;
                        build_force.physicMesh_inY = obj.physicMesh_inY;
                        build_force.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_force.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_force.continuityParameter = obj.continuityParameter;
                        build_force.domainProfile = obj.domainProfile;
                        build_force.domainProfileDer = obj.domainProfileDer;
                        build_force.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildIGAForce(build_force); 
                        
                        forceHistory(:,iteration) = b;

                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);              

                    end
                    
                    %-----------------------------------------------------------------%
                    %                  	      PROBLEM SOLUTION                        %
                    %-----------------------------------------------------------------%
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\b;

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                end
                                
                solveTime = toc;
                
                display = sprintf('Finished Solving the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,solveTime);
                disp(display);
                
                %% MEASUREMENT PARAMETERS

                % Selection matrix
                                
                M = obj.numbVerQuadNodes;
                Nh = obj.numbControlPoints;
                
                % Hs = obj.measureMatrix;
                
                % SELECT EVERY POINT IN THE DOMAIN
                Hs = eye(obj.numbControlPoints * obj.numbVerQuadNodes);
                
                % SELECT A LINE IN THE DOMAIN
                % Hs = zeros(Nh,M * Nh);
                
                % for ii = 1:Nh
                %     Hs(ii, ceil(M/2) * Nh + ii) = 1;
                %  end
                
                % IGA basis and modal basis matrices
                
                [Hb,Hm] = observStates(obj.dimModalBasis,obj.domainLimit_inX,obj.domainLimit_inY,...
                                       obj.stepMeshX,obj.label_upBoundDomain,obj.label_downBoundDomain,...
                                       obj.coefficientForm,obj.degreePolySplineBasis,obj.continuityParameter);
                                   
                % Observation matrix
                
                H = Hs * Hm * Hb;
                
                % Covariance matrices
                
                [numbStates,numbSamples] = size(solutionMatrix);
                
                Q1 = 10^-4;
                R1 = 10^-3;
                
                Q = Q1 * eye(numbStates);
                R = R1 * eye(size(Hs,1));
                
                %% ENSEMBLE PROPERTIES
                %---------------------------------------------------------%
                % Note: Definition of the properties of the random field to
                % be determined to create the samples of the ensemble.
                %---------------------------------------------------------%
                
                % TYPE
                %---------------------------------------------------------%
                % Note: The type of the random field can be set as Gaussian
                % ('gauss'), Exponential ('exp') or Turbulent
                % ('turbulent').
                %---------------------------------------------------------%
                
                corrProperties.name = 'gauss';
                
                % CORRELATION PARAMETER
                %---------------------------------------------------------%
                % Note: The scaling parameters for the correlation
                % function. 'c0' c0 may be a scalar for isotropic
                % correlation or a vector for an anisotropic correlation.
                % In the anisotropic case, the vector must have
                % 'n'elements, where 'n' is the dimension of the state
                % vector.
                %
                % (1) 'c0' =~ 0 (Uncorrelated samples)
                % (2) 'c0' =~ 1 (Completly correlated samples)
                %---------------------------------------------------------%
                
                corrProperties.c0 = 0.0001;
                
                % SAMPLE VARIANCE
                %---------------------------------------------------------%
                % Note: The variance 'sigma' of the samples. May be a 
                % scalar or a vector the size of the state vector.
                %---------------------------------------------------------%

                corrProperties.sigma = 10^-1;
                
                % DATA POINTS
                %---------------------------------------------------------%
                % Note: Data points where the samples will be evaluated.
                % Linear vector going form one to the number of states in
                % the state vector.
                %---------------------------------------------------------%
                
                dataPts = linspace(1,numbStates,numbStates)';
                
                %% EnKF - FILTER ALGORITHM
                %---------------------------------------------------------%
                % 1. Initialization
                % 2. Create state measurements
                % 3. Create initial ensemble set for state
                % 4. Create initial ensemble set for parameters
                % 5. Filtering loop
                %   5.1. Loop over time iterations
                %   5.2. Loop over ensemble set
                %       5.2.1. Compute the system matrix with previous
                %       assmilated state
                %       5.2.2. Compute the prediction for the current
                %       emsemble element
                %       5.2.3. Assimilate the measurements with the
                %       prediction for the current emsemble element
                %---------------------------------------------------------%
                
                tic;
                
                % 1. INITIALIZATION

                Z = zeros(size(Hs,1),numbSamples);
                ensObservation = [];
                ensObsNoise = [];
                currStateMean = zeros(numbStates,1);
                ensObsNoiseMean = zeros(size(Hs,1),1);
                currObsMean = zeros(size(Hs,1),1);
                normalizedState = [];
                normalizedObs = [];
                ensAssmilatedState = [];
                ensPredictionState = [];
                assStateHistory = [];
                condGain = [];
                
                % 2. CREATE STATE MEASUREMENTS
                
                for jj = 1:numbSamples
                    measureStd = 0.05 * max(H * solutionMatrix(:,jj));
                    v = measureStd * randn(size(Hs,1),1);
                    Z(:,jj) = H * solutionMatrix(:,jj) + v;
                    
                end
                
                % 3. CREATE INITIAL ENSEMBLE OF STATES 
                %---------------------------------------------------------%
                % Note: I need to create the object property for the number
                % of ensemble points to be created.
                %---------------------------------------------------------%
                
                [stateEnsembleSet,~] = randomfield(corrProperties,dataPts,...
                                                   'nsamples',obj.numbEnsembles,...
                                                   'mean',obj.initialState);
                
                % 4. CREATE INITIAL ENSEMBLE SET OF PARAMETERS
                %---------------------------------------------------------%
                % Note: This step is necessary when there is parameter
                % identification of the system.
                %---------------------------------------------------------%
                
                % 5. FILTERING LOOP

                for jj = 1:numbSamples-1
                    
                    % DEBUG
                    disp(['FILTER ITERATION : ',num2str(jj)]);
                    
                    %-----------------------------------------------------%
                    % Step 1 : Draw a statistically consistent observation
                    % set with the current observation
                    %-----------------------------------------------------%
                    tic;
                    for ii = 1:obj.numbEnsembles
                        
                        if (jj == 1)                            
                            measureStd = 0.05 * max(H * solutionMatrix(:,jj+1));
                            v = measureStd * randn(size(Hs,1),1);

                            ensObsNoise(:,ii) = v;

                            ensObservation(:,ii) = Z(:,jj) + v;
                        else                        
                            measureStd = 0.05 * max(H * solutionMatrix(:,jj));
                            v = measureStd * randn(size(Hs,1),1);

                            ensObsNoise(:,ii) = v;

                            ensObservation(:,ii) = Z(:,jj) + v;
                        end
                    end
                    tStep1 = toc;
                    disp(['TIME FOR STEP 1 : ',num2str(tStep1)]);
                    
                    %-----------------------------------------------------%
                    % Step 2 : Compute the ensemble means for the state,
                    % expected observation and noise used to generate the
                    % observation set
                    %-----------------------------------------------------%
                    
                    %-----------------------------------------------------%
                    % Note: Do not forget that the matrix containing the
                    % current set of state ensemble is update at each
                    % iteration of the filter.
                    %-----------------------------------------------------%
                    tic;
                    for ii = 1:obj.numbEnsembles
                        currStateMean = currStateMean + stateEnsembleSet(:,ii);
                        ensObsNoiseMean = ensObsNoiseMean + ensObsNoise(:,ii);
                        currObsMean = currObsMean + H * stateEnsembleSet(:,ii);
                    end
                    
                    currStateMean = currStateMean/obj.numbEnsembles;
                    ensObsNoiseMean = ensObsNoiseMean/obj.numbEnsembles;
                    currObsMean = currObsMean/obj.numbEnsembles;
                    
                    tStep2 = toc;
                    disp(['TIME FOR STEP 2 : ',num2str(tStep2)]);
                    
                    %-----------------------------------------------------%
                    % Step 3 : Compute the normalized anomalies
                    %-----------------------------------------------------%
                    tic;
                    for ii = 1:obj.numbEnsembles
                        
                        normalizedState(:,ii) = (stateEnsembleSet(:,ii) - currStateMean)/(sqrt(obj.numbEnsembles - 1));
                        normalizedObs(:,ii) = (H * stateEnsembleSet(:,ii) - ensObsNoise(:,ii) - currObsMean + ensObsNoiseMean)/(sqrt(obj.numbEnsembles - 1));
                        
                    end
                    tStep3 = toc;
                    disp(['TIME FOR STEP 3 : ',num2str(tStep3)]);
                    
                    %-----------------------------------------------------%
                    % Note : It is necessary to store the condition number
                    % of the anomaly matrices to analyse the convergence of
                    % the filter to the actual solution.
                    %-----------------------------------------------------%
                    
                    %-----------------------------------------------------%
                    % Step 4 : Compute the filter gain
                    %-----------------------------------------------------%
                    tic;
                    gainInv = (normalizedObs * normalizedObs')\eye(size(normalizedObs * normalizedObs'));
                    gainEnKF = (normalizedState * normalizedObs') * gainInv;
                    tStep4 = toc;
                    disp(['TIME FOR STEP 4 : ',num2str(tStep4)]);
                    
                    condGain = [condGain cond(gainEnKF,2)];
                    
                    %-----------------------------------------------------%
                    % Step 5 : Update the ensemble set using the weighted
                    % innovation
                    %-----------------------------------------------------%
                    tic;
                    for ii = 1:obj.numbEnsembles
                        
                        ensAssmilatedState(:,ii) = stateEnsembleSet(:,ii) + gainEnKF * (ensObservation(:,ii) - H * stateEnsembleSet(:,ii));
                        
                    end
                    tStep5 = toc;
                    disp(['TIME FOR STEP 5 : ',num2str(tStep5)]);
                    
                    %-----------------------------------------------------%
                    % Step 6 : Compute the new system matrix using the
                    % current state of the system. Necessary to be
                    % performed if a parameter of the system is used as
                    % additional state variable.
                    %-----------------------------------------------------%
                    
                    %-----------------------------------------------------%
                    % Step 7 : Compute the ensemble forecast using the
                    % current dynamic law for the system
                    %-----------------------------------------------------%
                    tic;
                    for ii = 1:obj.numbEnsembles
                        
                        ensPredictionState(:,ii) = stateMatrix * ensAssmilatedState(:,ii) + inputMatrix\forceHistory(:,jj);
                        
                    end
                    
                    stateEnsembleSet = ensPredictionState;
                    assStateHistory = [assStateHistory currStateMean];
                    
                    tStep6 = toc;
                    disp(['TIME FOR STEP 6 : ',num2str(tStep6)]);
                    
                end
                
                filterTime = toc;
                
                display = sprintf('Finished Filtering (EnKF without Parameter Identification)');
                disp(display);
               
                %% DATA EXPORT
                
                %-----------------------------------------------------------------%
                %                  	       DATA EXPORT                            %
                %-----------------------------------------------------------------%
                
                tic;
                
%                 for ii = 1:1000
%         
%                 checkFolder = exist(['Paraview_',num2str(ii)]);
%         
%                     if (checkFolder == 7)
%                         disp(['The folder Paraview_',num2str(ii),' already exists!'])
%                     else
%                         directoryPath = ['Paraview_',num2str(ii)];
%                         mkdir(directoryPath)
%                         break
%                     end
%                 end
                
%                 oldFolder = cd(directoryPath);
%                 
%                 dataFileApprox = cell(numbSamples-1,1);
%                 dataFileNoisy = cell(numbSamples-1,1);
%                 dataFileEstimate = cell(numbSamples-1,1);
% 
%                 errorEvolH1 = zeros(numbSamples-1,1);
%                 errorEvolL2 = zeros(numbSamples-1,1);
%                 errorNoisyH1 = zeros(numbSamples-1,1);
%                 errorNoisyL2 = zeros(numbSamples-1,1);
%                 errorEstimatedH1 = zeros(numbSamples-1,1);
%                 errorEstimatedL2 = zeros(numbSamples-1,1);
                
                %% CREATE EXPORT FOLDER
    
                for ii = 1:1000
        
                checkFolder = exist(['MatlabPlots',num2str(ii)]);
        
                    if (checkFolder == 7)
                        disp(['The folder MatlabPlots',num2str(ii),' already exists!'])
                    else
                        fileNameF = ['MatlabPlots',num2str(ii)];
                        mkdir(fileNameF)
                        break
                    end
                end
                
                %% PRINT SYSTEM EVOLUTION
                
                 for ii = 1:numbIterations-1
 
                 [errL2,errH1,~,~] = plot_KF( ...
                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.domainLimit_inY,obj.stepMeshX, ...
                 solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                 obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer, ...
                 obj.jacAtQuadNodes, obj.degreePolySplineBasis,obj.exactSolution_dX,obj.exactSolution_dY,obj.continuityParameter, ...
                 ii,obj.timeDomain,Z(:,ii),assStateHistory(:,ii),Hs,Hm,Hb);
             
                 errorHistoryL2(ii) = errL2;
                 errorHistoryH1(ii) = errH1;
             
                 end
                 
                 disp('Finished Plot Operation');
                 
                %% CREATE NOISY MEASUREMENT VIDEO
                
                for ii = 1:numbIterations-1
                    fileName = ['Plot_At_it=',num2str(ii),'.png'];
                    imageNames{ii} = fileName;
                end
                 
                currentFolder = pwd;
                workingDir = [currentFolder,'\',fileNameF];
                 
                outputVideo = VideoWriter(fullfile(workingDir,'Assimilated State.avi'));
                outputVideo.FrameRate = 10;
                open(outputVideo)

                for ii = 1:length(imageNames)
                    img = imread(fullfile(workingDir,imageNames{ii}));
                    writeVideo(outputVideo,img)
                end

                close(outputVideo)

                errorNormH1 = 0;
                errorNormL2 = 0;
                errorNoiseNormH1 = 0;
                errorNoiseNormL2 = 0;
                errorEstimateNormH1 = 0;
                errorEstimateNormL2 = 0;

                disp('Finished Computing Error for the IGA Method');

            end % End function
            
            %% Method 'solverEnKFID'
            
            function [solutionMatrix,Z,assStateHistory, ...
                      condGain,filterTime] = solverEnKFID(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                %% IMPORT CLASSES
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% SOLUTION EVOLUTION

                numbIterations = length(obj.timeDomain);
                numbStates = (round(1/obj.stepMeshX(1)) + 1) * obj.dimModalBasis(1);
                
                solutionMatrix = zeros(numbStates,numbIterations);
                solutionMatrix(:,1) = obj.initialState;
                
                forceHistory = zeros(numbStates,numbIterations - 1);
                
                tic;
                
                for iteration = 1 : numbIterations - 1
                    
                    if (iteration == 1)
                
                        % Definition of the Object from the AssemblerADRHandler class

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = obj.coefficientForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildSystemIGATransient(build_IGA);             
                        
                        % Computation of the State Matrix
                        
                        Mdt = M * (obj.timeStep^-1);

                        stateMatrix = (Mdt + A)\Mdt;

                        % Computation of the Input Matrix

                        inputMatrix = (Mdt + A);
                        
                        forceHistory(:,1) = b;
                        
                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);    
                    
                    else
                        
                        % Definition of the Object from the AssemblerADRHandler class

                        build_force = AssemblerADRHandler();

                        % Properties Assignment

                        build_force.dimModalBasis = obj.dimModalBasis(1);
                        build_force.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_force.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_force.stepMeshX = obj.stepMeshX(1);
                        build_force.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_force.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_force.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_force.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_force.coefficientForm = obj.coefficientForm;
                        build_force.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_force.geometricInfo = obj.geometricInfo;
                        build_force.robinCondStruct = obj.robinCondStruct;
                        build_force.couplingCond_DD = obj.couplingCond_DD;
                        build_force.physicMesh_inX = obj.physicMesh_inX;
                        build_force.physicMesh_inY = obj.physicMesh_inY;
                        build_force.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_force.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_force.continuityParameter = obj.continuityParameter;
                        build_force.domainProfile = obj.domainProfile;
                        build_force.domainProfileDer = obj.domainProfileDer;
                        build_force.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildIGAForce(build_force); 
                        
                        forceHistory(:,iteration) = b;

                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);              

                    end
                    
                    %-----------------------------------------------------------------%
                    %                  	      PROBLEM SOLUTION                        %
                    %-----------------------------------------------------------------%
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\b;

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                end
                                
                solveTime = toc;
                
                display = sprintf('Finished Solving the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,solveTime);
                disp(display);
                
                %% MEASUREMENT PARAMETERS

                % Selection matrix
                                
                M = obj.numbVerQuadNodes;
                Nh = obj.numbControlPoints;
                
                % Hs = obj.measureMatrix;
                
                % SELECT EVERY POINT IN THE DOMAIN
                Hs = eye(obj.numbControlPoints * obj.numbVerQuadNodes);
                
                % SELECT A LINE IN THE DOMAIN
                % Hs = zeros(Nh,M * Nh);
                
                % for ii = 1:Nh
                %     Hs(ii, ceil(M/2) * Nh + ii) = 1;
                %  end
                
                % IGA basis and modal basis matrices
                
                [Hb,Hm] = observStates(obj.dimModalBasis,obj.domainLimit_inX,obj.domainLimit_inY,...
                                       obj.stepMeshX,obj.label_upBoundDomain,obj.label_downBoundDomain,...
                                       obj.coefficientForm,obj.degreePolySplineBasis,obj.continuityParameter);
                                   
                % Observation matrix
                
                H = Hs * Hm * Hb;
                
                % Covariance matrices
                
                [numbStates,numbSamples] = size(solutionMatrix);
                
                Q1 = 10^-4;
                R1 = 10^-3;
                
                Q = Q1 * eye(numbStates);
                R = R1 * eye(size(Hs,1));
                
                %% ENSEMBLE PROPERTIES
                %---------------------------------------------------------%
                % Note: Definition of the properties of the random field to
                % be determined to create the samples of the ensemble.
                %---------------------------------------------------------%
                
                % TYPE
                %---------------------------------------------------------%
                % Note: The type of the random field can be set as Gaussian
                % ('gauss'), Exponential ('exp') or Turbulent
                % ('turbulent').
                %---------------------------------------------------------%
                
                corrProperties.name = 'gauss';
                
                % CORRELATION PARAMETER
                %---------------------------------------------------------%
                % Note: The scaling parameters for the correlation
                % function. 'c0' c0 may be a scalar for isotropic
                % correlation or a vector for an anisotropic correlation.
                % In the anisotropic case, the vector must have
                % 'n'elements, where 'n' is the dimension of the state
                % vector.
                %
                % (1) 'c0' =~ 0 (Uncorrelated samples)
                % (2) 'c0' =~ 1 (Completly correlated samples)
                %---------------------------------------------------------%
                
                corrProperties.c0 = 0.0001;
                
                % SAMPLE VARIANCE
                %---------------------------------------------------------%
                % Note: The variance 'sigma' of the samples. May be a 
                % scalar or a vector the size of the state vector.
                %---------------------------------------------------------%

                corrProperties.sigma = 10^-1;
                
                % DATA POINTS
                %---------------------------------------------------------%
                % Note: Data points where the samples will be evaluated.
                % Linear vector going form one to the number of states in
                % the state vector.
                %---------------------------------------------------------%
                
                dataPts = linspace(1,numbStates+obj.numbParam,numbStates+obj.numbParam)';
                
                %% EnKF - FILTER ALGORITHM
                %---------------------------------------------------------%
                % 1. Initialization
                % 2. Create state measurements
                % 3. Create initial ensemble set for state
                % 4. Create initial ensemble set for parameters
                % 5. Filtering loop
                %   5.1. Loop over time iterations
                %   5.2. Loop over ensemble set
                %       5.2.1. Compute the system matrix with previous
                %       assmilated state
                %       5.2.2. Compute the prediction for the current
                %       emsemble element
                %       5.2.3. Assimilate the measurements with the
                %       prediction for the current emsemble element
                %---------------------------------------------------------%
                
                tic;
                
                % 0. AUGMENTED PROPERTIES
                
                numbParameters = obj.numbParam;
                
                Haug = [H, zeros(size(H,1),numbParameters)];
                
                % 1. INITIALIZATION

                Z = zeros(size(Hs,1),numbSamples);
                ensObservation = [];
                ensObsNoise = [];
                currStateMean = zeros(numbStates + numbParameters,1);
                ensObsNoiseMean = zeros(size(Hs,1),1);
                currObsMean = zeros(size(Hs,1),1);
                normalizedState = [];
                normalizedObs = [];
                ensAssmilatedState = [];
                ensPredictionState = [];
                assStateHistory = [];
                condGain = [];
                
                % 2. CREATE STATE MEASUREMENTS
                
                for jj = 1:numbSamples
                    measureStd = 0.05 * max(H * solutionMatrix(:,jj));
                    v = measureStd * randn(size(Hs,1),1);
                    Z(:,jj) = H * solutionMatrix(:,jj) + v;
                    
                end
                
                % 3. CREATE INITIAL ENSEMBLE OF STATES 
                %---------------------------------------------------------%
                % Note: I need to create the object property for the number
                % of ensemble points to be created.
                %---------------------------------------------------------%
                
                initParameters = 0.1 * ones(numbParameters,1);
                initState = [obj.initialState; initParameters];
                assStateHistory = initState;
                
                [stateEnsembleSet,~] = randomfield(corrProperties,dataPts,...
                                                   'nsamples',obj.numbEnsembles,...
                                                   'mean',initState);
                
                % 4. CREATE INITIAL ENSEMBLE SET OF PARAMETERS
                %---------------------------------------------------------%
                % Note: This step is necessary when there is parameter
                % identification of the system.
                %---------------------------------------------------------%
                
                % 5. FILTERING LOOP

                for jj = 1:numbSamples-1
                    
                    % Debug
                    disp(['Simulation in sample : ',num2str(jj)]);
                    
                    %-----------------------------------------------------%
                    % Step 1 : Draw a statistically consistent observation
                    % set with the current observation
                    %-----------------------------------------------------%
                    
                    for ii = 1:obj.numbEnsembles
                        
                        if (jj == 1)                            
                            measureStd = 0.05 * max(H * solutionMatrix(:,jj+1));
                            v = measureStd * randn(size(Hs,1),1);

                            ensObsNoise(:,ii) = v;

                            ensObservation(:,ii) = Z(:,jj) + v;
                        else                        
                            measureStd = 0.05 * max(H * solutionMatrix(:,jj));
                            v = measureStd * randn(size(Hs,1),1);

                            ensObsNoise(:,ii) = v;

                            ensObservation(:,ii) = Z(:,jj) + v;
                        end
                    end
                        
                    %-----------------------------------------------------%
                    % Step 2 : Compute the ensemble means for the state,
                    % expected observation and noise used to generate the
                    % observation set
                    %-----------------------------------------------------%
                    
                    %-----------------------------------------------------%
                    % Note: Do not forget that the matrix containing the
                    % current set of state ensemble is update at each
                    % iteration of the filter.
                    %-----------------------------------------------------%
                    
                    for ii = 1:obj.numbEnsembles
                        currStateMean = currStateMean + stateEnsembleSet(:,ii);
                        ensObsNoiseMean = ensObsNoiseMean + ensObsNoise(:,ii);
                        currObsMean = currObsMean + Haug * stateEnsembleSet(:,ii);
                    end
                    
                    currStateMean = currStateMean/obj.numbEnsembles;
                    ensObsNoiseMean = ensObsNoiseMean/obj.numbEnsembles;
                    currObsMean = currObsMean/obj.numbEnsembles;
                    
                    %-----------------------------------------------------%
                    % Step 3 : Compute the normalized anomalies
                    %-----------------------------------------------------%
                    
                    for ii = 1:obj.numbEnsembles
                        
                        normalizedState(:,ii) = (stateEnsembleSet(:,ii) - currStateMean)/(sqrt(obj.numbEnsembles - 1));
                        normalizedObs(:,ii) = (Haug * stateEnsembleSet(:,ii) - ensObsNoise(:,ii) - currObsMean + ensObsNoiseMean)/(sqrt(obj.numbEnsembles - 1));
                        
                    end
                    
                    %-----------------------------------------------------%
                    % Note : It is necessary to store the condition number
                    % of the anomaly matrices to analyse the convergence of
                    % the filter to the actual solution.
                    %-----------------------------------------------------%
                    
                    %condNormStateAnomaly = [condNormStateAnomaly cond(normalizedState,2)];
                    %condNormObsAnomaly = [condNormObsAnomaly cond(normalizedObs,2)];
                    
                    %-----------------------------------------------------%
                    % Step 4 : Compute the filter gain
                    %-----------------------------------------------------%
                    
                    gainInv = (normalizedObs * normalizedObs')\eye(size(normalizedObs * normalizedObs'));
                    gainEnKF = (normalizedState * normalizedObs') * gainInv;
                    
                    condGain = [condGain cond(normalizedObs * normalizedObs',2)];
                    
                    %-----------------------------------------------------%
                    % Step 5 : Update the ensemble set using the weighted
                    % innovation
                    %-----------------------------------------------------%
                    
                    for ii = 1:obj.numbEnsembles
                        
                        ensAssmilatedState(:,ii) = stateEnsembleSet(:,ii) + gainEnKF * (ensObservation(:,ii) - Haug * stateEnsembleSet(:,ii));
                        
                    end
                    
                    %-----------------------------------------------------%
                    % Step 6 : Compute the new system matrix using the
                    % current state of the system. Necessary to be
                    % performed if a parameter of the system is used as
                    % additional state variable.
                    %-----------------------------------------------------%
                    
                    
                    
                    %-----------------------------------------------------%
                    % Step 7 : Compute the ensemble forecast using the
                    % current dynamic law for the system
                    %-----------------------------------------------------%
                    
                    parfor ii = 1:obj.numbEnsembles
                        
                        disp(['Entered analysis in sample : ',num2str(jj), ' Ensemble (',num2str(ii),')']);
                        tic;
                        
                        % Current parameters
                    
                        mu = @(x,y) ensAssmilatedState(end,ii);
                        beta1 = @(x,y) ensAssmilatedState(end-1,ii);
                        beta2 = @(x,y) ensAssmilatedState(end-2,ii);
                        sigma = @(x,y) ensAssmilatedState(end-3,ii);
                        chi = obj.coefficientForm.coeffrobin;

                        % Pass new parameters to the state matrix

                        coeffForm = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
                                             'sigma',sigma,'coeffrobin',chi);

                        % Compute new state matrix

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = coeffForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = jj-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,~,~,~,~,~,~] = buildSystemIGATransient(build_IGA);             

                        % Computation of the State Matrix

                        Mdt = M * (obj.timeStep^-1);

                        stateMatrix = (Mdt + A)\Mdt;

                        % Computation of the Input Matrix

                        inputMatrix = (Mdt + A);

                        % New Augmented System Matrices

                        augStateMatrix = [stateMatrix zeros(size(stateMatrix,1),numbParameters);...
                                          zeros(numbParameters,size(stateMatrix,2)) eye(numbParameters)];

                        augInputVect = [inputMatrix\forceHistory(:,jj); zeros(numbParameters,1)];
                        
                        ensPredictionState(:,ii) = augStateMatrix * ensAssmilatedState(:,ii) + augInputVect;
                        
                        AssEnsemTime = toc;
                        
                        disp(['Time for ensemble : ',num2str(AssEnsemTime)]);
                        
                    end
                    
                    for ii = 1:obj.numbEnsembles
                        currStateMean = currStateMean + ensPredictionState(:,ii);
                    end
                    
                    currStateMean = currStateMean/obj.numbEnsembles;
                    
                    stateEnsembleSet = ensPredictionState;
                    assStateHistory = [assStateHistory currStateMean];
                    
                end
                
                filterTime = toc;
                
                disp(['TIME REQUIRED TO PERFORM THE EnKF : ',num2str(filterTime)]);
                
                %% TESTS FOR THE PARAMETER ESTIMATION
                
                muReal = obj.coefficientForm.mu(0,0);
                beta1Real = obj.coefficientForm.beta1(0,0);
                beta2Real = obj.coefficientForm.beta2(0,0);
                sigmaReal = obj.coefficientForm.sigma(0,0);
                
                figure;
                for ii = 1:numbParameters
                    plot(assStateHistory(end - ii + 1,:),'-o');
                    hold on;
                end
                plot(muReal*ones(1,numbSamples),'--');
                hold on;
                plot(beta1Real*ones(1,numbSamples),'--');
                hold on;
                plot(beta2Real*ones(1,numbSamples),'--');
                hold on;
                plot(sigmaReal*ones(1,numbSamples),'--');
                hold on;
                
                % legend('$\bar\mu$','$\bar\beta_{1}$','$\bar\beta_{2}$','$\bar\sigma$',...
                %        '$\mu$','$\beta_{1}$','$\beta_{2}$','$\sigma$');
                
                disp(' ');
                disp('FINISHED FILTERING ALGORITHM (EnKF with PARAMETER IDENTIFICATION)');
                
                %% COMPUTATION OF THE ERROR
                
                mu      = assStateHistory(end-0,:);
                beta1   = assStateHistory(end-1,:);
                beta2   = assStateHistory(end-2,:);
                sigma   = assStateHistory(end-3,:);
                
                errorMu = mu - muReal;
                errorBeta1 = beta1 - beta1Real;
                errorBeta2 = beta2 - beta1Real;
                errorSigma = sigma - sigmaReal;
                
                exportTxt(mu, obj.numbEnsembles, filterTime);
                exportTxt(beta1, obj.numbEnsembles, filterTime);
                exportTxt(beta2, obj.numbEnsembles, filterTime);
                exportTxt(sigma, obj.numbEnsembles, filterTime);
                exportTxt(errorMu, obj.numbEnsembles, filterTime);
                exportTxt(errorBeta1, obj.numbEnsembles, filterTime);
                exportTxt(errorBeta2, obj.numbEnsembles, filterTime);
                exportTxt(errorSigma, obj.numbEnsembles, filterTime);

                errorNormH1 = 0;
                errorNormL2 = 0;
                errorNoiseNormH1 = 0;
                errorNoiseNormL2 = 0;
                errorEstimateNormH1 = 0;
                errorEstimateNormL2 = 0;

                disp('FINISHED COMPUTATION OF THE ERROR');

                return
                
                %% DATA EXPORT
                
                %-----------------------------------------------------------------%
                %                  	       DATA EXPORT                            %
                %-----------------------------------------------------------------%
                
                tic;
                
%                 for ii = 1:1000
%         
%                 checkFolder = exist(['Paraview_',num2str(ii)]);
%         
%                     if (checkFolder == 7)
%                         disp(['The folder Paraview_',num2str(ii),' already exists!'])
%                     else
%                         directoryPath = ['Paraview_',num2str(ii)];
%                         mkdir(directoryPath)
%                         break
%                     end
%                 end
                
%                 oldFolder = cd(directoryPath);
%                 
%                 dataFileApprox = cell(numbSamples-1,1);
%                 dataFileNoisy = cell(numbSamples-1,1);
%                 dataFileEstimate = cell(numbSamples-1,1);
% 
%                 errorEvolH1 = zeros(numbSamples-1,1);
%                 errorEvolL2 = zeros(numbSamples-1,1);
%                 errorNoisyH1 = zeros(numbSamples-1,1);
%                 errorNoisyL2 = zeros(numbSamples-1,1);
%                 errorEstimatedH1 = zeros(numbSamples-1,1);
%                 errorEstimatedL2 = zeros(numbSamples-1,1);
                
                %% CREATE EXPORT FOLDER
    
                for ii = 1:1000
        
                checkFolder = exist(['MatlabPlots',num2str(ii)]);
        
                    if (checkFolder == 7)
                        disp(['The folder MatlabPlots',num2str(ii),' already exists!'])
                    else
                        fileNameF = ['MatlabPlots',num2str(ii)];
                        mkdir(fileNameF)
                        break
                    end
                end
                
                %% PRINT SYSTEM EVOLUTION
                
                 for ii = 1:numbIterations-1
 
                 [errL2,errH1,~,~] = plot_KF( ...
                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.domainLimit_inY,obj.stepMeshX, ...
                 solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                 obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer, ...
                 obj.jacAtQuadNodes, obj.degreePolySplineBasis,obj.exactSolution_dX,obj.exactSolution_dY,obj.continuityParameter, ...
                 ii,obj.timeDomain,Z(:,ii),assStateHistory(1:end-4,ii),Hs,Hm,Hb);
             
                 errorHistoryL2(ii) = errL2;
                 errorHistoryH1(ii) = errH1;
             
                 end
                 
                 disp('Finished Plot Operation');
                 
                %% CREATE NOISY MEASUREMENT VIDEO
                
                for ii = 1:numbIterations-1
                    fileName = ['Plot_At_it=',num2str(ii),'.png'];
                    imageNames{ii} = fileName;
                end
                 
                currentFolder = pwd;
                workingDir = [currentFolder,'\',fileNameF];
                 
                outputVideo = VideoWriter(fullfile(workingDir,'Assimilated State.avi'));
                outputVideo.FrameRate = 10;
                open(outputVideo)

                for ii = 1:length(imageNames)
                    img = imread(fullfile(workingDir,imageNames{ii}));
                    writeVideo(outputVideo,img)
                end

                close(outputVideo)
            end % End function
            
            %% Method 'solverUKF'
            
            function [solutionMatrix,Z,Xa,filterTime] = solverUKF(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                %% IMPORT CLASSES
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% SOLUTION EVOLUTION

                numbIterations = length(obj.timeDomain);
                numbStates = (round(1/obj.stepMeshX(1)) + 1) * obj.dimModalBasis(1);
                
                solutionMatrix = zeros(numbStates,numbIterations);
                solutionMatrix(:,1) = obj.initialState;
                
                forceHistory = zeros(numbStates,numbIterations - 1);
                
                tic;
                
                for iteration = 1 : numbIterations - 1
                    
                    if (iteration == 1)
                
                        % Definition of the Object from the AssemblerADRHandler class

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = obj.coefficientForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildSystemIGATransient(build_IGA);             
                        
                        % Computation of the State Matrix
                        
                        Mdt = M * (obj.timeStep^-1);

                        stateMatrix = (Mdt + A)\Mdt;

                        % Computation of the Input Matrix

                        inputMatrix = (Mdt + A);
                        
                        forceHistory(:,1) = b;
                        
                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);    
                    
                    else
                        
                        % Definition of the Object from the AssemblerADRHandler class

                        build_force = AssemblerADRHandler();

                        % Properties Assignment

                        build_force.dimModalBasis = obj.dimModalBasis(1);
                        build_force.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_force.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_force.stepMeshX = obj.stepMeshX(1);
                        build_force.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_force.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_force.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_force.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_force.coefficientForm = obj.coefficientForm;
                        build_force.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_force.geometricInfo = obj.geometricInfo;
                        build_force.robinCondStruct = obj.robinCondStruct;
                        build_force.couplingCond_DD = obj.couplingCond_DD;
                        build_force.physicMesh_inX = obj.physicMesh_inX;
                        build_force.physicMesh_inY = obj.physicMesh_inY;
                        build_force.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_force.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_force.continuityParameter = obj.continuityParameter;
                        build_force.domainProfile = obj.domainProfile;
                        build_force.domainProfileDer = obj.domainProfileDer;
                        build_force.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildIGAForce(build_force); 
                        
                        forceHistory(:,iteration) = b;

                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);              

                    end
                    
                    %-----------------------------------------------------------------%
                    %                  	      PROBLEM SOLUTION                        %
                    %-----------------------------------------------------------------%
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\b;

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                end
                                
                solveTime = toc;
                
                display = sprintf('Finished Solving the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,solveTime);
                disp(display);
                
                %% MEASUREMENT PARAMETERS

                % Selection matrix
                                
                M = obj.numbVerQuadNodes;
                Nh = obj.numbControlPoints;
                
                % Hs = obj.measureMatrix;
                
                % SELECT EVERY POINT IN THE DOMAIN
                Hs = eye(obj.numbControlPoints * obj.numbVerQuadNodes);
                
                % SELECT A LINE IN THE DOMAIN
                % Hs = zeros(Nh,M * Nh);
                
                % for ii = 1:Nh
                %     Hs(ii, ceil(M/2) * Nh + ii) = 1;
                %  end
                
                % IGA basis and modal basis matrices
                
                [Hb,Hm] = observStates(obj.dimModalBasis,obj.domainLimit_inX,obj.domainLimit_inY,...
                                       obj.stepMeshX,obj.label_upBoundDomain,obj.label_downBoundDomain,...
                                       obj.coefficientForm,obj.degreePolySplineBasis,obj.continuityParameter);
                                   
                % Observation matrix
                
                H = Hs * Hm * Hb;
                
                %% SIGMA POINTS PROPERTIES
                %---------------------------------------------------------%
                % Note: Definition of the properties of the random field to
                % be determined to create the samples of the ensemble.
                %---------------------------------------------------------%
                
                % PRIMARY SCALLING PARAMETER
                %---------------------------------------------------------%
                % Note: This parameter defines how the sigma points will
                % spread from the mean vector.
                %---------------------------------------------------------%
                
                alpha = 1e-3;
                
                % SECONDARY SCALLING PARAMETER
                %---------------------------------------------------------%
                % Note: This parameter gives an extra degree of freedom to
                % define the scalling of sigma points. However, it is
                % usually set to zero.
                %---------------------------------------------------------%
                
                kapa = 0;
                
                % PRIOR STATISTICS
                %---------------------------------------------------------%
                % Note: In the case of the state vector being consider a
                % set of Gaussian random variable, the optimal parameter
                % beta is 2.
                %---------------------------------------------------------%
                
                beta = 2; 
                
                % DIMENSION OF THE STATE VECTOR
                %---------------------------------------------------------%
                % Note: Uses the dimension of the state vector to define
                % the number of sigma points to be computed.
                %---------------------------------------------------------%
                
                L = numbStates;
                
                % NUMBER OF SIGMA POINTS
                
                numbSigmaPts = 2*L + 1;
                
                %% UKF - FILTER ALGORITHM
                %---------------------------------------------------------%
                % 1. Initialization
                % 2. Create state measurements
                % 3. Create initial ensemble set for state
                % 4. Create initial ensemble set for parameters
                % 5. Filtering loop
                %   5.1. Loop over time iterations
                %   5.2. Loop over ensemble set
                %       5.2.1. Compute the system matrix with previous
                %       assmilated state
                %       5.2.2. Compute the prediction for the current
                %       emsemble element
                %       5.2.3. Assimilate the measurements with the
                %       prediction for the current emsemble element
                %---------------------------------------------------------%
                
                % 1. INITIALIZATION
 
                [numbStates,numbSamples] = size(solutionMatrix);
                numbMeasurements = size(Hs,1);

                Z = zeros(numbMeasurements,numbSamples);
                Wm = zeros(numbSigmaPts,1);
                Wc = zeros(numbSigmaPts,1);
                
                % 2. CREATE STATE MEASUREMENTS
                
                for jj = 1:numbSamples
                    measureStd = 0.1 * max(H * solutionMatrix(:,jj));
                    v = measureStd * randn(numbMeasurements,1);
                    Z(:,jj) = H * solutionMatrix(:,jj) + v;
                    
                end
                
                % 3. CREATE INITIAL SIGMA POINTS
                %---------------------------------------------------------%
                % Note: Initial sigma points are created for the augmented
                % state. In this case, augmented states mean state vector,
                % process noise vector and measurement noise vector.
                %---------------------------------------------------------%
                
                % PROCCESS COVARIANCE MATRIX
                
                Q1 = 10^-6;  
                Q = Q1 * eye(numbStates);
                
                % MEASUREMENT COVARIANCE MATRIX
                
                R1 = 10^-5;
                R = R1 * eye(numbMeasurements);
                
                % INITIAL STATE COVARIANCE MATRIX
                
                P1 = 10^-2;
                P = P1 * eye(numbStates);
                
                % DEFINITION OF THE AUGMENTED INITIAL STATE AND INITIAL
                % COVARIANCE MATRIX
                
                % augInitState = [obj.initialState; zeros(numbStates,1); zeros(size(Hs,1),1)];
                
                % augInitCovMat = [P, zeros(numbStates,numbStates), zeros(numbStates,numbMeasurements);
                %                  zeros(numbStates,numbStates), Q, zeros(numbStates,numbMeasurements);
                %                  zeros(numbMeasurements,numbStates), zeros(numbMeasurements,numbStates), R];
                
                % 4. COMPUTE THE SIGMA WEIGHTS
                %---------------------------------------------------------%
                % Note: We compute the weights associated with the sigma
                % points just computed. Wm is the set of weights used to
                % compute the mean of the sigma points, Wc is the set of
                % weights used to compute the covariance matrices.
                %---------------------------------------------------------%
                
                lambda = (alpha^2) * (L + kapa) - L;
                
                Wm(1) = lambda/(L + lambda);
                Wc(1) = lambda/(L + lambda) + (1 - alpha^2 + beta);
                
                for ii = 2:numbSigmaPts
                    Wm(ii) = 1/(2 * (L + lambda));
                    Wc(ii) = 1/(2 * (L + lambda));
                end
                
                % 5. FILTERING LOOP
                
                L = size(obj.initialState,1);
                
                sigmaMatrix = zeros(L,numbSigmaPts);
                Xf = zeros(numbStates,numbSigmaPts);
                Yf = zeros(numbMeasurements,numbSigmaPts);
                
                xa = zeros(numbStates,1);
                Pa = P;
                Xa = [xa];
                PaMat = [Pa];
                
                tic;

                for jj = 1:numbSamples-1
                    
                    % DEBUG
                    disp(['FILTER ITERATION : ',num2str(jj)]);
                    
                    %-----------------------------------------------------%
                    % Step 1 : Compute the new sigma points based on the
                    % optimal estimate of the previous step.
                    %-----------------------------------------------------%
                    tic;
                    
                    auxMatrix = chol((L + lambda) * Pa);
                    
                    for ii = 1:numbSigmaPts
                        if(ii == 1)
                            sigmaMatrix(:,ii) = xa;
                        elseif(ii < L+2)
                            sigmaMatrix(:,ii) = xa + auxMatrix(:,ii-1);
                        else
                            sigmaMatrix(:,ii) = xa - auxMatrix(:,ii - L - 1);
                        end
                    end
                    
                    tStep1 = toc;
                    disp(['TIME FOR STEP 1 : ',num2str(tStep1)]);
                    
                    %-----------------------------------------------------%
                    % Note: We have to extract each part of the augmented
                    % state before applying the next steps of the filter.
                    %-----------------------------------------------------%
                    
                    % stateSigma = sigmaMatrix(1:numbStates,:);
                    % processNoiseSigma = sigmaMatrix(numbStates + 1 : 2 * numbStates, :);
                    % measureNoiseSigma = sigmaMatrix(2 * numbStates + 1 : end,:);
                    
                    %-----------------------------------------------------%
                    % Step 2 : Compute the augmented state forecast (Xf,xf)
                    % and statistics (Pf) for each sigma point.
                    %-----------------------------------------------------%
                    
                    % Apply the Unscented Transform to the sigma points
                    
                    tic;
                    for ii = 1:numbSigmaPts
                        
                        Xf(:,ii) = stateMatrix * sigmaMatrix(:,ii) ...
                                   + inputMatrix\forceHistory(:,jj);
                        
                    end
                    tStep2 = toc;
                    disp(['TIME FOR STEP 2 : ',num2str(tStep2)]);
                    
                    % Compute the forecast using the transformed sigma
                    % points
                    
                    tic;
                    sum = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum = sum + Wm(ii) * Xf(:,ii);
                    end
                    
                    xf = sum;
                    
                    tStep3 = toc;
                    disp(['TIME FOR STEP 3 : ',num2str(tStep3)]);
                    
                    % Compute the covariance matrix of the forecast
                    
                    tic;
                    sum = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum = sum + Wc(ii) * (Xf(:,ii) - xf) * (Xf(:,ii) - xf)';
                    end
                    
                    Pf = sum + Q;
                    
                    tStep4 = toc;
                    disp(['TIME FOR STEP 4 : ',num2str(tStep4)]);
                    
                    % DEBUG
                    % disp(xf);
                    % disp(Pf);
               
                    % DEBUG
                    % disp('FORECASTED SIGMA POINTS : ');
                    % disp(Xf(:,1:3));
                    % disp('STATE SIGMA POINTS : ');
                    % disp(sigmaMatrix(:,1:3));
                    % disp('PROCESS NOISE POINTS : ');
                    % disp(processNoiseSigma(:,1:3));
                    % disp('MEASUREMENT NOISE POINTS : ');
                    % disp(measureNoiseSigma(:,1:3));
                    
                    %-----------------------------------------------------%
                    % Step 3 : Compute the expected observations using the
                    % forecast and the law for the observation  system.
                    %-----------------------------------------------------%
                    
                    tic;
                    % Compute the expected observation for each sigma
                    % forecast
                    
                    for ii = 1:numbSigmaPts
                        Yf(:,ii) = H * Xf(:,ii);
                    end
                    
                    % Compute the weighted mean of the expected
                    % observations
                    
                    sum = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum = sum + Wm(ii) * Yf(:,ii);
                    end
                    
                    yf = sum;
                    
                    tStep5 = toc;
                    disp(['TIME FOR STEP 5 : ',num2str(tStep5)]);
                    
                    %-----------------------------------------------------%
                    % Step 4: Compute the assimilated state using the
                    % measurements available.
                    %-----------------------------------------------------%
                    
                    % Compute the residual covariance matrices Pyy and Pxy
                    % that will be used to compute the gain
                    
                    tic;
                    sum1 = 0;
                    sum2 = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum1 = sum1 + Wc(ii) * (Yf(:,ii) - yf) * (Yf(:,ii) - yf)';
                        sum2 = sum2 + Wc(ii) * (Xf(:,ii) - xf) * (Yf(:,ii) - yf)';
                    end
                    
                    Pyy = sum1 + R;
                    Pxy = sum2;
                    tStep6 = toc;
                    disp(['TIME FOR STEP 6 : ',num2str(tStep6)]);
                    
                    % Compute the gain K
                    tic;
                    K = Pxy * (Pyy\eye(size(Pyy)));
                    tStep7 = toc;
                    disp(['TIME FOR STEP 7 : ',num2str(tStep7)]);
                    
                    % Compute the state assimilation from the state
                    % forecast and weighted innovation
                    tic;
                    xa = xf + K * ( Z(:,jj) - yf );
                    Pa = Pf - K * Pyy * K';
                    tStep8 = toc;
                    disp(['TIME FOR STEP 8 : ',num2str(tStep8)]);
                    %-----------------------------------------------------%
                    % Step 5: Store the assimilated state and covariance
                    % matrix.
                    %-----------------------------------------------------%
                    
                    Xa    = [ Xa xa ];
                    PaMat = [ PaMat Pa ]; 
                    
                    %-----------------------------------------------------%
                    % Step 6: Prepare the data for the next time iteration
                    % of the filter.
                    %-----------------------------------------------------%
                    
                end
                
                filterTime = toc;
                
                disp(['TIME REQUIRED TO PERFORM THE UKF : ',num2str(filterTime)]);
                disp(' ');
                disp('FINISHED FILTERING ALGORITHM (UKF without PARAMETER IDENTIFICATION)');

                %% CREATE EXPORT FOLDER
    
                for ii = 1:1000
        
                checkFolder = exist(['MatlabPlots',num2str(ii)]);
        
                    if (checkFolder == 7)
                        disp(['The folder MatlabPlots',num2str(ii),' already exists!'])
                    else
                        fileNameF = ['MatlabPlots',num2str(ii)];
                        mkdir(fileNameF)
                        break
                    end
                end
                
                disp('FINISHED CREATION OF EXPORT FOLDER');
                
                %% PRINT SYSTEM EVOLUTION
                
                 for ii = 1:numbIterations-1
 
                 [errL2,errH1,~,~] = plot_KF( ...
                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.domainLimit_inY,obj.stepMeshX, ...
                 solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                 obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer, ...
                 obj.jacAtQuadNodes, obj.degreePolySplineBasis,obj.exactSolution_dX,obj.exactSolution_dY,obj.continuityParameter, ...
                 ii,obj.timeDomain,Z(:,ii),Xa(:,ii),Hs,Hm,Hb);
             
                 errorHistoryL2(ii) = errL2;
                 errorHistoryH1(ii) = errH1;
             
                 end
                 
                 disp('FINISHED PLOT OPERATION');
                 
                %% CREATE NOISY MEASUREMENT VIDEO
                
                for ii = 1:numbIterations-1
                    fileName = ['Plot_At_it=',num2str(ii),'.png'];
                    imageNames{ii} = fileName;
                end
                
                currentFolder = pwd;
                workingDir = [currentFolder,'\',fileNameF];
                 
                outputVideo = VideoWriter(fullfile(workingDir,'Assimilated State.avi'));
                outputVideo.FrameRate = 10;
                open(outputVideo)

                for ii = 1:length(imageNames)
                    img = imread(fullfile(workingDir,imageNames{ii}));
                    writeVideo(outputVideo,img)
                end

                close(outputVideo)
                
                disp('FINISHED THE CREATION OF SOLUTION VIDEO');
                
                %% COMPUTATION OF THE ERROR

                errorNormH1 = 0;
                errorNormL2 = 0;
                errorNoiseNormH1 = 0;
                errorNoiseNormL2 = 0;
                errorEstimateNormH1 = 0;
                errorEstimateNormL2 = 0;

                disp('FINISHED COMPUTING THE ERROR FOR THE UKF');

            end % End function
            
            %% Method 'solverUKFID'
            
            function [solutionMatrix,Z,Xa,filterTime] = solverUKFID(obj)

                %%
                % solverIGA     - This function solves an elliptic problem definied
                %                 with a dominant dynamic in the X direction. We 
                %                 we assume that there exists Dirichlet boundary 
                %                 conditions at the INFLOW and Neumann boundary
                %                 conditions at the OUTFLOW. The solver receives in
                %                 the input all of the system data and returns the
                %                 numerical solution.
                %                 The solution of the problem is approximated using
                %                 an algorithm of type Domain Decomposition
                %                 (Robin-Robin). In each domain the problem is
                %                 solved using a Hi-Mod method with the
                %                 implementation of the specific assigned basis.
                %
                % The inputs are:
                %%
                %   (1)  domainLimit_inX        : Vector Containing the Extremes of the Domains
                %                                 in the X Direction
                %   (2)  domainLimit_inY        : Vector Containing the Extremes of the Domains
                %                                 in the Y Direction
                %   (3)  dimModalBasis          : Dimension of the Modal Basis in Each Domain
                %   (4)  stepMeshX              : Vector Containing the Step of the Finite
                %                                 Element Mesh
                %   (5)  label_upBoundDomain    : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Upper Limit of
                %                                 he Domain
                %   (6)  label_downBoundDomain  : Contains the Label Identifying the Nature of
                %                                 the Boundary Conditions on the Lower Limit of
                %                                 the Domain
                %   (7)  data_upBoundDomain     : Contains the Values of the Boundary Conditions
                %                                 on the Upper Limit of the Domain
                %   (8)  data_downBoundDomain   : Contains the Values of the Boundary Conditions
                %                                 on the Lower Limit of the Domain
                %   (9)  coefficientForm        : Data Strusture Containing All the @-Functions
                %                                 and the Constants Relative to the Bilinear Form
                %   (10) dirCondFuncStruct      : Data Structure Containing All the @-Functions
                %                                 for the Dirichlet Conditions at the Inflow and
                %                                 for the Exciting Force
                %   (11) geometricInfo          : Data Structure Containing All the
                %                                 Geometric Information regarding the
                %                                 Domain. The current version of the code
                %                                 works only for the specific condition of:
                %                                 (L = 1, a = 0, psi_x = 0)
                %   (12) robinCondStruct        : Data Structure Containing the Two Values of the
                %                                 Coefficients (R, L) for the Robin Condition Used
                %                                 in the Domain Decomposition
                %   (13) dataExportOption       : Labels the Kind of Plot Function that Will Be
                %                                 Used Once the Solution is Compluted
                %   (14) couplingCond_DD        : Contains the Label Adressing the Coupling Condition
                %                                 of the Problem
                %   (15) simulationCase         : Specify What Problem Case is Being Solved
                %   (16) exactSolution          : Exact Solution for the Differential Problem
                %   (17) exactSolution_dX       : Exact Solution for the Differential Problem
                %   (18) exactSolution_dY       : Exact Solution for the Differential Problem
                %   (19) physicMesh_inX         : Vector Containing the Physical Mesh in the X
                %                                 Direction
                %   (20) physicMesh_inY         : Vector Containing the Physical Mesh in the Y
                %                                 Direction
                %   (21) jacAtQuadNodes         : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (22) jacAtQuadNodes2        : Data Sturcture Used to Save the Value of the
                %                                 Jacobians Computed in the Quadratures Nodes 
                %   (23) degreePolySplineBasis  : Degree of the Polynomial B-Spline Basis
                %   (24) continuityParameter    : Degree of Continuity of the Basis 'C^(p-k)'
                %   (25) domainProfile          : Symbolic Function Defining the Profile of the
                %                                 Simulation Domain
                %   (26) domainProfileDer       : Symbolic Function Defining the Derivative of
                %                                 the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) u                : Solution of the Elliptic Problem
                %   (2) liftCoeffA       : Primo Coefficiente di Rilevamento
                %   (3) liftCoeffB       : Secondo Coefficiente di Rilevamento
                %   (4) errorNormL2      : Solution Error on the L2 Norm
                %   (5) errorNormH1      : Solution Error on the H1 Norm
                
                %% IMPORT CLASSES
                
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% SOLUTION EVOLUTION

                numbIterations = length(obj.timeDomain);
                numbStates = (round(1/obj.stepMeshX(1)) + 1) * obj.dimModalBasis(1);
                
                solutionMatrix = zeros(numbStates,numbIterations);
                solutionMatrix(:,1) = obj.initialState;
                
                forceHistory = zeros(numbStates,numbIterations - 1);
                
                tic;
                
                for iteration = 1 : numbIterations - 1
                    
                    if (iteration == 1)
                
                        % Definition of the Object from the AssemblerADRHandler class

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = obj.coefficientForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildSystemIGATransient(build_IGA);             
                        
                        % Computation of the State Matrix
                        
                        Mdt = M * (obj.timeStep^-1);

                        stateMatrix = (Mdt + A)\Mdt;

                        % Computation of the Input Matrix

                        inputMatrix = (Mdt + A);
                        
                        forceHistory(:,1) = b;
                        
                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);    
                    
                    else
                        
                        % Definition of the Object from the AssemblerADRHandler class

                        build_force = AssemblerADRHandler();

                        % Properties Assignment

                        build_force.dimModalBasis = obj.dimModalBasis(1);
                        build_force.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_force.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_force.stepMeshX = obj.stepMeshX(1);
                        build_force.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_force.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_force.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_force.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_force.coefficientForm = obj.coefficientForm;
                        build_force.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_force.geometricInfo = obj.geometricInfo;
                        build_force.robinCondStruct = obj.robinCondStruct;
                        build_force.couplingCond_DD = obj.couplingCond_DD;
                        build_force.physicMesh_inX = obj.physicMesh_inX;
                        build_force.physicMesh_inY = obj.physicMesh_inY;
                        build_force.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_force.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_force.continuityParameter = obj.continuityParameter;
                        build_force.domainProfile = obj.domainProfile;
                        build_force.domainProfileDer = obj.domainProfileDer;
                        build_force.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [b,~,liftCoeffA(1),liftCoeffB(1),~,~] = buildIGAForce(build_force); 
                        
                        forceHistory(:,iteration) = b;

                        % display = ['Finished Time Iteration ',num2str(iteration)];
                        % disp(display);              

                    end
                    
                    %-----------------------------------------------------------------%
                    %                  	      PROBLEM SOLUTION                        %
                    %-----------------------------------------------------------------%
                    
                    stateInovation = stateMatrix * solutionMatrix(:,iteration);
                    inputContrib = inputMatrix\b;

                    solutionMatrix(:,iteration + 1) = stateInovation + inputContrib;
                    
                end
                                
                solveTime = toc;
                
                display = sprintf('Finished Solving the Case %f (Filter %f) in %f',obj.caseSimul,obj.caseFilter,solveTime);
                disp(display);
                
                %% MEASUREMENT PARAMETERS

                % Selection matrix
                                
                M = obj.numbVerQuadNodes;
                Nh = obj.numbControlPoints;
                
                % Hs = obj.measureMatrix;
                
                % SELECT EVERY POINT IN THE DOMAIN
                Hs = eye(obj.numbControlPoints * obj.numbVerQuadNodes);
                
                % SELECT A LINE IN THE DOMAIN
                % Hs = zeros(Nh,M * Nh);
                
                % for ii = 1:Nh
                %     Hs(ii, ceil(M/2) * Nh + ii) = 1;
                %  end
                
                % IGA basis and modal basis matrices
                
                [Hb,Hm] = observStates(obj.dimModalBasis,obj.domainLimit_inX,obj.domainLimit_inY,...
                                       obj.stepMeshX,obj.label_upBoundDomain,obj.label_downBoundDomain,...
                                       obj.coefficientForm,obj.degreePolySplineBasis,obj.continuityParameter);
                                   
                % Observation matrix
                
                H = Hs * Hm * Hb;
                
                %% SIGMA POINTS PROPERTIES
                %---------------------------------------------------------%
                % Note: Definition of the properties of the random field to
                % be determined to create the samples of the ensemble.
                %---------------------------------------------------------%
                
                % PRIMARY SCALLING PARAMETER
                %---------------------------------------------------------%
                % Note: This parameter defines how the sigma points will
                % spread from the mean vector.
                %---------------------------------------------------------%
                
                alpha = 1e-3;
                
                % SECONDARY SCALLING PARAMETER
                %---------------------------------------------------------%
                % Note: This parameter gives an extra degree of freedom to
                % define the scalling of sigma points. However, it is
                % usually set to zero.
                %---------------------------------------------------------%
                
                kapa = 0;
                
                % PRIOR STATISTICS
                %---------------------------------------------------------%
                % Note: In the case of the state vector being consider a
                % set of Gaussian random variable, the optimal parameter
                % beta is 2.
                %---------------------------------------------------------%
                
                beta = 2;
                
                %% UKF - FILTER ALGORITHM
                %---------------------------------------------------------%
                % 1. Initialization
                % 2. Create state measurements
                % 3. Create initial ensemble set for state
                % 4. Create initial ensemble set for parameters
                % 5. Filtering loop
                %   5.1. Loop over time iterations
                %   5.2. Loop over ensemble set
                %       5.2.1. Compute the system matrix with previous
                %       assmilated state
                %       5.2.2. Compute the prediction for the current
                %       emsemble element
                %       5.2.3. Assimilate the measurements with the
                %       prediction for the current emsemble element
                %---------------------------------------------------------%
                
                tic;
                
                % 1. INITIALIZATION
 
                [numbStates,numbSamples] = size(solutionMatrix);
                numbMeasurements = size(Hs,1);
                
                evalMu = obj.coefficientForm.mu(0,0);
                evalBeta1 = obj.coefficientForm.beta1(0,0);
                evalBeta2 = obj.coefficientForm.beta2(0,0);
                evalSigma = obj.coefficientForm.sigma(0,0);
                
                if(not(evalMu == 0))
                    flagMu = 1;
                else
                    flagMu = 0;
                end
                if(not(evalBeta1 == 0))
                    flagBeta1 = 1;
                else
                    flagBeta1 = 0;
                end
                if(not(evalBeta2 == 0))
                    flagBeta2 = 1;
                else
                    flagBeta2 = 0;
                end
                if(not(evalSigma == 0))
                    flagSigma = 1;
                else
                    flagSigma = 0;
                end
                
                numbParameters = flagMu + flagBeta1 + flagBeta2 + flagSigma;
                numbParameters = 4;
                                
                % DEBUG
                disp(['NUMBER OF PARAMETERS : ',num2str(numbParameters)]);
                disp(['FLAG MU    : ',num2str(flagMu)]);
                disp(['FLAG BETA1 : ',num2str(flagBeta1)]);
                disp(['FLAG BETA2 : ',num2str(flagBeta2)]);
                disp(['FLAG SIGMA : ',num2str(flagSigma)]);

                L = numbStates + numbParameters;
                numbSigmaPts = 2*L + 1;

                Z = zeros(numbMeasurements,numbSamples);
                Wm = zeros(numbSigmaPts,1);
                Wc = zeros(numbSigmaPts,1);
                
                % 2. CREATE STATE MEASUREMENTS
                
                for jj = 1:numbSamples
                    measureStd = 0.05 * max(H * solutionMatrix(:,jj));
                    v = measureStd * randn(numbMeasurements,1);
                    Z(:,jj) = H * solutionMatrix(:,jj) + v;
                    
                end
                
                % 3. CREATE INITIAL SIGMA POINTS
                %---------------------------------------------------------%
                % Note: Initial sigma points are created for the augmented
                % state. In this case, augmented states mean state vector,
                % process noise vector and measurement noise vector.
                %---------------------------------------------------------%
                
                % PROCCESS COVARIANCE MATRIX
                
                stateQ1 = 10^-4;  
                stateQ = stateQ1 * eye(numbStates);
                
                paramQ1 = 1e-2;
                paramQ = paramQ1 * eye(numbParameters);
                
                Q = [stateQ, zeros(numbStates,numbParameters);
                     zeros(numbParameters,numbStates), paramQ];
                % MEASUREMENT COVARIANCE MATRIX
                
                R1 = 10^-3;
                R = R1 * eye(numbMeasurements);
                
                % INITIAL STATE COVARIANCE MATRIX
                
                stateP1 = 10^-4;
                stateP = stateP1 * eye(numbStates);
                
                paramP1 = 10^(-2.5);
                paramP = paramP1 * eye(numbParameters);
                
                P = [stateP, zeros(numbStates,numbParameters);
                     zeros(numbParameters,numbStates), paramP];
                
                % DEFINITION OF THE AUGMENTED INITIAL STATE AND INITIAL
                % COVARIANCE MATRIX
                
                % augInitState = [obj.initialState; zeros(numbStates,1); zeros(size(Hs,1),1)];
                
                % augInitCovMat = [P, zeros(numbStates,numbStates), zeros(numbStates,numbMeasurements);
                %                  zeros(numbStates,numbStates), Q, zeros(numbStates,numbMeasurements);
                %                  zeros(numbMeasurements,numbStates), zeros(numbMeasurements,numbStates), R];
                
                % 4. COMPUTE THE SIGMA WEIGHTS
                %---------------------------------------------------------%
                % Note: We compute the weights associated with the sigma
                % points just computed. Wm is the set of weights used to
                % compute the mean of the sigma points, Wc is the set of
                % weights used to compute the covariance matrices.
                %---------------------------------------------------------%
                
                % DEBUG
                % disp('L : ');
                % disp(L);
                % disp('NUMB STATES : ');
                % disp(numbStates);
                % disp('NUMB PARAMETERS : ');
                % disp(numbParameters);
                
                lambda = (alpha^2) * (L + kapa) - L;
                
                Wm(1) = lambda/(L + lambda);
                Wc(1) = lambda/(L + lambda) + (1 - alpha^2 + beta);
                
                for ii = 2:numbSigmaPts
                    Wm(ii) = 1/(2 * (L + lambda));
                    Wc(ii) = 1/(2 * (L + lambda));
                end
                
                % AUGMENTED SYSTEM PROPERTIES
                
                initState = ones(numbStates,1);
                initParameters = 0.1 * ones(numbParameters,1);
                augInitState = [initState; initParameters];
            
                % 5. FILTERING LOOP
                
                sigmaMatrix = zeros(L,numbSigmaPts);
                Xf = zeros(L,numbSigmaPts);
                Yf = zeros(numbMeasurements,numbSigmaPts);
                
                xa = augInitState;
                Pa = P;
                Xa = [xa];
                PaMat = [Pa];
                
                tic;

                for jj = 1:numbSamples-1
                    
                    
                    % DEBUG
                    disp(['FILTER ITERATION : ',num2str(jj)]);
                    
                    %-----------------------------------------------------%
                    % Step 1 : Compute the new sigma points based on the
                    % optimal estimate of the previous step.
                    %-----------------------------------------------------%
                    
                    auxMatrix = sqrt(L + lambda) * chol(Pa);
                    
                    for ii = 1:numbSigmaPts
                        if(ii == 1)
                            sigmaMatrix(:,ii) = xa;
                        elseif(ii < L+2)
                            sigmaMatrix(:,ii) = xa + auxMatrix(:,ii-1);
                        else
                            sigmaMatrix(:,ii) = xa - auxMatrix(:,ii - L - 1);
                        end
                    end
                    
                    %-----------------------------------------------------%
                    % Note: We have to extract each part of the augmented
                    % state before applying the next steps of the filter.
                    %-----------------------------------------------------%
                    
                    % stateSigma = sigmaMatrix(1:numbStates,:);
                    % processNoiseSigma = sigmaMatrix(numbStates + 1 : 2 * numbStates, :);
                    % measureNoiseSigma = sigmaMatrix(2 * numbStates + 1 : end,:);
                    
                    %-----------------------------------------------------%
                    % Step 2 : Compute the augmented state forecast (Xf,xf)
                    % and statistics (Pf) for each sigma point.
                    %-----------------------------------------------------%
                    
                    % Apply the Unscented Transform to the sigma points
                    
                    parfor ii = 1:numbSigmaPts
                        
                        % Current parameters
                        
                        mu = @(x,y) sigmaMatrix(end,ii);
                        beta1 = @(x,y) sigmaMatrix(end-1,ii);
                        beta2 = @(x,y) sigmaMatrix(end-2,ii);
                        sigma = @(x,y) sigmaMatrix(end-3,ii);
                        
                        chi = obj.coefficientForm.coeffrobin;

                        % Pass new parameters to the state matrix

                        coeffForm = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
                                             'sigma',sigma,'coeffrobin',chi);
                                         
                        % Definition of the Object from the AssemblerADRHandler class

                        build_IGA = AssemblerADRHandler();

                        % Properties Assignment

                        build_IGA.dimModalBasis = obj.dimModalBasis(1);
                        build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                        build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                        build_IGA.stepMeshX = obj.stepMeshX(1);
                        build_IGA.label_upBoundDomain = obj.label_upBoundDomain{1};
                        build_IGA.label_downBoundDomain = obj.label_downBoundDomain{1};
                        build_IGA.localdata_upBDomain = obj.data_upBoundDomain{1};
                        build_IGA.localdata_downBDomain = obj.data_downBoundDomain{1};
                        build_IGA.coefficientForm = coeffForm;
                        build_IGA.dirCondFuncStruct = obj.dirCondFuncStruct;
                        build_IGA.geometricInfo = obj.geometricInfo;
                        build_IGA.robinCondStruct = obj.robinCondStruct;
                        build_IGA.couplingCond_DD = obj.couplingCond_DD;
                        build_IGA.physicMesh_inX = obj.physicMesh_inX;
                        build_IGA.physicMesh_inY = obj.physicMesh_inY;
                        build_IGA.jacAtQuadNodes = obj.jacAtQuadNodes;
                        build_IGA.degreePolySplineBasis = obj.degreePolySplineBasis;
                        build_IGA.continuityParameter = obj.continuityParameter;
                        build_IGA.domainProfile = obj.domainProfile;
                        build_IGA.domainProfileDer = obj.domainProfileDer;
                        build_IGA.timeInstant = iteration-1;

                        % Call of the 'buildSystemIGA' Method

                        [A,M,~,~,~,~,~,~] = buildSystemIGATransient(build_IGA);             
                        
                        % Computation of the State Matrix
                        
                        Mdt = M * (obj.timeStep^-1);
                        stateMatrix = (Mdt + A)\Mdt;       
                        inputMatrix = (Mdt + A);
                        
                        % Compute the augmented system components
                        
                        augStateMat  = [stateMatrix, zeros(numbStates,numbParameters);
                                        zeros(numbParameters,numbStates), eye(numbParameters)];
                        augForceVect = [inputMatrix\forceHistory(:,jj); zeros(numbParameters,1)];
                                         
                        % Compute the state forecast for each sigma point
                        
                        Xf(:,ii) = augStateMat * sigmaMatrix(:,ii) + augForceVect;
                        
                    end
                    
                    % Compute the forecast mean using the transformed sigma
                    % points
                    
                    sum = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum = sum + Wm(ii) * Xf(:,ii);
                    end
                    
                    xf = sum;
                    
                    % Compute the covariance matrix of the forecast
                    
                    sum = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum = sum + Wc(ii) * (Xf(:,ii) - xf) * (Xf(:,ii) - xf)';
                    end
                    
                    Pf = sum + Q;
                    
                    % DEBUG
                    % disp(xf);
                    % disp(Pf);
               
                    % DEBUG
                    % disp('FORECASTED SIGMA POINTS : ');
                    % disp(Xf(:,1:3));
                    % disp('STATE SIGMA POINTS : ');
                    % disp(sigmaMatrix(:,1:3));
                    % disp('PROCESS NOISE POINTS : ');
                    % disp(processNoiseSigma(:,1:3));
                    % disp('MEASUREMENT NOISE POINTS : ');
                    % disp(measureNoiseSigma(:,1:3));
                    
                    %-----------------------------------------------------%
                    % Step 3 : Compute the expected observations using the
                    % forecast and the law for the observation  system.
                    %-----------------------------------------------------%
                    
                    % Compute the expected observation for each sigma
                    % forecast
                    
                    for ii = 1:numbSigmaPts
                        Yf(:,ii) = [H , zeros(numbMeasurements,numbParameters)] * Xf(:,ii);
                    end
                    
                    % Compute the weighted mean of the expected
                    % observations
                    
                    sum = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum = sum + Wm(ii) * Yf(:,ii);
                    end
                    
                    yf = sum;
                    
                    %-----------------------------------------------------%
                    % Step 4: Compute the assimilated state using the
                    % measurements available.
                    %-----------------------------------------------------%
                    
                    % Compute the residual covariance matrices Pyy and Pxy
                    % that will be used to compute the gain
                    
                    sum1 = 0;
                    sum2 = 0;
                    
                    for ii = 1:numbSigmaPts
                        sum1 = sum1 + Wc(ii) * (Yf(:,ii) - yf) * (Yf(:,ii) - yf)';
                        sum2 = sum2 + Wc(ii) * (Xf(:,ii) - xf) * (Yf(:,ii) - yf)';
                    end
                    
                    Pyy = sum1 + R;
                    Pxy = sum2;
                    
                    % Compute the gain K
                    
                    K = Pxy * (Pyy\eye(size(Pyy)));
                    
                    % Compute the state assimilation from the state
                    % forecast and weighted innovation
                    
                    xa = xf + K * ( Z(:,jj) - yf );
                    Pa = Pf - K * Pyy * K';
                    
                    % Conoditioning the covariance matrix of the
                    % assimilated state
                    
                    % epsilon = 1e-10;
                    % aux1 = 0.5 * Pa + 0.5 * Pa';
                    % aux2 = aux1 + epsilon * eye(size(aux1));
                    % Pa = aux2;
                    % disp('EIGENVALUES OF Pa');
                    % disp(eig(Pa));
                    
                    %-----------------------------------------------------%
                    % Step 5: Store the assimilated state and covariance
                    % matrix.
                    %-----------------------------------------------------%
                    
                    Xa    = [ Xa xa ];
                    PaMat = [ PaMat Pa ]; 
                    
                    %-----------------------------------------------------%
                    % Step 6: Prepare the data for the next time iteration
                    % of the filter.
                    %-----------------------------------------------------%
                    
                end
                
                filterTime = toc;
                
                %% TESTS FOR THE PARAMETER ESTIMATION
                
                muReal = obj.coefficientForm.mu(0,0);
                beta1Real = obj.coefficientForm.beta1(0,0);
                beta2Real = obj.coefficientForm.beta2(0,0);
                sigmaReal = obj.coefficientForm.sigma(0,0);
                
                figure;
                for ii = 1:numbParameters
                    plot(Xa(end - ii + 1,:),'-o');
                    hold on;
                end
                if (flagMu == 1)
                    plot(muReal*ones(1,numbSamples),'--');
                    hold on;
                end
                if (flagBeta1 == 1)
                    plot(beta1Real*ones(1,numbSamples),'--');
                    hold on;
                end
                if (flagBeta2 == 1)
                    plot(beta2Real*ones(1,numbSamples),'--');
                    hold on;
                end
                if (flagSigma == 1)
                    plot(sigmaReal*ones(1,numbSamples),'--');
                    hold on;
                end
                %legend('$\bar\mu$','$\bar\beta_{1}$','$\bar\beta_{2}$','$\bar\sigma$',...
                %       '$\mu$','$\beta_{1}$','$\beta_{2}$','$\sigma$');
                
                disp(['TIME REQUIRED TO PERFORM THE UKF : ',num2str(filterTime)]);
                disp(' ');
                disp('FINISHED FILTERING ALGORITHM (UKF without PARAMETER IDENTIFICATION)');
                
                %% COMPUTATION OF THE ERROR
                
                mu      = Xa(end-0,:);
                beta1   = Xa(end-1,:);
                beta2   = Xa(end-2,:);
                sigma   = Xa(end-3,:);
                
                errorMu = mu - muReal;
                errorBeta1 = beta1 - beta1Real;
                errorBeta2 = beta2 - beta1Real;
                errorSigma = sigma - sigmaReal;
                
                exportTxt(mu, numbSigmaPts, filterTime);
                exportTxt(beta1, numbSigmaPts, filterTime);
                exportTxt(beta2, numbSigmaPts, filterTime);
                exportTxt(sigma, numbSigmaPts, filterTime);
                exportTxt(errorMu, numbSigmaPts, filterTime);
                exportTxt(errorBeta1, numbSigmaPts, filterTime);
                exportTxt(errorBeta2, numbSigmaPts, filterTime);
                exportTxt(errorSigma, numbSigmaPts, filterTime);

                errorNormH1 = 0;
                errorNormL2 = 0;
                errorNoiseNormH1 = 0;
                errorNoiseNormL2 = 0;
                errorEstimateNormH1 = 0;
                errorEstimateNormL2 = 0;

                disp('FINISHED COMPUTING THE ERROR FOR THE UKF');
                
                return
                
                %% CREATE EXPORT FOLDER
    
                for ii = 1:1000
        
                checkFolder = exist(['MatlabPlots',num2str(ii)]);
        
                    if (checkFolder == 7)
                        disp(['The folder MatlabPlots',num2str(ii),' already exists!'])
                    else
                        fileNameF = ['MatlabPlots',num2str(ii)];
                        mkdir(fileNameF)
                        break
                    end
                end
                
                disp('FINISHED CREATION OF EXPORT FOLDER');
                
                %% PRINT SYSTEM EVOLUTION
                
                 for ii = 1:numbIterations-1
 
                 [errL2,errH1,~,~] = plot_KF( ...
                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.domainLimit_inY,obj.stepMeshX, ...
                 solutionMatrix(:,ii),obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                 obj.coefficientForm,obj.simulationCase,obj.exactSolution,obj.domainProfile,obj.domainProfileDer, ...
                 obj.jacAtQuadNodes, obj.degreePolySplineBasis,obj.exactSolution_dX,obj.exactSolution_dY,obj.continuityParameter, ...
                 ii,obj.timeDomain,Z(:,ii),Xa(1:end-numbParameters,ii),Hs,Hm,Hb);
             
                 errorHistoryL2(ii) = errL2;
                 errorHistoryH1(ii) = errH1;
             
                 end
                 
                 disp('FINISHED PLOT OPERATION');
                 
                %% CREATE NOISY MEASUREMENT VIDEO
                
                for ii = 1:numbIterations-1
                    fileName = ['Plot_At_it=',num2str(ii),'.png'];
                    imageNames{ii} = fileName;
                end
                
                currentFolder = pwd;
                workingDir = [currentFolder,'\',fileNameF];
                 
                outputVideo = VideoWriter(fullfile(workingDir,'Assimilated State.avi'));
                outputVideo.FrameRate = 10;
                open(outputVideo)

                for ii = 1:length(imageNames)
                    img = imread(fullfile(workingDir,imageNames{ii}));
                    writeVideo(outputVideo,img)
                end

                close(outputVideo)
                
                disp('FINISHED THE CREATION OF SOLUTION VIDEO');

            end % End function

    end
end