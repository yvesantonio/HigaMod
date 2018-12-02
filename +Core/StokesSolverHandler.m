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
            
            discStruct      % Structure containing the discretization parameters
            
            boundCondStruct % Structure containing the isogeometric basis
                            % paramters
                            
            probParameters  % Structure containing the problem parameters
            
            timeStruct      % Structure containing the time domain and time
                            % simulation parameters
                            
            quadProperties  % Structure containing the quadrature parameters
            
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

                [errL2,errH1] = plot_solution_IGA_scatter( ...
                obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
                obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map);

                errorNormH1 = errL2;
                errorNormL2 = errH1;
                
                disp('Finished PLOT OPERATION / ERROR with EXACT SOLUTION')
            
%                 [errL2,errH1] = computeErrorIGA_scatter( ...
%                 obj.dimModalBasis,liftCoeffA,liftCoeffB,obj.domainLimit_inX,obj.stepMeshX, u,obj.label_upBoundDomain,obj.label_downBoundDomain, ...
%                 obj.coefficientForm,obj.simulationCase,obj.degreePolySplineBasis,obj.continuityParameter,space,refDomain1D,obj.geometricInfo.map);
%              
%                 disp('Finished ERROR with FREEFEM++ SOLUTION')
% 
%                 errorNormH1 = errH1;
%                 errorNormL2 = errL2;

                disp('Finished Method SOLVER IGA')

            end % End function

    end
end