classdef PODSolverHandler < handle
    
    %% POD SOLVER HANDLER CLASS
    % The PODSolverHandler is a class that contains all of the scripts
    % responsible for the HiPOD model reduction of the differential
    % problems proposed in the examples. It calls all of the other basic
    % classes and functions already incorporated.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% POD SOLVER HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the SolverHandler
        % encapsulate all of the variables needed to run the methods 
        % associated with the computation of the solution for the
        % differential problems proposed. They basically represent all of
        % the variables required to define completely the equations, domain
        % and boundary conditions of the differential problem.
        % The properties are the same of SolverHandler, execept for
        %  * coefficientFormExpansion, which replaces SolverHandler.coefficientForm,
        %    which stores the affine expansion of the problem coefficients as
        %     {{@(param) function of param, @(x,y) function of x and y},
        %      {@(param) function of param, @(x,y) function of x and y},
        %      ...
        %      {@(param) function of param, @(x,y) function of x and y}}
        %  * dirCondFuncStructExpansion, which replaces SolverHandler.dirCondFuncStruct
        %  * exactSolution, exactSolution_dX and exactSolution_dY are not present,
        %    as the error analysis of the reduced model is carried out with respect
        %    to the HiMod solution, rather than the exact one.
        % and stores . These
        % properties will then be internally employed to create a SolverHandler class.
        
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
                          
            coefficientFormExpansion; % Data Strusture Containing All the @-Functions
                                      % and the Constants Relative to the Bilinear Form
                          
            dirCondFuncStructExpansion; % Data Structure Containing All the @-Functions
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
            parameterRange % range of variation for parameters
                           % [interval_1, interval_2, ..., interval_p]
                           % where p is the length of the parameter vector and
                           % interval_i = [a_i, b_i] is the range of variation for the
                           % i-th parameter
            Z % basis functions matrix
            reducedLHSExpansion % left-hand side expansion after projection on the reduced basis
            reducedRHSExpansion % left-hand side expansion after projection on the reduced basis
    end
    
    methods (Access = public)
        
        %% POD SOLVER HANDLER - CONSTRUCT METHOD
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
        
        %% POD SOLVER HANDLER - FUNCTIONAL METHODS
        % The following methods correspond to the construction (offline)
        % and evaluation (online) stages of the ROM.
        
            %% Method 'offline'.
            
            function [D] = offline(obj, nSnapshots, N)

                %%
                % offline
                % Carry out the offline stage of the HiPOD reduced order model.
                %
                % Class properties employed in this method are the same as 
                % SolverHandler.solverIGAScatter(), with the following exceptions:
                %   (*) coefficientFormExpansion replaces coefficientForm;
                %   (*) dirCondFuncStructExpansion replaces dirCondFuncStruct;
                % and are contained in the first input `obj'.
                % The remaining inputs are as follows:
                %   (2) nSnapshots: number of snapshots to build the correlation matrix
                %   (3) N: dimension of the reduced basis
                %
                % The outputs are:
                %%
                %   (1) D                : eigenvalues of the correlatio matrix
                
                
                % Generate deterministic training set
                p = size(obj.parameterRange, 1);
                if p == 1
                    trainingSet = num2cell(linspace(obj.parameterRange(1, 1), obj.parameterRange(1, 2), nSnapshots));
                else
                    error('Equispaced generation not implemented yet for p > 1')
                end
                
                % Compute snapshots matrix
                snpashotsMatrix = [];
                for i = 1:length(trainingSet)
                    ui = truthSolve(obj, trainingSet{i});
                    snpashotsMatrix = [snpashotsMatrix, ui];
                end
                
                % Compress by POD
                correlationMatrix = snpashotsMatrix'*snpashotsMatrix;
                [V, D] = eig(correlationMatrix, 'vector');
                [D, reorder] = sort(D, 'descend');
                V = V(:, reorder);
                V = V(:, 1:N);
                obj.Z = snpashotsMatrix*V;
                for i = 1:N
                    obj.Z(:, i) = 1/D(i)*obj.Z(:, i);
                end                
                % Precompute affine expansion storage
                [LHSExpansion, RHSExpansion] = truthAssemble(obj);
                obj.reducedLHSExpansion = {};
                for q = 1:length(LHSExpansion)
                    obj.reducedLHSExpansion{q} = {...
                        LHSExpansion{q}{1}, ...
                        obj.Z'*LHSExpansion{q}{2}*obj.Z ...
                    };
                end
                obj.reducedRHSExpansion = {};
                for q = 1:length(RHSExpansion)
                    obj.reducedRHSExpansion{q} = {...
                        RHSExpansion{q}{1}, ...
                        obj.Z'*RHSExpansion{q}{2} ...
                    };
                end                
            end % End function
            
            %% Method 'offline'.
            
            function [errors, speedups] = online(obj, nTest)

                %%
                % online
                % Carry out the online stage (with error analysis) of the HiPOD reduced order model.
                %
                % Class properties employed in this method are the same as 
                % SolverHandler.solverIGAScatter(), with the following exceptions:
                %   (*) coefficientFormExpansion replaces coefficientForm;
                %   (*) dirCondFuncStructExpansion replaces dirCondFuncStruct;
                % and are contained in the first input `obj'.
                % The remaining inputs are as follows:
                %   (2) nTests: number of elements in the testing set
                %
                % The outputs are:
                %%
                %   (1) errors                : average error on the testing set
                %   (1) speedups                : average speedup on the testing set
                
                % Generate random testing set
                p = size(obj.parameterRange, 1);
                if p == 1
                    testingSet = num2cell(obj.parameterRange(1, 1) + (obj.parameterRange(1, 2) - obj.parameterRange(1, 1))*rand(nTest, 1));
                else
                    error('Random generation not implemented yet for p > 1')
                end
                
                % Reduced basis dimension
                N = size(obj.Z, 2);
                
                % Storage for errors and speedups
                all_errors = zeros(length(testingSet), N);
                all_speedups = zeros(length(testingSet), N);
                
                % Carry out error analysis
                for i = 1:length(testingSet)
                    tic;
                    ui = truthSolve(obj, testingSet{i});
                    truthTime = toc;
                    for n = 1:N
                        tic;
                        uiReduced = reducedSolve(obj, testingSet{i}, n);
                        reducedTime = toc;
                        all_speedups(i, n) = truthTime/reducedTime;
                        uiReduced = obj.Z(:, 1:n)*uiReduced;
                        all_errors(i, n) = norm(ui - uiReduced)/norm(ui);
                    end
                end
                
                % Average errors and speedups over testing set
                errors = zeros(N, 1);
                speedups = zeros(N, 1);
                for n = 1:N
                    errors(n) = mean(all_errors(:, n));
                    speedups(n) = mean(all_speedups(:, n));
                end
                
            end % End function
    end
    
    methods (Access = private)
        
        %% POD SOLVER HANDLER - PRIVATE FUNCTIONAL METHODS
        % The following internal methods are called by the public methods.
        
            %% Method 'truthSolve'.
            
            function [u] = truthSolve(obj, param)

                %%
                % truthSolve
                % Perform an HiMod truth solve by properly
                % setting up a SolverHandler and assigning all require properties
                % based on the current value of the parameter mu
                %
                % Class properties employed in this method are the same as 
                % SolverHandler.solverIGAScatter(), with the following exceptions:
                %   (*) coefficientFormExpansion replaces coefficientForm;
                %   (*) dirCondFuncStructExpansion replaces dirCondFuncStruct;
                % and are contained in the first input `obj'.
                % Furthermore, the second input `param' contains the current value of the
                % parameters, and is a vector of real numbers.
                %
                % The outputs are:
                %%
                %   (1) u                : HiMod solution of the Elliptic Problem
                
                import Core.AssemblerADRHandler
                
                disp(['Truth solve for param=', num2str(param)])

                %---------------------------------------------------------------------%
                %               ASSEMBLING OF THE COMPLETE LINEAR SYSTEM              %
                %---------------------------------------------------------------------%

                % Definition of the Object from the AssemblerADRHandler class

                build_IGA = AssemblerADRHandler();

                % Properties Assignment
                
                % 1) Evaluate the affine expansion for the current value of the parameter mu
                % PDE coefficients on the left-hand or right-hand side might be parametrized,
                % and the corresponding function handles are summed after prescribing
                % the current value of the parameter.
                % In contrast, boundary conditions are not supposed to be parametrized,
                % so simply copy them.
                mu = sumHandles(obj.coefficientFormExpansion.mu, param);
                beta1 = sumHandles(obj.coefficientFormExpansion.beta1, param);
                beta2 = sumHandles(obj.coefficientFormExpansion.beta2, param);
                sigma = sumHandles(obj.coefficientFormExpansion.sigma, param);
                coeffrobin = obj.coefficientFormExpansion.coeffrobin;
                build_IGA.coefficientForm = struct('mu', mu, 'beta1', beta1, 'beta2', beta2, 'sigma', sigma, 'coeffrobin', coeffrobin);
                igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                force = sumHandles(obj.dirCondFuncStructExpansion.force, param);
                build_IGA.dirCondFuncStruct = struct('igaBoundCond', igaBoundCond, 'force', force);
                
                % 2) Copying the rest of the data from obj.

                build_IGA.dimModalBasis = obj.dimModalBasis(1);
                build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                build_IGA.stepMeshX = obj.stepMeshX(1);
                build_IGA.label_upBoundDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_UP_TAG;
                build_IGA.label_downBoundDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_DOWN_TAG;
                build_IGA.localdata_upBDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_UP_DATA;
                build_IGA.localdata_downBDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_DOWN_DATA;
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
                build_IGA.igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                
                % Call of the 'buildSystemIGA' Method

                [A,b,~,~,~,~,~,~] = buildSystemIGAScatter(build_IGA);
                
                %-----------------------------------------------------------------%
                %                  	      PROBLEM SOLUTION                        %
                %-----------------------------------------------------------------%

                u = A\b;
                
                % Note that we do not rebuild the solution to restore vector entries
                % related to boundary conditions, as we will need snapshots to be
                % of the same size of assembled tensors A and b.
                
            end % End function
            
            %% Method 'truthAssemble'.
            
            function [LHSExpansion, RHSExpansion] = truthAssemble(obj)

                %%
                % truthAssemble
                % Pre-assemble parameter independent matrices
                % on the left-hand side, as well as as parameter independent vectors
                % on the right-hand side.
                % TODO Dirichlet boundary conditions (if any) are always assumed
                %      to be homogeneous
                %
                % Class properties employed in this method are the same as 
                % SolverHandler.solverIGAScatter(), with the following exceptions:
                %   (*) coefficientFormExpansion replaces coefficientForm;
                %   (*) dirCondFuncStructExpansion replaces dirCondFuncStruct;
                % and are contained in the first input `obj'.
                %
                % The outputs are:
                %%
                %   (1) LHSExpansion: expansion of parameter independent matrices
                %   (2) RHSExpansion: expansion of parameter independent vectors
                
                import Core.AssemblerADRHandler
                
                % Definition of the Object from the AssemblerADRHandler class

                build_IGA = AssemblerADRHandler();

                % Properties Assignment, with the exception of coefficientForm and
                % dirCondFuncStruct

                build_IGA.dimModalBasis = obj.dimModalBasis(1);
                build_IGA.leftBDomain_inX = obj.domainLimit_inX(1);
                build_IGA.rightBDomain_inX = obj.domainLimit_inX(2);
                build_IGA.downBDomain_inY = obj.domainLimit_inY(1);
                build_IGA.upBDomain_inY = obj.domainLimit_inY(2);
                build_IGA.stepMeshX = obj.stepMeshX(1);
                build_IGA.label_upBoundDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_UP_TAG;
                build_IGA.label_downBoundDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_DOWN_TAG;
                build_IGA.localdata_upBDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_UP_DATA;
                build_IGA.localdata_downBDomain = obj.dirCondFuncStructExpansion.igaBoundCond.BC_DOWN_DATA;
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
                build_IGA.igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                
                % Prepare storage
                LHSExpansion = {};
                RHSExpansion = {};
                
                % Preassemble all terms related to mu
                for q = 1:length(obj.coefficientFormExpansion.mu)
                    mu = obj.coefficientFormExpansion.mu{q}{2};
                    beta1 = @(x, y) 0;
                    beta2 = @(x, y) 0;
                    sigma = @(x, y) 0;
                    coeffrobin = obj.coefficientFormExpansion.coeffrobin;
                    build_IGA.coefficientForm = struct('mu', mu, 'beta1', beta1, 'beta2', beta2, 'sigma', sigma, 'coeffrobin', coeffrobin);
                    igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                    force = @(x,y) 0;
                    build_IGA.dirCondFuncStruct = struct('igaBoundCond', igaBoundCond, 'force', force);
                    [A_q,~,~,~,~,~,~,~] = buildSystemIGAScatter(build_IGA);
                    LHSExpansion{length(LHSExpansion) + 1} = {...
                        obj.coefficientFormExpansion.mu{q}{1}, A_q ...
                    };
                end
                
                % Preassemble all terms related to beta1
                for q = 1:length(obj.coefficientFormExpansion.beta1)
                    mu = @(x, y) 0;
                    beta1 = obj.coefficientFormExpansion.beta1{q}{2};
                    beta2 = @(x, y) 0;
                    sigma = @(x, y) 0;
                    coeffrobin = obj.coefficientFormExpansion.coeffrobin;
                    build_IGA.coefficientForm = struct('mu', mu, 'beta1', beta1, 'beta2', beta2, 'sigma', sigma, 'coeffrobin', coeffrobin);
                    igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                    force = @(x,y) 0;
                    build_IGA.dirCondFuncStruct = struct('igaBoundCond', igaBoundCond, 'force', force);
                    [A_q,~,~,~,~,~,~,~] = buildSystemIGAScatter(build_IGA);
                    LHSExpansion{length(LHSExpansion) + 1} = {...
                        obj.coefficientFormExpansion.beta1{q}{1}, A_q ...
                    };
                end
                
                % Preassemble all terms related to beta2
                for q = 1:length(obj.coefficientFormExpansion.beta2)
                    mu = @(x, y) 0;
                    beta1 = @(x, y) 0;
                    beta2 = obj.coefficientFormExpansion.beta2{q}{2};
                    sigma = @(x, y) 0;
                    coeffrobin = obj.coefficientFormExpansion.coeffrobin;
                    build_IGA.coefficientForm = struct('mu', mu, 'beta1', beta1, 'beta2', beta2, 'sigma', sigma, 'coeffrobin', coeffrobin);
                    igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                    force = @(x,y) 0;
                    build_IGA.dirCondFuncStruct = struct('igaBoundCond', igaBoundCond, 'force', force);
                    [A_q,~,~,~,~,~,~,~] = buildSystemIGAScatter(build_IGA);
                    LHSExpansion{length(LHSExpansion) + 1} = {...
                        obj.coefficientFormExpansion.beta2{q}{1}, A_q ...
                    };
                end
                
                % Preassemble all terms related to sigma
                for q = 1:length(obj.coefficientFormExpansion.sigma)
                    mu = @(x, y) 0;
                    beta1 = @(x, y) 0;
                    beta2 = @(x, y) 0;
                    sigma = obj.coefficientFormExpansion.sigma{q}{2};
                    coeffrobin = obj.coefficientFormExpansion.coeffrobin;
                    build_IGA.coefficientForm = struct('mu', mu, 'beta1', beta1, 'beta2', beta2, 'sigma', sigma, 'coeffrobin', coeffrobin);
                    igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                    force = @(x,y) 0;
                    build_IGA.dirCondFuncStruct = struct('igaBoundCond', igaBoundCond, 'force', force);
                    [A_q,~,~,~,~,~,~,~] = buildSystemIGAScatter(build_IGA);
                    LHSExpansion{length(LHSExpansion) + 1} = {...
                        obj.coefficientFormExpansion.sigma{q}{1}, A_q ...
                    };
                end
                
                % Preassemble all terms related to force
                for q = 1:length(obj.dirCondFuncStructExpansion.force)
                    mu = @(x, y) 0;
                    beta1 = @(x, y) 0;
                    beta2 = @(x, y) 0;
                    sigma = @(x, y) 0;
                    coeffrobin = obj.coefficientFormExpansion.coeffrobin;
                    build_IGA.coefficientForm = struct('mu', mu, 'beta1', beta1, 'beta2', beta2, 'sigma', sigma, 'coeffrobin', coeffrobin);
                    igaBoundCond = obj.dirCondFuncStructExpansion.igaBoundCond;
                    force = obj.dirCondFuncStructExpansion.force{q}{2};
                    build_IGA.dirCondFuncStruct = struct('igaBoundCond', igaBoundCond, 'force', force);
                    [~,f_q,~,~,~,~,~,~] = buildSystemIGAScatter(build_IGA);
                    RHSExpansion{length(RHSExpansion) + 1} = {...
                        obj.dirCondFuncStructExpansion.force{q}{1}, f_q ...
                    };
                end
            end % End function
            
            %% Method 'reducedSolve'.
            
            function [u] = reducedSolve(obj, param, n)

                %%
                % reducedSolve
                % Perform an HiPOD reduced solve.
                %
                % Class properties employed in this method are:
                %   (1) lhs_expansion
                %   (2) rhs_expansion
                % and are contained in the first input `obj'.
                % Furthermore, the second input `param' contains the current value of the
                % parameters, and is a vector of real numbers.
                % The third input `n' contains the basis size, and should be less than
                % ore equal to obj.N
                %
                % The outputs are:
                %%
                %   (1) u                : HiPOD reduced solution of the Elliptic Problem
                
                disp(['Reduced solve for param=', num2str(param), ', n=', num2str(n)])
                A = 0;
                for q = 1:length(obj.reducedLHSExpansion)
                    A = A + obj.reducedLHSExpansion{q}{1}(param)*obj.reducedLHSExpansion{q}{2}(1:n, 1:n);
                end
                f = 0;
                for q = 1:length(obj.reducedRHSExpansion)
                    f = f + obj.reducedRHSExpansion{q}{1}(param)*obj.reducedRHSExpansion{q}{2}(1:n);
                end
                
                u = A\f;
                
            end % End function
                
    end
end

%% Auxiliary function to sum function handles appearing
% in the affine expansion
function [result] = sumHandles(expansion, param)
    result = @(x, y) 0;
    for q = 1:length(expansion)
        result = @(x, y) result(x, y) + expansion{q}{1}(param)*expansion{q}{2}(x, y);
    end
end

