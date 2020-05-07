classdef AssemblerADRHandler
    
    
    %% ASSEMBLER ADR HANDLER CLASS
    % the AssemblerADRHandler is a class that contain all the scripts
    % responsible for the assembling and building of the block matrix that
    % describe the discretized differential problem. On one hand, the class
    % properties, the properties of the objects were defined based on the
    % variables needed in the implementation of each function. On the other 
    % hand, all of the previous functions are concentrated and organized in
    % the class methods.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% ASSEMBLER ADR HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the AssemblerADRhandler
        % encapsulate all of the variables needed to run the methods
        % bellow, from the support methods to the final building methods.
        
        % BUILDING PROPERTIES
                
        dimModalBasis;              % Dimension of the Modal Basis
        
        leftBDomain_inX;            % Left Limit of the Domain in the X Direction
        
        rightBDomain_inX;           % Right Limit of the Domain in the X Direction
        
        downBDomain_inY;
        
        upBDomain_inY;
        
        stepMeshX;                  % Vector Containing the Step of the Finite
                                    % Element Mesh
                      
        label_upBoundDomain;        % Contains the Label Identifying the Nature of
                                    % the Boundary Conditions on the Upper Limit of
                                    % the Domain
                      
        label_downBoundDomain;      % Contains the Label Identifying the Nature of
                                    % the Boundary Conditions on the Lower Limit of
                                    % the Domain
                      
        localdata_upBDomain;        % Contains the Values of the Boundary Conditions
                                    % on the Upper Limir of the Domain
                      
        localdata_downBDomain;      % Contains the Values of the Boundary Conditions
                                    % on the Lower Limir of the Domain
                      
        domainPosition;             % Current domain used in the domain decomposition
                                    % Computations
        
        coefficientForm;            % Data Structure Containing All the @-Functions
                                    % and the Constants Relative to the Bilinear Form
                            
        dirCondFuncStruct;          % Data Structure Containing All the @-Functions
                                    % for the Dirichlet Conditions at the Inflow and
                                    % for the Exciting Force
                      
        geometricInfo;              % Data Structure Containing All the
                                    % Geometric Information regarding the
                                    % Domain. The current version of the code
                                    % works only for the specific condition of:
                                    % (L = 1, a = 0, psi_x = 0)
    
        robinCondStruct;            % Data Structure Containing the Two Values of the
                                    % Coefficients (R, L) for the Robin Condition Used
                                    % in the Domain Decomposition
                      
        numbDimMBEachDomain;        % Number of Elements in the Vector Containing the
                                    % Dimensions of the Modal Basis in Each Domain
                      
        couplingCond_DD;            % Contains the Label Adressing the Coupling Condition
                                    % of the Problem        
        
        physicMesh_inX;             % Vector Containing the Physical Mesh in the X
                                    % Direction
                      
        physicMesh_inY;             % Vector Containing the Physical Mesh in the Y
                                    % Direction
                      
        jacAtQuadNodes;             % Data Sturcture Used to Save the Value of the
                                    % Jacobians Computed in the Quadratures Nodes 
                      
        degreePolySplineBasis;      % Degree of the Polynomial B-Spline Basis
        
        continuityParameter;        % Degree of Continuity of the Basis 'C^(p-k)'
        
        domainProfile;              % Symbolic Function Defining the Profile of the
                                    % Simulation Domain
                      
        domainProfileDer;           % Symbolic Function Defining the Derivative of
                                    % the Profile of the Simulation Domain   
                                    
        timeInstant                 % Current time instant of the simulation
        
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
        
        delta;              % Parameter to tune for stabailization
        
        centreline;         % NURBS curve defining the centreline of 
                            % the physical domain
                            
        stabMethod;         % Specification of the stabilization method to
                            % be used in the solution
                            
        stabDelta;          % Scalling factor that may be used to control the
                            % amount of stabilization to be added in the
                            % problem
                            
        PeConfig            % Choice of the definition of the Pechlet number 
                            % to be used
                            
        timeDomain
        
        igaBoundCond
        
        currTimeSol
        
        assemblerStruct

    end
    
    methods (Access = public)
        
        %% ASSEMBLER ADR HANDLER - CONSTRUCT METHOD
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
        
        %% ASSEMBLER ADR HANDLER - BUILDING METHODS
        % The following functions correspond to the individual scripts
        % previously developed of the functions 'buildSystem',
        % 'buildSystemIGA' and 'buildSystemP2'.
        
            %% Method 'buildSystemFEM'
            
                function [A,b,modalBasis,liftCoeffA,liftCoeffB,verGLNodes,verGLWeights] = buildSystemFEM(obj)

                %%
                % buildSystemIGA - This function computes the assembled matrices
                %                    relative to the variational problem considering 
                %                    IGA modal basis.
                %
                % Note:
                % All of the following inputs are encapsulated in the
                % object properties.
                %
                % The inputs are:
                %%
                %   (1)  dimModalBasis              : Dimension of the Modal Basis
                %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
                %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
                %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
                %                                     Element Mesh
                %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
                %                                     the Boundary Conditions on the Upper Limit of
                %                                     the Domain
                %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
                %                                     the Boundary Conditions on the Lower Limit of
                %                                     the Domain
                %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
                %                                     on the Upper Limir of the Domain
                %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
                %                                     on the Lower Limir of the Domain
                %   (9) domainPosition              : Current domain used in the domain decomposition
                %                                     Computations
                %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
                %                                     and the Constants Relative to the Bilinear Form
                %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
                %                                     for the Dirichlet Conditions at the Inflow and
                %                                     for the Exciting Force
                %   (12) geometricInfo              : Data Structure Containing All the
                %                                     Geometric Information regarding the
                %                                     Domain. The current version of the code
                %                                     works only for the specific condition of:
                %                                     (L = 1, a = 0, psi_x = 0)
                %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
                %                                     Coefficients (R, L) for the Robin Condition Used
                %                                     in the Domain Decomposition
                %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
                %                                     Dimensions of the Modal Basis in Each Domain
                %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
                %                                     of the Problem
                %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
                %                                     Direction
                %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
                %                                     Direction
                %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
                %                                     Jacobians Computed in the Quadratures Nodes 
                %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
                %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
                %   (21) domainProfile              : Symbolic Function Defining the Profile of the
                %                                     Simulation Domain
                %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
                %                                     the Profile of the Simulation Domain
                % 
                % The outputs are:
                %%
                %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
                %   (2) b                   : Final Assembled Block Vector Using IGA Basis
                %   (3) modalBasis          : Modal Basis Obtained
                %   (3) liftCoeffA          : First Offset Adjustment Coefficient
                %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
                %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
                %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes

                %% IMPORT CLASSES

                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %% NUMBER OF NODES IN THE QUADRATURE FORMULA
                %-------------------------------------------------------------%
                % The values of 'nqnx' and 'nqny' does not change during the
                % computation of the Gauss-Legendre Integration Points.
                % Attention, if you increase the number of nodes in the
                % discretization, it is necessary to refine the quadrature
                % rule.
                %-------------------------------------------------------------%

                % Horizontal direction

                numbHorNodes = obj.numbHorQuadNodes;

                % Vertical Direction

                numbVerNodes = obj.numbVerQuadNodes;

                %% NUMBER OF ELEMENTS - FINITE ELEMENTS
                %-------------------------------------------------------------%
                % Total number of elements used in the finite element mesh.
                % Also called total number of intervals along the code.
                %-------------------------------------------------------------%

                numbElements  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);

                %% NUMBER OF NODES - FINITE ELEMENTS
                %-------------------------------------------------------------%
                % Total number of nodes used in the finite element mesh.
                %-------------------------------------------------------------%

                numbNodes  = numbElements+1;   

                %% GAUSS-LEGENDRE INTEGRATION NODES
                %-------------------------------------------------------------%
                % The following method reveives the number of integration nodes
                % in the horizontal and vertical direction and retruns the
                % respective Gauss-Legendre nodes and weigths for the
                % integration interval [0,1].
                %
                % Mesh FEM in X: Equispaced Nodes                            
                %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
                %                    Direction
                %  (2) horGLWeights : Vector of the Standard Weights on the X 
                %                     Direction
                %
                % Mesh FEM in Y: Equispaced Nodes                              
                %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
                %                    Direction  
                %  (2) verGLWeights : Vector of the Standard Weights on the Y 
                %                     Direction
                %-------------------------------------------------------------%

                % Horizontal direction

                obj_gaussLegendre_1 = IntegrateHandler();
                obj_gaussLegendre_1.numbQuadNodes = numbHorNodes;
                [~, horGLNodes, horGLWeights] = gaussLegendre(obj_gaussLegendre_1); 

                % Vertical direction

                obj_gaussLegendre_2 = IntegrateHandler();
                obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
                [~, verGLNodes, ~] = gaussLegendre(obj_gaussLegendre_2);             

                %% FEM MESH IN THE X DIRECTION
                %-------------------------------------------------------------%
                % Creation of the finite element mesh in the X direction
                % considering equispaced nodes. The mesh is created using the
                % total number of nodes and the left limit of the domain. Both
                % information come from the demo file and are passed here as a
                % property of the object.
                %-------------------------------------------------------------%

                meshFEM     = zeros(numbNodes,1);
                meshFEM(1)  = obj.leftBDomain_inX;

                for i=2:numbNodes
                    meshFEM(i) = meshFEM(i-1)+obj.stepMeshX;
                end

                %% AUGMENTED FEM MESH + WEIGHTS 
                %-------------------------------------------------------------%
                % Augmented finite element mesh containing all of the 
                % discretization nodes and the quadrature nodes contained in 
                % each element.
                %-------------------------------------------------------------%

                augMeshFEM       = zeros( numbElements*numbHorNodes, 1);

                %-------------------------------------------------------------%
                % Note: The loop allocates the correponding quadrature nodes
                % and weights corresponding to each element of the
                % discretization.
                %-------------------------------------------------------------%

                for i = 1:numbElements

                    % STEP 1
                    %---------------------------------------------------------%
                    % In the first step, the Gauss-Legendre nodes computed
                    % previously are rescaled to fit the interval corresponding
                    % to the current element.
                    %---------------------------------------------------------%                

                    obj_quadratureRule = IntegrateHandler();

                    obj_quadratureRule.leftBoundInterval = meshFEM(i);
                    obj_quadratureRule.rightBoundInterval = meshFEM(i+1);
                    obj_quadratureRule.inputNodes = horGLNodes;
                    obj_quadratureRule.inputWeights = horGLWeights;

                    [augMeshFEM((i-1)*numbHorNodes+1 : i*numbHorNodes), ...
                     augMeshFEMWeights((i-1)*numbHorNodes+1 : i*numbHorNodes), ] = ...
                                                     quadratureRule(obj_quadratureRule);

                    % STEP 2
                    %---------------------------------------------------------%
                    % In the second step, the nodes and weights are again
                    % rescaled to take in consideration the geometry of the
                    % domain, more specifically the profile of the centerline.
                    %---------------------------------------------------------%

                    auxMesh = augMeshFEM((i-1) * numbHorNodes + 1 : i * numbHorNodes);

                    for hp = 1: numbHorNodes

                        obj_gaussLegendre = IntegrateHandler();
                        obj_gaussLegendre.numbQuadNodes = 16;

                        [~,auxPoints,auxWeights] = gaussLegendre(obj_gaussLegendre);

                        auxPoints = auxMesh(hp) * auxPoints;
                        auxWeights = auxMesh(hp) * auxWeights;
                        auxMesh(hp) = sum(sqrt(1+(obj.domainProfileDer(auxPoints)).^2).*auxWeights); 

                    end
                    augMeshFEM((i-1) * numbHorNodes + 1 : i * numbHorNodes) = auxMesh;
                end

                %% MODAL BASIS IN THE Y DIRECTION
                %-------------------------------------------------------------%
                % The next method takes the coefficients of the bilinear form,
                % the dimension of the modal basis (assigned in the demo) and
                % the vertical evaluation points and computes the coefficients
                % of the modal basis when evaluated in the points of the
                % transverse fiber.
                %-------------------------------------------------------------%

                obj_newModalBasis = BasisHandler();

                obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
                obj_newModalBasis.evalNodesY = verGLNodes;
                obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
                obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
                obj_newModalBasis.coeffForm = obj.coefficientForm;

                [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);
                
                %% FEM BASIS IN THE X DIRECTION
                
                obj_newFiniteElementBasis = BasisHandler();
                
                obj_newFiniteElementBasis.degreeFiniteElement = 1;
                obj_newFiniteElementBasis.meshNodesX = meshFEM;
                obj_newFiniteElementBasis.meshQuadratureNodesX = augMeshFEM;

                [basisFEM, basisFEMDer] = newFiniteElementBasis(obj_newFiniteElementBasis);
                
                %% MEMORY ALLOCATION FOR SYSTEM MATRICES

                A = sparse( numbNodes*(obj.dimModalBasis), numbNodes*(obj.dimModalBasis) );
                b = zeros ( numbNodes*(obj.dimModalBasis), 1);

                %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
                %-------------------------------------------------------------%
                % We use the horizontal mesh created previously to evaluate the
                % thickness and the vertical coordinate of the centerline of
                % the channel.
                %-------------------------------------------------------------%

                [horEvalNodes,verEvalNodes]   = meshgrid( augMeshFEM, verGLNodes-0.5 );%!!!!!!

                % THICKNESS

                evalL     = obj.geometricInfo.L(horEvalNodes)';

                %% Y COORDINATE OF THE CENTERLINE
                %-------------------------------------------------------------%
                % In the current code this value is not used to perform any
                % computation.
                %-------------------------------------------------------------%

                % a_computed     = obj.geometricInfo.a(horEvalNodes)';        

                %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
                %-------------------------------------------------------------%
                % We use the horizontal and vertical meshes to evaluate the
                % coefficients of the bilinear form of the original equation in
                % the entire domain.
                %-------------------------------------------------------------%

                % DIFFUSION

                evalMu    = obj.coefficientForm.mu(horEvalNodes,verEvalNodes)';

                % ADVECTION

                evalBeta1 = obj.coefficientForm.beta1(horEvalNodes,verEvalNodes)';
                evalBeta2 = obj.coefficientForm.beta2(horEvalNodes,verEvalNodes)';

                % REACTION

                evalSigma = obj.coefficientForm.sigma(horEvalNodes,verEvalNodes)';

                %% EVALUATION OF THE EXCITING FORCE 
                %-------------------------------------------------------------%
                % Finally, we use the horizontal and vertical meshes to
                % evaluate the exciting force acting on the system in the whole
                % domain.
                %-------------------------------------------------------------%

                evalForce = obj.dirCondFuncStruct.force(horEvalNodes,verEvalNodes)';

                %-------------------------------------------------------------%
                % Note: All of the coefficients of the bilinear form and the
                % exciting term evaluated in the domain, as well as the mesh of
                % vertical coordinates are stored in a data structure called
                % "Computed" to be treated in the assembler function.
                %-------------------------------------------------------------%

                Computed = struct('mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2, ...
                                  'sigma_c',evalSigma,'force_c',evalForce,'y',verEvalNodes);

                %% LIFTING
                %-------------------------------------------------------------%
                % Compute the lifting contribution in the force vector due to
                % the non homogeneous Dirichlet boundary conditions in the
                % lower and upper boundary.
                %-------------------------------------------------------------%

                obj_liftBoundCond = BoundaryConditionHandler();

                obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
                obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
                obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
                obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
                obj_liftBoundCond.coeffForm = obj.coefficientForm;

                [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
                liftFunc = @(x,y) aLift * y + bLift;
                lifting = liftFunc(0,Computed.y)';

                %% JACOBIAN
                %-------------------------------------------------------------%
                % This section computes the curved domain and extract the
                % Jacobian vector that maps the centerline in the physical
                % domain to the computational domain.
                %-------------------------------------------------------------%

                % Mesh IGA of the curved domain

                meshFEMCurved    = zeros(numbNodes,1);
                [verGLNodes,verGLWeights] = gauss(numbVerNodes);

                verGLWeights = (obj.rightBDomain_inX - obj.leftBDomain_inX) * (obj.stepMeshX) * (verGLWeights*0.5);

                for j=1:numbElements

                    scalingVec = zeros(1,numbVerNodes);
                    auxNodes = (obj.rightBDomain_inX - obj.leftBDomain_inX) * obj.stepMeshX * (verGLNodes * 0.5 + 0.5) ...
                               + obj.leftBDomain_inX + (j-1) * obj.stepMeshX * (obj.rightBDomain_inX - obj.leftBDomain_inX);

                    for i=1:numbVerNodes
                        scalingVec(i) = sqrt(1+(obj.domainProfileDer(auxNodes(i)))^2);
                    end

                    meshFEMCurved(j+1) = meshFEMCurved(j) + sum(scalingVec * verGLWeights);

                end

                % Jacobian loop

                Jac = [];

                for iel = 1:numbElements

                    leftBound = meshFEMCurved(iel);     % Left Extreme
                    rightBound = meshFEMCurved(iel+1);  % Right Extreme

                    point = (rightBound - leftBound)*(0.5*verGLNodes + 0.5) + leftBound;

                    for j = 1 : numbVerNodes
                        Jac = [Jac sqrt(1+(obj.domainProfileDer(point(j)))^2)];
                    end
                end

                %% ASSEMBLING LOOP
                %-------------------------------------------------------------%
                % The assembling loop creates each one of the submatrices
                % correponding to the microstructure of the linear system. The
                % loop is repeated m^2 times to complete the macrostructure of
                % the system.
                %-------------------------------------------------------------%

                for imb = 1:obj.dimModalBasis
                    for kmb = 1:obj.dimModalBasis

                        [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerFEM( imb, kmb, numbNodes, ...
                                                verGLWeights,modalBasis(:,imb), ...
                                                modalBasisDer(:,imb),modalBasis(:,kmb),...
                                                modalBasisDer(:,kmb), evalL, ...
                                                Computed,obj.stepMeshX,...
                                                obj.degreePolySplineBasis,obj.continuityParameter,...
                                                lifting,aLift,bLift,Jac,numbHorNodes,numbVerNodes,...
                                                numbElements,augMeshFEMWeights,basisFEM,basisFEMDer,obj.dimModalBasis);

                        % Assignment of the Block Matrix Just Assembled

                        A(1+(imb-1)*numbNodes : imb*numbNodes , 1+(kmb-1)*numbNodes : kmb*numbNodes) = Amb;

                    end

                    % Assignment of the Block Vector Just Assembled
                    b( 1+(imb-1)*numbNodes : imb*numbNodes ) = bmb;

                end

                end
        
            %% Method 'buildSystemIGA'
            
            function [A,b,modalBasis,liftCoeffA,liftCoeffB,Jac] = buildSystemIGA(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler

            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*obj.continuityParameter + obj.degreePolySplineBasis + 1 - obj.continuityParameter;
            
            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%

            % Horizontal direction
            
            obj_gaussLegendre_1 = IntegrateHandler();
            obj_gaussLegendre_1.numbQuadNodes = numbHorNodes;
            [~, horGLNodes, horGLWeights] = gaussLegendre(obj_gaussLegendre_1); 
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);        
            
            %% ISOGEOMETRIC MESH IN THE X DIRECTION
            %-------------------------------------------------------------%
            % Creation of the isogeometric mesh in the X direction
            % considering the control point assigned to the desired spline
            % curve. The mesh is created using the total number of control
            % point. The number of control points depend of the degree of
            % the spline basis, the continuity parameter 'k' and the number
            % of knots with which the physical mesh was divided.
            %-------------------------------------------------------------%
            
            stepKnot = (obj.rightBDomain_inX - obj.leftBDomain_inX)/(numbControlPts-1);
            meshIGA = zeros(1,numbControlPts);
            
            for i = 2:numbControlPts
                 meshIGA(i) = meshIGA(i-1) + stepKnot;
            end
            
            %% AUGMENTED ISOGEOMETRIC MESH + WEIGHTS 
            %-------------------------------------------------------------%
            % Creation of the finite element mesh in the X direction
            % considering equispaced nodes. The mesh is created using the
            % total number of nodes and the left limit of the domain. Both
            % information come from the demo file and are passed here as a
            % property of the object.
            %-------------------------------------------------------------%

            augMeshIGA = zeros( (numbControlPts-1)*numbHorNodes, 1);

            %-------------------------------------------------------------%
            % Note: The loop allocates the correponding quadrature nodes
            % and weights corresponding to each knot of the physical mesh
            % in the isogeometric analysis.
            %-------------------------------------------------------------%
             
            for i = 1:numbKnots
                
                % STEP 1
                %---------------------------------------------------------%
                % In the first step, the Gauss-Legendre nodes computed
                % previously are rescaled to fit the interval corresponding
                % to the current element.
                %---------------------------------------------------------%
                
                obj_quadratureRule = IntegrateHandler();
                
                obj_quadratureRule.leftBoundInterval = meshIGA(i);
                obj_quadratureRule.rightBoundInterval = meshIGA(i+1);
                obj_quadratureRule.inputNodes = horGLNodes;
                obj_quadratureRule.inputWeights = horGLWeights;
                
                [augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes), ~] = ...
                                                 quadratureRule(obj_quadratureRule);
             
                % STEP 2
                %---------------------------------------------------------%
                % In the second step, the nodes and weights are again
                % rescaled to take in consideration the geometry of the
                % domain, more specifically the profile of the centerline.
                %---------------------------------------------------------%

%                 auxMesh = augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes);
% 
%                 for hp = 1: numbHorNodes
%                     
%                     obj_gaussLegendre = IntegrateHandler();
%                     obj_gaussLegendre.numbQuadNodes = 16;
%                     
%                     [~,auxPoints,auxWeights] = gaussLegendre(obj_gaussLegendre);
%                     
%                     auxPoints = auxMesh(hp) * auxPoints;
%                     auxWeights = auxMesh(hp) * auxWeights;
%                     auxMesh(hp) = sum(sqrt(1+(obj.domainProfileDer(auxPoints)).^2).*auxWeights); 
%                                        
%                 end
%                 augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes) = auxMesh;
                                                    
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Debug
            % figure;
            % plot(meshIGA);
            % figure;
            % plot(augMeshIGA);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% AUGMENTED VERTICAL POINTS (RESCALE)
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = obj.downBDomain_inY;
            objVertQuadRule.rightBoundInterval = obj.upBDomain_inY;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            % DEBUG
            % plot(augVerNodes)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Debug
            % figure;
            % plot(augVerNodes);
            % figure;
            % plot(verGLNodes);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);
            
            % DEBUG
            
            %size(modalBasis)
            %size(modalBasisDer)
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%

            [horEvalNodes,verEvalNodes]   = meshgrid( augMeshIGA, augVerNodes );
            
            % THICKNESS
            
            evalL     = obj.geometricInfo.L(horEvalNodes)';
            
            %% Y COORDINATE OF THE CENTERLINE
            %-------------------------------------------------------------%
            % In the current code this value is not used to perform any
            % computation.
            %-------------------------------------------------------------%
            
            % a_computed     = obj.geometricInfo.a(horEvalNodes)';        
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal and vertical meshes to evaluate the
            % coefficients of the bilinear form of the original equation in
            % the entire domain.
            %-------------------------------------------------------------%

            % DIFFUSION
            
            evalMu    = obj.coefficientForm.mu(horEvalNodes,verEvalNodes)';
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(horEvalNodes,verEvalNodes)';
            evalBeta2 = obj.coefficientForm.beta2(horEvalNodes,verEvalNodes)';
            
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(horEvalNodes,verEvalNodes)';
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(horEvalNodes,verEvalNodes)';
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed = struct('mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2, ...
                              'sigma_c',evalSigma,'force_c',evalForce,'y',verEvalNodes);
                          
            % Debug
            %-------------------------------------------------------------%
            %[~,~] = contourf(horEvalNodes,verEvalNodes,evalForce',20);    
            %colormap(jet); title('Force in Assembler!')
            %-------------------------------------------------------------%
                                                    
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y)';
            
            %% JACOBIAN
            %-------------------------------------------------------------%
            % This section computes the curved domain and extract the
            % Jacobian vector that maps the centerline in the physical
            % domain to the computational domain.
            %-------------------------------------------------------------%
            
            Jac = [];
            refCoordinate = linspace(obj.leftBDomain_inX,obj.rightBDomain_inX,numbKnots+1);
            
            for ii = 1:numbKnots
                
                x_left = refCoordinate(ii);                         % Right Extreme
                x_right  = refCoordinate(ii+1);                       % Left Extreme
                
%                 [evalPoints,~] = gauss(numbHorNodes);
                evalPoints = horGLNodes;
                
                evalPoints = (x_right-x_left)*(0.5*evalPoints + 0.5) + x_left;
                
                for jj = 1 : numbHorNodes
                    Jac = [Jac sqrt(1+(obj.domainProfileDer(evalPoints(jj)))^2)];
                end

            end
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis
                for kmb = 1:obj.dimModalBasis

                    [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerIGA( imb, kmb, numbControlPts, ...
                                            augVerWeights,modalBasis(:,imb), ...
                                            modalBasisDer(:,imb),modalBasis(:,kmb),...
                                            modalBasisDer(:,kmb), evalL, ...
                                            Computed,obj.stepMeshX,...
                                            obj.degreePolySplineBasis,obj.continuityParameter,obj.leftBDomain_inX,...
                                            obj.rightBDomain_inX,lifting,aLift,bLift,Jac,numbHorNodes,numbVerNodes,...
                                            numbKnots,horGLWeights,obj.D1,obj.D2);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;

            end

            end
                        
            %% Method 'buildIGAScatter'
            
            function [AA,bb,modalBasis,liftCoeffA,liftCoeffB,space,refDomain1D,bcStruct] = buildSystemIGAScatter(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler
            
            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve and that the
            % formula used bellow corresponds to use a straight line as
            % reference supporting fiber.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + ...
                             1 + obj.continuityParameter;

            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);        
            
            %% EXTRACT GEOMETRIC INFORMATION
            
            geometry  = obj.geometricInfo.geometry;
            map       = @(x,y) obj.geometricInfo.map(x,y);
            Jac       = @(x,y) obj.geometricInfo.Jac(x,y);
            Hes       = @(x,y) obj.geometricInfo.Hes(x,y);
            
            %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
            % Note: Create the knots and control points of the reference
            % supporting fiber where the coupled 1D problems will be
            % solved.
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsub = numbKnots;
            degree = obj.degreePolySplineBasis;
            regularity = obj.continuityParameter;
            
            [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);
            
            %% GENERATE ISOGEOMETRIC MESH FUNCTION
            % Note: Use the number of quadrature nodes to generate the
            % quadrature rule using the gauss nodes. Then, use this
            % information with the computed knots and control points to
            % generate the isogeometric mesh of the centreline.
            
            rule     = msh_gauss_nodes (numbHorNodes);
            [qn, qw] = msh_set_quad_nodes (zeta, rule);
            msh      = msh_cartesian (zeta, qn, qw, refDomain1D);
            
            %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
            % Note: Construct the functional space used to discretize the
            % problem along the centreline and perform the isogeometric
            % analysis. 
            % Here, we also evaluate the functional spaces in the specific
            % element. This will be used in the assembling loop to compute
            % the operators using the basis functions.
            
            space    = sp_bspline (knots, degree, msh);
            
            numbKnots = msh.nel_dir;
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                shapFunc = sp_evaluate_col (space, msh_col);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFunc{iel,1} = shapFunc;
                spaceFunc{iel,2} = gradFunc;
                spaceFunc{iel,3} = msh_col;

            end
            
            %% AUGMENTED HORIZONTAL POINTS 
            % Note: The FOR loop runs over the knots in the centerline ans
            % compute the required quadrature nodes. The quadrature nodes
            % will be necessary to evaluate the bilinear coefficients and
            % the forcing component.

            horEvalNodes = zeros(numbKnots * numbHorNodes,1);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                horEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  
                
            end
            
            %% AUGMENTED VERTICAL POINTS
            % Note: Since we are using the full map from the physical
            % domain to the reference domain, then the whole modal analysis
            % has to be performed in the vertical direction in the interval
            % [0,1];
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = 0;
            objVertQuadRule.rightBoundInterval = 1;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            obj_newModalBasis = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            obj_newModalBasis.dimLegendreBase = obj.dimModalBasis;
            obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%
            
            verEvalNodes = augVerNodes;
            
            X = mapOutParallel(horEvalNodes,verEvalNodes,map,1);
            Y = mapOutParallel(horEvalNodes,verEvalNodes,map,2);
            
            geoData.X = X;
            geoData.Y = Y;
            geoData.horNodes = horEvalNodes;
            geoData.verNodes = verEvalNodes;
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We need to evaluate all the coefficients of the bilinear form
            % in the quadrature nodes along the vertical direction to
            % perform the first integral (along the transverse fiber) to
            % obtain a the coefficients of the 1D coupled problem. The
            % result will be coefficients as a function of 'x'.
            %-------------------------------------------------------------%
            

            % DIFFUSION

            evalMu    = obj.coefficientForm.mu(X,Y);
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(X,Y);
            evalBeta2 = obj.coefficientForm.beta2(X,Y);
            
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(X,Y);
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(X,Y);
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed = struct('mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2, ...
                              'sigma_c',evalSigma,'force_c',evalForce,'y',verEvalNodes);
                          
            % Debug
            %-------------------------------------------------------------%
            %[~,~] = contourf(horEvalNodes,verEvalNodes,evalForce',20);    
            %colormap(jet); title('Force in Assembler!')
            %-------------------------------------------------------------%
            
            %% EVALUATION OF THE GEOMETRY PROPERTIES
            % Note: Since there is a transformation from the physical domain to the
            % physical domain, we must compute the map contribution used in the
            % computation of the coefficients and the Jacobian contribution to be
            % inserted in the integral (quadrature formula).
            
            evalJac     = jacOutParallel(horEvalNodes,verEvalNodes,Jac);
            evalDetJac  = detJacParallel(horEvalNodes,verEvalNodes,Jac);
            
            Phi1_dx = invJacParallel(1,1,evalJac);
            Phi1_dy = invJacParallel(1,2,evalJac);
            Phi2_dx = invJacParallel(2,1,evalJac);
            Phi2_dy = invJacParallel(2,2,evalJac);
            
            jacFunc.evalJac     = evalJac;
            jacFunc.evalDetJac  = evalDetJac;
            jacFunc.Phi1_dx     = Phi1_dx;
            jacFunc.Phi1_dy     = Phi1_dy;
            jacFunc.Phi2_dx     = Phi2_dx;
            jacFunc.Phi2_dy     = Phi2_dy;
                                                    
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y);
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis
                for kmb = 1:obj.dimModalBasis

                    [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerIGAScatter( imb, kmb, ...
                                            augVerWeights,modalBasis(:,imb),modalBasisDer(:,imb),modalBasis(:,kmb),...
                                            modalBasisDer(:,kmb),geoData,Computed,lifting,aLift,bLift,msh,space,jacFunc,...
                                            spaceFunc,obj.dirCondFuncStruct.igaBoundCond);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;
                    
                    disp(['FINISHED ASSEMBLING LOOP (',num2str(imb),' , ',num2str(kmb),')']);

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;
            end
            
            %% IMPOSE BOUNDARY CONDITIONS
            
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            BC_l = obj.igaBoundCond.BC_INF_TAG;
            BC_r = obj.igaBoundCond.BC_OUT_TAG;
            infBoundCond = obj.igaBoundCond.BC_INF_DATA;
            outBoundCond = obj.igaBoundCond.BC_OUT_DATA;
            
            obj_bcCoeff = BoundaryConditionHandler();
    
            obj_bcCoeff.infBoundCond  = infBoundCond;
            obj_bcCoeff.outBoundCond  = outBoundCond;
            obj_bcCoeff.augVerNodes   = augVerNodes;
            obj_bcCoeff.augVerWeights = augVerWeights;
            obj_bcCoeff.modalBasis = modalBasis;
            obj_bcCoeff.dimModalBasis = obj.dimModalBasis;
            obj_bcCoeff.coefficientForm = obj.coefficientForm;

            [infStruct,outStruct] = computeFourierCoeff(obj_bcCoeff);
            
            bcStruct.bcInfTag = BC_l;
            bcStruct.bcOutTag = BC_r;
            bcStruct.infStruct = infStruct;
            bcStruct.outStruct = outStruct;
            bcStruct.numbControlPts = numbControlPts;
            bcStruct.dimModalBasis = obj.dimModalBasis;
            
            [AA,bb] = impose_boundary(obj.dimModalBasis,BC_l, infStruct, BC_r, outStruct,A,b,numbControlPts);

            end
            
            %% Method 'buildIGAScatterWeak'
            
            function [AA,bb,modalBasis,liftCoeffA,liftCoeffB,space,refDomain1D,bcStruct] = buildSystemIGAScatterWeak(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler
            
            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve and that the
            % formula used bellow corresponds to use a straight line as
            % reference supporting fiber.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + ...
                             1 + obj.continuityParameter;

            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);        
            
            %% EXTRACT GEOMETRIC INFORMATION
            
            geometry  = obj.geometricInfo.geometry;
            map       = @(x,y) obj.geometricInfo.map(x,y);
            Jac       = @(x,y) obj.geometricInfo.Jac(x,y);
            Hes       = @(x,y) obj.geometricInfo.Hes(x,y);
            
            %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
            % Note: Create the knots and control points of the reference
            % supporting fiber where the coupled 1D problems will be
            % solved.
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsub = numbKnots;
            degree = obj.degreePolySplineBasis;
            regularity = obj.continuityParameter;
            
            [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);
            
            %% GENERATE ISOGEOMETRIC MESH FUNCTION
            % Note: Use the number of quadrature nodes to generate the
            % quadrature rule using the gauss nodes. Then, use this
            % information with the computed knots and control points to
            % generate the isogeometric mesh of the centreline.
            
            rule     = msh_gauss_nodes (numbHorNodes);
            [qn, qw] = msh_set_quad_nodes (zeta, rule);
            msh      = msh_cartesian (zeta, qn, qw, refDomain1D);
            
            %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
            % Note: Construct the functional space used to discretize the
            % problem along the centreline and perform the isogeometric
            % analysis. 
            % Here, we also evaluate the functional spaces in the specific
            % element. This will be used in the assembling loop to compute
            % the operators using the basis functions.
            
            space    = sp_bspline (knots, degree, msh);
            
            numbKnots = msh.nel_dir;
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                shapFunc = sp_evaluate_col (space, msh_col);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFunc{iel,1} = shapFunc;
                spaceFunc{iel,2} = gradFunc;
                spaceFunc{iel,3} = msh_col;

            end
            
            %% AUGMENTED HORIZONTAL POINTS 
            % Note: The FOR loop runs over the knots in the centerline ans
            % compute the required quadrature nodes. The quadrature nodes
            % will be necessary to evaluate the bilinear coefficients and
            % the forcing component.

            horEvalNodes = zeros(numbKnots * numbHorNodes,1);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                horEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  
                
            end
            
            %% AUGMENTED VERTICAL POINTS
            % Note: Since we are using the full map from the physical
            % domain to the reference domain, then the whole modal analysis
            % has to be performed in the vertical direction in the interval
            % [0,1];
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = 0;
            objVertQuadRule.rightBoundInterval = 1;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            obj_newModalBasis = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            obj_newModalBasis.dimLegendreBase = obj.dimModalBasis;
            obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis, modalBasisDer] = newModalBasisLegendre(obj_newModalBasis);
                
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%
            
            verEvalNodes = augVerNodes;
            
            X = mapOut(horEvalNodes,verEvalNodes,map,1);
            Y = mapOut(horEvalNodes,verEvalNodes,map,2);
            
            geoData.X = X;
            geoData.Y = Y;
            geoData.horNodes = horEvalNodes;
            geoData.verNodes = verEvalNodes;
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We need to evaluate all the coefficients of the bilinear form
            % in the quadrature nodes along the vertical direction to
            % perform the first integral (along the transverse fiber) to
            % obtain a the coefficients of the 1D coupled problem. The
            % result will be coefficients as a function of 'x'.
            %-------------------------------------------------------------%
            

            % DIFFUSION

            evalMu    = obj.coefficientForm.mu(X,Y);
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(X,Y);
            evalBeta2 = obj.coefficientForm.beta2(X,Y);
            
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(X,Y);
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(X,Y);
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed = struct('mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2, ...
                              'sigma_c',evalSigma,'force_c',evalForce,'y',verEvalNodes);
                          
            % Debug
            %-------------------------------------------------------------%
            %[~,~] = contourf(horEvalNodes,verEvalNodes,evalForce',20);    
            %colormap(jet); title('Force in Assembler!')
            %-------------------------------------------------------------%
            
            %% EVALUATION OF THE GEOMETRY PROPERTIES
            % Note: Since there is a transformation from the physical domain to the
            % physical domain, we must compute the map contribution used in the
            % computation of the coefficients and the Jacobian contribution to be
            % inserted in the integral (quadrature formula).
            
            evalJac     = jacOut(horEvalNodes,verEvalNodes,Jac);
            evalDetJac  = detJac(horEvalNodes,verEvalNodes,Jac);
            
            Phi1_dx = invJac(1,1,evalJac);
            Phi1_dy = invJac(1,2,evalJac);
            Phi2_dx = invJac(2,1,evalJac);
            Phi2_dy = invJac(2,2,evalJac);
            
            jacFunc.evalJac     = evalJac;
            jacFunc.evalDetJac  = evalDetJac;
            jacFunc.Phi1_dx     = Phi1_dx;
            jacFunc.Phi1_dy     = Phi1_dy;
            jacFunc.Phi2_dx     = Phi2_dx;
            jacFunc.Phi2_dy     = Phi2_dy;
                                                    
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y);
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis
                for kmb = 1:obj.dimModalBasis

                    [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerIGAScatterWeak( imb, kmb, ...
                                            augVerWeights,modalBasis(:,imb),modalBasisDer(:,imb),modalBasis(:,kmb),...
                                            modalBasisDer(:,kmb),geoData,Computed,lifting,aLift,bLift,msh,space,jacFunc,...
                                            spaceFunc,obj.dirCondFuncStruct.igaBoundCond);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;
                    
                    disp(['FINISHED ASSEMBLING LOOP (',num2str(imb),' , ',num2str(kmb),')']);

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;
            end
            
            %% IMPOSE BOUNDARY CONDITIONS
            
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            BC_l = obj.igaBoundCond.BC_INF_TAG;
            BC_r = obj.igaBoundCond.BC_OUT_TAG;
            infBoundCond = obj.igaBoundCond.BC_INF_DATA;
            outBoundCond = obj.igaBoundCond.BC_OUT_DATA;
            
            obj_bcCoeff = BoundaryConditionHandler();
    
            obj_bcCoeff.infBoundCond  = infBoundCond;
            obj_bcCoeff.outBoundCond  = outBoundCond;
            obj_bcCoeff.augVerNodes   = augVerNodes;
            obj_bcCoeff.augVerWeights = augVerWeights;
            obj_bcCoeff.modalBasis = modalBasis;
            obj_bcCoeff.dimModalBasis = obj.dimModalBasis;
            obj_bcCoeff.coefficientForm = obj.coefficientForm;

            [infStruct,outStruct] = computeFourierCoeff(obj_bcCoeff);
            
            bcStruct.bcInfTag = BC_l;
            bcStruct.bcOutTag = BC_r;
            bcStruct.infStruct = infStruct;
            bcStruct.outStruct = outStruct;
            bcStruct.numbControlPts = numbControlPts;
            bcStruct.dimModalBasis = obj.dimModalBasis;
            
            [AA,bb] = impose_boundary(obj.dimModalBasis,BC_l, infStruct, BC_r, outStruct,A,b,numbControlPts);

            end
            
            %% Method 'buildIGAScatterTransient'
            
            function [AA,MM,bb,sTS,modalBasis,space,refDomain1D,bcStruct] = buildSystemIGAScatterTransient(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler
            
            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve and that the
            % formula used bellow corresponds to use a straight line as
            % reference supporting fiber.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + ...
                             1 + obj.continuityParameter;

            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);        
            
            %% EXTRACT GEOMETRIC INFORMATION
            
            geometry  = obj.geometricInfo.geometry;
            map       = @(x,y) obj.geometricInfo.map(x,y);
            Jac       = @(x,y) obj.geometricInfo.Jac(x,y);
            Hes       = @(x,y) obj.geometricInfo.Hes(x,y);
            
            %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
            % Note: Create the knots and control points of the reference
            % supporting fiber where the coupled 1D problems will be
            % solved.
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsub = numbKnots;
            degree = obj.degreePolySplineBasis;
            regularity = obj.continuityParameter;
            
            [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);
            
            %% GENERATE ISOGEOMETRIC MESH FUNCTION
            % Note: Use the number of quadrature nodes to generate the
            % quadrature rule using the gauss nodes. Then, use this
            % information with the computed knots and control points to
            % generate the isogeometric mesh of the centreline.
            
            rule     = msh_gauss_nodes (numbHorNodes);
            [qn, qw] = msh_set_quad_nodes (zeta, rule);
            msh      = msh_cartesian (zeta, qn, qw, refDomain1D);
            
            %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
            % Note: Construct the functional space used to discretize the
            % problem along the centreline and perform the isogeometric
            % analysis. 
            % Here, we also evaluate the functional spaces in the specific
            % element. This will be used in the assembling loop to compute
            % the operators using the basis functions.
            
            space    = sp_bspline (knots, degree, msh);
            
            numbKnots = msh.nel_dir;
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                shapFunc = sp_evaluate_col (space, msh_col);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFunc{iel,1} = shapFunc;
                spaceFunc{iel,2} = gradFunc;
                spaceFunc{iel,3} = msh_col;

            end
            
            %% AUGMENTED HORIZONTAL POINTS 
            % Note: The FOR loop runs over the knots in the centerline ans
            % compute the required quadrature nodes. The quadrature nodes
            % will be necessary to evaluate the bilinear coefficients and
            % the forcing component.

            horEvalNodes = zeros(numbKnots * numbHorNodes,1);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                horEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  
                
            end
            
            %% AUGMENTED VERTICAL POINTS
            % Note: Since we are using the full map from the physical
            % domain to the reference domain, then the whole modal analysis
            % has to be performed in the vertical direction in the interval
            % [0,1];
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = 0;
            objVertQuadRule.rightBoundInterval = 1;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            M = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%
            
            verEvalNodes = augVerNodes;
            
            X = mapOut(horEvalNodes,verEvalNodes,map,1);
            Y = mapOut(horEvalNodes,verEvalNodes,map,2);
            
            geoData.X = X;
            geoData.Y = Y;
            geoData.horNodes = horEvalNodes;
            geoData.verNodes = verEvalNodes;
            
            %% TIME PROPERTIES
            
            it = obj.timeInstant;
            tDomain = obj.timeDomain;
            
            time = tDomain(it);
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We need to evaluate all the coefficients of the bilinear form
            % in the quadrature nodes along the vertical direction to
            % perform the first integral (along the transverse fiber) to
            % obtain a the coefficients of the 1D coupled problem. The
            % result will be coefficients as a function of 'x'.
            %-------------------------------------------------------------%
            

            % DIFFUSION

            evalMu    = obj.coefficientForm.mu(X,Y,time);
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(X,Y,time);
            evalBeta2 = obj.coefficientForm.beta2(X,Y,time);
            
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(X,Y,time);
            
            % MASS
            
            evalMass = ones(size(X));
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(X,Y,time);
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed = struct('mass_c',evalMass,'mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2, ...
                              'sigma_c',evalSigma,'force_c',evalForce,'y',verEvalNodes);
                          
            % Debug
            %-------------------------------------------------------------%
            %[~,~] = contourf(horEvalNodes,verEvalNodes,evalForce',20);    
            %colormap(jet); title('Force in Assembler!')
            %-------------------------------------------------------------%
            
            %% EVALUATION OF THE GEOMETRY PROPERTIES
            % Note: Since there is a transformation from the physical domain to the
            % physical domain, we must compute the map contribution used in the
            % computation of the coefficients and the Jacobian contribution to be
            % inserted in the integral (quadrature formula).
            
            evalJac     = jacOut(horEvalNodes,verEvalNodes,Jac);
            evalDetJac  = detJac(horEvalNodes,verEvalNodes,Jac);
            
            Phi1_dx = invJac(1,1,evalJac);
            Phi1_dy = invJac(1,2,evalJac);
            Phi2_dx = invJac(2,1,evalJac);
            Phi2_dy = invJac(2,2,evalJac);
            
            jacFunc.evalJac     = evalJac;
            jacFunc.evalDetJac  = evalDetJac;
            jacFunc.Phi1_dx     = Phi1_dx;
            jacFunc.Phi1_dy     = Phi1_dy;
            jacFunc.Phi2_dx     = Phi2_dx;
            jacFunc.Phi2_dy     = Phi2_dy;
                                                    
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y);
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis
                for kmb = 1:obj.dimModalBasis

                    [Amb,Mmb,bmb] = assemblerIGAScatterTransient( imb, kmb, ...
                                            augVerWeights,modalBasis(:,imb),modalBasisDer(:,imb),modalBasis(:,kmb),...
                                            modalBasisDer(:,kmb),geoData,Computed,msh,space,jacFunc,spaceFunc);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;
                    M(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Mmb;
                    
                    disp(['FINISHED ASSEMBLING LOOP (',num2str(imb),' , ',num2str(kmb),')']);

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;
            end
            
            %% IMPOSE BOUNDARY CONDITIONS
            
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            BC_l = obj.igaBoundCond.BC_INF_TAG;
            BC_r = obj.igaBoundCond.BC_OUT_TAG;
            infBoundCond = obj.igaBoundCond.BC_INF_DATA;
            outBoundCond = obj.igaBoundCond.BC_OUT_DATA;
            
            obj_bcCoeff = BoundaryConditionHandler();
    
            obj_bcCoeff.infBoundCond  = infBoundCond;
            obj_bcCoeff.outBoundCond  = outBoundCond;
            obj_bcCoeff.augVerNodes   = augVerNodes;
            obj_bcCoeff.augVerWeights = augVerWeights;
            obj_bcCoeff.modalBasis = modalBasis;
            obj_bcCoeff.dimModalBasis = obj.dimModalBasis;
            obj_bcCoeff.coefficientForm = obj.coefficientForm;

            [infStruct,outStruct] = computeFourierCoeff(obj_bcCoeff);
            
            bcStruct.bcInfTag = BC_l;
            bcStruct.bcOutTag = BC_r;
            bcStruct.infStruct = infStruct;
            bcStruct.outStruct = outStruct;
            bcStruct.numbControlPts = numbControlPts;
            bcStruct.dimModalBasis = obj.dimModalBasis;
            
            [AA,MM,bb,sTS] = impose_boundary_transient(obj.dimModalBasis,BC_l, infStruct, BC_r, outStruct,A,M,b,numbControlPts,obj.currTimeSol);

            end
            
            %% Method 'buildTransientStruct'
            
            function [assemblerStruct] = buildTransientStruct(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
                %% IMPORT CLASSES

                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %% NUMBER OF NODES IN THE QUADRATURE FORMULA
                %-------------------------------------------------------------%
                % The values of 'nqnx' and 'nqny' does not change during the
                % computation of the Gauss-Legendre Integration Points.
                % Attention, if you increase the number of nodes in the
                % discretization, it is necessary to refine the quadrature
                % rule.
                %-------------------------------------------------------------%

                % Horizontal direction

                numbHorNodes = obj.numbHorQuadNodes;
                assemblerStruct.numbHorNodes = numbHorNodes;

                % Vertical Direction

                numbVerNodes = obj.numbVerQuadNodes;
                assemblerStruct.numbVerNodes = numbVerNodes;

                %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
                %-------------------------------------------------------------%
                % Total number of knots used to divide the spline curve in the
                % isogeometric analysis.
                %-------------------------------------------------------------%

                numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
                assemblerStruct.numbKnots = numbKnots;

                %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
                %-------------------------------------------------------------%
                % Total number of control points used in the isogeometric
                % approach. Note that this number is equal to the dimension of
                % the basis used to define the Spline curve and that the
                % formula used bellow corresponds to use a straight line as
                % reference supporting fiber.
                %-------------------------------------------------------------%

                numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + 1 + obj.continuityParameter;
                assemblerStruct.numbControlPts = numbControlPts;

                %% GAUSS-LEGENDRE INTEGRATION NODES
                %-------------------------------------------------------------%
                % The following method reveives the number of integration nodes
                % in the horizontal and vertical direction and retruns the
                % respective Gauss-Legendre nodes and weigths for the
                % integration interval [0,1].   
                %
                % Mesh FEM in X: Equispaced Nodes                            
                %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
                %                    Direction
                %  (2) horGLWeights : Vector of the Standard Weights on the X 
                %                     Direction
                %
                % Mesh FEM in Y: Equispaced Nodes                              
                %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
                %                    Direction  
                %  (2) verGLWeights : Vector of the Standard Weights on the Y 
                %                     Direction
                %-------------------------------------------------------------%

                % Vertical direction

                obj_gaussLegendre_2 = IntegrateHandler();
                obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
                [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);   

                assemblerStruct.verGLNodes = verGLNodes;
                assemblerStruct.verWeights = verWeights;

                %% EXTRACT GEOMETRIC INFORMATION

                geometry  = obj.geometricInfo.geometry;
                map       = @(x,y) obj.geometricInfo.map(x,y);
                Jac       = @(x,y) obj.geometricInfo.Jac(x,y);
                Hes       = @(x,y) obj.geometricInfo.Hes(x,y);

                assemblerStruct.geometry = geometry;
                assemblerStruct.map      = map;
                assemblerStruct.Jac      = Jac;
                assemblerStruct.Hes      = Hes;

                %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
                % Note: Create the knots and control points of the reference
                % supporting fiber where the coupled 1D problems will be
                % solved.

                refDomain1D = geo_load(nrbline ([0 0], [1 0]));

                nsub = numbKnots;
                degree = obj.degreePolySplineBasis;
                regularity = obj.continuityParameter;

                [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);

                assemblerStruct.refDomain1D = refDomain1D;
                assemblerStruct.nsub = nsub;
                assemblerStruct.degree = degree;
                assemblerStruct.regularity = regularity;
                assemblerStruct.knots = knots;
                assemblerStruct.zeta = zeta;

                %% GENERATE ISOGEOMETRIC MESH FUNCTION
                % Note: Use the number of quadrature nodes to generate the
                % quadrature rule using the gauss nodes. Then, use this
                % information with the computed knots and control points to
                % generate the isogeometric mesh of the centreline.

                rule     = msh_gauss_nodes (numbHorNodes);
                [qn, qw] = msh_set_quad_nodes (zeta, rule);
                msh      = msh_cartesian (zeta, qn, qw, refDomain1D);

                assemblerStruct.rule = rule;
                assemblerStruct.qn   = qn;
                assemblerStruct.qw   = qw;
                assemblerStruct.msh  = msh;

                %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
                % Note: Construct the functional space used to discretize the
                % problem along the centreline and perform the isogeometric
                % analysis. 
                % Here, we also evaluate the functional spaces in the specific
                % element. This will be used in the assembling loop to compute
                % the operators using the basis functions.

                space    = sp_bspline (knots, degree, msh);

                numbKnots = msh.nel_dir;

                for iel = 1:numbKnots

                    msh_col = msh_evaluate_col (msh, iel);

                    %---------------------------------------------------------%
                    % Evaluated space to compute the operators
                    %---------------------------------------------------------%

                    gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                    shapFunc = sp_evaluate_col (space, msh_col);

                    %---------------------------------------------------------%
                    % Store the spaces in the cell matrix
                    %---------------------------------------------------------%

                    spaceFunc{iel,1} = shapFunc;
                    spaceFunc{iel,2} = gradFunc;
                    spaceFunc{iel,3} = msh_col;

                end

                assemblerStruct.space     = space;
                assemblerStruct.numbKnots = numbKnots;
                assemblerStruct.spaceFunc = spaceFunc;

                %% AUGMENTED HORIZONTAL POINTS 
                % Note: The FOR loop runs over the knots in the centerline ans
                % compute the required quadrature nodes. The quadrature nodes
                % will be necessary to evaluate the bilinear coefficients and
                % the forcing component.

                horEvalNodes = zeros(numbKnots * numbHorNodes,1);

                for iel = 1:numbKnots

                    msh_col = msh_evaluate_col (msh, iel);

                    localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                    horEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  

                end

                assemblerStruct.horEvalNodes = horEvalNodes;

                %% AUGMENTED VERTICAL POINTS
                % Note: Since we are using the full map from the physical
                % domain to the reference domain, then the whole modal analysis
                % has to be performed in the vertical direction in the interval
                % [0,1];

                objVertQuadRule = IntegrateHandler();

                objVertQuadRule.leftBoundInterval = 0;
                objVertQuadRule.rightBoundInterval = 1;
                objVertQuadRule.inputNodes = verGLNodes;
                objVertQuadRule.inputWeights = verWeights;

                [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);

                assemblerStruct.augVerNodes = augVerNodes;
                assemblerStruct.augVerWeights = augVerWeights;

                %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
                %-------------------------------------------------------------%
                % The next method takes the coefficients of the bilinear form,
                % the dimension of the modal basis (assigned in the demo) and
                % the vertical evaluation points and computes the coefficients
                % of the modal basis when evaluated in the points of the
                % transverse fiber.
                %-------------------------------------------------------------%

                obj_newModalBasis = BasisHandler();

                obj_newModalBasis.dimModalBasis      = obj.dimModalBasis;
                obj_newModalBasis.evalNodesY         = verGLNodes;
                obj_newModalBasis.labelUpBoundCond   = obj.label_upBoundDomain;
                obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
                obj_newModalBasis.coeffForm          = obj.coefficientForm;

                [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);

                assemblerStruct.modalBasis          = modalBasis;
                assemblerStruct.modalBasisDer       = modalBasisDer;
                assemblerStruct.dimModalBasis       = obj.dimModalBasis;
                assemblerStruct.labelUpBoundCond    = obj.label_upBoundDomain;
                assemblerStruct.labelDownBoundCond  = obj.label_downBoundDomain;
                
                %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
                %-------------------------------------------------------------%
                % We use the horizontal mesh created previously to evaluate the
                % thickness and the vertical coordinate of the centerline of
                % the channel.
                %-------------------------------------------------------------%

                verEvalNodes = augVerNodes;

                X = mapOutParallel(horEvalNodes,verEvalNodes,map,1);
                Y = mapOutParallel(horEvalNodes,verEvalNodes,map,2);

                geoData.X = X;
                geoData.Y = Y;
                geoData.horNodes = horEvalNodes;
                geoData.verNodes = verEvalNodes;

                assemblerStruct.geoData = geoData;

                %% EVALUATION OF THE GEOMETRY PROPERTIES
                % Note: Since there is a transformation from the physical domain to the
                % physical domain, we must compute the map contribution used in the
                % computation of the coefficients and the Jacobian contribution to be
                % inserted in the integral (quadrature formula).

                evalJac     = jacOutParallel(horEvalNodes,verEvalNodes,Jac);
                evalDetJac  = detJacParallel(horEvalNodes,verEvalNodes,Jac);

                Phi1_dx = invJacParallel(1,1,evalJac);
                Phi1_dy = invJacParallel(1,2,evalJac);
                Phi2_dx = invJacParallel(2,1,evalJac);
                Phi2_dy = invJacParallel(2,2,evalJac);

                jacFunc.evalJac     = evalJac;
                jacFunc.evalDetJac  = evalDetJac;
                jacFunc.Phi1_dx     = Phi1_dx;
                jacFunc.Phi1_dy     = Phi1_dy;
                jacFunc.Phi2_dx     = Phi2_dx;
                jacFunc.Phi2_dy     = Phi2_dy;

                assemblerStruct.jacFunc = jacFunc;
                
                %% BOUNDARY CONDITIONS PROPERTIES
                
                assemblerStruct.igaBoundCond = obj.dirCondFuncStruct.igaBoundCond;
                
                %% MODEL PROPERTIES
                
                assemblerStruct.coefficientForm   = obj.coefficientForm;
                assemblerStruct.dirCondFuncStruct = obj.dirCondFuncStruct;
                
                %% GENERATION OF PLOTTING STRUCTURE
                
                assemblerStruct.plotResolutionX = 100;
                assemblerStruct.plotResolutionY = 100;

            end
            
            %% Method 'buildTransientSystem'
            
            function [AA,MM,bb,sTS,bcStruct] = buildTransientSystem(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
                %% IMPORT CLASSES

                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %% NUMBER OF NODES IN THE QUADRATURE FORMULA
                %-------------------------------------------------------------%
                % The values of 'nqnx' and 'nqny' does not change during the
                % computation of the Gauss-Legendre Integration Points.
                % Attention, if you increase the number of nodes in the
                % discretization, it is necessary to refine the quadrature
                % rule.
                %-------------------------------------------------------------%

                % Horizontal direction

                numbHorNodes = obj.assemblerStruct.numbHorNodes;

                % Vertical Direction

                numbVerNodes = obj.assemblerStruct.numbVerNodes;

                %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
                %-------------------------------------------------------------%
                % Total number of knots used to divide the spline curve in the
                % isogeometric analysis.
                %-------------------------------------------------------------%

                numbKnots  = obj.assemblerStruct.numbKnots;

                %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
                %-------------------------------------------------------------%
                % Total number of control points used in the isogeometric
                % approach. Note that this number is equal to the dimension of
                % the basis used to define the Spline curve and that the
                % formula used bellow corresponds to use a straight line as
                % reference supporting fiber.
                %-------------------------------------------------------------%

                numbControlPts = obj.assemblerStruct.numbControlPts;

                %% GAUSS-LEGENDRE INTEGRATION NODES
                %-------------------------------------------------------------%
                % The following method reveives the number of integration nodes
                % in the horizontal and vertical direction and retruns the
                % respective Gauss-Legendre nodes and weigths for the
                % integration interval [0,1].   
                %
                % Mesh FEM in X: Equispaced Nodes                            
                %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
                %                    Direction
                %  (2) horGLWeights : Vector of the Standard Weights on the X 
                %                     Direction
                %
                % Mesh FEM in Y: Equispaced Nodes                              
                %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
                %                    Direction  
                %  (2) verGLWeights : Vector of the Standard Weights on the Y 
                %                     Direction
                %-------------------------------------------------------------%

                verGLNodes = obj.assemblerStruct.verGLNodes;
                verWeights = obj.assemblerStruct.verWeights;

                %% EXTRACT GEOMETRIC INFORMATION

                geometry  = obj.assemblerStruct.geometry;
                map       = obj.assemblerStruct.map;
                Jac       = obj.assemblerStruct.Jac;
                Hes       = obj.assemblerStruct.Hes;

                %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
                % Note: Create the knots and control points of the reference
                % supporting fiber where the coupled 1D problems will be
                % solved.

                refDomain1D = obj.assemblerStruct.refDomain1D;
                nsub        = obj.assemblerStruct.nsub;
                degree      = obj.assemblerStruct.degree;
                regularity  = obj.assemblerStruct.regularity;
                knots       = obj.assemblerStruct.knots;
                zeta        = obj.assemblerStruct.zeta;

                %% GENERATE ISOGEOMETRIC MESH FUNCTION
                % Note: Use the number of quadrature nodes to generate the
                % quadrature rule using the gauss nodes. Then, use this
                % information with the computed knots and control points to
                % generate the isogeometric mesh of the centreline.

                rule = obj.assemblerStruct.rule;
                qn   = obj.assemblerStruct.qn;
                qw   = obj.assemblerStruct.qw;
                msh  = obj.assemblerStruct.msh;

                %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
                % Note: Construct the functional space used to discretize the
                % problem along the centreline and perform the isogeometric
                % analysis. 
                % Here, we also evaluate the functional spaces in the specific
                % element. This will be used in the assembling loop to compute
                % the operators using the basis functions.

                space     = obj.assemblerStruct.space;
                numbKnots = obj.assemblerStruct.numbKnots;
                spaceFunc = obj.assemblerStruct.spaceFunc;

                %% AUGMENTED HORIZONTAL POINTS 
                % Note: The FOR loop runs over the knots in the centerline ans
                % compute the required quadrature nodes. The quadrature nodes
                % will be necessary to evaluate the bilinear coefficients and
                % the forcing component.

                horEvalNodes = obj.assemblerStruct.horEvalNodes;

                %% AUGMENTED VERTICAL POINTS
                % Note: Since we are using the full map from the physical
                % domain to the reference domain, then the whole modal analysis
                % has to be performed in the vertical direction in the interval
                % [0,1];
                
                augVerNodes   = obj.assemblerStruct.augVerNodes;
                augVerWeights = obj.assemblerStruct.augVerWeights;

                %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
                %-------------------------------------------------------------%
                % The next method takes the coefficients of the bilinear form,
                % the dimension of the modal basis (assigned in the demo) and
                % the vertical evaluation points and computes the coefficients
                % of the modal basis when evaluated in the points of the
                % transverse fiber.
                %-------------------------------------------------------------%
                
                modalBasis    = obj.assemblerStruct.modalBasis;
                modalBasisDer = obj.assemblerStruct.modalBasisDer;
                dimModBasis   = obj.assemblerStruct.dimModalBasis;

                %% MEMORY ALLOCATION FOR SYSTEM MATRICES

                A = sparse( numbControlPts * dimModBasis, numbControlPts * dimModBasis );
                M = sparse( numbControlPts * dimModBasis, numbControlPts * dimModBasis );
                b = zeros ( numbControlPts * dimModBasis, 1);

                %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
                %-------------------------------------------------------------%
                % We use the horizontal mesh created previously to evaluate the
                % thickness and the vertical coordinate of the centerline of
                % the channel.
                %-------------------------------------------------------------%

                geoData = obj.assemblerStruct.geoData;
                X       = obj.assemblerStruct.geoData.X;
                Y       = obj.assemblerStruct.geoData.Y;
                
                %% EVALUATION OF THE GEOMETRY PROPERTIES
                % Note: Since there is a transformation from the physical domain to the
                % physical domain, we must compute the map contribution used in the
                % computation of the coefficients and the Jacobian contribution to be
                % inserted in the integral (quadrature formula).

                jacFunc = obj.assemblerStruct.jacFunc;
                
                %% EVALUATION OF TIME PROPERTIES
                
                time = obj.assemblerStruct.time;

                %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
                %-------------------------------------------------------------%
                % We need to evaluate all the coefficients of the bilinear form
                % in the quadrature nodes along the vertical direction to
                % perform the first integral (along the transverse fiber) to
                % obtain a the coefficients of the 1D coupled problem. The
                % result will be coefficients as a function of 'x'.
                %-------------------------------------------------------------%


                % DIFFUSION

                evalMu    = obj.assemblerStruct.coefficientForm.mu(X,Y,time);

                % ADVECTION

                evalBeta1 = obj.assemblerStruct.coefficientForm.beta1(X,Y,time);
                evalBeta2 = obj.assemblerStruct.coefficientForm.beta2(X,Y,time);

                % REACTION

                evalSigma = obj.assemblerStruct.coefficientForm.sigma(X,Y,time);

                % MASS

                evalMass = ones(size(X));

                %% EVALUATION OF THE EXCITING FORCE 
                %-------------------------------------------------------------%
                % Finally, we use the horizontal and vertical meshes to
                % evaluate the exciting force acting on the system in the whole
                % domain.
                %-------------------------------------------------------------%

                evalForce = obj.assemblerStruct.dirCondFuncStruct.force(X,Y,time);

                %-------------------------------------------------------------%
                % Note: All of the coefficients of the bilinear form and the
                % exciting term evaluated in the domain, as well as the mesh of
                % vertical coordinates are stored in a data structure called
                % "Computed" to be treated in the assembler function.
                %-------------------------------------------------------------%
                
                boundCond = obj.assemblerStruct.igaBoundCond;

                Computed = struct('mass_c',evalMass,'mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2,'sigma_c',evalSigma,'force_c',evalForce);

                % Debug
                %-------------------------------------------------------------%
                %[~,~] = contourf(horEvalNodes,verEvalNodes,evalForce',20);    
                %colormap(jet); title('Force in Assembler!')
                %-------------------------------------------------------------%

                %% ASSEMBLING LOOP
                %-------------------------------------------------------------%
                % The assembling loop creates each one of the submatrices
                % correponding to the microstructure of the linear system. The
                % loop is repeated m^2 times to complete the macrostructure of
                % the system.
                %-------------------------------------------------------------%

%                 for imb = 1:dimModBasis
%                     for kmb = 1:dimModBasis
% 
%                         [Amb,Mmb,bmb] = assemblerIGAScatterTransient(imb,kmb,augVerWeights,modalBasis(:,imb),modalBasisDer(:,imb),modalBasis(:,kmb),...
%                                                                      modalBasisDer(:,kmb),geoData,Computed,msh,space,jacFunc,spaceFunc);
% 
%                         % Assignment of the Block Matrix Just Assembled
% 
%                         A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;
%                         M(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Mmb;
% 
%                     end
% 
%                     % Assignment of the Block Vector Just Assembled
%                     b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;
%                 end

                warning('off','all')

                Amb = cell(dimModBasis,dimModBasis);
                Mmb = cell(dimModBasis,dimModBasis);
                bmb = cell(dimModBasis,dimModBasis);
                
                modalBasis2 = modalBasis;
                modalBasisDer2 = modalBasisDer;
                
                parfor imb = 1:dimModBasis
                    for kmb = 1:dimModBasis

                        [Amb{imb,kmb},Mmb{imb,kmb},bmb{imb,kmb}] = assemblerIGAScatterTransient(imb,kmb,augVerWeights,modalBasis(:,imb),modalBasisDer(:,imb),modalBasis2(:,kmb),...
                                                                     modalBasisDer2(:,kmb),geoData,Computed,msh,space,jacFunc,spaceFunc);
                                                                 
                    end
                end
                
                for imb = 1:dimModBasis
                    for kmb = 1:dimModBasis

                        % Assignment of the Block Matrix Just Assembled

                        A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb{imb,kmb};
                        M(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Mmb{imb,kmb};
                        
                    end
                    
                    % Assignment of the Block Vector Just Assembled
                    
                    b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb{imb,imb};
                    
                end

                %% IMPOSE BOUNDARY CONDITIONS

                %-------------------------------------------------------------%
                % Compute the lifting contribution in the force vector due to
                % the non homogeneous Dirichlet boundary conditions in the
                % lower and upper boundary.
                %-------------------------------------------------------------%

                BC_l         = boundCond.BC_INF_TAG;
                BC_r         = boundCond.BC_OUT_TAG;
                infBoundCond = boundCond.BC_INF_DATA;
                outBoundCond = boundCond.BC_OUT_DATA;

                obj_bcCoeff = BoundaryConditionHandler();

                obj_bcCoeff.infBoundCond    = infBoundCond;
                obj_bcCoeff.outBoundCond    = outBoundCond;
                obj_bcCoeff.augVerNodes     = augVerNodes;
                obj_bcCoeff.augVerWeights   = augVerWeights;
                obj_bcCoeff.modalBasis      = modalBasis;
                obj_bcCoeff.dimModalBasis   = dimModBasis;
                obj_bcCoeff.coefficientForm = obj.assemblerStruct.coefficientForm;

                [infStruct,outStruct] = computeFourierCoeff(obj_bcCoeff);

                bcStruct.bcInfTag       = BC_l;
                bcStruct.bcOutTag       = BC_r;
                bcStruct.infStruct      = infStruct;
                bcStruct.outStruct      = outStruct;
                bcStruct.numbControlPts = numbControlPts;
                bcStruct.dimModalBasis  = dimModBasis;

                [AA,MM,bb,sTS] = impose_boundary_transient(dimModBasis,BC_l, infStruct, BC_r, outStruct,A,M,b,numbControlPts,obj.assemblerStruct.currTimeSol);

            end
            
            %% Method 'buildConvergenceStruct'
            
            function [convStruct] = buildConvergenceStruct(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
                %% IMPORT CLASSES

                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %% NUMBER OF NODES IN THE QUADRATURE FORMULA
                %-------------------------------------------------------------%
                % The values of 'nqnx' and 'nqny' does not change during the
                % computation of the Gauss-Legendre Integration Points.
                % Attention, if you increase the number of nodes in the
                % discretization, it is necessary to refine the quadrature
                % rule.
                %-------------------------------------------------------------%

                % Horizontal direction

                numbHorNodes = obj.numbHorQuadNodes;
                assemblerStruct.numbHorNodes = numbHorNodes;

                % Vertical Direction

                numbVerNodes = obj.numbVerQuadNodes;
                assemblerStruct.numbVerNodes = numbVerNodes;

                %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
                %-------------------------------------------------------------%
                % Total number of knots used to divide the spline curve in the
                % isogeometric analysis.
                %-------------------------------------------------------------%

                numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
                assemblerStruct.numbKnots = numbKnots;

                %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
                %-------------------------------------------------------------%
                % Total number of control points used in the isogeometric
                % approach. Note that this number is equal to the dimension of
                % the basis used to define the Spline curve and that the
                % formula used bellow corresponds to use a straight line as
                % reference supporting fiber.
                %-------------------------------------------------------------%

                numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + 1 + obj.continuityParameter;
                assemblerStruct.numbControlPts = numbControlPts;

                %% GAUSS-LEGENDRE INTEGRATION NODES
                %-------------------------------------------------------------%
                % The following method reveives the number of integration nodes
                % in the horizontal and vertical direction and retruns the
                % respective Gauss-Legendre nodes and weigths for the
                % integration interval [0,1].   
                %
                % Mesh FEM in X: Equispaced Nodes                            
                %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
                %                    Direction
                %  (2) horGLWeights : Vector of the Standard Weights on the X 
                %                     Direction
                %
                % Mesh FEM in Y: Equispaced Nodes                              
                %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
                %                    Direction  
                %  (2) verGLWeights : Vector of the Standard Weights on the Y 
                %                     Direction
                %-------------------------------------------------------------%

                % Vertical direction

                obj_gaussLegendre_2 = IntegrateHandler();
                obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
                [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);   

                assemblerStruct.verGLNodes = verGLNodes;
                assemblerStruct.verWeights = verWeights;

                %% EXTRACT GEOMETRIC INFORMATION

                geometry  = obj.geometricInfo.geometry;
                map       = @(x,y) obj.geometricInfo.map(x,y);
                Jac       = @(x,y) obj.geometricInfo.Jac(x,y);
                Hes       = @(x,y) obj.geometricInfo.Hes(x,y);

                assemblerStruct.geometry = geometry;
                assemblerStruct.map      = map;
                assemblerStruct.Jac      = Jac;
                assemblerStruct.Hes      = Hes;

                %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
                % Note: Create the knots and control points of the reference
                % supporting fiber where the coupled 1D problems will be
                % solved.

                refDomain1D = geo_load(nrbline ([0 0], [1 0]));

                nsub = numbKnots;
                degree = obj.degreePolySplineBasis;
                regularity = obj.continuityParameter;

                [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);

                assemblerStruct.refDomain1D = refDomain1D;
                assemblerStruct.nsub = nsub;
                assemblerStruct.degree = degree;
                assemblerStruct.regularity = regularity;
                assemblerStruct.knots = knots;
                assemblerStruct.zeta = zeta;

                %% GENERATE ISOGEOMETRIC MESH FUNCTION
                % Note: Use the number of quadrature nodes to generate the
                % quadrature rule using the gauss nodes. Then, use this
                % information with the computed knots and control points to
                % generate the isogeometric mesh of the centreline.

                rule     = msh_gauss_nodes (numbHorNodes);
                [qn, qw] = msh_set_quad_nodes (zeta, rule);
                msh      = msh_cartesian (zeta, qn, qw, refDomain1D);

                assemblerStruct.rule = rule;
                assemblerStruct.qn   = qn;
                assemblerStruct.qw   = qw;
                assemblerStruct.msh  = msh;

                %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
                % Note: Construct the functional space used to discretize the
                % problem along the centreline and perform the isogeometric
                % analysis. 
                % Here, we also evaluate the functional spaces in the specific
                % element. This will be used in the assembling loop to compute
                % the operators using the basis functions.

                space    = sp_bspline (knots, degree, msh);

                numbKnots = msh.nel_dir;

                for iel = 1:numbKnots

                    msh_col = msh_evaluate_col (msh, iel);

                    %---------------------------------------------------------%
                    % Evaluated space to compute the operators
                    %---------------------------------------------------------%

                    gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                    shapFunc = sp_evaluate_col (space, msh_col);

                    %---------------------------------------------------------%
                    % Store the spaces in the cell matrix
                    %---------------------------------------------------------%

                    spaceFunc{iel,1} = shapFunc;
                    spaceFunc{iel,2} = gradFunc;
                    spaceFunc{iel,3} = msh_col;

                end

                assemblerStruct.space     = space;
                assemblerStruct.numbKnots = numbKnots;
                assemblerStruct.spaceFunc = spaceFunc;

                %% AUGMENTED HORIZONTAL POINTS 
                % Note: The FOR loop runs over the knots in the centerline ans
                % compute the required quadrature nodes. The quadrature nodes
                % will be necessary to evaluate the bilinear coefficients and
                % the forcing component.

                horEvalNodes = zeros(numbKnots * numbHorNodes,1);

                for iel = 1:numbKnots

                    msh_col = msh_evaluate_col (msh, iel);

                    localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                    horEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  

                end

                assemblerStruct.horEvalNodes = horEvalNodes;

                %% AUGMENTED VERTICAL POINTS
                % Note: Since we are using the full map from the physical
                % domain to the reference domain, then the whole modal analysis
                % has to be performed in the vertical direction in the interval
                % [0,1];

                objVertQuadRule = IntegrateHandler();

                objVertQuadRule.leftBoundInterval = 0;
                objVertQuadRule.rightBoundInterval = 1;
                objVertQuadRule.inputNodes = verGLNodes;
                objVertQuadRule.inputWeights = verWeights;

                [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);

                assemblerStruct.augVerNodes = augVerNodes;
                assemblerStruct.augVerWeights = augVerWeights;

                %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
                %-------------------------------------------------------------%
                % The next method takes the coefficients of the bilinear form,
                % the dimension of the modal basis (assigned in the demo) and
                % the vertical evaluation points and computes the coefficients
                % of the modal basis when evaluated in the points of the
                % transverse fiber.
                %-------------------------------------------------------------%

                obj_newModalBasis = BasisHandler();

                obj_newModalBasis.dimModalBasis      = obj.dimModalBasis;
                obj_newModalBasis.evalNodesY         = verGLNodes;
                obj_newModalBasis.labelUpBoundCond   = obj.label_upBoundDomain;
                obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
                obj_newModalBasis.coeffForm          = obj.coefficientForm;

                [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);

                assemblerStruct.modalBasis          = modalBasis;
                assemblerStruct.modalBasisDer       = modalBasisDer;
                assemblerStruct.dimModalBasis       = obj.dimModalBasis;
                assemblerStruct.labelUpBoundCond    = obj.label_upBoundDomain;
                assemblerStruct.labelDownBoundCond  = obj.label_downBoundDomain;
                
                %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
                %-------------------------------------------------------------%
                % We use the horizontal mesh created previously to evaluate the
                % thickness and the vertical coordinate of the centerline of
                % the channel.
                %-------------------------------------------------------------%

                verEvalNodes = augVerNodes;

                X = mapOut(horEvalNodes,verEvalNodes,map,1);
                Y = mapOut(horEvalNodes,verEvalNodes,map,2);

                geoData.X = X;
                geoData.Y = Y;
                geoData.horNodes = horEvalNodes;
                geoData.verNodes = verEvalNodes;

                assemblerStruct.geoData = geoData;

                %% EVALUATION OF THE GEOMETRY PROPERTIES
                % Note: Since there is a transformation from the physical domain to the
                % physical domain, we must compute the map contribution used in the
                % computation of the coefficients and the Jacobian contribution to be
                % inserted in the integral (quadrature formula).

                evalJac     = jacOut(horEvalNodes,verEvalNodes,Jac);
                evalDetJac  = detJac(horEvalNodes,verEvalNodes,Jac);

                Phi1_dx = invJac(1,1,evalJac);
                Phi1_dy = invJac(1,2,evalJac);
                Phi2_dx = invJac(2,1,evalJac);
                Phi2_dy = invJac(2,2,evalJac);

                jacFunc.evalJac     = evalJac;
                jacFunc.evalDetJac  = evalDetJac;
                jacFunc.Phi1_dx     = Phi1_dx;
                jacFunc.Phi1_dy     = Phi1_dy;
                jacFunc.Phi2_dx     = Phi2_dx;
                jacFunc.Phi2_dy     = Phi2_dy;

                assemblerStruct.jacFunc = jacFunc;
                
                %% BOUNDARY CONDITIONS PROPERTIES
                
                assemblerStruct.igaBoundCond = obj.dirCondFuncStruct.igaBoundCond;
                
                %% MODEL PROPERTIES
                
                assemblerStruct.coefficientForm   = obj.coefficientForm;
                assemblerStruct.dirCondFuncStruct = obj.dirCondFuncStruct;
                
                %% GENERATION OF PLOTTING STRUCTURE
                
                assemblerStruct.plotResolutionX = 100;
                assemblerStruct.plotResolutionY = 100;

            end
            
            %% Method 'buildIGAScatterStabilized'
            
            function [A,b,modalBasis,liftCoeffA,liftCoeffB,space,refDomain1D] = buildSystemIGAScatterStabilized(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler
            
            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve and that the
            % formula used bellow corresponds to use a straight line as
            % reference supporting fiber.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + ...
                             1 + obj.continuityParameter;

            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);        
            
            %% EXTRACT GEOMETRIC INFORMATION
            
            geometry  = obj.geometricInfo.geometry;
            map       = @(x,y) obj.geometricInfo.map(x,y);
            Jac       = @(x,y) obj.geometricInfo.Jac(x,y);
            Hes       = @(x,y) obj.geometricInfo.Hes(x,y);
            
            disp('Finished EXTRACT GEOMETRIC INFORMATION');
            
            %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
            % Note: Create the knots and control points of the reference
            % supporting fiber where the coupled 1D problems will be
            % solved.
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsub = numbKnots;
            degree = obj.degreePolySplineBasis;
            regularity = obj.continuityParameter;
            
            [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);
            
            disp('Finished COMPUTE REFERENCE KNOT AND CONTROL POINTS');
            
            %% GENERATE ISOGEOMETRIC MESH FUNCTION
            % Note: Use the number of quadrature nodes to generate the
            % quadrature rule using the gauss nodes. Then, use this
            % information with the computed knots and control points to
            % generate the isogeometric mesh of the centreline.
            
            rule     = msh_gauss_nodes (numbHorNodes);
            [qn, qw] = msh_set_quad_nodes (zeta, rule);
            msh      = msh_cartesian (zeta, qn, qw, refDomain1D);
            
            disp('Finished GENERATE ISOGEOMETRIC MESH FUNCTION');
            
            %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
            % Note: Construct the functional space used to discretize the
            % problem along the centreline and perform the isogeometric
            % analysis. 
            % Here, we also evaluate the functional spaces in the specific
            % element. This will be used in the assembling loop to compute
            % the operators using the basis functions.
            
            space    = sp_bspline (knots, degree, msh);
            
            numbKnots = msh.nel_dir;
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                shapFunc = sp_evaluate_col (space, msh_col);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFunc{iel,1} = shapFunc;
                spaceFunc{iel,2} = gradFunc;
                spaceFunc{iel,3} = msh_col;

            end
            
            disp('Finished CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE');
            
            %% AUGMENTED HORIZONTAL POINTS 
            % Note: The FOR loop runs over the knots in the centerline ans
            % compute the required quadrature nodes. The quadrature nodes
            % will be necessary to evaluate the bilinear coefficients and
            % the forcing component.

            horEvalNodes = zeros(numbKnots * numbHorNodes,1);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                horEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  
                
            end
            
            disp('Finished AUGMENTED HORIZONTAL POINTS');
            
            %% AUGMENTED VERTICAL POINTS
            % Note: Since we are using the full map from the physical
            % domain to the reference domain, then the whole modal analysis
            % has to be performed in the vertical direction in the interval
            % [0,1];
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = 0;
            objVertQuadRule.rightBoundInterval = 1;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            disp('Finished AUGMENTED VERTICAL POINTS');
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);
            
            disp('Finished COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION');
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%
            
            verEvalNodes = augVerNodes;
            
            X = mapOut(horEvalNodes,verEvalNodes,map,1);
            Y = mapOut(horEvalNodes,verEvalNodes,map,2);
            
            geoData.X = X;
            geoData.Y = Y;
            geoData.horNodes = horEvalNodes;
            geoData.verNodes = verEvalNodes;
            geoData.numbKnots = numbKnots;
            geoData.numbInt = nsub;
            geoData.numbModes = obj.dimModalBasis;
            
            disp('Finished EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN');
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We need to evaluate all the coefficients of the bilinear form
            % in the quadrature nodes along the vertical direction to
            % perform the first integral (along the transverse fiber) to
            % obtain a the coefficients of the 1D coupled problem. The
            % result will be coefficients as a function of 'x'.
            %-------------------------------------------------------------%
            

            % DIFFUSION

            evalMu    = obj.coefficientForm.mu(X,Y);
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(X,Y);
            evalBeta2 = obj.coefficientForm.beta2(X,Y);
            
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(X,Y);
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(X,Y);
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%
            
            Computed.mu_c = evalMu;
            Computed.beta1_c = evalBeta1;
            Computed.beta2_c = evalBeta2;
            Computed.sigma_c = evalSigma;
            Computed.force_c = evalForce;
            Computed.y = verEvalNodes;
%             Computed.normL2beta = obj.coefficientForm.normL2beta;
%             Computed.normL2mu = obj.coefficientForm.normL2mu;
                          
            disp('Finished EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN');
            
            %% EVALUATION OF THE GEOMETRY PROPERTIES
            % Note: Since there is a transformation from the physical domain to the
            % physical domain, we must compute the map contribution used in the
            % computation of the coefficients and the Jacobian contribution to be
            % inserted in the integral (quadrature formula).
            
            evalJac     = jacOut(horEvalNodes,verEvalNodes,Jac);
            evalDetJac  = detJac(horEvalNodes,verEvalNodes,Jac);
            
            Phi1_dx = invJac(1,1,evalJac);
            Phi1_dy = invJac(1,2,evalJac);
            Phi2_dx = invJac(2,1,evalJac);
            Phi2_dy = invJac(2,2,evalJac);
            
            jacFunc.evalJac     = evalJac;
            jacFunc.evalDetJac  = evalDetJac;
            jacFunc.Phi1_dx     = Phi1_dx;
            jacFunc.Phi1_dy     = Phi1_dy;
            jacFunc.Phi2_dx     = Phi2_dx;
            jacFunc.Phi2_dy     = Phi2_dy;
            
            disp('Finished EVALUATION OF THE GEOMETRY PROPERTIES');
                                                    
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y);
            
            disp('Finished LIFTING');
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis
                for kmb = 1:obj.dimModalBasis

                    [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerIGAScatterStabilized( imb, kmb, ...
                                            augVerWeights,modalBasis(:,imb),modalBasisDer(:,imb),modalBasis(:,kmb),...
                                            modalBasisDer(:,kmb),geoData,Computed,lifting,aLift,bLift,msh,space,jacFunc,...
                                            spaceFunc,obj.dirCondFuncStruct.igaBoundCond,obj.stabMethod,obj.stabDelta,obj.PeConfig);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;
                    
                    disp(['FINISHED ASSEMBLING LOOP (',num2str(imb),' , ',num2str(kmb),')']);

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;
            end

            end
            
            %% Method 'buildIGAScatter3D'
            
            function [AA,bb,bcStruct,liftCoeffA,liftCoeffB,space,refDomain1D] = buildSystemIGAScatter3D(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler
            
            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve and that the
            % formula used bellow corresponds to use a straight line as
            % reference supporting fiber.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + ...
                             1 + obj.continuityParameter;

            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);   
            
            disp('Finished SETTING PARAMETERS')
            
            %% EXTRACT GEOMETRIC INFORMATION
            % Note: This section was modified to account for the third
            % coordinate of the domain.
            
            geometry  = obj.geometricInfo.geometry;
            map       = @(x,y,z) geometry.map([x;y;z]);
            Jac       = @(x,y,z) geometry.map_der([x;y;z]);
            Hes       = @(x,y,z) geometry.map_der2([x;y;z]);
            type      = obj.geometricInfo.Type;
            
            disp('Finished EXTRACT GEOMETRIC INFORMATION')
            
            %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
            % Note: Create the knots and control points of the reference
            % supporting fiber where the coupled 1D problems will be
            % solved.
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsub = numbKnots;
            degree = obj.degreePolySplineBasis;
            regularity = obj.continuityParameter;
            
            [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);
            
            disp('Finished COMPUTE REFERENCE KNOT AND CONTROL POINTS')
            
            %% GENERATE ISOGEOMETRIC MESH FUNCTION
            % Note: Use the number of quadrature nodes to generate the
            % quadrature rule using the gauss nodes. Then, use this
            % information with the computed knots and control points to
            % generate the isogeometric mesh of the centreline.
            
            rule     = msh_gauss_nodes (numbHorNodes);
            [qn, qw] = msh_set_quad_nodes (zeta, rule);
            msh      = msh_cartesian (zeta, qn, qw, refDomain1D);
            
            disp('Finished GENERATE ISOGEOMETRIC MESH FUNCTION')
            
            %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
            % Note: Construct the functional space used to discretize the
            % problem along the centreline and perform the isogeometric
            % analysis. 
            % Here, we also evaluate the functional spaces in the specific
            % element. This will be used in the assembling loop to compute
            % the operators using the basis functions.
            
            space    = sp_bspline (knots, degree, msh);
            
            numbKnots = msh.nel_dir;
            spaceFunc = cell(numbKnots,3);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                shapFunc = sp_evaluate_col (space, msh_col);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFunc{iel,1} = shapFunc;
                spaceFunc{iel,2} = gradFunc;
                spaceFunc{iel,3} = msh_col;

            end
            
            disp('Finished CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE')
            
            %% AUGMENTED HORIZONTAL POINTS 
            % Note: The FOR loop runs over the knots in the centerline ans
            % compute the required quadrature nodes. The quadrature nodes
            % will be necessary to evaluate the bilinear coefficients and
            % the forcing component.

            suppEvalNodes = zeros(numbKnots * numbHorNodes,1);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                suppEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  
                
            end
            
            %% AUGMENTED VERTICAL POINTS
            % Note: This section had to be modified to the change in the
            % transverse fiber. The transverse fiber now is a square cross
            % section of the reference domain.
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = 0;
            objVertQuadRule.rightBoundInterval = 1;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            disp('Finished AUGMENTED POINTS')
            
            %% COMPUTATION OF THE MODAL BASIS IN TRANSVERSE DOMAIN
            % Note: To deal with two dimensional problems the modal basis
            % now have to be computed as a function of (y,z) and not only
            % in the direction y. However, since the reference domain is a
            % cube, the transverse polygon will be an unit square and to
            % define the modal basis is only necessary to multiply the
            % modal basis along "y" and the one along "z".
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;
            
            [modalBasis, modalBasisDer1, modalBasisDer2] = newModalBasis3D(obj_newModalBasis);
            
            modalBasisStruct.modalBasis = modalBasis;
            modalBasisStruct.modalBasisDer1 = modalBasisDer1;
            modalBasisStruct.modalBasisDer2 = modalBasisDer2;
            modalBasisStruct.numbModes = obj.dimModalBasis;

            x = verGLNodes;
            y = verGLNodes;
            [X,Y] = meshgrid(x,y);
            
            figure;
            for ii = 1:obj.dimModalBasis^2
                subplot(obj.dimModalBasis,obj.dimModalBasis,ii)
                surf(X,Y,modalBasis(:,:,ii));
                az = 45;
                el = 30;
                view(az, el);
            end

            disp('Finished COMPUTATION OF THE MODAL BASIS IN TRANSVERSE DOMAIN')
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % Note: In the 3D case, we have to change all the avaluation
            % functions and the map and Jacobian contributions to consider
            % the extra dimension along the transverse direction.
            %-------------------------------------------------------------%
            
            sigmaEvalNodes = augVerNodes;
            
            [X,Y,Z] = mapOut3DHiMod(suppEvalNodes,sigmaEvalNodes,sigmaEvalNodes,obj.geometricInfo,type);
            
            geoData.X = X;
            geoData.Y = Y;
            geoData.Z = Z;
            geoData.horNodes = suppEvalNodes;
            geoData.verNodes = sigmaEvalNodes;
            
            disp('Finished EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN')
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We need to evaluate all the coefficients of the bilinear form
            % in the quadrature nodes along the vertical direction to
            % perform the first integral (along the transverse fiber) to
            % obtain a the coefficients of the 1D coupled problem. The
            % result will be coefficients as a function of 'x'.
            %-------------------------------------------------------------%
            

            % DIFFUSION

            evalMu    = obj.coefficientForm.mu(X,Y,Z);
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(X,Y,Z);
            evalBeta2 = obj.coefficientForm.beta2(X,Y,Z);
            evalBeta3 = obj.coefficientForm.beta3(X,Y,Z);
                         
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(X,Y,Z);
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(X,Y,Z);
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed.mu_c    = evalMu;
            Computed.beta1_c = evalBeta1;
            Computed.beta2_c = evalBeta2;
            Computed.beta3_c = evalBeta3;
            Computed.sigma_c = evalSigma;
            Computed.force_c = evalForce;
            Computed.y       = sigmaEvalNodes;
            
            disp('Finished EVALUATION OF COEFFICIENTS')
            
            %% EVALUATION OF THE GEOMETRY PROPERTIES
            % Note: Since there is a transformation from the physical domain to the
            % physical domain, we must compute the map contribution used in the
            % computation of the coefficients and the Jacobian contribution to be
            % inserted in the integral (quadrature formula).
            
            tic;
            [evalJac,structPsi,evalDetJac] = jacOut3DHiMod(suppEvalNodes,sigmaEvalNodes, ...
                                                      sigmaEvalNodes,obj.geometricInfo,type);   

            Psi1_dx = structPsi.Psi1_dx;
            Psi1_dy = structPsi.Psi1_dy;
            Psi1_dz = structPsi.Psi1_dz;
            Psi2_dx = structPsi.Psi2_dx;
            Psi2_dy = structPsi.Psi2_dy;
            Psi2_dz = structPsi.Psi2_dz;
            Psi3_dx = structPsi.Psi3_dx;
            Psi3_dy = structPsi.Psi3_dy;
            Psi3_dz = structPsi.Psi3_dz;
            
            T1 = toc;
            
            disp(['FINISHED EVALUATING JACOBIAN in t = ',num2str(T1),' s'])
            
            jacFunc.evalJac     = evalJac;
            jacFunc.evalDetJac  = evalDetJac;
            jacFunc.Psi1_dx     = Psi1_dx;
            jacFunc.Psi1_dy     = Psi1_dy;
            jacFunc.Psi1_dz     = Psi1_dz;
            jacFunc.Psi2_dx     = Psi2_dx;
            jacFunc.Psi2_dy     = Psi2_dy;
            jacFunc.Psi2_dz     = Psi2_dz;
            jacFunc.Psi3_dx     = Psi3_dx;
            jacFunc.Psi3_dy     = Psi3_dy;
            jacFunc.Psi3_dz     = Psi3_dz;
            
            disp('Finished EVALUATION OF THE JACOBIAN PROPERTIES')
                                                    
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary. The lifting function must change to
            % account for the extra dimension along the transverse
            % direction.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y);
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis^2
                for kmb = 1:obj.dimModalBasis^2

                    [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerIGAScatter3D( imb, kmb, ...
                                            augVerWeights,modalBasis(:,:,imb),modalBasisDer1(:,:,imb),modalBasisDer2(:,:,imb),modalBasis(:,:,kmb),...
                                            modalBasisDer1(:,:,kmb),modalBasisDer2(:,:,kmb),geoData,Computed,lifting,aLift,bLift,msh,space,jacFunc,...
                                            spaceFunc);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;
                    
                    if (imb == kmb)
                        disp(['FINISHED ASSEMBLING LOOP (',num2str(imb),' , ',num2str(kmb),')']);
                    end

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;
            end
            
            %% IMPOSE BOUNDARY CONDITIONS
            
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            BC_l = obj.igaBoundCond.BC_INF_TAG;
            BC_r = obj.igaBoundCond.BC_OUT_TAG;
            infBoundCond = obj.igaBoundCond.BC_INF_DATA;
            outBoundCond = obj.igaBoundCond.BC_OUT_DATA;
            
            obj_bcCoeff = BoundaryConditionHandler();
    
            obj_bcCoeff.infBoundCond    = infBoundCond;
            obj_bcCoeff.outBoundCond    = outBoundCond;
            obj_bcCoeff.augVerNodes     = augVerNodes;
            obj_bcCoeff.augVerWeights   = augVerWeights;
            obj_bcCoeff.modalBasis      = modalBasis;
            obj_bcCoeff.dimModalBasis   = obj.dimModalBasis;
            obj_bcCoeff.coefficientForm = obj.coefficientForm;

            [infStruct,outStruct] = computeFourierCoeff3D(obj_bcCoeff);
            
            bcStruct.bcInfTag = BC_l;
            bcStruct.bcOutTag = BC_r;
            bcStruct.infStruct = infStruct;
            bcStruct.outStruct = outStruct;
            bcStruct.numbControlPts = numbControlPts;
            bcStruct.dimModalBasis = obj.dimModalBasis;
            
            [AA,bb] = impose_boundary(obj.dimModalBasis,BC_l, infStruct, BC_r, outStruct,A,b,numbControlPts);
            
            end
            
            %% Method 'buildIGAScatter3DTransient'
            
            function [A,M,b,modalBasisStruct,liftCoeffA,liftCoeffB,space,refDomain1D] = buildSystemIGAScatter3DTransient(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler
            
            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve and that the
            % formula used bellow corresponds to use a straight line as
            % reference supporting fiber.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*(obj.degreePolySplineBasis - obj.continuityParameter) + ...
                             1 + obj.continuityParameter;

            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);   
            
            disp('Finished SETTING PARAMETERS')
            
            %% EXTRACT GEOMETRIC INFORMATION
            % Note: This section was modified to account for the third
            % coordinate of the domain.
            
            geometry  = obj.geometricInfo.geometry;
            map       = @(x,y,z) geometry.map([x;y;z]);
            Jac       = @(x,y,z) geometry.map_der([x;y;z]);
            Hes       = @(x,y,z) geometry.map_der2([x;y;z]);
            type      = obj.geometricInfo.Type;
            
            disp('Finished EXTRACT GEOMETRIC INFORMATION')
            
            %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
            % Note: Create the knots and control points of the reference
            % supporting fiber where the coupled 1D problems will be
            % solved.
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsub = numbKnots;
            degree = obj.degreePolySplineBasis;
            regularity = obj.continuityParameter;
            
            [knots, zeta] = kntrefine (refDomain1D.nurbs.knots, nsub-1, degree, regularity);
            
            disp('Finished COMPUTE REFERENCE KNOT AND CONTROL POINTS')
            
            %% GENERATE ISOGEOMETRIC MESH FUNCTION
            % Note: Use the number of quadrature nodes to generate the
            % quadrature rule using the gauss nodes. Then, use this
            % information with the computed knots and control points to
            % generate the isogeometric mesh of the centreline.
            
            rule     = msh_gauss_nodes (numbHorNodes);
            [qn, qw] = msh_set_quad_nodes (zeta, rule);
            msh      = msh_cartesian (zeta, qn, qw, refDomain1D);
            
            disp('Finished GENERATE ISOGEOMETRIC MESH FUNCTION')
            
            %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
            % Note: Construct the functional space used to discretize the
            % problem along the centreline and perform the isogeometric
            % analysis. 
            % Here, we also evaluate the functional spaces in the specific
            % element. This will be used in the assembling loop to compute
            % the operators using the basis functions.
            
            space    = sp_bspline (knots, degree, msh);
            
            numbKnots = msh.nel_dir;
            spaceFunc = cell(numbKnots,3);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFunc = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
                shapFunc = sp_evaluate_col (space, msh_col);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFunc{iel,1} = shapFunc;
                spaceFunc{iel,2} = gradFunc;
                spaceFunc{iel,3} = msh_col;

            end
            
            disp('Finished CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE')
            
            %% AUGMENTED HORIZONTAL POINTS 
            % Note: The FOR loop runs over the knots in the centerline ans
            % compute the required quadrature nodes. The quadrature nodes
            % will be necessary to evaluate the bilinear coefficients and
            % the forcing component.

            suppEvalNodes = zeros(numbKnots * numbHorNodes,1);
            
            for iel = 1:numbKnots
                
                msh_col = msh_evaluate_col (msh, iel);
                
                localNodes = reshape (msh_col.geo_map(1,:,:), msh_col.nqn, msh_col.nel);  
                suppEvalNodes((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes) = localNodes;  
                
            end
            
            %% AUGMENTED VERTICAL POINTS
            % Note: This section had to be modified to the change in the
            % transverse fiber. The transverse fiber now is a square cross
            % section of the reference domain.
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = 0;
            objVertQuadRule.rightBoundInterval = 1;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            disp('Finished AUGMENTED POINTS')
            
            %% COMPUTATION OF THE MODAL BASIS IN TRANSVERSE DOMAIN
            % Note: To deal with two dimensional problems the modal basis
            % now have to be computed as a function of (y,z) and not only
            % in the direction y. However, since the reference domain is a
            % cube, the transverse polygon will be an unit square and to
            % define the modal basis is only necessary to multiply the
            % modal basis along "y" and the one along "z".
            
            % (1) CHEBCHEV POLYNOMIALS
            
            % obj_newModalBasis = BasisHandler();
            
            % obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            % obj_newModalBasis.evalNodesY = verGLNodes;
            
            % [modalBasis, modalBasisDer1, modalBasisDer2] = newModalBasisCheb3D(obj_newModalBasis);
            
            % (2) LEGENDRE POLYNOMIALS
            
%             obj_newModalBasis = BasisHandler();
%             
%             obj_newModalBasis.dimLegendreBase = obj.dimModalBasis;
%             obj_newModalBasis.evalLegendreNodes = verGLNodes;
%             obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
%             obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
%             obj_newModalBasis.coeffForm = obj.coefficientForm;
%             
%             [modalBasis, modalBasisDer1, modalBasisDer2] = newModalBasisLegendre3D(obj_newModalBasis);
%             
%             modalBasisStruct.modalBasis = modalBasis;
%             modalBasisStruct.modalBasisDer1 = modalBasisDer1;
%             modalBasisStruct.modalBasisDer2 = modalBasisDer2;
%             modalBasisStruct.numbModes = obj.dimModalBasis;
            
%             size(modalBasis)
%             
%             x = verGLNodes;
%             y = verGLNodes;
%             [X,Y] = meshgrid(x,y);
%             
%             for ii = 1:obj.dimModalBasis^2
%                 figure;
%                 surf(X,Y,modalBasis(:,:,ii));
%                 az = 45;
%                 el = 30;
%                 view(az, el);
%             end
%             
%             return
            
            % (3) EDUCATED BASIS
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;
            
            [modalBasis, modalBasisDer1, modalBasisDer2] = newModalBasis3D(obj_newModalBasis);
            
            modalBasisStruct.modalBasis = modalBasis;
            modalBasisStruct.modalBasisDer1 = modalBasisDer1;
            modalBasisStruct.modalBasisDer2 = modalBasisDer2;
            modalBasisStruct.numbModes = obj.dimModalBasis;

            x = verGLNodes;
            y = verGLNodes;
            [X,Y] = meshgrid(x,y);
            
%             figure;
%             for ii = 1:obj.dimModalBasis^2
%                 subplot(obj.dimModalBasis,obj.dimModalBasis,ii)
%                 surf(X,Y,modalBasis(:,:,ii));
%                 az = 45;
%                 el = 30;
%                 view(az, el);
%             end

            disp('Finished COMPUTATION OF THE MODAL BASIS IN TRANSVERSE DOMAIN')
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % Note: In the 3D case, we have to change all the avaluation
            % functions and the map and Jacobian contributions to consider
            % the extra dimension along the transverse direction.
            %-------------------------------------------------------------%
            
            sigmaEvalNodes = augVerNodes;
            
            [X,Y,Z] = mapOut3DHiMod(suppEvalNodes,sigmaEvalNodes,sigmaEvalNodes,obj.geometricInfo,type);
            
            geoData.X = X;
            geoData.Y = Y;
            geoData.Z = Z;
            geoData.horNodes = suppEvalNodes;
            geoData.verNodes = sigmaEvalNodes;
            
            disp('Finished EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN')
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We need to evaluate all the coefficients of the bilinear form
            % in the quadrature nodes along the vertical direction to
            % perform the first integral (along the transverse fiber) to
            % obtain a the coefficients of the 1D coupled problem. The
            % result will be coefficients as a function of 'x'.
            %-------------------------------------------------------------%
            

            % DIFFUSION

            evalMu    = obj.coefficientForm.mu(X,Y,Z,obj.timeInstant);
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(X,Y,Z,obj.timeInstant);
            evalBeta2 = obj.coefficientForm.beta2(X,Y,Z,obj.timeInstant);
            evalBeta3 = obj.coefficientForm.beta3(X,Y,Z,obj.timeInstant);
                         
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(X,Y,Z,obj.timeInstant);
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(X,Y,Z,obj.timeInstant);
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed.mu_c    = evalMu;
            Computed.beta1_c = evalBeta1;
            Computed.beta2_c = evalBeta2;
            Computed.beta3_c = evalBeta3;
            Computed.sigma_c = evalSigma;
            Computed.force_c = evalForce;
            Computed.y       = sigmaEvalNodes;
            
            disp('Finished EVALUATION OF COEFFICIENTS')
            
            %% EVALUATION OF THE GEOMETRY PROPERTIES
            % Note: Since there is a transformation from the physical domain to the
            % physical domain, we must compute the map contribution used in the
            % computation of the coefficients and the Jacobian contribution to be
            % inserted in the integral (quadrature formula).
            
            tic;
            [evalJac,structPsi,evalDetJac] = jacOut3DHiMod(suppEvalNodes,sigmaEvalNodes, ...
                                                      sigmaEvalNodes,obj.geometricInfo,type);   

            Psi1_dx = structPsi.Psi1_dx;
            Psi1_dy = structPsi.Psi1_dy;
            Psi1_dz = structPsi.Psi1_dz;
            Psi2_dx = structPsi.Psi2_dx;
            Psi2_dy = structPsi.Psi2_dy;
            Psi2_dz = structPsi.Psi2_dz;
            Psi3_dx = structPsi.Psi3_dx;
            Psi3_dy = structPsi.Psi3_dy;
            Psi3_dz = structPsi.Psi3_dz;
            
            T1 = toc;
            
            disp(['FINISHED EVALUATING JACOBIAN in t = ',num2str(T1),' s'])
            
            jacFunc.evalJac     = evalJac;
            jacFunc.evalDetJac  = evalDetJac;
            jacFunc.Psi1_dx     = Psi1_dx;
            jacFunc.Psi1_dy     = Psi1_dy;
            jacFunc.Psi1_dz     = Psi1_dz;
            jacFunc.Psi2_dx     = Psi2_dx;
            jacFunc.Psi2_dy     = Psi2_dy;
            jacFunc.Psi2_dz     = Psi2_dz;
            jacFunc.Psi3_dx     = Psi3_dx;
            jacFunc.Psi3_dy     = Psi3_dy;
            jacFunc.Psi3_dz     = Psi3_dz;
            
            disp('Finished EVALUATION OF THE JACOBIAN PROPERTIES')
                                                    
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary. The lifting function must change to
            % account for the extra dimension along the transverse
            % direction.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y);
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis^2
                for kmb = 1:obj.dimModalBasis^2

                    [Amb,Mmb,bmb,liftCoeffA,liftCoeffB] = assemblerIGAScatter3DTransient( imb, kmb, ...
                                            augVerWeights,modalBasis(:,:,imb),modalBasisDer1(:,:,imb),modalBasisDer2(:,:,imb),modalBasis(:,:,kmb),...
                                            modalBasisDer1(:,:,kmb),modalBasisDer2(:,:,kmb),geoData,Computed,lifting,aLift,bLift,msh,space,jacFunc,...
                                            spaceFunc);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;
                    M(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Mmb;
                    
                    if (imb == kmb)
                        disp(['FINISHED ASSEMBLING LOOP (',num2str(imb),' , ',num2str(kmb),')']);
                    end

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;
            end
            
            end
        
            %% Method 'buildIGAForce'
            
            function [b,modalBasis,liftCoeffA,liftCoeffB,intNodesGaussLeg,intNodesWeights] = buildIGAForce(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Stiffness Matrix Using IGA Basis
            %   (2) M                   : Final Assembled Mass Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler

            % Number of Nodes for the Quadrature Formula

            %---------------------------------------------------------------------%
            % Note:
            % Attention, if you increase the number nodes it is also necessary to
            % refine the quadrature mesh.
            %---------------------------------------------------------------------%

            nqnx = 8;   
            nqny = 64;

            ne  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);   % Number of Intervals
            nx  = ne+1;                                                              % Number of Nodes     
            ncp = ne*obj.continuityParameter + obj.degreePolySplineBasis + 1 - obj.continuityParameter;   % Spline Basis Dimension (Number of Control Points) 

            % Gauss-Legendre Integration Nodes

            %---------------------------------------------------------------------%
            % Note:
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            %---------------------------------------------------------------------%

            obj_gaussLegendre_1 = IntegrateHandler();
            obj_gaussLegendre_2 = IntegrateHandler();
                
            obj_gaussLegendre_1.numbQuadNodes = nqnx;
            obj_gaussLegendre_2.numbQuadNodes = nqny;
            
            [~, xq, wxq] = gaussLegendre(obj_gaussLegendre_1); 
            [~, intNodesGaussLeg, intNodesWeights] = gaussLegendre(obj_gaussLegendre_2); 

            %---------------------------------------------------------------------%
            % Description of Obtained Data:                                       %
            %---------------------------------------------------------------------%                    
            % Mesh FEM in X: Equispaced Nodes                                     %
            %   (1) xq  : Vector of the 8 Standard Nodes on the X Direction       %
            %   (2) wxq : Vector of the 8 Standard Weights on the X Direction     %
            % Mesh FEM in Y: Equispaced Nodes                                     %
            %   (1) yq  : Vector of the 8 Standard Nodes on the Y Direction       %
            %   (2) wyq : Vector of the 8 Standard Weights on the Y Direction     %
            %---------------------------------------------------------------------%

            % FEM MESH IN THE X DIRECTION
            %---------------------------------------------------------------------%
            % Note: In this case we are considering equispaced nodes in the mesh.
            %---------------------------------------------------------------------%

            meshx     = zeros(nx,1);
            meshx(1)  = obj.leftBDomain_inX;
            meshx2    = zeros(nx,1);
            meshx2(1) = obj.leftBDomain_inX;

            for i=2:nx
                meshx(i) = meshx(i-1)+obj.stepMeshX;
            end

            for i=2:nx
                meshx2(i) = meshx2(i-1)+obj.jacAtQuadNodes(i-1); 
            end

            %---------------------------------------------------------------------%
            %                  QUADRATURE NODES ALONG THE MESH                    %
            %---------------------------------------------------------------------%
            % Notes:                                                              %
            %                                                                     %
            % Each interval has its nodes and weights. Both meshes, based on P1   %
            % and NURBS, are considered. Additionally, we have to evaluate also   %
            % the weights for both cases.                                         %
            %---------------------------------------------------------------------%

            mesh_xx    = zeros( ne*nqnx, 1);        % NODES CLASSIC LINEAR
            mesh_wx    = zeros( ne*nqnx, 1);        % WEIGHTS CLASSIC LINEAR

            mesh_xxIGA = zeros( (ncp-1)*nqnx, 1);   % NODES IGA
            mesh_wxIGA = zeros( (ncp-1)*nqnx, 1); 	% WEIGHTS IGA

            for i = 1:ne
                
                obj_quadratureRule = IntegrateHandler();
                
                obj_quadratureRule.leftBoundInterval = meshx(i);
                obj_quadratureRule.rightBoundInterval = meshx(i+1);
                obj_quadratureRule.inputNodes = xq;
                obj_quadratureRule.inputWeights = wxq;
                
                [mesh_xx((i-1)*nqnx+1 : i*nqnx), mesh_wx((i-1)*nqnx+1 : i*nqnx)] = ...
                                                 quadratureRule(obj_quadratureRule);

                %-----------------------------------------------------------------%
                % Note:                                                           %
                % We now have the mesh in the desired interval. However, it is    %
                % still necessary to rescale it in the most convenient way taking %
                % in consideration the shape of the domain.                       %
                %-----------------------------------------------------------------%

                xgauss = mesh_xx((i-1)*nqnx+1 : i*nqnx);

                for hp = 1: nqnx
                    
                    obj_gaussLegendre = IntegrateHandler();
                    obj_gaussLegendre.numbQuadNodes = 10;
                    
                    [~,gpos,gpes] = gaussLegendre(obj_gaussLegendre);
                    
                    gpos = xgauss(hp) * gpos;
                    gpes = xgauss(hp) * gpes;
                    xgauss(hp) = sum(sqrt(1+(obj.domainProfileDer(gpos)).^2).*gpes); 
                    
                end
                mesh_xx((i-1)*nqnx+1 : i*nqnx) = xgauss;
            end

            %---------------------------------------------------------------------%
            % Note:                                                               %
            % The vector 'mesh_xx' contains all the quadrature nodes divided into %
            % the respective elements. The nodes used are the classic nodes used  %
            % with linear elements.                                               %
            %---------------------------------------------------------------------%

            higa = 1/(ncp-1);
            meshIGA = zeros(1,ncp);
            
            for i = 2:ncp
                 meshIGA(i) = meshIGA(i-1) + higa;
            end


            for i = 1:ncp-1
                
                obj_quadratureRule = IntegrateHandler();
                
                obj_quadratureRule.leftBoundInterval = meshIGA(i);
                obj_quadratureRule.rightBoundInterval = meshIGA(i+1);
                obj_quadratureRule.inputNodes = xq;
                obj_quadratureRule.inputWeights = wxq;
                
                [mesh_xxIGA((i-1)*nqnx+1 : i*nqnx), mesh_wxIGA((i-1)*nqnx+1 : i*nqnx)] = ...
                                                        quadratureRule(obj_quadratureRule);
                                                    
            end

            size_fb  = ncp;

            % Computation of the Modal Basis in the Y Direction
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = intNodesGaussLeg;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis,~] = newModalBasis(obj_newModalBasis);

            % Block Matrices of the System with the Known Term

            b = zeros ( size_fb*(obj.dimModalBasis), 1);

            [x,~]          = meshgrid( mesh_xx, intNodesGaussLeg );
            L_computed     = obj.geometricInfo.L(x)';   % Thickness
            a_computed     = obj.geometricInfo.a(x)';   % Inferior Y Coordinate of the Channel

            y = a_computed';                          % Insertion of the Dimentionality

            [x,y]          = meshgrid( mesh_xx, intNodesGaussLeg-0.5 );

            mu_c    = obj.coefficientForm.mu(x,y)';
            beta1_c = obj.coefficientForm.beta1(x,y)';
            beta2_c = obj.coefficientForm.beta2(x,y)';
            sigma_c = obj.coefficientForm.sigma(x,y)';

            % Definition of the Exciting Force Acting on the System

            force_c = obj.dirCondFuncStruct.force(x,y,obj.timeInstant + 1)';

            Computed = struct('mu_c',mu_c,'beta1_c',beta1_c,'beta2_c',beta2_c, ...
                              'sigma_c',sigma_c,'force_c',force_c,'y',y);

            %---------------------------------------------------------------------%
            %                          ASSEMBLING LOOP                            %
            %---------------------------------------------------------------------%
            % Each frequency pair is considered. Then, the respective block of    %
            % the matrix and block of the known term are computed. In the         %
            % function both variables 'degree_p' and 'degree_k' are used.         %
            %---------------------------------------------------------------------%
            
            for imb = 1:obj.dimModalBasis

                [bmb,liftCoeffA,liftCoeffB] = assemblerIGAForce(intNodesWeights,modalBasis(:,imb),L_computed, ...
                                            obj.label_upBoundDomain,obj.label_downBoundDomain,obj.localdata_upBDomain,...
                                            obj.localdata_downBDomain,obj.coefficientForm,Computed,obj.stepMeshX,...
                                            obj.degreePolySplineBasis,obj.continuityParameter,obj.domainProfileDer,obj.leftBDomain_inX,...
                                            obj.rightBDomain_inX);

                b( 1+(imb-1)*size_fb : imb*size_fb ) = bmb;

            end
            end
        
            %% Method 'buildSystemIGASTR'
            
            function [A,b,modalBasis,liftCoeffA,liftCoeffB,verGLNodes,verGLWeights] = buildSystemIGASTR(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            %   (23) delta                      : Parameter to tune for
            %                                     stabilization
            %
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler

            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*obj.continuityParameter + obj.degreePolySplineBasis + 1 - obj.continuityParameter;
            
            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%

            % Horizontal direction
            
            obj_gaussLegendre_1 = IntegrateHandler();
            obj_gaussLegendre_1.numbQuadNodes = numbHorNodes;
            [~, horGLNodes, horGLWeights] = gaussLegendre(obj_gaussLegendre_1); 
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);        
            
            %% ISOGEOMETRIC MESH IN THE X DIRECTION
            %-------------------------------------------------------------%
            % Creation of the isogeometric mesh in the X direction
            % considering the control point assigned to the desired spline
            % curve. The mesh is created using the total number of control
            % point. The number of control points depend of the degree of
            % the spline basis, the continuity parameter 'k' and the number
            % of knots with which the physical mesh was divided.
            %-------------------------------------------------------------%
            
            stepKnot = (obj.rightBDomain_inX - obj.leftBDomain_inX)/(numbControlPts-1);
            meshIGA = zeros(1,numbControlPts);
            
            for i = 2:numbControlPts
                 meshIGA(i) = meshIGA(i-1) + stepKnot;
            end
            
            %% AUGMENTED ISOGEOMETRIC MESH + WEIGHTS 
            %-------------------------------------------------------------%
            % Creation of the finite element mesh in the X direction
            % considering equispaced nodes. The mesh is created using the
            % total number of nodes and the left limit of the domain. Both
            % information come from the demo file and are passed here as a
            % property of the object.
            %-------------------------------------------------------------%

            augMeshIGA = zeros( (numbControlPts-1)*numbHorNodes, 1);

            %-------------------------------------------------------------%
            % Note: The loop allocates the correponding quadrature nodes
            % and weights corresponding to each knot of the physical mesh
            % in the isogeometric analysis.
            %-------------------------------------------------------------%
             
            for i = 1:numbKnots
                
                % STEP 1
                %---------------------------------------------------------%
                % In the first step, the Gauss-Legendre nodes computed
                % previously are rescaled to fit the interval corresponding
                % to the current element.
                %---------------------------------------------------------%
                
                obj_quadratureRule = IntegrateHandler();
                
                obj_quadratureRule.leftBoundInterval = meshIGA(i);
                obj_quadratureRule.rightBoundInterval = meshIGA(i+1);
                obj_quadratureRule.inputNodes = horGLNodes;
                obj_quadratureRule.inputWeights = horGLWeights;
                
                [augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes), ~] = ...
                                                 quadratureRule(obj_quadratureRule);
             
                % STEP 2
                %---------------------------------------------------------%
                % In the second step, the nodes and weights are again
                % rescaled to take in consideration the geometry of the
                % domain, more specifically the profile of the centerline.
                %---------------------------------------------------------%

                auxMesh = augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes);

                for hp = 1: numbHorNodes
                    
                    obj_gaussLegendre = IntegrateHandler();
                    obj_gaussLegendre.numbQuadNodes = 16;
                    
                    [~,auxPoints,auxWeights] = gaussLegendre(obj_gaussLegendre);
                    
                    auxPoints = auxMesh(hp) * auxPoints;
                    auxWeights = auxMesh(hp) * auxWeights;
                    auxMesh(hp) = sum(sqrt(1+(obj.domainProfileDer(auxPoints)).^2).*auxWeights); 
                                       
                end
                augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes) = auxMesh;
                                                    
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Debug
            % figure;
            % plot(meshIGA);
            % figure;
            % plot(augMeshIGA);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% AUGMENTED VERTICAL POINTS (RESCALE)
            
            objVertQuadRule = IntegrateHandler();

            objVertQuadRule.leftBoundInterval = obj.downBDomain_inY;
            objVertQuadRule.rightBoundInterval = obj.upBDomain_inY;
            objVertQuadRule.inputNodes = verGLNodes;
            objVertQuadRule.inputWeights = verWeights;

            [augVerNodes, augVerWeights] = quadratureRule(objVertQuadRule);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Debug
            % figure;
            % plot(augVerNodes);
            % figure;
            % plot(verGLNodes);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = augVerNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%

            [horEvalNodes,verEvalNodes]   = meshgrid( augMeshIGA, augVerNodes );
            
            % THICKNESS
            
            evalL     = obj.geometricInfo.L(horEvalNodes)';
            
            %% Y COORDINATE OF THE CENTERLINE
            %-------------------------------------------------------------%
            % In the current code this value is not used to perform any
            % computation.
            %-------------------------------------------------------------%
            
            % a_computed     = obj.geometricInfo.a(horEvalNodes)';        
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal and vertical meshes to evaluate the
            % coefficients of the bilinear form of the original equation in
            % the entire domain.
            %-------------------------------------------------------------%

            % DIFFUSION
            
            evalMu    = obj.coefficientForm.mu(horEvalNodes,verEvalNodes)';
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(horEvalNodes,verEvalNodes)';
            evalBeta2 = obj.coefficientForm.beta2(horEvalNodes,verEvalNodes)';
            
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(horEvalNodes,verEvalNodes)';
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(horEvalNodes,verEvalNodes)';
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed = struct('mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2, ...
                              'sigma_c',evalSigma,'force_c',evalForce,'y',verEvalNodes);
                          
            % Debug
            %-------------------------------------------------------------%
            %[~,~] = contourf(horEvalNodes,verEvalNodes,evalForce',20);    
            %colormap(jet); title('Force in Assembler!')
            %-------------------------------------------------------------%
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            %% LOCAL PECHLET
            %---------------------------------------------------------------------%
            % Note: Computation of the local Pechlet number using the L2-norm and
            % the values of the convective field in each point of the domain.
            %---------------------------------------------------------------------%

            normBetaL2 = sqrt(evalBeta1.^2 + evalBeta2.^2);
    
            Pechlet = normBetaL2 .* obj.stepMeshX ./ (2*evalMu);
            
            % EVALUATION OF PECHLET
            
            figure;
            
            [~,~] = contourf(horEvalNodes,verEvalNodes,Pechlet',20);    
            colormap(jet);
            cmin = min(min(Pechlet));
            cmax = max(max(Pechlet));
            caxis([cmin cmax])
            colorbar();
            % axis([minX maxX minY maxY]);
            axis equal
            set(gca, 'FontSize', 14)
            
            % STABILIZED PECHLET
            
            stabPechlet = normBetaL2 .* obj.stepMeshX ./ (2 .* (1 + obj.delta .* Pechlet) .* evalMu);
            
            % EVALUATION OF STABILIZED PECHLET
            
            figure;
            
            [~,~] = contourf(horEvalNodes,verEvalNodes,stabPechlet',20);    
            colormap(jet);
            cmin = min(min(stabPechlet));
            cmax = max(max(stabPechlet));
            caxis([cmin cmax])
            colorbar();
            % axis([minX maxX minY maxY]);
            axis equal
            set(gca, 'FontSize', 14)
            
                          
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y)';
            
            %% JACOBIAN
            %-------------------------------------------------------------%
            % This section computes the curved domain and extract the
            % Jacobian vector that maps the centerline in the physical
            % domain to the computational domain.
            %-------------------------------------------------------------%
            
            % Mesh IGA of the curved domain
            
            meshIGACurved    = zeros(numbKnots + 1,1);
            [verGLNodes,verGLWeights] = gauss(numbVerNodes);
            
            verGLWeights = (obj.rightBDomain_inX - obj.leftBDomain_inX) * (obj.stepMeshX) * (verGLWeights + 1)*((obj.upBDomain_inY - obj.downBDomain_inY)/2);
                        
            for j=1:numbKnots
        
                scalingVec = zeros(1,numbVerNodes);
                
                % Pay attention to this code and what it does with the
                % vertical nodes. It needs modification.
                
                auxNodes = (obj.rightBDomain_inX - obj.leftBDomain_inX) * obj.stepMeshX * (verGLWeights + 1)*((obj.upBDomain_inY - obj.downBDomain_inY)/2) ...
                           + obj.leftBDomain_inX + (j-1) * obj.stepMeshX * (obj.rightBDomain_inX - obj.leftBDomain_inX);

                for i=1:numbVerNodes
                    scalingVec(i) = sqrt(1+(obj.domainProfileDer(auxNodes(i)))^2);
                end

                % Also pay attetion to verGLWeights here. It makes no
                % sense. May be the sorce of error when computing cases
                % with the deformed supporting fiber.
                
                meshIGACurved(j+1) = meshIGACurved(j) + sum(scalingVec * (verGLWeights + 1)*((obj.upBDomain_inY - obj.downBDomain_inY)/2));

            end
            
            % Jacobian loop
            
            Jac = [];
            refCoordinate = linspace(obj.leftBDomain_inX,obj.rightBDomain_inX,length(augMeshIGA));
            
            for ii = 1:length(augMeshIGA)
                
                Jac = [Jac sqrt(1+(obj.domainProfileDer(refCoordinate(ii)))^2)];
                
            end
   
%             for iel = 1:numbKnots
% 
%                 x_dx = meshIGACurved(iel);                         % Right Extreme
%                 x_sx = meshIGACurved(iel+1);                       % Left Extreme
% 
%                 point = (x_sx-x_dx)*(verGLWeights + 1)*((obj.upBDomain_inY - obj.downBDomain_inY)/2) + x_dx;
% 
%                 for j = 1 : numbVerNodes
%                     Jac = [Jac sqrt(1+(obj.domainProfileDer(point(j)))^2)];
%                 end
%             end
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis
                for kmb = 1:obj.dimModalBasis

                    [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerIGASTR( imb, kmb, numbControlPts, ...
                                            verGLWeights,modalBasis(:,imb), ...
                                            modalBasisDer(:,imb),modalBasis(:,kmb),...
                                            modalBasisDer(:,kmb), evalL, ...
                                            Computed,obj.stepMeshX,...
                                            obj.degreePolySplineBasis,obj.continuityParameter,obj.leftBDomain_inX,...
                                            obj.rightBDomain_inX,lifting,aLift,bLift,Jac,numbHorNodes,numbVerNodes,...
                                            numbKnots,horGLWeights,obj.delta,Pechlet);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;

            end

            end
            
            %% Method 'buildSystemIGARTS'
            
            function [A,b,modalBasis,liftCoeffA,liftCoeffB,verGLNodes,verGLWeights] = buildSystemIGARTS(obj)
                                                         
            %%
            % buildSystemIGA - This function computes the assembled matrices
            %                    relative to the variational problem considering 
            %                    IGA modal basis.
            %
            % Note:
            % All of the following inputs are encapsulated in the
            % object properties.
            %
            % The inputs are:
            %%
            %   (1)  dimModalBasis              : Dimension of the Modal Basis
            %   (2)  leftBDomain_inX            : Left Limit of the Domain in the X Direction
            %   (3)  rightBDomain_inX           : Right Limit of the Domain in the X Direction
            %   (4)  stepMeshX                  : Vector Containing the Step of the Finite
            %                                     Element Mesh
            %   (5)  label_upBoundDomain        : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Upper Limit of
            %                                     the Domain
            %   (6)  label_downBoundDomain      : Contains the Label Identifying the Nature of
            %                                     the Boundary Conditions on the Lower Limit of
            %                                     the Domain
            %   (7)  localdata_upBDomain        : Contains the Values of the Boundary Conditions
            %                                     on the Upper Limir of the Domain
            %   (8) localdata_downBDomain       : Contains the Values of the Boundary Conditions
            %                                     on the Lower Limir of the Domain
            %   (9) domainPosition              : Current domain used in the domain decomposition
            %                                     Computations
            %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                     and the Constants Relative to the Bilinear Form
            %   (11) dirCondFuncStruct          : Data Structure Containing All the @-Functions
            %                                     for the Dirichlet Conditions at the Inflow and
            %                                     for the Exciting Force
            %   (12) geometricInfo              : Data Structure Containing All the
            %                                     Geometric Information regarding the
            %                                     Domain. The current version of the code
            %                                     works only for the specific condition of:
            %                                     (L = 1, a = 0, psi_x = 0)
            %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
            %                                     Coefficients (R, L) for the Robin Condition Used
            %                                     in the Domain Decomposition
            %   (14) numbDimMBEachDomain        : Number of Elements in the Vector Containing the
            %                                     Dimensions of the Modal Basis in Each Domain
            %   (15) couplingCond_DD            : Contains the Label Adressing the Coupling Condition
            %                                     of the Problem
            %   (16) physicMesh_inX             : Vector Containing the Physical Mesh in the X
            %                                     Direction
            %   (17) physicMesh_inY             : Vector Containing the Physical Mesh in the Y
            %                                     Direction
            %   (18) jacAtQuadNodes             : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (19) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
            %   (20) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
            %   (21) domainProfile              : Symbolic Function Defining the Profile of the
            %                                     Simulation Domain
            %   (22) domainProfileDer           : Symbolic Function Defining the Derivative of
            %                                     the Profile of the Simulation Domain
            %   (23) delta                      : Parameter to tune for
            %                                     stabilization
            %
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) liftCoeffA          : First Offset Adjustment Coefficient
            %   (4) liftCoeffB          : Second Offset Adjustment Coefficient
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler

            %% NUMBER OF NODES IN THE QUADRATURE FORMULA
            %-------------------------------------------------------------%
            % The values of 'nqnx' and 'nqny' does not change during the
            % computation of the Gauss-Legendre Integration Points.
            % Attention, if you increase the number of nodes in the
            % discretization, it is necessary to refine the quadrature
            % rule.
            %-------------------------------------------------------------%

            % Horizontal direction
            
            numbHorNodes = obj.numbHorQuadNodes;
            
            % Vertical Direction
            
            numbVerNodes = obj.numbVerQuadNodes;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            numbKnots  = round((obj.rightBDomain_inX-obj.leftBDomain_inX)/obj.stepMeshX);
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve.
            %-------------------------------------------------------------%
            
            numbControlPts = numbKnots*obj.continuityParameter + obj.degreePolySplineBasis + 1 - obj.continuityParameter;
            
            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].
            %
            % Mesh FEM in X: Equispaced Nodes                            
            %  (1) horGLNodes  : Vector of the Standard Nodes on the X 
            %                    Direction
            %  (2) horGLWeights : Vector of the Standard Weights on the X 
            %                     Direction
            %
            % Mesh FEM in Y: Equispaced Nodes                              
            %  (1) verGLNodes  : Vector of the Standard Nodes on the Y 
            %                    Direction  
            %  (2) verGLWeights : Vector of the Standard Weights on the Y 
            %                     Direction
            %-------------------------------------------------------------%

            % Horizontal direction
            
            obj_gaussLegendre_1 = IntegrateHandler();
            obj_gaussLegendre_1.numbQuadNodes = numbHorNodes;
            [~, horGLNodes, horGLWeights] = gaussLegendre(obj_gaussLegendre_1); 
            
            % Vertical direction
            
            obj_gaussLegendre_2 = IntegrateHandler();
            obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
            [~, verGLNodes, ~] = gaussLegendre(obj_gaussLegendre_2);             
            
            %% ISOGEOMETRIC MESH IN THE X DIRECTION
            %-------------------------------------------------------------%
            % Creation of the isogeometric mesh in the X direction
            % considering the control point assigned to the desired spline
            % curve. The mesh is created using the total number of control
            % point. The number of control points depend of the degree of
            % the spline basis, the continuity parameter 'k' and the number
            % of knots with which the physical mesh was divided.
            %-------------------------------------------------------------%
            
            stepKnot = 1/(numbControlPts-1);
            meshIGA = zeros(1,numbControlPts);
            
            for i = 2:numbControlPts
                 meshIGA(i) = meshIGA(i-1) + stepKnot;
            end
            
            %% AUGMENTED ISOGEOMETRIC MESH + WEIGHTS 
            %-------------------------------------------------------------%
            % Creation of the finite element mesh in the X direction
            % considering equispaced nodes. The mesh is created using the
            % total number of nodes and the left limit of the domain. Both
            % information come from the demo file and are passed here as a
            % property of the object.
            %-------------------------------------------------------------%

            augMeshIGA = zeros( (numbControlPts-1)*numbHorNodes, 1);

            %-------------------------------------------------------------%
            % Note: The loop allocates the correponding quadrature nodes
            % and weights corresponding to each knot of the physical mesh
            % in the isogeometric analysis.
            %-------------------------------------------------------------%
             
            for i = 1:numbKnots
                
                % STEP 1
                %---------------------------------------------------------%
                % In the first step, the Gauss-Legendre nodes computed
                % previously are rescaled to fit the interval corresponding
                % to the current element.
                %---------------------------------------------------------%
                
                obj_quadratureRule = IntegrateHandler();
                
                obj_quadratureRule.leftBoundInterval = meshIGA(i);
                obj_quadratureRule.rightBoundInterval = meshIGA(i+1);
                obj_quadratureRule.inputNodes = horGLNodes;
                obj_quadratureRule.inputWeights = horGLWeights;
                
                [augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes), ~] = ...
                                                 quadratureRule(obj_quadratureRule);
             
                % STEP 2
                %---------------------------------------------------------%
                % In the second step, the nodes and weights are again
                % rescaled to take in consideration the geometry of the
                % domain, more specifically the profile of the centerline.
                %---------------------------------------------------------%

                auxMesh = augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes);

                for hp = 1: numbHorNodes
                    
                    obj_gaussLegendre = IntegrateHandler();
                    obj_gaussLegendre.numbQuadNodes = 16;
                    
                    [~,auxPoints,auxWeights] = gaussLegendre(obj_gaussLegendre);
                    
                    auxPoints = auxMesh(hp) * auxPoints;
                    auxWeights = auxMesh(hp) * auxWeights;
                    auxMesh(hp) = sum(sqrt(1+(obj.domainProfileDer(auxPoints)).^2).*auxWeights); 
                                       
                end
                augMeshIGA((i-1)*numbHorNodes+1 : i*numbHorNodes) = auxMesh;
                                                    
            end
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasis, modalBasisDer] = newModalBasis(obj_newModalBasis);
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( numbControlPts*(obj.dimModalBasis), numbControlPts*(obj.dimModalBasis) );
            b = zeros ( numbControlPts*(obj.dimModalBasis), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%

            [horEvalNodes,verEvalNodes]   = meshgrid( augMeshIGA, verGLNodes-0.5 );%!!!!!!
            
            % THICKNESS
            
            evalL     = obj.geometricInfo.L(horEvalNodes)';
            
            %% Y COORDINATE OF THE CENTERLINE
            %-------------------------------------------------------------%
            % In the current code this value is not used to perform any
            % computation.
            %-------------------------------------------------------------%
            
            % a_computed     = obj.geometricInfo.a(horEvalNodes)';        
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal and vertical meshes to evaluate the
            % coefficients of the bilinear form of the original equation in
            % the entire domain.
            %-------------------------------------------------------------%

            % DIFFUSION
            
            evalMu    = obj.coefficientForm.mu(horEvalNodes,verEvalNodes)';
            
            % ADVECTION
            
            evalBeta1 = obj.coefficientForm.beta1(horEvalNodes,verEvalNodes)';
            evalBeta2 = obj.coefficientForm.beta2(horEvalNodes,verEvalNodes)';
            
            % REACTION
            
            evalSigma = obj.coefficientForm.sigma(horEvalNodes,verEvalNodes)';
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForce = obj.dirCondFuncStruct.force(horEvalNodes,verEvalNodes)';
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            Computed = struct('mu_c',evalMu,'beta1_c',evalBeta1,'beta2_c',evalBeta2, ...
                              'sigma_c',evalSigma,'force_c',evalForce,'y',verEvalNodes);
                          
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomain;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomain;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLift,bLift] = liftBoundCond(obj_liftBoundCond);
            liftFunc = @(x,y) aLift * y + bLift;
            lifting = liftFunc(0,Computed.y)';
            
            %% JACOBIAN
            %-------------------------------------------------------------%
            % This section computes the curved domain and extract the
            % Jacobian vector that maps the centerline in the physical
            % domain to the computational domain.
            %-------------------------------------------------------------%
            
            % Mesh IGA of the curved domain
            
            meshIGACurved    = zeros(numbKnots + 1,1);
            [verGLNodes,verGLWeights] = gauss(numbVerNodes);
            
            verGLWeights = (obj.rightBDomain_inX - obj.leftBDomain_inX) * (obj.stepMeshX) * (verGLWeights*0.5);
                        
            for j=1:numbKnots
        
                scalingVec = zeros(1,numbVerNodes);
                auxNodes = (obj.rightBDomain_inX - obj.leftBDomain_inX) * obj.stepMeshX * (verGLNodes * 0.5 + 0.5) ...
                           + obj.leftBDomain_inX + (j-1) * obj.stepMeshX * (obj.rightBDomain_inX - obj.leftBDomain_inX);

                for i=1:numbVerNodes
                    scalingVec(i) = sqrt(1+(obj.domainProfileDer(auxNodes(i)))^2);
                end

                meshIGACurved(j+1) = meshIGACurved(j) + sum(scalingVec * verGLWeights);

            end
            
            % Jacobian loop
            
            Jac = [];
   
            for iel = 1:numbKnots

                x_dx = meshIGACurved(iel);                         % Right Extreme
                x_sx = meshIGACurved(iel+1);                       % Left Extreme

                point = (x_sx-x_dx)*(0.5*verGLNodes + 0.5) + x_dx;

                for j = 1 : numbVerNodes
                    Jac = [Jac sqrt(1+(obj.domainProfileDer(point(j)))^2)];
                end
            end
            
            %% ASSEMBLING LOOP
            %-------------------------------------------------------------%
            % The assembling loop creates each one of the submatrices
            % correponding to the microstructure of the linear system. The
            % loop is repeated m^2 times to complete the macrostructure of
            % the system.
            %-------------------------------------------------------------%

            for imb = 1:obj.dimModalBasis
                for kmb = 1:obj.dimModalBasis

                    [Amb,bmb,liftCoeffA,liftCoeffB] = assemblerIGARTS( imb, kmb, numbControlPts, ...
                                            verGLWeights,modalBasis(:,imb), ...
                                            modalBasisDer(:,imb),modalBasis(:,kmb),...
                                            modalBasisDer(:,kmb), evalL, ...
                                            Computed,obj.stepMeshX,...
                                            obj.degreePolySplineBasis,obj.continuityParameter,obj.leftBDomain_inX,...
                                            obj.rightBDomain_inX,lifting,aLift,bLift,Jac,numbHorNodes,numbVerNodes,...
                                            numbKnots,horGLWeights,obj.delta);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPts : imb*numbControlPts , 1+(kmb-1)*numbControlPts : kmb*numbControlPts) = Amb;

                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPts : imb*numbControlPts ) = bmb;

            end

            end

    end
end

%% ASSEMBLER ADR HANDLER - ASSEMBLING METHODS
% The following functions correspond to the individual scripts
% previously developed of the functions 'assembla_x',
% 'assembla_legendre' and 'assembla_x_IGA'.

%% Method 'assemblerFEM'
            
function [Al,bl,aLift,bLift] = assemblerFEM(imb,kmb,numbNodes, ...
                                verGLWeights,mb_i,mb_yi,mb_k,mb_yk, ...
                                L, Computed, ...
                                horStep,p,k, ...
                                lifting,aLift,bLift,jacFEM,...
                                numbHorQuadNodes,~,...
                                numbElements,horWeights,basisFEM, ...
                                basisFEMDer,numbModes)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% JACOBIAN IMPOSITIONS
    %---------------------------------------------------------------------%
    % This values will come from the demo when the code is changed to
    % incorporate also cases where L(x) is not constant. Note that the name
    % of the variables were chosen to match the same name in the reference
    % paper used to create this function.
    %---------------------------------------------------------------------%
    
    D1 = 0;
    D2 = 1;
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%

    objForceInt = IntegrateHandler();
    
    objForceInt.funcToIntegrate =   L.*Computed.force_c - ...
                                    L.*(aLift).*Computed.beta2_c - ...
                                    Computed.sigma_c.*L.*lifting;
                                
    objForceInt.funcWeight = mb_i.*verGLWeights;
    
    forceVec  = integrate(objForceInt);

    %---------------------------------------------------------------------%
    %                          DISCONTINUOUS FORCES                       %
    %---------------------------------------------------------------------%
    % Note: The following script can be used when the force exciting the
    %       system in discontinuous in the domain.
    %
    % forc = integrate(  - L.*a_ril.*Computed.beta2_c - ... 
    %                   Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);
    % nqny=64;
    % [~, adhocnodes, adhocw ] = gaussLegendre( nqny );
    % adhocnodes=adhocnodes/5+0.4;
    % adhocw=adhocw/5;
    % mb_iaus = newModalBasis( imb, adhocnodes, bc_up,bc_down,Coeff_forma);
    % forc = forc + integrate( L.*Computed.force_c, mb_iaus(:,imb).*adhocw);
    %---------------------------------------------------------------------%

    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
        
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    obj_integrate_3 = IntegrateHandler();
    obj_integrate_4 = IntegrateHandler();
    obj_integrate_5 = IntegrateHandler();
    obj_integrate_6 = IntegrateHandler();
    obj_integrate_7 = IntegrateHandler();
    obj_integrate_8 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_2.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_3.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_4.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_5.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_6.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_7.funcToIntegrate = L.*Computed.beta2_c;
    obj_integrate_8.funcToIntegrate = L.*Computed.sigma_c;
    
    obj_integrate_1.funcWeight = mb_k  .*mb_i  .*verGLWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_yi .*verGLWeights;
    obj_integrate_3.funcWeight = mb_k  .*mb_i  .*verGLWeights;
    obj_integrate_4.funcWeight = mb_yk .*mb_i  .*verGLWeights;
    obj_integrate_5.funcWeight = mb_yk .*mb_yi .*verGLWeights;
    obj_integrate_6.funcWeight = mb_yk .*mb_i  .*verGLWeights;
    obj_integrate_7.funcWeight = mb_yk .*mb_i  .*verGLWeights;
    obj_integrate_8.funcWeight = mb_k  .*mb_i  .*verGLWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    lambda_3   = integrate(obj_integrate_3);
    lambda_4   = integrate(obj_integrate_4);
    lambda_5   = integrate(obj_integrate_5);
    lambda_6   = integrate(obj_integrate_6);
    lambda_7   = integrate(obj_integrate_7);
    lambda_8   = integrate(obj_integrate_8);
    
    r00 = lambda_5 * (D1^2 + D2^2) + lambda_6 * D1 + lambda_7 * D2 + lambda_8;
    r10 = lambda_2 * D1 + lambda_3;
    r01 = lambda_4 * D1;
    r11 = lambda_1;

    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. This version of assemblerFEM works only for basis of
    % degree equal to 1.
    %---------------------------------------------------------------------%
    bl       = zeros( numbElements + 1, 1);
    maindiag = zeros( numbElements + 1, 1);
    diagsup  = zeros( numbElements + 1, 1); % Spdiags taglier?? il primo elemento
    diaginf  = zeros( numbElements + 1, 1); % Spdiags invece taglier?? l'ultimo

    for ie = 1 : numbElements % Loop Runs in Each Interval

        % Computation of the Weights in the Given Interval

        %-----------------------------------------------------------------%
        % Note: The index inside the function 'mesh_wx()' selects the
        %       correct points interval to evaluate the weights.
        %-----------------------------------------------------------------%

        localWeight = horWeights((ie - 1) * numbHorQuadNodes + 1 : ie * numbHorQuadNodes)';

        % Selection of the Coefficients Relative to the Given Interval

        r11Local = r11 ( (ie-1) * numbHorQuadNodes + 1 : ie * numbHorQuadNodes)';
        r01Local = r01 ( (ie-1) * numbHorQuadNodes + 1 : ie * numbHorQuadNodes)';
        r10Local = r10 ( (ie-1) * numbHorQuadNodes + 1 : ie * numbHorQuadNodes)';
        r00Local = r00 ( (ie-1) * numbHorQuadNodes + 1 : ie * numbHorQuadNodes)';

        % Loop Runs the Following Computations for each Function in the 
        % Shape of EF

        for i = 1:p+1

            % Selection of the Correct Pieces for the Basis in the FEM

            shape_i  = basisFEM    ( (ie-1) * numbHorQuadNodes + 1:ie * numbHorQuadNodes, i)';
            dshape_i = basisFEMDer ( (ie-1) * numbHorQuadNodes + 1:ie * numbHorQuadNodes, i)';

            % Loop Runs Again These Computations for each Function in the 
            % Shape of EF

            for j = 1:p+1

                shape_j  = basisFEM    ( (ie-1) * numbHorQuadNodes + 1:ie * numbHorQuadNodes, j)';
                dshape_j = basisFEMDer ( (ie-1) * numbHorQuadNodes + 1:ie * numbHorQuadNodes, j)';

                if(i==j)

                    maindiag(ie+i-1)=maindiag(ie+i-1)+...
                        (  r11Local.*dshape_i.*dshape_j ...
                        + r01Local.*dshape_i.*shape_j...
                        + r10Local.*shape_i.*dshape_j ...
                        + r00Local.*shape_i.*shape_j  ...
                        )*localWeight;

                elseif(i>j)
                    diaginf(ie)=diaginf(ie)+(...
                        r11Local.*dshape_i.*dshape_j ...
                        + r01Local.*dshape_i.*shape_j...
                        + r10Local.*shape_i.*dshape_j ...
                        + r00Local.*shape_i.*shape_j  ...
                        )*localWeight;
                else
                    diagsup(ie+1)=diagsup(ie+1)+( ...
                        r11Local.*dshape_i.*dshape_j ...
                        + r01Local.*dshape_i.*shape_j...
                        + r10Local.*shape_i.*dshape_j ...
                        + r00Local.*shape_i.*shape_j  ...
                        )*localWeight;
                end

                %--------------------------------------------------------------------------------------------------------------%
                % Note: This compact version is more readable but less
                %       effective even considering the 30% computational
                %       cost reduction during the assembling operation.
                % 
                % Code Example: Al(ie+i-1,ie+j-1) =  Al(ie+i-1,ie+j-1)+integrate( m.*femb_xi.*femb_xj+b.*femb_i.*femb_xj ...
                %                                                                 + g.*femb_xi.*femb_j+s.*femb_i.*femb_j   ...
                %                                                                 + funzionale_robin_up.*femb_i.*femb_j    ...
                %                                                                 + funzionale_robin_down.*femb_i.*femb_j  ...
                %                                                                 , w);
                %--------------------------------------------------------------------------------------------------------------%

            end
            
            bl(ie+i-1) = bl(ie+i-1) + ( forceVec((ie-1)*numbHorQuadNodes+1:ie*numbHorQuadNodes)'.*shape_i )*(localWeight );
            
        end
    end
    
    % Assignment of the Diagonal Elements of Matrix A

    Al = spdiags( [diaginf,maindiag,diagsup], -1:1, numbNodes, numbNodes);

    disp('Finished ASSEMBLING LOOP');

end
            
%% Method 'assemblerIGA'
            
function [Al,bl,aLift,bLift] = assemblerIGA(imb,kmb,numbControlPts, ...
                                augVerWeights,mb_i,mb_yi,mb_k,mb_yk, ...
                                L, Computed, ...
                                horStep,p,k, ...
                                domainLeftLimit,domainRightLimit, ...
                                lifting,aLift,bLift,jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights,D1,D2)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%

    objForceInt = IntegrateHandler();
    
    objForceInt.funcToIntegrate =   L.*Computed.force_c - ...
                                    L.*(aLift).*Computed.beta2_c - ...
                                    Computed.sigma_c.*L.*lifting;
                                
    objForceInt.funcWeight = mb_i.*augVerWeights;
    
    forceVec  = integrate(objForceInt);
    
%     % Debug
%     %---------------------------------------------------------------------%
%     
%     disp(size(forceVec))
%     disp(size(Computed.force_c))
%     disp(size(mb_i))
    
%     disp(max(max(lifting)))
%     disp(min(min(lifting)))
%     disp(max(max(aLift)))
%     disp(min(min(aLift)))
%     disp(max(max(bLift)))
%     disp(min(min(bLift)))
    %---------------------------------------------------------------------%
    
    % Debug 
    %---------------------------------------------------------------------%
    %disp(['Max L : ',num2str(max(max(L)))]);
    %disp(['Min L : ',num2str(min(min(L)))]);
    %plot(linspace(0,length(forceVec),length(forceVec)),forceVec);
    %---------------------------------------------------------------------%

    %---------------------------------------------------------------------%
    %                          DISCONTINUOUS FORCES                       %
    %---------------------------------------------------------------------%
    % Note: The following script can be used when the force exciting the
    %       system in discontinuous in the domain.
    %
    % forc = integrate(  - L.*a_ril.*Computed.beta2_c - ... 
    %                   Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);
    % nqny=64;
    % [~, adhocnodes, adhocw ] = gaussLegendre( nqny );
    % adhocnodes=adhocnodes/5+0.4;
    % adhocw=adhocw/5;
    % mb_iaus = newModalBasis( imb, adhocnodes, bc_up,bc_down,Coeff_forma);
    % forc = forc + integrate( L.*Computed.force_c, mb_iaus(:,imb).*adhocw);
    %---------------------------------------------------------------------%
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    obj_integrate_3 = IntegrateHandler();
    obj_integrate_4 = IntegrateHandler();
    obj_integrate_5 = IntegrateHandler();
    obj_integrate_6 = IntegrateHandler();
    obj_integrate_7 = IntegrateHandler();
    obj_integrate_8 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_2.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_3.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_4.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_5.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_6.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_7.funcToIntegrate = L.*Computed.beta2_c;
    obj_integrate_8.funcToIntegrate = L.*Computed.sigma_c;
    
    obj_integrate_1.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_yi .*augVerWeights;
    obj_integrate_3.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    obj_integrate_4.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_5.funcWeight = mb_yk .*mb_yi .*augVerWeights;
    obj_integrate_6.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_7.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_8.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    lambda_3   = integrate(obj_integrate_3);
    lambda_4   = integrate(obj_integrate_4);
    lambda_5   = integrate(obj_integrate_5);
    lambda_6   = integrate(obj_integrate_6);
    lambda_7   = integrate(obj_integrate_7);
    lambda_8   = integrate(obj_integrate_8);
    
    r00 = lambda_5 * (D1^2 + D2^2) + lambda_6 * D1 + lambda_7 * D2 + lambda_8;
    r10 = lambda_2 * D1 + lambda_3;
    r01 = lambda_4 * D1;
    r11 = lambda_1;

    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = p;
    objBasisIGA.continuityParameter   = k;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGA,derBasisIGA,controlPts] = newIsoGeoBasis(objBasisIGA);  
    
    %% LOCAL JACOBIAN
    
%     Jac = [];
%     meshx_p1 = linspace(0,1,nel+1);
% 
%     for iel = 1:nel
%         
%         x_dx = meshx_p1(iel);                         % Right Extreme
%         x_sx = meshx_p1(iel+1);                       % Left Extreme
%         [point,~] = gauss(nqnx);                      % Creation of the 'p+1' Quadrature Points
%         
%         point = (x_sx-x_dx)*(0.5*point + 0.5) + x_dx;
%         
%         for j = 1 : nqnx
%            Jac = [Jac sqrt(1+(dforma(point(j)))^2)];
%         end
%     end
    
    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    globalForce  = zeros(numbControlPts,1);
    row   = zeros(1,numbControlPts^2);
    col   = zeros(1,numbControlPts^2);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        forceElem  = zeros(p+1,1);
        
        deltaElem1 = zeros(p+1,p+1);
        deltaElem2 = zeros(p+1,p+1);
        deltaElem3 = zeros(p+1,p+1);
        deltaElem4 = zeros(p+1,p+1);
        
        % INTEGRAL OVER THE SUPPORTING FIBER
        %-----------------------------------------------------------------%
        % We compute the coefficients referring to the integrals along the
        % supporting fiber.
        %-----------------------------------------------------------------%
        
        for igauss = 1:numbHorQuadNodes

            % SHAPE FUNCTIONS
            %-------------------------------------------------------------%
            % Compute the shape functions and their derivatives with
            % respect to the parametric variable, evaluated over the
            % Gauss-Legendre integration points.
            %-------------------------------------------------------------%

            shape  = basisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            dshape = derBasisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            cpi    = controlPts((k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            
            % COORDINATE OF THE GAUSS POINT IN THE PHYSICAL DOMAIN
            %-------------------------------------------------------------%
            % This componponent is not being used now because the thickness
            % of the channel is constant and equal to '1'. However, it is
            % extremely necessary for more complex geometries because it
            % gives the correct horizontal component in the physical domain
            % to evaluate the function L(x).
            %-------------------------------------------------------------%
            
            x      = shape*cpi'*jacIGA((iel-1)*numbHorQuadNodes+igauss);
            
            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshape*cpi';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

%             shape = shape/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));
            dshape = dshape/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;
            
%             gwt = horWeights(igauss);
            
            % LOCAL RHS VECTOR
            %-------------------------------------------------------------%
            % The current code uses 'numbHorQuadNodes' instead of
            % 'numbVerQuadNodes' in the next 5 lines of code. Check if
            % the solution is consistent.
            %-------------------------------------------------------------%

            forceElem       = forceElem  + shape'*forceVec((iel-1)*numbHorQuadNodes+igauss)*gwt;    % Force
            
            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shape' * shape) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (shape' *dshape) * gwt;
            deltaElem3 = deltaElem3 + r01((iel-1)*numbHorQuadNodes+igauss) * (dshape'* shape) * gwt;
            deltaElem4 = deltaElem4 + r11((iel-1)*numbHorQuadNodes+igauss) * (dshape'*dshape) * gwt;
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:p+1
            
            i = (k-1)*(iel-1) + iel + p - (a - 1);

            if (i>1)

                globalForce(i) = globalForce(i) + forceElem(a);

                for b = 1:p+1
                    
                    j = (k-1)*(iel-1) + iel + p - (b - 1);
                    icount = icount + 1;
                    
                    row(icount)   = i;
                    col(icount)   = j;
                    
                    deltaLocal1(icount) = deltaElem1(a,b);
                    deltaLocal2(icount) = deltaElem2(a,b);
                    deltaLocal3(icount) = deltaElem3(a,b);
                    deltaLocal4(icount) = deltaElem4(a,b);
                    
                end
            else 
                icount = icount + 1;
                row(icount)  = i;
                col(icount)  = i;
                
                deltaLocal1(icount) = 0;
                deltaLocal2(icount) = 0;
                deltaLocal3(icount) = 0;
                deltaLocal4(icount) = 1;
                
            end
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPts,numbControlPts);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPts,numbControlPts);
    Delta3 = sparse(row,col,deltaLocal3,numbControlPts,numbControlPts);
    Delta4 = sparse(row,col,deltaLocal4,numbControlPts,numbControlPts);
    
    Al   = Delta1 + Delta2 + Delta3 + Delta4;
    bl   = globalForce;
    
    % Debug
    %---------------------------------------------------------------------%
    % figure; plot(bl);
    %---------------------------------------------------------------------%
    
    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%
    
%     if (imb == kmb)
%         Al(1,1)  = 1;
%         Al(1,2)  = 0;
%         Al(1,3)  = 0;
%     else
%         Al(1,1)  = 1e-15;
%         Al(1,2)  = 0;
%         Al(1,3)  = 0;
%     end
    
%     if (imb == kmb)
%         Al(end,end)  = 1;
%         Al(end,end-1)  = 0;
%         Al(end,end-2)  = 0;
%     else
%         Al(end,end)  = 1e-15;
%         Al(end,end-1)  = 0;
%         Al(end,end-2)  = 0;
%     end
    
    if (imb == kmb)
        Al(1,1)  = 1;
        Al(1,2)  = 0;
        Al(1,3)  = 0;
    else
        Al(1,1)  = 1e-15;
        Al(1,2)  = 0;
        Al(1,3)  = 0;
    end
    
    bl(1) = 0;
%     bl(end) = 0;
end

%% Method 'assemblerIGAScatter'
            
function [Al,bl,aLift,bLift] = assemblerIGAScatter(imb,kmb,augVerWeights,mb_i,mb_yi,mb_k,mb_yk,geoData,...
                                Computed,lifting,aLift,bLift,msh,space,jacFunc,spaceFunc,BoundCond)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%
    
    funcToIntegrale = (Computed.force_c - ...
                      (aLift).* Computed.beta2_c - ...
                      Computed.sigma_c .* lifting) .* ...
                      jacFunc.evalDetJac;
                                
    funcWeight = mb_i .* augVerWeights;
    
    forceVec  = sum(funcToIntegrale .* funcWeight , 1);
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;
    alpha2      = jacFunc.Phi2_dx.^2 + jacFunc.Phi2_dy.^2;
    varbeta1    = Computed.beta1_c .* jacFunc.Phi1_dx + ...
                  Computed.beta2_c .* jacFunc.Phi1_dy;
    varbeta2    = Computed.beta1_c .* jacFunc.Phi2_dx + ...
                  Computed.beta2_c .* jacFunc.Phi2_dy;
    delta       = jacFunc.Phi1_dx .* jacFunc.Phi2_dx + ...
                  jacFunc.Phi1_dy .* jacFunc.Phi2_dy;
    
    funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* alpha1;
    funcToIntegrate_2 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_3 = jacFunc.evalDetJac .* varbeta1;
    funcToIntegrate_4 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_5 = jacFunc.evalDetJac .* Computed.mu_c .* alpha2;
    funcToIntegrate_6 = jacFunc.evalDetJac .* varbeta2;
    funcToIntegrate_7 = jacFunc.evalDetJac .* Computed.sigma_c;
    
    funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_2 = mb_k  .* mb_yi .* augVerWeights;
    funcWeight_3 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_4 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_5 = mb_yk .* mb_yi .* augVerWeights;
    funcWeight_6 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_7 = mb_k  .* mb_i  .* augVerWeights;
    
    aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);
    aux2   = sum(funcToIntegrate_2 .* funcWeight_2 , 1);
    aux3   = sum(funcToIntegrate_3 .* funcWeight_3 , 1);
    aux4   = sum(funcToIntegrate_4 .* funcWeight_4 , 1);
    aux5   = sum(funcToIntegrate_5 .* funcWeight_5 , 1);
    aux6   = sum(funcToIntegrate_6 .* funcWeight_6 , 1);
    aux7   = sum(funcToIntegrate_7 .* funcWeight_7 , 1);
    
    r00 = aux5 + aux6 + aux7;
    r10 = aux2 + aux3;
    r01 = aux4;
    r11 = aux1;
    
    coeff.r00 = r00;
    coeff.r10 = r10;
    coeff.r01 = r01;
    coeff.r11 = r11;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = msh.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_10 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_01 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_11 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    rhs      = zeros (space.ndof, 1);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00 = Local_00 + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, r00Local);
        Local_10 = Local_10 + op_gradu_v(spaceFunc{iel,2}, spaceFunc{iel,1}, spaceFunc{iel,3}, r10Local);
        Local_01 = Local_01 + op_u_gradv(spaceFunc{iel,1}, spaceFunc{iel,2}, spaceFunc{iel,3}, r01Local);
        Local_11 = Local_11 + op_gradu_gradv(spaceFunc{iel,2}, spaceFunc{iel,2}, spaceFunc{iel,3}, r11Local);
        
        rhs      = rhs + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, fLocal);

    end

    Al       = Local_11 + Local_10 + Local_01 + Local_00;
    bl       = rhs;

    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%
    
%     [Al,bl] = imposeBoundIGA(Al,bl,BoundCond,space,msh,coeff,geoData,spaceFunc,imb,kmb);
% 
%         % IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX
%         
%         if (imb == kmb)
%             AModified = Al;
%             AModified(1,:) = 0;
%             AModified(1,1) = 1e10;
%             AModified(end,:) = 0;
%             AModified(end,end) = 1e10;
%         else
%             AModified = Al;
%             AModified(1,:) = 0;
%             AModified(end,:) = 0;
%         end
%         
%         Al = AModified;
%         
%         % IMPLEMENT MODIFICATIONS IN THE FORCING TERM
%         
%         if (imb == kmb)
%             bModified = bl;
%             bModified(1) = 0;
%             bModified(end) = 0;
%         else
%             bModified = bl;
%         end   
%         
%         bl = bModified;
%     
%     % NEUMANN BOUNDARY CONDITIONS
%     % Note: Modify the orinial stiffness matrix and right-handside term to
%     % include the Neumann boundary condtions, if any.
%     
%     if(imb == kmb)
%         for iside = BoundCond.neuSides
%             
%             if (iside == 1)
%               x = msh.breaks{1}(1);
%             else
%               x = msh.breaks{1}(end);
%             end
%             
%             sp_side = space.boundary(iside);
%             bl(sp_side.dofs) = bl(sp_side.dofs) + BoundCond.neu(x,iside);
%             
%         end
%     end  
    
end

%% Method 'assemblerIGAScatterWeak'
            
function [Al,bl,aLift,bLift] = assemblerIGAScatterWeak(imb,kmb,augVerWeights,mb_i,mb_yi,mb_k,mb_yk,geoData,...
                                Computed,lifting,aLift,bLift,msh,space,jacFunc,spaceFunc,BoundCond)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% PENALTY TERMS
    %---------------------------------------------------------------------%
    % Here we computes the penalty contributions divided in up/down and
    % inflow/outflow contributions
    %---------------------------------------------------------------------%
    
    %penalty up and down contributions
    
    [r00Up, penalty_forceUp] = penalty_upBC(imb, kmb, mb_i(end), mb_yi(end), mb_k(end), mb_yk(end), jacFunc.Phi2_dy(end,:), LateralConditions);
    [r00Down, penalty_forceDown] = penalty_downBC(imb, kmb, mb_i(1), mb_yi(1), mb_k(1), mb_yk(1), jacFunc.Phi2_dy(1,:), LateralConditions);
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%
    
    funcToIntegrale = (Computed.force_c - ...
                      (aLift).* Computed.beta2_c - ...
                      Computed.sigma_c .* lifting) .* ...
                      jacFunc.evalDetJac;
                                
    funcWeight = mb_i .* augVerWeights;
    
    forceVec  = sum(funcToIntegrale .* funcWeight , 1);
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;
    alpha2      = jacFunc.Phi2_dx.^2 + jacFunc.Phi2_dy.^2;
    varbeta1    = Computed.beta1_c .* jacFunc.Phi1_dx + ...
                  Computed.beta2_c .* jacFunc.Phi1_dy;
    varbeta2    = Computed.beta1_c .* jacFunc.Phi2_dx + ...
                  Computed.beta2_c .* jacFunc.Phi2_dy;
    delta       = jacFunc.Phi1_dx .* jacFunc.Phi2_dx + ...
                  jacFunc.Phi1_dy .* jacFunc.Phi2_dy;
    
    funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* alpha1;
    funcToIntegrate_2 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_3 = jacFunc.evalDetJac .* varbeta1;
    funcToIntegrate_4 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_5 = jacFunc.evalDetJac .* Computed.mu_c .* alpha2;
    funcToIntegrate_6 = jacFunc.evalDetJac .* varbeta2;
    funcToIntegrate_7 = jacFunc.evalDetJac .* Computed.sigma_c;
    
    funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_2 = mb_k  .* mb_yi .* augVerWeights;
    funcWeight_3 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_4 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_5 = mb_yk .* mb_yi .* augVerWeights;
    funcWeight_6 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_7 = mb_k  .* mb_i  .* augVerWeights;
    
    aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);
    aux2   = sum(funcToIntegrate_2 .* funcWeight_2 , 1);
    aux3   = sum(funcToIntegrate_3 .* funcWeight_3 , 1);
    aux4   = sum(funcToIntegrate_4 .* funcWeight_4 , 1);
    aux5   = sum(funcToIntegrate_5 .* funcWeight_5 , 1);
    aux6   = sum(funcToIntegrate_6 .* funcWeight_6 , 1);
    aux7   = sum(funcToIntegrate_7 .* funcWeight_7 , 1);
    
    r00 = aux5 + aux6 + aux7;
    r10 = aux2 + aux3;
    r01 = aux4;
    r11 = aux1;
    
    coeff.r00 = r00;
    coeff.r10 = r10;
    coeff.r01 = r01;
    coeff.r11 = r11;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = msh.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_10 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_01 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_11 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    rhs      = zeros (space.ndof, 1);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00 = Local_00 + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, r00Local);
        Local_10 = Local_10 + op_gradu_v(spaceFunc{iel,2}, spaceFunc{iel,1}, spaceFunc{iel,3}, r10Local);
        Local_01 = Local_01 + op_u_gradv(spaceFunc{iel,1}, spaceFunc{iel,2}, spaceFunc{iel,3}, r01Local);
        Local_11 = Local_11 + op_gradu_gradv(spaceFunc{iel,2}, spaceFunc{iel,2}, spaceFunc{iel,3}, r11Local);
        
        rhs      = rhs + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, fLocal);

    end

    
    Al       = Local_11 + Local_10 + Local_01 + Local_00;
    bl       = rhs;
    
end

%% Method 'assemblerIGAScatterTransient'
            
function [Al,Ml,bl] = assemblerIGAScatterTransient(imb,kmb,augVerWeights,mb_i,mb_yi,mb_k,mb_yk,geoData,...
                                Computed,msh,space,jacFunc,spaceFunc)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    warning('off','all')
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%
    
    funcToIntegrale = Computed.force_c .* jacFunc.evalDetJac;
                                
    funcWeight = mb_i .* augVerWeights;
    
    forceVec  = sum(funcToIntegrale .* funcWeight , 1);
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;
    alpha2      = jacFunc.Phi2_dx.^2 + jacFunc.Phi2_dy.^2;
    varbeta1    = Computed.beta1_c .* jacFunc.Phi1_dx + ...
                  Computed.beta2_c .* jacFunc.Phi1_dy;
    varbeta2    = Computed.beta1_c .* jacFunc.Phi2_dx + ...
                  Computed.beta2_c .* jacFunc.Phi2_dy;
    delta       = jacFunc.Phi1_dx .* jacFunc.Phi2_dx + ...
                  jacFunc.Phi1_dy .* jacFunc.Phi2_dy;
    
    funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* alpha1;
    funcToIntegrate_2 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_3 = jacFunc.evalDetJac .* varbeta1;
    funcToIntegrate_4 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_5 = jacFunc.evalDetJac .* Computed.mu_c .* alpha2;
    funcToIntegrate_6 = jacFunc.evalDetJac .* varbeta2;
    funcToIntegrate_7 = jacFunc.evalDetJac .* Computed.sigma_c;
    funcToIntegrate_8 = jacFunc.evalDetJac .* Computed.mass_c;
    
    funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_2 = mb_k  .* mb_yi .* augVerWeights;
    funcWeight_3 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_4 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_5 = mb_yk .* mb_yi .* augVerWeights;
    funcWeight_6 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_7 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_8 = mb_k  .* mb_i  .* augVerWeights;
    
    aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);
    aux2   = sum(funcToIntegrate_2 .* funcWeight_2 , 1);
    aux3   = sum(funcToIntegrate_3 .* funcWeight_3 , 1);
    aux4   = sum(funcToIntegrate_4 .* funcWeight_4 , 1);
    aux5   = sum(funcToIntegrate_5 .* funcWeight_5 , 1);
    aux6   = sum(funcToIntegrate_6 .* funcWeight_6 , 1);
    aux7   = sum(funcToIntegrate_7 .* funcWeight_7 , 1);
    aux8   = sum(funcToIntegrate_8 .* funcWeight_8 , 1);
    
    r00 = aux5 + aux6 + aux7;
    r10 = aux2 + aux3;
    r01 = aux4;
    r11 = aux1;
    m00 = aux8;
    
    coeff.r00 = r00;
    coeff.r10 = r10;
    coeff.r01 = r01;
    coeff.r11 = r11;
    coeff.m00 = m00;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = msh.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Mass_00 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_00 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_10 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_01 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_11 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    rhs      = zeros (space.ndof, 1);
    
    for iel = 1:numbKnots
        
        m00Local = m00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Mass_00  = Mass_00 + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, m00Local);
        Local_00 = Local_00 + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, r00Local);
        Local_10 = Local_10 + op_gradu_v(spaceFunc{iel,2}, spaceFunc{iel,1}, spaceFunc{iel,3}, r10Local);
        Local_01 = Local_01 + op_u_gradv(spaceFunc{iel,1}, spaceFunc{iel,2}, spaceFunc{iel,3}, r01Local);
        Local_11 = Local_11 + op_gradu_gradv(spaceFunc{iel,2}, spaceFunc{iel,2}, spaceFunc{iel,3}, r11Local);
        
        rhs      = rhs + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, fLocal);

    end

    Ml = Mass_00;
    Al = Local_11 + Local_10 + Local_01 + Local_00;
    bl = rhs;

%     % IMPOSITION OF DIRICHLET BOUNDARY CONDITION
%     ---------------------------------------------------------------------%
%     This step is necessary to guarantee that the resulting linear system
%     is not singular.
%     ---------------------------------------------------------------------%
%     
%     [Al,bl] = imposeBoundIGA(Al,bl,BoundCond,space,msh,coeff,geoData,spaceFunc,imb,kmb);
% 
%         IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX
%         
%         if (imb == kmb)
%             AModified = Al;
%             AModified(1,:) = 0;
%             AModified(1,1) = 1e10;
%             AModified(end,:) = 0;
%             AModified(end,end) = 1e10;
%         else
%             AModified = Al;
%             AModified(1,:) = 0;
%             AModified(end,:) = 0;
%         end
%         
%         Al = AModified;
%         
%         IMPLEMENT MODIFICATIONS IN THE FORCING TERM
%         
%         if (imb == kmb)
%             bModified = bl;
%             bModified(1) = 0;
%             bModified(end) = 0;
%         else
%             bModified = bl;
%         end   
%         
%         bl = bModified;
%     
%     NEUMANN BOUNDARY CONDITIONS
%     Note: Modify the orinial stiffness matrix and right-handside term to
%     include the Neumann boundary condtions, if any.
%     
%     if(imb == kmb)
%         for iside = BoundCond.neuSides
%             
%             if (iside == 1)
%               x = msh.breaks{1}(1);
%             else
%               x = msh.breaks{1}(end);
%             end
%             
%             sp_side = space.boundary(iside);
%             bl(sp_side.dofs) = bl(sp_side.dofs) + BoundCond.neu(x,iside);
%             
%         end
%     end  
    
end

%% Method 'assemblerIGAScatterStabilized'
            
function [Al,bl,aLift,bLift] = assemblerIGAScatterStabilized(imb,kmb,augVerWeights,mb_i,mb_yi,mb_k,mb_yk,geoData,...
                                Computed,lifting,aLift,bLift,msh,space,jacFunc,spaceFunc,BoundCond,stabMethod,stabDelta,PeConfig)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%
    
    funcToIntegrale = (Computed.force_c - ...
                      (aLift).* Computed.beta2_c - ...
                      Computed.sigma_c .* lifting) .* ...
                      jacFunc.evalDetJac;
                                
    funcWeight = mb_i .* augVerWeights;
    
    forceVec  = sum(funcToIntegrale .* funcWeight , 1);
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;
    alpha2      = jacFunc.Phi2_dx.^2 + jacFunc.Phi2_dy.^2;
    varbeta1    = Computed.beta1_c .* jacFunc.Phi1_dx + ...
                  Computed.beta2_c .* jacFunc.Phi1_dy;
    varbeta2    = Computed.beta1_c .* jacFunc.Phi2_dx + ...
                  Computed.beta2_c .* jacFunc.Phi2_dy;
    delta       = jacFunc.Phi1_dx .* jacFunc.Phi2_dx + ...
                  jacFunc.Phi1_dy .* jacFunc.Phi2_dy;
    
    funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* alpha1;
    funcToIntegrate_2 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_3 = jacFunc.evalDetJac .* varbeta1;
    funcToIntegrate_4 = jacFunc.evalDetJac .* Computed.mu_c .* delta;
    funcToIntegrate_5 = jacFunc.evalDetJac .* Computed.mu_c .* alpha2;
    funcToIntegrate_6 = jacFunc.evalDetJac .* varbeta2;
    funcToIntegrate_7 = jacFunc.evalDetJac .* Computed.sigma_c;
    
    funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_2 = mb_k  .* mb_yi .* augVerWeights;
    funcWeight_3 = mb_k  .* mb_i  .* augVerWeights;
    funcWeight_4 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_5 = mb_yk .* mb_yi .* augVerWeights;
    funcWeight_6 = mb_yk .* mb_i  .* augVerWeights;
    funcWeight_7 = mb_k  .* mb_i  .* augVerWeights;
    
    aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);
    aux2   = sum(funcToIntegrate_2 .* funcWeight_2 , 1);
    aux3   = sum(funcToIntegrate_3 .* funcWeight_3 , 1);
    aux4   = sum(funcToIntegrate_4 .* funcWeight_4 , 1);
    aux5   = sum(funcToIntegrate_5 .* funcWeight_5 , 1);
    aux6   = sum(funcToIntegrate_6 .* funcWeight_6 , 1);
    aux7   = sum(funcToIntegrate_7 .* funcWeight_7 , 1);
    
    r00 = aux5 + aux6 + aux7;
    r10 = aux2 + aux3;
    r01 = aux4;
    r11 = aux1;
    
    %% COMPUTE THE PECHLET NUMBER
    
    switch PeConfig
                
        case 1

            % Compute the step of the mesh used in the supporting
            % fiber
            
            xMax = max(max(geoData.X));
            h = sqrt((xMax/geoData.numbInt)^2 + (1/geoData.numbModes)^2);

            % Compute the Pechlet matrix
            
            beta = sqrt(Computed.beta1_c.^2 + Computed.beta2_c.^2);
            Pechlet = (h/2) * abs(beta)./abs(Computed.mu_c);

        case 2

            % Compute the step of the mesh used in the supporting
            % fiber
            
            h = 1/geoData.numbInt;

            % Compute the Pechlet matrix

            Pechlet = (h/2) * abs(r10)./abs(r11);

    end
        
    %% COMPUTE THE STABILIZATION COEFFICIENTS
    
    switch stabMethod
        
    case{'USTR'}
        
        % Compute the step of the mesh used in the supporting
        % fiber

        xMax = max(max(geoData.X));
        h = sqrt((xMax/geoData.numbInt)^2 + (1/geoData.numbModes)^2);

        % Compute the Pechlet matrix

        beta = sqrt(Computed.beta1_c.^2 + Computed.beta2_c.^2);
        Pechlet = (h/2) * abs(beta)./abs(Computed.mu_c);

        alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;
        alpha2      = jacFunc.Phi2_dx.^2 + jacFunc.Phi2_dy.^2;
        delta       = jacFunc.Phi1_dx .* jacFunc.Phi2_dx + ...
                      jacFunc.Phi1_dy .* jacFunc.Phi2_dy;

        funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* Pechlet .* stabDelta .* alpha1;
        funcToIntegrate_2 = jacFunc.evalDetJac .* Computed.mu_c .* Pechlet .* stabDelta .* delta;
        funcToIntegrate_3 = jacFunc.evalDetJac .* Computed.mu_c .* Pechlet .* stabDelta .* delta;
        funcToIntegrate_4 = jacFunc.evalDetJac .* Computed.mu_c .* Pechlet .* stabDelta .* alpha2;

        funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;
        funcWeight_2 = mb_k  .* mb_yi .* augVerWeights;
        funcWeight_3 = mb_yk .* mb_i  .* augVerWeights;
        funcWeight_4 = mb_yk .* mb_yi .* augVerWeights;

        aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);
        aux2   = sum(funcToIntegrate_2 .* funcWeight_2 , 1);
        aux3   = sum(funcToIntegrate_3 .* funcWeight_3 , 1);
        aux4   = sum(funcToIntegrate_4 .* funcWeight_4 , 1);

        s00 = aux4;
        s10 = aux2;
        s01 = aux3;
        s11 = aux1;
        
    case {'URTS'}
        
        % Compute the step of the mesh used in the supporting
        % fiber

        h = 1/geoData.numbInt;

        % Compute the Pechlet matrix

        Pechlet = (h/2) * abs(r10)./abs(r11);
        
        alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;

        funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* alpha1;

        funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;

        aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);

        s00 = 0;
        s10 = 0;
        s01 = 0;
        s11 = aux1 .* Pechlet .* stabDelta;
        
    case {'SDSTR'}
        
        % Compute the step of the mesh used in the supporting
        % fiber

        xMax = max(max(geoData.X));
        h = sqrt((xMax/geoData.numbInt)^2 + (1/geoData.numbModes)^2);

        % Compute the Pechlet matrix

        beta = sqrt(Computed.beta1_c.^2 + Computed.beta2_c.^2);
        Pechlet = (h/2) * abs(beta)./abs(Computed.mu_c);
        
        varbeta1    = Computed.beta1_c .* jacFunc.Phi1_dx + ...
                      Computed.beta2_c .* jacFunc.Phi1_dy;
        varbeta2    = Computed.beta1_c .* jacFunc.Phi2_dx + ...
                      Computed.beta2_c .* jacFunc.Phi2_dy;
                  
        PHIsd = Pechlet./sqrt(Computed.beta1_c.^2 + Computed.beta2_c.^2);
        
        PHIsd(1)
        
        funcToIntegrate_1 = jacFunc.evalDetJac .* PHIsd .* varbeta1.^2;
        funcToIntegrate_2 = jacFunc.evalDetJac .* PHIsd .* varbeta1 .* varbeta2;
        funcToIntegrate_3 = jacFunc.evalDetJac .* PHIsd .* varbeta1 .* varbeta2;
        funcToIntegrate_4 = jacFunc.evalDetJac .* PHIsd .* varbeta2.^2;

        funcWeight_1 = mb_k   .* mb_i   .* augVerWeights;
        funcWeight_2 = mb_k   .* mb_yi  .* augVerWeights;
        funcWeight_3 = mb_yk  .* mb_i   .* augVerWeights;
        funcWeight_4 = mb_yk  .* mb_yi  .* augVerWeights;
        
        aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);
        aux2   = sum(funcToIntegrate_2 .* funcWeight_2 , 1);
        aux3   = sum(funcToIntegrate_3 .* funcWeight_3 , 1);
        aux4   = sum(funcToIntegrate_4 .* funcWeight_4 , 1);
        
        s00 = stabDelta .* aux4;
        s10 = stabDelta .* aux2;
        s01 = stabDelta .* aux3;
        s11 = stabDelta .* aux1;
        
    case {'SDRTS'}
        
        % Compute the step of the mesh used in the supporting
        % fiber

        h = 1/geoData.numbInt;

        % Compute the Pechlet matrix

        Pechlet = (h/2) * abs(r10)./abs(r11);
        
        PHIsd = (Pechlet.^2)./r10;
        
        s00 = 0;
        s10 = 0;
        s01 = 0;
        s11 = stabDelta .* PHIsd;
        
        PHIsd(1)
        
    case {'SGSTR'}
        
        % Compute the step of the mesh used in the supporting
        % fiber

        xMax = max(max(geoData.X));
        h = sqrt((xMax/geoData.numbInt)^2 + (1/geoData.numbModes)^2);

        % Compute the Pechlet matrix

        beta = sqrt(Computed.beta1_c.^2 + Computed.beta2_c.^2);
        Pechlet = (h/2) * abs(beta)./abs(Computed.mu_c);
        
        Bern = @(x) 1.*(x == 0) + (x./(exp(x) - 1)).*(x>0);
        PHIsg = @(x) x - 1 + Bern(2 .* x);
        
        alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;
        alpha2      = jacFunc.Phi2_dx.^2 + jacFunc.Phi2_dy.^2;
        delta       = jacFunc.Phi1_dx .* jacFunc.Phi2_dx + ...
                      jacFunc.Phi1_dy .* jacFunc.Phi2_dy;

        funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* PHIsg(Pechlet) .* stabDelta .* alpha1;
        funcToIntegrate_2 = jacFunc.evalDetJac .* Computed.mu_c .* PHIsg(Pechlet) .* stabDelta .* delta;
        funcToIntegrate_3 = jacFunc.evalDetJac .* Computed.mu_c .* PHIsg(Pechlet) .* stabDelta .* delta;
        funcToIntegrate_4 = jacFunc.evalDetJac .* Computed.mu_c .* PHIsg(Pechlet) .* stabDelta .* alpha2;

        funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;
        funcWeight_2 = mb_k  .* mb_yi .* augVerWeights;
        funcWeight_3 = mb_yk .* mb_i  .* augVerWeights;
        funcWeight_4 = mb_yk .* mb_yi .* augVerWeights;

        aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);
        aux2   = sum(funcToIntegrate_2 .* funcWeight_2 , 1);
        aux3   = sum(funcToIntegrate_3 .* funcWeight_3 , 1);
        aux4   = sum(funcToIntegrate_4 .* funcWeight_4 , 1);

        s00 = aux4;
        s10 = aux2;
        s01 = aux3;
        s11 = aux1;
        
    case {'SGRTS'}
        
        % Compute the step of the mesh used in the supporting
        % fiber

        h = 1/geoData.numbInt;

        % Compute the Pechlet matrix

        Pechlet = (h/2) * abs(r10)./abs(r11);
        
        Bern = @(x) 1.*(x == 0) + (x./(exp(x) - 1)).*(x>0);
        PHIsg = @(x) x - 1 + Bern(2 .* x);
        
        alpha1      = jacFunc.Phi1_dx.^2 + jacFunc.Phi1_dy.^2;

        funcToIntegrate_1 = jacFunc.evalDetJac .* Computed.mu_c .* alpha1;

        funcWeight_1 = mb_k  .* mb_i  .* augVerWeights;

        aux1   = sum(funcToIntegrate_1 .* funcWeight_1 , 1);

        s00 = 0;
        s10 = 0;
        s01 = 0;
        s11 = aux1 .* PHIsg(Pechlet) .* stabDelta;
    
    case {'None'}
            
        s00 = 0;
        s10 = 0;
        s01 = 0;
        s11 = 0;
        
    end
    
    r00 = r00 + s00;
    r10 = r10 + s10;
    r01 = r01 + s01;
    r11 = r11 + s11;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = msh.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_10 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_01 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_11 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    rhs      = zeros (space.ndof, 1);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00 = Local_00 + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, r00Local);
        Local_10 = Local_10 + op_gradu_v(spaceFunc{iel,2}, spaceFunc{iel,1}, spaceFunc{iel,3}, r10Local);
        Local_01 = Local_01 + op_u_gradv(spaceFunc{iel,1}, spaceFunc{iel,2}, spaceFunc{iel,3}, r01Local);
        Local_11 = Local_11 + op_gradu_gradv(spaceFunc{iel,2}, spaceFunc{iel,2}, spaceFunc{iel,3}, r11Local);
        
        rhs      = rhs + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, fLocal);

    end

    
    Al       = Local_11 + Local_10 + Local_01 + Local_00;
    bl       = rhs;

    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%
    
%     [Al,bl] = imposeBoundIGA(Al,bl,BoundCond,space,msh,coeff,geoData,spaceFunc,imb,kmb);

        % IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX
        
        if (imb == kmb)
            AModified = Al;
            AModified(1,:) = 0;
            AModified(1,1) = 1e10;
            AModified(end,:) = 0;
            AModified(end,end) = 1e10;
        else
            AModified = Al;
            AModified(1,:) = 0;
            AModified(end,:) = 0;
        end
        
        Al = AModified;
        
        % IMPLEMENT MODIFICATIONS IN THE FORCING TERM
        
        if (imb == kmb)
            bModified = bl;
            bModified(1) = 0;
            bModified(end) = 0;
        else
            bModified = bl;
        end   
        
        bl = bModified;
    
    % NEUMANN BOUNDARY CONDITIONS
    % Note: Modify the orinial stiffness matrix and right-handside term to
    % include the Neumann boundary condtions, if any.
    
    if(imb == kmb)
        for iside = BoundCond.neuSides
            
            if (iside == 1)
              x = msh.breaks{1}(1);
            else
              x = msh.breaks{1}(end);
            end
            
            sp_side = space.boundary(iside);
            bl(sp_side.dofs) = bl(sp_side.dofs) + BoundCond.neu(x,iside);
            
        end
    end  
    
end

%% Method 'assemblerIGAScatter3D'
            
function [Al,bl,aLift,bLift] = assemblerIGAScatter3D(imb,kmb,augVerWeights,mb_i,mb_yi,mb_zi,mb_k,mb_yk,mb_zk,geoData,...
                                Computed,lifting,aLift,bLift,msh,space,jacFunc,spaceFunc,BoundCond)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%
    
    % funcToIntegrate = (Computed.force_c - ...
    %                   (aLift).* Computed.beta2_c - ...
    %                   Computed.sigma_c .* lifting) .* ...
    %                   jacFunc.evalDetJac;

    funcToIntegrate = jacFunc.evalDetJac .* Computed.force_c;
    
    weightMat = augVerWeights * augVerWeights';
    funcWeight = mb_i .* weightMat;

    aux = funcToIntegrate .* funcWeight;
    auxx = sum(sum(aux));
    forceVec(:) = auxx(1,1,:);

    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    Psi11 = jacFunc.Psi1_dx.^2 + jacFunc.Psi1_dy.^2 + jacFunc.Psi1_dz.^2;
    Psi22 = jacFunc.Psi2_dx.^2 + jacFunc.Psi2_dy.^2 + jacFunc.Psi2_dz.^2;
    Psi33 = jacFunc.Psi3_dx.^2 + jacFunc.Psi3_dy.^2 + jacFunc.Psi3_dz.^2;
    
    Psi12 = jacFunc.Psi1_dx .* jacFunc.Psi2_dx + ...
            jacFunc.Psi1_dy .* jacFunc.Psi2_dy + ...
            jacFunc.Psi1_dz .* jacFunc.Psi2_dz;
    Psi21 = jacFunc.Psi2_dx .* jacFunc.Psi1_dx + ...
            jacFunc.Psi2_dy .* jacFunc.Psi1_dy + ...
            jacFunc.Psi2_dz .* jacFunc.Psi1_dz;
    Psi13 = jacFunc.Psi1_dx .* jacFunc.Psi3_dx + ...
            jacFunc.Psi1_dy .* jacFunc.Psi3_dy + ...
            jacFunc.Psi1_dz .* jacFunc.Psi3_dz;
    Psi31 = jacFunc.Psi3_dx .* jacFunc.Psi1_dx + ...
            jacFunc.Psi3_dy .* jacFunc.Psi1_dy + ...
            jacFunc.Psi3_dz .* jacFunc.Psi1_dz;
    Psi23 = jacFunc.Psi2_dx .* jacFunc.Psi3_dx + ...
            jacFunc.Psi2_dy .* jacFunc.Psi3_dy + ...
            jacFunc.Psi2_dz .* jacFunc.Psi3_dz;
    Psi32 = jacFunc.Psi3_dx .* jacFunc.Psi2_dx + ...
            jacFunc.Psi3_dy .* jacFunc.Psi2_dy + ...
            jacFunc.Psi3_dz .* jacFunc.Psi2_dz;
    
    BetaPsi1 = Computed.beta1_c .* jacFunc.Psi1_dx + ...
               Computed.beta2_c .* jacFunc.Psi1_dy + ...
               Computed.beta3_c .* jacFunc.Psi1_dz;
    BetaPsi2 = Computed.beta1_c .* jacFunc.Psi2_dx + ...
               Computed.beta2_c .* jacFunc.Psi2_dy + ...
               Computed.beta3_c .* jacFunc.Psi2_dz;
    BetaPsi3 = Computed.beta1_c .* jacFunc.Psi3_dx + ...
               Computed.beta2_c .* jacFunc.Psi3_dy + ...
               Computed.beta3_c .* jacFunc.Psi3_dz;

    funcToIntegrate_1  = jacFunc.evalDetJac .* Computed.mu_c .* Psi11;
    funcToIntegrate_2  = jacFunc.evalDetJac .* Computed.mu_c .* Psi12;
    funcToIntegrate_3  = jacFunc.evalDetJac .* Computed.mu_c .* Psi13;
    funcToIntegrate_4  = jacFunc.evalDetJac .* Computed.mu_c .* Psi21;
    funcToIntegrate_5  = jacFunc.evalDetJac .* Computed.mu_c .* Psi22;
    funcToIntegrate_6  = jacFunc.evalDetJac .* Computed.mu_c .* Psi23;
    funcToIntegrate_7  = jacFunc.evalDetJac .* Computed.mu_c .* Psi31;
    funcToIntegrate_8  = jacFunc.evalDetJac .* Computed.mu_c .* Psi32;
    funcToIntegrate_9  = jacFunc.evalDetJac .* Computed.mu_c .* Psi33;
    funcToIntegrate_10 = jacFunc.evalDetJac .* BetaPsi1;
    funcToIntegrate_11 = jacFunc.evalDetJac .* BetaPsi2;
    funcToIntegrate_12 = jacFunc.evalDetJac .* BetaPsi3;
    
    % DEBUG
    
    
    
    funcWeight_1  = mb_k  .* mb_i  .* weightMat;
    funcWeight_2  = mb_k  .* mb_yi .* weightMat;
    funcWeight_3  = mb_k  .* mb_zi .* weightMat;
    funcWeight_4  = mb_yk .* mb_i  .* weightMat;
    funcWeight_5  = mb_yk .* mb_yi .* weightMat;
    funcWeight_6  = mb_yk .* mb_zi .* weightMat;
    funcWeight_7  = mb_zk .* mb_i  .* weightMat;
    funcWeight_8  = mb_zk .* mb_yi .* weightMat;
    funcWeight_9  = mb_zk .* mb_zi .* weightMat;
    funcWeight_10 = mb_k  .* mb_i  .* weightMat;
    funcWeight_11 = mb_yk .* mb_i  .* weightMat;
    funcWeight_12 = mb_zk .* mb_i  .* weightMat;

    aux1  = funcToIntegrate_1  .* funcWeight_1;
    aux2  = funcToIntegrate_2  .* funcWeight_2;
    aux3  = funcToIntegrate_3  .* funcWeight_3;
    aux4  = funcToIntegrate_4  .* funcWeight_4;
    aux5  = funcToIntegrate_5  .* funcWeight_5;
    aux6  = funcToIntegrate_6  .* funcWeight_6;
    aux7  = funcToIntegrate_7  .* funcWeight_7;
    aux8  = funcToIntegrate_8  .* funcWeight_8;
    aux9  = funcToIntegrate_9  .* funcWeight_9;
    aux10 = funcToIntegrate_10 .* funcWeight_10;
    aux11 = funcToIntegrate_11 .* funcWeight_11;
    aux12 = funcToIntegrate_12 .* funcWeight_12;

    auxVec1  = sum(sum(aux1));
    auxVec2  = sum(sum(aux2));
    auxVec3  = sum(sum(aux3));
    auxVec4  = sum(sum(aux4));
    auxVec5  = sum(sum(aux5));
    auxVec6  = sum(sum(aux6));
    auxVec7  = sum(sum(aux7));
    auxVec8  = sum(sum(aux8));
    auxVec9  = sum(sum(aux9));
    auxVec10 = sum(sum(aux10));
    auxVec11 = sum(sum(aux11));
    auxVec12 = sum(sum(aux12));
    
    r00 = auxVec5 + auxVec6 + auxVec8 + auxVec9 + auxVec11 + auxVec12;
    r10 = auxVec2 + auxVec3 + auxVec10;
    r01 = auxVec4 + auxVec7;
    r11 = auxVec1;
    
    auxx1(:) = r00(1,1,:);
    auxx2(:) = r01(1,1,:);
    auxx3(:) = r10(1,1,:);
    auxx4(:) = r11(1,1,:);
    r00 = auxx1;
    r01 = auxx2;
    r10 = auxx3;
    r11 = auxx4;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = msh.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_10 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_01 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_11 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    rhs      = zeros (space.ndof, 1);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00 = Local_00 + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, r00Local);
        Local_10 = Local_10 + op_gradu_v(spaceFunc{iel,2}, spaceFunc{iel,1}, spaceFunc{iel,3}, r10Local);
        Local_01 = Local_01 + op_u_gradv(spaceFunc{iel,1}, spaceFunc{iel,2}, spaceFunc{iel,3}, r01Local);
        Local_11 = Local_11 + op_gradu_gradv(spaceFunc{iel,2}, spaceFunc{iel,2}, spaceFunc{iel,3}, r11Local);
        
        rhs      = rhs + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, fLocal);

    end
    
    Al       = Local_11 + Local_10 + Local_01 + Local_00;
    bl       = rhs;

end

%% Method 'assemblerIGAScatter3DTransient'
            
function [Al,Ml,bl,aLift,bLift] = assemblerIGAScatter3DTransient(imb,kmb,augVerWeights,mb_i,mb_yi,mb_zi,mb_k,mb_yk,mb_zk,geoData,...
                                Computed,lifting,aLift,bLift,msh,space,jacFunc,spaceFunc,BoundCond)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%
    
    % funcToIntegrate = (Computed.force_c - ...
    %                   (aLift).* Computed.beta2_c - ...
    %                   Computed.sigma_c .* lifting) .* ...
    %                   jacFunc.evalDetJac;

    funcToIntegrate = jacFunc.evalDetJac .* Computed.force_c;
    
    weightMat = augVerWeights * augVerWeights';
    funcWeight = mb_i .* weightMat;

    aux = funcToIntegrate .* funcWeight;
    auxx = sum(sum(aux));
    forceVec(:) = auxx(1,1,:);

    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    Psi11 = jacFunc.Psi1_dx.^2 + jacFunc.Psi1_dy.^2 + jacFunc.Psi1_dz.^2;
    Psi22 = jacFunc.Psi2_dx.^2 + jacFunc.Psi2_dy.^2 + jacFunc.Psi2_dz.^2;
    Psi33 = jacFunc.Psi3_dx.^2 + jacFunc.Psi3_dy.^2 + jacFunc.Psi3_dz.^2;
    
    Psi12 = jacFunc.Psi1_dx .* jacFunc.Psi2_dx + ...
            jacFunc.Psi1_dy .* jacFunc.Psi2_dy + ...
            jacFunc.Psi1_dz .* jacFunc.Psi2_dz;
    Psi21 = jacFunc.Psi2_dx .* jacFunc.Psi1_dx + ...
            jacFunc.Psi2_dy .* jacFunc.Psi1_dy + ...
            jacFunc.Psi2_dz .* jacFunc.Psi1_dz;
    Psi13 = jacFunc.Psi1_dx .* jacFunc.Psi3_dx + ...
            jacFunc.Psi1_dy .* jacFunc.Psi3_dy + ...
            jacFunc.Psi1_dz .* jacFunc.Psi3_dz;
    Psi31 = jacFunc.Psi3_dx .* jacFunc.Psi1_dx + ...
            jacFunc.Psi3_dy .* jacFunc.Psi1_dy + ...
            jacFunc.Psi3_dz .* jacFunc.Psi1_dz;
    Psi23 = jacFunc.Psi2_dx .* jacFunc.Psi3_dx + ...
            jacFunc.Psi2_dy .* jacFunc.Psi3_dy + ...
            jacFunc.Psi2_dz .* jacFunc.Psi3_dz;
    Psi32 = jacFunc.Psi3_dx .* jacFunc.Psi2_dx + ...
            jacFunc.Psi3_dy .* jacFunc.Psi2_dy + ...
            jacFunc.Psi3_dz .* jacFunc.Psi2_dz;
    
    BetaPsi1 = Computed.beta1_c .* jacFunc.Psi1_dx + ...
               Computed.beta2_c .* jacFunc.Psi1_dy + ...
               Computed.beta3_c .* jacFunc.Psi1_dz;
    BetaPsi2 = Computed.beta1_c .* jacFunc.Psi2_dx + ...
               Computed.beta2_c .* jacFunc.Psi2_dy + ...
               Computed.beta3_c .* jacFunc.Psi2_dz;
    BetaPsi3 = Computed.beta1_c .* jacFunc.Psi3_dx + ...
               Computed.beta2_c .* jacFunc.Psi3_dy + ...
               Computed.beta3_c .* jacFunc.Psi3_dz;
           
    Sigma = Computed.sigma_c;
    Mass  = ones(size(Sigma));

    funcToIntegrate_1  = jacFunc.evalDetJac .* Computed.mu_c .* Psi11;
    funcToIntegrate_2  = jacFunc.evalDetJac .* Computed.mu_c .* Psi12;
    funcToIntegrate_3  = jacFunc.evalDetJac .* Computed.mu_c .* Psi13;
    funcToIntegrate_4  = jacFunc.evalDetJac .* Computed.mu_c .* Psi21;
    funcToIntegrate_5  = jacFunc.evalDetJac .* Computed.mu_c .* Psi22;
    funcToIntegrate_6  = jacFunc.evalDetJac .* Computed.mu_c .* Psi23;
    funcToIntegrate_7  = jacFunc.evalDetJac .* Computed.mu_c .* Psi31;
    funcToIntegrate_8  = jacFunc.evalDetJac .* Computed.mu_c .* Psi32;
    funcToIntegrate_9  = jacFunc.evalDetJac .* Computed.mu_c .* Psi33;
    funcToIntegrate_10 = jacFunc.evalDetJac .* BetaPsi1;
    funcToIntegrate_11 = jacFunc.evalDetJac .* BetaPsi2;
    funcToIntegrate_12 = jacFunc.evalDetJac .* BetaPsi3;
    funcToIntegrate_13 = jacFunc.evalDetJac .* Sigma;
    funcToIntegrate_14 = jacFunc.evalDetJac .* Mass;
    
    funcWeight_1  = mb_k  .* mb_i  .* weightMat;
    funcWeight_2  = mb_k  .* mb_yi .* weightMat;
    funcWeight_3  = mb_k  .* mb_zi .* weightMat;
    funcWeight_4  = mb_yk .* mb_i  .* weightMat;
    funcWeight_5  = mb_yk .* mb_yi .* weightMat;
    funcWeight_6  = mb_yk .* mb_zi .* weightMat;
    funcWeight_7  = mb_zk .* mb_i  .* weightMat;
    funcWeight_8  = mb_zk .* mb_yi .* weightMat;
    funcWeight_9  = mb_zk .* mb_zi .* weightMat;
    funcWeight_10 = mb_k  .* mb_i  .* weightMat;
    funcWeight_11 = mb_yk .* mb_i  .* weightMat;
    funcWeight_12 = mb_zk .* mb_i  .* weightMat;
    funcWeight_13 = mb_k  .* mb_i  .* weightMat;
    funcWeight_14 = mb_k  .* mb_i  .* weightMat;

    aux1  = funcToIntegrate_1  .* funcWeight_1;
    aux2  = funcToIntegrate_2  .* funcWeight_2;
    aux3  = funcToIntegrate_3  .* funcWeight_3;
    aux4  = funcToIntegrate_4  .* funcWeight_4;
    aux5  = funcToIntegrate_5  .* funcWeight_5;
    aux6  = funcToIntegrate_6  .* funcWeight_6;
    aux7  = funcToIntegrate_7  .* funcWeight_7;
    aux8  = funcToIntegrate_8  .* funcWeight_8;
    aux9  = funcToIntegrate_9  .* funcWeight_9;
    aux10 = funcToIntegrate_10 .* funcWeight_10;
    aux11 = funcToIntegrate_11 .* funcWeight_11;
    aux12 = funcToIntegrate_12 .* funcWeight_12;
    aux13 = funcToIntegrate_13 .* funcWeight_13;
    aux14 = funcToIntegrate_14 .* funcWeight_14;
    
    auxVec1  = sum(sum(aux1));
    auxVec2  = sum(sum(aux2));
    auxVec3  = sum(sum(aux3));
    auxVec4  = sum(sum(aux4));
    auxVec5  = sum(sum(aux5));
    auxVec6  = sum(sum(aux6));
    auxVec7  = sum(sum(aux7));
    auxVec8  = sum(sum(aux8));
    auxVec9  = sum(sum(aux9));
    auxVec10 = sum(sum(aux10));
    auxVec11 = sum(sum(aux11));
    auxVec12 = sum(sum(aux12));
    auxVec13 = sum(sum(aux13));
    auxVec14 = sum(sum(aux14));
    
    r00 = auxVec5 + auxVec6 + auxVec8 + auxVec9 + auxVec11 + auxVec12 + auxVec13;
    r10 = auxVec2 + auxVec3 + auxVec10;
    r01 = auxVec4 + auxVec7;
    r11 = auxVec1;
    mass = auxVec14;
    
    auxx1(:) = r00(1,1,:);
    auxx2(:) = r01(1,1,:);
    auxx3(:) = r10(1,1,:);
    auxx4(:) = r11(1,1,:);
    auxx5(:) = mass(1,1,:);
    r00  = auxx1;
    r01  = auxx2;
    r10  = auxx3;
    r11  = auxx4;
    mass = auxx5;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = msh.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_10 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_01 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_11 = spalloc (space.ndof, space.ndof, 3*space.ndof);
    Local_mass = spalloc (space.ndof, space.ndof, 3*space.ndof);
    rhs      = zeros (space.ndof, 1);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        MassLocal = mass((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00 = Local_00 + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, r00Local);
        Local_10 = Local_10 + op_gradu_v(spaceFunc{iel,2}, spaceFunc{iel,1}, spaceFunc{iel,3}, r10Local);
        Local_01 = Local_01 + op_u_gradv(spaceFunc{iel,1}, spaceFunc{iel,2}, spaceFunc{iel,3}, r01Local);
        Local_11 = Local_11 + op_gradu_gradv(spaceFunc{iel,2}, spaceFunc{iel,2}, spaceFunc{iel,3}, r11Local);
        Local_mass = Local_mass + op_u_v(spaceFunc{iel,1}, spaceFunc{iel,1}, spaceFunc{iel,3}, MassLocal);
        
        rhs      = rhs + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, fLocal);

    end
    
    Ml       = Local_mass;
    Al       = Local_11 + Local_10 + Local_01 + Local_00;
    bl       = rhs;
    
    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%

        % IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX
        
        if (imb == kmb)
            AModified = Al;
            AModified(1,:) = 0;
            AModified(1,1) = 1e30;
            AModified(end,:) = 0;
            AModified(end,end) = 1e30;
        else
            AModified = Al;
            AModified(1,:) = 0;
            AModified(end,:) = 0;
        end
        
        Al = AModified;
        
        % IMPLEMENT MODIFICATIONS IN THE FORCING TERM
        
        if (imb == kmb)
            bModified = bl;
            bModified(1) = 0;
            bModified(end) = 0;
        else
            bModified = bl;
        end   
        
        bl = bModified;

end

%% Method 'assemblerIGATransient'
            
function [Al,Ml,bl,a_ril,b_ril] = assemblerIGATransient(imb,kmb,wyq,mb_i,mb_yi,mb_k,mb_yk, ...
                                L,bc_up,bc_down,dato_uploc, dato_downloc,Coeff_forma,...
                                Computed,hx,p,k,dforma,in,fi)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Stiffness Matrix Using IGA Basis
    %   (2) Ml    : Final Assembled Mass Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler

    % Computation of the Offset Adjustment

    obj_liftBoundCond = BoundaryConditionHandler();
    
    obj_liftBoundCond.labelUpBoundCond = bc_up;
    obj_liftBoundCond.labelDownBoundCond = bc_down;
    obj_liftBoundCond.dataUpBoundCond = dato_uploc;
    obj_liftBoundCond.dataDownBoundCond = dato_downloc;
    obj_liftBoundCond.coeffForm = Coeff_forma;
    
    [a_ril,b_ril] = liftBoundCond(obj_liftBoundCond);
    aus = @(x,y) a_ril * y + b_ril;
    rilevamento_c = aus(0,Computed.y)';

    % Computation of the Force Exciting the System

    obj_integrate_1 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = L.*Computed.force_c - L.*(a_ril).*Computed.beta2_c - Computed.sigma_c.*L.*rilevamento_c;
    obj_integrate_1.funcWeight = mb_i.*wyq;
    
    forc  = integrate(obj_integrate_1);

    %---------------------------------------------------------------------%
    %                          DISCONTINUOUS FORCES                       %
    %---------------------------------------------------------------------%
    % Note: The following script can be used when the force exciting the
    %       system in discontinuous in the domain.
    %
    % forc = integrate(  - L.*a_ril.*Computed.beta2_c - ... 
    %                   Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);
    % nqny=64;
    % [~, adhocnodes, adhocw ] = gaussLegendre( nqny );
    % adhocnodes=adhocnodes/5+0.4;
    % adhocw=adhocw/5;
    % mb_iaus = newModalBasis( imb, adhocnodes, bc_up,bc_down,Coeff_forma);
    % forc = forc + integrate( L.*Computed.force_c, mb_iaus(:,imb).*adhocw);
    %---------------------------------------------------------------------%

    % Computation of the Coefficients 'r_{ik}^{st}'
    
    obj_integrate_2 = IntegrateHandler();
    obj_integrate_3 = IntegrateHandler();
    obj_integrate_4 = IntegrateHandler();
    obj_integrate_5 = IntegrateHandler();
    obj_integrate_6 = IntegrateHandler();
    obj_integrate_7 = IntegrateHandler();
    
    obj_integrate_2.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_3.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_4.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_5.funcToIntegrate = L.*Computed.beta2_c;
    obj_integrate_6.funcToIntegrate = L.*Computed.sigma_c;
    obj_integrate_7.funcToIntegrate = L.*ones(size(Computed.sigma_c));
    
    obj_integrate_2.funcWeight = mb_k .*mb_i .*wyq;
    obj_integrate_3.funcWeight = mb_yk.*mb_yi.*wyq;
    obj_integrate_4.funcWeight = mb_k .*mb_i .*wyq;
    obj_integrate_5.funcWeight = mb_yk.*mb_i .*wyq;
    obj_integrate_6.funcWeight = mb_k .*mb_i .*wyq;
    obj_integrate_7.funcWeight = mb_k .*mb_i .*wyq;

    muv             = integrate(obj_integrate_2);
    muv_y           = integrate(obj_integrate_3);
    betav           = integrate(obj_integrate_4);
    betav_y         = integrate(obj_integrate_5);
    sigmav          = integrate(obj_integrate_6);
    auxCoeffMass    = integrate(obj_integrate_7);

    %---------------------------------------------------------------------%
    % Note: To perform other cases in simulation we can also define the   %
    % following variables:                                                %
    % in = 0;                                                             %
    % fi = 1;                                                             %
    %---------------------------------------------------------------------%
    
    ngauss = 8;
    nqnx = 8;
    nel = 1/hx;   
    
    if p > 9
        ngauss = 15;
        display('ngauss = 15... check Gauss rule for very high degrees!')
    end

    % Set Mesh [Linear Parameterization, Continuity C^(p-k)]
    
    %---------------------------------------------------------------------%
    % Note:
    % The following functions, 'augknt' and 'bspkntins', are used to set
    % the mesh used in the solution of the problem with isogeometric
    % approach.
    %---------------------------------------------------------------------%
    
    knot      = augknt([in fi],p+1);
    h         = 1/(nel);
    ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
    ins       = (fi-in)*ins+in;
    
    %---------------------------------------------------------------------%
    % Note:
    % The function 'bspkntins' inserts the knots points into a univariate
    % B-Spline and is defined by the following parameters:
    %
    % (1ST INPUT) Degree of the B-Spline
    % (2ND INPUT) Control points, matrix of size (dim,nc)
    % (3RD INPUT) Knot sequence, row vector of size nk
    % (4TH INPUT) Row vector of knots to be inserted, size nu
    %
    % (1ST OUTPUT) Control points of the new B-Spline, of size (dim,nc+nu) 
    % (2ND OUTPUT) Knot vector of the new B-Spline, of size (nk+nu)
    %---------------------------------------------------------------------%
    
    [cp,knot] = bspkntins(p,in:(fi-in)*1/p:fi,knot,ins);

    % cp    : Vector Containing the Mesh Control Points 
    % knot  : Vector Containing the Mesh Nodes

    ncp   = length(cp);
    nspan = ncp - p;

    Jac = [];
    meshx_p1 = linspace(0,1,nel+1);

    for iel = 1:nel
        
        x_dx = meshx_p1(iel);                         % Right Extreme
        x_sx = meshx_p1(iel+1);                       % Left Extreme
        [point,~] = gauss(nqnx);                      % Creation of the 'p+1' Quadrature Points
        
        point = (x_sx-x_dx)*(0.5*point + 0.5) + x_dx;
        
        for j = 1 : nqnx
           Jac = [Jac sqrt(1+(dforma(point(j)))^2)];
        end
    end

    % This Loop Computes the 'meshS' (Curvilinear Abscissa)

    fi = 1; in = 0;
    [nodi,pesi] = gauss(nqnx);        % Nodes and Weights for the Quadrature Formula
    pesi = (fi-in)*(h)*(pesi*0.5);    % Rescale of the Weights Considering the Domain (Rescaled to [0,1])
    meshS    = zeros(nel+1,1);
    meshS(1) = 0;

    for j=1:nel
        
        ll = zeros(1,nqnx);
        nodi1 = (fi-in)*h*(nodi*0.5+0.5)+ in +(j-1)*h*(fi-in);
        
        for i=1:nqnx
            ll(i) = sqrt(1+(dforma(nodi1(i)))^2);
        end
        
        % debug
%         disp(size(ll))
%         disp(size(pesi))
        
        meshS(j+1) = meshS(j) + sum(ll*pesi);
        
    end

    % Vector Initialization
    
    f_gl  = zeros(ncp,1);
    row   = zeros(1,ncp^2);
    col   = zeros(1,ncp^2);
    val   = zeros(1,ncp^2);

    % Precompute the Values of Shape Functions and Derivatives at Gauss Points

    [gp,gw] = gauss(ngauss);

    gausspt = zeros(1,ngauss*nel);
    n = 0;
    
    for n0 = p+1:p+nspan
        
        if knot(n0) ~= knot(n0+1)
            gausspt(n*ngauss+1:(n+1)*ngauss) =...
                ((knot(n0+1) - knot(n0))*gp + knot(n0+1) + knot(n0))/2;
            n = n + 1;
        end
        
    end

    N  = spcol(knot,p+1,sort([gausspt,gausspt]));
    
    dN = N(2:2:end,:);     % IGA Basis
    N  = N(1:2:end,:);     % IGA Basis Derivative     
    
    clear gausspt

    % Loop Over Elements for Local Array Construction

    icount = 0;
    for iel = 1:nel
        
        dl = h/2;
        f_el  = zeros(p+1,1);      
        k_el  = zeros(p+1,p+1);    
        a_el  = zeros(p+1,p+1);
        a1_el = zeros(p+1,p+1);
        c_el  = zeros(p+1,p+1);
        c1_el = zeros(p+1,p+1);
        coeffMass = zeros(p+1,p+1);

        % Loop Over Gauss Points 

        [xgauss,~] = gauss(p+1);

        % Rescale of the Nodes Considering the Specific Element

        xgauss = (fi-in)*h*Jac(iel)*(0.5*xgauss+0.5)+meshS(iel);
        
        for igauss = 1:ngauss

            % Compute Shape Functions and Their Derivatives with Respect to
            % the Parametric Variable, Evaluated at Current Gauss Point

            shape  = N((iel-1) * ngauss + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            dshape = dN((iel-1) * ngauss + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            cpi    = cp((k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel); 

            % Compute the Coordinates of the Gauss Point in the Physical 
            % Space and the Jacobian at the Gauss Point

            x      = shape*cpi'*Jac((iel-1)*ngauss+igauss);       
            J      = dshape*cpi';

            % Compute the Shape Function Derivatives with respect to the
            % Physical Variable

            dshape = dshape/(J*Jac((iel-1)*ngauss+igauss));

            % Compute the Integration Weights

            gwt = gw(igauss)*(J*Jac((iel-1)*ngauss+igauss))*dl;

            % Manufactured Body Load
            % Compute Local R.H.S. Vector and Stiffness

            f_el  = f_el  + shape'*forc((iel-1)*nqnx+igauss)*gwt;                         % Force
            k_el  = k_el  + muv((iel-1)*nqnx+igauss)*(dshape'*dshape)*gwt;                % Diffusion
            a_el  = a_el  + betav((iel-1)*nqnx+igauss)*(shape'*dshape)*gwt;               % Transport
            c_el  = c_el  + sigmav((iel-1)*nqnx+igauss)*(shape'*shape)*gwt;               % Reaction
            a1_el = a1_el + muv_y((iel-1)*nqnx+igauss)*(shape'*shape)*gwt;                % Diffusion on Y
            c1_el = c1_el + betav_y((iel-1)*nqnx+igauss)*(shape'*shape)*gwt;              % Transport on Y
            coeffMass = coeffMass + auxCoeffMass((iel-1)*nqnx+igauss)*(shape'*shape)*gwt; % Mass
        end

        % Assemble Global R.H.S. Vector and Prepare Vectors for Sparse Stiffness Assembly

        for a = 1:p+1
            
            i = (k-1)*(iel-1) + iel + p - (a - 1);

            if (i>1)

                f_gl(i) = f_gl(i) + f_el(a);

                for b = 1:p+1
                    
                    j = (k-1)*(iel-1) + iel + p - (b - 1);
                    icount = icount + 1;
                    
                    row(icount)   = i;
                    col(icount)   = j;
                    
                    val(icount)   = k_el(a,b);
                    val1(icount)  = a_el(a,b);
                    val2(icount)  = c_el(a,b);
                    val3(icount)  = a1_el(a,b);
                    val4(icount)  = c1_el(a,b);
                    val5(icount)  = coeffMass(a,b);
                end
            else 
                icount = icount + 1;
                row(icount)  = i;
                col(icount)  = i;
                val(icount)  = 1;
                val1(icount) = 0;
                val2(icount) = 0;
                val3(icount) = 0;
                val4(icount) = 0;
                val5(icount) = 0;
            end
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    val = val(1:icount);

    % Assemble Stiffness and the Other Matrices

    k_gl      = sparse(row,col,val,ncp,ncp);
    a_gl      = sparse(row,col,val1,ncp,ncp);
    c_gl      = sparse(row,col,val2,ncp,ncp);
    a1_gl     = sparse(row,col,val3,ncp,ncp);
    c1_gl     = sparse(row,col,val4,ncp,ncp);
    coeffMass = sparse(row,col,val5,ncp,ncp);
    
    Al = k_gl + a_gl + c_gl + a1_gl + c1_gl;
    Ml = coeffMass;
    bl = f_gl;
    
    % IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX

    if (imb == kmb)
        AModified = Al;
        AModified(1,:) = 0;
        AModified(1,1) = 1e30;
        AModified(end,:) = 0;
        AModified(end,end) = 1e30;
    else
        AModified = Al;
        AModified(1,:) = 0;
        AModified(end,:) = 0;
    end

    Al = AModified;

    % IMPLEMENT MODIFICATIONS IN THE FORCING TERM

    if (imb == kmb)
        bModified = bl;
        bModified(1) = 0;
        bModified(end) = 0;
    else
        bModified = bl;
    end   

    bl = bModified;

end

%% Method 'assemblerIGATransientFull'

function [Al,Ml,bl,aLift,bLift] = assemblerIGATransientFull(imb,kmb,numbControlPts, ...
                                augVerWeights,mb_i,mb_yi,mb_k,mb_yk, ...
                                L, Computed, ...
                                horStep,p,k, ...
                                domainLeftLimit,domainRightLimit, ...
                                lifting,aLift,bLift,jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights,D1,D2)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%

    objForceInt = IntegrateHandler();
    
    objForceInt.funcToIntegrate =   L.*Computed.force_c - ...
                                    L.*(aLift).*Computed.beta2_c - ...
                                    Computed.sigma_c.*L.*lifting;
                                
    objForceInt.funcWeight = mb_i.*augVerWeights;
    
    forceVec  = integrate(objForceInt);
    
%     % Debug
%     %---------------------------------------------------------------------%
%     
%     disp(size(forceVec))
%     disp(size(Computed.force_c))
%     disp(size(mb_i))
    
%     disp(max(max(lifting)))
%     disp(min(min(lifting)))
%     disp(max(max(aLift)))
%     disp(min(min(aLift)))
%     disp(max(max(bLift)))
%     disp(min(min(bLift)))
    %---------------------------------------------------------------------%
    
    % Debug 
    %---------------------------------------------------------------------%
    %disp(['Max L : ',num2str(max(max(L)))]);
    %disp(['Min L : ',num2str(min(min(L)))]);
    %plot(linspace(0,length(forceVec),length(forceVec)),forceVec);
    %---------------------------------------------------------------------%

    %---------------------------------------------------------------------%
    %                          DISCONTINUOUS FORCES                       %
    %---------------------------------------------------------------------%
    % Note: The following script can be used when the force exciting the
    %       system in discontinuous in the domain.
    %
    % forc = integrate(  - L.*a_ril.*Computed.beta2_c - ... 
    %                   Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);
    % nqny=64;
    % [~, adhocnodes, adhocw ] = gaussLegendre( nqny );
    % adhocnodes=adhocnodes/5+0.4;
    % adhocw=adhocw/5;
    % mb_iaus = newModalBasis( imb, adhocnodes, bc_up,bc_down,Coeff_forma);
    % forc = forc + integrate( L.*Computed.force_c, mb_iaus(:,imb).*adhocw);
    %---------------------------------------------------------------------%
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    obj_integrate_3 = IntegrateHandler();
    obj_integrate_4 = IntegrateHandler();
    obj_integrate_5 = IntegrateHandler();
    obj_integrate_6 = IntegrateHandler();
    obj_integrate_7 = IntegrateHandler();
    obj_integrate_8 = IntegrateHandler();
    obj_integrate_9 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_2.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_3.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_4.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_5.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_6.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_7.funcToIntegrate = L.*Computed.beta2_c;
    obj_integrate_8.funcToIntegrate = L.*Computed.sigma_c;
    obj_integrate_9.funcToIntegrate = L.*ones(size(Computed.sigma_c));
    
    obj_integrate_1.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_yi .*augVerWeights;
    obj_integrate_3.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    obj_integrate_4.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_5.funcWeight = mb_yk .*mb_yi .*augVerWeights;
    obj_integrate_6.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_7.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_8.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    obj_integrate_9.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    lambda_3   = integrate(obj_integrate_3);
    lambda_4   = integrate(obj_integrate_4);
    lambda_5   = integrate(obj_integrate_5);
    lambda_6   = integrate(obj_integrate_6);
    lambda_7   = integrate(obj_integrate_7);
    lambda_8   = integrate(obj_integrate_8);
    lambda_9   = integrate(obj_integrate_9);
    
    r00 = lambda_5 * (D1^2 + D2^2) + lambda_6 * D1 + lambda_7 * D2 + lambda_8;
    r10 = lambda_2 * D1 + lambda_3;
    r01 = lambda_4 * D1;
    r11 = lambda_1;
    
    auxCoeffMass = lambda_9;

    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = p;
    objBasisIGA.continuityParameter   = k;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGA,derBasisIGA,controlPts] = newIsoGeoBasis(objBasisIGA);  
    
    %% LOCAL JACOBIAN
    
%     Jac = [];
%     meshx_p1 = linspace(0,1,nel+1);
% 
%     for iel = 1:nel
%         
%         x_dx = meshx_p1(iel);                         % Right Extreme
%         x_sx = meshx_p1(iel+1);                       % Left Extreme
%         [point,~] = gauss(nqnx);                      % Creation of the 'p+1' Quadrature Points
%         
%         point = (x_sx-x_dx)*(0.5*point + 0.5) + x_dx;
%         
%         for j = 1 : nqnx
%            Jac = [Jac sqrt(1+(dforma(point(j)))^2)];
%         end
%     end
    
    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    globalForce  = zeros(numbControlPts,1);
    row   = zeros(1,numbControlPts^2);
    col   = zeros(1,numbControlPts^2);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        forceElem  = zeros(p+1,1);
        
        deltaElem1 = zeros(p+1,p+1);
        deltaElem2 = zeros(p+1,p+1);
        deltaElem3 = zeros(p+1,p+1);
        deltaElem4 = zeros(p+1,p+1);
        coeffMass  = zeros(p+1,p+1);
        
        % INTEGRAL OVER THE SUPPORTING FIBER
        %-----------------------------------------------------------------%
        % We compute the coefficients referring to the integrals along the
        % supporting fiber.
        %-----------------------------------------------------------------%
        
        for igauss = 1:numbHorQuadNodes

            % SHAPE FUNCTIONS
            %-------------------------------------------------------------%
            % Compute the shape functions and their derivatives with
            % respect to the parametric variable, evaluated over the
            % Gauss-Legendre integration points.
            %-------------------------------------------------------------%

            shape  = basisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            dshape = derBasisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            cpi    = controlPts((k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            
            % COORDINATE OF THE GAUSS POINT IN THE PHYSICAL DOMAIN
            %-------------------------------------------------------------%
            % This componponent is not being used now because the thickness
            % of the channel is constant and equal to '1'. However, it is
            % extremely necessary for more complex geometries because it
            % gives the correct horizontal component in the physical domain
            % to evaluate the function L(x).
            %-------------------------------------------------------------%
            
            x      = shape*cpi'*jacIGA((iel-1)*numbHorQuadNodes+igauss);
            
            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshape*cpi';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshape = dshape/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL RHS VECTOR
            %-------------------------------------------------------------%
            % The current code uses 'numbHorQuadNodes' instead of
            % 'numbVerQuadNodes' in the next 5 lines of code. Check if
            % the solution is consistent.
            %-------------------------------------------------------------%

            forceElem       = forceElem  + shape'*forceVec((iel-1)*numbHorQuadNodes+igauss)*gwt;    % Force
            
            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shape' * shape) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (shape' *dshape) * gwt;
            deltaElem3 = deltaElem3 + r01((iel-1)*numbHorQuadNodes+igauss) * (dshape'* shape) * gwt;
            deltaElem4 = deltaElem4 + r11((iel-1)*numbHorQuadNodes+igauss) * (dshape'*dshape) * gwt;
            
            % LOCAL COMPONENTS OF THE MASS MATRIX
            
            coeffMass = coeffMass + auxCoeffMass((iel-1)*numbHorQuadNodes+igauss) * (shape' * shape) * gwt;
            
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:p+1
            
            i = (k-1)*(iel-1) + iel + p - (a - 1);

            if (i>1)

                globalForce(i) = globalForce(i) + forceElem(a);

                for b = 1:p+1
                    
                    j = (k-1)*(iel-1) + iel + p - (b - 1);
                    icount = icount + 1;
                    
                    row(icount)   = i;
                    col(icount)   = j;
                    
                    deltaLocal1(icount) = deltaElem1(a,b);
                    deltaLocal2(icount) = deltaElem2(a,b);
                    deltaLocal3(icount) = deltaElem3(a,b);
                    deltaLocal4(icount) = deltaElem4(a,b);
                    deltaLocal5(icount) = coeffMass(a,b);
                    
                end
            else 
                icount = icount + 1;
                row(icount)  = i;
                col(icount)  = i;
                
                deltaLocal1(icount) = 0;
                deltaLocal2(icount) = 0;
                deltaLocal3(icount) = 0;
                deltaLocal4(icount) = 1;
                deltaLocal5(icount) = 0;
                
            end
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPts,numbControlPts);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPts,numbControlPts);
    Delta3 = sparse(row,col,deltaLocal3,numbControlPts,numbControlPts);
    Delta4 = sparse(row,col,deltaLocal4,numbControlPts,numbControlPts);
    CoeffMass = sparse(row,col,deltaLocal4,numbControlPts,numbControlPts);
    
    Al   = Delta1 + Delta2 + Delta3 + Delta4;
    Ml   = CoeffMass;
    bl   = globalForce;
    
    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%

        % IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX
        
        if (imb == kmb)
            AModified = Al;
            AModified(1,:) = 0;
            AModified(1,1) = 1e30;
            AModified(end,:) = 0;
            AModified(end,end) = 1e30;
        else
            AModified = Al;
            AModified(1,:) = 0;
            AModified(end,:) = 0;
        end
        
        Al = AModified;
        
        % IMPLEMENT MODIFICATIONS IN THE FORCING TERM
        
        if (imb == kmb)
            bModified = bl;
            bModified(1) = 0;
            bModified(end) = 0;
        else
            bModified = bl;
        end   
        
        bl = bModified;

    disp('Finished ASSEMBLING LOOP');
end

%% Method 'assemblerIGAForce'
            
function [bl,a_ril,b_ril] = assemblerIGAForce(wyq,mb_i,L,bc_up,bc_down, ...
                                dato_uploc,dato_downloc,Coeff_forma,Computed, ...
                                hx,p,k,dforma,in,fi)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Stiffness Matrix Using IGA Basis
    %   (2) Ml    : Final Assembled Mass Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler

    % Computation of the Offset Adjustment

    obj_liftBoundCond = BoundaryConditionHandler();
    
    obj_liftBoundCond.labelUpBoundCond = bc_up;
    obj_liftBoundCond.labelDownBoundCond = bc_down;
    obj_liftBoundCond.dataUpBoundCond = dato_uploc;
    obj_liftBoundCond.dataDownBoundCond = dato_downloc;
    obj_liftBoundCond.coeffForm = Coeff_forma;
    
    [a_ril,b_ril] = liftBoundCond(obj_liftBoundCond);
    aus = @(x,y) a_ril * y + b_ril;
    rilevamento_c = aus(0,Computed.y)';

    % Computation of the Force Exciting the System

    obj_integrate_1 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = L.*Computed.force_c - L.*(a_ril).*Computed.beta2_c - Computed.sigma_c.*L.*rilevamento_c;
    obj_integrate_1.funcWeight = mb_i.*wyq;
    
    forc  = integrate(obj_integrate_1);

    %---------------------------------------------------------------------%
    %                          DISCONTINUOUS FORCES                       %
    %---------------------------------------------------------------------%
    % Note: The following script can be used when the force exciting the
    %       system in discontinuous in the domain.
    %
    % forc = integrate(  - L.*a_ril.*Computed.beta2_c - ... 
    %                   Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);
    % nqny=64;
    % [~, adhocnodes, adhocw ] = gaussLegendre( nqny );
    % adhocnodes=adhocnodes/5+0.4;
    % adhocw=adhocw/5;
    % mb_iaus = newModalBasis( imb, adhocnodes, bc_up,bc_down,Coeff_forma);
    % forc = forc + integrate( L.*Computed.force_c, mb_iaus(:,imb).*adhocw);
    %---------------------------------------------------------------------%    

    %---------------------------------------------------------------------%
    % Note: To perform other cases in simulation we can also define the   %
    % following variables:                                                %
    % in = 0;                                                             %
    % fi = 1;                                                             %
    %---------------------------------------------------------------------%
    
    ngauss = 8;
    nqnx = 8;
    nel = 1/hx;   
    
    if p > 9
        ngauss = 15;
        display('ngauss = 15... check Gauss rule for very high degrees!')
    end

    % Set Mesh [Linear Parameterization, Continuity C^(p-k)]
    
    %---------------------------------------------------------------------%
    % Note:
    % The following functions, 'augknt' and 'bspkntins', are used to set
    % the mesh used in the solution of the problem with isogeometric
    % approach.
    %---------------------------------------------------------------------%
    
    knot      = augknt([in fi],p+1);
    h         = 1/(nel);
    ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
    ins       = (fi-in)*ins+in;
    
    %---------------------------------------------------------------------%
    % Note:
    % The function 'bspkntins' inserts the knots points into a univariate
    % B-Spline and is defined by the following parameters:
    %
    % (1ST INPUT) Degree of the B-Spline
    % (2ND INPUT) Control points, matrix of size (dim,nc)
    % (3RD INPUT) Knot sequence, row vector of size nk
    % (4TH INPUT) Row vector of knots to be inserted, size nu
    %
    % (1ST OUTPUT) Control points of the new B-Spline, of size (dim,nc+nu) 
    % (2ND OUTPUT) Knot vector of the new B-Spline, of size (nk+nu)
    %---------------------------------------------------------------------%
    
    [cp,knot] = bspkntins(p,in:(fi-in)*1/p:fi,knot,ins);

    % cp    : Vector Containing the Mesh Control Points 
    % knot  : Vector Containing the Mesh Nodes

    ncp   = length(cp);
    nspan = ncp - p;

    Jac = [];
    meshx_p1 = linspace(0,1,nel+1);

    for iel = 1:nel
        
        x_dx = meshx_p1(iel);                         % Right Extreme
        x_sx = meshx_p1(iel+1);                       % Left Extreme
        [point,~] = gauss(nqnx);                      % Creation of the 'p+1' Quadrature Points
        
        point = (x_sx-x_dx)*(0.5*point + 0.5) + x_dx;
        
        for j = 1 : nqnx
           Jac = [Jac sqrt(1+(dforma(point(j)))^2)];
        end
    end

    % This Loop Computes the 'meshS' (Curvilinear Abscissa)

    fi = 1; in = 0;
    [nodi,pesi] = gauss(nqnx);        % Nodes and Weights for the Quadrature Formula
    pesi = (fi-in)*(h)*(pesi*0.5);    % Rescale of the Weights Considering the Domain (Rescaled to [0,1])
    meshS    = zeros(nel+1,1);
    meshS(1) = 0;

    for j=1:nel
        
        ll = zeros(1,nqnx);
        nodi1 = (fi-in)*h*(nodi*0.5+0.5)+ in +(j-1)*h*(fi-in);
        
        for i=1:nqnx
            ll(i) = sqrt(1+(dforma(nodi1(i)))^2);
        end
        
        meshS(j+1) = meshS(j) + sum(ll*pesi);
        
    end

    % Vector Initialization
    
    f_gl  = zeros(ncp,1);
    
    % Precompute the Values of Shape Functions and Derivatives at Gauss Points

    [gp,gw] = gauss(ngauss);

    gausspt = zeros(1,ngauss*nel);
    n = 0;
    
    for n0 = p+1:p+nspan;
        
        if knot(n0) ~= knot(n0+1)
            gausspt(n*ngauss+1:(n+1)*ngauss) =...
                ((knot(n0+1) - knot(n0))*gp + knot(n0+1) + knot(n0))/2;
            n = n + 1;
        end
        
    end

    N  = spcol(knot,p+1,sort([gausspt,gausspt]));
    
    dN = N(2:2:end,:);     % IGA Basis
    N  = N(1:2:end,:);     % IGA Basis Derivative     
    
    clear gausspt

    % Loop Over Elements for Local Array Construction

    for iel = 1:nel
        
        dl = h/2;
        f_el  = zeros(p+1,1);      
        
        for igauss = 1:ngauss

            % Compute Shape Functions and Their Derivatives with Respect to
            % the Parametric Variable, Evaluated at Current Gauss Point

            shape  = N((iel-1) * ngauss + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            dshape = dN((iel-1) * ngauss + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            cpi    = cp((k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel); 

            % Compute the Coordinates of the Gauss Point in the Physical 
            % Space and the Jacobian at the Gauss Point
       
            J      = dshape*cpi';
            
            % Compute the Integration Weights

            gwt = gw(igauss)*(J*Jac((iel-1)*ngauss+igauss))*dl;

            % Manufactured Body Load
            % Compute Local R.H.S. Vector and Stiffness

            f_el  = f_el  + shape'*forc((iel-1)*nqnx+igauss)*gwt;                         % Force
            
        end

        % Assemble Global R.H.S. Vector and Prepare Vectors for Sparse Stiffness Assembly

        for a = 1:p+1
            
            i = (k-1)*(iel-1) + iel + p - (a - 1);

            if (i>1)

                f_gl(i) = f_gl(i) + f_el(a);

            end
        end
    end

    % Assemble Force Vector
    
    bl = f_gl;

end

%% Method 'assemblerIGASTR'
            
function [Al,bl,aLift,bLift] = assemblerIGASTR(imb,kmb,numbControlPts, ...
                                augVerWeights,mb_i,mb_yi,mb_k,mb_yk, ...
                                L, Computed, ...
                                horStep,p,k, ...
                                domainLeftLimit,domainRightLimit, ...
                                lifting,aLift,bLift,jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights,delta,Pechlet)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%

    objForceInt = IntegrateHandler();
    
    objForceInt.funcToIntegrate =   L.*Computed.force_c - ...
                                    L.*(aLift).*Computed.beta2_c - ...
                                    Computed.sigma_c.*L.*lifting;
                                
    objForceInt.funcWeight = mb_i.*augVerWeights;
    
    forceVec  = integrate(objForceInt);
    
%     % Debug
%     %---------------------------------------------------------------------%
%     
%     disp(size(forceVec))
%     disp(size(Computed.force_c))
%     disp(size(mb_i))
    
%     disp(max(max(lifting)))
%     disp(min(min(lifting)))
%     disp(max(max(aLift)))
%     disp(min(min(aLift)))
%     disp(max(max(bLift)))
%     disp(min(min(bLift)))
    %---------------------------------------------------------------------%
    
    % Debug 
    %---------------------------------------------------------------------%
    %disp(['Max L : ',num2str(max(max(L)))]);
    %disp(['Min L : ',num2str(min(min(L)))]);
    %plot(linspace(0,length(forceVec),length(forceVec)),forceVec);
    %---------------------------------------------------------------------%

    %---------------------------------------------------------------------%
    %                          DISCONTINUOUS FORCES                       %
    %---------------------------------------------------------------------%
    % Note: The following script can be used when the force exciting the
    %       system in discontinuous in the domain.
    %
    % forc = integrate(  - L.*a_ril.*Computed.beta2_c - ... 
    %                   Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);
    % nqny=64;
    % [~, adhocnodes, adhocw ] = gaussLegendre( nqny );
    % adhocnodes=adhocnodes/5+0.4;
    % adhocw=adhocw/5;
    % mb_iaus = newModalBasis( imb, adhocnodes, bc_up,bc_down,Coeff_forma);
    % forc = forc + integrate( L.*Computed.force_c, mb_iaus(:,imb).*adhocw);
    %---------------------------------------------------------------------%
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    % NOTE:  METHOD STABILIZE THEN REDUC?
    % STABILIZATION with mu->  :
    % mu*(1+beta1*sqrt(hx^2+1)/(2*mu))
    %---------------------------------------------------------------------%

    D1 = 0;
    D2 = 1;
    
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    obj_integrate_3 = IntegrateHandler();
    obj_integrate_4 = IntegrateHandler();
    obj_integrate_5 = IntegrateHandler();
    obj_integrate_6 = IntegrateHandler();
    obj_integrate_7 = IntegrateHandler();
    obj_integrate_8 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = L.*Computed.mu_c.*(1+delta.*Pechlet);
    obj_integrate_2.funcToIntegrate = L.*Computed.mu_c.*(1+delta.*Pechlet);
    obj_integrate_3.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_4.funcToIntegrate = L.*Computed.mu_c.*(1+delta.*Pechlet);
    obj_integrate_5.funcToIntegrate = L.*Computed.mu_c.*(1+delta.*Pechlet);
    obj_integrate_6.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_7.funcToIntegrate = L.*Computed.beta2_c;
    obj_integrate_8.funcToIntegrate = L.*Computed.sigma_c;
    
    obj_integrate_1.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_yi .*augVerWeights;
    obj_integrate_3.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    obj_integrate_4.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_5.funcWeight = mb_yk .*mb_yi .*augVerWeights;
    obj_integrate_6.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_7.funcWeight = mb_yk .*mb_i  .*augVerWeights;
    obj_integrate_8.funcWeight = mb_k  .*mb_i  .*augVerWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    lambda_3   = integrate(obj_integrate_3);
    lambda_4   = integrate(obj_integrate_4);
    lambda_5   = integrate(obj_integrate_5);
    lambda_6   = integrate(obj_integrate_6);
    lambda_7   = integrate(obj_integrate_7);
    lambda_8   = integrate(obj_integrate_8);
    
    r00 = lambda_5 * (D1^2 + D2^2) + lambda_6 * D1 + lambda_7 * D2 + lambda_8;
    r10 = lambda_2 * D1 + lambda_3;
    r01 = lambda_4 * D1;
    r11 = lambda_1;

    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = p;
    objBasisIGA.continuityParameter   = k;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGA,derBasisIGA,controlPts] = newIsoGeoBasis(objBasisIGA);  
    
    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    globalForce  = zeros(numbControlPts,1);
    row   = zeros(1,numbControlPts^2);
    col   = zeros(1,numbControlPts^2);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        forceElem  = zeros(p+1,1);
        
        deltaElem1 = zeros(p+1,p+1);
        deltaElem2 = zeros(p+1,p+1);
        deltaElem3 = zeros(p+1,p+1);
        deltaElem4 = zeros(p+1,p+1);
        
        % INTEGRAL OVER THE SUPPORTING FIBER
        %-----------------------------------------------------------------%
        % We compute the coefficients referring to the integrals along the
        % supporting fiber.
        %-----------------------------------------------------------------%
        
        for igauss = 1:numbHorQuadNodes

            % SHAPE FUNCTIONS
            %-------------------------------------------------------------%
            % Compute the shape functions and their derivatives with
            % respect to the parametric variable, evaluated over the
            % Gauss-Legendre integration points.
            %-------------------------------------------------------------%

            shape  = basisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            dshape = derBasisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            cpi    = controlPts((k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            
            % COORDINATE OF THE GAUSS POINT IN THE PHYSICAL DOMAIN
            %-------------------------------------------------------------%
            % This componponent is not being used now because the thickness
            % of the channel is constant and equal to '1'. However, it is
            % extremely necessary for more complex geometries because it
            % gives the correct horizontal component in the physical domain
            % to evaluate the function L(x).
            %-------------------------------------------------------------%
            
            x      = shape*cpi'*jacIGA((iel-1)*numbHorQuadNodes+igauss);
            
            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshape*cpi';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshape = dshape/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL RHS VECTOR
            %-------------------------------------------------------------%
            % The current code uses 'numbHorQuadNodes' instead of
            % 'numbVerQuadNodes' in the next 5 lines of code. Check if
            % the solution is consistent.
            %-------------------------------------------------------------%

            forceElem       = forceElem  + shape'*forceVec((iel-1)*numbHorQuadNodes+igauss)*gwt;    % Force
            
            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shape' * shape) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (shape' *dshape) * gwt;
            deltaElem3 = deltaElem3 + r01((iel-1)*numbHorQuadNodes+igauss) * (dshape'* shape) * gwt;
            deltaElem4 = deltaElem4 + r11((iel-1)*numbHorQuadNodes+igauss) * (dshape'*dshape) * gwt;
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:p+1
            
            i = (k-1)*(iel-1) + iel + p - (a - 1);

            if (i>1)

                globalForce(i) = globalForce(i) + forceElem(a);

                for b = 1:p+1
                    
                    j = (k-1)*(iel-1) + iel + p - (b - 1);
                    icount = icount + 1;
                    
                    row(icount)   = i;
                    col(icount)   = j;
                    
                    deltaLocal1(icount) = deltaElem1(a,b);
                    deltaLocal2(icount) = deltaElem2(a,b);
                    deltaLocal3(icount) = deltaElem3(a,b);
                    deltaLocal4(icount) = deltaElem4(a,b);
                    
                end
            else 
                icount = icount + 1;
                row(icount)  = i;
                col(icount)  = i;
                
                deltaLocal1(icount) = 0;
                deltaLocal2(icount) = 0;
                deltaLocal3(icount) = 0;
                deltaLocal4(icount) = 1;
                
            end
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPts,numbControlPts);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPts,numbControlPts);
    Delta3 = sparse(row,col,deltaLocal3,numbControlPts,numbControlPts);
    Delta4 = sparse(row,col,deltaLocal4,numbControlPts,numbControlPts);
    
    Al   = Delta1 + Delta2 + Delta3 + Delta4;
    bl   = globalForce;
    
    % Debug
    %---------------------------------------------------------------------%
    % figure; plot(bl);
    %---------------------------------------------------------------------%
    
    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%
    
    if (imb == kmb)
        Al(1,1)  = 1;
        Al(1,2)  = 0;
        Al(1,3)  = 0;
    else
        Al(1,1)  = 1e-15;
        Al(1,2)  = 0;
        Al(1,3)  = 0;
    end

    disp('Finished ASSEMBLING LOOP');
end

%% Method 'assemblerIGARTS'
            
function [Al,bl,aLift,bLift] = assemblerIGARTS(imb,kmb,numbControlPts, ...
                                verGLWeights,mb_i,mb_yi,mb_k,mb_yk, ...
                                L, Computed, ...
                                horStep,p,k, ...
                                domainLeftLimit,domainRightLimit, ...
                                lifting,aLift,bLift,jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights,delta)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (2)  kmb          : Coefficients Used in the Assembling Operation
    %                       Loop
    %   (3)  size_fb      : Size of the Functional Basis. Equal to the 
    %                       Dimention of the Spline Basis (Number of
    %                       Control Points)
    %   (4)  mesh_wx      : Vector Containing the Weights of the Quadrature 
    %                       Nodes in the Whole Mesh 
    %   (5)  wyq          : Vector Containing the Weights of the
    %                       Gauss-Legendre Integration Nodes
    %   (6)  mb_i         : Modal Basis referring for coeff. (imb,kmb)
    %   (7)  mb_yi        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (8)  mb_k         : Modal Basis referring to coeff. (imb,kmb)
    %   (9)  mb_yk        : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L            : Thickness of the Channel Used in the Example
    %   (11) bc_up        : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Upper Limit of
    %                       the Domain
    %   (11) bc_down      : Contains the Label Identifying the Nature of
    %                       the Boundary Conditions on the Lower Limit of
    %                       the Domain
    %   (12) dato_uploc         : Contains the Values of the Boundary Conditions
    %                             on the Upper Limir of the Domain
    %   (13) dato_downloc       : Contains the Values of the Boundary Conditions
    %                             on the Lower Limir of the Domain
    %   (14) posizione_dominio (domain position)  : ??
    %   (15) Coeff_forma        : Data Strusture Containing All the @-Functions
    %                             and the Constants Relative to the Bilinear Form
    %   (16) gamma        : Data Structure Containing the Two Values of the
    %                       Coefficients (R, L) for the Robin Condition Used
    %                       in the Domain Decomposition
    %   (17) Computed     : Structure containing the coefficiets mu, beta
    %                       1, beta 2, sigma and the force exciting the system 
    %   (18) nd           : Number of Elements in the Vector Containing the
    %                       Dimensions of the Modal Basis in Each Domain
    %   (19) coupling     : Contains the Label Adressing the Coupling Condition
    %                       of the Problem
    %   (20) hx           : Vector Containing the Step of the Finite
    %                       Element Mesh
    %   (20) p            : Degree of the Polynomial B-Spline Basis
    %   (21) k            : Degree of Continuity of the Basis 'C^(p-k)'
    %   (23) dforma       : Symbolic Function Defining the Derivative of
    %                       the Profile of the Simulation Domain
    %   (24) in, fi       : Values Defining the Extremes of the Domains 
    %                       Used to Compute the Error of the Solution 
    %   (25) delta        : Parameter to tune for stabilization
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) a_ril : Primo Coefficiente di Rilevamento
    %   (4) b_ril : Secondo Coefficiente di Rilevamento
    %
    % Notes: (1) The input variables 'degree', 'degree_k', 'n', 'nqnx',
    %            'mesh_wx', 'IGA_basis', 'IGA_basis_x', 'mesh_fisx', 
    %            'mesh_fisy' and 'forma' apparently are not used in this
    %            function.
    %        (2) The inputs 'bc_up' and 'bc_down' currently are able to
    %            receive the labels 'rob' and 'dir', corresponding
    %            respectively to the Robin and Dirichlet boundary
    %            condditions. The Neumann boundary conditions are not yet
    %            added to the code, but should not be hard to implement.
    %
    %            'rob':  mu * du/dnu + coeffrobin * u = data
    % 		     'dir':  u = data

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% JACOBIAN IMPOSITIONS
    %---------------------------------------------------------------------%
    % This values will come from the demo when the code is changed to
    % incorporate also cases where L(x) is not constant. Note that the name
    % of the variables were chosen to match the same name in the reference
    % paper used to create this function.
    %---------------------------------------------------------------------%
    
    D1 = 0;
    D2 = 1;
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%

    objForceInt = IntegrateHandler();
    
    objForceInt.funcToIntegrate =   L.*Computed.force_c - ...
                                    L.*(aLift).*Computed.beta2_c - ...
                                    Computed.sigma_c.*L.*lifting;
                                
    objForceInt.funcWeight = mb_i.*verGLWeights;
    
    forceVec  = integrate(objForceInt);

    %---------------------------------------------------------------------%
    %                          DISCONTINUOUS FORCES                       %
    %---------------------------------------------------------------------%
    % Note: The following script can be used when the force exciting the
    %       system in discontinuous in the domain.
    %
    % forc = integrate(  - L.*a_ril.*Computed.beta2_c - ... 
    %                   Computed.sigma_c.*L.*rilevamento_c, mb_i.*wyq);
    % nqny=64;
    % [~, adhocnodes, adhocw ] = gaussLegendre( nqny );
    % adhocnodes=adhocnodes/5+0.4;
    % adhocw=adhocw/5;
    % mb_iaus = newModalBasis( imb, adhocnodes, bc_up,bc_down,Coeff_forma);
    % forc = forc + integrate( L.*Computed.force_c, mb_iaus(:,imb).*adhocw);
    %---------------------------------------------------------------------%
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    % NOTE:  METHOD REDUCE THEN STABILIZE
    % STABILIZATION with mu->  :
    % mu*(1+delta*betainf.*h./(2*mu)  ) 
    %---------------------------------------------------------------------%
    betainf=max(max(sqrt(Computed.beta1_c.^2+Computed.beta2_c.^2)));
       
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    obj_integrate_3 = IntegrateHandler();
    obj_integrate_4 = IntegrateHandler();
    obj_integrate_5 = IntegrateHandler();
    obj_integrate_6 = IntegrateHandler();
    obj_integrate_7 = IntegrateHandler();
    obj_integrate_8 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = L.*Computed.mu_c.*(1+delta*betainf.*horStep./(2*Computed.mu_c)  );
    obj_integrate_2.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_3.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_4.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_5.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_6.funcToIntegrate = L.*Computed.beta1_c;
    obj_integrate_7.funcToIntegrate = L.*Computed.beta2_c;
    obj_integrate_8.funcToIntegrate = L.*Computed.sigma_c;
    
    obj_integrate_1.funcWeight = mb_k  .*mb_i  .*verGLWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_yi .*verGLWeights;
    obj_integrate_3.funcWeight = mb_k  .*mb_i  .*verGLWeights;
    obj_integrate_4.funcWeight = mb_yk .*mb_i  .*verGLWeights;
    obj_integrate_5.funcWeight = mb_yk .*mb_yi .*verGLWeights;
    obj_integrate_6.funcWeight = mb_yk .*mb_i  .*verGLWeights;
    obj_integrate_7.funcWeight = mb_yk .*mb_i  .*verGLWeights;
    obj_integrate_8.funcWeight = mb_k  .*mb_i  .*verGLWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    lambda_3   = integrate(obj_integrate_3);
    lambda_4   = integrate(obj_integrate_4);
    lambda_5   = integrate(obj_integrate_5);
    lambda_6   = integrate(obj_integrate_6);
    lambda_7   = integrate(obj_integrate_7);
    lambda_8   = integrate(obj_integrate_8);
    
    r00 = lambda_5 * (D1^2 + D2^2) + lambda_6 * D1 + lambda_7 * D2 + lambda_8;
    r10 = lambda_2 * D1 + lambda_3;
    r01 = lambda_4 * D1;
    r11 = lambda_1;

    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = p;
    objBasisIGA.continuityParameter   = k;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGA,derBasisIGA,controlPts] = newIsoGeoBasis(objBasisIGA);
    
    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    globalForce  = zeros(numbControlPts,1);
    row   = zeros(1,numbControlPts^2);
    col   = zeros(1,numbControlPts^2);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        forceElem  = zeros(p+1,1);
        
        deltaElem1 = zeros(p+1,p+1);
        deltaElem2 = zeros(p+1,p+1);
        deltaElem3 = zeros(p+1,p+1);
        deltaElem4 = zeros(p+1,p+1);
        
        % INTEGRAL OVER THE SUPPORTING FIBER
        %-----------------------------------------------------------------%
        % We compute the coefficients referring to the integrals along the
        % supporting fiber.
        %-----------------------------------------------------------------%
        
        for igauss = 1:numbHorQuadNodes

            % SHAPE FUNCTIONS
            %-------------------------------------------------------------%
            % Compute the shape functions and their derivatives with
            % respect to the parametric variable, evaluated over the
            % Gauss-Legendre integration points.
            %-------------------------------------------------------------%

            shape  = basisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            dshape = derBasisIGA((iel-1) * numbHorQuadNodes + igauss,(k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            cpi    = controlPts((k-1) * (iel-1) + iel + p:-1:(k-1) * (iel-1) + iel);
            
            % COORDINATE OF THE GAUSS POINT IN THE PHYSICAL DOMAIN
            %-------------------------------------------------------------%
            % This componponent is not being used now because the thickness
            % of the channel is constant and equal to '1'. However, it is
            % extremely necessary for more complex geometries because it
            % gives the correct horizontal component in the physical domain
            % to evaluate the function L(x).
            %-------------------------------------------------------------%
            
            x      = shape*cpi'*jacIGA((iel-1)*numbHorQuadNodes+igauss);
            
            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshape*cpi';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshape = dshape/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL RHS VECTOR
            %-------------------------------------------------------------%
            % The current code uses 'numbHorQuadNodes' instead of
            % 'numbVerQuadNodes' in the next 5 lines of code. Check if
            % the solution is consistent.
            %-------------------------------------------------------------%

            forceElem       = forceElem  + shape'*forceVec((iel-1)*numbHorQuadNodes+igauss)*gwt;    % Force
            
            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shape' * shape) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (shape' *dshape) * gwt;
            deltaElem3 = deltaElem3 + r01((iel-1)*numbHorQuadNodes+igauss) * (dshape'* shape) * gwt;
            deltaElem4 = deltaElem4 + r11((iel-1)*numbHorQuadNodes+igauss) * (dshape'*dshape) * gwt;
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:p+1
            
            i = (k-1)*(iel-1) + iel + p - (a - 1);

            if (i>1)

                globalForce(i) = globalForce(i) + forceElem(a);

                for b = 1:p+1
                    
                    j = (k-1)*(iel-1) + iel + p - (b - 1);
                    icount = icount + 1;
                    
                    row(icount)   = i;
                    col(icount)   = j;
                    
                    deltaLocal1(icount) = deltaElem1(a,b);
                    deltaLocal2(icount) = deltaElem2(a,b);
                    deltaLocal3(icount) = deltaElem3(a,b);
                    deltaLocal4(icount) = deltaElem4(a,b);
                    
                end
            else 
                icount = icount + 1;
                row(icount)  = i;
                col(icount)  = i;
                
                deltaLocal1(icount) = 0;
                deltaLocal2(icount) = 0;
                deltaLocal3(icount) = 0;
                deltaLocal4(icount) = 1;
                
            end
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPts,numbControlPts);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPts,numbControlPts);
    Delta3 = sparse(row,col,deltaLocal3,numbControlPts,numbControlPts);
    Delta4 = sparse(row,col,deltaLocal4,numbControlPts,numbControlPts);
    
    Al   = Delta1 + Delta2 + Delta3 + Delta4;
    bl   = globalForce;
    
    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%
    
    if (imb == kmb)
        Al(1,1)  = 1;
        Al(1,2)  = 0;
        Al(1,3)  = 0;
    else
        Al(1,1)  = 1e-15;
        Al(1,2)  = 0;
        Al(1,3)  = 0;
    end

    disp('Finished ASSEMBLING LOOP');
end