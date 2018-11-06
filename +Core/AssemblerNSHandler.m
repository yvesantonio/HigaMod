classdef AssemblerNSHandler   
    %% ASSEMBLER NS HANDLER CLASS
    % the AssemblerNSHandler is a class that contain all the scripts
    % responsible for the assembling and building of the block matrix that
    % describe the discretized Navier Stokes differential problem. On one hand, 
    % the classproperties, the properties of the objects were defined based on 
    % the variables needed in the implementation of each function. On the other 
    % hand, all of the previous functions are concentrated and organized in
    % the class methods.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% ASSEMBLER NS HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the AssemblerNSHandler
        % encapsulate all of the variables needed to run the methods
        % bellow, from the support methods to the final building methods.
        
        % BUILDING PROPERTIES
                
        dimModalBasis;              % Dimension of the Modal Basis
        
        leftBDomain_inX;            % Left Limit of the Domain in the X Direction
        
        rightBDomain_inX;           % Right Limit of the Domain in the X Direction
        
        wallDistance;               % Distance between Domain Walls
                      
        numbHorQuadNodes;           % Number of horizontal nodes to apply the quadrature
                                    % formula
                                    
        numbVerQuadNodes;           % Number of vertical nodes to apply the quadrature
                                    % formula
                                    
        numbFiniteElements;         % Number finite elements on the central axis                           
                                        
        inOutCondStruct;            % Data Structure used to save the Labels Identifying the Nature of
                                     % the Inflow and Outflow Boundary Conditions('dir'/'neu') 
                                    % and their Values. 

        timeStep;                    % Length of each Time Step (delta t)
        
        physicalProp;
        
        discrStruct;
        
        mapStruct;
        
        forceVector;
        
        oldFixedPointVel;
        
        oldTimeDisp;
        
        jInOut;
        
        w;
        
        currentTime;
        
        oldTimeVelWallUp;
        
        oldTimeVelWallDown;
        
        aleMap;
        
    end
    
    methods (Access = public)
        
        %% ASSEMBLER NS HANDLER - CONSTRUCT METHOD
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
        
        %% ASSEMBLER NS HANDLER - BUILDING METHODS
            
            %% Method 'buildNSwithFSI'
            
            function [DiscrStruct]= buildNSwithFSI(obj)

                %% DEBUG TUTTO DA RISCRIVERE
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

                import Core.AssemblerNSHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %% NUMBER OF NODES IN THE QUADRATURE FORMULA
                %-------------------------------------------------------------%
                % The values of 'numbHorQuadNodes' and 'numbVerQuadNodes' does 
                % not change during the computation of the Gauss-Legendre Integration 
                % Points. Attention, if you increase the number of nodes in the
                % discretization, it is necessary to refine the quadrature
                % rule.
                %-------------------------------------------------------------%

                % Horizontal direction

                NumbHorQuadNodes = obj.numbHorQuadNodes;

                % Vertical Direction

                NumbVerQuadNodes = obj.numbVerQuadNodes;

                %% NUMBER OF ELEMENTS - FINITE ELEMENTS
                %-------------------------------------------------------------%
                % Total number of elements used in the finite element mesh.
                % Also called total number of intervals along the code.
                %-------------------------------------------------------------%

                numbFE  = obj.numbFiniteElements;

                %% NUMBER OF NODES - FINITE ELEMENTS
                %-------------------------------------------------------------%
                % Total number of nodes used in the finite element mesh.
                % Pressure: we use a linear finite element discretization so 
                % nodes are the edge of the elements (P1).
                % Velocity: for each component of the velocity (u1, u2) we
                % use a discretization of second order (P2).
                %-------------------------------------------------------------%

                numbNodesP1  = numbFE+1;   
                numbNodesP2  = numbFE+numbNodesP1;
                
                %% GAUSS-LEGENDRE INTEGRATION NODES
                %-------------------------------------------------------------%
                % The following method reveives the number of integration nodes
                % in the horizontal and vertical direction and retruns the
                % respective Gauss-Legendre nodes and weigths for the
                % integration interval [0,1].
                %-------------------------------------------------------------%

                % Horizontal direction

                obj_gaussLegendre_1 = IntegrateHandler();
                obj_gaussLegendre_1.numbQuadNodes = NumbHorQuadNodes;
                [~, horGLNodesOn01, horGLWeightsOn01] = gaussLegendre(obj_gaussLegendre_1); 
                
                
                % Vertical direction

                obj_gaussLegendre_2 = IntegrateHandler();
                obj_gaussLegendre_2.numbQuadNodes = NumbVerQuadNodes;
                [~, verGLNodesOn01, verGLWeightsOn01] = gaussLegendre(obj_gaussLegendre_2);             
                
                %% RISCALED GAUSS-LEGENDRE VERTICAL NODES
                %-------------------------------------------------------------%
                % The quadrature rule riscales the Gauss-Legendre nodes and
                % weight throughout the transvelsal fiber.
                %-------------------------------------------------------------%
                
                obj_gaussLegendre_1.inputNodes = verGLNodesOn01;
                obj_gaussLegendre_1.inputWeights = verGLWeightsOn01;
                boundInterval = obj.wallDistance /2;
                obj_gaussLegendre_1.leftBoundInterval = - boundInterval;
                obj_gaussLegendre_1.rightBoundInterval = boundInterval;
                [verGLNodes, verGLWeights] = quadratureRule(obj_gaussLegendre_1);
                
                %% FEM MESH IN THE X DIRECTION
                %-------------------------------------------------------------%
                % Creation of the finite element mesh in the X direction
                % considering equispaced nodes. The mesh is created using the
                % total number of nodes and the left limit of the domain. Both
                % information come from the demo file and are passed here as a
                % property of the object.
                %-------------------------------------------------------------%

                meshFEM         = linspace( obj.leftBDomain_inX, obj.rightBDomain_inX, numbNodesP1);
                meshFEMorderP2  = linspace( obj.leftBDomain_inX, obj.rightBDomain_inX, numbNodesP2);
                
                %% AUGMENTED FEM MESH + WEIGHTS 
                %-------------------------------------------------------------%
                % Augmented finite element mesh containing all of the 
                % discretization nodes and the quadrature nodes contained in 
                % each element.
                %-------------------------------------------------------------%
                              
                augMeshFEM        = zeros( numbFE*NumbHorQuadNodes, 1);
                augMeshFEMWeights = zeros( numbFE*NumbHorQuadNodes, 1);

                %-------------------------------------------------------------%
                % Note: The loop allocates the correponding quadrature nodes
                % and weights corresponding to each finite element of the
                % discretization.
                %-------------------------------------------------------------%

                for i = 1:numbFE

                    %---------------------------------------------------------%
                    % The Gauss-Legendre nodes computed previously
                    % are rescaled to fit the interval corresponding
                    % to the current element.
                    %---------------------------------------------------------%                

                    obj_quadratureRule = IntegrateHandler();

                    obj_quadratureRule.leftBoundInterval = meshFEM(i);
                    obj_quadratureRule.rightBoundInterval = meshFEM(i+1);
                    obj_quadratureRule.inputNodes = horGLNodesOn01;
                    obj_quadratureRule.inputWeights = horGLWeightsOn01;

                    [augMeshFEM((i-1)*NumbHorQuadNodes+1 : i*NumbHorQuadNodes), ...
                     augMeshFEMWeights((i-1)*NumbHorQuadNodes+1 : i*NumbHorQuadNodes), ] = ...
                                                     quadratureRule(obj_quadratureRule);

                end
                
                %% ALL DOMAIN AUGMENTED MESH
                %-------------------------------------------------------------%
                % Command meshgrid creates two vectors Xhat and Yhat in
                % which there are the coordinates (x,y) of each mesh points
                %-------------------------------------------------------------%
                
                [InitialGridX,InitialGridY] = meshgrid( augMeshFEM,verGLNodes );
                    
                %% MODAL BASIS IN THE Y DIRECTION FOR VELOCITY AND PRESSURE
                %-------------------------------------------------------------%
                % The next method takes the coefficients of the bilinear form,
                % the dimension of the modal basis (assigned in the demo) and
                % the vertical evaluation points and computes the coefficients
                % of the modal basis when evaluated in the points of the
                % transverse fiber.
                %-------------------------------------------------------------%
                % NOTE:
                % Selecting modal basis for the transversal fiber we have
                % to distinguish different cases:
                % - Y component of velocity has robin boundary condition
                % due to FSI, so we set the number of modal basis equal to
                % modal basis dimension set in demo file;
                % - X component of velocity has homogeneous dirichlet boundary
                % condition, so we set the number of modal basis equal to
                % the modal basis dimension decreased by two;
                % - pressure has robin boundary condition, but we need to
                % descrease by two the modal basis dimension to
                % ensure numerical convergence.
                %-------------------------------------------------------------%
                
                numbModalBasisVelX = obj.dimModalBasis-2;
                numbModalBasisVelY = obj.dimModalBasis;
                numbModalBasisPres = obj.dimModalBasis-2;
                
                % Set basis for X component of the velocity
                obj_newModalBasis_1 = BasisHandler();
                obj_newModalBasis_1.evalLegendreNodes = verGLNodes;
                obj_newModalBasis_1.dimLegendreBase = numbModalBasisVelX;
                obj_newModalBasis_1.labelUpBoundCondLegendre = 'dir';
                obj_newModalBasis_1.labelDownBoundCondLegendre = 'dir';
                
                [modalBasisVelX, modalBasisDerVelX, ~] = modalBasisLegendre(obj_newModalBasis_1);
                
                % DEBUG 
                % --------------------------------- %
%                 display (modalBasisVelX);
                % --------------------------------- %
                
                % Set basis for Y component of the velocity
                obj_newModalBasis_2 = BasisHandler();
                obj_newModalBasis_2.evalLegendreNodes = verGLNodes;
                obj_newModalBasis_2.dimLegendreBase = numbModalBasisVelY;
                obj_newModalBasis_2.labelUpBoundCondLegendre = 'rob';
                obj_newModalBasis_2.labelDownBoundCondLegendre = 'rob';
                
                [modalBasisVelY, modalBasisDerVelY, ~] = modalBasisLegendre(obj_newModalBasis_2);
                
                % Set basis for pressure
                obj_newModalBasis_3 = BasisHandler();
                obj_newModalBasis_3.evalLegendreNodes = verGLNodes;
                obj_newModalBasis_3.dimLegendreBase = numbModalBasisPres;
                obj_newModalBasis_3.labelUpBoundCondLegendre = 'rob';
                obj_newModalBasis_3.labelDownBoundCondLegendre = 'rob';
                
                [modalBasisPres, ~, ~] = modalBasisLegendre(obj_newModalBasis_3);
                
                %% MODAL BASIS ON DOMAIN WALL
                %-------------------------------------------------------------%
                % Buondary modal basis is needed in order to imposed robin
                % boundary condition on the domain wall. 
                %-------------------------------------------------------------%
                
                obj_newModalBasis_4 = BasisHandler();
                boundInterval = obj.wallDistance /2;
                obj_newModalBasis_4.evalLegendreNodes = [ -boundInterval, boundInterval ];
                obj_newModalBasis_4.dimLegendreBase = numbModalBasisVelY;
                obj_newModalBasis_4.labelUpBoundCondLegendre = 'rob';
                obj_newModalBasis_4.labelDownBoundCondLegendre = 'rob';
                
                [modalBasisWall, modalBasisDerWall, modalBasisDer2Wall] = modalBasisLegendre(obj_newModalBasis_4);
                                    % new
                                    
                %% FEM BASIS IN THE X DIRECTION FOR VELOCITY AND PRESSURE
                %-------------------------------------------------------------%
                % We impose Neumann condition in inflow and outflow
                % boundary, for both velocity in x direction and velocity
                % in y direction.
                %-------------------------------------------------------------%
                obj_newFiniteElementBasis_1 = BasisHandler();
                obj_newFiniteElementBasis_1.meshNodesX = meshFEM(1:2);
                obj_newFiniteElementBasis_1.meshQuadratureNodesX = augMeshFEM(1:NumbHorQuadNodes);
                obj_newFiniteElementBasis_1.degreeFiniteElement = 2;
                [basisFEMVel, basisFEMDerVel, basisFEMDer2Vel] = finiteElementBasis(obj_newFiniteElementBasis_1);
                
                obj_newFiniteElementBasis_2 = BasisHandler();
                obj_newFiniteElementBasis_2.meshNodesX = meshFEM(1:2);
                obj_newFiniteElementBasis_2.meshQuadratureNodesX = augMeshFEM(1:NumbHorQuadNodes);
                obj_newFiniteElementBasis_2.degreeFiniteElement = 1;
                [basisFEMPres, ~, ~] = finiteElementBasis(obj_newFiniteElementBasis_2);
                
                %% DISCRETIZATION STRUCTURE
              
                DiscrStruct = struct( ...
                    'numbNodesP1',       numbNodesP1,        ...
                    'numbNodesP2',       numbNodesP2,        ...
                    'verGLNodes',        verGLNodes,         ...
                    'verGLWeights',      verGLWeights,       ...
                    'meshFEM',           meshFEM,            ...
                    'meshFEMorderP2',    meshFEMorderP2,     ...
                    'augMeshFEM',        augMeshFEM,         ...
                    'augMeshFEMWeights', augMeshFEMWeights,  ...
                    'initialGridX',      InitialGridX,       ...
                    'initialGridY',      InitialGridY,       ...
                    'numbModalBasisVelX',numbModalBasisVelX, ...
                    'modalBasisVelX',    modalBasisVelX,     ...
                    'modalBasisDerVelX', modalBasisDerVelX,  ...
                    'numbModalBasisVelY',numbModalBasisVelY, ...
                    'modalBasisVelY',    modalBasisVelY,     ...
                    'modalBasisDerVelY', modalBasisDerVelY,  ...
                    'numbModalBasisPres',numbModalBasisPres, ...
                    'modalBasisPres',    modalBasisPres,     ...
                    'modalBasisWall',    modalBasisWall,     ...
                    'modalBasisDerWall', modalBasisDerWall,  ... % new
                    'modalBasisDer2Wall',modalBasisDer2Wall, ... % new
                    'basisFEMVel',       basisFEMVel,        ...
                    'basisFEMDerVel',    basisFEMDerVel,     ...
                    'basisFEMDer2Vel',   basisFEMDer2Vel,    ... % new
                    'basisFEMPres',      basisFEMPres );
                                
            end

            %% Method 'assemblerNSwithFSI'
            
            function [A,b]= assemblerNSwithFSI(obj)

                %% DEBUG TUTTO DA RISCRIVERE
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

                import Core.AssemblerNSHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler

                %% NUMBER OF NODES IN THE QUADRATURE FORMULA
                %-------------------------------------------------------------%
                % The values of 'numbHorQuadNodes' and 'numbVerQuadNodes' does 
                % not change during the computation of the Gauss-Legendre Integration 
                % Points. Attention, if you increase the number of nodes in the
                % discretization, it is necessary to refine the quadrature
                % rule.
                %-------------------------------------------------------------%

                % Horizontal direction

                NumbHorQuadNodes = obj.numbHorQuadNodes;

                %% NUMBER OF ELEMENTS - FINITE ELEMENTS
                %-------------------------------------------------------------%
                % Total number of elements used in the finite element mesh.
                % Also called total number of intervals along the code.
                %-------------------------------------------------------------%

                numbFE  = obj.numbFiniteElements;

                %% MESH AND BASIS
                
                numbNodesP1         = obj.discrStruct.numbNodesP1;
                numbNodesP2         = obj.discrStruct.numbNodesP2;               
                verGLWeights        = obj.discrStruct.verGLWeights;
                augMeshFEM          = obj.discrStruct.augMeshFEM; % new
                augMeshFEMWeights   = obj.discrStruct.augMeshFEMWeights;   
                numbModalBasisVelX  = obj.discrStruct.numbModalBasisVelX;
                modalBasisVelX      = obj.discrStruct.modalBasisVelX;         
                modalBasisDerVelX   = obj.discrStruct.modalBasisDerVelX;
                numbModalBasisVelY  = obj.discrStruct.numbModalBasisVelY;
                modalBasisVelY      = obj.discrStruct.modalBasisVelY;         
                modalBasisDerVelY   = obj.discrStruct.modalBasisDerVelY;
                numbModalBasisPres  = obj.discrStruct.numbModalBasisPres;
                modalBasisPres      = obj.discrStruct.modalBasisPres;         
                modalBasisWall      = obj.discrStruct.modalBasisWall;
                modalBasisDerWall   = obj.discrStruct.modalBasisDerWall; % new
                modalBasisDer2Wall  = obj.discrStruct.modalBasisDer2Wall; % new
                basisFEMVel         = obj.discrStruct.basisFEMVel;               
                basisFEMDerVel      = obj.discrStruct.basisFEMDerVel;         
                basisFEMDer2Vel     = obj.discrStruct.basisFEMDerVel; % new
                basisFEMPres        = obj.discrStruct.basisFEMPres;
                
                %% INITILIZATION PARAMETERS
                
                Map = obj.mapStruct;
                Fx = obj.forceVector.Fx;
                Fy = obj.forceVector.Fy;
                U_kmeno1_x = obj.oldFixedPointVel.Ux;
                U_kmeno1_y = obj.oldFixedPointVel.Uy;
                Etaup_old = obj.oldTimeDisp.EtaUp;
                Etadown_old = obj.oldTimeDisp.EtaDown;
                % new
                EtaupDer_old = obj.oldTimeDisp.EtaUpDer_old;
                EtadownDer_old = obj.oldTimeDisp.EtaDownDer_old;
                EtaupDer2_old = obj.oldTimeDisp.EtaUpDer2_old;
                EtadownDer2_old = obj.oldTimeDisp.EtaDownDer2_old;
                % end new
                W = obj.w;
                
                %% MEMORY ALLOCATION FOR SYSTEM MATRICES
                %-------------------------------------------------------------%
                % The final algebrical system A*solution = b is
                % characterized by a NxN matrix A, with: 
                % N = numbDofVelX + numbDofVelY + numbDofPres.
                %
                % This matrix is built using submatrices such as:
                % - Ax(Ay) ( represent the operator that associates the x(y)
                %       component of the velocity to the x(y) component of the 
                %       velocity test function );
                % - Ayx(Axy) ( represent the operator that associates the
                %       y(x) component of the velocity to the x(y) component 
                %       of the velocity test function );
                % - Bx(By) ( represent the operator that associates the
                %       x(y) component of the velocity to the pressure test 
                %       function );
                %-------------------------------------------------------------%
                
                numbDofVelX=numbModalBasisVelX*numbNodesP2;
                numbDofVelY=numbModalBasisVelY*numbNodesP2;
                numbDofPres=numbModalBasisPres*numbNodesP1;
                
                Ax = sparse( numbDofVelX, numbDofVelX );
                Ay = sparse( numbDofVelY, numbDofVelY );
        
                Axy = sparse( numbDofVelY, numbDofVelX );
                Ayx = sparse( numbDofVelX, numbDofVelY );
        
                Bx = sparse( numbDofPres,numbDofVelX );
                By = sparse( numbDofPres,numbDofVelY );
        
                bx = zeros( numbDofVelX ,1);
                by = zeros( numbDofVelY ,1);
                
                %% ASSEMBLER LOOPS
                
                % Matrix Ax
                for k = 1:numbModalBasisVelX
                    MbKX=modalBasisVelX(:,k);
                    DMbKX=modalBasisDerVelX(:,k);
                    for r = 1:numbModalBasisVelX
                        MbRX=modalBasisVelX(:,r);
                        DMbRX=modalBasisDerVelX(:,r);
                        
                        [ Ax(1+(r-1)*numbNodesP2 : r*numbNodesP2 , 1+(k-1)*numbNodesP2 : k*numbNodesP2), ...
                            bx(1+(r-1)*numbNodesP2 : r*numbNodesP2 ) ]= ...
                            AssembleVel1DX(Map,verGLWeights, MbKX, MbRX,DMbKX,DMbRX,numbFE,...
                            augMeshFEMWeights,NumbHorQuadNodes,basisFEMVel',basisFEMDerVel',1/obj.timeStep,...
                            obj.physicalProp.nu,Fx,U_kmeno1_x,U_kmeno1_y-W);
                    end
                end
                
                % DEBUG 
                % ------------------------------ %
%                 display (Ax);
                % ------------------------------ %

                % Matrix Axy
                for k = 1:numbModalBasisVelX
                    DMbKX=modalBasisDerVelX(:,k);
                    for r = 1:numbModalBasisVelY
                        MbRY=modalBasisVelY(:,r);
                        DMbRY=modalBasisDerVelY(:,r);

                        Axy(1+(r-1)*numbNodesP2 : r*numbNodesP2 , 1+(k-1)*numbNodesP2 : k*numbNodesP2)=...
                            AssembleVelXY(Map,verGLWeights,...
                            DMbKX,MbRY,DMbRY,...
                            basisFEMVel',basisFEMDerVel',numbFE,...
                            augMeshFEMWeights,NumbHorQuadNodes,obj.physicalProp.nu);
                    end
                end

                % DEBUG 
                % ------------------------------ %
%                 display (Axy);
                % ------------------------------ %

                % Matrix Ay
                for k = 1:numbModalBasisVelY
                    MbKY=modalBasisVelY(:,k);
                    DMbKY=modalBasisDerVelY(:,k);

                    MBrobKY=modalBasisWall(:,k);
                    DMBrobKY=modalBasisDerWall(:,k); % new
                    D2MBrobKY=modalBasisDer2Wall(:,k);

                    for r = 1:numbModalBasisVelY
                        MBrobRY=modalBasisWall(:,r);
                        DMBrobRY=modalBasisDerWall(:,r); % new
                        D2MBrobRY=modalBasisDer2Wall(:,r);

                        MbRY=modalBasisVelY(:,r);
                        DMbRY=modalBasisDerVelY(:,r);

                        [ Ay(1+(r-1)*numbNodesP2 : r*numbNodesP2 , 1+(k-1)*numbNodesP2 : k*numbNodesP2), ...
                            by(1+(r-1)*numbNodesP2 : r*numbNodesP2 ) ]= ...
                            FSIAssembleVel1DY(Map,verGLWeights, MbKY, MbRY,DMbKY,DMbRY,numbFE,...
                            augMeshFEMWeights,NumbHorQuadNodes,basisFEMVel',basisFEMDerVel',basisFEMDer2Vel',obj.timeStep,obj.physicalProp.nu,...
                            Fy,U_kmeno1_x,U_kmeno1_y-W,obj.physicalProp.K,MBrobKY,DMBrobKY,D2MBrobKY,MBrobRY,DMBrobRY,D2MBrobRY,Etadown_old,Etaup_old,EtadownDer_old,EtaupDer_old,EtadownDer2_old,EtaupDer2_old, ... % new DMBrobKY, DMBrobRY, EtadownDer_old, EtaupDer_old
                            obj.physicalProp.rhoStructure,obj.physicalProp.hStructure,obj.oldTimeVelWallUp,obj.oldTimeVelWallDown,...
                            obj.physicalProp.a,obj.physicalProp.b,obj.physicalProp.c,obj.aleMap);
                        
%                         [ Ay(1+(r-1)*numbNodesP2 : r*numbNodesP2 , 1+(k-1)*numbNodesP2 : k*numbNodesP2), ...
%                             by(1+(r-1)*numbNodesP2 : r*numbNodesP2 ) ]= ...
%                             FSIAssembleVel1DYold(Map,verGLWeights, MbKY, MbRY,DMbKY,DMbRY,numbFE,...
%                             augMeshFEMWeights,NumbHorQuadNodes,basisFEMVel',basisFEMDerVel',obj.timeStep,obj.physicalProp.nu,...
%                             Fy,U_kmeno1_x,U_kmeno1_y-W,obj.physicalProp.K,MBrobKY,MBrobRY,Etadown_old,Etaup_old, ...
%                             obj.physicalProp.rhoStructure,obj.physicalProp.hStructure,obj.oldTimeVelWallUp,obj.oldTimeVelWallDown,...
%                             obj.physicalProp.b,obj.physicalProp.c,obj.aleMap);
                    end
                end

                % DEBUG 
                % ------------------------------ %
%                 display (Ay);
                % ------------------------------ %

                % Matrix Ayx
                for k = 1:numbModalBasisVelY
                    MbKY=modalBasisVelY(:,k);
                    DMbKY=modalBasisDerVelY(:,k);

                    for r = 1:numbModalBasisVelX
                        DMbRX=modalBasisDerVelX(:,r);

                        Ayx(1+(r-1)*numbNodesP2 : r*numbNodesP2 , 1+(k-1)*numbNodesP2 : k*numbNodesP2)=...
                            AssembleVelYX(Map,verGLWeights,...
                            MbKY,DMbKY,DMbRX,...
                            basisFEMVel',basisFEMDerVel',numbFE,...
                            augMeshFEMWeights,NumbHorQuadNodes,obj.physicalProp.nu);
                    end
                end

                % DEBUG 
                % ------------------------------ %
%                 display (Ayx);
                % ------------------------------ %

                % Matrix Bx
                for k = 1:numbModalBasisVelX
                    MbKX=modalBasisVelX(:,k);
                    DMbKX=modalBasisDerVelX(:,k);

                    for s=1:numbModalBasisPres

                        MbSP=modalBasisPres(:,s);

                        Bx(1+(s-1)*numbNodesP1:s*numbNodesP1, 1+(k-1)*numbNodesP2 : k*numbNodesP2) = ...
                            AssembleBX(Map,verGLWeights,MbKX,MbSP,DMbKX,numbFE,augMeshFEMWeights,NumbHorQuadNodes,basisFEMVel',basisFEMDerVel',basisFEMPres');
                    end
                end
                
                % DEBUG 
                % ------------------------------ %
%                 display (Bx);
                % ------------------------------ %

                % Matrix By
                for k = 1:numbModalBasisVelY

                    DMbKY=modalBasisDerVelY(:,k);

                    for s=1:numbModalBasisPres

                        MbSP=modalBasisPres(:,s);

                        By(1+(s-1)*numbNodesP1:s*numbNodesP1, 1+(k-1)*numbNodesP2 : k*numbNodesP2) = ...
                            AssembleBY(Map,verGLWeights,MbSP,DMbKY,numbFE,augMeshFEMWeights,NumbHorQuadNodes,basisFEMVel',basisFEMPres');
                    end
                end

                % DEBUG 
                % ------------------------------ %
%                 display (By);
                % ------------------------------ %
                
                %% BOUNDARY CONDITION IMPOSITION
                
                % Definition of the Object from the BoundaryConditionHandler class
                
                obj_BC = BoundaryConditionHandler();
                
                % Properties Assignment
                
                obj_BC.matAx = Ax;
                obj_BC.matAy = Ay;
                obj_BC.vectBx = bx;
                obj_BC.vectBy = by;
                obj_BC.boundaryCondStruct = obj.inOutCondStruct;
                obj_BC.numbModalBasisVelX = numbModalBasisVelX;
                obj_BC.numbModalBasisVelY = numbModalBasisVelY;
                obj_BC.numbNodesFEM = numbNodesP2;
                obj_BC.verGLWeights = verGLWeights;
                obj_BC.modalBasisVelX = modalBasisVelX;
                obj_BC.jInOut = (obj.jInOut)';
                obj_BC.currentTime = obj.currentTime;
                
                % Call of the 'imposeBCInflowOutflowNS' Method

                [Ax,Ay,bx,by] = imposeBCInflowOutflowNS (obj_BC);

                %% ASSEMBLING END 
                        
                A = [Ax,Ayx,Bx';Axy,Ay,By';Bx,By,zeros(numbDofPres,numbDofPres)];
                b = [bx;by;zeros(numbDofPres,1)];

            end
            
    end
end

%% Method 'AssembleVel1DX'

function [Al,bl]=AssembleVel1DX(Map,wyq,mb_k,mb_r,mb_yk,mb_yr,ne,mesh_wx,nqnx,FEM_basis,FEM_basis_x,usudt,nu,Force,Bx,By)

% per ogni x (vettorialmente)
% calcola r^st_kr
r11v =(2*nu*Map.Jac)*( mb_k.*mb_r.*wyq);

r10v = ( 2*nu*Map.Jac.*Map.D)*( mb_k.*mb_yr .*wyq ) ...
     + (    Map.Jac.*Bx)*( mb_k  .*mb_r .*wyq );

r01v = ( 2*nu*Map.Jac.*Map.D)*( mb_yk.*mb_r .*wyq );

r00v = ( nu*Map.Jac.*(2*Map.D.^2+Map.J.^2) ) *(mb_yk.*mb_yr.*wyq )...
     + (    Map.Jac*usudt ) *(mb_k.*mb_r.*wyq ) ...
     + (    Map.Jac.*(Map.D.*Bx+Map.J.*By))*(mb_yk .*mb_r .*wyq);
 
forcev = (  Map.Jac.*Force)*( mb_r.*wyq);
nx=ne+1;

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag =   zeros( ne + nx, 1);
diagsup  =   zeros( ne + nx, 1); % Spdiags taglierà il primo elemento
diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglierà l'ultimo
diagsuper  = zeros( ne + nx, 1); % Spdiags taglierà i primi due elementi
diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglierà gli ultimi due
bl         = zeros( ne + nx, 1);
for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficienti relativi all'intervallo in questione
    r11 = r11v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r01 = r01v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:3                  % Per ogni funzione di forma EF test
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basis  ( i, :);
        femb_xi = FEM_basis_x( i, :);
        for j = 1:3             % Per ogni funzione di forma EF sol
            
            femb_j  = FEM_basis  ( j, :);
            femb_xj = FEM_basis_x( j, :);
            
            if(i==j)
                maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+...
                    (  r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j+1 && i==2)
                diaginf(2*ie-1)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+1 && i==3)
                diaginf(2*ie)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+2)
                diaginfer(2*ie-1)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-2)
                diagsuper(2*ie+1)=(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-1 && i==2)
                diagsup(2*ie+1)=( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j-1 && i==1)
                diagsup(2*ie)=( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            end
        end
        bl(2*(ie-1)+i) = bl(2*(ie-1)+i) + ( forcev((ie-1)*nqnx+1:ie*nqnx)'.*femb_i )*w;
    end
end
Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A
end

%% Method 'FSIAssembleVel1DY'

function [Al,bl]=FSIAssembleVel1DY(Map,wyq,mb_k,mb_r,mb_yk,mb_yr,ne,...
    mesh_wx,nqnx,FEM_basis,FEM_basis_x,FEM_basis_xx,dt,nu,Force,Bx,By,...
    E,MBrobKY,DMBrobKY,D2MBrobKY,MBrobRY,DMBrobRY,D2MBrobRY,etadown_old,etaup_old,etadownDer_old,etaupDer_old,etadownDer2_old,etaupDer2_old, ... % new DMBrobKY, DMBrobRY, etadownDer_old, etaupDer_old
    rho_s,h_s,uup_old,udown_old,a,b,c,AleMap)

% per ogni x (vettorialmente)
% calcola r^st_kr       I SEGNI SONO CORRETTI ????????????
r11v = (2*nu*Map.Jac)*(mb_k.*mb_r.*wyq)... % mancava 2*
     + (b*dt+c)*AleMap.JacUp.*Map.Aup*MBrobKY(2)*MBrobRY(2) ...
     + (b*dt+c)*AleMap.JacDown.*Map.Adown*MBrobKY(1)*MBrobRY(1) ...
     + (4*a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^2)*DMBrobKY(2)*DMBrobRY(2) ...
     + (4*a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^2)*DMBrobKY(1)*DMBrobRY(1);
r10v = ( 2*nu*Map.Jac.*Map.D)*( mb_k.*mb_yr .*wyq ) ... % mancava 2*
     + (    Map.Jac.*Bx)*( mb_k  .*mb_r .*wyq ) ...
     + (b*dt+c)*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp*MBrobKY(2)*DMBrobRY(2) ... % new
     + (b*dt+c)*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown*MBrobKY(1)*DMBrobRY(1) ... % new
     + (2*a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^3)*DMBrobKY(2)*D2MBrobRY(2) ...
     + (2*a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^3)*DMBrobKY(1)*D2MBrobRY(1) ...
     + (2*a*dt)*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp.*Map.LambdaDer2Up*DMBrobKY(2)*DMBrobRY(2) ...
     + (2*a*dt)*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown.*Map.LambdaDer2Down*DMBrobKY(1)*DMBrobRY(1);
r01v = ( 2*nu*Map.Jac.*Map.D)*( mb_yk.*mb_r .*wyq ) ... % mancava 2*
     + (b*dt+c)*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp*DMBrobKY(2)*MBrobRY(2) ... % new
     + (b*dt+c)*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown*DMBrobKY(1)*MBrobRY(1) ... % new
     + (2*a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^3)*D2MBrobKY(2)*DMBrobRY(2) ...
     + (2*a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^3)*D2MBrobKY(1)*DMBrobRY(1) ...
     + (2*a*dt)*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp.*Map.LambdaDer2Up*DMBrobKY(2)*DMBrobRY(2) ...
     + (2*a*dt)*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown.*Map.LambdaDer2Down*DMBrobKY(1)*DMBrobRY(1);
r00v = (nu*Map.Jac.*(Map.D.^2+2*Map.J.^2) ) *(mb_yk.*mb_yr.*wyq ) ...
     + (    Map.Jac/dt ) *(mb_k.*mb_r.*wyq ) ...
     + (    Map.Jac.*(Map.D.*Bx+Map.J.*By))*(mb_yk .*mb_r .*wyq) ...
	 + (rho_s*h_s/dt+E*dt)*AleMap.JacUp.*Map.Aup*MBrobKY(2)*MBrobRY(2) ...
     + (rho_s*h_s/dt+E*dt)*AleMap.JacDown.*Map.Adown*MBrobKY(1)*MBrobRY(1) ... 
     + (b*dt+c)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^2)*DMBrobKY(2)*DMBrobRY(2) ... % new
     + (b*dt+c)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^2)*DMBrobKY(1)*DMBrobRY(1) ... % new
     + (a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^4)*D2MBrobKY(2)*D2MBrobRY(2) ...
     + (a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^4)*D2MBrobKY(1)*D2MBrobRY(1) ...
     + (a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^2).*Map.LambdaDer2Up*D2MBrobKY(2)*DMBrobRY(2) ...
     + (a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^2).*Map.LambdaDer2Down*D2MBrobKY(1)*DMBrobRY(1) ...
     + (a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^2).*Map.LambdaDer2Up*DMBrobKY(2)*D2MBrobRY(2) ...
     + (a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^2).*Map.LambdaDer2Down*DMBrobKY(1)*D2MBrobRY(1) ...
     + (a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDer2Up.^2)*DMBrobKY(2)*DMBrobRY(2) ...
     + (a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDer2Down.^2)*DMBrobKY(1)*DMBrobRY(1);
% new
r22v = + (a*dt)*AleMap.JacUp.*Map.Aup*MBrobKY(2)*MBrobRY(2) ...
       + (a*dt)*AleMap.JacDown.*Map.Adown*MBrobKY(1)*MBrobRY(1);
r20v = + (a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^2)*MBrobKY(2)*D2MBrobRY(2) ...
       + (a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^2)*MBrobKY(1)*D2MBrobRY(1) ...
       + (a*dt)*AleMap.JacUp.*Map.Aup.*Map.LambdaDer2Up*MBrobKY(2)*DMBrobRY(2) ...
       + (a*dt)*AleMap.JacDown.*Map.Adown.*Map.LambdaDer2Down*MBrobKY(1)*DMBrobRY(1);
r21v = + (2*a*dt)*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp*MBrobKY(2)*DMBrobRY(2) ...
       + (2*a*dt)*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown*MBrobKY(1)*DMBrobRY(1);
r12v = + (2*a*dt)*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp*DMBrobKY(2)*MBrobRY(2) ...
       + (2*a*dt)*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown*DMBrobKY(1)*MBrobRY(1);
r02v = + (a*dt)*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^2)*D2MBrobKY(2)*MBrobRY(2) ...
       + (a*dt)*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^2)*D2MBrobKY(1)*MBrobRY(1) ...
       + (a*dt)*AleMap.JacUp.*Map.Aup.*Map.LambdaDer2Up*DMBrobKY(2)*MBrobRY(2) ...
       + (a*dt)*AleMap.JacDown.*Map.Adown.*Map.LambdaDer2Down*DMBrobKY(1)*MBrobRY(1);
b21inv = - (a*dt)*MBrobKY(2)*MBrobRY(2) ...
         - (a*dt)*MBrobKY(1)*MBrobRY(1) *0;
b21outv = - (a*dt)*MBrobKY(2)*MBrobRY(2) ...
          - (a*dt)*MBrobKY(1)*MBrobRY(1) *0;
b20inv = - (a*dt)*Map.LambdaDerUp(1)*MBrobKY(2)*DMBrobRY(2) ...
         - (a*dt)*Map.LambdaDerDown(1)*MBrobKY(1)*DMBrobRY(1);
b20outv = - (a*dt)*Map.LambdaDerUp(end)*MBrobKY(2)*DMBrobRY(2) ...
          - (a*dt)*Map.LambdaDerDown(end)*MBrobKY(1)*DMBrobRY(1);
b11inv = - (2*a*dt)*Map.LambdaDerUp(1)*DMBrobKY(2)*MBrobRY(2) ...
         - (2*a*dt)*Map.LambdaDerDown(1)*DMBrobKY(1)*MBrobRY(1);
b11outv = - (2*a*dt)*Map.LambdaDerUp(end)*DMBrobKY(2)*MBrobRY(2) ...
          - (2*a*dt)*Map.LambdaDerDown(end)*DMBrobKY(1)*MBrobRY(1);
b10inv = - (2*a*dt)*(Map.LambdaDerUp(1)^2)*DMBrobKY(2)*DMBrobRY(2) ...
         - (2*a*dt)*(Map.LambdaDerDown(1)^2)*DMBrobKY(1)*DMBrobRY(1);
b10outv = - (2*a*dt)*(Map.LambdaDerUp(end)^2)*DMBrobKY(2)*DMBrobRY(2) ...
          - (2*a*dt)*(Map.LambdaDerDown(end)^2)*DMBrobKY(1)*DMBrobRY(1);
b01inv = - (a*dt)*(Map.LambdaDerUp(1)^2)*D2MBrobKY(2)*MBrobRY(2) ...
         - (a*dt)*(Map.LambdaDerDown(1)^2)*D2MBrobKY(1)*MBrobRY(1) ...
         - (a*dt)*Map.LambdaDer2Up(1)*DMBrobKY(2)*MBrobRY(2) ...
         - (a*dt)*Map.LambdaDer2Down(1)*DMBrobKY(1)*MBrobRY(1);
b01outv = - (a*dt)*(Map.LambdaDerUp(end)^2)*D2MBrobKY(2)*MBrobRY(2) ...
          - (a*dt)*(Map.LambdaDerDown(end)^2)*D2MBrobKY(1)*MBrobRY(1) ...
          - (a*dt)*Map.LambdaDer2Up(end)*DMBrobKY(2)*MBrobRY(2) ...
          - (a*dt)*Map.LambdaDer2Down(end)*DMBrobKY(1)*MBrobRY(1);
b00inv = - (a*dt)*(Map.LambdaDerUp(1)^3)*D2MBrobKY(2)*DMBrobRY(2) ...
         - (a*dt)*(Map.LambdaDerDown(1)^3)*D2MBrobKY(1)*DMBrobRY(1) ...
         - (a*dt)*Map.LambdaDer2Up(1)*Map.LambdaDerUp(1)*DMBrobKY(2)*DMBrobRY(2) ...
         - (a*dt)*Map.LambdaDer2Down(1)*Map.LambdaDerDown(1)*DMBrobKY(1)*DMBrobRY(1);
b00outv = - (a*dt)*(Map.LambdaDerUp(end)^3)*D2MBrobKY(2)*DMBrobRY(2) ...
          - (a*dt)*(Map.LambdaDerDown(end)^3)*D2MBrobKY(1)*DMBrobRY(1) ...
          - (a*dt)*Map.LambdaDer2Up(end)*Map.LambdaDerUp(end)*DMBrobKY(2)*DMBrobRY(2) ...
          - (a*dt)*Map.LambdaDer2Down(end)*Map.LambdaDerDown(end)*DMBrobKY(1)*DMBrobRY(1);
 
force0v = ( Map.Jac.*Force )*( mb_r.*wyq)...
        - E*etaup_old.*AleMap.JacUp.*Map.Aup*MBrobRY(2) ...
        - E*etadown_old.*AleMap.JacDown.*Map.Adown*MBrobRY(1)...
        + rho_s*h_s/dt*uup_old.*AleMap.JacUp.*Map.Aup*MBrobRY(2) ...
        + rho_s*h_s/dt*udown_old.*AleMap.JacDown.*Map.Adown*MBrobRY(1)...
        - b*etaupDer_old.*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp*DMBrobRY(2) ... % new
        - b*etadownDer_old.*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown*DMBrobRY(1) ... % new
        - a*etaupDer2_old.*AleMap.JacUp.*Map.Aup.*(Map.LambdaDerUp.^2)*D2MBrobRY(2) ...
        - a*etadownDer2_old.*AleMap.JacDown.*Map.Adown.*(Map.LambdaDerDown.^2)*D2MBrobRY(1) ...
        - a*etaupDer2_old.*AleMap.JacUp.*Map.Aup.*Map.LambdaDer2Up*DMBrobRY(2) ...
        - a*etadownDer2_old.*AleMap.JacDown.*Map.Adown.*Map.LambdaDer2Down*DMBrobRY(1);
force1v = - b*etaupDer_old.*AleMap.JacUp.*Map.Aup*MBrobRY(2) ... % new
          - b*etadownDer_old.*AleMap.JacDown.*Map.Adown*MBrobRY(1) ... % new
          - 2*a*etaupDer2_old.*AleMap.JacUp.*Map.Aup.*Map.LambdaDerUp*DMBrobRY(2) ...
          - 2*a*etadownDer2_old.*AleMap.JacDown.*Map.Adown.*Map.LambdaDerDown*DMBrobRY(1);
% new
force2v = - a*etaupDer2_old.*AleMap.JacUp.*Map.Aup*DMBrobRY(2) ...
          - a*etadownDer2_old.*AleMap.JacDown.*Map.Adown*DMBrobRY(1);
force1inv = + a*etaupDer2_old(1)*MBrobRY(2) ...
            + a*etadownDer2_old(1)*MBrobRY(1) *0; % new
force1outv = + a*etaupDer2_old(end)*MBrobRY(2) ...
             + a*etadownDer2_old(end)*MBrobRY(1) *0; % new
force0inv = + a*etaupDer2_old(1)*Map.LambdaDerUp(1)*DMBrobRY(2) ...
            + a*etadownDer2_old(1)*Map.LambdaDerDown(1)*DMBrobRY(1);
force0outv = + a*etaupDer2_old(end)*Map.LambdaDerUp(end)*DMBrobRY(2) ...
             + a*etadownDer2_old(end)*Map.LambdaDerDown(end)*DMBrobRY(1);
    
nx=ne+1;

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag =   zeros( ne + nx, 1);
diagsup  =   zeros( ne + nx, 1); % Spdiags taglierà il primo elemento
diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglierà l'ultimo
diagsuper  = zeros( ne + nx, 1); % Spdiags taglierà i primi due elementi
diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglierà gli ultimi due
bl         = zeros( ne + nx, 1);
for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r11 = r11v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r01 = r01v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r22 = r22v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r20 = r20v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r21 = r21v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r12 = r12v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r02 = r02v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:3                  % Per ogni funzione di forma EF test
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i   = FEM_basis  ( i, :);
        femb_xi  = FEM_basis_x( i, :);
        femb_xxi = FEM_basis_xx( i, :);
        for j = 1:3             % Per ogni funzione di forma EF sol
            
            femb_j   = FEM_basis  ( j, :);
            femb_xj  = FEM_basis_x( j, :);
            femb_xxj = FEM_basis_xx( j, :);
            
            if(i==j)
                maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    + r22.*femb_xxi.*femb_xxj ...
                    + r20.*femb_i.*femb_xxj ...
                    + r21.*femb_xi.*femb_xxj ...
                    + r12.*femb_xxi.*femb_xj ...
                    + r02.*femb_xxi.*femb_j ...
                    )*w;
                if(ie==1)
                    maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+(...
                      b21inv*femb_xi(1)*femb_xxj(1) ...
                    + b20inv*femb_i(1)*femb_xxj(1) ...
                    + b11inv*femb_xi(1)*femb_xj(1) ...
                    + b10inv*femb_i(1)*femb_xj(1) ...
                    + b01inv*femb_xi(1)*femb_j(1) ...
                    + b00inv*femb_i(1)*femb_j(1) );
                end
                if(ie==ne)
                    maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+(...
                      b21outv*femb_xi(end)*femb_xxj(end) ...
                    + b20outv*femb_i(end)*femb_xxj(end) ...
                    + b11outv*femb_xi(end)*femb_xj(end) ...
                    + b10outv*femb_i(end)*femb_xj(end) ...
                    + b01outv*femb_xi(end)*femb_j(end) ...
                    + b00outv*femb_i(end)*femb_j(end) );
                end
            elseif(i==j+1 && i==2)
                diaginf(2*ie-1)=diaginf(2*ie-1)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    + r22.*femb_xxi.*femb_xxj ...
                    + r20.*femb_i.*femb_xxj ...
                    + r21.*femb_xi.*femb_xxj ...
                    + r12.*femb_xxi.*femb_xj ...
                    + r02.*femb_xxi.*femb_j ...
                    )*w;
                if(ie==1)
                    diaginf(2*ie-1)=diaginf(2*ie-1)+(...
                      b21inv*femb_xi(1)*femb_xxj(1) ...
                    + b20inv*femb_i(1)*femb_xxj(1) ...
                    + b11inv*femb_xi(1)*femb_xj(1) ...
                    + b10inv*femb_i(1)*femb_xj(1) ...
                    + b01inv*femb_xi(1)*femb_j(1) ...
                    + b00inv*femb_i(1)*femb_j(1) );
                end
                if(ie==ne)
                    diaginf(2*ie-1)=diaginf(2*ie-1)+(...
                      b21outv*femb_xi(end)*femb_xxj(end) ...
                    + b20outv*femb_i(end)*femb_xxj(end) ...
                    + b11outv*femb_xi(end)*femb_xj(end) ...
                    + b10outv*femb_i(end)*femb_xj(end) ...
                    + b01outv*femb_xi(end)*femb_j(end) ...
                    + b00outv*femb_i(end)*femb_j(end) );
                end
            elseif(i==j+1 && i==3)
                diaginf(2*ie)=diaginf(2*ie)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    + r22.*femb_xxi.*femb_xxj ...
                    + r20.*femb_i.*femb_xxj ...
                    + r21.*femb_xi.*femb_xxj ...
                    + r12.*femb_xxi.*femb_xj ...
                    + r02.*femb_xxi.*femb_j ...
                    )*w;
                if(ie==1)
                    diaginf(2*ie)=diaginf(2*ie)+(...
                      b21inv*femb_xi(1)*femb_xxj(1) ...
                    + b20inv*femb_i(1)*femb_xxj(1) ...
                    + b11inv*femb_xi(1)*femb_xj(1) ...
                    + b10inv*femb_i(1)*femb_xj(1) ...
                    + b01inv*femb_xi(1)*femb_j(1) ...
                    + b00inv*femb_i(1)*femb_j(1) );
                end
                if(ie==ne)
                    diaginf(2*ie)=diaginf(2*ie)+(...
                      b21outv*femb_xi(end)*femb_xxj(end) ...
                    + b20outv*femb_i(end)*femb_xxj(end) ...
                    + b11outv*femb_xi(end)*femb_xj(end) ...
                    + b10outv*femb_i(end)*femb_xj(end) ...
                    + b01outv*femb_xi(end)*femb_j(end) ...
                    + b00outv*femb_i(end)*femb_j(end) );
                end
            elseif(i==j+2)
                diaginfer(2*ie-1)=diaginfer(2*ie-1)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    + r22.*femb_xxi.*femb_xxj ...
                    + r20.*femb_i.*femb_xxj ...
                    + r21.*femb_xi.*femb_xxj ...
                    + r12.*femb_xxi.*femb_xj ...
                    + r02.*femb_xxi.*femb_j ...
                    )*w;
                if(ie==1)
                    diaginfer(2*ie-1)=diaginfer(2*ie-1)+(...
                      b21inv*femb_xi(1)*femb_xxj(1) ...
                    + b20inv*femb_i(1)*femb_xxj(1) ...
                    + b11inv*femb_xi(1)*femb_xj(1) ...
                    + b10inv*femb_i(1)*femb_xj(1) ...
                    + b01inv*femb_xi(1)*femb_j(1) ...
                    + b00inv*femb_i(1)*femb_j(1) );
                end
                if(ie==ne)
                    diaginfer(2*ie-1)=diaginfer(2*ie-1)+(...
                      b21outv*femb_xi(end)*femb_xxj(end) ...
                    + b20outv*femb_i(end)*femb_xxj(end) ...
                    + b11outv*femb_xi(end)*femb_xj(end) ...
                    + b10outv*femb_i(end)*femb_xj(end) ...
                    + b01outv*femb_xi(end)*femb_j(end) ...
                    + b00outv*femb_i(end)*femb_j(end) );
                end
            elseif(i==j-2)
                diagsuper(2*ie+1)=diagsuper(2*ie+1)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    + r22.*femb_xxi.*femb_xxj ...
                    + r20.*femb_i.*femb_xxj ...
                    + r21.*femb_xi.*femb_xxj ...
                    + r12.*femb_xxi.*femb_xj ...
                    + r02.*femb_xxi.*femb_j ...
                    )*w;
                if(ie==1)
                    diagsuper(2*ie+1)=diagsuper(2*ie+1)+(...
                      b21inv*femb_xi(1)*femb_xxj(1) ...
                    + b20inv*femb_i(1)*femb_xxj(1) ...
                    + b11inv*femb_xi(1)*femb_xj(1) ...
                    + b10inv*femb_i(1)*femb_xj(1) ...
                    + b01inv*femb_xi(1)*femb_j(1) ...
                    + b00inv*femb_i(1)*femb_j(1) );
                end
                if(ie==ne)
                    diagsuper(2*ie+1)=diagsuper(2*ie+1)+(...
                      b21outv*femb_xi(end)*femb_xxj(end) ...
                    + b20outv*femb_i(end)*femb_xxj(end) ...
                    + b11outv*femb_xi(end)*femb_xj(end) ...
                    + b10outv*femb_i(end)*femb_xj(end) ...
                    + b01outv*femb_xi(end)*femb_j(end) ...
                    + b00outv*femb_i(end)*femb_j(end) );
                end
            elseif(i==j-1 && i==2)
                diagsup(2*ie+1)=diagsup(2*ie+1)+( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    + r22.*femb_xxi.*femb_xxj ...
                    + r20.*femb_i.*femb_xxj ...
                    + r21.*femb_xi.*femb_xxj ...
                    + r12.*femb_xxi.*femb_xj ...
                    + r02.*femb_xxi.*femb_j ...
                    )*w;
                if(ie==1)
                    diagsup(2*ie+1)=diagsup(2*ie+1)+( ...
                      b21inv*femb_xi(1)*femb_xxj(1) ...
                    + b20inv*femb_i(1)*femb_xxj(1) ...
                    + b11inv*femb_xi(1)*femb_xj(1) ...
                    + b10inv*femb_i(1)*femb_xj(1) ...
                    + b01inv*femb_xi(1)*femb_j(1) ...
                    + b00inv*femb_i(1)*femb_j(1) );
                end
                if(ie==ne)
                    diagsup(2*ie+1)=diagsup(2*ie+1)+( ...
                      b21outv*femb_xi(end)*femb_xxj(end) ...
                    + b20outv*femb_i(end)*femb_xxj(end) ...
                    + b11outv*femb_xi(end)*femb_xj(end) ...
                    + b10outv*femb_i(end)*femb_xj(end) ...
                    + b01outv*femb_xi(end)*femb_j(end) ...
                    + b00outv*femb_i(end)*femb_j(end) );
                end
            elseif(i==j-1 && i==1)
                diagsup(2*ie)=diagsup(2*ie)+( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    + r22.*femb_xxi.*femb_xxj ...
                    + r20.*femb_i.*femb_xxj ...
                    + r21.*femb_xi.*femb_xxj ...
                    + r12.*femb_xxi.*femb_xj ...
                    + r02.*femb_xxi.*femb_j ...
                    )*w;
                if(ie==1)
                    diagsup(2*ie)=diagsup(2*ie)+( ...
                      b21inv*femb_xi(1)*femb_xxj(1) ...
                    + b20inv*femb_i(1)*femb_xxj(1) ...
                    + b11inv*femb_xi(1)*femb_xj(1) ...
                    + b10inv*femb_i(1)*femb_xj(1) ...
                    + b01inv*femb_xi(1)*femb_j(1) ...
                    + b00inv*femb_i(1)*femb_j(1) );
                end
                if(ie==ne)
                    diagsup(2*ie)=diagsup(2*ie)+( ...
                      b21outv*femb_xi(end)*femb_xxj(end) ...
                    + b20outv*femb_i(end)*femb_xxj(end) ...
                    + b11outv*femb_xi(end)*femb_xj(end) ...
                    + b10outv*femb_i(end)*femb_xj(end) ...
                    + b01outv*femb_xi(end)*femb_j(end) ...
                    + b00outv*femb_i(end)*femb_j(end) );
                end
            end
        end
        bl(2*(ie-1)+i) = bl(2*(ie-1)+i) ...
                       + ( force0v((ie-1)*nqnx+1:ie*nqnx)'.*femb_i )*w ...
                       + ( force1v((ie-1)*nqnx+1:ie*nqnx)'.*femb_xi )*w ...
                       + ( force2v((ie-1)*nqnx+1:ie*nqnx)'.*femb_xxi )*w;
        if(ie==1)
            bl(2*(ie-1)+i) = bl(2*(ie-1)+i) ...
                           + force1inv*femb_xi(1) ...
                           + force0inv*femb_i(1);
        end
        if(ie==ne)
            bl(2*(ie-1)+i) = bl(2*(ie-1)+i) ...
                           + force1outv*femb_xi(end) ...
                           + force0outv*femb_i(end);
        end
    end
        
end
Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A
end

function [Al,bl]=FSIAssembleVel1DYold(Map,wyq,mb_k,mb_r,mb_yk,mb_yr,ne,...
    mesh_wx,nqnx,FEM_basis,FEM_basis_x,dt,nu,Force,Bx,By,...
    E,MBrobKY,MBrobRY,etadown_old,etaup_old,...
    rho_s,h_s,uup_old,udown_old,b,c,AleMap)

% per ogni x (vettorialmente)
% calcola r^st_kr    
r11v = ( nu*Map.Jac)*( mb_k.*mb_r.*wyq)...
     - (b*dt+c)*AleMap.JacUp.*Map.Aup*MBrobKY(1)*MBrobRY(1) ...
     - (b*dt+c)*AleMap.JacDown.*Map.Adown*MBrobKY(2)*MBrobRY(2);
r10v = ( nu*Map.Jac.*Map.D)*( mb_k.*mb_yr .*wyq ) ...
     + (    Map.Jac.*Bx)*( mb_k  .*mb_r .*wyq );
r01v = ( nu*Map.Jac.*Map.D)*( mb_yk.*mb_r .*wyq );
r00v =(nu*Map.Jac.*(Map.D.^2+2*Map.J.^2) ) *(mb_yk.*mb_yr.*wyq ) ...
     + (    Map.Jac/dt ) *(mb_k.*mb_r.*wyq ) ...
     + (    Map.Jac.*(Map.D.*Bx+Map.J.*By))*(mb_yk .*mb_r .*wyq) ...
	 + (rho_s*h_s/dt+E*dt)*AleMap.JacUp.*Map.Aup*MBrobKY(1)*MBrobRY(1) ...
     + (rho_s*h_s/dt+E*dt)*AleMap.JacDown.*Map.Adown*MBrobKY(2)*MBrobRY(2);
forcev = ( Map.Jac.*Force )*( mb_r.*wyq)...
        - E*etaup_old.*AleMap.JacUp.*Map.Aup*MBrobRY(1) ...
        - E*etadown_old.*AleMap.JacDown.*Map.Adown*MBrobRY(2)...
        + rho_s*h_s/dt*uup_old.*AleMap.JacUp.*Map.Aup*MBrobRY(1) ...
        + rho_s*h_s/dt*udown_old.*AleMap.JacDown.*Map.Adown*MBrobRY(2)...
        - b*etaup_old.*AleMap.JacUp.*Map.Aup*MBrobRY(1) ...
        - b*etadown_old.*AleMap.JacDown.*Map.Adown*MBrobRY(2);
    
nx=ne+1;

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag =   zeros( ne + nx, 1);
diagsup  =   zeros( ne + nx, 1); % Spdiags taglierà il primo elemento
diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglierà l'ultimo
diagsuper  = zeros( ne + nx, 1); % Spdiags taglierà i primi due elementi
diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglierà gli ultimi due
bl         = zeros( ne + nx, 1);
for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r11 = r11v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r01 = r01v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:3                  % Per ogni funzione di forma EF test
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basis  ( i, :);
        femb_xi = FEM_basis_x( i, :);
        for j = 1:3             % Per ogni funzione di forma EF sol
            
            femb_j  = FEM_basis  ( j, :);
            femb_xj = FEM_basis_x( j, :);
            
            if(i==j)
                maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+...
                    (  r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;                
            elseif(i==j+1 && i==2)
                diaginf(2*ie-1)=diaginf(2*ie-1)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+1 && i==3)
                diaginf(2*ie)=diaginf(2*ie)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+2)
                diaginfer(2*ie-1)=diaginfer(2*ie-1)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-2)
                diagsuper(2*ie+1)=diagsuper(2*ie+1)+(...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-1 && i==2)
                diagsup(2*ie+1)=diagsup(2*ie+1)+( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j-1 && i==1)
                diagsup(2*ie)=diagsup(2*ie)+( ...
                    r11.*femb_xi.*femb_xj ...
                    + r01.*femb_xi.*femb_j...
                    + r10.*femb_i.*femb_xj ...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            end
        end
        bl(2*(ie-1)+i) = bl(2*(ie-1)+i) + ( forcev((ie-1)*nqnx+1:ie*nqnx)'.*femb_i )*w;
    end
end
Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A
end

%% Method 'AssembleVelYX'

function Al = AssembleVelYX(Map,wyq,MbKY,DMbKY,DMbRX,FEM_basis,FEM_basis_x,ne,...
    mesh_wx,nqnx,nu)

% per ogni x (vettorialmente)
% calcola r^st_kr
r10v = ( nu*Map.Jac.*Map.J)*( MbKY.*DMbRX .*wyq );
r00v = ( nu*Map.Jac.*Map.D.*Map.J ) *(DMbKY.*DMbRX.*wyq );

nx=ne+1;

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag =   zeros( ne + nx, 1);
diagsup  =   zeros( ne + nx, 1); % Spdiags taglier� il primo elemento
diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglier� l'ultimo
diagsuper  = zeros( ne + nx, 1); % Spdiags taglier� i primi due elementi
diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglier� gli ultimi due

for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:3                  % Per ogni funzione di forma EF test
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basis  ( i , : );
        for j = 1:3             % Per ogni funzione di forma EF sol
            
            femb_j  = FEM_basis  ( j, :);
            femb_xj  = FEM_basis_x  ( j, :);
            
            if(i==j)
                maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+...
                    (  r10.*femb_i.*femb_xj...
                    +r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j+1 && i==2)
                diaginf(2*ie-1)=(...
                    + r10.*femb_i.*femb_xj...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+1 && i==3)
                diaginf(2*ie)=(...
                    + r10.*femb_i.*femb_xj...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+2)
                diaginfer(2*ie-1)=(...
                    + r10.*femb_i.*femb_xj...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-2)
                diagsuper(2*ie+1)=(...
                    + r10.*femb_i.*femb_xj...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-1 && i==2)
                diagsup(2*ie+1)=( ...
                    + r10.*femb_i.*femb_xj...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j-1 && i==1)
                diagsup(2*ie)=( ...
                    + r10.*femb_i.*femb_xj...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            end
        end
    end
end
Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A

end

%% Method 'AssembleVelXY'

function Al = AssembleVelXY(Map,wyq,DMbKX,MbRY,DMbRY,FEM_basis,FEM_basis_x,ne,...
    mesh_wx,nqnx,nu)

% per ogni x (vettorialmente)
% calcola r^st_kr
r01v = ( nu*Map.Jac.*Map.J)*( DMbKX.*MbRY .*wyq );
r00v = ( nu*Map.Jac.*Map.D.*Map.J ) *(DMbKX.*DMbRY.*wyq );

nx=ne+1;

% Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
maindiag =   zeros( ne + nx, 1);
diagsup  =   zeros( ne + nx, 1); % Spdiags taglierà il primo elemento
diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglierà l'ultimo
diagsuper  = zeros( ne + nx, 1); % Spdiags taglierà i primi due elementi
diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglierà gli ultimi due

for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r01 = r01v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for i = 1:3                  % Per ogni funzione di forma EF test
        
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basis  (i,:);
        femb_xi = FEM_basis_x(i,:);
        for j = 1:3             % Per ogni funzione di forma EF sol
            
            femb_j  = FEM_basis  ( j,:);
            
            if(i==j)
                maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+...
                    (  r01.*femb_xi.*femb_j...
                    +r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j+1 && i==2)
                diaginf(2*ie-1)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+1 && i==3)
                diaginf(2*ie)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j+2)
                diaginfer(2*ie-1)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-2)
                diagsuper(2*ie+1)=(...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            elseif(i==j-1 && i==2)
                diagsup(2*ie+1)=( ...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
                
            elseif(i==j-1 && i==1)
                diagsup(2*ie)=( ...
                    + r01.*femb_xi.*femb_j...
                    + r00.*femb_i.*femb_j  ...
                    )*w;
            end
        end
    end
end
Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A
end

%% Method 'AssembleBX'

function Al=AssembleBX(Map,wyq,mb_k,mb_s,mb_yk,ne,mesh_wx,nqnx,FEM_basis,FEM_basis_x,FEM_basisP)

r00v = -( Map.Jac.*Map.D) *(mb_yk.*mb_s.*wyq );
r10v = -( Map.Jac )*(mb_k.*mb_s.*wyq );

nx=ne+1;
Al=zeros(nx,nx+ne);
for ie = 1 : ne  % Per ogni intervallo
    
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    
    for j = 1:3     % Per ogni funzione di forma EF vel
        femb_j  = FEM_basis  ( j, :);
        femb_xj = FEM_basis_x( j, :);
        for i = 1:2  % Per ogni funzione di forma EF press
            %selezione dei pezzi corretti per le basi FEM
            femb_i  = FEM_basisP  ( i, :);
            Al(ie-1+i,2*ie-2+j) =  Al(ie-1+i,2*ie-2+j)+(r00.*femb_i.*femb_j+r10.*femb_i.*femb_xj)*w;
        end
    end
end
end

%% Method 'AssembleBY'

function Al=AssembleBY(Map,wyq,mb_s,mb_yk,ne,mesh_wx,nqnx,FEM_basis,FEM_basisP)

r00v = -( Map.Jac.*Map.J) *(mb_yk.*mb_s.*wyq );

nx=ne+1;
Al=zeros(nx,nx+ne);
for ie = 1 : ne  % Per ogni intervallo
    %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
    w = mesh_wx( (ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
    % Seleziono i coefficiweenti relativi all'intervallo in questione
    r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
    for i = 1:2                  % Per ogni funzione di forma EF pressione
        %selezione dei pezzi corretti per le basi FEM
        femb_i  = FEM_basisP  ( i, :);
        for j = 1:3             % Per ogni funzione di forma EF velocità
            femb_j  = FEM_basis  ( j, :);
            Al(ie-1+i,2*ie-2+j) = Al(ie-1+i,2*ie-2+j) + (r00.*femb_i.*femb_j)*w;
        end
    end
end
end

%% Extract
% %% Method 'FSIAssembleVel1DYold'
% 
% function [Al,bl]=FSIAssembleVel1DYold(Map,wyq,mb_k,mb_r,mb_yk,mb_yr,ne,...
%     mesh_wx,nqnx,FEM_basis,FEM_basis_x,dt,nu,Force,Bx,By,...
%     E,MBrobKY,MBrobRY,etadown_old,etaup_old)
% 
% % per ogni x (vettorialmente)
% % calcola r^st_kr
% r11v = ( nu*Map.Jac)*( mb_k.*mb_r.*wyq);
% r10v = ( nu*Map.Jac.*Map.D)*( mb_k.*mb_yr .*wyq ) ...
%      + (    Map.Jac.*Bx)*( mb_k  .*mb_r .*wyq );
% r01v = ( nu*Map.Jac.*Map.D)*( mb_yk.*mb_r .*wyq );
% r00v =(nu*Map.Jac.*(Map.D.^2+2*Map.J.^2) ) *(mb_yk.*mb_yr.*wyq ) ...
%      + (    Map.Jac/dt ) *(mb_k.*mb_r.*wyq ) ...
%      + (    Map.Jac.*(Map.D.*Bx+Map.J.*By))*(mb_yk .*mb_r .*wyq) ...
%      + E*dt*Map.Aup*MBrobKY(1)*MBrobRY(1) ...
%      + E*dt*Map.Adown*MBrobKY(2)*MBrobRY(2);
% forcev = ( Map.Jac.*Force )*( mb_r.*wyq)...
%         - E*etaup_old.*Map.Aup*MBrobRY(1) ...
%         - E*etadown_old.*Map.Adown*MBrobRY(2);
%     
% nx=ne+1;
% 
% % Inizializzo le diagonali del blocco, si aumenta l'efficienza nell'assemblaggio di matrici sparse
% maindiag =   zeros( ne + nx, 1);
% diagsup  =   zeros( ne + nx, 1); % Spdiags taglierà il primo elemento
% diaginf  =   zeros( ne + nx, 1); % Spdiags invece taglierà l'ultimo
% diagsuper  = zeros( ne + nx, 1); % Spdiags taglierà i primi due elementi
% diaginfer  = zeros( ne + nx, 1); % Spdiags invece taglierà gli ultimi due
% bl         = zeros( ne + nx, 1);
% for ie = 1 : ne  % Per ogni intervallo
%     
%     %(ie-1)*nqnx+1:ie*nqnx seleziona l'intervallo di punti corretto
%     w = mesh_wx((ie-1)*nqnx+1:ie*nqnx); % Calcolo i pesi dell'intervallo in questione
%     
%     % Seleziono i coefficiweenti relativi all'intervallo in questione
%     r11 = r11v ( (ie-1)*nqnx+1 : ie*nqnx)';
%     r10 = r10v ( (ie-1)*nqnx+1 : ie*nqnx)';
%     r01 = r01v ( (ie-1)*nqnx+1 : ie*nqnx)';
%     r00 = r00v ( (ie-1)*nqnx+1 : ie*nqnx)';
%     
%     for i = 1:3                  % Per ogni funzione di forma EF test
%         
%         %selezione dei pezzi corretti per le basi FEM
%         femb_i  = FEM_basis  ( i, :);
%         femb_xi = FEM_basis_x( i, :);
%         for j = 1:3             % Per ogni funzione di forma EF sol
%             
%             femb_j  = FEM_basis  ( j, :);
%             femb_xj = FEM_basis_x( j, :);
%             
%             if(i==j)
%                 maindiag(2*ie+i-2)=maindiag(2*ie+i-2)+...
%                     (  r11.*femb_xi.*femb_xj ...
%                     + r01.*femb_xi.*femb_j...
%                     + r10.*femb_i.*femb_xj ...
%                     + r00.*femb_i.*femb_j  ...
%                     )*w;                
%             elseif(i==j+1 && i==2)
%                 diaginf(2*ie-1)=diaginf(2*ie-1)+(...
%                     r11.*femb_xi.*femb_xj ...
%                     + r01.*femb_xi.*femb_j...
%                     + r10.*femb_i.*femb_xj ...
%                     + r00.*femb_i.*femb_j  ...
%                     )*w;
%             elseif(i==j+1 && i==3)
%                 diaginf(2*ie)=diaginf(2*ie)+(...
%                     r11.*femb_xi.*femb_xj ...
%                     + r01.*femb_xi.*femb_j...
%                     + r10.*femb_i.*femb_xj ...
%                     + r00.*femb_i.*femb_j  ...
%                     )*w;
%             elseif(i==j+2)
%                 diaginfer(2*ie-1)=diaginfer(2*ie-1)+(...
%                     r11.*femb_xi.*femb_xj ...
%                     + r01.*femb_xi.*femb_j...
%                     + r10.*femb_i.*femb_xj ...
%                     + r00.*femb_i.*femb_j  ...
%                     )*w;
%             elseif(i==j-2)
%                 diagsuper(2*ie+1)=diagsuper(2*ie+1)+(...
%                     r11.*femb_xi.*femb_xj ...
%                     + r01.*femb_xi.*femb_j...
%                     + r10.*femb_i.*femb_xj ...
%                     + r00.*femb_i.*femb_j  ...
%                     )*w;
%             elseif(i==j-1 && i==2)
%                 diagsup(2*ie+1)=diagsup(2*ie+1)+( ...
%                     r11.*femb_xi.*femb_xj ...
%                     + r01.*femb_xi.*femb_j...
%                     + r10.*femb_i.*femb_xj ...
%                     + r00.*femb_i.*femb_j  ...
%                     )*w;
%                 
%             elseif(i==j-1 && i==1)
%                 diagsup(2*ie)=diagsup(2*ie)+( ...
%                     r11.*femb_xi.*femb_xj ...
%                     + r01.*femb_xi.*femb_j...
%                     + r10.*femb_i.*femb_xj ...
%                     + r00.*femb_i.*femb_j  ...
%                     )*w;
%             end
%         end
%         bl(2*(ie-1)+i) = bl(2*(ie-1)+i) + ( forcev((ie-1)*nqnx+1:ie*nqnx)'.*femb_i )*w;
%     end
% end
% Al=spdiags( [diaginfer,diaginf,maindiag,diagsup,diagsuper], -2:2, (ne+nx), (nx+ne) ); % Assegno le diagonali alla matrice A
% end