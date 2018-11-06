classdef AssemblerStokesHandler
    
    %% ASSEMBLER STOKES HANDLER CLASS
    % the AssemblerStokesHandler is a class that contain all the scripts
    % responsible for the assembling and building of the block matrix that
    % describe the discretized differential problem. On one hand, the class
    % properties, the properties of the objects were defined based on the
    % variables needed in the implementation of each function. On the other 
    % hand, all of the previous functions are concentrated and organized in
    % the class methods.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% ASSEMBLER STOKES HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the AssemblerStokeshandler
        % encapsulate all of the variables needed to run the methods
        % bellow, from the support methods to the final building methods.
        
        % BUILDING PROPERTIES
                
        dimModalBasisU;              % Dimension of the Modal Basis for the velocity         
        dimModalBasisP;              % Dimension of the Modal Basis for the pressure
        
        leftBDomain_inX;            % Left Limit of the Domain in the X Direction
   
        rightBDomain_inX;           % Right Limit of the Domain in the X Direction
        
        stepMeshX;                  % Vector Containing the Step of the Finite
                                    % Element Mesh
                      
        label_upBoundDomain;        % Contains the Label Identifying the Nature of
                                    % the Boundary Conditions on the Upper Limit of
                                    % the Domain                 
        label_downBoundDomain;     % Contains the Label Identifying the Nature of
                                    % the Boundary Conditions on the Lower Limit of
                                    % the Domain
                      
        localdata_upBDomainX;       % Contains the Values of the Boundary Conditions
        localdata_upBDomainY;       % on the Upper Limir of the Domain
                      
        localdata_downBDomainX;      % Contains the Values of the Boundary Conditions
        localdata_downBDomainY;      % on the Lower Limir of the Domain
                      
        
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
    
        robinCondStruct;            % NO Data Structure Containing the Two Values of the
                                    % Coefficients (R, L) for the Robin Condition Used
                                    % in the Domain Decomposition
                      
        physicMesh_inX;             % Vector Containing the Physical Mesh in the X
                                    % Direction
                      
        physicMesh_inY;             % Vector Containing the Physical Mesh in the Y
                                    % Direction
                      
        jacAtQuadNodes;             % Data Sturcture Used to Save the Value of the
                                    % Jacobians Computed in the Quadratures Nodes 
                      
        degreePolySplineBasisU;      % Degree of the Polynomial B-Spline Basis
        degreePolySplineBasisP;
        
        continuityParameterU;        % Degree of Continuity of the Basis 'C^(p-k)'
        continuityParameterP;
        
        domainProfile;              % Symbolic Function Defining the Profile of the
                                    % Simulation Domain
                      
        domainProfileDer;           % Symbolic Function Defining the Derivative of
                                    % the Profile of the Simulation Domain   
        
        numbHorQuadNodes;           % Number of horizontal nodes to apply the quadrature
                                    % formula
                                    
        numbVerQuadNodes;           % Number of vertical nodes to apply the quadrature
                                    % formula
        
    end
    
    methods (Access = public)
        
        %% ASSEMBLER STOKES HANDLER - CONSTRUCT METHOD
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
        
        %% ASSEMBLER STOKES HANDLER - BUILDING METHODS
        
        
            %% Method 'buildSystemIGA'
            
            function [A,b,modalBasisU,modalBasisP,liftCoeffAX,liftCoeffBX,liftCoeffAY,liftCoeffBY,...
                      verGLNodes,verGLWeights] = buildSystemIGA(obj)
                                                         
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
            %   (1)  dimModalBasisU          : Dimension of the Modal basis for velocity
            %   (2)  dimModalBasisP            : Dimension of the Modal basis for pressure             
            %   (3)  leftBDomain_inX           : Left Limit of the Domain in the X Direction
            %   (4)  rightBDomain_inX          : Right Limit of the Domain in the X Direction
            %   (5)  stepMeshX                 : Vector Containing the Step of the Finite
            %                                    Element Mesh
            %   (6)  label_upBoundDomainX      : Contains the Label Identifying the Nature of
            %                                    the Boundary Conditions on the Upper Limit of
            %                                    the Domain for u_x
            %   (7)  label_upBoundDomainY      : Contains the Label Identifying the Nature of
            %                                    the Boundary Conditions on the Upper Limit of
            %                                    the Domain for u_y
            %   (8)  label_downBoundDomainX    : Contains the Label Identifying the Nature of
            %                                    the Boundary Conditions on the Lower Limit of
            %                                    the Domain for u_x
            %   (9)  label_downBoundDomainY    : Contains the Label Identifying the Nature of
            %                                    the Boundary Conditions on the Lower Limit of
            %                                    the Domain for u_y
            %   (10)  localdata_upBDomainX     : Contains the Values of the Boundary Conditions
            %                                    on the Upper Limir of the
            %                                    Domain for u_x
            %   (11)  localdata_upBDomainY     : Contains the Values of the Boundary Conditions
            %                                    on the Upper Limir of the
            %                                    Domain for u_y
            %   (12) localdata_downBDomainX    : Contains the Values of the Boundary Conditions
            %                                    on the Lower Limir of the
            %                                    Domain for u_x
            %   (13) localdata_downBDomainY    : Contains the Values of the Boundary Conditions
            %                                    on the Lower Limir of the
            %                                    Domain for u_y
            %   (14) coefficientForm            : Data Strusture Containing All the @-Functions
            %                                    and the Constants Relative to the Bilinear Form
            %   (15) dirCondFuncStruct         : Data Structure Containing All the @-Functions
            %                                    for the Dirichlet Conditions at the Inflow and
            %                                    for the Exciting Forces
            %   (16) geometricInfo             : Data Structure Containing All the
            %                                    Geometric Information regarding the
            %                                    Domain. The current version of the code
            %                                    works only for the specific condition of:
            %                                    (L = 1, a = 0, psi_x = 0)
            %   (17) physicMesh_inX            : Vector Containing the Physical Mesh in the X
            %                                    Direction
            %   (18) physicMesh_inY            : Vector Containing the Physical Mesh in the Y
            %                                    Direction
            %   (19) jacAtQuadNodes            : Data Sturcture Used to Save the Value of the
            %                                     Jacobians Computed in the Quadratures Nodes 
            %   (20) degreePolySplineBasisU    : Degree of the Polynomial
            %                                    B-Spline Basis for the velocity
            %   (21) degreePolySplineBasisP    : Degree of the Polynomial B-Spline Basis
            %   (22) continuityParameterU      : Degree of Continuity of
            %                                    the Basis 'C^(p-k)' for the velocity
            %   (23) continuityParameterP      : Degree of Continuity of
            %                                     the Basis 'C^(p-k)' for the pressure
            %   (24) domainProfile             : Symbolic Function Defining the Profile of the
            %                                    Simulation Domain
            %   (25) domainProfileDer          : Symbolic Function Defining the Derivative of
            %                                    the Profile of the Simulation Domain
            % 
            % The outputs are:
            %%
            %   (1) A                   : Final Assembled Block Matrix Using IGA Basis
            %   (2) b                   : Final Assembled Block Vector Using IGA Basis
            %   (3) modalBasis          : Modal Basis Obtained
            %   (3) aLiftX              : First Offset Adjustment Coefficient for u_x
            %   (4) bLiftX              : Second Offset Adjustment Coefficient for u_x
            %   (5) aLiftY              : First Offset Adjustment Coefficient for u_y
            %   (6) bLiftY              : Second Offset Adjustment Coefficient for u_y
            %   (5) intNodesGaussLeg    : Gauss-Legendre Integration Nodes
            %   (6) intNodesWeights     : Weights of the Gauss-Legendre Integration Nodes
            
            %% IMPORT CLASSES
            
            import Core.AssemblerStokesHandler
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
            
            numbControlPtsU = numbKnots*obj.continuityParameterU + obj.degreePolySplineBasisU +...
                              1 - obj.continuityParameterU;
            numbControlPtsP = numbKnots*obj.continuityParameterP + obj.degreePolySplineBasisP +...
                              1 - obj.continuityParameterP;
            
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
            
            % Isogeometric mesh for u
            stepKnotU = 1/(numbControlPtsU-1);
            meshIGAU = zeros(1,numbControlPtsU);
            
            for i = 2:numbControlPtsU
                 meshIGAU(i) = meshIGAU(i-1) + stepKnotU;
            end
            
            % Isogeometric mesh for p
            stepKnotP = 1/(numbControlPtsP-1);
            meshIGAP = zeros(1,numbControlPtsP);
            
            for i = 2:numbControlPtsP
                 meshIGAP(i) = meshIGAP(i-1) + stepKnotP;
            end
            
                        %% AUGMENTED ISOGEOMETRIC MESH + WEIGHTS
            %-------------------------------------------------------------%
            % Creation of the finite element mesh in the X direction
            % considering equispaced nodes. The mesh is created using the
            % total number of nodes and the left limit of the domain. Both
            % information come from the demo file and are passed here as a
            % property of the object.
            %-------------------------------------------------------------%

            augMeshIGAU = zeros( (numbControlPtsU-1)*numbHorNodes, 1);
            augMeshIGAP = zeros( (numbControlPtsP-1)*numbHorNodes, 1);

            %-------------------------------------------------------------%
            % Note: The loop allocates the correponding quadrature nodes
            % and weights corresponding to each knot of the physical mesh
            % in the isogeometric analysis.
            %-------------------------------------------------------------%
            
            % Velocity
            for i = 1:numbKnots 
                
                % STEP 1
                %---------------------------------------------------------%
                % In the first step, the Gauss-Legendre nodes computed
                % previously are rescaled to fit the interval corresponding
                % to the current element.
                %---------------------------------------------------------%
                
                obj_quadratureRule = IntegrateHandler();
                
                obj_quadratureRule.leftBoundInterval = meshIGAU(i);
                obj_quadratureRule.rightBoundInterval = meshIGAU(i+1);
                obj_quadratureRule.inputNodes = horGLNodes;
                obj_quadratureRule.inputWeights = horGLWeights;
                
                [augMeshIGAU((i-1)*numbHorNodes+1 : i*numbHorNodes), ~] = ...
                                                 quadratureRule(obj_quadratureRule);
             
                % STEP 2
                %---------------------------------------------------------%
                % In the second step, the nodes and weights are again
                % rescaled to take in consideration the geometry of the
                % domain, more specifically the profile of the centerline.
                %---------------------------------------------------------%

                auxMesh = augMeshIGAU((i-1)*numbHorNodes+1 : i*numbHorNodes);

                for hp = 1: numbHorNodes
                    
                    obj_gaussLegendre = IntegrateHandler();
                    obj_gaussLegendre.numbQuadNodes = 16;
                    
                    [~,auxPoints,auxWeights] = gaussLegendre(obj_gaussLegendre);
                    
                    auxPoints = auxMesh(hp) * auxPoints;
                    auxWeights = auxMesh(hp) * auxWeights;
                    auxMesh(hp) = sum(sqrt(1+(obj.domainProfileDer(auxPoints)).^2).*auxWeights); 
                                       
                end
                augMeshIGAU((i-1)*numbHorNodes+1 : i*numbHorNodes) = auxMesh;
                                                    
            end
            
            % Pressure
            for i = 1:numbKnots
                
                % STEP 1
                %---------------------------------------------------------%
                % In the first step, the Gauss-Legendre nodes computed
                % previously are rescaled to fit the interval corresponding
                % to the current element.
                %---------------------------------------------------------%
                
                obj_quadratureRule = IntegrateHandler();
                
                obj_quadratureRule.leftBoundInterval = meshIGAP(i);
                obj_quadratureRule.rightBoundInterval = meshIGAP(i+1);
                obj_quadratureRule.inputNodes = horGLNodes;
                obj_quadratureRule.inputWeights = horGLWeights;
                
                [augMeshIGAP((i-1)*numbHorNodes+1 : i*numbHorNodes), ~] = ...
                                                 quadratureRule(obj_quadratureRule);
             
                % STEP 2
                %---------------------------------------------------------%
                % In the second step, the nodes and weights are again
                % rescaled to take in consideration the geometry of the
                % domain, more specifically the profile of the centerline.
                %---------------------------------------------------------%

                auxMesh = augMeshIGAP((i-1)*numbHorNodes+1 : i*numbHorNodes);

                for hp = 1: numbHorNodes
                    
                    obj_gaussLegendre = IntegrateHandler();
                    obj_gaussLegendre.numbQuadNodes = 16;
                    
                    [~,auxPoints,auxWeights] = gaussLegendre(obj_gaussLegendre);
                    
                    auxPoints = auxMesh(hp) * auxPoints;
                    auxWeights = auxMesh(hp) * auxWeights;
                    auxMesh(hp) = sum(sqrt(1+(obj.domainProfileDer(auxPoints)).^2).*auxWeights); 
                                       
                end
                augMeshIGAP((i-1)*numbHorNodes+1 : i*numbHorNodes) = auxMesh;
                                                    
            end
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            % Velocity
            obj_newModalBasis = BasisHandler();
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasisU;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomain;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomain;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasisU, modalBasisDerU] = newModalBasis(obj_newModalBasis);
                      
            
            % Pressure: I choose Neumann boundary conditions since no
            % Dirichlet conditions are assigned to the pressure.
            obj_newModalBasis2 = BasisHandler();
            
            obj_newModalBasis2.dimModalBasis = obj.dimModalBasisP;
            obj_newModalBasis2.evalNodesY = verGLNodes;
            obj_newModalBasis2.labelUpBoundCond = 'neu';
            obj_newModalBasis2.labelDownBoundCond = 'neu';
            obj_newModalBasis2.coeffForm = obj.coefficientForm;

            [modalBasisP, modalBasisDerP] = newModalBasis(obj_newModalBasis2);
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES

            A = sparse( 2*numbControlPtsU*obj.dimModalBasisU + numbControlPtsP*obj.dimModalBasisP, ...
                        2*numbControlPtsU*obj.dimModalBasisU + numbControlPtsP*obj.dimModalBasisP );
            b = zeros ( 2*numbControlPtsU*obj.dimModalBasisU + numbControlPtsP*obj.dimModalBasisP, 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%

            [horEvalNodesU,verEvalNodesU]   = meshgrid( augMeshIGAU, verGLNodes-0.5 );%!!!!!!
            [horEvalNodesP,~]   = meshgrid( augMeshIGAP, verGLNodes-0.5 );%!!!!!!
            
            % THICKNESS
            
            evalLU     = obj.geometricInfo.L(horEvalNodesU)';
            evalLP     = obj.geometricInfo.L(horEvalNodesP)';
               
            
            %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal and vertical meshes to evaluate the
            % coefficients of the bilinear form of the original equation in
            % the entire domain.
            %-------------------------------------------------------------%

            % DIFFUSION
            evalMu    = obj.coefficientForm.mu(horEvalNodesU,verEvalNodesU)';
            
            
            %% EVALUATION OF THE EXCITING FORCE 
            %-------------------------------------------------------------%
            % Finally, we use the horizontal and vertical meshes to
            % evaluate the exciting force acting on the system in the whole
            % domain.
            %-------------------------------------------------------------%

            evalForceX = obj.dirCondFuncStruct.forceX(horEvalNodesU,verEvalNodesU)';
            evalForceY = obj.dirCondFuncStruct.forceY(horEvalNodesU,verEvalNodesU)';
            
            %-------------------------------------------------------------%
            % Note: All of the coefficients of the bilinear form and the
            % exciting term evaluated in the domain, as well as the mesh of
            % vertical coordinates are stored in a data structure called
            % "Computed" to be treated in the assembler function.
            %-------------------------------------------------------------%

            ComputedX = struct('mu_c',evalMu,'force_c',evalForceX,'y',verEvalNodesU);
            ComputedY = struct('mu_c',evalMu,'force_c',evalForceY,'y',verEvalNodesU);

                          
            %% LIFTING
            %-------------------------------------------------------------%
            % Compute the lifting contribution in the force vector due to
            % the non homogeneous Dirichlet boundary conditions in the
            % lower and upper boundary.
            % Actually, it is not used at the moment, since homogeneous
            % Dirichlet conditions are considered.
            %-------------------------------------------------------------%
            
            obj_liftBoundCond = BoundaryConditionHandler();
    
            obj_liftBoundCond.labelUpBoundCond = obj.label_upBoundDomain;
            obj_liftBoundCond.labelDownBoundCond = obj.label_downBoundDomain;
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomainX;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomainX;
            obj_liftBoundCond.coeffForm = obj.coefficientForm;

            [aLiftX,bLiftX] = liftBoundCond(obj_liftBoundCond);
            
            obj_liftBoundCond.dataUpBoundCond = obj.localdata_upBDomainY;
            obj_liftBoundCond.dataDownBoundCond = obj.localdata_downBDomainY;
            
            [aLiftY,bLiftY] = liftBoundCond(obj_liftBoundCond);
            liftFuncX = @(x,y) aLiftX * y + bLiftX;
            liftFuncY = @(x,y) aLiftY*y + bLiftY;
            liftingX = liftFuncX(0,ComputedX.y)';
            liftingY = liftFuncY(0,ComputedY.y)';
            
            
            
            %% JACOBIAN
            %-------------------------------------------------------------%
            % This section computes the curved domain and extract the
            % Jacobian vector that maps the centerline in the physical
            % domain to the computational domain.
            %-------------------------------------------------------------%
            
            % Mesh IGA of the curved domain
            
            meshIGACurved    = zeros(numbKnots + 1,1);
            [verGLNodes,verGLWeights] = gauss(numbVerNodes);
            
            verGLWeights = (obj.rightBDomain_inX - obj.leftBDomain_inX)*(obj.stepMeshX)*(verGLWeights*0.5);
                        
            for j=1:numbKnots
        
                scalingVec = zeros(1,numbVerNodes);
                auxNodes = (obj.rightBDomain_inX - obj.leftBDomain_inX)*obj.stepMeshX*(verGLNodes * 0.5 + 0.5) ...
                           + obj.leftBDomain_inX + (j-1)*obj.stepMeshX*(obj.rightBDomain_inX - obj.leftBDomain_inX);

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
            
            % Assembling loop for A_x
            for imb = 1:obj.dimModalBasisU
                for kmb = 1:obj.dimModalBasisU

                    [Amb,bmb,liftCoeffAX,liftCoeffBX] = assemblerIGA( imb, kmb, numbControlPtsU, ...
                                            verGLWeights,modalBasisU(:,imb), ...
                                            modalBasisDerU(:,imb),modalBasisU(:,kmb),...
                                            modalBasisDerU(:,kmb), evalLU, ...
                                            ComputedX,obj.stepMeshX,...
                                            obj.degreePolySplineBasisU,obj.continuityParameterU,...
                                            obj.leftBDomain_inX,...
                                            obj.rightBDomain_inX,liftingX,aLiftX,bLiftX,Jac,...
                                            numbHorNodes,numbVerNodes,...
                                            numbKnots,horGLWeights);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPtsU : imb*numbControlPtsU , ...
                      1+(kmb-1)*numbControlPtsU : kmb*numbControlPtsU) = Amb;      
                end

                % Assignment of the Block Vector Just Assembled
                b( 1+(imb-1)*numbControlPtsU : imb*numbControlPtsU ) = bmb;               
            end
            disp('FINISHED ASSEMBLING A_X');
            
            %assembling loop for A_y
            for imb = 1:obj.dimModalBasisU
                for kmb = 1:obj.dimModalBasisU

                    [Amb,bmb,liftCoeffAY,liftCoeffBY] = assemblerIGA( imb, kmb, numbControlPtsU, ...
                                            verGLWeights,modalBasisU(:,imb), ...
                                            modalBasisDerU(:,imb),modalBasisU(:,kmb),...
                                            modalBasisDerU(:,kmb), evalLU, ...
                                            ComputedY,obj.stepMeshX,...
                                            obj.degreePolySplineBasisU,obj.continuityParameterU,...
                                            obj.leftBDomain_inX,...
                                            obj.rightBDomain_inX,liftingY,aLiftY,bLiftY,Jac,...
                                            numbHorNodes,numbVerNodes,...
                                            numbKnots,horGLWeights);

                    % Assignment of the Block Matrix Just Assembled           
                    A(1+(obj.dimModalBasisU+imb-1)*numbControlPtsU : ...
                      (imb+ obj.dimModalBasisU)*numbControlPtsU ...
                      ,1+(obj.dimModalBasisU+kmb-1)*numbControlPtsU :...
                      (kmb+obj.dimModalBasisU)*numbControlPtsU) = Amb;
                end

                % Assignment of the Block Vector Just Assembled                
                b( 1+(obj.dimModalBasisU+imb-1)*numbControlPtsU :...
                    (imb+obj.dimModalBasisU)*numbControlPtsU ) = bmb;

            end
            disp('FINISHED ASSEMBLING A_Y');         
            
            % Assembling loop for P_x
            for imb = 1:obj.dimModalBasisU
                for kmb = 1:obj.dimModalBasisP

                    Amb   = assemblerIGA_PX( imb, kmb,...
                                            numbControlPtsU, numbControlPtsP, ...
                                            verGLWeights,modalBasisU(:,imb), modalBasisDerU(:,imb), ...
                                            modalBasisP(:,kmb), ...
                                            evalLU, obj.stepMeshX,...
                                            obj.degreePolySplineBasisU,obj.degreePolySplineBasisP,...
                                            obj.continuityParameterU,obj.continuityParameterP, ...
                                            obj.leftBDomain_inX, obj.rightBDomain_inX,...
                                            Jac,numbHorNodes,numbVerNodes, ...
                                            numbKnots,horGLWeights);

                    % Assignment of the Block Matrix Just Assembled

                    A(1+(imb-1)*numbControlPtsU : imb*numbControlPtsU , ...
                      2*numbControlPtsU*obj.dimModalBasisU+1+(kmb-1)*numbControlPtsP : ...
                      2*numbControlPtsU*obj.dimModalBasisU+kmb*numbControlPtsP) = -Amb;

                end

            end
            disp('FINISHED ASSEMBLING P_X');    
            
            % Assembling loop for P_y
            for imb = 1:obj.dimModalBasisU
                for kmb = 1:obj.dimModalBasisP

                    Amb   = assemblerIGA_PY( imb, kmb, ...
                                            numbControlPtsU, numbControlPtsP, ...
                                            verGLWeights,modalBasisU(:,imb), modalBasisDerU(:,imb), ...
                                            modalBasisP(:,kmb),  ...
                                            evalLU, obj.stepMeshX,...
                                            obj.degreePolySplineBasisU,obj.degreePolySplineBasisP,...
                                            obj.continuityParameterU,obj.continuityParameterP, ...
                                            obj.leftBDomain_inX, obj.rightBDomain_inX,...
                                            Jac,numbHorNodes,numbVerNodes, ...
                                            numbKnots,horGLWeights);

                    % Assignment of the Block Matrix Just Assembled

                    A(numbControlPtsU*obj.dimModalBasisU+1+(imb-1)*numbControlPtsU : ...
                      numbControlPtsU*obj.dimModalBasisU+imb*numbControlPtsU , ...
                      2*numbControlPtsU*obj.dimModalBasisU+1+(kmb-1)*numbControlPtsP : ...
                      2*numbControlPtsU*obj.dimModalBasisU+kmb*numbControlPtsP) = -Amb;

                end

            end
            disp('FINISHED ASSEMBLING P_Y'); 
            
            % Assembling loop for Q_x
            for imb = 1:obj.dimModalBasisP
                for kmb = 1:obj.dimModalBasisU

                    [Amb,bmb]  = assemblerIGA_QX( imb, kmb, ...
                                            numbControlPtsP, numbControlPtsU, ...
                                            verGLWeights,modalBasisP(:,imb), ...
                                            modalBasisU(:,kmb), modalBasisDerU(:,kmb), ...
                                            evalLP, obj.stepMeshX,...
                                            obj.degreePolySplineBasisP,obj.continuityParameterP,...
                                            obj.degreePolySplineBasisU,obj.continuityParameterU, ...
                                            obj.leftBDomain_inX, obj.rightBDomain_inX,...
                                            liftingX,aLiftX,bLiftX,Jac,numbHorNodes,numbVerNodes, ...
                                            numbKnots,horGLWeights);

                    % Assignment of the Block Matrix Just Assembled

                    A(2*numbControlPtsU*obj.dimModalBasisU+1+(imb-1)*numbControlPtsP :...
                      2*numbControlPtsU*obj.dimModalBasisU+imb*numbControlPtsP , ...
                      1+(kmb-1)*numbControlPtsU : ...
                      kmb*numbControlPtsU) = -Amb;


                end
                
                b(2*numbControlPtsU*obj.dimModalBasisU+1+(imb-1)*numbControlPtsP : ...
                2*numbControlPtsU*obj.dimModalBasisU+imb*numbControlPtsP) = bmb;

            end
            disp('FINISHED ASSEMBLING Q_X');
            
            % Assembling loop for Q_y
            for imb = 1:obj.dimModalBasisP
                for kmb = 1:obj.dimModalBasisU

                    [Amb,bmb] = assemblerIGA_QY( imb, kmb,...
                                            numbControlPtsP, numbControlPtsU, ...
                                            verGLWeights,modalBasisP(:,imb),  ...
                                            modalBasisU(:,kmb), modalBasisDerU(:,kmb), ...
                                            evalLP, obj.stepMeshX,...
                                            obj.degreePolySplineBasisP,obj.continuityParameterP,...
                                            obj.degreePolySplineBasisU,obj.continuityParameterU, ...
                                            obj.leftBDomain_inX, obj.rightBDomain_inX,...
                                            liftingY,aLiftY,bLiftY,Jac,numbHorNodes,numbVerNodes, ...
                                            numbKnots,horGLWeights);

                    % Assignment of the Block Matrix Just Assembled

                    A(2*numbControlPtsU*obj.dimModalBasisU+1+(imb-1)*numbControlPtsP :...
                      2*numbControlPtsU*obj.dimModalBasisU+imb*numbControlPtsP , ...
                      numbControlPtsU*obj.dimModalBasisU+1+(kmb-1)*numbControlPtsU : ...
                      numbControlPtsU*obj.dimModalBasisU+kmb*numbControlPtsU) = -Amb;

                end
                b(2*numbControlPtsU*obj.dimModalBasisU+1+(imb-1)*numbControlPtsP : ...
                2*numbControlPtsU*obj.dimModalBasisU+imb*numbControlPtsP) = bmb;
            end
            disp('FINISHED ASSEMBLING Q_Y');

            end
    end
end

%% ASSEMBLER STOKES HANDLER - ASSEMBLING METHODS

   
%% Method 'assemblerIGA'
            
function [Al,bl,aLift,bLift] = assemblerIGA(imb,kmb,numbControlPts, ...
                                verGLWeights,mb_i,mb_yi,mb_k,mb_yk, ...
                                L, Computed, ...
                                horStep,p,k, ...
                                domainLeftLimit,domainRightLimit, ...
                                ~,aLift,bLift,jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights)

    %%
    % assemblerIGA   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
    % The inputs are:
    %%
    %   (1)  imb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (2)  kmb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (3)  numbControlPts    : Size of the Functional Basis. Equal to the 
    %                            Dimention of the Spline Basis (Number of
    %                            Control Points)
    %   (4)  verGLWeights      : Vector Containing the Weights of the Quadrature 
    %                            Nodes in the Whole Mesh 
    %   (5)  mb_i              : Modal Basis referring for coeff. (imb,kmb)
    %   (6)  mb_yi             : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (7)  mb_k              : Modal Basis referring to coeff. (imb,kmb)
    %   (8)  mb_yk             : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %   (10) L                 : Thickness of the Channel Used in the Example
    %   (11) Computed          : Structure containing the coefficiets mu  and the force exciting the system 
    %   (12) horStep           : Step of the Finite Element Mesh
    %   (13) p                 : Degree of the Polynomial B-Spline Basis
    %   (14) k                 : Degree of Continuity of the Basis 'C^(p-k)'
    %   (15) domainLeftLimit,
    %        domainRightLimit  : Extremes of the domain
    %   (16) aLift             : Coefficient of lifting
    %   (17) bLift             : Coefficient of lifting
    %   (18) jacIGA            : Jacobian computed at the quadrature nodes
    %   (19) numbHordQuadNodes : Number of quadrature nodes in the
    %                            horizontal direction
    %   (20) nKnots            : Number of knots for the isogeometric mesh
    %   (21) horGLWeights      : Horizontal Gauss Legendre integration
    %                            nodes
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final Assembled Block Vector Using IGA Basis
    %   (3) aLift : Primo Coefficiente di Rilevamento
    %   (4) bLift : Secondo Coefficiente di Rilevamento
    

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
    
    objForceInt.funcToIntegrate =   L.*Computed.force_c;
                                
    objForceInt.funcWeight = mb_i.*verGLWeights;
    
    forceVec  = integrate(objForceInt);

    
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

    
    obj_integrate_1.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_2.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_3.funcToIntegrate = L.*Computed.mu_c;
    obj_integrate_4.funcToIntegrate = L.*Computed.mu_c;

    
    obj_integrate_1.funcWeight = mb_k  .*mb_i  .*verGLWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_yi .*verGLWeights;
    obj_integrate_3.funcWeight = mb_yk .*mb_i  .*verGLWeights;
    obj_integrate_4.funcWeight = mb_yk .*mb_yi .*verGLWeights;

    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    lambda_3   = integrate(obj_integrate_3);
    lambda_4   = integrate(obj_integrate_4);

    
    r00 = lambda_4 * (D1^2 + D2^2);
    r10 = lambda_2 * D1;
    r01 = lambda_3 * D1;
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
            
            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshape*cpi';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshape = dshape/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL RHS VECTOR

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
        Al(1,1)  = 0;
        Al(1,2)  = 0;
        Al(1,3)  = 0;
    end

    %disp('Finished ASSEMBLING LOOP');
end


%% Method 'assemblerIGA_PX'
            
function Al = assemblerIGA_PX(~,~,numbControlPtsU,numbControlPtsP, ...
                                verGLWeights,mb_i,mb_yi,mb_k, ...
                                LU, ...
                                horStep,pU,pP,kU,kP, ...
                                domainLeftLimit,domainRightLimit, ...
                                jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights)

    %%
    % assemblerIGA_PX   - This function computes the assembled matrices
    %                     relative to the coupling of pressure and x component of velocity considering 
    %                     IGA modal basis.
    %
        % The inputs are:
    %%
    %   (1)  imb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (2)  kmb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (3)  numbControlPts    : Size of the Functional Basis. Equal to the 
    %                            Dimention of the Spline Basis (Number of
    %                            Control Points)
    %   (4)  verGLWeights      : Vector Containing the Weights of the Quadrature 
    %                            Nodes in the Whole Mesh 
    %   (5)  mb_i              : Modal Basis referring for coeff. (imb,kmb)
    %                            for velocity
    %   (6)  mb_yi             : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %                            for velocity
    %   (7)  mb_k              : Modal Basis referring to coeff. (imb,kmb)
    %                            for pressure
    %   (8) LU                 : Thickness of the Channel Used in the Example
    %   (9) Computed           : Structure containing the coefficiets mu  and the force exciting the system 
    %   (10) horStep           : Step of the Finite Element Mesh
    %   (11) pU                : Degree of the Polynomial B-Spline Basis
    %                            for the velocity
    %   (12) pP                : Degree of the Polynomial B-Spline Basis
    %                            for the pressure
    %   (13) kU                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the velocity
    %   (14) kP                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the pressure
    %   (15) domainLeftLimit,
    %        domainRightLimit  : Extremes of the domain
    %   (16) jacIGA            : Jacobian computed at the quadrature nodes
    %   (17) numbHordQuadNodes : Number of quadrature nodes in the
    %                            horizontal direction
    %   (18) nKnots            : Number of knots for the isogeometric mesh
    %   (19) horWeights        : Horizontal Gauss Legendre integration
    %                            nodes
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis

    %% IMPORT CLASS
    
    import Core.AssemblerStokesHandler
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
    
    D1 = 1;
    D2 = 0;

    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'rho_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = ones(size(LU));
    obj_integrate_2.funcToIntegrate = ones(size(LU));
    
    obj_integrate_1.funcWeight = mb_k  .*mb_yi  .*verGLWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_i .*verGLWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    
    r00 = D2*lambda_1;
    r10 = D1*lambda_2;


    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    %Isogeometric basis for the velocity
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = pU;
    objBasisIGA.continuityParameter   = kU;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGAU,derBasisIGAU,controlPtsU] = newIsoGeoBasis(objBasisIGA);
    
    %Isogeometric basis for the pressure
    objBasisIGA.degreePolySplineBasis = pP;
    objBasisIGA.continuityParameter   = kP;
    
    [basisIGAP,~,~] = newIsoGeoBasis(objBasisIGA);
    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    row   = zeros(1,numbControlPtsU*numbControlPtsP);
    col   = zeros(1,numbControlPtsU*numbControlPtsP);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        deltaElem1 = zeros(pU+1,pP+1);
        deltaElem2 = zeros(pU+1,pP+1);
        
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

            shapeU  = basisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            dshapeU = derBasisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            cpiU    = controlPtsU((kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            
            shapeP  = basisIGAP((iel-1) * numbHorQuadNodes + igauss,(kP-1) * (iel-1) + iel + pP:-1:(kP-1) * (iel-1) + iel);

            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshapeU*cpiU';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshapeU = dshapeU/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shapeU' * shapeP) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (dshapeU'* shapeP) * gwt;
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:pU+1
            
            i = (kU-1)*(iel-1) + iel + pU - (a - 1);
            
                for b = 1:pP+1
                    
                    j = (kP-1)*(iel-1) + iel + pP - (b - 1);
                    icount = icount + 1;
                    
                    row(icount)   = i;
                    col(icount)   = j;
                    
                    deltaLocal1(icount) = deltaElem1(a,b);
                    deltaLocal2(icount) = deltaElem2(a,b);                    
                end
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPtsU,numbControlPtsP);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPtsU,numbControlPtsP);
    
    Al   = Delta1 + Delta2;
    
    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%
    
    Al(1,1) = 0;
    Al(1,2) = 0;
    Al(1,3) = 0;

    %disp('Finished ASSEMBLING LOOP');
end


%% Method 'assemblerIGA_PY'
            
function Al = assemblerIGA_PY(~,~,numbControlPtsU,numbControlPtsP, ...
                                verGLWeights,mb_i,mb_yi,mb_k, ...
                                LU, ...
                                horStep,pU,pP,kU,kP, ...
                                domainLeftLimit,domainRightLimit, ...
                                jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights)

    %%
    % assemblerIGA_PY   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
            % The inputs are:
    %%
    %   (1)  imb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (2)  kmb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (3)  numbControlPts    : Size of the Functional Basis. Equal to the 
    %                            Dimention of the Spline Basis (Number of
    %                            Control Points)
    %   (4)  verGLWeights      : Vector Containing the Weights of the Quadrature 
    %                            Nodes in the Whole Mesh 
    %   (5)  mb_i              : Modal Basis referring for coeff. (imb,kmb)
    %                            for velocity
    %   (6)  mb_yi             : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %                            for velocity
    %   (7)  mb_k              : Modal Basis referring to coeff. (imb,kmb)
    %                            for pressure
    %   (8) LU                 : Thickness of the Channel Used in the Example
    %   (9) Computed           : Structure containing the coefficiets mu  and the force exciting the system 
    %   (10) horStep           : Step of the Finite Element Mesh
    %   (11) pU                : Degree of the Polynomial B-Spline Basis
    %                            for the velocity
    %   (12) pP                : Degree of the Polynomial B-Spline Basis
    %                            for the pressure
    %   (13) kU                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the velocity
    %   (14) kP                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the pressure
    %   (15) domainLeftLimit,
    %        domainRightLimit  : Extremes of the domain
    %   (16) jacIGA            : Jacobian computed at the quadrature nodes
    %   (17) numbHordQuadNodes : Number of quadrature nodes in the
    %                            horizontal direction
    %   (18) nKnots            : Number of knots for the isogeometric mesh
    %   (19) horWeights        : Horizontal Gauss Legendre integration
    %                            nodes
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis

    %% IMPORT CLASS
    
    import Core.AssemblerStokesHandler
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

    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = ones(size(LU));
    obj_integrate_2.funcToIntegrate = ones(size(LU));
    
    obj_integrate_1.funcWeight = mb_k  .*mb_yi  .*verGLWeights;
    obj_integrate_2.funcWeight = mb_k  .*mb_i .*verGLWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    
    r00 = D2*lambda_1;
    r10 = D1*lambda_2;


    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    %Isogeometric basis for the velocity
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = pU;
    objBasisIGA.continuityParameter   = kU;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGAU,derBasisIGAU,controlPtsU] = newIsoGeoBasis(objBasisIGA);
    
    %Isogeometric basis for the pressure
    objBasisIGA.degreePolySplineBasis = pP;
    objBasisIGA.continuityParameter   = kP;
    
    [basisIGAP,~,~] = newIsoGeoBasis(objBasisIGA);
    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    row   = zeros(1,numbControlPtsU*numbControlPtsP);
    col   = zeros(1,numbControlPtsU*numbControlPtsP);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        deltaElem1 = zeros(pU+1,pP+1);
        deltaElem2 = zeros(pU+1,pP+1);
        
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

            shapeU  = basisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            dshapeU = derBasisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            cpiU    = controlPtsU((kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            
            shapeP  = basisIGAP((iel-1) * numbHorQuadNodes + igauss,(kP-1) * (iel-1) + iel + pP:-1:(kP-1) * (iel-1) + iel);

            
            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshapeU*cpiU';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshapeU = dshapeU/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shapeU' * shapeP) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (dshapeU'* shapeP) * gwt;
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:pU+1
            
            i = (kU-1)*(iel-1) + iel + pU - (a - 1);

                for b = 1:pP+1
                    
                    j = (kP-1)*(iel-1) + iel + pP - (b - 1);
                    icount = icount + 1;
                    
                    row(icount)   = i;
                    col(icount)   = j;
                    
                    deltaLocal1(icount) = deltaElem1(a,b);
                    deltaLocal2(icount) = deltaElem2(a,b);                    
                end
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPtsU,numbControlPtsP);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPtsU,numbControlPtsP);
    
    Al   = Delta1 + Delta2;
    
    %% IMPOSITION OF DIRICHLET BOUNDARY CONDITION
    %---------------------------------------------------------------------%
    % This step is necessary to guarantee that the resulting linear system
    % is not singular.
    %---------------------------------------------------------------------%
    
    Al(1,1)  = 0;
    Al(1,2)  = 0;
    Al(1,3)  = 0;

    %disp('Finished ASSEMBLING LOOP');
end


%% Method 'assemblerIGA_QX'
            
function [Al,bl,aLift,bLift] = assemblerIGA_QX(~,~,numbControlPtsP,numbControlPtsU, ...
                                verGLWeights,mb_i, mb_k,mb_yk, ...
                                L, ...
                                horStep,pP,kP,pU,kU, ...
                                domainLeftLimit,domainRightLimit, ...
                                ~,aLift,bLift,jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights)

    %%
    % assemblerIGA_QX   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
            % The inputs are:
    %%
    %   (1)  imb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (2)  kmb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (3)  numbControlPts    : Size of the Functional Basis. Equal to the 
    %                            Dimention of the Spline Basis (Number of
    %                            Control Points)
    %   (4)  verGLWeights      : Vector Containing the Weights of the Quadrature 
    %                            Nodes in the Whole Mesh 
    %   (5)  mb_i              : Modal Basis referring for coeff. (imb,kmb)
    %                            for pressure
    %   (6)  mb_k              : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %                            for velocity
    %   (7)  mb_yk              : Modal Basis referring to coeff. (imb,kmb)
    %                            for velocity
    %   (8) L                 : Thickness of the Channel Used in the Example
    %   (9) Computed           : Structure containing the coefficiets mu  and the force exciting the system 
    %   (10) horStep           : Step of the Finite Element Mesh
    %   (11) pU                : Degree of the Polynomial B-Spline Basis
    %                            for the velocity
    %   (12) pP                : Degree of the Polynomial B-Spline Basis
    %                            for the pressure
    %   (13) kU                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the velocity
    %   (14) kP                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the pressure
    %   (15) domainLeftLimit,
    %        domainRightLimit  : Extremes of the domain
    %   (16) aLift, bLift      : Coefficients of lifting
    %   (17) jacIGA            : Jacobian computed at the quadrature nodes
    %   (18) numbHordQuadNodes : Number of quadrature nodes in the
    %                            horizontal direction
    %   (19) nKnots            : Number of knots for the isogeometric mesh
    %   (20) horWeights        : Horizontal Gauss Legendre integration
    %                            nodes
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final rhs vector
    

    %% IMPORT CLASS
    
    import Core.AssemblerStokesHandler
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
    
    D1 = 1;
    D2 = 0;
    
    %% EXCITATION FORCE
    %---------------------------------------------------------------------%
    % Computation of the excitation force condiering the boundary
    % conditions acting on the system. The weights used to perform the
    % integration are the Gauss-Legendre nodes used in the horizontal
    % mesh.
    %---------------------------------------------------------------------%

    objForceInt = IntegrateHandler();
    
    objForceInt.funcToIntegrate =   -L.*(aLift);
                                
    objForceInt.funcWeight = mb_i.*verGLWeights;
    
    forceVec  = integrate(objForceInt);
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = ones(size(L));
    obj_integrate_2.funcToIntegrate = ones(size(L));
    
    obj_integrate_1.funcWeight = mb_i  .*mb_yk  .*verGLWeights;
    obj_integrate_2.funcWeight = mb_i  .*mb_k .*verGLWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    
    r00 = D2*lambda_1;
    r10 = D1*lambda_2;


    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    %Isogeometric basis for the velocity
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = pU;
    objBasisIGA.continuityParameter   = kU;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGAU,derBasisIGAU,controlPtsU] = newIsoGeoBasis(objBasisIGA);
    
    %Isogeometric basis for the pressure
    objBasisIGA.degreePolySplineBasis = pP;
    objBasisIGA.continuityParameter   = kP;
    
    [basisIGAP,~,~] = newIsoGeoBasis(objBasisIGA);
    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    globalForce  = zeros(numbControlPtsP,1);
    row   = zeros(1,numbControlPtsU*numbControlPtsP);
    col   = zeros(1,numbControlPtsU*numbControlPtsP);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        forceElem  = zeros(pP+1,1);
        
        deltaElem1 = zeros(pP+1,pU+1);
        deltaElem2 = zeros(pP+1,pU+1);
        
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

            shapeU  = basisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            dshapeU = derBasisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            cpiU    = controlPtsU((kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            
            shapeP  = basisIGAP((iel-1) * numbHorQuadNodes + igauss,(kP-1) * (iel-1) + iel + pP:-1:(kP-1) * (iel-1) + iel);
            
            
            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshapeU*cpiU';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshapeU = dshapeU/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL RHS VECTOR
            %-------------------------------------------------------------%
            % The current code uses 'numbHorQuadNodes' instead of
            % 'numbVerQuadNodes' in the next 5 lines of code. Check if
            % the solution is consistent.
            %-------------------------------------------------------------%

            forceElem       = forceElem  + forceVec((iel-1)*numbHorQuadNodes+igauss)*gwt;    % Force
            
            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shapeP' * shapeU) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (shapeP'* dshapeU) * gwt;
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:pP+1
            
            i = (kP-1)*(iel-1) + iel + pP - (a - 1);
            globalForce(i) = globalForce(i) + forceElem(a);

            for b = 1:pU+1

                j = (kU-1)*(iel-1) + iel + pU - (b - 1);
                icount = icount + 1;

                row(icount)   = i;
                col(icount)   = j;

                deltaLocal1(icount) = deltaElem1(a,b);
                deltaLocal2(icount) = deltaElem2(a,b);                    
            end

        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPtsP,numbControlPtsU);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPtsP,numbControlPtsU);
    
    Al   = Delta1 + Delta2;
    bl   = globalForce;
end 

%% Method 'assemblerIGA_QY'
            
function [Al,bl,aLift,bLift] = assemblerIGA_QY(~,~,numbControlPtsP,numbControlPtsU, ...
                                verGLWeights,mb_i,mb_k,mb_yk, ...
                                L, ...
                                horStep,pP,kP,pU,kU, ...
                                domainLeftLimit,domainRightLimit, ...
                                ~,aLift,bLift,jacIGA,...
                                numbHorQuadNodes,~,...
                                nknots,horWeights)

    %%
    % assemblerIGA_QY   - This function computes the assembled matrices
    %                  relative to the variational problem considering 
    %                  IGA modal basis.
    %
        % The inputs are:
    %%
    %   (1)  imb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (2)  kmb               : Coefficients Used in the Assembling Operation
    %                            Loop
    %   (3)  numbControlPts    : Size of the Functional Basis. Equal to the 
    %                            Dimention of the Spline Basis (Number of
    %                            Control Points)
    %   (4)  verGLWeights      : Vector Containing the Weights of the Quadrature 
    %                            Nodes in the Whole Mesh 
    %   (5)  mb_i              : Modal Basis referring for coeff. (imb,kmb)
    %                            for pressure
    %   (6)  mb_k              : Modal Basis in the Y Direction for coeff. (imb,kmb)
    %                            for velocity
    %   (7)  mb_yk             : Modal Basis referring to coeff. (imb,kmb)
    %                            for velocity
    %   (8) L                  : Thickness of the Channel Used in the Example
    %   (9) Computed           : Structure containing the coefficiets mu  and the force exciting the system 
    %   (10) horStep           : Step of the Finite Element Mesh
    %   (11) pU                : Degree of the Polynomial B-Spline Basis
    %                            for the velocity
    %   (12) pP                : Degree of the Polynomial B-Spline Basis
    %                            for the pressure
    %   (13) kU                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the velocity
    %   (14) kP                : Degree of Continuity of the Basis
    %                            'C^(p-k)' for the pressure
    %   (15) domainLeftLimit,
    %        domainRightLimit  : Extremes of the domain
    %   (16) aLift, bLift      : Coefficients of lifting
    %   (17) jacIGA            : Jacobian computed at the quadrature nodes
    %   (18) numbHordQuadNodes : Number of quadrature nodes in the
    %                            horizontal direction
    %   (19) nKnots            : Number of knots for the isogeometric mesh
    %   (20) horWeights        : Horizontal Gauss Legendre integration
    %                            nodes
    % 
    %
    % The outputs are:
    %%
    %   (1) Al    : Final Assembled Block Matrix Using IGA Basis
    %   (2) bl    : Final rhs vector

    %% IMPORT CLASS
    
    import Core.AssemblerStokesHandler
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
    
    objForceInt.funcToIntegrate =   -L.*(aLift);
                                
    objForceInt.funcWeight = mb_i.*verGLWeights;
    
    forceVec  = integrate(objForceInt);

    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    obj_integrate_1 = IntegrateHandler();
    obj_integrate_2 = IntegrateHandler();
    
    obj_integrate_1.funcToIntegrate = ones(size(L));
    obj_integrate_2.funcToIntegrate = ones(size(L));
    
    obj_integrate_1.funcWeight = mb_i  .*mb_yk  .*verGLWeights;
    obj_integrate_2.funcWeight = mb_i  .*mb_k .*verGLWeights;
    
    lambda_1   = integrate(obj_integrate_1);
    lambda_2   = integrate(obj_integrate_2);
    
    r00 = D2*lambda_1;
    r10 = D1*lambda_2;


    %% ISOGEOMETRIC BASIS AND DERIVATIVES
    %---------------------------------------------------------------------%
    % Precomputes the value of the shape function and its first derivative
    % evaluated in the Gauss points. Note that the shape functions
    % mentioned here are the functions named by in the reference
    % papers.
    %---------------------------------------------------------------------%
    
    %Isogeometric basis for the velocity
    objBasisIGA = BasisHandler();
    
    objBasisIGA.degreePolySplineBasis = pU;
    objBasisIGA.continuityParameter   = kU;
    objBasisIGA.numbElements          = nknots;
    objBasisIGA.numbQuadPointPerElem  = numbHorQuadNodes;
    objBasisIGA.leftBoundDomainX      = domainLeftLimit;
    objBasisIGA.rightBoundDomainX     = domainRightLimit;
    
    [basisIGAU,derBasisIGAU,controlPtsU] = newIsoGeoBasis(objBasisIGA);
    %Isogeometric basis for the pressure
    objBasisIGA.degreePolySplineBasis = pP;
    objBasisIGA.continuityParameter   = kP;
    
    [basisIGAP,~,~] = newIsoGeoBasis(objBasisIGA);

    %% LOOP OF ASSEMBLING LOCAL ELEMENTS
    %---------------------------------------------------------------------%
    % This loop computes the local coefficients (diffusion, advection and
    % reaction) and the local force components necessary to assemble the
    % global matrices. 
    %---------------------------------------------------------------------%
    
    % Vector Initialization
    
    globalForce  = zeros(numbControlPtsP,1);
    row   = zeros(1,numbControlPtsU*numbControlPtsP);
    col   = zeros(1,numbControlPtsU*numbControlPtsP);
    icount = 0;
    
    for iel = 1:nknots
        
        dl = horStep/2;
        
        % MEMORY ALLOCATION FOR BILINEAR COEFFICIENTS
        
        forceElem  = zeros(pP+1,1);
        
        deltaElem1 = zeros(pP+1,pU+1);
        deltaElem2 = zeros(pP+1,pU+1);
        
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

            shapeU  = basisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            dshapeU = derBasisIGAU((iel-1) * numbHorQuadNodes + igauss,(kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            cpiU    = controlPtsU((kU-1) * (iel-1) + iel + pU:-1:(kU-1) * (iel-1) + iel);
            
            shapeP  = basisIGAP((iel-1) * numbHorQuadNodes + igauss,(kP-1) * (iel-1) + iel + pP:-1:(kP-1) * (iel-1) + iel);

            % JACOBIAN EVALUATED AT THE GAUSS POINT
            
            J      = dshapeU*cpiU';

            % SHAPE FUNCTION DERIVATIVE WITH RESPECT TO THE PHYSICAL DOMAIN

            dshapeU = dshapeU/(J*jacIGA((iel-1)*numbHorQuadNodes+igauss));

            % INTEGRATION WEIGHTS IN THE PHYSICAL DOMAIN

            gwt = horWeights(igauss)*(J*jacIGA((iel-1)*numbHorQuadNodes+igauss))*dl;

            % LOCAL RHS VECTOR
            %-------------------------------------------------------------%
            % The current code uses 'numbHorQuadNodes' instead of
            % 'numbVerQuadNodes' in the next 5 lines of code. Check if
            % the solution is consistent.
            %-------------------------------------------------------------%

            forceElem       = forceElem  + forceVec((iel-1)*numbHorQuadNodes+igauss)*gwt;    % Force
            
            % LOCAL COMPONENTS OF THE STIFFNESS MATRIX
            
            deltaElem1 = deltaElem1 + r00((iel-1)*numbHorQuadNodes+igauss) * (shapeP' * shapeU) * gwt;
            deltaElem2 = deltaElem2 + r10((iel-1)*numbHorQuadNodes+igauss) * (shapeP'* dshapeU) * gwt;
        end

        % LOOP OF ASSEMBLING GLOBAL ELEMENTS
        %-----------------------------------------------------------------%
        % The coeeficients computed in the previous loops are assigned to
        % the correct position in the global matrix to assemble the global
        % stiffeness and the gloabl rhs matrix.
        %-----------------------------------------------------------------%

        for a = 1:pP+1
            
            i = (kP-1)*(iel-1) + iel + pP - (a - 1);
            
            globalForce(i) = globalForce(i) + forceElem(a);

            for b = 1:pU+1

                j = (kU-1)*(iel-1) + iel + pU - (b - 1);
                icount = icount + 1;

                row(icount)   = i;
                col(icount)   = j;

                deltaLocal1(icount) = deltaElem1(a,b);
                deltaLocal2(icount) = deltaElem2(a,b);                    
            end
            
        end
    end

    row = row(1:icount);
    col = col(1:icount);
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    Delta1 = sparse(row,col,deltaLocal1,numbControlPtsP,numbControlPtsU);
    Delta2 = sparse(row,col,deltaLocal2,numbControlPtsP,numbControlPtsU);
    
    Al   = Delta1 + Delta2;
    bl   = globalForce;
    
end 