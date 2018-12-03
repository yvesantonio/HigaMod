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
        
        domainProfile;              % Symbolic Function Defining the Profile of the
                                    % Simulation Domain
                      
        domainProfileDer;           % Symbolic Function Defining the Derivative of
                                    % the Profile of the Simulation Domain   
        
        numbHorQuadNodes;           % Number of horizontal nodes to apply the quadrature
                                    % formula
                                    
        numbVerQuadNodes;           % Number of vertical nodes to apply the quadrature
                                    % formula
                                    
        dimModalBasisP              % Dimension of the modal basis for the pressure
        
        dimModalBasisUx             % Dimension of the modal basis for the X component
                                    % of the velocity field
                                    
        dimModalBasisUy             % Dimension of the modal basis for the Y component
                                    % of the velocity field
                                    
        label_upBoundDomainP        % Tag of the boundary condition on the top lateral boundary
                                    % for the pressure
        
        label_downBoundDomainP      % Tag of the boundary condition on the bottom lateral boundary
                                    % for the pressure
                                    
        label_upBoundDomainUx       % Tag of the boundary condition on the top lateral boundary
                                    % for the pressure
        
        label_downBoundDomainUx     % Tag of the boundary condition on the bottom lateral boundary
                                    % for the pressure
                                    
        label_upBoundDomainUy       % Tag of the boundary condition on the top lateral boundary
                                    % for the pressure
        
        label_downBoundDomainUy     % Tag of the boundary condition on the bottom lateral boundary
                                    % for the pressure
                                    
        discStruct      % Structure containing the discretization parameters
        
        boundCondStruct % Structure containing the boundary condition 
                        % information

        igaBasisStruct  % Structure containing the isogeometric basis
                        % paramters

        probParameters  % Structure containing the problem parameters

        timeStruct      % Structure containing the time domain and time
                        % simulation parameters

        quadProperties  % Structure containing the quadrature parameters
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
        
            %% Method 'buildIGAScatter'
            
            function [stiffStruct,massStruct,forceStruct] = buildSystemIGAScatter(obj)
                                                         
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE QUADRATURE NODES %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Horizontal Direction
            
            numbHorNodesP = obj.quadProperties.numbHorNodesP;
            
            % Vertical Direction
            
            numbVerNodesP = obj.quadProperties.numbVerNodesP;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X QUADRATURE NODES %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Horizontal Direction
            
            numbHorNodesUx = obj.quadProperties.numbHorNodesUx;
            
            % Vertical Direction
            
            numbVerNodesUx = obj.quadProperties.numbVerNodesUx;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y QUADRATURE NODES %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Horizontal Direction
            
            numbHorNodesUy = obj.quadProperties.numbHorNodesUx;
            
            % Vertical Direction
            
            numbVerNodesUy = obj.quadProperties.numbVerNodesUx;
            
            %% NUMBER OF KNOTS - IDOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of knots used to divide the spline curve in the
            % isogeometric analysis.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE NUMBER KNOTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            numbKnotsP  = obj.discStruct.numbElementsP;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X NUMBER KNOTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            numbKnotsUx  = obj.discStruct.numbElementsUx;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y NUMBER KNOTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            numbKnotsUy  = obj.discStruct.numbElementsUy;
            
            %% NUMBER OF CONTROL POINTS - ISOGEOMETRIC ANALYSIS
            %-------------------------------------------------------------%
            % Total number of control points used in the isogeometric
            % approach. Note that this number is equal to the dimension of
            % the basis used to define the Spline curve and that the
            % formula used bellow corresponds to use a straight line as
            % reference supporting fiber.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE NUMBER OF CONTROL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            degreePolySplineBasisP  = obj.igaBasisStruct.degreeSplineBasisP;
            continuityParameterP    = obj.igaBasisStruct.continuityParameterP; 
            numbControlPtsP         = numbKnotsP*(degreePolySplineBasisP - continuityParameterP) + ...
                                      1 + continuityParameterP;
                          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X NUMBER OF CONTROL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            degreePolySplineBasisUx = obj.igaBasisStruct.degreeSplineBasisUx;
            continuityParameterUx   = obj.igaBasisStruct.continuityParameterUx; 
            numbControlPtsUx        = numbKnotsUx * (degreePolySplineBasisUx - continuityParameterUx) + ...
                                      1 + continuityParameterUx;
                          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y NUMBER OF CONTROL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            degreePolySplineBasisUy = obj.igaBasisStruct.degreeSplineBasisUy;
            continuityParameterUy   = obj.igaBasisStruct.continuityParameterUy; 
            numbControlPtsUy        = numbKnotsUy * (degreePolySplineBasisUy - continuityParameterUy) + ...
                                      1 + continuityParameterUy;

            %% GAUSS-LEGENDRE INTEGRATION NODES
            %-------------------------------------------------------------%
            % The following method reveives the number of integration nodes
            % in the horizontal and vertical direction and retruns the
            % respective Gauss-Legendre nodes and weigths for the
            % integration interval [0,1].   
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE MODAL GL INTEGRATION NODES %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_gaussLegendreP = IntegrateHandler();
            obj_gaussLegendreP.numbQuadNodes = numbVerNodesP;
            [~, verGLNodesP, verWeightsP] = gaussLegendre(obj_gaussLegendreP);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X MODAL GL INTEGRATION NODES %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_gaussLegendreUx = IntegrateHandler();
            obj_gaussLegendreUx.numbQuadNodes = numbVerNodesUx;
            [~, verGLNodesUx, verWeightsUx] = gaussLegendre(obj_gaussLegendreUx);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y MODAL GL INTEGRATION NODES %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_gaussLegendreUy = IntegrateHandler();
            obj_gaussLegendreUy.numbQuadNodes = numbVerNodesUy;
            [~, verGLNodesUy, verWeightsUy] = gaussLegendre(obj_gaussLegendreUy);
            
            %% EXTRACT GEOMETRIC INFORMATION
            
            geometry  = obj.geometricInfo.geometry;
            map       = @(x,y) obj.geometricInfo.map(x,y);
            Jac       = @(x,y) obj.geometricInfo.Jac(x,y);
            Hes       = @(x,y) obj.geometricInfo.Hes(x,y);
            
            %% COMPUTE REFERENCE KNOT AND CONTROL POINTS
            %-------------------------------------------------------------%
            % Note: Create the knots and control points of the reference
            % supporting fiber where the coupled 1D problems will be
            % solved.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE REFERENCE KNOTS AND CONTROL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsubP = numbKnotsP;
            degreeP = degreePolySplineBasisP;
            regularityP = continuityParameterP;
            
            [knotsP, zetaP] = kntrefine (refDomain1D.nurbs.knots, nsubP-1, degreeP, regularityP);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X REFERENCE KNOTS AND CONTROL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsubUx = numbKnotsUx;
            degreeUx = degreePolySplineBasisUx;
            regularityUx = continuityParameterUx;
            
            [knotsUx, zetaUx] = kntrefine (refDomain1D.nurbs.knots, nsubUx-1, degreeUx, regularityUx);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y REFERENCE KNOTS AND CONTROL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            refDomain1D = geo_load(nrbline ([0 0], [1 0]));
            
            nsubUy = numbKnotsUy;
            degreeUy = degreePolySplineBasisUy;
            regularityUy = continuityParameterUy;
            
            [knotsUy, zetaUy] = kntrefine (refDomain1D.nurbs.knots, nsubUy-1, degreeUy, regularityUy);
            
            %% GENERATE ISOGEOMETRIC MESH FUNCTION
            %-------------------------------------------------------------%
            % Note: Use the number of quadrature nodes to generate the
            % quadrature rule using the gauss nodes. Then, use this
            % information with the computed knots and control points to
            % generate the isogeometric mesh of the centreline.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%
            % PRESSURE MESH %
            %%%%%%%%%%%%%%%%%
            
            ruleP       = msh_gauss_nodes (numbHorNodesP);
            [qnP, qwP]  = msh_set_quad_nodes (zetaP, ruleP);
            mshP        = msh_cartesian (zetaP, qnP, qwP, refDomain1D);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X MESH %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ruleUx       = msh_gauss_nodes (numbHorNodesUx);
            [qnUx, qwUx] = msh_set_quad_nodes (zetaUx, ruleUx);
            mshUx        = msh_cartesian (zetaUx, qnUx, qwUx, refDomain1D);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y  MESH %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ruleUy       = msh_gauss_nodes (numbHorNodesUy);
            [qnUy, qwUy] = msh_set_quad_nodes (zetaUy, ruleUy);
            mshUy        = msh_cartesian (zetaUy, qnUy, qwUy, refDomain1D);
            
            %% CONSTRUCT THE ISOGEOMETRIC FUNCTIONAL SPACE
            %-------------------------------------------------------------%
            % Note: Construct the functional space used to discretize the
            % problem along the centreline and perform the isogeometric
            % analysis. 
            % Here, we also evaluate the functional spaces in the specific
            % element. This will be used in the assembling loop to compute
            % the operators using the basis functions.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE ISOGEOMETRIC FUNCTIONAL SPACE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            spaceP    = sp_bspline (knotsP, degreeP, mshP);
            
            numbKnotsP = mshP.nel_dir;
            
            for iel = 1:numbKnotsP
                
                msh_colP = msh_evaluate_col (mshP, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFuncP = sp_evaluate_col (spaceP, msh_colP, 'value', false, 'gradient', true);
                shapFuncP = sp_evaluate_col (spaceP, msh_colP);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFuncP{iel,1} = shapFuncP;
                spaceFuncP{iel,2} = gradFuncP;
                spaceFuncP{iel,3} = msh_colP;

            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X ISOGEOMETRIC FUNCTIONAL SPACE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            spaceUx    = sp_bspline (knotsUx, degreeUx, mshUx);
            
            numbKnotsUx = mshUx.nel_dir;
            
            for iel = 1:numbKnotsUx
                
                msh_colUx = msh_evaluate_col (mshUx, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFuncUx = sp_evaluate_col (spaceUx, msh_colUx, 'value', false, 'gradient', true);
                shapFuncUx = sp_evaluate_col (spaceUx, msh_colUx);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFuncUx{iel,1} = shapFuncUx;
                spaceFuncUx{iel,2} = gradFuncUx;
                spaceFuncUx{iel,3} = msh_colUx;

            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y ISOGEOMETRIC FUNCTIONAL SPACE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            spaceUy    = sp_bspline (knotsUy, degreeUy, mshUy);
            
            numbKnotsUy = mshUy.nel_dir;
            
            for iel = 1:numbKnotsUy
                
                msh_colUy = msh_evaluate_col (mshUy, iel);
                
                %---------------------------------------------------------%
                % Evaluated space to compute the operators
                %---------------------------------------------------------%
                
                gradFuncUy = sp_evaluate_col (spaceUy, msh_colUy, 'value', false, 'gradient', true);
                shapFuncUy = sp_evaluate_col (spaceUy, msh_colUy);
                
                %---------------------------------------------------------%
                % Store the spaces in the cell matrix
                %---------------------------------------------------------%
                
                spaceFuncUy{iel,1} = shapFuncUy;
                spaceFuncUy{iel,2} = gradFuncUy;
                spaceFuncUy{iel,3} = msh_colUy;

            end
            
            %% AUGMENTED HORIZONTAL POINTS 
            %-------------------------------------------------------------%
            % Note: The FOR loop runs over the knots in the centerline ans
            % compute the required quadrature nodes. The quadrature nodes
            % will be necessary to evaluate the bilinear coefficients and
            % the forcing component.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE AUGMENTED HORIZONTAL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            horEvalNodesP = zeros(numbKnotsP * numbHorNodesP,1);
            
            for iel = 1:numbKnotsP
                
                msh_colP = msh_evaluate_col (mshP, iel);
                
                localNodesP = reshape (msh_colP.geo_map(1,:,:), msh_colP.nqn, msh_colP.nel);  
                horEvalNodesP((iel - 1)*numbHorNodesP + 1 : iel*numbHorNodesP) = localNodesP;  
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X AUGMENTED HORIZONTAL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            horEvalNodesUx = zeros(numbKnotsUx * numbHorNodesUx,1);
            
            for iel = 1:numbKnotsUx
                
                msh_colUx = msh_evaluate_col (mshUx, iel);
                
                localNodesUx = reshape (msh_colUx.geo_map(1,:,:), msh_colUx.nqn, msh_colUx.nel);  
                horEvalNodesUx((iel - 1)*numbHorNodesUx + 1 : iel*numbHorNodesUx) = localNodesUx;  
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y AUGMENTED HORIZONTAL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            horEvalNodesUy = zeros(numbKnotsUy * numbHorNodesUy,1);
            
            for iel = 1:numbKnotsUy
                
                msh_colUy = msh_evaluate_col (mshUy, iel);
                
                localNodesUy = reshape (msh_colUy.geo_map(1,:,:), msh_colUy.nqn, msh_colUy.nel);  
                horEvalNodesUy((iel - 1)*numbHorNodesUy + 1 : iel*numbHorNodesUy) = localNodesUy;  
                
            end
            
            %% AUGMENTED VERTICAL POINTS
            %-------------------------------------------------------------%
            % Note: Since we are using the full map from the physical
            % domain to the reference domain, then the whole modal analysis
            % has to be performed in the vertical direction in the interval
            % [0,1];
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE AUGMENTED VERTICAL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            objVertQuadRuleP = IntegrateHandler();

            objVertQuadRuleP.leftBoundInterval = 0;
            objVertQuadRuleP.rightBoundInterval = 1;
            objVertQuadRuleP.inputNodes = verGLNodesP;
            objVertQuadRuleP.inputWeights = verWeightsP;

            [augVerNodesP, augVerWeightsP] = quadratureRule(objVertQuadRuleP);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X AUGMENTED VERTICAL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            objVertQuadRuleUx = IntegrateHandler();

            objVertQuadRuleUx.leftBoundInterval = 0;
            objVertQuadRuleUx.rightBoundInterval = 1;
            objVertQuadRuleUx.inputNodes = verGLNodesUx;
            objVertQuadRuleUx.inputWeights = verWeightsUx;

            [augVerNodesUx, augVerWeightsUx] = quadratureRule(objVertQuadRuleUx);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y AUGMENTED VERTICAL POINTS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            objVertQuadRuleUy = IntegrateHandler();

            objVertQuadRuleUy.leftBoundInterval = 0;
            objVertQuadRuleUy.rightBoundInterval = 1;
            objVertQuadRuleUy.inputNodes = verGLNodesUy;
            objVertQuadRuleUy.inputWeights = verWeightsUy;

            [augVerNodesUy, augVerWeightsUy] = quadratureRule(objVertQuadRuleUy);
            
            %% COMPUTATION OF THE MODAL BASIS IN THE Y DIRECTION
            %-------------------------------------------------------------%
            % The next method takes the coefficients of the bilinear form,
            % the dimension of the modal basis (assigned in the demo) and
            % the vertical evaluation points and computes the coefficients
            % of the modal basis when evaluated in the points of the
            % transverse fiber.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE MODAL BASIS %
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_newModalBasisP = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            % obj_newModalBasis.dimLegendreBase = obj.dimModalBasisP;
            % obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasisP.dimModalBasis = obj.discStruct.numbModesP;
            obj_newModalBasisP.evalNodesY = verGLNodesP;
            obj_newModalBasisP.labelUpBoundCond = obj.boundCondStruct.bc_up_tag_P;
            obj_newModalBasisP.labelDownBoundCond = obj.boundCondStruct.bc_down_tag_P;
            obj_newModalBasisP.coeffForm = obj.probParameters;

            [modalBasisP, modalBasisDerP] = newModalBasisStokes(obj_newModalBasisP);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X MODAL BASIS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_newModalBasisUx = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            % obj_newModalBasis.dimLegendreBase = obj.dimModalBasisP;
            % obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasisUx.dimModalBasis = obj.discStruct.numbModesUx;
            obj_newModalBasisUx.evalNodesY = verGLNodesUx;
            obj_newModalBasisUx.labelUpBoundCond = obj.boundCondStruct.bc_up_tag_Ux;
            obj_newModalBasisUx.labelDownBoundCond = obj.boundCondStruct.bc_down_tag_Ux;
            obj_newModalBasisUx.coeffForm = obj.probParameters;

            [modalBasisUx, modalBasisDerUx] = newModalBasisStokes(obj_newModalBasisUx);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y MODAL BASIS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_newModalBasisUy = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            % obj_newModalBasis.dimLegendreBase = obj.dimModalBasisP;
            % obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasisUy.dimModalBasis = obj.discStruct.numbModesUy;
            obj_newModalBasisUy.evalNodesY = verGLNodesUy;
            obj_newModalBasisUy.labelUpBoundCond = obj.boundCondStruct.bc_up_tag_Uy;
            obj_newModalBasisUy.labelDownBoundCond = obj.boundCondStruct.bc_down_tag_Uy;
            obj_newModalBasisUy.coeffForm = obj.probParameters;

            [modalBasisUy, modalBasisDerUy] = newModalBasisStokes(obj_newModalBasisUy);
            
            %% MEMORY ALLOCATION FOR SYSTEM MATRICES
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BLOCK COMPONENTS OF THE MASS MATRIX %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Mxx = sparse( numbControlPts*(obj.dimModalBasisUx), numbControlPts*(obj.dimModalBasisUx) );
            Myy = sparse( numbControlPts*(obj.dimModalBasisUy), numbControlPts*(obj.dimModalBasisUy) );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BLOCK COMPONENTS OF THE STIFFNESS MATRIX %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Axx = sparse( numbControlPts*(obj.dimModalBasisUx), numbControlPts*(obj.dimModalBasisUx) );
            Bxy = sparse( numbControlPts*(obj.dimModalBasisUx), numbControlPts*(obj.dimModalBasisUy) );
            Byx = sparse( numbControlPts*(obj.dimModalBasisUy), numbControlPts*(obj.dimModalBasisUx) );
            Ayy = sparse( numbControlPts*(obj.dimModalBasisUy), numbControlPts*(obj.dimModalBasisUy) );
            Px  = sparse( numbControlPts*(obj.dimModalBasisUx), numbControlPts*(obj.dimModalBasisP) );
            Py  = sparse( numbControlPts*(obj.dimModalBasisUy), numbControlPts*(obj.dimModalBasisP) );
            Qx  = sparse( numbControlPts*(obj.dimModalBasisP), numbControlPts*(obj.dimModalBasisUx) );
            Qy  = sparse( numbControlPts*(obj.dimModalBasisP), numbControlPts*(obj.dimModalBasisUy) );            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BLOCK COMPONENTS OF THE SOURCE TERM %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Fx = zeros ( numbControlPts*(obj.dimModalBasisUx), 1);
            Fy = zeros ( numbControlPts*(obj.dimModalBasisUx), 1);
            Fp = zeros ( numbControlPts*(obj.dimModalBasisUx), 1);
            
            %% EVALUATION OF THE GEOMETRIC PROPERTIES OF THE DOMAIN
            %-------------------------------------------------------------%
            % We use the horizontal mesh created previously to evaluate the
            % thickness and the vertical coordinate of the centerline of
            % the channel.
            %-------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE PHYSICAL DOMAIN %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            verEvalNodesP = augVerNodesP;
            XP = mapOut(horEvalNodesP,verEvalNodesP,map,1);
            YP = mapOut(horEvalNodesP,verEvalNodesP,map,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG X PHYSICAL DOMAIN %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            verEvalNodesUx = augVerNodesUx;
            XUx = mapOut(horEvalNodesUx,verEvalNodesUx,map,1);
            YUx = mapOut(horEvalNodesUx,verEvalNodesUx,map,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VELOCITY FIELD ALONG Y PHYSICAL DOMAIN %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            verEvalNodesUy = augVerNodesUy;
            XUy = mapOut(horEvalNodesUy,verEvalNodesUy,map,1);
            YUy = mapOut(horEvalNodesUy,verEvalNodesUy,map,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CREATE GEOMETRIC STRUCTURE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            geoData.XP = XP;
            geoData.YP = YP;
            geoData.XUx = XUx;
            geoData.YUx = YUx;
            geoData.XUy = XUy;
            geoData.YUy = YUy;
            geoData.horNodesP = horEvalNodesP;
            geoData.verNodesP = verEvalNodesP;
            geoData.horNodesUx = horEvalNodesUx;
            geoData.verNodesUx = verEvalNodesUx;
            geoData.horNodesUy = horEvalNodesUy;
            geoData.verNodesUy = verEvalNodesUy;
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % JACOBIAN PROPERTIES %
            %%%%%%%%%%%%%%%%%%%%%%%
            
            evalJac     = jacOut(horEvalNodesP,verEvalNodesP,Jac);
            evalDetJac  = detJac(horEvalNodesP,verEvalNodesP,Jac);
            
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
            
            %% ASSEMBLE TIME DEPENDENT STRUCTURES
            
            for time = obj.timeStruct.timeDomain
                
                %% EVALUATION OF THE BILINEAR COEFFICIENTS OF THE DOMAIN
                %-------------------------------------------------------------%
                % We need to evaluate all the coefficients of the bilinear form
                % in the quadrature nodes along the vertical direction to
                % perform the first integral (along the transverse fiber) to
                % obtain a the coefficients of the 1D coupled problem. The
                % result will be coefficients as a function of 'x'.
                %-------------------------------------------------------------%

                % KINETIC VISCOSITY OF THE FLUID

                evalNu    = obj.probParameters.nu(XP,YP,time);

                %% EVALUATION OF THE EXCITING FORCE 
                %-------------------------------------------------------------%
                % Finally, we use the horizontal and vertical meshes to
                % evaluate the exciting force acting on the system in the whole
                % domain.
                %-------------------------------------------------------------%

                evalForce = obj.probParameters.force(XP,YP,time);

                %-------------------------------------------------------------%
                % Note: All of the coefficients of the bilinear form and the
                % exciting term evaluated in the domain, as well as the mesh of
                % vertical coordinates are stored in a data structure called
                % "Computed" to be treated in the assembler function.
                %-------------------------------------------------------------%

                Computed = struct('nu_c',evalNu,'force_c',evalForce,'y',verEvalNodesP);            

                %% AXX ASSEMBLING LOOP

                for kmb = 1:obj.discStruct.numbModesUx
                    for jmb = 1:obj.discStruct.numbModesUx

                        [Amb,bmb] = assemblerIGAScatterAxx( kmb, jmb, ...
                                                augVerWeights,modalBasis(:,kmb),modalBasisDer(:,kmb),modalBasis(:,jmb),...
                                                modalBasisDer(:,jmb),geoData,Computed,lifting,aLift,bLift,msh,space,jacFunc,...
                                                spaceFuncP,obj.dirCondFuncStruct.igaBoundCond);

                        % Assignment of the Block Matrix Just Assembled

                        A(1+(kmb-1)*numbControlPts : kmb*numbControlPts , 1+(jmb-1)*numbControlPts : jmb*numbControlPts) = Amb;

                        disp(['FINISHED ASSEMBLING LOOP (',num2str(kmb),' , ',num2str(jmb),')']);

                    end

                    % Assignment of the Block Vector Just Assembled
                    b( 1+(kmb-1)*numbControlPts : kmb*numbControlPts ) = bmb;
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
                
            end
            
    end
end

%% ASSEMBLER STOKES HANDLER - ASSEMBLING METHODS

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


    
end