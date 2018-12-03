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
            
            jacFunc.modGradPhi1         = Phi1_dx.^2 + jacFunc.Phi1_dy.^2;
            jacFunc.modGradPhi2         = Phi2_dx.^2 + jacFunc.Phi2_dy.^2;
            jacFunc.modGradPhi1_Proj1   = Phi1_dx.^2;
            jacFunc.modGradPhi1_Proj2   = Phi1_dy.^2;
            jacFunc.modGradPhi2_Proj1   = Phi2_dx.^2;
            jacFunc.modGradPhi2_Proj2   = Phi2_dy.^2;

            jacFunc.GPhi1DotGPhi2       = Phi1_dx .* Phi2_dx + Phi1_dy .* Phi2_dy;
            jacFunc.GPhi1DotGPhi2_Proj1 = Phi1_dx .* Phi2_dx;
            jacFunc.GPhi1DotGPhi2_Proj2 = Phi1_dy .* Phi2_dy;

            jacFunc.GPhi1Proj1DotGPhi1Proj2 = Phi1_dx .* Phi1_dy;
            jacFunc.GPhi2Proj1DotGPhi1Proj2 = Phi2_dx .* Phi1_dy;
            jacFunc.GPhi1Proj1DotGPhi2Proj2 = Phi1_dx .* Phi2_dy;
            jacFunc.GPhi2Proj1DotGPhi2Proj2 = Phi2_dx .* Phi2_dy;
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % STRCUT CONTAINING TIME DEPENDENT INFORMATION %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            stiffTimeStruct = cell(size(obj.timeStruct.timeDomain));
            forceTimeStruct = cell(size(obj.timeStruct.timeDomain));
            massTimeStruct  = cell(size(obj.timeStruct.timeDomain));
            
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

                evalForceX = obj.probParameters.force_x(XP,YP,time);
                evalForceY = obj.probParameters.force_y(XP,YP,time);

                %-------------------------------------------------------------%
                % Note: All of the coefficients of the bilinear form and the
                % exciting term evaluated in the domain, as well as the mesh of
                % vertical coordinates are stored in a data structure called
                % "Computed" to be treated in the assembler function.
                %-------------------------------------------------------------%

                Computed = struct('nu',evalNu,'forceX',evalForceX,'forceY',evalForceY,'y',verEvalNodesP);            

                %% Axx ASSEMBLING LOOP
                
                assemblerStruct = [];

                for kmb = 1:obj.discStruct.numbModesUx
                    for jmb = 1:obj.discStruct.numbModesUx
                        
                        assemblerStruct.index1     = kmb;
                        assemblerStruct.index2     = jmb;
                        assemblerStruct.wgh1       = augVerWeightsUx;
                        assemblerStruct.wgh2       = augVerWeightsUx;
                        assemblerStruct.mb1        = modalBasisUx;
                        assemblerStruct.mb2        = modalBasisUx;
                        assemblerStruct.dmb1       = modalBasisDerUx;
                        assemblerStruct.dmb2       = modalBasisDerUx;
                        assemblerStruct.param      = Computed;
                        assemblerStruct.geodata    = geoData;
                        assemblerStruct.jacFunc    = jacFunc;
                        assemblerStruct.msh1       = mshUx;
                        assemblerStruct.msh2       = mshUx;
                        assemblerStruct.space1     = spaceUx;
                        assemblerStruct.space2     = spaceUx;
                        assemblerStruct.spaceFunc1 = spaceFuncUx;
                        assemblerStruct.spaceFunc1 = spaceFuncUx;
                        assemblerStruct.aLift      = aLift;
                        assemblerStruct.bLift      = bLift;
                        assemblerStruct.lifting    = lifting;
                        
                        [Mx_mb,Ax_mb,Fx_mb] = assemblerIGAScatterAx(assemblerStruct);

                        % Assignment of the Block Matrix Just Assembled

                        Axx(1 + (kmb-1) * numbControlPtsUx : kmb * numbControlPtsUx , 1 + (jmb-1) * numbControlPtsUx : jmb * numbControlPtsUx) = Ax_mb;
                        Mxx(1 + (kmb-1) * numbControlPtsUx : kmb * numbControlPtsUx , 1 + (jmb-1) * numbControlPtsUx : jmb * numbControlPtsUx) = Mx_mb;

                        disp(['Ax - FINISHED ASSEMBLING LOOP (',num2str(kmb),' , ',num2str(jmb),')']);

                    end

                    % Assignment of the Block Vector Just Assembled
                    Fx( 1 + (kmb-1) * numbControlPtsUx : kmb * numbControlPtsUx ) = Fx_mb;
                end
                
                %% Ayy ASSEMBLING LOOP
                
                assemblerStruct = [];

                for kmb = 1:obj.discStruct.numbModesUy
                    for jmb = 1:obj.discStruct.numbModesUy
                        
                        assemblerStruct.index1     = kmb;
                        assemblerStruct.index2     = jmb;
                        assemblerStruct.wgh1       = augVerWeightsUy;
                        assemblerStruct.wgh2       = augVerWeightsUy;
                        assemblerStruct.mb1        = modalBasisUy;
                        assemblerStruct.mb2        = modalBasisUy;
                        assemblerStruct.dmb1       = modalBasisDerUy;
                        assemblerStruct.dmb2       = modalBasisDerUy;
                        assemblerStruct.param      = Computed;
                        assemblerStruct.geodata    = geoData;
                        assemblerStruct.jacFunc    = jacFunc;
                        assemblerStruct.msh1       = mshUy;
                        assemblerStruct.msh2       = mshUy;
                        assemblerStruct.space1     = spaceUy;
                        assemblerStruct.space2     = spaceUy;
                        assemblerStruct.spaceFunc1 = spaceFuncUy;
                        assemblerStruct.spaceFunc1 = spaceFuncUy;
                        assemblerStruct.aLift      = aLift;
                        assemblerStruct.bLift      = bLift;
                        assemblerStruct.lifting    = lifting;
                        
                        [My_mb,Ay_mb,Fy_mb] = assemblerIGAScatterAy(assemblerStruct);

                        % Assignment of the Block Matrix Just Assembled

                        Ayy(1 + (kmb-1) * numbControlPtsUy : kmb * numbControlPtsUy , 1 + (jmb-1) * numbControlPtsUy : jmb * numbControlPtsUy) = Ay_mb;
                        Myy(1 + (kmb-1) * numbControlPtsUy : kmb * numbControlPtsUy , 1 + (jmb-1) * numbControlPtsUy : jmb * numbControlPtsUy) = My_mb;

                        disp(['Ay - FINISHED ASSEMBLING LOOP (',num2str(kmb),' , ',num2str(jmb),')']);

                    end

                    % Assignment of the Block Vector Just Assembled
                    Fy( 1 + (kmb-1) * numbControlPtsUy : kmb * numbControlPtsUy ) = Fy_mb;
                end
                
                %% Bxy ASSEMBLING LOOP
                
                assemblerStruct = [];

                for kmb = 1:obj.discStruct.numbModesUx
                    for jmb = 1:obj.discStruct.numbModesUy
                        
                        assemblerStruct.index1     = kmb;
                        assemblerStruct.index2     = jmb;
                        assemblerStruct.wgh1       = augVerWeightsUx;
                        assemblerStruct.wgh2       = augVerWeightsUy;
                        assemblerStruct.mb1        = modalBasisUx;
                        assemblerStruct.mb2        = modalBasisUy;
                        assemblerStruct.dmb1       = modalBasisDerUx;
                        assemblerStruct.dmb2       = modalBasisDerUy;
                        assemblerStruct.param      = Computed;
                        assemblerStruct.geodata    = geoData;
                        assemblerStruct.jacFunc    = jacFunc;
                        assemblerStruct.msh1       = mshUx;
                        assemblerStruct.msh2       = mshUy;
                        assemblerStruct.space1     = spaceUx;
                        assemblerStruct.space2     = spaceUy;
                        assemblerStruct.spaceFunc1 = spaceFuncUx;
                        assemblerStruct.spaceFunc1 = spaceFuncUy;
                        assemblerStruct.aLift      = aLift;
                        assemblerStruct.bLift      = bLift;
                        assemblerStruct.lifting    = lifting;
                        
                        [Bxy_mb] = assemblerIGAScatterBxy(assemblerStruct);

                        % Assignment of the Block Matrix Just Assembled

                        Bxy(1 + (kmb-1) * numbControlPtsUx : kmb * numbControlPtsUx , 1 + (jmb-1) * numbControlPtsUy : jmb * numbControlPtsUy) = Bxy_mb;

                        disp(['Bxy - FINISHED ASSEMBLING LOOP (',num2str(kmb),' , ',num2str(jmb),')']);

                    end
                end
                
                %% Byx ASSEMBLING LOOP
                
                assemblerStruct = [];

                for kmb = 1:obj.discStruct.numbModesUy
                    for jmb = 1:obj.discStruct.numbModesUx
                        
                        assemblerStruct.index1     = kmb;
                        assemblerStruct.index2     = jmb;
                        assemblerStruct.wgh1       = augVerWeightsUy;
                        assemblerStruct.wgh2       = augVerWeightsUx;
                        assemblerStruct.mb1        = modalBasisUy;
                        assemblerStruct.mb2        = modalBasisUx;
                        assemblerStruct.dmb1       = modalBasisDerUy;
                        assemblerStruct.dmb2       = modalBasisDerUx;
                        assemblerStruct.param      = Computed;
                        assemblerStruct.geodata    = geoData;
                        assemblerStruct.jacFunc    = jacFunc;
                        assemblerStruct.msh1       = mshUy;
                        assemblerStruct.msh2       = mshUx;
                        assemblerStruct.space1     = spaceUy;
                        assemblerStruct.space2     = spaceUx;
                        assemblerStruct.spaceFunc1 = spaceFuncUy;
                        assemblerStruct.spaceFunc1 = spaceFuncUx;
                        assemblerStruct.aLift      = aLift;
                        assemblerStruct.bLift      = bLift;
                        assemblerStruct.lifting    = lifting;
                        
                        [Byx_mb] = assemblerIGAScatterByx(assemblerStruct);

                        % Assignment of the Block Matrix Just Assembled

                        Byx(1 + (kmb-1) * numbControlPtsUy : kmb * numbControlPtsUy , 1 + (jmb-1) * numbControlPtsUx : jmb * numbControlPtsUx) = Byx_mb;

                        disp(['Byx - FINISHED ASSEMBLING LOOP (',num2str(kmb),' , ',num2str(jmb),')']);

                    end
                end
                
                %% Px ASSEMBLING LOOP
                
                assemblerStruct = [];

                for kmb = 1:obj.discStruct.numbModesP
                    for jmb = 1:obj.discStruct.numbModesUx
                        
                        assemblerStruct.index1     = kmb;
                        assemblerStruct.index2     = jmb;
                        assemblerStruct.wgh1       = augVerWeightsP;
                        assemblerStruct.wgh2       = augVerWeightsUx;
                        assemblerStruct.mb1        = modalBasisP;
                        assemblerStruct.mb2        = modalBasisUx;
                        assemblerStruct.dmb1       = modalBasisDerP;
                        assemblerStruct.dmb2       = modalBasisDerUx;
                        assemblerStruct.param      = Computed;
                        assemblerStruct.geodata    = geoData;
                        assemblerStruct.jacFunc    = jacFunc;
                        assemblerStruct.msh1       = mshP;
                        assemblerStruct.msh2       = mshUx;
                        assemblerStruct.space1     = spaceP;
                        assemblerStruct.space2     = spaceUx;
                        assemblerStruct.spaceFunc1 = spaceFuncP;
                        assemblerStruct.spaceFunc1 = spaceFuncUx;
                        assemblerStruct.aLift      = aLift;
                        assemblerStruct.bLift      = bLift;
                        assemblerStruct.lifting    = lifting;
                        
                        [Px_mb] = assemblerIGAScatterPx(assemblerStruct);

                        % Assignment of the Block Matrix Just Assembled

                        Px(1 + (kmb-1) * numbControlPtsP : kmb * numbControlPtsP , 1 + (jmb-1) * numbControlPtsUx : jmb * numbControlPtsUx) = Px_mb;

                        disp(['Px - FINISHED ASSEMBLING LOOP (',num2str(kmb),' , ',num2str(jmb),')']);

                    end
                end
                
                %% Py ASSEMBLING LOOP
                
                assemblerStruct = [];

                for kmb = 1:obj.discStruct.numbModesP
                    for jmb = 1:obj.discStruct.numbModesUy
                        
                        assemblerStruct.index1     = kmb;
                        assemblerStruct.index2     = jmb;
                        assemblerStruct.wgh1       = augVerWeightsP;
                        assemblerStruct.wgh2       = augVerWeightsUy;
                        assemblerStruct.mb1        = modalBasisP;
                        assemblerStruct.mb2        = modalBasisUy;
                        assemblerStruct.dmb1       = modalBasisDerP;
                        assemblerStruct.dmb2       = modalBasisDerUy;
                        assemblerStruct.param      = Computed;
                        assemblerStruct.geodata    = geoData;
                        assemblerStruct.jacFunc    = jacFunc;
                        assemblerStruct.msh1       = mshP;
                        assemblerStruct.msh2       = mshUy;
                        assemblerStruct.space1     = spaceP;
                        assemblerStruct.space2     = spaceUy;
                        assemblerStruct.spaceFunc1 = spaceFuncP;
                        assemblerStruct.spaceFunc1 = spaceFuncUy;
                        assemblerStruct.aLift      = aLift;
                        assemblerStruct.bLift      = bLift;
                        assemblerStruct.lifting    = lifting;
                        
                        [Py_mb] = assemblerIGAScatterPy(assemblerStruct);

                        % Assignment of the Block Matrix Just Assembled

                        Py(1 + (kmb-1) * numbControlPtsP : kmb * numbControlPtsP , 1 + (jmb-1) * numbControlPtsUy : jmb * numbControlPtsUy) = Py_mb;

                        disp(['Py - FINISHED ASSEMBLING LOOP (',num2str(kmb),' , ',num2str(jmb),')']);

                    end
                end
                
                %% Qx ASSEMBLING LOOP
                
                Qx(:) = Px(:)';
                
                %% Qy ASSEMBLING LOOP
                
                Qy(:) = Py(:)';

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

%% Method 'assemblerIGAScatterAx'
            
function [Mx,Ax,Fx] = assemblerIGAScatterAx(assemblerStruct)

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
    
    funcToIntegrale = (assemblerStruct.param.forceX) .* assemblerStruct.jacFunc.evalDetJac;
    funcWeight      = assemblerStruct.mb2 .* assemblerStruct.wgh2;
    forceVec        = sum(funcToIntegrale .* funcWeight , 1);
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    % Computation of coefficients of the modal expansion
    
    funcToIntegrate_01 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi1;
    funcToIntegrate_02 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi1_Proj1;
    funcToIntegrate_03 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2;
    funcToIntegrate_04 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2_Proj1;
    funcToIntegrate_05 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2;
    funcToIntegrate_06 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2_Proj1;
    funcToIntegrate_07 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi2;
    funcToIntegrate_08 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi2_Proj1;
    funcToIntegrate_09 = jacFunc.evalDetJac .* ones(size(assemblerStruct.param.nu));
    
    % Computation of quadrature weights for the modal expansion
    
    funcWeight_01 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_02 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_03 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_04 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_05 = assemblerStruct.dmb1 .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_06 = assemblerStruct.dmb1 .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_07 = assemblerStruct.dmb1 .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_08 = assemblerStruct.dmb1 .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_09 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    
    % Numerical integration of the modal coefficients
    
    aux1   = sum(funcToIntegrate_01 .* funcWeight_01 , 1);
    aux2   = sum(funcToIntegrate_02 .* funcWeight_02 , 1);
    aux3   = sum(funcToIntegrate_03 .* funcWeight_03 , 1);
    aux4   = sum(funcToIntegrate_04 .* funcWeight_04 , 1);
    aux5   = sum(funcToIntegrate_05 .* funcWeight_05 , 1);
    aux6   = sum(funcToIntegrate_06 .* funcWeight_06 , 1);
    aux7   = sum(funcToIntegrate_07 .* funcWeight_07 , 1);
    aux8   = sum(funcToIntegrate_08 .* funcWeight_08 , 1);
    aux9   = sum(funcToIntegrate_09 .* funcWeight_09 , 1);
    
    % Definition of the modal coefficients
    
    m00 = aux9;
    r00 = aux7 + aux8;
    r10 = aux5 + aux6;
    r01 = aux3 + aux4;
    r11 = aux1 + aux2;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = assemblerStruct.msh1.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_Mass  = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_00    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_10    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_01    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_11    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    rhs         = zeros (assemblerStruct.space1.ndof, 1);
    
    for iel = 1:numbKnots
        
        m00Local = m00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_Mass  = Local_Mass + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, m00Local);
        
        Local_00    = Local_00   + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r00Local);
        Local_10    = Local_10   + op_gradu_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r10Local);
        Local_01    = Local_01   + op_u_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r01Local);
        Local_11    = Local_11   + op_gradu_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r11Local);
        
        rhs         = rhs + op_f_v (assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, fLocal);

    end

    Mx  = Local_Mass;
    Ax  = Local_11 + Local_10 + Local_01 + Local_00;
    Fx  = rhs;
    
end

%% Method 'assemblerIGAScatterAy'
            
function [My,Ay,Fy] = assemblerIGAScatterAy(assemblerStruct)

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
    
    funcToIntegrale = (assemblerStruct.param.forceY) .* assemblerStruct.jacFunc.evalDetJac;
    funcWeight      = assemblerStruct.mb2 .* assemblerStruct.wgh2;
    forceVec        = sum(funcToIntegrale .* funcWeight , 1);
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    % Computation of coefficients of the modal expansion
    
    funcToIntegrate_01 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi1;
    funcToIntegrate_02 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi1_Proj2;
    funcToIntegrate_03 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2;
    funcToIntegrate_04 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2_Proj2;
    funcToIntegrate_05 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2;
    funcToIntegrate_06 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1DotGPhi2_Proj2;
    funcToIntegrate_07 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi2;
    funcToIntegrate_08 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.modGradPhi2_Proj2;
    funcToIntegrate_09 = jacFunc.evalDetJac .* ones(size(assemblerStruct.param.nu));
    
    % Computation of quadrature weights for the modal expansion
    
    funcWeight_01 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_02 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_03 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_04 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_05 = assemblerStruct.dmb1 .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_06 = assemblerStruct.dmb1 .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_07 = assemblerStruct.dmb1 .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_08 = assemblerStruct.dmb1 .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_09 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    
    % Numerical integration of the modal coefficients
    
    aux1   = sum(funcToIntegrate_01 .* funcWeight_01 , 1);
    aux2   = sum(funcToIntegrate_02 .* funcWeight_02 , 1);
    aux3   = sum(funcToIntegrate_03 .* funcWeight_03 , 1);
    aux4   = sum(funcToIntegrate_04 .* funcWeight_04 , 1);
    aux5   = sum(funcToIntegrate_05 .* funcWeight_05 , 1);
    aux6   = sum(funcToIntegrate_06 .* funcWeight_06 , 1);
    aux7   = sum(funcToIntegrate_07 .* funcWeight_07 , 1);
    aux8   = sum(funcToIntegrate_08 .* funcWeight_08 , 1);
    aux9   = sum(funcToIntegrate_09 .* funcWeight_09 , 1);
    
    % Definition of the modal coefficients
    
    m00 = aux9;
    r00 = aux7 + aux8;
    r10 = aux5 + aux6;
    r01 = aux3 + aux4;
    r11 = aux1 + aux2;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = assemblerStruct.msh1.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_Mass  = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_00    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_10    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_01    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_11    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    rhs         = zeros (assemblerStruct.space1.ndof, 1);
    
    for iel = 1:numbKnots
        
        m00Local = m00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        fLocal   = forceVec((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_Mass  = Local_Mass + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, m00Local);
        
        Local_00    = Local_00   + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r00Local);
        Local_10    = Local_10   + op_gradu_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r10Local);
        Local_01    = Local_01   + op_u_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r01Local);
        Local_11    = Local_11   + op_gradu_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r11Local);
        
        rhs         = rhs + op_f_v (assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, fLocal);

    end

    My  = Local_Mass;
    Ay  = Local_11 + Local_10 + Local_01 + Local_00;
    Fy  = rhs;
    
end

%% Method 'assemblerIGAScatterBxy'
            
function [Bxy] = assemblerIGAScatterBxy(assemblerStruct)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    % Computation of coefficients of the modal expansion
    
    funcToIntegrate_01 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1Proj1DotGPhi1Proj2;
    funcToIntegrate_02 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi2Proj1DotGPhi1Proj2;
    funcToIntegrate_03 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1Proj1DotGPhi2Proj2;
    funcToIntegrate_04 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi2Proj1DotGPhi2Proj2;
    
    % Computation of quadrature weights for the modal expansion
    
    funcWeight_01 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_02 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_03 = assemblerStruct.dmb1 .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_04 = assemblerStruct.dmb1 .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    
    % Numerical integration of the modal coefficients
    
    aux1   = sum(funcToIntegrate_01 .* funcWeight_01 , 1);
    aux2   = sum(funcToIntegrate_02 .* funcWeight_02 , 1);
    aux3   = sum(funcToIntegrate_03 .* funcWeight_03 , 1);
    aux4   = sum(funcToIntegrate_04 .* funcWeight_04 , 1);
    
    % Definition of the modal coefficients
    
    r00 = aux4;
    r10 = aux3;
    r01 = aux2;
    r11 = aux1;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = assemblerStruct.msh1.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_10    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_01    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_11    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00    = Local_00   + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r00Local);
        Local_10    = Local_10   + op_gradu_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r10Local);
        Local_01    = Local_01   + op_u_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r01Local);
        Local_11    = Local_11   + op_gradu_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r11Local);
        
    end

    Bxy = Local_11 + Local_10 + Local_01 + Local_00;
    
end

%% Method 'assemblerIGAScatterByx'
            
function [Byx] = assemblerIGAScatterByx(assemblerStruct)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    % Computation of coefficients of the modal expansion
    
    funcToIntegrate_01 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1Proj1DotGPhi1Proj2;
    funcToIntegrate_02 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi1Proj1DotGPhi2Proj2;
    funcToIntegrate_03 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi2Proj1DotGPhi1Proj2;
    funcToIntegrate_04 = jacFunc.evalDetJac .* assemblerStruct.param.nu .* assemblerStruct.jacFunc.GPhi2Proj1DotGPhi2Proj2;
    
    % Computation of quadrature weights for the modal expansion
    
    funcWeight_01 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_02 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    funcWeight_03 = assemblerStruct.dmb1 .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_04 = assemblerStruct.dmb1 .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    
    % Numerical integration of the modal coefficients
    
    aux1   = sum(funcToIntegrate_01 .* funcWeight_01 , 1);
    aux2   = sum(funcToIntegrate_02 .* funcWeight_02 , 1);
    aux3   = sum(funcToIntegrate_03 .* funcWeight_03 , 1);
    aux4   = sum(funcToIntegrate_04 .* funcWeight_04 , 1);
    
    % Definition of the modal coefficients
    
    r00 = aux4;
    r10 = aux3;
    r01 = aux2;
    r11 = aux1;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = assemblerStruct.msh1.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_10    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_01    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_11    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r10Local = r10((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r11Local = r11((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00    = Local_00   + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r00Local);
        Local_10    = Local_10   + op_gradu_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r10Local);
        Local_01    = Local_01   + op_u_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r01Local);
        Local_11    = Local_11   + op_gradu_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r11Local);
        
    end

    Byx = Local_11 + Local_10 + Local_01 + Local_00;
    
end

%% Method 'assemblerIGAScatterPx'
            
function [Px] = assemblerIGAScatterPx(assemblerStruct)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    % Computation of coefficients of the modal expansion
    
    funcToIntegrate_01 = jacFunc.evalDetJac .* assemblerStruct.jacFunc.Phi1_dx;
    funcToIntegrate_02 = jacFunc.evalDetJac .* assemblerStruct.jacFunc.Phi2_dx;
    
    % Computation of quadrature weights for the modal expansion
    
    funcWeight_01 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_02 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    
    % Numerical integration of the modal coefficients
    
    aux1   = sum(funcToIntegrate_01 .* funcWeight_01 , 1);
    aux2   = sum(funcToIntegrate_02 .* funcWeight_02 , 1);
    
    % Definition of the modal coefficients
    
    r00 = aux2;
    r01 = aux1;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = assemblerStruct.msh1.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_01    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00    = Local_00   + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r00Local);
        Local_01    = Local_01   + op_u_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r01Local);
        
    end

    Px = Local_01 + Local_00;
    
end

%% Method 'assemblerIGAScatterPy'
            
function [Py] = assemblerIGAScatterPy(assemblerStruct)

    %% IMPORT CLASS
    
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% MODAL BASIS INTEGRATION
    %---------------------------------------------------------------------%
    % Computation of the auxiliary coefficients, presented in the original
    % work as 'r_{ik}^{st}'. Those coeffients simplify the computation of 
    % the parts to be assembled.
    %---------------------------------------------------------------------%
    
    % Computation of coefficients of the modal expansion
    
    funcToIntegrate_01 = jacFunc.evalDetJac .* assemblerStruct.jacFunc.Phi1_dy;
    funcToIntegrate_02 = jacFunc.evalDetJac .* assemblerStruct.jacFunc.Phi2_dy;
    
    % Computation of quadrature weights for the modal expansion
    
    funcWeight_01 = assemblerStruct.mb1  .* assemblerStruct.mb2  .* assemblerStruct.wgh2;
    funcWeight_02 = assemblerStruct.mb1  .* assemblerStruct.dmb2 .* assemblerStruct.wgh2;
    
    % Numerical integration of the modal coefficients
    
    aux1   = sum(funcToIntegrate_01 .* funcWeight_01 , 1);
    aux2   = sum(funcToIntegrate_02 .* funcWeight_02 , 1);
    
    % Definition of the modal coefficients
    
    r00 = aux2;
    r01 = aux1;
    
    %% ASSEMBLE OF STIFFNESS MATRIX AND RHS VECTOR
    
    numbKnots = assemblerStruct.msh1.nel_dir;
    numbHorNodes = length(r00)/numbKnots;
    
    Local_00    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    Local_01    = spalloc (assemblerStruct.space1.ndof, assemblerStruct.space2.ndof, 3 * assemblerStruct.space2.ndof);
    
    for iel = 1:numbKnots
        
        r00Local = r00((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        r01Local = r01((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
        
        Local_00    = Local_00   + op_u_v(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r00Local);
        Local_01    = Local_01   + op_u_gradv(assemblerStruct.spaceFunc1{iel,1}, assemblerStruct.spaceFunc2{iel,1}, assemblerStruct.spaceFunc1{iel,3}, r01Local);
        
    end

    Py = Local_01 + Local_00;
    
end