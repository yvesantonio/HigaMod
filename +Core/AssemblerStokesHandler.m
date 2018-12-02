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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MODAL BASIS FOR THE PRESSURE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_newModalBasis = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            % obj_newModalBasis.dimLegendreBase = obj.dimModalBasisP;
            % obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasisP;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomainP;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomainP;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasisP, modalBasisDerP] = newModalBasis(obj_newModalBasis);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MODAL BASIS FOR THE VELOCITY IN X %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_newModalBasis = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            % obj_newModalBasis.dimLegendreBase = obj.dimModalBasisP;
            % obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasisUx;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomainUx;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomainUx;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasisUx, modalBasisDerUx] = newModalBasis(obj_newModalBasis);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MODAL BASIS FOR THE VELOCITY IN Y %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj_newModalBasis = BasisHandler();
            
            % Set variables to use a Legendre Modal Base
            
            % obj_newModalBasis.dimLegendreBase = obj.dimModalBasisP;
            % obj_newModalBasis.evalLegendreNodes = verGLNodes;
            
            % Set variables to use a Educated Modal Basis
            
            obj_newModalBasis.dimModalBasis = obj.dimModalBasisUy;
            obj_newModalBasis.evalNodesY = verGLNodes;
            obj_newModalBasis.labelUpBoundCond = obj.label_upBoundDomainUy;
            obj_newModalBasis.labelDownBoundCond = obj.label_downBoundDomainUy;
            obj_newModalBasis.coeffForm = obj.coefficientForm;

            [modalBasisUy, modalBasisDerUy] = newModalBasis(obj_newModalBasis);
            
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