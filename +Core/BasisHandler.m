classdef BasisHandler
    
    %% BASIS HANDLER CLASS
    % The BasisHandler is a class that contain all the scripts responsible 
    % for the generation of the Finite Element basis, IGA basis and the 
    % modal basis used in the Hierarchical Reduction of the problem. On one
    % hand, the class properties, i.e. the properties of the objects were 
    % defined based on the variables needed to implement each function. On
    % the other hand, all of the previous functions are concentrated and
    % organized in the class methods.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% BASIS HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the BasisHandler
        % encapsulate all of the variables needed to run the methods
        % generating the base functions needed through out the solution
        % of the problem.
        
        % NEW FEM BASIS PROPERTIES
        
        degreeFiniteElement;     % Degree of the Finite Elements Used (In Our Case
                                 % We Consider It Equal to '2'
                      
        meshNodesX;              % Original Mesh
        
        meshQuadratureNodesX;    % Mesh Containing the Quadrature Nodes
        
        % NEW MODAL BASIS PROPERTIES
        
        dimModalBasis;           % Dimension of the Modal Basis
        
        evalNodesY;              % Nodes to Evaluate the Modal Basis and their
                                 % Derivatives
        
        labelUpBoundCond;        % Contains the Label Identifying the Nature of
                                 % the Boundary Conditions on the Upper Limit of
                                 % the Domain
                      
        labelDownBoundCond;      % Contains the Label Identifying the Nature of
                                 % the Boundary Conditions on the Lower Limit of
                                 % the Domain
                      
        coeffForm;               % Data Structure Containing All the @-Functions
                                 % and the Constants Relative to the Bilinear Form
                      
        % NEW IGA BASIS PROPERTIES
        
        degreePolySplineBasis;   % Degree of the Polynomial B-Spline Basis
        
        continuityParameter;     % Degree of Continuity of the Basis 'C^(p-k)'
        
        numbElements;            % Number of Elements
        
        numbQuadPointPerElem;    % Number of Quadrature Points per Element
        
        leftBoundDomainX;        % Left Bound of the Domain in the X Direction
        
        rightBoundDomainX;       % Right Bound of the Domain in the X Direction
        
        % NEW FOURIER MODAL BASIS PROPERTIES
        
        dimFourierBase;          % Number of Functions in the Base (Dirichlet). If we
                                 % consider the Robin Case, the number of functions
                                 % goes to 'm+1' 
                      
        evalFourierNodes;        % Point in which the Base Functions will be Evaluated.
                                 % They often coincide with the Quadrature Nodes
                      
        % NEW LEGENDRE MODAL BASIS PROPERTIES
        
        dimLegendreBase;         % Number of Functions in the Base (Dirichlet). If we
                                 % consider the Robin Case, the number of functions
                                 % goes to 'm+1' 
                      
        evalLegendreNodes;       % Point in which the Base Functions will be Evaluated.
                                 % They often coincide with the Quadrature Nodes
        
    end
    
    methods (Access = public)
        
        %% BASIS HANDLER - CONSTRUCT METHOD
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

        %% BASIS HANDLER - METHODS
        % The following functions correspond to the individual scripts
        % previously developed of the functions 'newFiniteElementBasis',
        % 'newModalBasis', 'newModalBasisFourier', 'newIsoGeoBasis'
        % and 'newModalBasisLegendre'.
        
            %% Method 'newFiniteElementBasis'
            
            function [shapeFunction,shapeFunctionDerivative] = newFiniteElementBasis(obj)

                %%
                % newFiniteElementBasis - This function evaluates the standard finite element
                %                         basis and their at the nodes of the mesh.
                %
                % The inputs are:
                %%
                %   (1)  degreeFiniteElement    : Degree of the Finite Elements Used (In Our Case
                %                                 We Consider It Equal to '2'
                %   (2)  meshNodesX             : Original Mesh
                %   (3)  meshQuadratureNodesX   : Mesh Containing the Quadrature Nodes
                %   
                % The outputs are:
                %%
                %   (1) shapeFunction           : Matrix Containing the Shape Functions
                %                                 Evaluated in the Nodes Defined in 'meshQuadratureNodesX' 
                %   (2) shapeFunctionDerivative : Matrix Conatining the Derivative of the Shape
                %                                 Functions Evaluated in the Nodes Defined in
                %                                 'meshQuadratureNodesX'

                % Number of Intervals

                  numbIntervals = length(obj.meshNodesX) - 1;

                % Number of Quadrature Nodes in Each Interval

                  numbQuadNodesPerInterval = length(obj.meshQuadratureNodesX) / numbIntervals;

                % Basis Function of the Elements

                  shapeFunction = zeros(length(obj.meshQuadratureNodesX),obj.degreeFiniteElement + 1);

                % Coordinates of the Elements

                  shapeFunctionDerivative = zeros(length(obj.meshQuadratureNodesX),obj.degreeFiniteElement + 1);

                  for ie = 1:numbIntervals

                    xl = obj.meshNodesX(ie);             % Left Bound of the Interval
                    xr = obj.meshNodesX(ie+1);           % Right Bound of the Interval

                    for iqn = 1:numbQuadNodesPerInterval

                      i_x = (ie-1) * numbQuadNodesPerInterval + iqn;

                      % Extraction of the Integration Nodes for the Current Interval

                      x = obj.meshQuadratureNodesX(i_x);
                      
                      %---------------------------------------------------%
                      % Note: The function fem_bas is in the folder Util
                      % and computes the coefficients of the finite
                      % element basis based on the degree, current
                      % integration node and the left and right limits of
                      % the current element. 
                      % This means that the resulting matrices 'shape' and
                      % 'shapeDer' contains the finite element basis
                      % evaluated in all discretization and integration
                      % nodes.
                      %---------------------------------------------------%

                      [shapeFunction(i_x,:),shapeFunctionDerivative(i_x,:)] = fem_bas(obj.degreeFiniteElement,x,xl,xr);

                    end
                  end
                return
            end
            
            %% Method 'newIsoGeoBasis'
            
            function [coeffIsoGeoBase,coeffIsoGeoBaseDer,controlPts] = newIsoGeoBasis(obj)

                %%
                % newIsoGeoBasis  - This function evaluates the
                % isogeometric basis in the desired points.
                %
                % The inputs are:
                %%
                %   (1)  degreePolySplineBasis  : Degree Used in for the Polynomial Basis
                %   (2)  continuityParameter    : Parameter Defining the the Continuity Properties
                %                                 among the Elements of the Mesh
                %   (3)  numbElements           : Number of Elements
                %   (4)  numbQuadPointPerElem   : Number of Quadrature Points per Element
                %   (5)  leftBoundDomainX       : Left Bound of the Domain in the X Direction
                %   (6)  rightBoundDomainX      : Right Bound of the Domain in the X Direction
                %   
                % The outputs are:
                %%
                %   (1) coeffIsoGeoBase         : Matrix Containing the Base Functions
                %                                 Evaluated at the Nodes
                %   (2) coeffIsoGeoBaseDer      : Matrix Conatining the Derivative of the Base
                %                                 Functions Evaluated at the Nodes

                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                % Number of Integration Points

                ngauss = obj.numbQuadPointPerElem;

                % Construction of the Mesh to Generate the Base Functions

                %---------------------------------------------------------------------%
                % Note:
                % We constructed the mesh considering a linear parametrization
                % procedure, with continuity C^(p-k).
                %---------------------------------------------------------------------%

                knot      = augknt([obj.leftBoundDomainX obj.rightBoundDomainX],obj.degreePolySplineBasis + 1);
                
                h         = 1/obj.numbElements;
                
                internKnot = linspace(0,1,obj.numbElements+1)';
                internKnot = internKnot(2:end-1);
                
                ins       = sort(reshape(internKnot*ones(1,obj.continuityParameter),1,[]));      %,obj.continuityParameter*(obj.numbElements - 1)));
                
                ins       = (obj.rightBoundDomainX - obj.leftBoundDomainX) * ins + obj.leftBoundDomainX;
                
                [controlPts,knot] = bspkntins(obj.degreePolySplineBasis, ...
                            obj.leftBoundDomainX:(obj.rightBoundDomainX - obj.leftBoundDomainX)*1/obj.degreePolySplineBasis:obj.rightBoundDomainX,knot,ins);
                
                numbControlPts   = length(controlPts);
                
                nspan = numbControlPts - obj.degreePolySplineBasis;
                
                % Definition of the Nodes and Weights inside the Reference Interval
                
                objGaussLeg = IntegrateHandler();
                objGaussLeg.numbQuadNodes = ngauss;
                [~,gp,~] = gaussLegendre(objGaussLeg); 
                
                % [gp,~] = gauss(ngauss);

                %---------------------------------------------------------------------%
                % Note:
                % The 'gauss' function is responsable for generating the nodes and
                % respective weights for a given interval from the number of desired
                % integration points.
                %---------------------------------------------------------------------%

                gausspt = zeros(1,ngauss * obj.numbElements);

                n = 0;

                for n0 = obj.degreePolySplineBasis + 1:obj.degreePolySplineBasis + nspan

                    if knot(n0) ~= knot(n0+1)

                        gausspt(n*ngauss+1:(n+1)*ngauss) =...
                            ((knot(n0+1) - knot(n0))*gp + knot(n0+1) + knot(n0))/2;
                        n = n + 1;

                    end

                end
                
                coeffIsoGeoBase  = spcol(knot,obj.degreePolySplineBasis + 1,sort([gausspt,gausspt]));
                
                % Derivative of the Base Functions Evaluated

                coeffIsoGeoBaseDer = coeffIsoGeoBase(2:2:end,:);

                % Base Functions Evaluated

                coeffIsoGeoBase  = coeffIsoGeoBase(1:2:end,:); 

                %---------------------------------------------------------------------%
                % Note:
                % One possible idea to improve the code is to build Jacobians for each
                % element already inside this script.
                %---------------------------------------------------------------------%

                return 

            end
            
            %% Method 'newModalBasis'
            
            function [coeffModalBase,coeffModalBaseDer] = newModalBasis(obj)

                %%
                % newModalBasis   - This function evaluates the modal basis and their
                %                   derivatives at the nodes 'yq', solving the problem
                %                   of SL.
                %
                % The inputs are:
                %%
                %   (1)  dimModalBasis         : Dimension of the Modal Basis
                %   (2)  evalNodesY            : Nodes to Evaluate the Modal Basis and their
                %                                Derivatives
                %   (3)  labelUpBoundCond      : Contains the Label Identifying the Nature of
                %                                the Boundary Conditions on the Upper Limit of
                %                                the Domain
                %   (4)  labelDownBoundCond    : Contains the Label Identifying the Nature of
                %                                the Boundary Conditions on the Lower Limit of
                %                                the Domain
                %   (5)  coeffForm             : Data Strusture Containing All the @-Functions
                %                                and the Constants Relative to the Bilinear Form
                %
                % The outputs are:
                %%
                %   (1) coeffModalBase         : Coefficients of the Modal Basis at the Nodes
                %   (2) coeffModalBaseDer      : Coefficients of the Derivative of the Modal
                %                                Basis at the Nodes

                %% IMPORT CLASSES
            
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
            
                %% Evaluation of the Problem Coefficients
                
                if(nargin(obj.coeffForm.mu) == 3) muval = obj.coeffForm.mu(1,1,1);
                else muval = obj.coeffForm.mu(1,1);
                end
                
                sigma = obj.coeffForm.coeffrobin;

                %% Initialization of the Vector Containing the Cefficients in the New
                % Basis

                coeffModalBase   = zeros( length(obj.evalNodesY), obj.dimModalBasis);
                coeffModalBaseDer = zeros( length(obj.evalNodesY), obj.dimModalBasis);

                % Initialization of the Vector Containing the Eigenvalues of the New
                % Basis

                lambda = zeros(obj.dimModalBasis,1);

                %% Computation of the Eigenvalues of the New Basis

                %---------------------------------------------------------------------%
                % Note:
                % The eigenvalues of the new basis are computed considering the
                % different possible cases fot the boundary conditions in the Upper 
                % and Lower limit of the domain.
                %---------------------------------------------------------------------%    

                switch [obj.labelUpBoundCond,obj.labelDownBoundCond]

                    case 'dirdir'   % Dirichlet Dirichlet

                        for i = 1:obj.dimModalBasis
                            lambda(i) = i * pi;
                        end

                        B = @(lambda) 0;

                    case 'dirrob'   % Dirichlet Robin

                        stlv  = @(lmb)  tan(lmb) + muval*lmb/sigma;
                        B     = @(lmb) -tan(lmb);
                        
                        obj_computeEigenvalues_1 = IntegrateHandler();
                        
                        obj_computeEigenvalues_1.funcStourmLiou = stlv;
                        obj_computeEigenvalues_1.numbEigenvalues = obj.dimModalBasis;
                        obj_computeEigenvalues_1.labelUpBoundCond = obj.labelUpBoundCond;
                        obj_computeEigenvalues_1.labelDownBoundCond = obj.labelDownBoundCond;
                        obj_computeEigenvalues_1.coeffForm = obj.coeffForm;

                        lambda = computeEigenvalues(obj_computeEigenvalues_1);

                    case 'robrob'   % Robin Robin

                        stlv  = @(lmb) 2*muval*lmb + tan(lmb).*(sigma - muval.^2*lmb.^2/sigma);
                        B     = @(lmb) muval*lmb/sigma;

                        obj_computeEigenvalues_2 = IntegrateHandler();
                        
                        obj_computeEigenvalues_2.funcStourmLiou = stlv;
                        obj_computeEigenvalues_2.numbEigenvalues = obj.dimModalBasis;
                        obj_computeEigenvalues_2.labelUpBoundCond = obj.labelUpBoundCond;
                        obj_computeEigenvalues_2.labelDownBoundCond = obj.labelDownBoundCond;
                        obj_computeEigenvalues_2.coeffForm = obj.coeffForm;

                        lambda = computeEigenvalues(obj_computeEigenvalues_2);
                        
                    case 'robdir'   % Robin Dirichlet

                        stlv  = @(lmb) muval*lmb/sigma + tan(lmb);
                        B     = @(lmb) 0;

                        obj_computeEigenvalues_3 = IntegrateHandler();
                        
                        obj_computeEigenvalues_3.funcStourmLiou = stlv;
                        obj_computeEigenvalues_3.numbEigenvalues = obj.dimModalBasis;
                        obj_computeEigenvalues_3.labelUpBoundCond = obj.labelUpBoundCond;
                        obj_computeEigenvalues_3.labelDownBoundCond = obj.labelDownBoundCond;
                        obj_computeEigenvalues_3.coeffForm = obj.coeffForm;

                        lambda = computeEigenvalues(obj_computeEigenvalues_3);

                    case 'neuneu'

                        for i = 1:obj.dimModalBasis
                            lambda(i) = (i-1) * pi;
                        end

                        B = @(lmb) 1;

                    otherwise
                        disp('In newModalBasis: Boundary Conditions Not Recognized or Not Yet Available')
                        
                end
                
                % Evaluation of the Coefficient in the New Basis

                %---------------------------------------------------------------------%
                % Note:
                % The coeffients at the new basis are also computed considering the
                % different possible cases fot the boundary conditions in the Upper 
                % and Lower limit of the domain.
                %---------------------------------------------------------------------%

                switch [obj.labelUpBoundCond,obj.labelDownBoundCond]

                    case {'dirdir','dirrob','robrob','robdir'}

                        for n = 1 : obj.dimModalBasis
                            b         = B( lambda(n) );
                            L2norm    = sqrt( ( 1 + b^2 )/2 + ( b^2 - 1 )*sin( 2*lambda(n) )/( 4*lambda(n) ) + b*( sin( lambda(n) ) )^2/lambda(n) );

                            % Normalization for the Coefficient A

                            for i = 1 : length(obj.evalNodesY)
                                coeffModalBase( i, n) = sin( lambda(n) * obj.evalNodesY(i) )/L2norm + b * cos( lambda(n)*obj.evalNodesY(i) )/L2norm;
                                coeffModalBaseDer( i, n) = lambda(n)*cos( lambda(n)*obj.evalNodesY(i) )/L2norm - b * lambda(n) * sin( lambda(n)*obj.evalNodesY(i) )/L2norm;
                            end

                        end

                    case 'neuneu'

                        coeffModalBase  ( :, 1) = ones ( length(obj.evalNodesY), 1 );
                        coeffModalBaseDer( :, 1) = zeros( length(obj.evalNodesY), 1 );

                        for n = 2 : obj.dimModalBasis
                            for i = 1 : length(obj.evalNodesY)
                                coeffModalBase  ( i, n) =             sqrt(2)*cos( lambda(n)*obj.evalNodesY(i) );
                                coeffModalBaseDer( i, n) = - lambda(n)*sqrt(2)*sin( lambda(n)*obj.evalNodesY(i) );
                            end
                        end
                end
                
            end            
                    
            %% Method 'newModalBasisFourier'
            
            function [coeffModalBaseFourier,coeffModalBaseFourierDer] = newModalBasisFourier(obj)

                %%
                % newModalBasisFourier    - This function evaluates the modal basis and
                %                           their derivatives at the nodes. In this
                %                           case the we consider the functions of the
                %                           modal basis as a fourier decomposition.
                %
                % The inputs are:
                %%
                %   (1)  dimFourierBase       : Number of Functions in the Base (Dirichlet). If we
                %                               consider the Robin Case, the number of functions
                %                               goes to 'm+1' 
                %   (2)  evalFourierNodes     : Point in which the Base Functions will be Evaluated.
                %                               They often coincide with the Quadrature Nodes
                %   
                % The outputs are:
                %%
                %   (1) coeffModalBaseFourier       : Matrix Containing the Base Functions (Columns)
                %                                     Evaluated at the Nodes (Rows)
                %   (2) coeffModalBaseFourierDer    : Matrix Conatining the Derivative of the Base
                %                                     Functions (Columns) Evaluated at the Nodes (Rows)

                % Initialization of the Vector Containing the Cefficients in the New
                % Basis

                 coeffModalBaseFourier   = zeros( length(obj.evalFourierNodes), obj.dimFourierBase);
                 coeffModalBaseFourierDer = zeros( length(obj.evalFourierNodes), obj.dimFourierBase);

                 for i = 1:length(obj.evalFourierNodes)
                     coeffModalBaseFourier  ( i , 1 ) = 1;
                     coeffModalBaseFourierDer( i , 1 ) = 0;
                 end

                 % Loop on the Base Functions

                 for n = 2:2:obj.dimFourierBase

                     % Loop on the Section

                     for i = 1:length(obj.evalFourierNodes)
                         coeffModalBaseFourier  ( i, n ) = sqrt(2) * sin(n * pi * obj.evalFourierNodes(i));
                         coeffModalBaseFourierDer( i, n ) = n * pi * sqrt(2) * cos(n * pi * obj.evalFourierNodes(i));
                     end

                 end

                 for n = 3:2:obj.dimFourierBase

                     % Loop on the Section

                     for i = 1:length(obj.evalFourierNodes)
                         coeffModalBaseFourier  ( i, n ) = sqrt(2)*cos((n-1) * pi * obj.evalFourierNodes(i));
                         coeffModalBaseFourierDer( i, n ) = -(n-1) * pi * sqrt(2) * sin((n-1) * pi * obj.evalFourierNodes(i));
                     end

                 end
            end
            
            %% Method 'newModalBasisLegendre'
            
            function [coeffModalBaseLegendre,coeffModalBaseLegendreDer] = newModalBasisLegendre(obj)

                %%
                % newModalBasisLegendre    - This function evaluates the modal basis 
                %                            and their derivatives at the nodes. In 
                %                            this case the we consider the functions 
                %                            of the modal basis Standard Legendre.
                %
                % The inputs are:
                %%
                %   (1)  dimLegendreBase    : Number of Functions in the Base (Dirichlet). If we
                %                             consider the Robin Case, the number of functions
                %                             goes to 'm+1' 
                %   (2)  evalLegendreNodes  : Point in which the Base Functions will be Evaluated.
                %                             They often coincide with the Quadrature Nodes
                %   
                % The outputs are:
                %%
                %   (1) coeffModalBaseLegendre    : Matrix Containing the Base Functions (Columns)
                %                                   Evaluated at the Nodes (Rows)
                %   (2) coeffModalBaseLegendreDer : Matrix Conatining the Derivative of the Base
                %                                   Functions (Columns)
                %                                   Evaluated at the Nodes (Rows) 

                if(strcmp(obj.labelUpBoundCond,'rob') && strcmp(obj.labelDownBoundCond,'rob'))

                    coeffModalBaseLegendre = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);
                    coeffModalBaseLegendreDer = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);

                    obj_polyLegendre = IntegrateHandler();
                    obj_polyLegendre.degreePolyLegendre = obj.dimLegendreBase;
                    [P,Pd] = polyLegendre(obj_polyLegendre);

                    % Loop on the Base Functions

                    for n = 1:obj.dimLegendreBase

                        % Loop on the Section

                        coeffModalBaseLegendre(:,n) = polyval(P{n},obj.evalLegendreNodes);
                        coeffModalBaseLegendreDer(:,n) = polyval(Pd{n},obj.evalLegendreNodes);

                    end

                elseif(strcmp(obj.labelUpBoundCond,'dir') && strcmp(obj.labelDownBoundCond,'dir'))

                    coeffModalBaseLegendre   = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);
                    coeffModalBaseLegendreDer = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);

                    % Loop on the Base Functions

                    for n = 1:obj.dimLegendreBase

                        % Loop on the Section

                        for i = 1:length(obj.evalLegendreNodes)
                            coeffModalBaseLegendre  ( i, n ) = obj.evalLegendreNodes(i).^(n-1)*(1 - obj.evalLegendreNodes(i).^2);
                            coeffModalBaseLegendreDer( i, n ) = (n-1) * obj.evalLegendreNodes(i).^(n-2) * (1 - obj.evalLegendreNodes(i).^2) ...
                                            + obj.evalLegendreNodes(i).^(n-1) * (-2 * obj.evalLegendreNodes(i));
                        end

                    end

                end
            end
            
            %% Method 'newModalBasis3D'
            
            function [coeffs,coeffsDer1,coeffsDer2] = newModalBasis3D(obj)

                %%
                % newModalBasis   - This function evaluates the modal basis and their
                %                   derivatives at the nodes 'yq', solving the problem
                %                   of SL.
                %
                % The inputs are:
                %%
                %   (1)  dimModalBasis         : Dimension of the Modal Basis
                %   (2)  evalNodesY            : Nodes to Evaluate the Modal Basis and their
                %                                Derivatives
                %   (3)  labelUpBoundCond      : Contains the Label Identifying the Nature of
                %                                the Boundary Conditions on the Upper Limit of
                %                                the Domain
                %   (4)  labelDownBoundCond    : Contains the Label Identifying the Nature of
                %                                the Boundary Conditions on the Lower Limit of
                %                                the Domain
                %   (5)  coeffForm             : Data Strusture Containing All the @-Functions
                %                                and the Constants Relative to the Bilinear Form
                %
                % The outputs are:
                %%
                %   (1) coeffModalBase         : Coefficients of the Modal Basis at the Nodes
                %   (2) coeffModalBaseDer      : Coefficients of the Derivative of the Modal
                %                                Basis at the Nodes

                %% IMPORT CLASSES
            
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
            
                %% Evaluation of the Problem Coefficients
                
                if(nargin(obj.coeffForm.mu) == 4) muval = obj.coeffForm.mu(0,0,0,0);
                elseif(nargin(obj.coeffForm.mu) == 3) muval = obj.coeffForm.mu(0,0,0);
                else muval = obj.coeffForm.mu(0,0);
                end
                
                sigma = obj.coeffForm.coeffrobin;

                %% Initialization of the Vector Containing the Cefficients in the New
                % Basis

                coeffModalBase   = zeros( length(obj.evalNodesY), obj.dimModalBasis);
                coeffModalBaseDer = zeros( length(obj.evalNodesY), obj.dimModalBasis);

                % Initialization of the Vector Containing the Eigenvalues of the New
                % Basis

                lambda = zeros(obj.dimModalBasis,1);

                %% Computation of the Eigenvalues of the New Basis

                %---------------------------------------------------------------------%
                % Note:
                % The eigenvalues of the new basis are computed considering the
                % different possible cases fot the boundary conditions in the Upper 
                % and Lower limit of the domain.
                %---------------------------------------------------------------------%    

                switch [obj.labelUpBoundCond,obj.labelDownBoundCond]

                    case 'dirdir'   % Dirichlet Dirichlet

                        for i = 1:obj.dimModalBasis
                            lambda(i) = i * pi;
                        end

                        B = @(lambda) 0;

                    case 'dirrob'   % Dirichlet Robin

                        stlv  = @(lambda)  tan(lambda) + muval*lambda/sigma;
                        B     = @(lambda) -tan(lambda);
                        
                        obj_computeEigenvalues_1 = IntegrateHandler();
                        
                        obj_computeEigenvalues_1.funcStourmLiou = stlv;
                        obj_computeEigenvalues_1.numbEigenvalues = obj.dimModalBasis;
                        obj_computeEigenvalues_1.labelUpBoundCond = obj.labelUpBoundCond;
                        obj_computeEigenvalues_1.labelDownBoundCond = obj.labelDownBoundCond;
                        obj_computeEigenvalues_1.coeffForm = obj.coeffForm;

                        lambda = computeEigenvalues(obj_computeEigenvalues_1);

                    case 'robrob'   % Robin Robin

                        stlv  = @(lambda) 2*muval*lambda + tan(lambda).*(sigma - muval.^2*lambda.^2/sigma);
                        B     = @(lambda) muval*lambda/sigma;

                        obj_computeEigenvalues_2 = IntegrateHandler();
                        
                        obj_computeEigenvalues_2.funcStourmLiou = stlv;
                        obj_computeEigenvalues_2.numbEigenvalues = obj.dimModalBasis;
                        obj_computeEigenvalues_2.labelUpBoundCond = obj.labelUpBoundCond;
                        obj_computeEigenvalues_2.labelDownBoundCond = obj.labelDownBoundCond;
                        obj_computeEigenvalues_2.coeffForm = obj.coeffForm;

                        lambda = computeEigenvalues(obj_computeEigenvalues_2);
                        
                    case 'robdir'   % Robin Dirichlet

                        stlv  = @(lambda) muval*lambda/sigma + tan(lambda);
                        B     = @(lambda) 0;

                        obj_computeEigenvalues_3 = IntegrateHandler();
                        
                        obj_computeEigenvalues_3.funcStourmLiou = stlv;
                        obj_computeEigenvalues_3.numbEigenvalues = obj.dimModalBasis;
                        obj_computeEigenvalues_3.labelUpBoundCond = obj.labelUpBoundCond;
                        obj_computeEigenvalues_3.labelDownBoundCond = obj.labelDownBoundCond;
                        obj_computeEigenvalues_3.coeffForm = obj.coeffForm;

                        lambda = computeEigenvalues(obj_computeEigenvalues_3);

                    case 'neuneu'

                        for i = 1:obj.dimModalBasis
                            lambda(i) = (i-1) * pi;
                        end

                        B = @(lambda) 1;

                    otherwise
                        disp('In newModalBasis: Boundary Conditions Not Recognized or Not Yet Available')
                end

                % Evaluation of the Coefficient in the New Basis

                %---------------------------------------------------------------------%
                % Note:
                % The coeffients at the new basis are also computed considering the
                % different possible cases fot the boundary conditions in the Upper 
                % and Lower limit of the domain.
                %---------------------------------------------------------------------%

                switch [obj.labelUpBoundCond,obj.labelDownBoundCond]

                    case {'dirdir','dirrob','robrob','robdir'}

                        for n = 1 : obj.dimModalBasis
                            b         = B( lambda(n) );
                            L2norm    = sqrt( ( 1 + b^2 )/2 + ( b^2 - 1 )*sin( 2*lambda(n) )/( 4*lambda(n) ) + b*( sin( lambda(n) ) )^2/lambda(n) );

                            % Normalization for the Coefficient A

                            for i = 1 : length(obj.evalNodesY)
                                coeffModalBase( i, n) = sin( lambda(n) * obj.evalNodesY(i) )/L2norm + b * cos( lambda(n)*obj.evalNodesY(i) )/L2norm;
                                coeffModalBaseDer( i, n) = lambda(n)*cos( lambda(n)*obj.evalNodesY(i) )/L2norm - b * lambda(n) * sin( lambda(n)*obj.evalNodesY(i) )/L2norm;
                            end

                        end

                    case 'neuneu'

                        coeffModalBase  ( :, 1) = ones ( length(obj.evalNodesY), 1 );
                        coeffModalBaseDer( :, 1) = zeros( length(obj.evalNodesY), 1 );

                        for n = 2 : obj.dimModalBasis
                            for i = 1 : length(obj.evalNodesY)
                                coeffModalBase  ( i, n) =             sqrt(2)*cos( lambda(n)*obj.evalNodesY(i) );
                                coeffModalBaseDer( i, n) = - lambda(n)*sqrt(2)*sin( lambda(n)*obj.evalNodesY(i) );
                            end
                        end
                end
                
                N = length(obj.evalNodesY);
                m = obj.dimModalBasis;
                
                coeffs = zeros(N,N,m^2);
                coeffsDer1 = zeros(N,N,m^2);
                coeffsDer2 = zeros(N,N,m^2);
                for kk = 1:m
                    for ll = 1:m
                        for ii = 1:N
                            for jj = 1:N
                                aux1(ii,jj,ll) = coeffModalBase(ii,kk) * coeffModalBase(jj,ll);
                                aux2(ii,jj,ll) = coeffModalBaseDer(ii,kk) * coeffModalBase(jj,ll);
                                aux3(ii,jj,ll) = coeffModalBase(ii,kk) * coeffModalBaseDer(jj,ll);
                            end
                        end
                    end
                    coeffs(:,:,(kk-1)*m+1:kk*m) = aux1;
                    coeffsDer1(:,:,(kk-1)*m+1:kk*m) = aux2;
                    coeffsDer2(:,:,(kk-1)*m+1:kk*m) = aux3;
                end
            end
            
            %% Method 'newModalBasisLegendre3D'
            
            function [coeffs,coeffsDer1,coeffsDer2] = newModalBasisLegendre3D(obj)

                %%
                % newModalBasisLegendre    - This function evaluates the modal basis 
                %                            and their derivatives at the nodes. In 
                %                            this case the we consider the functions 
                %                            of the modal basis Standard Legendre.
                %
                % The inputs are:
                %%
                %   (1)  dimLegendreBase    : Number of Functions in the Base (Dirichlet). If we
                %                             consider the Robin Case, the number of functions
                %                             goes to 'm+1' 
                %   (2)  evalLegendreNodes  : Point in which the Base Functions will be Evaluated.
                %                             They often coincide with the Quadrature Nodes
                %   
                % The outputs are:
                %%
                %   (1) coeffModalBaseLegendre    : Matrix Containing the Base Functions (Columns)
                %                                   Evaluated at the Nodes (Rows)
                %   (2) coeffModalBaseLegendreDer : Matrix Conatining the Derivative of the Base
                %                                   Functions (Columns)
                %                                   Evaluated at the Nodes (Rows) 

                evalNodes = obj.evalLegendreNodes * 2 - 1;
                obj.evalLegendreNodes = evalNodes;
                
%                 if(strcmp(obj.labelUpBoundCond,'rob') && strcmp(obj.labelDownBoundCond,'rob'))
% 
%                     coeffModalBaseLegendre = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);
%                     coeffModalBaseLegendreDer = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);
% 
%                     obj_polyLegendre = IntegrateHandler();
%                     obj_polyLegendre.degreePolyLegendre = obj.dimLegendreBase;
%                     [P,Pd] = polyLegendre(obj_polyLegendre);
% 
%                     % Loop on the Base Functions
% 
%                     for m = 1:obj.dimLegendreBase
% 
%                         % Loop on the Section
% 
%                         coeffModalBaseLegendre(:,m) = polyval(P{m},obj.evalLegendreNodes);
%                         coeffModalBaseLegendreDer(:,m) = polyval(Pd{m},obj.evalLegendreNodes);
% 
%                     end
% 
%                 elseif(strcmp(obj.labelUpBoundCond,'dir') && strcmp(obj.labelDownBoundCond,'dir'))
% 
%                     coeffModalBaseLegendre   = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);
%                     coeffModalBaseLegendreDer = zeros( length(obj.evalLegendreNodes), obj.dimLegendreBase);
% 
%                     % Loop on the Base Functions
%                     
%                     for m = 1:obj.dimLegendreBase
% 
%                         % Loop on the Section
%                         
%                         syms f(x)
%                         f(x) = legendreP(m,x);
%                         df = diff(f,x);
%                         
%                         coeffModalBaseLegendre(:,m) = f(obj.evalLegendreNodes);
%                         coeffModalBaseLegendreDer(:,m) = df(obj.evalLegendreNodes);
% 
%                     end
%                 end

%                 syms f(y,z,modBasisJ,modBasisK)
%                 f(y,z,modBasisJ,modBasisK) = legendreP(modBasisJ,y) .* legendreP(modBasisK,z);
%                 dfy = diff(f,y);
%                 dfz = diff(f,z);
                
                [Y,Z] = meshgrid(evalNodes,evalNodes);
                
                N = length(evalNodes);
                m = obj.dimLegendreBase;
                
                coeffs = zeros(N,N,m^2);
                coeffsDer1 = zeros(N,N,m^2);
                coeffsDer2 = zeros(N,N,m^2);
                
                for kk = 1:m
                    for ll = 1:m
                        
                        syms f(y,z)
                        f(y,z) = legendreP(kk,y) .* legendreP(ll,z);
                        dfy = diff(f,y);
                        dfz = diff(f,z);
                        
                        aux1(:,:,ll) = f(Y,Z);
                        aux2(:,:,ll) = dfy(Y,Z);
                        aux3(:,:,ll) = dfz(Y,Z);
                    end
                    coeffs(:,:,(kk-1)*m+1:kk*m) = aux1;
                    coeffsDer1(:,:,(kk-1)*m+1:kk*m) = aux2;
                    coeffsDer2(:,:,(kk-1)*m+1:kk*m) = aux3;
                end
                
                
                
            end
            
            %% Method 'newModalBasisCheb3D'
            
            function [coeffs,coeffsDer1,coeffsDer2] = newModalBasisCheb3D(obj)
                
                % SET EIGENVALUE PRECISION
                
                numbEig = 36; 
                
                % COMPUTE CHEBCHEV DIFFERENTIATION MATRIX
                
                [chebMat,chebPts] = cheb(numbEig); 
                D2 = chebMat^2; 
                D2 = D2(2:numbEig,2:numbEig);
                
                % COMPUTE EIGENFUNCTIONS OF CHEBCHEV DIFFERENTIATION MATRIX
                
                [V,Lam] = eig(D2); 
                lam = diag(Lam);
                
                % SORT EIGENFUNCTIONS AND EIGENVALUES
                
                [~,ii] = sort(-lam);
                V = V(:,ii); 
                clf
                
                % EVALUATE THE EIGENFUNCTIONS IN THE INTEGRATION NODES
                
                evalNodes = obj.evalNodesY * 2 - 1;
                evalU = zeros(length(evalNodes),obj.dimModalBasis);
                
                figure
                for ii = 1:obj.dimModalBasis
                    
                    u = [0;V(:,ii);0];
                    fitU = polyfit(chebPts,u,numbEig);
                    evalU(:,ii) = polyval(fitU,evalNodes);
                    
                    % subplot(obj.dimModalBasis,1,ii)
                    % plot(chebPts,u,'.','markersize',12)
                    % grid on
                    % line(evalNodes,evalU(:,ii),'linewidth',.7), axis off

                end
                
                for ii = 1:length(evalNodes)
                    for jj = 1:length(evalNodes)
                        for kk = 1:obj.dimModalBasis
                            
                            coeffs(ii,jj,kk) = evalU(ii,kk) * evalU(jj,kk);
                            coeffsDer1(ii,jj,kk) = evalU(ii,kk) * evalU(jj,kk);
                            coeffsDer2(ii,jj,kk) = evalU(ii,kk) * evalU(jj,kk);
                            
                        end
                    end
                end
                
            end
    end
end