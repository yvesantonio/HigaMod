classdef BoundaryConditionHandler
    
    %% BOUNDARY CONDITION HANDLER CLASS
    % The BoundaryConditionHandler is a class that contains all of the
    % scripts responsible for the definition and handling of the properties
    % of the boundary conditions in the domain of the differential problem.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% BOUNDARY CONDITION HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the
        % BoundaryConditionHandler encapsulate all of the variables needed
        % to run the methods associated with the boundary conditons of the
        % problem. Differently from the AssemblerADRHandler class, here the
        % methods are completily independdent from each other and be
        % directly accessed.
        
            %% Offset Adjustment Properties

            labelUpBoundCond;   % Contains the Label Identifying the Nature of
                                % the Boundary Conditions on the Upper Limit of
                                % the Domain

            labelDownBoundCond; % Contains the Label Identifying the Nature of
                                % the Boundary Conditions on the Lower Limit of
                                % the Domain

            dataUpBoundCond;    % Contains the Values of the Boundary Conditions
                                % on the Upper Limir of the Domain

            dataDownBoundCond;  % Contains the Values of the Boundary Conditions
                                % on the Lower Limir of the Domain

            coeffForm;          % Data Strusture Containing All the @-Functions
                                % and the Constants Relative to the Bilinear Form

            %% Boundary Imposition Properties

            dimModalBasis;         % Dimension of the Modal Basis

            leftLimitDomainX;      % Left Limit of the Domain in the X Direction 

            rightLimitDomainX;     % Right Limit of the Domain in the X Direction

            stepMeshX;             % Vector Containing the Step of the Finite
                                   % Element Mesh

            leftLabelBoundCond;    % Boundary Conditions on the Left Bound of the
                                   % Domain

            rightLabelBoundCond    % Boundary Conditions on the Right Bound of the
                                   % Domain

            leftDomainBoundary;    % Left Boundary of the Domain

            rightDomainBoundary;   % Right Boundary of the Domain

            blockMatA;             % Bloc Matrix Defining the Problem

            blockVecB;             % Bloc Vector Defining the Problem   

            degreePolySplineBasis; % Degree of the Polynomial B-Spline Basis

            continuityParameter;   % Degree of Continuity of the Basis 'C^(p-k)'
            
            boundaryStruct;        % Structure containing all the information concerning
                                   % the boundary conditions of the
                                   % problems

            %% Inflow Interpolation Properties

            liftAdjustment      % Offset Adjustment to the Boundary Conditions

            modalBasis;         % Modal Basis Used to Represent the Vertical
                                % Direction

            nodeWeightInflow;   % Weight of the Node Used in the Current
                                % Computation

            dataDirBoundCond;   % Data Structure Containing the Dirichlet
                                % Boundary Conditions

            %% Interface Interpolation Properties

            boundModalBasis;        % Modal Basis Facing the Boundary

            inputLiftAdjust;        % Offset Adjustment Coefficient at the Input

            outputLiftAdjust;       % Offset Adjustment Coefficient at the Output

            inputCoeffModalBasis;   % Coefficients of the Modal Basis at the Input

            outputCoeffModalBasis;  % Coefficients of the Modal Basis at the Output

            labelInterBoundCond;    % Label of the Type of Boundary Condition
                                    % at the Interface

            nodeWeightInterface;    % Weight of the Node Used in the Current
                                    % Computation. It is commmented because was
                                    % previously defined
                                    
            timeInstant;
            
            infBoundCond;
            
            outBoundCond;
            
            augVerNodes;
            
            augVerWeights;
            
            bcStruct;
            
            bcInfTag;
            
            bcOutTag;
            
            uRid;
            
            numbControlPts;
            
            coefficientForm;

    end
    
    methods (Access = public)
        
        %% BOUNDARY CONDITION HANDLER - CONSTRUCT METHOD
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
        
        %% BOUNDARY CONDITION HANDLER - FUNCTIONAL METHODS
        % The following functions correspond to the individual scripts
        % previously developed of the functions 'liftBoundCond',
        % 'inflowInterpol', 'interfaceInterpol', 'imposeBoundaryP1'
        % and 'imposeBoundaryP2'.
        
            %% Method 'computeFourierCoeff'
            
            function [infStruct,outStruct] = computeFourierCoeff(obj)
                
                %% IMPORT CLASSES
            
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% IMPORT FEATURES

                infBC = obj.infBoundCond;
                outBC = obj.outBoundCond;
                nodes = obj.augVerNodes;
                wghts = obj.augVerWeights;
                modBasis = obj.modalBasis;
                
                %% COMPUTE REFINED QUADRATURE NODES
                
                obj_gaussLegendre_2 = IntegrateHandler();
                obj_gaussLegendre_2.numbQuadNodes = 128;
                [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);  
                
                %% COMPUTE REFINED MODAL BASIS
                
                obj_newModalBasis = BasisHandler();
            
                obj_newModalBasis.dimModalBasis = obj.dimModalBasis;
                obj_newModalBasis.evalNodesY = verGLNodes;
                obj_newModalBasis.labelUpBoundCond = 'dir';
                obj_newModalBasis.labelDownBoundCond = 'dir';
                obj_newModalBasis.coeffForm = obj.coefficientForm;

                [refModalBasis, ~] = newModalBasis(obj_newModalBasis);
                
                %% COMPUTE PROJECTION
                
                valueInfBC = infBC(verGLNodes);
                valueOutBC = outBC(verGLNodes);
                
                infStruct = zeros(obj.dimModalBasis,1);
                outStruct = zeros(obj.dimModalBasis,1);
                
                for ii = 1:obj.dimModalBasis
                    
                    infStruct(ii) = (valueInfBC .* refModalBasis(:,ii))' * verWeights;
                    outStruct(ii) = (valueOutBC .* refModalBasis(:,ii))' * verWeights;                                                                                                                          
                    
                end
                
                aux1 = infStruct;
                aux2 = outStruct;
                infStruct = fliplr(aux1);
                outStruct = fliplr(aux2);
                
            end
            
            %% Method 'computeFourierCoeff3D'
            
            function [infStruct,outStruct] = computeFourierCoeff3D(obj)
                
                %% IMPORT CLASSES
            
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% IMPORT FEATURES

                infBC = obj.infBoundCond;
                outBC = obj.outBoundCond;
                nodes = obj.augVerNodes;
                wghts = obj.augVerWeights;
                modBasis = obj.modalBasis;
                
                weightMat = wghts * wghts';
                [X,Y]   = meshgrid(nodes,nodes);
                
                %% COMPUTE PROJECTION
                
                valueInfBC = infBC(X,Y);
                valueOutBC = outBC(X,Y);
                
                infStruct = [];
                outStruct = [];
                
                for ii = 1:size(modBasis,3)
                    
                    auxInf = sum(sum(valueInfBC .* modBasis(:,:,ii) .* weightMat));
                    auxOut = sum(sum(valueOutBC .* modBasis(:,:,ii) .* weightMat));
                    
                    infStruct(ii) = auxInf;
                    outStruct(ii) = auxOut;                                                                                                                         
                    
                end
                
                aux1 = infStruct;
                aux2 = outStruct;
                infStruct = fliplr(aux1);
                outStruct = fliplr(aux2);                
                
            end
            
            %% Method 'buildBoundCond'
            
            function [uAug] = buildBoundCond(obj)
                
                BC_l    = obj.bcStruct.bcInfTag;
                BC_r    = obj.bcStruct.bcOutTag;
                nx      = obj.bcStruct.numbControlPts;
                size_mb = obj.bcStruct.dimModalBasis;
                uAug    = zeros(size_mb * nx);
                
                if (strcmp(BC_l,'dir') && strcmp(BC_r,'dir'))
    
                    for imb = 1 : size_mb     % Per ogni frequenza
                        
                        uAug((imb-1) * nx + 1)                = obj.bcStruct.infStruct(imb);
                        uAug((imb-1) * nx + 2 : imb * nx - 1) = obj.uRid((imb-1) * (nx-2) + 1 : imb * (nx-2));
                        uAug(imb * nx)                        = obj.bcStruct.outStruct(imb);

                    end

                elseif (strcmp(BC_l,'dir') && strcmp(BC_r,'neu'))
                    
                    for imb = 1 : size_mb

                        uAug((imb-1) * nx + 1)            = obj.bcStruct.infStruct(imb);
                        uAug((imb-1) * nx + 2 : imb * nx) = obj.uRid((imb-1) * (nx-1) + 1 : imb * (nx-1));
                        
                    end

                elseif (strcmp(BC_l,'neu') && strcmp(BC_r,'dir'))

                    for imb = 1 : size_mb

                        uAug(imb * nx)                        = obj.bcStruct.outStruct(imb);
                        uAug((imb-1) * nx + 1 : imb * nx - 1) = obj.uRid((imb-1) * (nx-1) + 1 : imb * (nx-1));
                        
                    end
                    
                elseif ( strcmp(BC_l,'rob') && (strcmp(BC_r,'rob') || (strcmp(BC_r,'neu') ) ) ) %robrob e robneu condividono lo stesso codice.

                    uAug = obj.uRid;
                    
                elseif ( strcmp(BC_l,'neu') && (strcmp(BC_r,'neu') ) )
                    
                    uAug = obj.uRid;
                    
                elseif ( strcmp(BC_l,'dir') && strcmp(BC_r,'rob') )

                    for imb = 1 : size_mb

                        uAug((imb-1) * nx + 1)            = obj.bcStruct.infStruct(imb);
                        uAug((imb-1) * nx + 2 : imb * nx) = obj.uRid((imb-1) * (nx-1) + 1 : imb * (nx-1));
                        
                    end
                    
                end
                
            end
        
            %% Method 'liftBoundCond'
            
            function [liftCoeffA,liftCoeffB] = liftBoundCond(obj)

                %%
                % liftBoundCond - This function computes the coefficients corresponding to
                %                 the linear offset adjustment of the boundary conditions.
                %
                % The inputs are:
                %%
                %   (1) labelUpBoundCond         : Contains the Label Identifying the Nature of
                %                                  the Boundary Conditions on the Upper Limit of
                %                                  the Domain
                %   (2) labelDownBoundCond       : Contains the Label Identifying the Nature of
                %                                  the Boundary Conditions on the Lower Limit of
                %                                  the Domain
                %   (3) dataUpBoundCond          : Numerical Value of the Upper Bound of the
                %                                  Boundary Conditions
                %   (4) dataDownBoundCond        : Numerical Value of the Lower Bound of the
                %                                  Boundary Conditions
                %   (5) coeffForm                : Data Strusture Containing All the @-Functions
                %                                  and the Constants Relative to the Bilinear Form
                %
                % The outputs are:
                %%
                %   (1) liftCoeffA   : Coefficient of the Line Used to Corret the Boundary 
                %                      Conditions
                %   (2) liftCoeffB   : Coefficient of the Line Used to Corret the Boundary 
                %                      Conditions
                %---------------------------------------------------------------------%

                % Evaluation of the Mu Coefficient

                if(nargin(obj.coeffForm.mu) == 4) muval = obj.coeffForm.mu(0,0,0,0);
                elseif(nargin(obj.coeffForm.mu) == 3) muval = obj.coeffForm.mu(0,0,0);
                else muval = obj.coeffForm.mu(0,0);
                end
                
                % Evaluation of the Robin Coefficient

                coeffrobin = obj.coeffForm.coeffrobin;

                switch [obj.labelUpBoundCond,obj.labelDownBoundCond]

                    % Dirichlet-Dirichlet Boundary Condition Case

                    case 'dirdir'
                        liftCoeffA = obj.dataUpBoundCond - obj.dataDownBoundCond;
                        liftCoeffB = obj.dataDownBoundCond;

                    % Dirichlet-Robin Boundary Condition Case

                    case 'dirrob'
                        if(coeffrobin + muval == 0)
                            display('In liftBoundCond: Non Homogeneous Dirichlet-Robin Conditions Requires a Quadratic Adjustment Not Implemented Yet.')
                        end

                        liftCoeffB = (muval * obj.dataUpBoundCond + obj.dataDownBoundCond)/(muval + coeffrobin);
                        liftCoeffA = (obj.dataUpBoundCond - liftCoeffB);

                    % Robin-Robin Boundary Condition Case

                    case 'robrob'
                        if(2*muval*coeffrobin+coeffrobin^2 == 0)
                            display('In liftBoundCond: Non Homogeneous Robin-Robin Conditions Requires a Quadratic Adjustment Not Implemented Yet.')
                        end
                     
                            liftCoeffA = ( obj.dataUpBoundCond - obj.dataDownBoundCond )/(coeffrobin + 2*muval);
                            liftCoeffB = ( obj.dataDownBoundCond + muval*liftCoeffA )/coeffrobin;

                    % Robin-Dirichlet Boundary Condition Case

                    case 'robdir'
                        if(coeffrobin + muval == 0)
                            display('In liftBoundCond: Non Homogeneous Robin-Dirichlet Conditions Requires a Quadratic Adjustment Not Implemented Yet.')
                        end
                        
                        liftCoeffB = obj.dataDownBoundCond;
                        liftCoeffA = ( obj.dataUpBoundCond - coeffrobin * obj.dataDownBoundCond )/( coeffrobin + muval );

                    % Neumann-Neumann Boundary Condition Case    

                    case 'neuneu'
                        if( obj.dataUpBoundCond == 0 && obj.dataDownBoundCond == 0)
                            liftCoeffA = 0;
                            liftCoeffB = 0;
                        else
                            display('In liftBoundCond: Non homogeneus Neumann BC not available.')
                        end

                    otherwise
                        disp('In liftBoundCond: Boundary Conditions Not Recognized')
                end
                
            end
        
            %% Method 'inflowInterpol'
            
            function newBoundProfile = inflowInterpol(obj)

                %%
                % inflowInterpol - This function computes the boudnary with the offset
                %                  adjustment in the boundary conditions.
                %
                % The inputs are:
                %%
                %   (1) liftAdjust        : Offset Adjustment to the Boundary Conditions
                %   (2) modalBasis        : Modal Basis Used to Represent the Vertical
                %                           Direction
                %   (3) nodeWeightInflow  : Weight of the Node Used in the Current
                %                           Computation
                %   (4) dataDirBoundCond  : Data Structure Containing the Dirichlet
                %                           Boundary Conditions
                %
                % The outputs are:
                %%
                %   (1) newBoundProfile   : Reshaped Boundary with New Boundary Conditions
                %
                % Notes:
                %%
                %   (1) The linear offset 'ril' is already present in the quadrature
                %       nodes.
                %---------------------------------------------------------------------%

                m = size(obj.modalBasis,2);
                
                newBoundProfile = zeros(m,1);

                for i = 1 : m

                    %------------------------------------------------------------------%
                    % Note:
                    % The offset adjustment 'liftAdjust' is taken at the inflow boundary
                    %------------------------------------------------------------------%

                        % Creation of an Object for the IntegrateHandler
                        % Class
                        
                        obj_int = IntegrateHandler();
                        
                        % Properties Assignment
                        
                        obj_int.funcToIntegrate = ((obj.dataDirBoundCond - obj.liftAdjustment).*obj.modalBasis(:,i))';
                        obj_int.funcWeight = obj.nodeWeightInflow;
                        
                        % Call of the Integrate Method
                    
                        newBoundProfile(i) = integrate(obj_int);

                end
                
                disp('Finished INFLOW INTERPOLATION');
                
            end
            
            %% Method 'interfaceInterpol'
            
            function coeffProjData = interfaceInterpol(obj)

                %%
                % interfaceInterpol - This function takes the data corrisponding to the
                %                     adjacent sub-domains (offset adjustment already 
                %                     considered) and projects it into the basis of each
                %                     i-th specif sub-domain. 
                %
                % The inputs are:
                %%
                %   (1) boundModalBasis         : Modal Basis Facing the Boundary
                %   (2) inputLiftAdjust         : Offset Adjustment Coefficient at the Input
                %   (3) outputLiftAdjust        : Offset Adjustment Coefficient at the Output
                %   (4) inputCoeffModalBasis    : Coefficients of the Modal Basis at the Input
                %   (5) outputCoeffModalBasis   : Coefficients of the Modal Basis at the Output
                %   (6) nodeWeightInterface     : Weights of the Node Used in the Current
                %                                 Computation
                %   (7) labelInterBoundCond     : Label of the Type of Boundary Condition
                %                                 at the Interface
                %
                % The outputs are:
                %%
                %   (1) coeffProjData  : Coefficients of the Projected Data in the
                %                        Basis of the i-th Sub-Domain
                %
                % Notes:
                %%
                %   (1) The linear offset 'liftAdjust' is already present in the quadrature
                %       nodes.
                %   (2) This code can be optimized if we consider that in the case of
                %       equal basis, we only have to project the offset adjustment. If
                %       also the offset adjustment is equal, no extra computation is
                %       required.
                %   (3) This optimization can be placed before the definition of the
                %       function output 'coeff_output'. This way, we can first check 
                %       the need for the computations before enterring the 
                %       'interfaceInterpol'.
                %---------------------------------------------------------------------%

                % Number of Elements based on the Domain of Interest (sx)

                m_out = size( obj.outputCoeffModalBasis , 2 ); 

                % Number of Elements based on the Faced Domain (dx)

                m_in  = size( obj.inputCoeffModalBasis  , 2 ); 

                % Initialization of the Output Vector

                coeffProjData = zeros(m_out, 1 );

                if (strcmp(obj.labelInterBoundCond,'dir')||strcmp(obj.labelInterBoundCond,'rob'))

                    % Difference between the Output and Input Offset Adjustments

                    rill   = obj.inputLiftAdjust - obj.outputLiftAdjust;

                end

                for i = 1 : m_out

                    if (strcmp(obj.labelInterBoundCond,'dir') || strcmp(obj.labelInterBoundCond,'rob'))
                        
                        % Creation of an Object for the IntegrateHandler
                        % Class
                        
                        obj_int_out = IntegrateHandler();
                        
                        % Properties Assignment
                        
                        obj_int_out.funcToIntegrate = (rill.*obj.outputCoeffModalBasis(:,i))';
                        obj_int_out.funcWeight = obj.nodeWeightInterface;
                        
                        % Call of the Integrate Method
                        
                        coeffProjData(i) = integrate(obj_int_out);
                        
                    end

                    for j = 1 : m_in
                        
                        % Creation of an Object for the IntegrateHandler
                        % Class
                        
                        obj_int_in = IntegrateHandler();
                        
                        % Properties Assignment
                        
                        obj_int_in.funcToIntegrate = (obj.outputCoeffModalBasis(:,i).*obj.inputCoeffModalBasis(:,j))';
                        obj_int_in.funcWeight = obj.nodeWeightInterface;
                        
                        % Call of the Integrate Method
                        
                        coeffProjData(i) = coeffProjData(i) + obj.boundModalBasis(j) * integrate(obj_int_in);
                    end    
                end
                
                disp('Finished INTERFACE INTERPOLATION');
                
            end
            
            %% Method 'weakImposeBoundaryConditions3D'
            
            function [stiffMat,rhsVect] = weakImposeBoundaryConditions3D(obj)
                
                % EXTRACT THE INFORMATION FROM THE BOUNDARY STRUCTURE
                % RECEIVED AS IMPUT TO THE METHOD
                
                inflowINFO  = obj.boundaryStruct.inflowINFO;
                outflowINFO = obj.boundaryStruct.inflowINFO;
                inflowDATAdir  = obj.boundaryStruct.inflowINFO;
                inflowDATAneu  = obj.boundaryStruct.inflowINFO;
                outflowDATAdir = obj.boundaryStruct.inflowINFO;
                outflowDATAneu = obj.boundaryStruct.inflowINFO;
                dimModes = obj.boundaryStruct.dimModes;
                dimNodes = obj.boundaryStruct.dimNodes;
                stiffMat = obj.boundaryStruct.stiffMat;
                rhsVect = obj.boundaryStruct.rhsVect;
                coeffModalBasis = obj.boundaryStruct.coeffModalBasis;
                weightMat = obj.boundaryStruct.weightMat;
                jacFunc = obj.boundaryStruct.jacFunc;
                
                % IMPOSE INFLOW BOUNDARY CONDITIONS
                
                switch inflowINFO
                    
                    case 'dir'
                        
                        % DEFINITION OF THE COEFFICIENT USED IN THE WEAK
                        % IMPOSITION OF THE BOUNDARY CONDITIONS
                        
                        gamma = 1e30 * max(max(max(abs(stiffMat)),max(abs(rhsVect))));
                        
                        % IMPOSITION OF THE DIRICHLET BOUNDARY CONDITIONS 

                        for kk = 1:dimModes
                            stiffMat( (kk-1)*dimNodes + 1, (kk-1)*dimNodes + 1) = gamma;
                            rhsVect( (kk-1)*dimNodes + 1 ) = gamma * inflowDATAdir(kk);
                        end

                    case 'neu'
                        
                        % IN THE CASE OF NEUMANN BOUNDARY CONDITION IN
                        % INFLOW WE JUST HAVE TO UPDATE THE RIGHT HAND SIDE
                        % TERM TO CONTEMPLATE THE MODAL BASIS INTEGRAL
                        
                        for ii = 1:size(coeffModalBasis,2)
                            rhsVect(1 + (ii-1) * dimNodes) = ...
                            rhsVect(1 +(ii-1) * dimNodes) + ... 
                            inflowDATAneu * jacDet(:,1)' * (coeffModalBasis(:,:,ii) .* weightMat);
                        end
                        
                end
                
                % IMPOSE OUTFLOW BOUNDARY CONDITIONS
                
                switch outflowINFO
                    
                    case 'dir'
                        
                        % DEFINITION OF THE COEFFICIENT USED IN THE WEAK
                        % IMPOSITION OF THE BOUNDARY CONDITIONS
                        
                        gamma = 1e30 * max(max(max(abs(stiffMat)),max(abs(rhsVect))));
                        
                        % IMPOSITION OF THE DIRICHLET BOUNDARY CONDITIONS 

                        for kk = 1:dimModes
                            stiffMat( kk * dimNodes , kk * dimNodes) = gamma;
                            rhsVect( kk * dimNodes ) = gamma * outflowDATAdir(kk);
                        end

                    case 'neu'
                        
                        % IN THE CASE OF NEUMANN BOUNDARY CONDITION IN
                        % INFLOW WE JUST HAVE TO UPDATE THE RIGHT HAND SIDE
                        % TERM TO CONTEMPLATE THE MODAL BASIS INTEGRAL
                        
                        for ii = 1:size(coeffModalBasis,3)
                            rhsVect(kk * dimNodes) = ...
                            rhsVect(kk * dimNodes) + ... 
                            outflowDATAneu * jacDet(:,1)' * (coeffModalBasis(:,ii) .* weightMat);
                        end
                        
                end
            
            end
            
            %% Method 'imposeBoundaryP1'
            
            function [reducedMatA, reducedVecB] = imposeBoundaryP1(obj)
                               %%
                % imposeBoundaryP1 - This function imposes all the boundary
                %                    conditions relative to the interface, inflow 
                %                    and outflow of the fluid. If necessary, it is 
                %                    possible to reduce the problem by eliminating 
                %                    the Dirichlet's Degrees of Freedom.
                %
                % The inputs are:
                %%
                %   (1)  dimModalBasis        : Dimension of the Modal Basis
                %   (2)  leftLimitDomainX     : Left Limit of the Domain in the X Direction
                %   (3)  rightLimitDomainX    : Right Limit of the Domain in the X Direction
                %                               Boundary Conditions
                %   (4)  stepMeshX            : Vector Containing the Step of the Finite
                %                                Element Mesh
                %   (5)  leftLabelBoundCond   : Boundary Conditions on the Left Bound of the
                %                               Domain
                %   (6)  rightLabelBoundCond  : Boundary Conditions on the Right Bound of the
                %                               Domain
                %   (7)  leftDomainBoundary   : Left Boundary of the Domain
                %   (8)  rightDomainBoundary  : Right Boundary of the Domain
                %   (9)  blockMatA, blockVecB : Block Matrices Defining the Problem
                %
                % The outputs are:
                %%
                %   (1) reducedMatA   : Block Matrix of the Problem with all the Boundary
                %                       Conditions Imposed
                %   (2) reducedVecB   : Known Block Vector of the Problem with all the
                %                       Boundary Conditions Imposed

                ne = round((obj.rightLimitDomainX - obj.leftLimitDomainX)/obj.stepMeshX);  % Number of Intervals
                nx = ne+1;                                                                 % Number of Nodes

                %---------------------------------------------------------------------%
                % FEM MESH IN THE X DIRECTION
                % Note: In this case we are considering equispaced nodes in the mesh.
                %---------------------------------------------------------------------%

                meshx = zeros(nx,1);
                meshx(1) = obj.leftLimitDomainX;

                for i=2:nx

                    meshx(i) = meshx(i-1) + obj.stepMeshX;

                end

                %---------------------------------------------------------------------%
                % INITIALIZATION OF THE 'brid' (Reduced b)
                % Note: 
                % The blocks of 'brid' will be added considering each frequency
                % of the decomposed problem.
                %---------------------------------------------------------------------%

                reducedVecB=[];

                %---------------------------------------------------------------------%
                % TREATMENT OF THE BOUNDARY
                % Note: 
                % The loop goes back into the orginal global matrix and takes the
                % values that are usefull based on the boundary conditions defined for
                % the problem.
                %---------------------------------------------------------------------%

                %-----------------------------------------%
                % DIRICHLET-DIRICHLET BOUNDARY CONDITIONS %
                %-----------------------------------------%

                if (strcmp(obj.leftLabelBoundCond,'dir') && strcmp(obj.rightLabelBoundCond,'dir'))

                    %-----------------------------------------------------------------%
                    % Note:
                    % The value 'size_mb' indicates the number of frequencies
                    % considered for the modal basis.
                    %-----------------------------------------------------------------%
                    
                    for imb = 1 : obj.dimModalBasis

                        buff = zeros( nx-2, 1);

                        %-------------------------------------------------------------%
                        % Note: 
                        % Initialization of the contribution to the Known Term of the 
                        % Dirichlet node, as large as the number of internal nodes.   
                        %-------------------------------------------------------------%

                        for jmb = 1 : obj.dimModalBasis

                            reducedMatA( (imb-1)*(nx-2) + 1 : imb*(nx-2) , (jmb-1)*(nx-2) + 1 : jmb*(nx-2) ) = ...
                                obj.blockMatA( (imb-1)*nx + 2     : imb*nx - 1 , (jmb-1)*nx + 2     : jmb*nx - 1 );

                            buff = buff + obj.blockMatA( (imb-1)*nx + 2 : imb*nx-1, (jmb-1)*nx + 1 ) * obj.leftDomainBoundary( jmb ) ...
                                + obj.blockMatA( (imb-1)*nx + 2 : imb*nx-1, jmb*nx         ) * obj.rightDomainBoundary( jmb );

                        end

                        reducedVecB = [ reducedVecB; obj.blockVecB( (imb-1)*nx + 2 : imb*nx - 1 ) - buff ];

                    end

                %---------------------------------------%
                % DIRICHLET-NEUMANN BOUNDARY CONDITIONS %
                %---------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'dir') && strcmp(obj.rightLabelBoundCond,'neu') )
                    
                    for imb = 1 : obj.dimModalBasis

                        buff = zeros( nx-1, 1);

                        for jmb = 1 : obj.dimModalBasis
                            reducedMatA( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                                obj.blockMatA( (imb-1)*nx + 2     : imb*nx     , (jmb-1)*nx + 2     : jmb*nx     );

                            buff = buff + obj.blockMatA( (imb-1)*nx + 2 : imb*nx, (jmb-1)*nx + 1 ) * obj.leftDomainBoundary(jmb);

                        end

                        reducedVecB = [ reducedVecB; obj.blockVecB( (imb-1)*nx + 2 : imb*nx ) - buff ];

                        reducedVecB(end) = reducedVecB( end ) + obj.rightDomainBoundary( imb ); % Neumann Data
                    end

                %---------------------------------------%
                % NEUMANN-DIRICHLET BOUNDARY CONDITIONS %
                %---------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'neu') && strcmp(obj.rightLabelBoundCond,'dir') )

                    for imb = 1 : obj.dimModalBasis

                        buff = zeros( nx-1, 1);

                        for jmb = 1 : obj.dimModalBasis
                            reducedMatA( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                                obj.blockMatA( (imb-1)*nx + 1     : imb*nx-1   , (jmb-1)*nx + 1     : jmb*nx - 1 );

                            buff = buff + obj.blockMatA( (imb-1)*nx + 1 : imb*nx - 1 , jmb*nx ) * obj.rightDomainBoundary( jmb );
                        end

                        reducedVecB = [ reducedVecB ; obj.blockVecB( (imb-1)*nx + 1 : imb*nx - 1 ) - buff ];

                        reducedVecB(end-nx+2) = reducedVecB(end-nx+2) + obj.leftDomainBoundary( imb );

                    end

                %------------------------------------------------%
                % ROBIN-ROBIN, ROBIN-NEUMANN BOUNDARY CONDITIONS %
                %------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % In the case of Robin-Robin and Robin-Neumann boundary conditions, the
                % code used to process the contribution of the known term is the same.
                %---------------------------------------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'rob') && (strcmp(obj.rightLabelBoundCond,'rob') || (strcmp(obj.rightLabelBoundCond,'neu') ) ) )

                    reducedMatA = obj.blockMatA;
                    reducedVecB = obj.blockVecB;
                    
                    for imb = 1 : obj.dimModalBasis
                        
                        reducedVecB((imb-1)*nx+1) = reducedVecB((imb-1)*nx+1) + obj.leftDomainBoundary ( imb );
                        reducedVecB(imb*nx) = reducedVecB(imb*nx) + obj.rightDomainBoundary( imb );
                        
                    end

                %-------------------------------------%
                % NEUMANN-NEUMANN BOUNDARY CONDITIONS %
                %-------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'neu') && (strcmp(obj.rightLabelBoundCond,'neu') ) )
                    
                    reducedMatA = obj.blockMatA;
                    reducedVecB = obj.blockVecB;
                    
                    for imb = 1 : obj.dimModalBasis
                        
                        reducedVecB((imb-1)*nx+1) = reducedVecB((imb-1)*nx+1) + obj.leftDomainBoundary( imb );
                        reducedVecB(imb*nx) = reducedVecB(imb*nx) + obj.rightDomainBoundary( imb );
                        
                    end

                %-------------------------------------%
                % DIRICHLET-ROBIN BOUNDARY CONDITIONS %
                %-------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'dir') && strcmp(obj.rightLabelBoundCond,'rob') )

                    for imb = 1 : obj.dimModalBasis

                        buff = zeros( nx-1, 1);

                        for jmb = 1 : obj.dimModalBasis
                            
                            reducedMatA( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                                obj.blockMatA( (imb-1)*nx + 2     : imb*nx     , (jmb-1)*nx + 2     : jmb*nx     );
                            
                            buff = buff + obj.blockMatA( (imb-1)*nx + 2 : imb*nx, (jmb-1)*nx + 1 ) * obj.leftDomainBoundary(jmb);

                        end

                        reducedVecB = [ reducedVecB; obj.blockVecB( (imb-1)*nx + 2 : imb*nx ) - buff ];
                        reducedVecB(end) = reducedVecB( end ) + obj.rightDomainBoundary( imb ); % Robin Robin Data
                    end
                end
            end
     
            %% Method 'imposeBoundaryP2'
            
            function [reducedMatA, reducedVecB] = imposeBoundaryP2(obj)

                %%
                % imposeBoundaryP2   - This function imposes all the boundary
                %                      conditions relative to the interface, inflow 
                %                      and outflow of the fluid. If necessary, it is 
                %                      possible to reduce the problem by eliminating 
                %                      the Dirichlet's Degrees of Freedom.
                %
                % The inputs are:
                %%
                %   (1)  dimModalBasis         : Dimension of the Modal Basis
                %   (2)  leftLimitDomainX      : Left Limit of the Domain in the X Direction
                %   (3)  rightLimitDomainX     : Right Limit of the Domain in the X Direction
                %                                Boundary Conditions
                %   (4)  stepMeshX             : Vector Containing the Step of the Finite
                %                                Element Mesh
                %   (5)  leftLabelBoundCond    : Boundary Conditions on the Left Bound of the
                %                                Domain
                %   (6)  rightLabelBoundCond   : Boundary Conditions on the Right Bound of the
                %                                Domain
                %   (7)  leftDomainBoundary    : Left Boundary of the Domain
                %   (8)  rightDomainBoundary   : Right Boundary of the Domain
                %   (9)  blockMatA, blockVecB  : Block Matrices Defining the Problem
                %   (10) degreePolySplineBasis : Degree of the Polynomial B-Spline Basis
                %   (11) continuityParameter   : Degree of Continuity of the Basis 'C^(p-k)'
                %
                % The outputs are:
                %%
                %   (1) reducedMatA   : Block Matrix of the Problem with all the Boundary
                %                       Conditions Imposed
                %   (2) reducedVecB   : Known Block Vector of the Problem with all the
                %                       Boundary Conditions Imposed
                %---------------------------------------------------------------------%

                ne = round((obj.rightLimitDomainX - obj.leftLimitDomainX)/obj.stepMeshX);  % Number of Intervals
                nx = ne+1;                                                                 % Number of Nodes

                %---------------------------------------------------------------------%
                % FEM MESH IN THE X DIRECTION
                % Note: In this case we are considering equispaced nodes in the mesh.
                %---------------------------------------------------------------------%

                meshx = zeros(nx,1);
                meshx(1) = obj.leftLimitDomainX;
                
                nx = ne * obj.continuityParameter + obj.degreePolySplineBasis + 1 - obj.continuityParameter;
                hxi = obj.stepMeshX/2;

                for i=2:nx
                    meshx(i) = meshx(i-1) + hxi;
                end

                %---------------------------------------------------------------------%
                % INITIALIZATION OF THE 'brid' (Reduced b)
                % Note: 
                % The blocks of 'brid' will be added considering each frequency
                % of the decomposed problem.
                %---------------------------------------------------------------------%

                reducedVecB = [];

                %---------------------------------------------------------------------%
                % TREATMENT OF THE BOUNDARY
                % Note: 
                % The loop goes back into the orginal global matrix and takes the
                % values that are usefull based on the boundary conditions defined for
                % the problem.
                %---------------------------------------------------------------------%

                %-----------------------------------------%
                % DIRICHLET-DIRICHLET BOUNDARY CONDITIONS %
                %-----------------------------------------%

                if (strcmp(obj.leftLabelBoundCond,'dir') && strcmp(obj.rightLabelBoundCond,'dir'))

                    %-----------------------------------------------------------------%
                    % Note:
                    % The value 'size_mb' indicates the number of frequencies
                    % considered for the modal basis.
                    %-----------------------------------------------------------------%

                    for imb = 1 : obj.dimModalBasis 

                        buff = zeros( nx-2, 1);

                        %-------------------------------------------------------------%
                        % Note: 
                        % Initialization of the contribution to the Known Term of the 
                        % Dirichlet node, as large as the number of internal nodes.   
                        %-------------------------------------------------------------%

                        for jmb = 1 : obj.dimModalBasis

                            reducedMatA( (imb-1)*(nx-2) + 1 : imb*(nx-2) , (jmb-1)*(nx-2)+1 : jmb*(nx-2) ) = ...
                               obj.blockMatA( (imb-1)*nx + 2     : imb*nx - 1 , (jmb-1)*nx + 2     : jmb*nx - 1 );

                            buff = buff + obj.blockMatA( (imb-1)*nx + 2 : imb*nx-1, (jmb-1)*nx + 1 )*obj.leftDomainBoundary( jmb ) ...
                                        + obj.blockMatA( (imb-1)*nx + 2 : imb*nx-1, jmb*nx         )*obj.rightDomainBoundary( jmb );

                        end

                        reducedVecB = [ reducedVecB; obj.blockVecB( (imb-1)*nx + 2 : imb*nx - 1 ) - buff ];

                    end

                %---------------------------------------%
                % DIRICHLET-NEUMANN BOUNDARY CONDITIONS %
                %---------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'dir') && strcmp(obj.rightLabelBoundCond,'neu') )
                    for imb = 1 : obj.dimModalBasis

                        buff = zeros( nx-1, 1);

                        for jmb = 1 : obj.dimModalBasis
                            
                            reducedMatA( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                               obj.blockMatA( (imb-1)*nx + 2     : imb*nx     , (jmb-1)*nx + 2     : jmb*nx     );

                            buff = buff + obj.blockMatA( (imb-1)*nx + 2 : imb*nx, (jmb-1)*nx + 1 ) * obj.leftDomainBoundary(jmb);

                        end
                        
                        reducedVecB = [ reducedVecB; obj.blockVecB( (imb-1)*nx + 2 : imb*nx ) - buff ];			
                        
                        reducedVecB(end) = reducedVecB( end ) + obj.rightDomainBoundary( imb );   % Neumann Data
                    end

                %---------------------------------------%
                % NEUMANN-DIRICHLET BOUNDARY CONDITIONS %
                %---------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'neu') && strcmp(obj.rightLabelBoundCond,'dir') )

                    for imb = 1 : obj.dimModalBasis

                        buff = zeros( nx-1, 1);

                        for jmb = 1 : obj.dimModalBasis
                            
                            reducedMatA( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                               obj.blockMatA( (imb-1)*nx + 1     : imb*nx-1   , (jmb-1)*nx + 1     : jmb*nx - 1 );

                            buff = buff + obj.blockMatA( (imb-1)*nx + 1 : imb*nx - 1 , jmb*nx ) * obj.rightDomainBoundary( jmb );
                        end

                        reducedVecB = [ reducedVecB ; obj.blockVecB( (imb-1)*nx + 1 : imb*nx - 1 ) - buff ];

                        reducedVecB(end-nx+2) = reducedVecB(end-nx+2) + obj.leftDomainBoundary( imb );

                    end

                %------------------------------------------------%
                % ROBIN-ROBIN, ROBIN-NEUMANN BOUNDARY CONDITIONS %
                %------------------------------------------------%

                %---------------------------------------------------------------------%
                % Note:
                % In the case of Robin-Robin and Robin-Neumann boundary conditions, the
                % code used to process the contribution of the known term is the same.
                %---------------------------------------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'rob') && (strcmp(obj.rightLabelBoundCond,'rob') || (strcmp(obj.rightLabelBoundCond,'neu') ) ) )

                    reducedMatA = obj.blockMatA;
                    reducedVecB = obj.blockVecB;
                    
                    for imb = 1 : obj.dimModalBasis
                        
                        reducedVecB((imb-1)*nx+1) = reducedVecB((imb-1)*nx+1) + obj.leftDomainBoundary( imb );
                        reducedVecB(imb*nx) = reducedVecB(imb*nx) + obj.rightDomainBoundary( imb );		
                    end

                %-------------------------------------%
                % NEUMANN-NEUMANN BOUNDARY CONDITIONS %
                %-------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'neu') && (strcmp(obj.rightLabelBoundCond,'neu') ) )
                    
                    reducedMatA = obj.blockMatA;
                    reducedVecB = obj.blockVecB;
                    
                    for imb = 1 : obj.dimModalBasis
                        
                        reducedVecB((imb-1)*nx+1) = reducedVecB((imb-1)*nx+1) + obj.leftDomainBoundary( imb );
                        reducedVecB(imb*nx) = reducedVecB(imb*nx) + obj.rightDomainBoundary( imb );		
                        
                    end

                %-------------------------------------%
                % DIRICHLET-ROBIN BOUNDARY CONDITIONS %
                %-------------------------------------%

                elseif ( strcmp(obj.leftLabelBoundCond,'dir') && strcmp(obj.rightLabelBoundCond,'rob') )

                    for imb = 1 : obj.dimModalBasis

                        buff = zeros( nx-1, 1);

                        for jmb = 1 : obj.dimModalBasis
                            
                            reducedMatA( (imb-1)*(nx-1) + 1 : imb*(nx-1) , (jmb-1)*(nx-1) + 1 : jmb*(nx-1) ) = ...
                            obj.blockMatA( (imb-1)*nx + 2     : imb*nx     , (jmb-1)*nx + 2     : jmb*nx     );
                            buff = buff + obj.blockMatA( (imb-1)*nx + 2 : imb*nx, (jmb-1)*nx + 1 ) * obj.leftDomainBoundary(jmb);

                        end

                        reducedVecB = [ reducedVecB; obj.blockVecB( (imb-1)*nx + 2 : imb*nx ) - buff ];		
                        reducedVecB(end) = reducedVecB( end ) + obj.rightDomainBoundary( imb );   % Robin-Robin Data
                        
                    end
                end
                
                disp('Finished BOUNDARY IMPOSITION');
                
            end
            
            %% Method 'computeFourierCoeffStokes'
            
            function [infStruct,outStruct] = computeFourierCoeffStokes(obj)
                
                %% IMPORT CLASSES
            
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% IMPORT FEATURES
                
                infBC    = obj.boundaryStruct.bc_inf_data;
                outBC    = obj.boundaryStruct.bc_out_data;
                nodes    = obj.boundaryStruct.augVerNodes;
                wghts    = obj.boundaryStruct.augVerWeights;
                modBasis = obj.boundaryStruct.modalBasis;
                
                %% COMPUTE PROJECTION
                
                valueInfBC = infBC(nodes);
                valueOutBC = outBC(nodes);
                
                infStruct = [];
                outStruct = [];
                
                for ii = 1:size(modBasis,2)
                    
                    infStruct(ii) = (valueInfBC .* modBasis(:,ii))' * wghts;
                    outStruct(ii) = (valueOutBC .* modBasis(:,ii))' * wghts;                                                                                                                          
                    
                end
                
            end
            
            %% Method 'computeFourierCoeffStokesTransient'
            
            function [infStruct,outStruct] = computeFourierCoeffStokesTransient(obj)
                
                %% IMPORT CLASSES
            
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% IMPORT FEATURES
                
                infBC    = obj.boundaryStruct.bc_inf_data;
                outBC    = obj.boundaryStruct.bc_out_data;
                nodes    = obj.boundaryStruct.augVerNodes;
                wghts    = obj.boundaryStruct.augVerWeights;
                modBasis = obj.boundaryStruct.modalBasis;
                
                %% COMPUTE PROJECTION
                
                valueInfBC = infBC(nodes,obj.boundaryStruct.time);
                valueOutBC = outBC(nodes,obj.boundaryStruct.time);
                
                infStruct = [];
                outStruct = [];
                
                for ii = 1:size(modBasis,2)
                    
                    infStruct(ii) = (valueInfBC .* modBasis(:,ii))' * wghts;
                    outStruct(ii) = (valueOutBC .* modBasis(:,ii))' * wghts;                                                                                                                          
                    
                end
                
            end
            
            %% Method 'imposeBoundaryStokes'
            
            function [A,B,P,f] = imposeBoundaryStokes(obj)
                
                bc_inf_tag  = obj.boundaryStruct.bc_inf_tag;
                bc_out_tag  = obj.boundaryStruct.bc_out_tag;
                numbModes   = obj.boundaryStruct.numbModes;
                numbCtrlPts = obj.boundaryStruct.numbCtrlPts;
                projInfBC   = obj.boundaryStruct.projInfBC;
                projOutBC   = obj.boundaryStruct.projOutBC;
                A           = obj.boundaryStruct.stiffMatrix;
                B           = obj.boundaryStruct.couplingMat;
                P           = obj.boundaryStruct.pressureMat;
                f           = obj.boundaryStruct.forceTerm;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DEFINITION OF THE WEAK BC IMPOSITION %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                delta = 1;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % IMPOSITION OF THE BC IN THE GLOBAL STIFFNESS MATRIX %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % (1) DIRICHLET - DIRICHLET
                
                if (strcmp(bc_inf_tag,'dir') && strcmp(bc_out_tag,'dir'))

                    for imb = 1 : numbModes
                        
                        A((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(A((imb - 1) * numbCtrlPts + 1 , :)));
                        A(imb * numbCtrlPts , :) = zeros(size(A(imb * numbCtrlPts , :)));
                        B((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(B((imb - 1) * numbCtrlPts + 1 , :)));
                        B(imb * numbCtrlPts , :) = zeros(size(B(imb * numbCtrlPts , :)));
                        P((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(P((imb - 1) * numbCtrlPts + 1 , :)));
                        P(imb * numbCtrlPts , :) = zeros(size(P(imb * numbCtrlPts , :)));
                        
                        A((imb - 1) * numbCtrlPts + 1 , (imb - 1) * numbCtrlPts + 1) = delta;
                        A(imb * numbCtrlPts , imb * numbCtrlPts) = delta;
                        
                        f((imb - 1) * numbCtrlPts + 1) = delta * projInfBC(imb);
                        f(imb * numbCtrlPts) = delta * projOutBC(imb);

                    end
                    
                % (2) DIRICHLET - NEUMANN

                elseif ( strcmp(bc_inf_tag,'dir') && strcmp(bc_out_tag,'neu') )
                    
                    for imb = 1 : numbModes
                        
                        A((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(A((imb - 1) * numbCtrlPts + 1 , :)));
                        B((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(B((imb - 1) * numbCtrlPts + 1 , :)));
                        P((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(P((imb - 1) * numbCtrlPts + 1 , :)));
                        
                        A((imb - 1) * numbCtrlPts + 1 , (imb - 1) * numbCtrlPts + 1) = delta;
                        
                        f((imb - 1) * numbCtrlPts + 1) = delta * projInfBC(imb);
                        f(imb * numbCtrlPts) = f(imb * numbCtrlPts) + projOutBC(imb);
                        
                    end

                % (3) NEUMANN - DIRICHLET
                
                elseif ( strcmp(bc_inf_tag,'neu') && strcmp(bc_out_tag,'dir') )

                    for imb = 1 : numbModes
                        
                        A(imb * numbCtrlPts , :) = zeros(size(A(imb * numbCtrlPts , :)));
                        B(imb * numbCtrlPts , :) = zeros(size(B(imb * numbCtrlPts , :)));
                        P(imb * numbCtrlPts , :) = zeros(size(P(imb * numbCtrlPts , :)));
                        
                        A(imb * numbCtrlPts , imb * numbCtrlPts) = delta;
                        
                        f((imb - 1) * numbCtrlPts + 1) = f((imb - 1) * numbCtrlPts + 1) + projInfBC(imb);
                        f(imb * numbCtrlPts) = delta * projOutBC(imb);
                        
                    end
                    
                % (4) ROBIN - ROBIAN & ROBIN - NEUMANN
                    
                elseif ( strcmp(bc_inf_tag,'rob') && (strcmp(bc_out_tag,'rob') || (strcmp(bc_out_tag,'neu') ) ) ) %robrob e robneu condividono lo stesso codice.

                    for imb = 1 : numbModes
                        
                        f((imb - 1) * numbCtrlPts + 1) = f((imb - 1) * numbCtrlPts + 1) + projInfBC(imb);
                        f(imb * numbCtrlPts) = f(imb * numbCtrlPts) + projOutBC(imb);
                        
                    end
                    
                % (5) NEUMANN - NEUMANN
                    
                elseif ( strcmp(bc_inf_tag,'neu') && (strcmp(bc_out_tag,'neu') ) )

                    for imb = 1 : numbModes
                        
                        f((imb - 1) * numbCtrlPts + 1) = f((imb - 1) * numbCtrlPts + 1) + projInfBC(imb);
                        f(imb * numbCtrlPts) = f(imb * numbCtrlPts) + projOutBC(imb);
                        
                    end
                    
                % DIRICHLET - ROBIN
                    
                elseif ( strcmp(bc_inf_tag,'dir') && strcmp(bc_out_tag,'rob') )

                    for imb = 1 : numbModes
                        
                        A((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(A((imb - 1) * numbCtrlPts + 1 , :)));
                        B((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(B((imb - 1) * numbCtrlPts + 1 , :)));
                        P((imb - 1) * numbCtrlPts + 1 , :) = zeros(size(P((imb - 1) * numbCtrlPts + 1 , :)));
                        
                        A((imb - 1) * numbCtrlPts + 1 , (jmb - 1) * numbCtrlPts + 1) = delta;
                        
                        f((imb - 1) * numbCtrlPts + 1) = delta * projInfBC(imb);
                        f(imb * numbCtrlPts) = f(imb * numbCtrlPts) + projOutBC(imb);
                        
                    end
                end
                
            end
            
            %% Method 'computeFourierCoeffStokes3D'
            
            function [bcDataStruct] = computeFourierCoeffStokes3D(obj)
                
                %% IMPORT CLASSES
            
                import Core.AssemblerADRHandler
                import Core.BoundaryConditionHandler
                import Core.IntegrateHandler
                import Core.EvaluationHandler
                import Core.BasisHandler
                import Core.SolverHandler
                
                %% IMPORT FEATURES
                
                infBCUx    = obj.boundaryStruct.bcInfDataUx;
                outBCUx    = obj.boundaryStruct.bcOutDataUx;
                infBCUy    = obj.boundaryStruct.bcInfDataUy;
                outBCUy    = obj.boundaryStruct.bcOutDataUy;
                infBCUz    = obj.boundaryStruct.bcInfDataUz;
                outBCUz    = obj.boundaryStruct.bcOutDataUz;
                
                modBasisUx = obj.boundaryStruct.modalBasisUx;
                modBasisUy = obj.boundaryStruct.modalBasisUy;
                modBasisUz = obj.boundaryStruct.modalBasisUz;
                
                nodesUx     = obj.boundaryStruct.augVerNodesUx;
                wghtsUx     = obj.boundaryStruct.augVerWeightsUx;
                weightMatUx = wghtsUx * wghtsUx';
                [XUx,YUx]   = meshgrid(nodesUx,nodesUx);
                
                nodesUy     = obj.boundaryStruct.augVerNodesUy;
                wghtsUy     = obj.boundaryStruct.augVerWeightsUy;
                weightMatUy = wghtsUy * wghtsUy';
                [XUy,YUy]   = meshgrid(nodesUy,nodesUy);
                
                nodesUz     = obj.boundaryStruct.augVerNodesUz;
                wghtsUz     = obj.boundaryStruct.augVerWeightsUz;
                weightMatUz = wghtsUz * wghtsUz';
                [XUz,YUz]   = meshgrid(nodesUz,nodesUz);
                
                %% COMPUTE PROJECTION Ux
                
                valueInfBCUx = infBCUx(XUx,YUx);
                valueOutBCUx = outBCUx(XUx,YUx);
                
                infStructUx = [];
                outStructUx = [];
                
                for ii = 1:size(modBasisUx,3)
                    
                    auxInfUx = sum(sum(valueInfBCUx .* modBasisUx(:,:,ii) .* weightMatUx));
                    auxOutUx = sum(sum(valueOutBCUx .* modBasisUx(:,:,ii) .* weightMatUx));
                    
                    infStructUx(ii) = auxInfUx;
                    outStructUx(ii) = auxOutUx;                                                                                                                         
                    
                end
                
                bcDataStruct.infStructUx = infStructUx;
                bcDataStruct.outStructUx = outStructUx;
                
                %% COMPUTE PROJECTION Uy
                
                valueInfBCUy = infBCUy(XUy,YUy);
                valueOutBCUy = outBCUy(XUy,YUy);
                
                infStructUy = [];
                outStructUy = [];
                
                for ii = 1:size(modBasisUy,3)
                    
                    auxInfUy = sum(sum(valueInfBCUy .* modBasisUy(:,:,ii) .* weightMatUy));
                    auxOutUy = sum(sum(valueOutBCUy .* modBasisUy(:,:,ii) .* weightMatUy));
                    
                    infStructUy(ii) = auxInfUy;
                    outStructUy(ii) = auxOutUy;                                                                                                                         
                    
                end
                
                bcDataStruct.infStructUy = infStructUy;
                bcDataStruct.outStructUy = outStructUy;
                
                %% COMPUTE PROJECTION Uz
                
                valueInfBCUz = infBCUz(XUz,YUz);
                valueOutBCUz = outBCUz(XUz,YUz);
                
                infStructUz = [];
                outStructUz = [];
                
                for ii = 1:size(modBasisUz,3)
                    
                    auxInfUz = sum(sum(valueInfBCUz .* modBasisUz(:,:,ii) .* weightMatUz));
                    auxOutUz = sum(sum(valueOutBCUz .* modBasisUz(:,:,ii) .* weightMatUz));
                    
                    infStructUz(ii) = auxInfUz;
                    outStructUz(ii) = auxOutUz;                                                                                                                         
                    
                end
                
                bcDataStruct.infStructUz = infStructUz;
                bcDataStruct.outStructUz = outStructUz;
                
            end
            
            %% Method 'imposeBoundaryStokes3D'
            
            function [probDataStruct] = imposeBoundaryStokes3D(obj)
                
                bcInfTagUx  = obj.boundaryStruct.bcInfTagUx;
                bcOutTagUx  = obj.boundaryStruct.bcOutTagUx;
                bcInfTagUy  = obj.boundaryStruct.bcInfTagUy;
                bcOutTagUy  = obj.boundaryStruct.bcOutTagUy;
                bcInfTagUz  = obj.boundaryStruct.bcInfTagUz;
                bcOutTagUz  = obj.boundaryStruct.bcOutTagUz;
                
                numbModesUx   = obj.boundaryStruct.numbModesUx;
                numbCtrlPtsUx = obj.boundaryStruct.numbCtrlPtsUx;
                numbModesUy   = obj.boundaryStruct.numbModesUy;
                numbCtrlPtsUy = obj.boundaryStruct.numbCtrlPtsUy;
                numbModesUz   = obj.boundaryStruct.numbModesUz;
                numbCtrlPtsUz = obj.boundaryStruct.numbCtrlPtsUz;
                
                projInfBCUx = obj.boundaryStruct.projInfBCUx;
                projOutBCUx = obj.boundaryStruct.projOutBCUx;
                projInfBCUy = obj.boundaryStruct.projInfBCUy;
                projOutBCUy = obj.boundaryStruct.projOutBCUy;
                projInfBCUz = obj.boundaryStruct.projInfBCUz;
                projOutBCUz = obj.boundaryStruct.projOutBCUz;
                
                Axx         = obj.boundaryStruct.Axx;
                Ayy         = obj.boundaryStruct.Ayy;
                Azz         = obj.boundaryStruct.Azz;
                Bxy         = obj.boundaryStruct.Bxy;
                Bxz         = obj.boundaryStruct.Bxz;
                Byz         = obj.boundaryStruct.Byz;
                Byx         = obj.boundaryStruct.Byx;
                Bzx         = obj.boundaryStruct.Bzx;
                Bzy         = obj.boundaryStruct.Bzy;
                Px          = obj.boundaryStruct.Px;
                Py          = obj.boundaryStruct.Py;
                Pz          = obj.boundaryStruct.Pz;
                Fx          = obj.boundaryStruct.Fx;
                Fy          = obj.boundaryStruct.Fy;
                Fz          = obj.boundaryStruct.Fz;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DEFINITION OF THE WEAK BC IMPOSITION %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                delta = 1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % IMPOSITION OF THE BC IN THE GLOBAL STIFFNESS MATRIX %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %% IMPOSE VELOCITY PROFILE ALONG X DIRECTION
                
                % (1) DIRICHLET - DIRICHLET
                
                if (strcmp(bcInfTagUx,'dir') && strcmp(bcOutTagUx,'dir'))

                    for imb = 1 : numbModesUx
                        
                        % VELOCITY COMPONENT IN X
                        
                        Axx((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Axx((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        Axx(imb * numbCtrlPtsUx , :) = zeros(size(Axx(imb * numbCtrlPtsUx , :)));
                        
                        Axx((imb - 1) * numbCtrlPtsUx + 1 , (imb - 1) * numbCtrlPtsUx + 1) = delta;
                        Axx(imb * numbCtrlPtsUx , imb * numbCtrlPtsUx) = delta;
                        
                        % VELOCITY COUPLING IN Y
                        
                        Bxy((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Bxy((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        Bxy(imb * numbCtrlPtsUx , :) = zeros(size(Bxy(imb * numbCtrlPtsUx , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Bxz((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Bxz((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        Bxz(imb * numbCtrlPtsUx , :) = zeros(size(Bxz(imb * numbCtrlPtsUx , :)));
                        
                        % PRESSURE COMPONENT IN X
                        
                        Px((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Px((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        Px(imb * numbCtrlPtsUx , :) = zeros(size(Px(imb * numbCtrlPtsUx , :)));
                        
                        % FORCING COMPONENT IN X
                        
                        Fx((imb - 1) * numbCtrlPtsUx + 1) = delta * projInfBCUx(imb);
                        Fx(imb * numbCtrlPtsUx) = delta * projOutBCUx(imb);

                    end
                    
                % (2) DIRICHLET - NEUMANN

                elseif ( strcmp(bcInfTagUx,'dir') && strcmp(bcOutTagUx,'neu') )
                    
                    for imb = 1 : numbModesUx
                        
                        % VELOCITY COMPONENT IN X
                        
                        Axx((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Axx((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        Axx((imb - 1) * numbCtrlPtsUx + 1 , (imb - 1) * numbCtrlPtsUx + 1) = delta;
                        
                        % VELOCITY COUPLING IN Y
                        
                        Bxy((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Bxy((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Bxz((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Bxz((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        
                        % PRESSURE COMPONENT IN X
                        
                        Px((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Px((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        
                        % FORCING COMPONENT IN X
                        
                        Fx((imb - 1) * numbCtrlPtsUx + 1) = delta * projInfBCUx(imb);
                        Fx(imb * numbCtrlPtsUx) = Fx(imb * numbCtrlPtsUx) + projOutBCUx(imb);
                        
                    end

                % (3) NEUMANN - DIRICHLET
                
                elseif ( strcmp(bcInfTagUx,'neu') && strcmp(bcOutTagUx,'dir') )

                    for imb = 1 : numbModesUx
                        
                        % VELOCITY COMPONENT IN X
                        
                        Axx(imb * numbCtrlPtsUx , :) = zeros(size(Axx(imb * numbCtrlPtsUx , :)));
                        Axx(imb * numbCtrlPtsUx , imb * numbCtrlPtsUx) = delta;

                        % VELOCITY COUPLING IN Y
                        
                        Bxy(imb * numbCtrlPtsUx , :) = zeros(size(Bxy(imb * numbCtrlPtsUx , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Bxz(imb * numbCtrlPtsUx , :) = zeros(size(Bxz(imb * numbCtrlPtsUx , :)));
                        
                        % PRESSURE COMPONENT IN X
                        
                        Px(imb * numbCtrlPtsUx , :) = zeros(size(Px(imb * numbCtrlPtsUx , :)));
                        
                        % FORCING COMPONENT IN X
                        
                        Fx((imb - 1) * numbCtrlPtsUx + 1) = Fx((imb - 1) * numbCtrlPtsUx + 1) + projInfBCUx(imb);
                        Fx(imb * numbCtrlPtsUx) = delta * projOutBCUx(imb);
                        
                    end
                    
                % (4) ROBIN - ROBIAN & ROBIN - NEUMANN
                    
                elseif ( strcmp(bcInfTagUx,'rob') && (strcmp(bcOutTagUx,'rob') || (strcmp(bcOutTagUx,'neu') ) ) )

                    for imb = 1 : numbModesUx
                        
                        % FORCING COMPONENT IN X
                        
                        Fx((imb - 1) * numbCtrlPtsUx + 1) = Fx((imb - 1) * numbCtrlPtsUx + 1) + projInfBCUx(imb);
                        Fx(imb * numbCtrlPtsUx) = Fx(imb * numbCtrlPtsUx) + projOutBCUx(imb);
                        
                    end
                    
                % (5) NEUMANN - NEUMANN
                    
                elseif ( strcmp(bcInfTagUx,'neu') && (strcmp(bcOutTagUx,'neu') ) )

                    for imb = 1 : numbModesUx
                        
                        % FORCING COMPONENT IN X
                        
                        Fx((imb - 1) * numbCtrlPtsUx + 1) = Fx((imb - 1) * numbCtrlPtsUx + 1) + projInfBCUx(imb);
                        Fx(imb * numbCtrlPtsUx) = Fx(imb * numbCtrlPtsUx) + projOutBCUx(imb);
                        
                    end
                    
                % DIRICHLET - ROBIN
                    
                elseif ( strcmp(bcInfTagUx,'dir') && strcmp(bcOutTagUx,'rob') )
                    
                    for imb = 1 : numbModesUx
                        
                        % VELOCITY COMPONENT IN X
                        
                        Axx((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Axx((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        Axx((imb - 1) * numbCtrlPtsUx + 1 , (imb - 1) * numbCtrlPtsUx + 1) = delta;
                        
                        % VELOCITY COUPLING IN Y
                        
                        Bxy((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Bxy((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Bxz((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Bxz((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        
                        % PRESSURE COMPONENT IN X
                        
                        Px((imb - 1) * numbCtrlPtsUx + 1 , :) = zeros(size(Px((imb - 1) * numbCtrlPtsUx + 1 , :)));
                        
                        % FORCING COMPONENT IN X
                        
                        Fx((imb - 1) * numbCtrlPtsUx + 1) = delta * projInfBCUx(imb);
                        Fx(imb * numbCtrlPtsUx) = Fx(imb * numbCtrlPtsUx) + projOutBCUx(imb);
                        
                    end
                
                end
                
                %% IMPOSE VELOCITY PROFILE ALONG Y DIRECTION
                
                % (1) DIRICHLET - DIRICHLET
                
                if (strcmp(bcInfTagUy,'dir') && strcmp(bcOutTagUy,'dir'))

                    for imb = 1 : numbModesUy
                        
                        % VELOCITY COMPONENT IN Y
                        
                        Ayy((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Ayy((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        Ayy(imb * numbCtrlPtsUy , :) = zeros(size(Ayy(imb * numbCtrlPtsUy , :)));
                        
                        Ayy((imb - 1) * numbCtrlPtsUy + 1 , (imb - 1) * numbCtrlPtsUy + 1) = delta;
                        Ayy(imb * numbCtrlPtsUy , imb * numbCtrlPtsUy) = delta;
                        
                        % VELOCITY COUPLING IN X
                        
                        Byx((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Byx((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        Byx(imb * numbCtrlPtsUy , :) = zeros(size(Byx(imb * numbCtrlPtsUy , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Byz((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Byz((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        Byz(imb * numbCtrlPtsUy , :) = zeros(size(Byz(imb * numbCtrlPtsUy , :)));
                        
                        % PRESSURE COMPONENT IN Y
                        
                        Py((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Py((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        Py(imb * numbCtrlPtsUy , :) = zeros(size(Py(imb * numbCtrlPtsUy , :)));
                        
                        % FORCING COMPONENT IN Y
                        
                        Fy((imb - 1) * numbCtrlPtsUy + 1) = delta * projInfBCUy(imb);
                        Fy(imb * numbCtrlPtsUy) = delta * projOutBCUy(imb);

                    end
                    
                % (2) DIRICHLET - NEUMANN

                elseif ( strcmp(bcInfTagUy,'dir') && strcmp(bcOutTagUy,'neu') )
                    
                    for imb = 1 : numbModesUy
                        
                        % VELOCITY COMPONENT IN Y
                        
                        Ayy((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Ayy((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        Ayy((imb - 1) * numbCtrlPtsUy + 1 , (imb - 1) * numbCtrlPtsUy + 1) = delta;
                        
                        % VELOCITY COUPLING IN X
                        
                        Byx((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Byx((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Byz((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Byz((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        
                        % PRESSURE COMPONENT IN Y
                        
                        Py((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Py((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        
                        % FORCING COMPONENT IN Y
                        
                        Fy((imb - 1) * numbCtrlPtsUy + 1) = delta * projInfBCUy(imb);
                        Fy(imb * numbCtrlPtsUy) = Fy(imb * numbCtrlPtsUy) + projOutBCUy(imb);
                        
                    end

                % (3) NEUMANN - DIRICHLET
                
                elseif ( strcmp(bcInfTagUy,'neu') && strcmp(bcOutTagUy,'dir') )

                    for imb = 1 : numbModesUy
                        
                        % VELOCITY COMPONENT IN Y
                        
                        Ayy(imb * numbCtrlPtsUy , :) = zeros(size(Ayy(imb * numbCtrlPtsUy , :)));
                        Ayy(imb * numbCtrlPtsUy , imb * numbCtrlPtsUy) = delta;

                        % VELOCITY COUPLING IN X
                        
                        Byx(imb * numbCtrlPtsUy , :) = zeros(size(Byx(imb * numbCtrlPtsUy , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Byz(imb * numbCtrlPtsUy , :) = zeros(size(Byz(imb * numbCtrlPtsUy , :)));
                        
                        % PRESSURE COMPONENT IN Y
                        
                        Py(imb * numbCtrlPtsUy , :) = zeros(size(Py(imb * numbCtrlPtsUy , :)));
                        
                        % FORCING COMPONENT IN Y
                        
                        Fy((imb - 1) * numbCtrlPtsUy + 1) = Fy((imb - 1) * numbCtrlPtsUy + 1) + projInfBCUy(imb);
                        Fy(imb * numbCtrlPtsUy) = delta * projOutBCUy(imb);
                        
                    end
                    
                % (4) ROBIN - ROBIAN & ROBIN - NEUMANN
                    
                elseif ( strcmp(bcInfTagUy,'rob') && (strcmp(bcOutTagUy,'rob') || (strcmp(bcOutTagUy,'neu') ) ) )

                    for imb = 1 : numbModesUy
                        
                        % FORCING COMPONENT IN Y
                        
                        Fy((imb - 1) * numbCtrlPtsUy + 1) = Fy((imb - 1) * numbCtrlPtsUy + 1) + projInfBCUy(imb);
                        Fy(imb * numbCtrlPtsUy) = Fy(imb * numbCtrlPtsUy) + projOutBCUy(imb);
                        
                    end
                    
                % (5) NEUMANN - NEUMANN
                    
                elseif ( strcmp(bcInfTagUy,'neu') && (strcmp(bcOutTagUy,'neu') ) )

                    for imb = 1 : numbModesUy
                        
                        % FORCING COMPONENT IN Y
                        
                        Fy((imb - 1) * numbCtrlPtsUy + 1) = Fy((imb - 1) * numbCtrlPtsUy + 1) + projInfBCUy(imb);
                        Fy(imb * numbCtrlPtsUy) = Fy(imb * numbCtrlPtsUy) + projOutBCUy(imb);
                        
                    end
                    
                % DIRICHLET - ROBIN
                    
                elseif ( strcmp(bcInfTagUy,'dir') && strcmp(bcOutTagUy,'rob') )
                    
                    for imb = 1 : numbModesUy
                        
                        % VELOCITY COMPONENT IN Y
                        
                        Ayy((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Ayy((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        Ayy((imb - 1) * numbCtrlPtsUy + 1 , (imb - 1) * numbCtrlPtsUy + 1) = delta;
                        
                        % VELOCITY COUPLING IN X
                        
                        Byx((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Byx((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        
                        % VELOCITY COUPLING IN Z
                        
                        Byz((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Byz((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        
                        % PRESSURE COMPONENT IN Y
                        
                        Py((imb - 1) * numbCtrlPtsUy + 1 , :) = zeros(size(Py((imb - 1) * numbCtrlPtsUy + 1 , :)));
                        
                        % FORCING COMPONENT IN Y
                        
                        Fy((imb - 1) * numbCtrlPtsUy + 1) = delta * projInfBCUy(imb);
                        Fy(imb * numbCtrlPtsUy) = Fy(imb * numbCtrlPtsUy) + projOutBCUy(imb);
                        
                    end
                
                end
                
                %% IMPOSE VELOCITY PROFILE ALONG Z DIRECTION
                
                % (1) DIRICHLET - DIRICHLET
                
                if (strcmp(bcInfTagUz,'dir') && strcmp(bcOutTagUz,'dir'))

                    for imb = 1 : numbModesUz
                        
                        % VELOCITY COMPONENT IN Z
                        
                        Azz((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Azz((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        Azz(imb * numbCtrlPtsUz , :) = zeros(size(Azz(imb * numbCtrlPtsUz , :)));
                        
                        Azz((imb - 1) * numbCtrlPtsUz + 1 , (imb - 1) * numbCtrlPtsUz + 1) = delta;
                        Azz(imb * numbCtrlPtsUz , imb * numbCtrlPtsUz) = delta;
                        
                        % VELOCITY COUPLING IN X
                        
                        Bzx((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Bzx((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        Bzx(imb * numbCtrlPtsUz , :) = zeros(size(Bzx(imb * numbCtrlPtsUz , :)));
                        
                        % VELOCITY COUPLING IN Y
                        
                        Bzy((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Bzy((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        Bzy(imb * numbCtrlPtsUz , :) = zeros(size(Bzy(imb * numbCtrlPtsUz , :)));
                        
                        % PRESSURE COMPONENT IN Z
                        
                        Pz((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Pz((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        Pz(imb * numbCtrlPtsUz , :) = zeros(size(Pz(imb * numbCtrlPtsUz , :)));
                        
                        % FORCING COMPONENT IN Z
                        
                        Fz((imb - 1) * numbCtrlPtsUz + 1) = delta * projInfBCUz(imb);
                        Fz(imb * numbCtrlPtsUz) = delta * projOutBCUz(imb);

                    end
                    
                % (2) DIRICHLET - NEUMANN

                elseif ( strcmp(bcInfTagUz,'dir') && strcmp(bcOutTagUz,'neu') )
                    
                    for imb = 1 : numbModesUz
                        
                        % VELOCITY COMPONENT IN Z
                        
                        Azz((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Azz((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        Azz((imb - 1) * numbCtrlPtsUz + 1 , (imb - 1) * numbCtrlPtsUz + 1) = delta;
                        
                        % VELOCITY COUPLING IN X
                        
                        Bzx((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Bzx((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        
                        % VELOCITY COUPLING IN Y
                        
                        Bzx((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Bzx((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        
                        % PRESSURE COMPONENT IN Z
                        
                        Pz((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Pz((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        
                        % FORCING COMPONENT IN Z
                        
                        Fz((imb - 1) * numbCtrlPtsUz + 1) = delta * projInfBCUz(imb);
                        Fz(imb * numbCtrlPtsUz) = Fz(imb * numbCtrlPtsUz) + projOutBCUz(imb);
                        
                    end

                % (3) NEUMANN - DIRICHLET
                
                elseif ( strcmp(bcInfTagUz,'neu') && strcmp(bcOutTagUz,'dir') )

                    for imb = 1 : numbModesUz
                        
                        % VELOCITY COMPONENT IN Z
                        
                        Azz(imb * numbCtrlPtsUz , :) = zeros(size(Azz(imb * numbCtrlPtsUz , :)));
                        Azz(imb * numbCtrlPtsUz , imb * numbCtrlPtsUz) = delta;

                        % VELOCITY COUPLING IN X
                        
                        Bzx(imb * numbCtrlPtsUz , :) = zeros(size(Bzx(imb * numbCtrlPtsUz , :)));
                        
                        % VELOCITY COUPLING IN Y
                        
                        Bzx(imb * numbCtrlPtsUz , :) = zeros(size(Bzx(imb * numbCtrlPtsUz , :)));
                        
                        % PRESSURE COMPONENT IN Z
                        
                        Pz(imb * numbCtrlPtsUz , :) = zeros(size(Pz(imb * numbCtrlPtsUz , :)));
                        
                        % FORCING COMPONENT IN Z
                        
                        Fz((imb - 1) * numbCtrlPtsUz + 1) = Fz((imb - 1) * numbCtrlPtsUz + 1) + projInfBCUz(imb);
                        Fz(imb * numbCtrlPtsUz) = delta * projOutBCUz(imb);
                        
                    end
                    
                % (4) ROBIN - ROBIAN & ROBIN - NEUMANN
                    
                elseif ( strcmp(bcInfTagUz,'rob') && (strcmp(bcOutTagUz,'rob') || (strcmp(bcOutTagUz,'neu') ) ) )

                    for imb = 1 : numbModesUz
                        
                        % FORCING COMPONENT IN Z
                        
                        Fz((imb - 1) * numbCtrlPtsUz + 1) = Fz((imb - 1) * numbCtrlPtsUz + 1) + projInfBCUz(imb);
                        Fz(imb * numbCtrlPtsUz) = Fz(imb * numbCtrlPtsUz) + projOutBCUz(imb);
                        
                    end
                    
                % (5) NEUMANN - NEUMANN
                    
                elseif ( strcmp(bcInfTagUz,'neu') && (strcmp(bcOutTagUz,'neu') ) )

                    for imb = 1 : numbModesUz
                        
                        % FORCING COMPONENT IN Z
                        
                        Fz((imb - 1) * numbCtrlPtsUz + 1) = Fz((imb - 1) * numbCtrlPtsUz + 1) + projInfBCUz(imb);
                        Fz(imb * numbCtrlPtsUz) = Fz(imb * numbCtrlPtsUz) + projOutBCUz(imb);
                        
                    end
                    
                % DIRICHLET - ROBIN
                    
                elseif ( strcmp(bcInfTagUz,'dir') && strcmp(bcOutTagUz,'rob') )
                    
                    for imb = 1 : numbModesUz
                        
                        % VELOCITY COMPONENT IN Z
                        
                        Azz((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Azz((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        Azz((imb - 1) * numbCtrlPtsUz + 1 , (imb - 1) * numbCtrlPtsUz + 1) = delta;
                        
                        % VELOCITY COUPLING IN X
                        
                        Bzx((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Bzx((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        
                        % VELOCITY COUPLING IN Y
                        
                        Bzx((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Bzx((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        
                        % PRESSURE COMPONENT IN Z
                        
                        Pz((imb - 1) * numbCtrlPtsUz + 1 , :) = zeros(size(Pz((imb - 1) * numbCtrlPtsUz + 1 , :)));
                        
                        % FORCING COMPONENT IN Z
                        
                        Fz((imb - 1) * numbCtrlPtsUz + 1) = delta * projInfBCUz(imb);
                        Fz(imb * numbCtrlPtsUz) = Fz(imb * numbCtrlPtsUz) + projOutBCUz(imb);
                        
                    end
                
                end
                
                %% EXPORT UPDATE DATA STRUCTURES
                
                probDataStruct.Axx = Axx;
                probDataStruct.Ayy = Ayy;
                probDataStruct.Azz = Azz;
                probDataStruct.Bxy = Bxy;
                probDataStruct.Bxz = Bxz;
                probDataStruct.Byz = Byz;
                probDataStruct.Byx = Byx;
                probDataStruct.Bzx = Bzx;
                probDataStruct.Bzy = Bzy;
                probDataStruct.Px  = Px;
                probDataStruct.Py  = Py;
                probDataStruct.Pz  = Pz;
                probDataStruct.Fx  = Fx;
                probDataStruct.Fy  = Fy;
                probDataStruct.Fz  = Fz;
                  
            end
    end
end