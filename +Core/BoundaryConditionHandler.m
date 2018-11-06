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
            
            function [liftCoeffA,liftCoeffB] = computeFourierCoeff(obj)

                infBC = obj.infBoundCond;
                outBC = obj.outBoundCond;
                
                valueInfBC = 
                    
                for ii = 1:obj.dimModalBasis
                    
                    projectInfCoeff = 
                    
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
    end
end