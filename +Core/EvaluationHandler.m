classdef EvaluationHandler
    
    %% EVALUATION HANDLER CLASS
    % The EvaluationHandler is a class that contains all of the scripts
    % that auxiliate the computation of the solution for the differential
    % problems proposed in the examples. It calls all of the other basic
    % classes and functions already incorporated.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% EVALUATION HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the EvaluationHandler
        % encapsulate all of the variables needed to run the methods 
        % that auxiliate the computation of the solution for the
        % differential problems proposed. They basically represent all of
        % the variables required to define complitily the equations, domain
        % and boundary conditions of the differential problem.
        
            %% Evaluation Properties

            leftLimitDomainX;       % Left Limit of the Domain in the X Direction

            rightLimitDomainX;       % Right Limit of the Domain in the X Direction

            leftLabelBoundCond;     % Boundary Conditions on the Left Bound of the
                                    % Domain

            rightLabelBoundCond;    % Boundary Conditions on the Right Bound of the
                                    % Domain

            leftDomainBoundary;     % Left Boundary of the Domain

            rightDomainBoundary;    % Right Boundary of the Domain

            matrixA;                % Bloc Matrix Defining the Problem

            vectorB;                % Bloc Vector Defining the Problem   

            solWithoutDirichlet;    % Solution without the Dirichlet Contributions
            
            basisProject;           % Basis in which the functions will be projected
            
            funcProject;            % Function to be projected
            
            coeffWeights;           % Relative Weight of the Basis Coefficients of the
                                    % Projected Function
                          
            solDiffProblem;         % Solution of the Differential Problem
            
            numbSolNodes;           % Vector of Nodes of the Mesh
            
            numbQuadNodes;          % Number of Quadrature Nodes
            
            modalBasis;             % Modal Basis HigaMod
            
            stepMeshX;              % Vector Containing the Step of the Finite
                                    % Element Mesh
                          
            degreePolySplineBasis;  % Degree of the Polynomial B-Spline Basis
            
            continuityParameter;    % Degree of Continuity of the Basis 'C^(p-k)'
            
            dimModalBasis;          % Dimension of the Modal Basis
            
            coefficientForm;        % Data Strusture Containing All the @-Functions
                                    % and the Constants Relative to the Bilinear Form
                          
            robinCondStruct;        % Data Structure Containing the Two Values of the
                                    % Coefficients (R, L) for the Robin Condition Used
                                    % in the Domain Decomposition

    end
    
    methods (Access = public)
        
        %% EVALUATION HANDLER - CONSTRUCT METHOD
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
        
        %% EVALUATION HANDLER - FUNCTIONAL METHODS
        % The following functions correspond to the individual scripts
        % previously developed of the functions 'evalHiMod',
        % 'extractSolutionP1', 'extractSolutionP2' and 'coeffBasisProject'.
        
            %% Method 'evalHiMod'
            
            function solVector = evalHiMod(obj)

                %%
                % evalHiMod - Evaluate the solution u of a differential probelm on HiMod 'grid' and get the 
                %             result on sol
                % The inputs are:
                %   (1) solDiffProblem  : Solution of the Differential Problem
                %   (2) modalBasis      : Modal Basis Used in the HiMod Solution
                %   (3) numbSolNodes    : Number of Nodes Used in the Solution
                %   (4) numbQuadNodes   : Number of Quadrature Nodes
                %   (5) stepMeshX       : Vector Containing the Step of the Finite
                %                         Element Mesh
                %
                % The outputs are:
                %   (1) solVector       : Vector where the value in the index are evaluated 
                %                         in accordance with HiMod nodes + modes ordering
                
                solVector = zeros(length(obj.numbSolNodes),size(obj.modalBasis,1));
                ne = length(obj.numbSolNodes)/obj.numbQuadNodes;
                nx = ne + 1;

                for h = 1:ne
                    xl = (h-1) * obj.stepMeshX;
                    xr = h * obj.stepMeshX;
                    x1 = obj.numbSolNodes((h-1) * obj.numbQuadNodes + 1);
                    x2 = obj.numbSolNodes((h-1) * obj.numbQuadNodes + 2);
                    x3 = obj.numbSolNodes((h-1) * obj.numbQuadNodes + 3);
                    x4 = obj.numbSolNodes((h-1) * obj.numbQuadNodes + 4);
                    
                    for ii = 1:size(obj.modalBasis,1)
                        for imb = 1:size(obj.modalBasis,2)

                            ul = obj.solDiffProblem(h+(imb-1)*nx);
                            uri = obj.solDiffProblem(h+1+(imb-1)*nx);

                            c1 = 1/obj.stepMeshX * (uri*(x1-xl)-ul*(x1-xr));
                            solVector((h-1) * obj.numbQuadNodes + 1,ii) = solVector((h-1) * obj.numbQuadNodes + 1,ii) + c1 * obj.modalBasis(ii,imb);

                            c2 = 1/obj.stepMeshX * (uri*(x2-xl)-ul*(x2-xr));
                            solVector((h-1) * obj.numbQuadNodes + 2,ii) = solVector((h-1) * obj.numbQuadNodes + 2,ii) + c2 * obj.modalBasis(ii,imb);

                            c3 = 1/obj.stepMeshX * (uri*(x3-xl)-ul*(x3-xr));
                            solVector((h-1) * obj.numbQuadNodes + 3,ii) = solVector((h-1) * obj.numbQuadNodes + 3,ii) + c3 * obj.modalBasis(ii,imb);

                            c4 = 1/obj.stepMeshX * (uri*(x4-xl)-ul*(x4-xr));
                            solVector((h-1) * obj.numbQuadNodes + 4,ii) = solVector((h-1) * obj.numbQuadNodes + 4,ii) + c4 * obj.modalBasis(ii,imb);
                        end
                    end
                end
            end
            
            %% Method 'extractSolutionP1'
            
            function [solWithDirichlet,leftDirContrib,rightDirContrib,leftNeuContrib,...
                      rightNeuContrib,leftRobContrib,rightRobContrib] = extractSolutionP1(obj)
    
                %%
                % extractSolutionP1  - This function adds the Dirichlet values at the
                %                      INFLOW to the final solution and computes the
                %                      quantities at the interface for the Domain
                %                      Decompositon cycle.
                %
                % The inputs are:
                %%
                %   (1)  dimModalBasis        : Dimension of the Modal Basis
                %   (2)  leftLimitDomainX     : Left Limit of the Domain in the X Direction
                %   (3)  rightLimitDomainX    : Right Limit of the Domain in the X Direction
                %   (4)  leftLabelBoundCond   : Label fo the Type of Boundary Condition at the
                %                               Left Bound of the Domain
                %   (5)  rightLabelBoundCond  : Label fo the Type of Boundary Condition at the
                %                               Right Bound of the Domain
                %   (6)  stepMeshX            : Vector Containing the Step of the Finite
                %                               Element Mesh
                %   (7)  leftDomainBoundary   : Profile of the Boundary Condition at the Left
                %                               Extreme
                %   (8)  rightDomainBoundary  : Profile of the Boundary Condition at the Right
                %                               Extreme
                %   (9)  solWithoutDirichlet  : Solution of the Differential Problem Without
                %                               the Dirichlet Contribution
                %   (10) coefficientForm            : Data Strusture Containing All the @-Functions
                %                               and the Constants Relative to the Bilinear Form
                %   (11) matrixA              : Bloc Matrix of the Differential Problem
                %   (12) vectorB              : Bloc Vaector with Known Terms
                %   (13) robinCondStruct      : Data Structure Containing the Two Values of the
                %                               Coefficients (R, L) for the Robin Condition Used
                %                               in the Domain Decomposition
                %
                % The outputs are:
                %%
                %   (1) solWithDirichlet      : Solution of the Differential Problem with the
                %                               Dirichlet Contribution
                %   (2) leftDirContrib        : Dirichlet Contribution at the Left Bound of the
                %                               Domain
                %   (3) rightDirContrib       : Dirichlet Contribution at the Right Bound of the
                %                               Domain
                %   (4) leftNeuContrib        : Neumann Contribution at the Left Bound of the
                %                               Domain
                %   (5) rightNeuContrib       : Neumann Contribution at the Right Bound of the
                %                               Domain
                %   (6) leftRobContrib        : Robin Contribution at the Left Bound of the
                %                               Domain
                %   (7) rightRobContrib       : Robin Contribution at the Right Bound of the
                %                               Domain
                %---------------------------------------------------------------------%

                % Adds the Dirichlet Values at INFLOW to the Solution and Computes the
                % Quantities at the Interface for the Domain Decomposition Loop

                leftDirContrib  = 0;
                rightDirContrib = 0;
                leftNeuContrib  = 0;
                rightNeuContrib = 0;
                leftRobContrib  = 0;
                rightRobContrib = 0;

                % Number of Intervals

                ne = round((obj.rightLimitDomainX - obj.leftLimitDomainX)/obj.stepMeshX);

                % Number of Extremes (Nodes in the X Direction)

                nx = ne+1;

                % Finite Element Mesh in the X Direction (Equispaced Nodes Considered)

                meshx    = zeros( nx, 1);
                meshx(1) = obj.leftLimitDomainX;

                for i = 2 : nx

                    meshx(i) = meshx(i-1) + obj.stepMeshX;

                end

                if ( strcmp(obj.leftLabelBoundCond, 'dir') && strcmp(obj.rightLabelBoundCond, 'dir'))

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.leftDomainBoundary( imb ); obj.solWithoutDirichlet( (imb-1)*(nx-2) + 1 : imb*(nx-2) ); obj.rightDomainBoundary( imb )];
                        
                    end

                elseif ( strcmp(obj.leftLabelBoundCond, 'dir') && strcmp(obj.rightLabelBoundCond, 'neu') )

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.leftDomainBoundary( imb ); obj.solWithoutDirichlet( (imb-1)*(nx-1) + 1 : imb*(nx-1) ) ];
                        
                    end

                elseif ( strcmp(obj.leftLabelBoundCond, 'neu') && strcmp(obj.rightLabelBoundCond, 'dir'))

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.solWithoutDirichlet( (imb-1)*(nx-1) + 1 : imb*(nx-1) ); obj.rightDomainBoundary( imb )];
                        
                    end

                elseif ( strcmp(obj.leftLabelBoundCond, 'dir') && strcmp(obj.rightLabelBoundCond, 'rob') )

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.leftDomainBoundary( imb ); obj.solWithoutDirichlet( (imb-1)*(nx-1) + 1 : imb*(nx-1) ) ];
                        
                    end

                elseif (strcmp(obj.leftLabelBoundCond, 'rob') && strcmp(obj.rightLabelBoundCond, 'rob') )
                    
                    solWithDirichlet = obj.solWithoutDirichlet;
                    
                elseif (strcmp(obj.leftLabelBoundCond, 'rob') && strcmp(obj.rightLabelBoundCond, 'neu') )
                    
                    solWithDirichlet = obj.solWithoutDirichlet;
                    
                elseif (strcmp(obj.leftLabelBoundCond, 'neu') && strcmp(obj.rightLabelBoundCond, 'neu') )
                    
                    solWithDirichlet = obj.solWithoutDirichlet;
                    
                else
                    
                    display('Error in extractSolutionP1.m');
                    
                end

                %---------------------------------------------------------------------%
                % Note:
                % The following lines of code corrispond to the previous
                % Dirichlet-Neumann Algorithm. In the current version of the code, this
                % section is never called. The code is still present in the event one
                % needs to prove the Dirichlet-Neumann Condition.
                %---------------------------------------------------------------------%

                if(~   ( strcmp(obj.leftLabelBoundCond,'rob') || strcmp(obj.leftLabelBoundCond,'rob') )    )

                    % Extraction of the Solution at the Boundary (Dirichlet and Neumann)

                    leftDirContrib = zeros( obj.dimModalBasis, 1 );
                    rightDirContrib = zeros( obj.dimModalBasis, 1 );
                    leftNeuContrib = zeros( obj.dimModalBasis, 1 );
                    rightNeuContrib = zeros( obj.dimModalBasis, 1 );

                    for imb = 1 : obj.dimModalBasis
                        
                        leftDirContrib( imb ) = solWithDirichlet ( (imb-1)*nx + 1  );
                        rightDirContrib( imb ) = solWithDirichlet ( (imb-1)*nx + nx );
                        
                    end

                    % Global Residual

                    residuo = obj.vectorB - obj.matrixA * solWithDirichlet;

                    for imb = 1 : obj.dimModalBasis
                        
                        leftNeuContrib( imb ) = - residuo( (imb-1)*nx + 1  );
                        rightNeuContrib( imb ) = - residuo( (imb-1)*nx + nx );
                        
                    end

                % Robin-Robin Case

                else

                    mu_valutato = obj.coefficientForm.mu(1,1);

                    %-----------------------------------------------------------------%
                    % Note:
                    % The previous line of code needs to be changed in the case of
                    % non constant 'mu', because otherwise the sentence loses the
                    % physical meaning.
                    %-----------------------------------------------------------------%

                    for imb= 1 : obj.dimModalBasis
                        
                        leftRobContrib(imb) = mu_valutato * (solWithDirichlet( (imb-1)*nx + 2 ) - solWithDirichlet( (imb-1)*nx + 1))/obj.stepMeshX + obj.robinCondStruct.R * solWithDirichlet( (imb-1)*nx + 1 );
                        rightRobContrib(imb) = -mu_valutato * (solWithDirichlet(imb*nx) - solWithDirichlet(imb*nx-1))/obj.stepMeshX + obj.robinCondStruct.L*solWithDirichlet(imb*nx);
                        
                    end
                end
            end
       
            %% Method 'extractSolutionP2'
           
            function [solWithDirichlet,leftDirContrib,rightDirContrib,leftNeuContrib,...
                      rightNeuContrib,leftRobContrib,rightRobContrib] = extractSolutionP2(obj)

                %%
                % extractSolutionP2 -   This function adds the Dirichlet values at the
                %                       INFLOW to the final solution and computes the
                %                       quantities at the interface for the Domain
                %                       Decompositon cycle.
                %
                % The inputs are:
                %%
                %   (1)  dimModalBasis              : Dimension of the Modal Basis
                %   (2)  leftLimitDomainX           : Left Limit of the Domain in the X Direction
                %   (3)  rightLimitDomainX          : Right Limit of the Domain in the X Direction
                %   (4)  leftLabelBoundCond         : Label fo the Type of Boundary Condition at the
                %                                     Left Bound of the Domain
                %   (5)  rightLabelBoundCond        : Label fo the Type of Boundary Condition at the
                %                                     Right Bound of the Domain
                %   (6)  stepMeshX                  : Vector Containing the Step of the Finite
                %                                     Element Mesh
                %   (7)  leftDomainBoundary         : Profile of the Boundary Condition at the Left
                %                                     Extreme
                %   (8)  rightDomainBoundary        : Profile of the Boundary Condition at the Right
                %                                     Extreme
                %   (9)  solWithoutDirichlet        : Solution of the Differential Problem Without
                %                                     the Dirichlet Contribution
                %   (10) coefficientForm                  : Data Strusture Containing All the @-Functions
                %                                     and the Constants Relative to the Bilinear Form
                %   (11) matrixA                    : Bloc Matrix of the Differential Problem
                %   (12) vectorB                    : Bloc Vaector with Known Terms
                %   (13) robinCondStruct            : Data Structure Containing the Two Values of the
                %                                     Coefficients (R, L) for the Robin Condition Used
                %                                     in the Domain Decomposition
                %   (14) degreePolySplineBasis      : Degree of the Polynomial B-Spline Basis
                %   (15) continuityParameter        : Degree of Continuity of the Basis 'C^(p-k)'
                %
                % The outputs are:
                %%
                %   (1) solWithDirichlet    : Solution of the Differential Problem with the
                %                             Dirichlet Contribution
                %   (2) leftDirContrib      : Dirichlet Contribution at the Left Bound of the
                %                             Domain
                %   (3) rightDirContrib     : Dirichlet Contribution at the Right Bound of the
                %                             Domain
                %   (4) leftNeuContrib      : Neumann Contribution at the Left Bound of the
                %                             Domain
                %   (5) rightNeuContrib     : Neumann Contribution at the Right Bound of the
                %                             Domain
                %   (6) leftRobContrib      : Robin Contribution at the Left Bound of the
                %                             Domain
                %   (7) rightRobContrib     : Robin Contribution at the Right Bound of the
                %                             Domain
                %---------------------------------------------------------------------%


                % Adds the Dirichlet Values at INFLOW to the Solution and Computes the
                % Quantities at the Interface for the Domain Decomposition Loop

                leftDirContrib  = 0;
                rightDirContrib = 0;
                leftNeuContrib  = 0;
                rightNeuContrib = 0;
                leftRobContrib  = 0;
                rightRobContrib = 0;

                % Number of Intervals

                ne = round((obj.rightLimitDomainX - obj.leftLimitDomainX)/obj.stepMeshX);

                % Number of Extremes (Nodes in the X Direction)

                nx = ne + 1;

                % Finite Element Mesh in the X Direction (Equispaced Nodes Considered)

                meshx    = zeros( nx, 1);
                meshx(1) = obj.leftLimitDomainX;

                for i = 2 : nx
                    meshx(i) = meshx(i-1) + obj.stepMeshX;
                end

                nx = ne * obj.continuityParameter + obj.degreePolySplineBasis + 1 - obj.continuityParameter;

                if ( strcmp(obj.leftLabelBoundCond, 'dir') && strcmp(obj.leftLabelBoundCond, 'dir'))

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.leftDomainBoundary( imb ); obj.solWithoutDirichlet( (imb-1)*(nx-2) + 1 : imb*(nx-2) ); obj.rightDomainBoundary( imb )];
                        
                    end

                elseif ( strcmp(obj.leftLabelBoundCond, 'dir') && strcmp(obj.rightLabelBoundCond, 'neu') )

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.leftDomainBoundary( imb ); obj.solWithoutDirichlet( (imb-1)*(nx-1) + 1 : imb*(nx-1) ) ];
                        
                    end

                elseif ( strcmp(obj.leftLabelBoundCond, 'neu') && strcmp(obj.rightLabelBoundCond, 'dir'))

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.solWithoutDirichlet( (imb-1)*(nx-1) + 1 : imb*(nx-1) ); obj.rightDomainBoundary( imb )];
                        
                    end

                elseif ( strcmp(obj.leftLabelBoundCond, 'dir') && strcmp(obj.rightLabelBoundCond, 'rob') )

                    solWithDirichlet=[];

                    for imb = 1 : obj.dimModalBasis
                        
                        solWithDirichlet = [solWithDirichlet; obj.leftDomainBoundary( imb ); obj.solWithoutDirichlet( (imb-1)*(nx-1) + 1 : imb*(nx-1) ) ];
                        
                    end

                elseif (strcmp(obj.leftLabelBoundCond, 'rob') && strcmp(obj.rightLabelBoundCond, 'rob') )
                    
                    solWithDirichlet = obj.solWithoutDirichlet;
                    
                elseif (strcmp(obj.leftLabelBoundCond, 'rob') && strcmp(obj.rightLabelBoundCond, 'neu') )
                    
                    solWithDirichlet = obj.solWithoutDirichlet;
                    
                elseif (strcmp(obj.leftLabelBoundCond, 'neu') && strcmp(obj.rightLabelBoundCond, 'neu') )
                    
                    solWithDirichlet = obj.solWithoutDirichlet;
                    
                else
                    
                    display('Error in extractSolutionP2.m');
                    
                end

                %---------------------------------------------------------------------%
                % Note:
                % The following lines of code corrispond to the previous
                % Dirichlet-Neumann Algorithm. In the current version of the code, this
                % section is never called. The code is still present in the event one
                % needs to prove the Dirichlet-Neumann Condition.
                %---------------------------------------------------------------------%

                if(~   ( strcmp(obj.leftLabelBoundCond,'rob') || strcmp(obj.rightLabelBoundCond,'rob') )    )

                    % Extraction of the Solution at the Boundary (Dirichlet and Neumann)
                    
                    leftDirContrib = zeros( obj.dimModalBasis, 1 );
                    rightDirContrib = zeros( obj.dimModalBasis, 1 );
                    leftNeuContrib = zeros( obj.dimModalBasis, 1 );
                    rightNeuContrib = zeros( obj.dimModalBasis, 1 );

                    for imb = 1 : obj.dimModalBasis
                        
                        leftDirContrib( imb ) = solWithDirichlet ( (imb-1)*nx + 1  );
                        rightDirContrib( imb ) = solWithDirichlet ( (imb-1)*nx + nx );
                        
                    end

                    % Global Residual

                    residuo = obj.vectorB - obj.matrixA * solWithDirichlet;

                    for imb = 1 : obj.dimModalBasis
                        
                        leftNeuContrib( imb ) = - residuo( (imb-1)*nx + 1  );
                        rightNeuContrib( imb ) = - residuo( (imb-1)*nx + nx );
                        
                    end

                % Robin-Robin Case

                else

                    mu_valutato = obj.coefficientForm.mu(1,1);

                    %-----------------------------------------------------------------%
                    % Note:
                    % The previous line of code needs to be changed in the case of
                    % non constant 'mu', because otherwise the sentence loses the
                    % physical meaning.
                    %-----------------------------------------------------------------%

                    for imb= 1 : obj.dimModalBasis
                        
                        leftRobContrib(imb) =  mu_valutato*(solWithDirichlet( (imb-1)*nx + 2 ) - solWithDirichlet( (imb-1)*nx + 1))/obj.stepMeshX + obj.robinCondStruct.R * solWithDirichlet( (imb-1)*nx + 1 );
                        rightRobContrib(imb) = -mu_valutato*(solWithDirichlet(imb*nx) - solWithDirichlet(imb*nx-1))/obj.stepMeshX + obj.robinCondStruct.L * solWithDirichlet(imb * nx);
                        
                    end
                end
            end

            %% Method 'coeffBasisProject'
            
            function coeffFourier = coeffBasisProject(obj)
                
                %%
                % coeffBasisProject - This function computes the Fourier Coefficients relative
                %                     to the projection of the functions into the orthogonal
                %                     basis.
                %
                % The inputs are:
                %%
                %   (1) basisProject       : Basis in which the functions will be projected
                %   (2) funcProject        : Function to be projected
                %   (3) coeffWeights       : Relative Weight of the Basis Coefficients of the
                %                            Projected Function
                %
                % The outputs are:
                %%
                %   (1) coeffFourier   : Vector Containing the 'm' coeffients relative
                %                        to the Projection of the Function into the New
                %                        Basis
                %---------------------------------------------------------------------%

                mii = size(obj.basisProject,2);

                coeffFourier = zeros(mii,1);

                for i = 1 : mii
                    
                    obj_integrate = IntegrateHandler();
                    
                    obj_integrate.funcToIntegrate = (obj.funcProject.*obj.basisProject(:,i))';
                    obj_integrate.funcWeight = obj.coeffWeights;
                    
                    coeffFourier(i) = integrate(obj_integrate);
                end
            end
            
    end
end
