classdef IntegrateHandler
    
    %% INTEGRATE HANDLER CLASS
    % The IntegrateHandler is a class that contains all of the
    % scripts responsible for the integration and projection functions used
    % in the script of the solver applied. This function are used through
    % out the other classes but were not added to the utility class because
    % they gather a unique functionality.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% INTEGRATE HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the
        % IntegrateHandler encapsulate all of the variables needed
        % to run the methods associated with the integration and projection
        % of values throught out the solution of the differential problem.
        % Differently from the AssemblerADRHandler class, here the
        % methods are completily independdent from each other and be
        % directly accessed.
        
            %% Integrate Properties

            funcToIntegrate;     % Function to be Integrated
            funcWeight;          % Weight

            %% Gauss-Legendre Properties

            numbQuadNodes;       % Number of Desired Quadrature nodes

            %% Compute Eigenvalues Properties

            funcStourmLiou;     % Function relative to the Stourm-Liouville
                                % Problem
                          
            numbEigenvalues;    % Desired Number of Eigenvalues to be Computed
            
            labelUpBoundCond;   % Contains the Label Identifying the Nature of
                                % the Boundary Conditions on the Upper Limit of
                                % the Domain
                          
            labelDownBoundCond; % Contains the Label Identifying the Nature of
                                % the Boundary Conditions on the Lower Limit of
                                % the Domain
                          
            coeffForm;          % Data Strusture Containing All the @-Functions
                                % and the Constants Relative to the Bilinear Form

            %% Legendre Polynomial Properties

            degreePolyLegendre; % Degree of the Legendre Polynomial

            %% Qudrature Rule Properties
            
            leftBoundInterval;  % Left Bound of the New Interval
            
            rightBoundInterval; % Right Bound of the New Interval
            
            inputNodes;         % Input Nodes to be Shifted
            
            inputWeights;       % Input Weights to be Scaled

    end
    
    methods (Access = public)
        
        %% INTEGRATION HANDLER - CONSTRUCT METHOD
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
        
        %% INTEGRATION HANDLER - FUNCTIONAL METHODS
        % The following functions correspond to the individual scripts
        % previously developed of the functions 'computeEigenvalues',
        % 'gaussLegendre', 'polyLegendre', 'quadratureRule' and
        % 'integrate'.
        
            %% Method 'computeEigenvalues'
            
            function lambda = computeEigenvalues(obj)

                %%
                % computeEigenvalues - This function computes the first 'm' eigenvalues
                %                      of the Stourm-Liouville problem associated with
                %                      the chosen boundary conditions.
                %
                % The inputs are:
                %%
                %   (1)  funcStourmLiou     : Function relative to the Stourm-Liouville
                %                            Problem
                %   (2)  numbEigenvalues    : Desired Number of Eigenvalues to be Computed
                %   (3)  labelUpBoundCond   : Contains the Label Identifying the Nature of
                %                             the Boundary Conditions on the Upper Limit of
                %                             the Domain
                %   (4)  labelDownBoundCond : Contains the Label Identifying the Nature of
                %                             the Boundary Conditions on the Lower Limit of
                %                             the Domain
                %   (5)  coeffForm          : Data Strusture Containing All the @-Functions
                %                             and the Constants Relative to the Bilinear Form
                %
                % The outputs are:
                %%
                %   (1) lambda        : Vector Containing thefirst 'm' eigenvalues
                %                       of the Stourm-Liouville problem associated with
                %                       the chosen boundary conditions

                % Initialization of the Eigenvalue Vector

                lambda = zeros( obj.numbEigenvalues, 1);

                % Comuptation of the Eigenvalues for the Implemented Boundary
                % Conditions

                switch [obj.labelUpBoundCond,obj.labelDownBoundCond]

                    %-----------------------------------------%
                    % Dirichlet-Dirichlet Boundary Conditions %
                    %-----------------------------------------%

                    case 'dirdir'
                        for i= 1 : obj.numbEigenvalues
                            lambda(i) = i*pi;
                        end

                    %-----------------------------------------------------%
                    % Robin-Dirichlet/Dirichlet-Robin Boundary Conditions %
                    %-----------------------------------------------------%

                    case {'dirrob','robdir'}
                        for i = 1 : obj.numbEigenvalues
                            a = -pi/2 + i*pi + 1e-2;
                            b =  -1e-2 + i*pi + pi/2;           % Attention!
                            lambda( i ) = fzero(obj.funcStourmLiou, [a, b]);

                        end

                    %---------------------------------%
                    % Robin-Robin Boundary Conditions %
                    %---------------------------------%

                    case 'robrob'
                        i=1;

                        %-------------------------------------------------------------%
                        % Attention:
                        % Why 'mu' was not considered inside the conditions for the
                        % 'if'?
                        %-------------------------------------------------------------%

                        if(obj.coeffForm.coeffrobin <= pi/2)

                        sv = 1e-10;
                        a = sv;
                        b = pi/2-sv;
                        found = false;
                        
                            while(~found && a < 1)
                                    if obj.funcStourmLiou(a) * obj.funcStourmLiou(b)<0

                                        [lambda(i)] = fzero(obj.funcStourmLiou,[a,b]);
                                        found = true;

                                    end
                                a = a * 10;
                                b = b * 10;
                            end			
                            i=i+1;
                        end
                        
                        j=1;
                        
                        while(i <= obj.numbEigenvalues)
                            
                            a = -pi/2 + j*pi + 1e-8;
                            b =  -1e-8 + j*pi + pi/2;
                            
                            if obj.funcStourmLiou(a) * obj.funcStourmLiou(b)>0
                                ac = (a + b)/2;
                                bc = ac;

                                lambda(i) = fzero(obj.funcStourmLiou,[a,bc]);               
                                i = i+1;				
                                lambda(i) = fzero(obj.funcStourmLiou,[ac,b]);
                                
                            else			
                                
                                lambda(i) = fzero(obj.funcStourmLiou,[a,b]);
                                
                            end
                            
                            i=i+1;
                            j=j+1;
                        end   

                    otherwise
                        disp('In computeEigenvalues: Boundary Conditions Not Recognized')
                end
            end
            
            %% Method 'gaussLegendre'
            
            function [numbObtainedNodes,quadNodes,wquadNodesWeight] = gaussLegendre(obj) 

                %%
                % gaussLegendre - This function implements the one-dimensional
                %                 Gauss-Legendre quadrature rule. It returns the 'NP'
                %                 quadrature nodes ('QNODES') and the 'NP' weights
                %                 ('w') of the quadrature rule in the interval (0,1).
                %
                % The inputs are:
                %%
                %   (1) numbQuadNodes   : Number of Desired Quadrature nodes
                %
                % The outputs are:
                %%
                %   (1) numbObtainedNodes   : Number of Quadrature Nodes
                %   (2) quadNodes           : Quadrature Nodes
                %   (3) quadNodesWeight     : Wights of the Quadratures Nodes

                % Definition of the Number of Quadrature Nodes

                numbObtainedNodes = obj.numbQuadNodes;

                % Computation of the Quadrature Nodes

                if obj.numbQuadNodes <= 1    
                   x = 1/2;
                   wquadNodesWeight = 1;    
                   return 
                end

                jac = zeros(numbObtainedNodes);   
                k = (1:obj.numbQuadNodes - 1);

                v   = (k)./(sqrt(4*(k.^2)-1)); 
                jac = jac+diag(v,1)+diag(v,-1); 
                [wquadNodesWeight,x] = eig(jac); 
                norm2 = sqrt(diag(wquadNodesWeight'*wquadNodesWeight));
                wquadNodesWeight = (2*wquadNodesWeight(1,:)'.^2)./norm2;
                x = diag(x); 

                for k = 1:obj.numbQuadNodes - 1 
                    
                   l = k;
                   value_a = x(k);    
                   
                   for j = k + 1:obj.numbQuadNodes   
                       
                      if x(j) < value_a          
                         l = j;          
                         value_a = x(j);          
                      end      
                      
                   end
                   
                   if l ~= k       
                      t = x(k);  x(k) = x(l);  x(l) = t;       
                      t = wquadNodesWeight(k);  wquadNodesWeight(k) = wquadNodesWeight(l);  wquadNodesWeight(l) = t;       
                   end    
                   
                end 

                quadNodes = 1/2*(x+1);
                wquadNodesWeight = wquadNodesWeight/2;
                
                return
            end
            
            %% Method 'integrate'
            
            function [integral] = integrate(obj) 
                
                %%
                % integrate - This function integrates the function 'F'
                %                 with weoght 'w'
                % The inputs are:
                %%
                %   (1) funcToIntegrate     : Function to be Integrated
                %   (2) funcWeight          : Weight
                %
                % The outputs are:
                %%
                %   (1) integral  : Integral of F with Weight 'w'
                %---------------------------------------------------------------------%

                  integral = obj.funcToIntegrate * obj.funcWeight;
            end
                
            %% Method 'polyLegendre'
            
            function [polyLegendre,polyLegendreDer] = polyLegendre(obj)

                %%
                % polyLegendre   - This function creates a Legendre Polynomial of
                %                  degree equal to 'm'.
                %
                % The inputs are:
                %%
                %   (1) degreePolyLegendre  : Degree of the Legendre Polynomial
                %
                % The outputs are:
                %%
                %   (1) polyLegendre        : Legendre Polynomial
                %   (2) polyLegendreDer     : Derivative of the Legendre Polynomial
                %---------------------------------------------------------------------%

                polyLegendre = cell(obj.degreePolyLegendre,1);
                polyLegendreDer = cell(obj.degreePolyLegendre,1);

                polyLegendre{1} = 1;
                polyLegendreDer{1} = polyder(polyLegendre{1});

                polyLegendre{2} = [1 0 ];
                polyLegendreDer{2} = polyder(polyLegendre{2});

                for ii = 1:obj.degreePolyLegendre - 2

                    polyLegendre{ii+2} = 1/(ii+1)*((2*ii+1)*[polyLegendre{ii+1} 0] - ii * [0 0 polyLegendre{ii}]);
                    polyLegendreDer{ii+2} = polyder(polyLegendre{ii+2});

                end

                for ii = 1:obj.degreePolyLegendre

                    C = sqrt((2*(ii-1)+1)/2);
                    polyLegendre{ii} = polyLegendre{ii}*C;
                    polyLegendreDer{ii} = polyLegendreDer{ii}*C;

                end
            end
            
            %% Method 'quadratureRule'
            
            function [shiftedNodes,scaledWeights] = quadratureRule(obj)

                %%
                % quadratureRule  - This function updates a 1D quadrature formula from
                %                   the interval (0,1) to the new interval (a,b), by 
                %                   shifting the nodes and scaling the weights 
                %                   accordingly.
                %
                % The inputs are:
                %%
                %   (1) leftBoundInterval   : Left Bound of the New Interval
                %   (2) rightBoundInterval  : Right Bound of the New Interval
                %   (3) inputNodes          : Input Nodes to be Shifted
                %   (4) inputWeights        : Input Weights to be Scaled
                %
                % The outputs are:
                %%
                %   (1) shiftedNodes        : Shifted Nodes from the Quadrature Rule
                %   (2) scaledWeights       : Scaled Weights from the Quadrature Rule
                %---------------------------------------------------------------------%

                h = obj.rightBoundInterval - obj.leftBoundInterval;
                shiftedNodes = h * obj.inputNodes + obj.leftBoundInterval;
                scaledWeights = h * obj.inputWeights;

            return
            end
    end
end