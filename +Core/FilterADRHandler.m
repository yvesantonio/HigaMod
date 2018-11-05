classdef FilterADRHandler
    
    %% FILTER ADR HANDLER CLASS
    % Class description and documentation.
    
    properties (SetAccess = public, GetAccess = public)
        
        %% FILTER ADR HANDLER - OBJECT PROPERTIES
        % The properties of the objects used in the FilterADRhandler
        % encapsulate all of the variables needed to run the methods
        % bellow, from the support methods to the final building methods.
        
        stateMatrix;    % State matrix of the linear system (A)
        
        inputMatrix;    % Input matrix of the linear system (B)
        
        outputMatrix;   % Output matrix of the linear system (C)
        
        solutionMatrix; % Solution of the differential problem in every timestep
        
        ratioSignalNoise;   % Ratio from the signal to the white gaussian noise
                            % as a function of the spectral power of both
                            % signals.
        
    end
    
    methods (Access = public)
        
        %% FILTER ADR HANDLER - CONSTRUCT METHOD
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
        
        %% FILTER ADR HANDLER - BUILDING METHODS
        % The following functions correspond to the individual scripts
        % previously developed of the functions '...'.
        
            %% Method 'observabilityCheck'

            function [decomposedA,decomposedB,decomposedC,transformMatrix,observedStates] = observabilityCheck(obj)

                %%
                % observabilityCheck - This function checks if the input
                %                      system is observable. If the system is 
                %                      not observable, then the fucntion
                %                      returns the system matrix decomposed
                %                      on the respective observable and non
                %                      -observable parts.
                % 
                % Note:
                % All of the following inputs are encapsulated in the
                % object properties.
                %
                % The inputs are:
                %%
                %   (1)  stateMatrix      : State matrix of the linear
                %                           system. (A)
                %   (2)  inputMatrix      : Input matrix of the linear
                %                           system. (B)
                %   (3)  outputMatrix     : Output matrix of the linear
                %                           system. (C)
                %
                % The outputs are:
                %%
                %   (1) decomposedA      : Decomposed State Matrix for the
                %                          linear system
                %   (2) decomposedB      : Decomposed Input Matrix for the
                %                          linear system
                %   (3) decomposedC      : Decomposed Output Matrix for the
                %                          linear system
                %   (3) transformMatrix  : Linear transformation matrix
                %                          used to find the decomposed
                %                          system matrices
                %   (4) observedStates   : Vector containing the map of
                %                          observable states

                % COMPUTES THE OBSERVABILITY MATRIX
                
                observMatrix = obsv(obj.stateMatrix,obj.outputMatrix);
                
                % COMPUTES THE NUMBER OF NON OBSERVABLE STATE
                
                %---------------------------------------------------------%
                % Note:
                % If the number of non observable states is greater than
                % zero, then the system matrices need to be decomposed into
                % its observable and non-observable components.
                %---------------------------------------------------------%
                
                nonObservStates = length(obj.stateMatrix) - rank(observMatrix);
                
                % COMPUTES THE DECOMPOSED MATRICES
                
                if(nonObservStates > 0)
                    
                    [decomposedA,decomposedB,decomposedC,transformMatrix,observedStates] = ...
                     obsvf(obj.stateMatrix,obj.inputMatrix,obj.outputMatrix);
                    
                    disp('The system is not observable and needed to be decomposed');
                    
                else
                    
                    decomposedA = obj.stateMatrix;
                    decomposedB = obj.inputMatrix;
                    decomposedC = obj.outputMatrix;
                    
                    disp('The system is observable and did not need to be decomposed');
                    
                end

            end
            
            %% Method 'noisyMeasurement'
            
            function [noisyMeasure] = noisyMeasurement(obj)
                
                [~,numbIterations] = size(obj.solutionMatrix);
                
                noisyMeasure = zeros(size(obj.solutionMatrix));
                noisyMeasure(:,1) = obj.solutionMatrix(:,1);
                
                for ii = 1:numbIterations - 1
                    
                    noisyMeasure(:,iteration + 1) = awgn(obj.solutionMatrix(:,iteration + 1),obj.ratioSignalNoise);
    
                end
                
            end
            
    end
end