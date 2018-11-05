clear all
close all

%% DEFINITION TEST CASE
dt = 2*pi/50;
T = 0 : dt : 2*pi;
x = sin(T);
y = cos(T);
param = ones(size(T));

numbSamples = length(T);
numbStates  = 3;
numbEnsembles = 10;
numbMeasState = 2;

alpha = 1;

M = [1 dt 0; -dt alpha 0; 0 0 1];
H = [1 0 0; 0 1 0];

solutionMatrix = [x;y;param];

initState = [0;0;10];

% Check the evolution matrix

checkEvol = zeros(numbStates,numbSamples);
checkEvol(:,1) = initState;

for ii = 1:numbSamples-1
    checkEvol(:,ii+1) = M * checkEvol(:,ii);
end

figure
plot(checkEvol(1,:),checkEvol(2,:));

%% ENSEMBLE PROPERTIES
%---------------------------------------------------------%
% Note: Definition of the properties of the random field to
% be determined to create the samples of the ensemble.
%---------------------------------------------------------%

% TYPE
%---------------------------------------------------------%
% Note: The type of the random field can be set as Gaussian
% ('gauss'), Exponential ('exp') or Turbulent
% ('turbulent').
%---------------------------------------------------------%

corrProperties.name = 'gauss';

% CORRELATION PARAMETER
%---------------------------------------------------------%
% Note: The scaling parameters for the correlation
% function. 'c0' c0 may be a scalar for isotropic
% correlation or a vector for an anisotropic correlation.
% In the anisotropic case, the vector must have
% 'n'elements, where 'n' is the dimension of the state
% vector.
%
% (1) 'c0' =~ 0 (Uncorrelated samples)
% (2) 'c0' =~ 1 (Completly correlated samples)
%---------------------------------------------------------%

corrProperties.c0 = 0.0001;

% SAMPLE VARIANCE
%---------------------------------------------------------%
% Note: The variance 'sigma' of the samples. May be a 
% scalar or a vector the size of the state vector.
%---------------------------------------------------------%

corrProperties.sigma = 10^-1.5;

% DATA POINTS
%---------------------------------------------------------%
% Note: Data points where the samples will be evaluated.
% Linear vector going form one to the number of states in
% the state vector.
%---------------------------------------------------------%

dataPts = linspace(1,numbStates,numbStates)';

%% EnKF - FILTER ALGORITHM
%---------------------------------------------------------%
% 1. Initialization
% 2. Create state measurements
% 3. Create initial ensemble set for state
% 4. Create initial ensemble set for parameters
% 5. Filtering loop
%   5.1. Loop over time iterations
%   5.2. Loop over ensemble set
%       5.2.1. Compute the system matrix with previous
%       assmilated state
%       5.2.2. Compute the prediction for the current
%       emsemble element
%       5.2.3. Assimilate the measurements with the
%       prediction for the current emsemble element
%---------------------------------------------------------%

tic;

% 1. INITIALIZATION

Z = zeros(numbMeasState,numbSamples);
ensObservation = [];
ensObsNoise = [];
currStateMean = zeros(numbStates,1);
ensObsNoiseMean = zeros(numbMeasState,1);
currObsMean = zeros(numbMeasState,1);
normalizedState = [];
normalizedObs = [];
ensAssmilatedState = [];
ensPredictionState = [];
assStateHistory = [];

% 2. CREATE STATE MEASUREMENTS

for jj = 1:numbSamples
    measureStd = 0.05 * max(H * solutionMatrix(:,jj));
    v = measureStd * randn(numbMeasState,1);
    Z(:,jj) = H * solutionMatrix(:,jj) + v;

end

% 3. CREATE INITIAL ENSEMBLE OF STATES 
%---------------------------------------------------------%
% Note: I need to create the object property for the number
% of ensemble points to be created.
%---------------------------------------------------------%

[stateEnsembleSet,~] = randomfield(corrProperties,dataPts,'nsamples',numbEnsembles,'mean',initState);

figure; plot(x,y,'--'); hold on;
plot(Z(1,:),Z(2,:),'+');

figure; plot(x,y,'--'); hold on;
for ii = 1:numbEnsembles
    plot(stateEnsembleSet(1,:),stateEnsembleSet(2,:),'*');
end

% 4. CREATE INITIAL ENSEMBLE SET OF PARAMETERS
%---------------------------------------------------------%
% Note: This step is necessary when there is parameter
% identification of the system.
%---------------------------------------------------------%

% 5. FILTERING LOOP

figure
plot(x,y,'--');
hold on;

for jj = 1:numbSamples
    
    plot(stateEnsembleSet(1,:),stateEnsembleSet(2,:),'*');

    %-----------------------------------------------------%
    % Step 1 : Draw a statistically consistent observation
    % set with the current observation
    %-----------------------------------------------------%

    for ii = 1:numbEnsembles

        if (jj == 1)                            
            measureStd = 0.05 * max(H * solutionMatrix(:,jj));
            v = measureStd * randn(numbMeasState,1);

            ensObsNoise(:,ii) = v;

            ensObservation(:,ii) = Z(:,jj) + v;
        else                        
            measureStd = 0.05 * max(H * solutionMatrix(:,jj));
            v = measureStd * randn(numbMeasState,1);

            ensObsNoise(:,ii) = v;

            ensObservation(:,ii) = Z(:,jj) + v;
        end
    end

    %-----------------------------------------------------%
    % Step 2 : Compute the ensemble means for the state,
    % expected observation and noise used to generate the
    % observation set
    %-----------------------------------------------------%

    %-----------------------------------------------------%
    % Note: Do not forget that the matrix containing the
    % current set of state ensemble is update at each
    % iteration of the filter.
    %-----------------------------------------------------%

    for ii = 1:numbEnsembles
        currStateMean = currStateMean + stateEnsembleSet(:,ii);
        ensObsNoiseMean = ensObsNoiseMean + ensObsNoise(:,ii);
        currObsMean = currObsMean + H * stateEnsembleSet(:,ii);
    end

    currStateMean = currStateMean/numbEnsembles;
    ensObsNoiseMean = ensObsNoiseMean/numbEnsembles;
    currObsMean = currObsMean/numbEnsembles;

    %-----------------------------------------------------%
    % Step 3 : Compute the normalized anomalies
    %-----------------------------------------------------%

    for ii = 1:numbEnsembles

        normalizedState(:,ii) = (stateEnsembleSet(:,ii) - currStateMean)/(sqrt(numbEnsembles - 1));
        normalizedObs(:,ii) = (H * stateEnsembleSet(:,ii) - ensObsNoise(:,ii) - currObsMean + ensObsNoiseMean)/(sqrt(numbEnsembles - 1));

    end

    %-----------------------------------------------------%
    % Step 4 : Compute the filter gain
    %-----------------------------------------------------%

    gainInv = (normalizedObs * normalizedObs')\eye(size(normalizedObs * normalizedObs'));
    gainEnKF = (normalizedState * normalizedObs') * gainInv;

    %-----------------------------------------------------%
    % Step 5 : Update the ensemble set using the weighted
    % innovation
    %-----------------------------------------------------%

    for ii = 1:numbEnsembles

        ensAssmilatedState(:,ii) = stateEnsembleSet(:,ii) + gainEnKF * (ensObservation(:,ii) - H * stateEnsembleSet(:,ii));

    end

    %-----------------------------------------------------%
    % Step 6 : Compute the new system matrix using the
    % current state of the system. Necessary to be
    % performed if a parameter of the system is used as
    % additional state variable.
    %-----------------------------------------------------%
    
    

    %-----------------------------------------------------%
    % Step 7 : Compute the ensemble forecast using the
    % current dynamic law for the system
    %-----------------------------------------------------%

    for ii = 1:numbEnsembles
        
        alpha = ensAssmilatedState(3,ii);
        M = [1 dt 0; -dt alpha 0; 0 0 1];

        ensPredictionState(:,ii) = M * ensAssmilatedState(:,ii);

    end

    stateEnsembleSet = ensPredictionState;
    assStateHistory = [assStateHistory currStateMean];
    
    

end

% Debug
figure; plot(assStateHistory(1,:),assStateHistory(2,:),'-o'); hold on;
plot(solutionMatrix(1,:),solutionMatrix(2,:),'--');

filterTime = toc;

display = sprintf('Finished Filtering (EnKF without Parameter Identification)');
disp(display);

return