function [Wom] = Womersley(yy,tt,womStruct)
    
    %%%%%%%%%%%%%
    % PARAMTERS %
    %%%%%%%%%%%%%
    
    % INFLOW PRESSURE
    
    Pinf = womStruct.Pinf;
    
    % OUTFLOW PRESSURE
    
    Pout = womStruct.Pout;
    
    % PRESSURE DROP
    
    dP = Pinf - Pout;
    
    % VISCOSITY
    
    mu = womStruct.mu;
    
    % FLUID DENSITY
    
    rho = womStruct.rho;
    
    % DYNAMIC VISCOSITY
    
    nu = mu/rho;
        
    % LENGTH OF THE CYLINDER
    
    L = womStruct.L;
    
    % RADIUS OF THE CYLINDER
    
    R = L/2;
    
    % PULSATING ANGULAR FREQUENCY
    
    w = (2 * pi)/womStruct.finalTime;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WOMERSLEY INFLOW PROFILE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % NUMBER OF FOURIER COMPONENTS
    
    numbModes = womStruct.numbModes;
    
    % WOMERSLEY SIDE PARAMETER
    
    sigma = nu / (L^2);
    
    % WOMERSLEY AUXILIARY FUNCTIONS
    
    dk = @(y,t,kk) 4 * dP ./ (rho .* pi .* (2 * kk + 1) .* (((2 .* kk + 1).^4) .* (sigma^2) * (pi^4) + w^2));
    Sk = @(y,t,kk) sin((pi/L) * (2 .* kk + 1) .* y);
    Pk = @(y,t,kk) ((2 .* kk + 1).^2) .* sigma .* (pi^2) .* sin(w .* t) - w .* cos(w .* t) + w .* exp(-((2 .* kk + 1).^2) .* sigma .* (pi^2) .* t);
    
    % WOMERSLEY INFLOW PROFILE
    
    Wom = zeros(length(yy),length(tt));
    
    for ii = 1:length(yy)
        for jj = 1:length(tt)
            Wom(ii,jj) = sum(dk(yy(ii),tt(jj),1:numbModes) .* Sk(yy(ii),tt(jj),1:numbModes) .* Pk(yy(ii),tt(jj),1:numbModes));
        end
    end

end