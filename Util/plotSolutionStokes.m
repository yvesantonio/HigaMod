function [] = plotSolutionStokes(plotStruct)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINTION HORIZONTAL NODES AND INTERVALS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numbKnotsP      = plotStruct.numbKnotsP;
    numbKnotsUx     = plotStruct.numbKnotsUx;
    numbKnotsUy     = plotStruct.numbKnotsUx;
    
    numbControlPtsP  = plotStruct.numbControlPtsP;
    numbControlPtsUx = plotStruct.numbControlPtsUx;
    numbControlPtsUy = plotStruct.numbControlPtsUy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINITION NUMBER OF VERTICAL EVALUATION POINTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M = 1000;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINITION VERTICAL EVALUATION POINTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    evalNodesY = linspace(0,1,M);
    numbEvalNodesY = length(evalNodesY);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALIZATION OF THE SOLUTION MATRICES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    solP   = zeros(numbControlPtsP ,numbEvalNodesY);
    solUx  = zeros(numbControlPtsUx,numbEvalNodesY);
    solUy  = zeros(numbControlPtsUy,numbEvalNodesY);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % PRESSURE MODAL BASIS %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    import Core.BasisHandler

    objModalBasis = BasisHandler();

    objModalBasis.dimModalBasis         = plotStruct.discStruct.numbModesP;
    objModalBasis.evalNodesY            = evalNodesY;
    objModalBasis.labelUpBoundCond      = plotStruct.boundaryStruct.boundCondStruc.bc_up_tag_P;
    objModalBasis.labelDownBoundCond    = plotStruct.boundaryStruct.boundCondStruc.bc_down_tag_P;
    objModalBasis.coeffForm             = plotStruct.probParameters;

    [modalBasisP,modalBasisDerP] = newModalBasis(objModalBasis);
    
    %%%%%%%%%%%%%%%%%%
    % Ux MODAL BASIS %
    %%%%%%%%%%%%%%%%%%
    
    import Core.BasisHandler

    objModalBasis = BasisHandler();

    objModalBasis.dimModalBasis         = plotStruct.discStruct.numbModesUx;
    objModalBasis.evalNodesY            = evalNodesY;
    objModalBasis.labelUpBoundCond      = plotStruct.boundaryStruct.boundCondStruc.bc_up_tag_Ux;
    objModalBasis.labelDownBoundCond    = plotStruct.boundaryStruct.boundCondStruc.bc_down_tag_Ux;
    objModalBasis.coeffForm             = plotStruct.probParameters;

    [modalBasisUx,modalBasisDerUx] = newModalBasis(objModalBasis);
    
    %%%%%%%%%%%%%%%%%%
    % Uy MODAL BASIS %
    %%%%%%%%%%%%%%%%%%
    
    import Core.BasisHandler

    objModalBasis = BasisHandler();

    objModalBasis.dimModalBasis         = plotStruct.discStruct.numbModesUy;
    objModalBasis.evalNodesY            = evalNodesY;
    objModalBasis.labelUpBoundCond      = plotStruct.boundaryStruct.boundCondStruc.bc_up_tag_Uy;
    objModalBasis.labelDownBoundCond    = plotStruct.boundaryStruct.boundCondStruc.bc_down_tag_Uy;
    objModalBasis.coeffForm             = plotStruct.probParameters;

    [modalBasisUy,modalBasisDerUy] = newModalBasis(objModalBasis);

    %---------------------------------------------------------------------%
    % Note:
    % The coefficients 'coeffModalBase' and 'coeffModalBaseDer' are
    % vectors conrresponding respectively to the coefficients of the
    % modal bases and derivative of the modal basis on the points
    % assigned in the Y direction.
    %---------------------------------------------------------------------%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRESSURE HORIZONTAL MESH %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    structEvalNodesXP = {linspace(0,1,numbControlPtsP)};
    evalNodesXP       = linspace(0,1,numbControlPtsP);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Ux HORIZONTAL MESH %
    %%%%%%%%%%%%%%%%%%%%%%
    
    structEvalNodesXUx = {linspace(0,1,numbControlPtsUx)};
    evalNodesXUx       = linspace(0,1,numbControlPtsUx);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Uy HORIZONTAL MESH %
    %%%%%%%%%%%%%%%%%%%%%%
    
    structEvalNodesXUy = {linspace(0,1,numbControlPtsUy)};
    evalNodesXUy       = linspace(0,1,numbControlPtsUy);

    %---------------------------------------------------------------------%
    % Note:
    % The function 'bspeval' used in the following loops is taken from
    % the NURBS toolbox and its primary objective is to evaluate a 
    % univariate B-Spline function in the desired nodes.
    %
    % The parameter involved in this functions are (in order):
    %
    % (1ST INPUT) Degree of the B-Spline;
    % (2ND INPUT) Control Points, matrix of size (dim,nc), where 'dim' 
    %             is the degree of the B-Spline function and 'nc' is
    %             the number of control points.
    % (3RD INPUT) Knot sequence, row vector of size nk, where 'nk' is
    %             the number of knot points.
    % (4TH INPUT) Parametric evaluation points, row vector of size nu,
    %             where 'nu' is the number of points used to represent
    %             the domain.
    %
    % (OUTPUT) Evaluated points, matrix of size (dim,nu).
    %---------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT SOLUTION VECTOR %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    uP  = cell(length(plotStruct.timeStruct.timeDomain),1);
    uUx = cell(length(plotStruct.timeStruct.timeDomain),1);
    uUy = cell(length(plotStruct.timeStruct.timeDomain),1);
    
    for tt = 1:length(plotStruct.timeStruct.timeDomain)
        
        AUX = plotStruct.solStruct{tt};
        
        index1 = numbControlPtsUx * plotStruct.discStruct.numbModesUx;
        index2 = numbControlPtsUy * plotStruct.discStruct.numbModesUy;
        index3 = numbControlPtsP  * plotStruct.discStruct.numbModesP;
        
        uUx{tt} = AUX(1:index1,1);
        uUy{tt} = AUX(index1 + 1 : index1 + index2,1);
        uP{tt}  = AUX(index1 + index2 + 1 : index1 + index2 + index3,1);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % PRESSURE PROJECTION %
    %%%%%%%%%%%%%%%%%%%%%%%
    
    switch caso
            
        case {1,2,3,4,5,6,7,8,9,10}
        
            for m = 1:plotStruct.discStruct.numbModesP
                
                approx( (m-1) * numbControlPtsP+1:m*numbControlPtsP)  = sp_eval(u((m-1)*numbControlPtsP+1:m*numbControlPtsP), space, geometry, structEvalNodesXP);
                dapprox((m-1)*numbControlPtsP+1:m*numbControlPtsP) = sp_eval(u((m-1)*numbControlPtsP+1:m*numbControlPtsP), space, geometry, structEvalNodesXP,'gradient');

            end

            % COMPUTATION OF THE NURBS BASE

            u  = approx';

            % COMPUTATION OF THE DERIVATIVE OF THE NURBS BASE

            ux = dapprox';
    end
    
    solVect = zeros(numbControlPtsP * M,1);
    solMat = zeros(M,numbControlPtsP);

    for h = 1:numbControlPtsP

        for k = 1:M

            for imb = 1:size_mb

                % COMPUTATION OF THE APPROXIMATED SOLUTION VECTOR FIELD
                % EVALUATED IN THE POINTS OF THE DOMAIN MESH

                solP(h,k)   = solP(h,k)   + u(h+(imb-1)*numbControlPtsP)*modalBasisP(k,imb);
                solVect(h + (k-1)*numbControlPtsP) = solVect(h + (k-1)*numbControlPtsP) + u(h+(imb-1)*numbControlPtsP)*modalBasisP(k,imb);

                % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
                % SOLUTION VECTOR FIELD EVALUATED IN THE POINTS OF THE
                % DOMAIN MESH

                sol_x(h,k) = sol_x(h,k) + ux(h+(imb-1)*numbControlPtsP)*modalBasisP(k,imb);
                sol_y(h,k) = sol_y(h,k) + u(h+(imb-1)*numbControlPtsP)*modalBasisDerP(k,imb);

            end

        solP(h,k) = solP(h,k)+ a_ril(1) * evalNodesY(k) + b_ril(1);

        end

    end
    
    finalSol = solP;
    for ii = 1:M
        solMat(ii,:) = finalSol((ii-1)*numbControlPtsP + 1 : ii * numbControlPtsP);
    end

    X = mapOut(evalNodesXP,evalNodesY,map,1);
    Y = mapOut(evalNodesXP,evalNodesY,map,2);
    
    Nx = 103;
    Ny = 103;
    
    x_eval = linspace(0,1,Nx);
    y_eval = linspace(0,1,Ny);
    Xeval  = mapOut(x_eval,y_eval,map,1);
    Yeval  = mapOut(x_eval,y_eval,map,2);
    
    higaSol = griddata(X,Y,solMat,Xeval,Yeval);
    [numbX,numbY] = size(higaSol);
    
    for ii = 1:numbX
        for jj = 1:numbY
            if(isnan(higaSol(ii,jj)))
                higaSol(ii,jj) = 0;
            end
        end
    end

    %--------------------------------------------------%
    % CONTOUR PLOT OF THE EXACT SOLUTION IN THE DOMAIN
    %--------------------------------------------------%
    
    minX = min(min(Xeval));
    minY = min(min(Yeval));
    maxX = max(max(Xeval));
    maxY = max(max(Yeval));
    
    scale = 0.01;
    
    % CONTOUR PLOT FOR THE SOLUTION
    
%     % SURF PLOT OF THE SOLUTION
%     
%     figure;
%     mesh(Xeval,Yeval,higaSol);
%     colormap(jet);
%     cmin = min(min(higaSol));
%     cmax = max(max(higaSol));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY cmin cmax]);
%     axis equal
%     daspect([1 1 scale])
%     pbaspect([1 1 scale])
%     set(gca, 'FontSize', 14)
%     
%     % SURF PLOT OF THE CUTS IN THE SOLUTION
%     
%     xCut1 = (Nx-1)/2 + 1;   % Position of the centerline
%     xCut2 = (Nx-3)/4 + 1;   % Position of the first quarter line
%     xCut3 = (Nx-3)*3/4 + 2; % Position of the third quarter line
%     
%     solX1 = higaSol(xCut1,:);   % Solution at the centerline
%     solX2 = higaSol(xCut2,:);   % Solution at the first quarter line
%     solX3 = higaSol(xCut3,:);   % Solution at the third quarter line
    
%     figure;
%     plot(Xeval(xCut1,:),solX1,'b','LineWidth',2); hold on;
% %     plot(Xeval(xCut2,:),solX2,'r','LineWidth',2); 
% %     plot(Xeval(xCut3,:),solX3,'g','LineWidth',2);
%     cmin = 0.02; % min(min(higaSol));
%     cmax = 0.02; % max(max(higaSol));
%     axis([minX maxX cmin cmax]);
%     axis equal
%     set(gca, 'FontSize', 14)
%     grid on
% %     daspect([1 scale scale])
% %     pbaspect([1 scale scale])

    %%%%%%%%%%%%%%%%%%%%%
    % CHECK EXPORT FOLDER
    %%%%%%%%%%%%%%%%%%%%%
    
    for ii = 1:1000

        checkFolder = exist(['MatlabPlots',num2str(ii)]);

        if (checkFolder == 7)
        else
            fileName = ['MatlabPlots',num2str(ii-1)];
            break
        end
    end

    %%%%%%%%%%%%%%%%%%%%%
    % OPEN CURRENT FOLDER 
    %%%%%%%%%%%%%%%%%%%%%
    
    cd(fileName);

    %%%%%%%%%%%%%%%%%%%%
    % CREATE FIGURE PLOT
    %%%%%%%%%%%%%%%%%%%%
    
    fig1 = figure('visible','off');
            
    [~,~] = contourf(Xeval,Yeval,higaSol,20);    
    colormap(jet);
    cmin = min(min(higaSol));
    cmax = max(max(higaSol));
    caxis([cmin cmax])
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    daspect([1 1 scale])
    pbaspect([1 1 scale])
    set(gca, 'FontSize', 14)
    hold on;

    %%%%%%%%%%%%%
    % SAVE FIGURE
    %%%%%%%%%%%%%
    
    print(fig1,['Plot_At_t=',num2str(iteration)],'-dpng')
    
    cd ..

    errL2 = 0;
    errH1 = 0;

end
