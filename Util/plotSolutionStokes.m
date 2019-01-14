function [] = plotSolutionStokes(plotStruct)

    import Core.AssemblerADRHandler
    import Core.AssemblerStokesHandler
    import Core.AssemblerNSHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler

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
    
    evalNodesY     = linspace(0,1,M);
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

    [modalBasisP,modalBasisDerP] = newModalBasisStokes(objModalBasis);
    
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

    [modalBasisUx,modalBasisDerUx] = newModalBasisStokes(objModalBasis);
    
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

    [modalBasisUy,modalBasisDerUy] = newModalBasisStokes(objModalBasis);

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
    
    AUX = plotStruct.solStruct;

    index1 = numbControlPtsUx * plotStruct.discStruct.numbModesUx;
    index2 = numbControlPtsUy * plotStruct.discStruct.numbModesUy;
    index3 = numbControlPtsP  * plotStruct.discStruct.numbModesP;

    uUx = AUX(1:index1,1);
    uUy = AUX(index1 + 1 : index1 + index2,1);
    uP  = AUX(index1 + index2 + 1 : index1 + index2 + index3,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % SOLUTION PROJECTION %
    %%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Project the solution vector of the pressure using the
    % information of the isogeometric basis and the modal basis.
    %---------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COORDINATES OF THE PROJECTION IN THE PHYSICAL DOMAIN %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PROJECTION POINTS
    
    XP  = mapOut(evalNodesXP,evalNodesY,plotStruct.map,1);
    YP  = mapOut(evalNodesXP,evalNodesY,plotStruct.map,2);

    XUx = mapOut(evalNodesXUx,evalNodesY,plotStruct.map,1);
    YUx = mapOut(evalNodesXUx,evalNodesY,plotStruct.map,2);

    XUy = mapOut(evalNodesXUy,evalNodesY,plotStruct.map,1);
    YUy = mapOut(evalNodesXUy,evalNodesY,plotStruct.map,2);

    % EXTENDED PROJECTION POINTS

    Nx = 100;
    Ny = 100;

    x_eval = linspace(0,1,Nx);
    y_eval = linspace(0,1,Ny);
    Xeval  = mapOut(x_eval,y_eval,plotStruct.map,1);
    Yeval  = mapOut(x_eval,y_eval,plotStruct.map,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROJECTION USING THE ISOGEOMETRIC BASIS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for m = 1:plotStruct.discStruct.numbModesP

        currUP = uP;
        spaceP    = plotStruct.spaceP;
        geoP = plotStruct.refDomain1D;
        approxP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP)  = sp_eval(currUP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP), spaceP, geoP, structEvalNodesXP);
        dapproxP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP) = sp_eval(currUP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP), spaceP, geoP, structEvalNodesXP,'gradient');

    end

    for m = 1:plotStruct.discStruct.numbModesUx

        currUUx = uUx;
        spaceUx = plotStruct.spaceUx;
        geoUx   = plotStruct.refDomain1D;
        approxUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx)  = sp_eval(currUUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx), spaceUx, geoUx, structEvalNodesXUx);
        dapproxUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx) = sp_eval(currUUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx), spaceUx, geoUx, structEvalNodesXUx,'gradient');

    end

    for m = 1:plotStruct.discStruct.numbModesUy

        currUUy = uUy;
        spaceUy = plotStruct.spaceUy;
        geoUy   = plotStruct.refDomain1D;
        approxUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy)  = sp_eval(currUUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy), spaceUy, geoUy, structEvalNodesXUy);
        dapproxUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy) = sp_eval(currUUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy), spaceUy, geoUy, structEvalNodesXUy,'gradient');

    end

    igaProjP  = approxP';
    igaProjUx = approxUx';
    igaProjUy = approxUy';

    igaProjP_dx  = dapproxP';
    igaProjUx_dx = dapproxUx';
    igaProjUy_dx = dapproxUy';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATRIX STRUCTURE FOR THE PROJECTED INFORMATION %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solMatP   = zeros(M,numbControlPtsP);
    solMatUx  = zeros(M,numbControlPtsUx);
    solMatUy  = zeros(M,numbControlPtsUy);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROJECTION USING THE MODAL BASIS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for kk = 1:M

        % PROJECTION OF THE PRESSURE

        for hh = 1:numbControlPtsP
            for mb = 1:plotStruct.discStruct.numbModesP
                solP(hh,kk) = solP(hh,kk) + igaProjP(hh + (mb-1) * numbControlPtsP) * modalBasisP(kk,mb);
            end
        end

        % PROJECTION OF Ux

        for hh = 1:numbControlPtsUx
            for mb = 1:plotStruct.discStruct.numbModesUx
                solUx(hh,kk) = solUx(hh,kk) + igaProjUx(hh + (mb-1) * numbControlPtsUx) * modalBasisUx(kk,mb);
            end
        end

        % PROJECTION OF Uy

        for hh = 1:numbControlPtsUy
            for mb = 1:plotStruct.discStruct.numbModesUy
                solUy(hh,kk) = solUy(hh,kk) + igaProjUy(hh + (mb-1) * numbControlPtsUy) * modalBasisUy(kk,mb);
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE THE SOLUTION MATRIX %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    finalSolP  = solP;
    finalSolUx = solUx;
    finalSolUy = solUy;

    for ii = 1:M
        solMatP(ii,:)  = finalSolP((ii-1)  * numbControlPtsP  + 1 : ii * numbControlPtsP );
        solMatUx(ii,:) = finalSolUx((ii-1) * numbControlPtsUx + 1 : ii * numbControlPtsUx);
        solMatUy(ii,:) = finalSolUy((ii-1) * numbControlPtsUy + 1 : ii * numbControlPtsUy);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROJECTION IN FINE MESH FOR THE PLOTTING %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : This step is necessary because during the previous projection
    % step we could not choose the number of points to plot the solution
    % along the centerline, only along the transverse direction.
    %---------------------------------------------------------------------%

    higaSolP  = griddata(XP ,YP ,solMatP ,Xeval,Yeval);
    higaSolUx = griddata(XUx,YUx,solMatUx,Xeval,Yeval);
    higaSolUy = griddata(XUy,YUy,solMatUy,Xeval,Yeval);

    [numbX,numbY] = size(higaSolP);

    for ii = 1:numbX
        for jj = 1:numbY
            if(isnan(higaSolP(ii,jj)))
                higaSolP(ii,jj) = 0;
            end
            if(isnan(higaSolUx(ii,jj)))
                higaSolUx(ii,jj) = 0;
            end
            if(isnan(higaSolUy(ii,jj)))
                higaSolUy(ii,jj) = 0;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STORE SOLUTION STRUCTURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solStruct.P  = higaSolP;
    solStruct.Ux = higaSolUx;
    solStruct.Uy = higaSolUy;
    solStruct.X  = Xeval;
    solStruct.Y  = Yeval;

    %%%%%%%%%%%%%%%%%%%%%%%
    % CHECK EXPORT FOLDER %
    %%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Create the export folder to save the plots and video of the
    % simulation.
    %---------------------------------------------------------------------%
    
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
    %---------------------------------------------------------------------%
    % Note : Open the current export folder.
    %---------------------------------------------------------------------%
    
    cd(fileName);

    %%%%%%%%%%%%%%%%%%%%
    % CREATE FIGURE PLOT
    %%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Create figure to store the plots and exporting format.
    %---------------------------------------------------------------------%
        
    X  = solStruct.X;
    Y  = solStruct.Y;
    P  = solStruct.P;
    Ux = solStruct.Ux;
    Uy = solStruct.Uy;

    fig1 = figure('visible','off');
    quiver(X,Y,Ux,Uy);

    fig2 = figure('visible','off');
    [~,~] = contourf(X,Y,P,20);    
    colormap(jet);
    cmin = min(min(P));
    cmax = max(max(P));
    caxis([cmin cmax]);
    colorbar();
    minX = min(min(X));
    maxX = max(max(X));
    minY = min(min(Y));
    maxY = max(max(Y));
    axis([minX maxX minY maxY]);
    axis equal
    scale = 0.1;
    daspect([1 1 scale]);
    pbaspect([1 1 scale]);
    set(gca, 'FontSize', 14);
    hold on;

    fig3 = figure('visible','off');
    [~,~] = contourf(X,Y,Ux,20);    
    colormap(jet);
    cmin = min(min(Ux));
    cmax = max(max(Ux));
    caxis([cmin cmax]);
    colorbar();
    minX = min(min(X));
    maxX = max(max(X));
    minY = min(min(Y));
    maxY = max(max(Y));
    axis([minX maxX minY maxY]);
    axis equal
    scale = 0.1;
    daspect([1 1 scale]);
    pbaspect([1 1 scale]);
    set(gca, 'FontSize', 14);
    hold on;

    fig4 = figure('visible','off');
    [~,~] = contourf(X,Y,Uy,20);    
    colormap(jet);
    cmin = min(min(Uy));
    cmax = max(max(Uy));
    caxis([cmin cmax]);
    colorbar();
    minX = min(min(X));
    maxX = max(max(X));
    minY = min(min(Y));
    maxY = max(max(Y));
    axis([minX maxX minY maxY]);
    axis equal
    scale = 0.1;
    daspect([1 1 scale]);
    pbaspect([1 1 scale]);
    set(gca, 'FontSize', 14);
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT INFLOW AND OUTFLOW PROFILES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numbPts = 12;
    
    fig5 = figure('visible','off');
    plot(Ux(linspace(1,Ny,numbPts),1),...
        Y(linspace(1,Ny,numbPts),1),...
        'o',...
        'LineWidth',2,...
        'MarkerSize',5,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor',[1.0,1.0,1.0]);
    hold on
    plot(plotStruct.boundCondStruct.bc_inf_data_Ux(Y(:,1)),...
        Y(:,1),...
        '--.',...
        'LineWidth',2,...
        'MarkerSize',8,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1.0,1.0,1.0]);
    l = legend({'$\tilde{u}_{x}|_{\Gamma_{out}}$','$u_{x}|_{\Gamma_{out}}$'},'FontSize',14);
    set(l, 'Interpreter', 'latex')
    
    
    fig6 = figure('visible','off');
    plot(Ux(linspace(1,Ny,numbPts),end),...
        Y(linspace(1,Ny,numbPts),end),...
        'o',...
        'LineWidth',2,...
        'MarkerSize',5,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor',[1.0,1.0,1.0]);
    hold on
    plot(plotStruct.boundCondStruct.bc_inf_data_Ux(Y(:,end)),...
        Y(:,end),...
        '--.',...
        'LineWidth',2,...
        'MarkerSize',8,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1.0,1.0,1.0]);
    l = legend({'$\tilde{u}_{x}|_{\Gamma_{out}}$','$u_{x}|_{\Gamma_{out}}$'},'FontSize',14);
    set(l, 'Interpreter', 'latex')

    %%%%%%%%%%%%%
    % SAVE FIGURE
    %%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Save the current figure.
    %---------------------------------------------------------------------%

    print(fig1,'Quiver','-dpng');
    print(fig2,'Contour_P','-dpng');
    print(fig3,'Contour_Ux','-dpng');
    print(fig4,'Contour_Uy','-dpng');
    print(fig5,'Inflow_Profile','-dpng');
    print(fig6,'Outflow_Profile','-dpng');
    
    cd ..
    
    %%%%%%%%%%%%%%%%%
    % COMPUTE ERROR %
    %%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : This function only plots the HigaMod solution and its
    % gradient. We do not compute the error with respect to a analytical or
    % reference solution, so the L2 and H1 norms of the error are set to
    % zero.
    %---------------------------------------------------------------------%

end
